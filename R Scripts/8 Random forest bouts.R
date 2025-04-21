#######################################################################################
################### BBA RANDOM FOREST: BOUT-BASED #####################################
#######################################################################################

# Fits random forest models to predict vessel following from immersion metrics on a bout-by-bout basis. 
# Implements custom k-fold cross-validation to assess model performance
# (for full dataset, k = 5; for demonstration dataset, k = 2). 
# Demonstration dataset will not produce same results.


################### CONTENTS ####################

# 1. Load packages and data
# 2. Function to get foraging bouts
# 3. Custom folding function for cross-validation
# 4. Get bouts and folds
# 5. Random forest 5: Regularity only
# 6. Random forest 6: All predictors
# 7. Random forest 7: Regularity only, uncorrected GLS
# 8. Random forest 8: Regularity only, uncorrected GLS


################## 1. LOAD PACKAGES AND DATA ##################################

library(tidyverse)
library(lubridate)
library(pbapply)
library(parallel)
library(caret)
library(randomForest)
theme_set(theme_bw())

# Set seed
set.seed(123)

# Set WD
setwd("___________")

# Read in data
imm_metrics_list1 <- readRDS("Output Dataset Files/imm_metrics_list1.RData")
imm_metrics_list1_gls <- readRDS("Output Dataset Files/imm_metrics_list1_gls.RData")







######################## 2. FUNCTION TO GET FORAGING BOUTS #################################

# Reconfigure imm_metrics_list1 dataframes to be 1 row per foraging bout
# Function takes in imm_metrics_list1 list element (dataframe) for 1 bird
# Also takes in rolling window length (in minutes) and if it should use corrected immersion (logical)
# Re-calculates immersion metrics within each bout
get.bouts <- function(imm_df, window, corrected){
  
  # Return NULL if track does not have necessary immersion vector
  if(corrected == FALSE & is.na(imm_df$leg_wet[1])){
    return(NULL)
  }
  
  # Set buffer of size window/2
  # buffer unit = rows, window unit = mins, dt = secs
  dt <- as.numeric(difftime(imm_df$DateTime[2], imm_df$DateTime[1], units = "secs"))
  buffer <- floor(window*60/(dt*2))
  
  ### STEP 1: Classify wet/dry/intermediate bouts based on window size
  
  # Identify wet and dry periods
  # Switch to using uncorrected GLS immersion if corrected == FALSE
  if(corrected == FALSE){
    imm_df$leg_wet_C <- imm_df$leg_wet
  }
  imm_df <- imm_df %>% mutate(transition = ifelse(leg_wet_C != lag(leg_wet_C), 1, 0))
  imm_df$transition[1] <- 1
  imm_df <- imm_df %>% mutate(period = cumsum(transition), P_wet = NA, imm_factor = "all dry")
  
  # Isolate wet/dry periods and iterate through them to classify as all wet, all dry, or intermediate
  imm_factor_list <- lapply(1:max(imm_df$period), function(p){
    
    # Get dataframe of a single wet or dry period
    wd_period <- imm_df[imm_df$period == p, c("leg_wet_C", "imm_factor")]
    
    # If wet/dry time period is shorter than the window length (buffer*2)
    if(nrow(wd_period) <= buffer*2){
      
      # Set whole period to "intermediate" 
      wd_period$imm_factor <- "intermediate"
      
      # Else if time period is long enough
    } else {
      
      # Identify the wet periods (imm_factor is "dry" by default)
      P_wet <- sum(wd_period$leg_wet_C == "wet")/nrow(wd_period)
      if(P_wet == 1){
        wd_period$imm_factor <- "all wet"
      }
      
      # Buffer the boundaries of each wet/dry period with "intermediate"
      wd_period$imm_factor[1:buffer] <- "intermediate"
      wd_period$imm_factor[(nrow(wd_period) - buffer + 1):nrow(wd_period)] <- "intermediate"
    }
    
    # Simplify to a vector
    return(wd_period$imm_factor)
  })
  
  # Flatten the list into single imm_factor vector, put back into imm_df
  imm_factor_vec <- unlist(imm_factor_list)
  imm_df$imm_factor <- imm_factor_vec
  
  
  ### STEP 2: Get metrics for each intermediate (foraging) bout
  
  # Add columns for start of bout (0 or 1) and bout number
  imm_df <- imm_df %>% mutate(start_bout = as.numeric(NA))
  
  # Get starts of bouts
  imm_df$start_bout <- ifelse(imm_df$imm_factor == "intermediate" & 
                             lag(imm_df$imm_factor != "intermediate"),
                           1, 0)
  if(imm_df$imm_factor[1] == "intermediate"){
    imm_df$start_bout[1] <- 1
  }
  
  # Remove non-intermediate timesteps
  # Label bouts by number
  # Recode follow vector to "N" and "F"
  imm_df2 <- imm_df %>%
    filter(imm_factor == "intermediate") %>%
    mutate(bout = cumsum(start_bout)) %>%
    mutate(follow = ifelse(is.na(follow), "N", "F"))
  
  # If using uncorrected leg_wet data, change column names to match those 
  if(corrected == FALSE){
    colnames(imm_df2) <- gsub("_gls", "", colnames(imm_df2))
  }
  
  # Iterate through each foraging bout to calculate metrics
  bouts_list <- lapply(unique(imm_df2$bout), function(bout_num){
    
    # Get df for single bout
    bout_df <- imm_df2 %>% filter(bout == bout_num)
    
    # Define landings and takeoffs (start and end of wet periods)
    bout_df <- mutate(bout_df, 
                        landing = ifelse(bout_df$leg_wet_C == "wet" & 
                                           (lag(bout_df$leg_wet_C) == "dry"), 
                                         1, 0), 
                        takeoff = ifelse(bout_df$leg_wet_C == "wet" & 
                                           (lead(bout_df$leg_wet_C) == "dry"), 
                                         1, 0))
    
    # Set first row of both to 0 (built in buffer means this has to be the case)
    #bout_df$landing[1] <- 0
    #bout_df$takeoff[1] <- 0
    
    # If any NAs left for takeoffs and landings, switch to 0 (includes first row)
    bout_df$takeoff <- ifelse(is.na(bout_df$takeoff), 0, bout_df$takeoff)
    bout_df$landing <- ifelse(is.na(bout_df$landing), 0, bout_df$landing)
    
    # Initiate vectors for durations of wet and dry
    wet_vec <- c()
    dry_vec <- c()
    dur_wet <- 0
    dur_dry <- 0
    
    # Calculate durations (in seconds) and add to vectors
    for(j in 1:nrow(bout_df)){
      
      # If end of wet period
      if (bout_df$takeoff[j] == 1){
        dur_wet <- dur_wet + dt
        wet_vec <- c(wet_vec, dur_wet)
        dur_wet <- 0
        
        # If start or middle of wet period  
      } else if (bout_df$leg_wet_C[j] == "wet"){
        dur_wet <- dur_wet + dt
        
        # If end of dry period
      } else if (bout_df$leg_wet_C[j] == "dry" & bout_df$landing[j+1] %in% c(1, NA)){
        dur_dry <- dur_dry + dt
        dry_vec <- c(dry_vec, dur_dry)
        dur_dry <- 0
        
        # If start or middle of dry period
      } else if (bout_df$leg_wet_C[j] == "dry"){
        dur_dry <- dur_dry + dt
      }
    }
    
    # Calculate mean of wet and dry durations
    mean_wet <- ifelse(!is.null(wet_vec), mean(wet_vec), 0)
    mean_dry <- ifelse(!is.null(dry_vec), mean(dry_vec), 0)
    
    # Calculate other important metrics for bouts dataframe
    P_wet <- sum(bout_df$leg_wet_C == "wet")/nrow(bout_df)
    N_landings <- sum(bout_df$landing)
    P_LWLD <- max(bout_df$P_LWLD)
    start <- min(bout_df$DateTime)
    end <- max(bout_df$DateTime)
    duration <- difftime(end, start, units = "mins")
    follow <- ifelse(sum(bout_df$follow == "F") > 0, "F", "N")
    
    # Return 1-row dataframe for this bout
    bout_row <- data.frame(Bird = bout_df$Bird[1],
                           bout = bout_num,
                           start = start,
                           end = end,
                           duration = duration,
                           mean_wet = mean_wet,
                           mean_dry = mean_dry,
                           N_landings = N_landings,
                           P_wet = P_wet,
                           P_LWLD = P_LWLD,
                           follow = follow)
    return(bout_row)
  })
  
  # Bind the rows together and return bouts dataframe
  bouts_df <- do.call(rbind, bouts_list) %>%
    mutate(follow = relevel(as.factor(follow), ref = "N"))
  return(bouts_df)
}







############### 3. CUSTOM FOLDING FUNCTION FOR CROSS-VALIDATION ###################
# Goal: Create k folds of dataset
# Ensuring that all observations from a single individual are kept in the same fold
# And each fold has roughly equal number of birds that follow vessels

# Data: Should take in all data in single dataframe
# k: Number of folds
# Output: List of folds, each containing a mutually exclusive vector of row numbers
# This output will be used for the indexOut argument of trainControl()
custom_fold <- function(data, k){
  
  # Get vectors of birds with (F) and without (N) vessel following
  birds <- unique(data$Bird)
  birds_F <- unique(data$Bird[data$follow == "F"])
  birds_N <- birds[!(birds %in% birds_F)]
  
  # Calculate how many birds in each fold
  N_birds_F <- round(length(birds_F)/k)
  N_birds_N <- round(length(birds_N)/k)
  
  # Initialise list to contain folds, and empty vector of birds already used
  folds <- vector("list", k)
  
  # Randomly allocate birds to folds
  # Ensure equal number of vessel-following birds in each fold
  for(i in 1:k){
    
    # If we're not on the last fold
    if(i < k){
      birds_in_fold <- c(sample(birds_F, N_birds_F), # With following
                         sample(birds_N, N_birds_N)) # Without following
      folds[[i]] <- which(data$Bird %in% birds_in_fold)
      birds_F <- birds_F[!(birds_F %in% birds_in_fold)]
      birds_N <- birds_N[!(birds_N %in% birds_in_fold)]
      
      # Else if we're on the last fold
    } else {
      birds_in_fold <- c(birds_F, birds_N)
      folds[[i]] <- which(data$Bird %in% birds_in_fold)
    }
  }
  
  # Return list of folds
  return(folds)
}








###################### 4. GET BOUTS AND FOLDS ###############################

# Create dataframe of bouts for whole dataset (corrected immersion)
imm_metrics_bouts1_list <- pblapply(imm_metrics_list1, get.bouts, window = 12, corrected = TRUE)
imm_metrics_bouts1 <- do.call(rbind, imm_metrics_bouts1_list)
saveRDS(imm_metrics_bouts1, "Output Dataset Files/imm_metrics_bouts1.RData")

# Create dataframe of bouts for whole dataset (uncorrected immersion)
imm_metrics_bouts1_gls_list <- pblapply(imm_metrics_list1, get.bouts, window = 12, corrected = FALSE)
imm_metrics_bouts1_gls <- do.call(rbind, imm_metrics_bouts1_gls_list)
saveRDS(imm_metrics_bouts1_gls, "Output Dataset Files/imm_metrics_bouts1_gls.RData")

# # Inspect the bouts
# nrow(imm_metrics_bouts1) # 1863
# sum(imm_metrics_bouts1$follow == "F") # 19
# imm_metrics_bouts1 %>% filter(follow == "F")
# imm_metrics_bouts1 %>% filter(follow == "F") %>% summarise(dur = mean(duration))
# imm_metrics_bouts1 %>% filter(follow == "N") %>% summarise(dur = mean(duration))


# Create folds for k-fold cross-validation
# k = 2 for demonstration; k = 5 for full dataset
folds <- custom_fold(imm_metrics_bouts1, k = 2)
folds_gls <- custom_fold(imm_metrics_bouts1_gls, k = 2)










############################### 5. RANDOM FOREST 5: Regularity only ####################################

# Set FP threshold
fp_threshold <- 0.1

# Create trainControl object
# 2-fold CV for demo; 5-fold CV for full dataset
# Custom summary function to assess model performance
train_control5 <- trainControl(
  method = "cv",
  number = 2,
  summaryFunction = custom_summary,
  classProbs = TRUE,  # To calculate probabilities
  savePredictions = "final",
  selectionFunction = "best",
  indexOut = folds # Stratify by individual birds
)

# Define parameters
mtry <- sqrt(1) # Only 1 predictor here
tunegrid <- expand.grid(.mtry = mtry)

# Train the random forest model
rf5 <- train(
  follow ~ P_LWLD,
  data = imm_metrics_bouts1,
  method = "rf",
  trControl = train_control5,
  tuneGrid = tunegrid,
  metric = "sensitivity_adj",
  maximize = TRUE
)

# Save model
saveRDS(rf5, file = "Output Dataset Files/rf5.RData")

# Results
rf5
summaries5 <- lapply(unique(rf5$pred$Resample), function(resample){
  custom_summary(rf5$pred[rf5$pred$Resample == resample,])
})
summaries5_df <- as.data.frame(do.call(rbind, summaries5))
summaries5_df











###################### 6. RANDOM FOREST 6: All predictors ###############################

# Create trainControl object
# 2-fold CV for demo; 5-fold CV for full dataset
# Custom summary function to assess model performance
train_control6 <- trainControl(
  method = "cv",
  number = 2,
  summaryFunction = custom_summary,
  classProbs = TRUE,  # To calculate probabilities
  savePredictions = "final",
  selectionFunction = "best",
  indexOut = folds
)

# Define parameters
mtry <- sqrt(5)
tunegrid <- expand.grid(.mtry = mtry)

# Train the random forest model
rf6 <- train(
  follow ~ P_LWLD + P_wet + mean_wet + mean_dry + N_landings,
  data = imm_metrics_bouts1,
  method = "rf",
  trControl = train_control6,
  tuneGrid = tunegrid,
  metric = "sensitivity_adj",
  maximize = TRUE
)

# Save model
saveRDS(rf6, file = "Output Dataset Files/rf6.RData")

# Results
rf6
summaries6 <- lapply(unique(rf6$pred$Resample), function(resample){
  custom_summary(rf6$pred[rf6$pred$Resample == resample,])
})
summaries6_df <- as.data.frame(do.call(rbind, summaries6))
summaries6_df









############ 7. RANDOM FOREST 7: Regularity only, uncorrected GLS ###############################

# Set FP threshold
fp_threshold <- 0.1

# Create trainControl object
# 2-fold CV for demo; 5-fold CV for full dataset
# Custom summary function to assess model performance
train_control7 <- trainControl(
  method = "cv",
  number = 2,
  summaryFunction = custom_summary,
  classProbs = TRUE,  # To calculate probabilities
  savePredictions = "final",
  selectionFunction = "best",
  indexOut = folds_gls # Stratify by individual birds
)

# Define parameters
mtry <- sqrt(1) # Only 1 predictor here
tunegrid <- expand.grid(.mtry = mtry)

# Train the random forest model
rf7 <- train(
  follow ~ P_LWLD,
  data = imm_metrics_bouts1_gls,
  method = "rf",
  trControl = train_control7,
  tuneGrid = tunegrid,
  metric = "sensitivity_adj",
  maximize = TRUE
)

# Save model
saveRDS(rf7, file = "Output Dataset Files/rf7.RData")

# Results
rf7
summaries7 <- lapply(unique(rf7$pred$Resample), function(resample){
  custom_summary(rf7$pred[rf7$pred$Resample == resample,])
})
summaries7_df <- as.data.frame(do.call(rbind, summaries7))
summaries7_df









###################### 8. RANDOM FOREST 8: All predictors, uncorrected GLS ######################

# Create trainControl object
# 2-fold CV for demo; 5-fold CV for full dataset
# Custom summary function to assess model performance
train_control8 <- trainControl(
  method = "cv",
  number = 2,
  summaryFunction = custom_summary,
  classProbs = TRUE,  # To calculate probabilities
  savePredictions = "final",
  selectionFunction = "best",
  indexOut = folds_gls
)

# Define parameters
mtry <- sqrt(5)
tunegrid <- expand.grid(.mtry = mtry)

# Train the random forest model
rf8 <- train(
  follow ~ P_LWLD + P_wet + mean_wet + mean_dry + N_landings,
  data = imm_metrics_bouts1_gls,
  method = "rf",
  trControl = train_control8,
  tuneGrid = tunegrid,
  metric = "sensitivity_adj",
  maximize = TRUE
)

# Save model
saveRDS(rf8, file = "Output Dataset Files/rf8.RData")

# Results
rf8
summaries8 <- lapply(unique(rf8$pred$Resample), function(resample){
  custom_summary(rf8$pred[rf8$pred$Resample == resample,])
})
summaries8_df <- as.data.frame(do.call(rbind, summaries4))
summaries8_df



