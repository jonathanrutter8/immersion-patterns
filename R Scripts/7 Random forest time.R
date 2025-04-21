#######################################################################################
################### BBA RANDOM FOREST: TIME-BASED #####################################
#######################################################################################

# Fits random forest models to predict vessel following from immersion metrics on a timestep-by-timestep basis. 
# Implements custom k-fold cross-validation to assess model performance
# (for full dataset, k = 5; for demonstration dataset, k = 2). 
# Demonstration dataset will not produce same results.


################### CONTENTS ####################

# 1. Load packages and data
# 2. Custom folding function for cross-validation
# 3. Custom summary function
# 4. Random forest 1: Regularity only
# 5. Random forest 2: All predictors
# 6. Random forest 3: Regularity only, uncorrected GLS
# 7. Random forest 4: Regularity only, uncorrected GLS
# 8. All summaries




################## 1. LOAD PACKAGES AND DATA ##################################

# Packages
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

# Combine into single dataframe
imm_metrics_df1 <- do.call(rbind, imm_metrics_list1)

# Add bout labels (optional)
imm_metrics_df1 <- imm_metrics_df1 %>% mutate(start_bout = as.numeric(NA))
imm_metrics_df1$start_bout <- ifelse(imm_metrics_df1$imm_factor == "intermediate" & 
                                       lag(imm_metrics_df1$imm_factor != "intermediate"),
                                     1, 0)
if(imm_metrics_df1$imm_factor[1] == "intermediate"){
  imm_metrics_df1$start_bout[1] <- 1
}
imm_metrics_df1 <- imm_metrics_df1 %>%
  filter(imm_factor == "intermediate") %>%
  mutate(bout = cumsum(start_bout)) %>%
  dplyr::select(-start_bout)

# Filter to only intermediate immersion timestamps
imm_metrics_df2 <- imm_metrics_df1 %>%
  filter(imm_factor == "intermediate") %>%
  mutate(follow = ifelse(!is.na(follow), "F", "N")) %>%
  mutate(follow = relevel(as.factor(follow), ref = "N"))


# Get dataframe with just birds with GLS
imm_metrics_list1_gls <- readRDS("Output Dataset Files/imm_metrics_list1_gls.RData")
imm_metrics_gls1 <- do.call(rbind, imm_metrics_list1_gls) %>%
  filter(imm_factor == "intermediate") %>%
  mutate(follow = ifelse(!is.na(follow), "F", "N")) %>%
  mutate(follow = relevel(as.factor(follow), ref = "N"))







############### 2. CUSTOM FOLDING FUNCTION FOR CROSS-VALIDATION ###################
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

# Optional: Custom summary function for stratifying by bouts
custom_fold_bouts <- function(data, k){
  
  # Get vectors of bouts with (F) and without (N) vessel following
  bouts <- unique(data$bout)
  bouts_F <- unique(data$bout[data$follow == "F"])
  bouts_N <- bouts[!(bouts %in% bouts_F)]
  
  # Calculate how many bouts in each fold
  N_bouts_F <- round(length(bouts_F)/k)
  N_bouts_N <- round(length(bouts_N)/k)
  
  # Initialise list to contain folds, and empty vector of bouts already used
  folds <- vector("list", k)
  
  # Randomly allocate bouts to folds
  # Ensure equal number of vessel-following bouts in each fold
  for(i in 1:k){
    
    # If we're not on the last fold
    if(i < k){
      bouts_in_fold <- c(sample(bouts_F, N_bouts_F), # With following
                         sample(bouts_N, N_bouts_N)) # Without following
      folds[[i]] <- which(data$bout %in% bouts_in_fold)
      bouts_F <- bouts_F[!(bouts_F %in% bouts_in_fold)]
      bouts_N <- bouts_N[!(bouts_N %in% bouts_in_fold)]
      
      # Else if we're on the last fold
    } else {
      bouts_in_fold <- c(bouts_F, bouts_N)
      folds[[i]] <- which(data$bout %in% bouts_in_fold)
    }
  }
  
  # Return list of folds
  return(folds)
}

# Create folds for k-fold cross-validation (P_LWLD only)
# k = 2 for demonstration; k = 5 for full dataset
folds <- custom_fold(imm_metrics_df2, k = 2)
folds_bouts <- custom_fold_bouts(imm_metrics_df2, k = 2)
folds_gls <- custom_fold(imm_metrics_gls1, k = 2)
folds_bouts_gls <- custom_fold_bouts(imm_metrics_gls1, k = 2)

# Create folds for k-fold CV (all predictors - need to narrow dataset to rows without NA)
# k = 2 for demonstration; k = 5 for full dataset
imm_metrics_df2.1 <- imm_metrics_df2[!is.na(imm_metrics_df2$mean_wet),]
folds2 <- custom_fold(imm_metrics_df2.1, k = 2)
folds2_bouts <- custom_fold_bouts(imm_metrics_df2.1, k = 2)
imm_metrics_gls1.1 <- imm_metrics_gls1[!is.na(imm_metrics_gls1$mean_wet),]
folds2_gls <- custom_fold(imm_metrics_gls1.1, k = 2)
folds2_bouts_gls <- custom_fold_bouts(imm_metrics_gls1.1, k = 2)






############## 3. CUSTOM SUMMARY FUNCTIONS ##########################

# Custom summary function
# Goal: Maximise sensitivity while keeping false pos:obs pos ratio < 0.1
custom_summary <- function(data, lev = NULL, model = NULL){
  
  # Extract observed and predicted probabilities
  obs <- data$obs
  pred <- data$pred
  #prob <- data$F
  
  # Get vectors of false/true pos/neg and observed pos/neg
  false_pos <- sum(ifelse(pred == "F" & obs == "N", 1, 0))
  false_neg <- sum(ifelse(pred == "N" & obs == "F", 1, 0))
  true_pos <- sum(ifelse(pred == "F" & obs == "F", 1, 0))
  true_neg <- sum(ifelse(pred == "N" & obs == "N", 1, 0))
  obs_pos <- sum(ifelse(obs == "F", 1, 0))
  obs_neg <- sum(ifelse(obs == "N", 1, 0))
  
  # Calculate sensitivity and specificity for this threshold
  sensitivity <- true_pos/obs_pos
  specificity <- true_neg/obs_neg
  
  # Calculate false pos percent (want this to be below fp_threshold)
  false_pos_percent <- false_pos/obs_pos
  
  # Prevent dividing by 0 (should not happen with stratified folds)
  if(obs_pos == 0){
    return(c(sensitivity = NA, 
             specificity = NA,
             false_pos = NA,
             obs_pos = NA,
             false_pos_percent = NA,
             sensitivity_adj = 0))
  } else {
    
    # If false positives sufficiently low
    if(false_pos_percent < fp_threshold){ # This is the FP threshold, from outside function
      message("This model met the FP threshold.")
      return(c(sensitivity = sensitivity, 
               specificity = specificity,
               false_pos = false_pos,
               obs_pos = obs_pos,
               false_pos_percent = false_pos_percent,
               sensitivity_adj = sensitivity))
      
    # Else if FP threshold
    } else {
      message("This model did NOT meet the FP threshold.")
      return(c(sensitivity = sensitivity, 
               specificity = specificity,
               false_pos = false_pos,
               obs_pos = obs_pos,
               false_pos_percent = false_pos_percent,
               sensitivity_adj = 0)) # If false pos percent too high, set sensitivity to 0 so it's not maxed
    }
  }
}









############ 4. RANDOM FOREST 1: Regularity only ###############################

# Set FP threshold
fp_threshold <- 0.1

# Create trainControl object
# 2-fold CV for demo; 5-fold CV for full dataset
# Custom summary function to assess model performance
train_control1 <- trainControl(
  method = "cv",
  number = 2, # Change to 5 for full dataset
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
rf1 <- train(
  follow ~ P_LWLD,
  data = imm_metrics_df2,
  method = "rf",
  trControl = train_control1,
  tuneGrid = tunegrid,
  metric = "sensitivity_adj",
  maximize = TRUE
)

# Save model
saveRDS(rf1, file = "Output Dataset Files/rf1.RData")

# Results
rf1
summaries1 <- lapply(unique(rf1$pred$Resample), function(resample){
  custom_summary(rf1$pred[rf1$pred$Resample == resample,])
})
summaries1_df <- as.data.frame(do.call(rbind, summaries1))
summaries1_df









###################### 5. RANDOM FOREST 2: All predictors ###############################

# Create trainControl object
# 2-fold CV for demo; 5-fold CV for full dataset
# Custom summary function to assess model performance
train_control2 <- trainControl(
  method = "cv",
  number = 2, # Change to 5 for full dataset
  summaryFunction = custom_summary,
  classProbs = TRUE,  # To calculate probabilities
  savePredictions = "final",
  selectionFunction = "best",
  indexOut = folds2
)

# Define parameters
mtry <- sqrt(5)
tunegrid <- expand.grid(.mtry = mtry)

# Train the random forest model
rf2 <- train(
  follow ~ P_LWLD + P_wet + mean_wet + mean_dry + N_landings,
  data = imm_metrics_df2.1,
  method = "rf",
  trControl = train_control2,
  tuneGrid = tunegrid,
  metric = "sensitivity_adj",
  maximize = TRUE
)

# Save model
saveRDS(rf2, file = "Output Dataset Files/rf2.RData")

# Results
rf2
summaries2 <- lapply(unique(rf2$pred$Resample), function(resample){
  custom_summary(rf2$pred[rf2$pred$Resample == resample,])
})
summaries2_df <- as.data.frame(do.call(rbind, summaries2))
summaries2_df











############ 6. RANDOM FOREST 3: Regularity only, uncorrected GLS ###############################

# Set FP threshold
fp_threshold <- 0.1

# Create trainControl object
# 2-fold CV for demo; 5-fold CV for full dataset
# Custom summary function to assess model performance
train_control3 <- trainControl(
  method = "cv",
  number = 2, # Change to 5 for full dataset
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
rf3 <- train(
  follow ~ P_LWLD_gls,
  data = imm_metrics_gls1,
  method = "rf",
  trControl = train_control3,
  tuneGrid = tunegrid,
  metric = "sensitivity_adj",
  maximize = TRUE
)

# Save model
saveRDS(rf3, file = "Output Dataset Files/rf3.RData")

# Results
rf3
summaries3 <- lapply(unique(rf3$pred$Resample), function(resample){
  custom_summary(rf3$pred[rf3$pred$Resample == resample,])
})
summaries3_df <- as.data.frame(do.call(rbind, summaries3))
summaries3_df









###################### 7. RANDOM FOREST 4: All predictors, uncorrected GLS ######################

# Create trainControl object
# 2-fold CV for demo; 5-fold CV for full dataset
# Custom summary function to assess model performance
train_control4 <- trainControl(
  method = "cv",
  number = 2, # Change to 5 for full dataset
  summaryFunction = custom_summary,
  classProbs = TRUE,  # To calculate probabilities
  savePredictions = "final",
  selectionFunction = "best",
  indexOut = folds2_gls
)

# Define parameters
mtry <- sqrt(5)
tunegrid <- expand.grid(.mtry = mtry)

# Train the random forest model
rf4 <- train(
  follow ~ P_LWLD_gls + P_wet_gls + mean_wet_gls + mean_dry_gls + N_landings_gls,
  data = imm_metrics_gls1.1,
  method = "rf",
  trControl = train_control4,
  tuneGrid = tunegrid,
  metric = "sensitivity_adj",
  maximize = TRUE
)

# Save model
saveRDS(rf4, file = "Output Dataset Files/rf4.RData")

# Results
rf4
summaries4 <- lapply(unique(rf4$pred$Resample), function(resample){
  custom_summary(rf4$pred[rf4$pred$Resample == resample,])
})
summaries4_df <- as.data.frame(do.call(rbind, summaries4))
summaries4_df










################ 8. All summaries #################

# Package to arrange plots
library(ggpubr)

# Read in data if needed
rf1 <- readRDS("Output Dataset Files/rf1.RData")
rf2 <- readRDS("Output Dataset Files/rf2.RData")
rf3 <- readRDS("Output Dataset Files/rf3.RData")
rf4 <- readRDS("Output Dataset Files/rf4.RData")

# Summarised model outputs
rf1
rf2
rf3
rf4

# Variable importance
varImp2 <- varImp(rf2)$importance
varImp4 <- varImp(rf4)$importance
vip2 <- ggplot(varImp2, aes(x = Overall, 
                            y = reorder(rownames(varImp2), Overall))) + 
  geom_col() +
  scale_y_discrete(labels = c("Number of landings", "Proportion wet", "Mean dry", "Mean wet", "Regularity")) +
  labs(x = "Importance (scaled mean decrease in accuracy)", y = "Predictor") +
  theme(panel.grid.major.y = element_blank())
vip4 <- ggplot(varImp4, aes(x = Overall, y = reorder(gsub("_gls", "", rownames(varImp4)), Overall))) + 
  geom_col() +
  scale_y_discrete(labels = c("Number of landings", "Proportion wet", "Mean dry", "Mean wet", "Regularity")) +
  labs(x = "Importance (scaled mean decrease in accuracy)", y = "Predictor") +
  theme(panel.grid.major.y = element_blank())

vip2
vip4






