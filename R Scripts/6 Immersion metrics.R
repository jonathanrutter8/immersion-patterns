###############################################################
################# Immersion Metrics Analysis ##################
###############################################################

# Calculates several rolling metrics, including immersion regularity, from GLS-immersion data.
# Tests the predictive performance of the regularity metric across parameter space using a binary threshold classifier model.

################ CONTENTS ##############################

# 1. Packages, WD, data
# 2. Function: Calculate regularity and other metrics in sliding window
# 3. Function: Calculate regularity and other metrics for uncorrected immersion
# 4. Functions to test regularity across parameter space
# 5. Overarching function to run it all together
# 6. Run the functions in parallel: wide parameter space
# 7. Run the functions in parallel: narrow parameter space
# 8. Test out variations on window2 and peak expansion
# 9. Produce final immersion dataframe with corrected immersion
# 10. Produce final immersion dataframe with uncorrected immersion




############# 1. PACKAGES, WD, DATA #############################

library(tidyverse)
library(lubridate)
library(pbapply)
library(parallel)
library(caret)
theme_set(theme_bw())

# Set seed
set.seed(88)

# Set WD
setwd("___________")

# Read in corrected regular immersion data
imm_list <- readRDS("Output Dataset Files/Wet_Data_Reg_Corr3.RData")
imm_list <- lapply(imm_list, function(x){
  mutate(x, leg_wet_C = as.character(leg_wet_C))
})

# Narrow down to birds with GLS or HMM immersion data to analyse
imm_list <- imm_list[!is.na(sapply(imm_list, function(x){x$leg_wet_C[1]}))]

# Read in visual interactions analysis data
int_summary <- read.csv("Visual Interactions Analysis/BBA Interactions Summary.csv") %>%
  mutate(Real_start = dmy_hm(Real_start), Real_end = dmy_hm(Real_end))
follows <- int_summary[int_summary$Int_Class == "Follow", c(1,2,4,6,7)]

# # If necessary, read in imm_metrics1 list (created from P_LWLD function below)
# imm_metrics1 <- readRDS("Output Dataset Files/imm_metrics.Rdata")


# May need to correct interaction columns
imm_list <- pblapply(imm_list, function(x){

  # Which bird is it
  bird <- x$Bird[1]

  # Clear the follow and stop vectors
  x$follow <- NA
  x$stop <- NA

  # Get interaction times for this bird
  int_bird <- int_summary[int_summary$Bird == bird, 1:7]

  # If the bird had interaction
  if(nrow(int_bird) > 0){

    for(i in 1:nrow(int_bird)){

      # If bird is following vessel, add Int_Num to track and immersion dataframes
      if(int_bird$Int_Class[i] == "Follow"){
        x_Fvec <- ifelse(x$DateTime >= int_bird$Real_start[i] &
                           x$DateTime <= int_bird$Real_end[i],
                         int_bird$Int_Num[i], NA)
        x$follow <- ifelse(is.na(x$follow), x_Fvec, x$follow)

        # Else if bird stops near vessel, add Int_Num to track and immersion dataframes
      } else if(int_bird$Int_Class[i] == "Stop"){
        x_Svec <- ifelse(x$DateTime >= int_bird$Real_start[i] &
                           x$DateTime <= int_bird$Real_end[i],
                         int_bird$Int_Num[i], NA)
        x$stop <- ifelse(is.na(x$stop), x_Svec, x$stop)
      }
    }
  }

  return(x)

})





############### 2. FUNCTION: CALCULATE REGULARITY AND OTHER METRICS IN SLIDING WINDOW #############################

# Function to calculate regularity (P_LWLD, proportion low wet low dry) and other metrics in sliding window
# window is to calculate most metrics. Should be in minutes, and should be cleanly divisible by imm_df interval
# SDthreshold is what defines a "low" SD of wet or dry duration. Should be in minutes.
# window2 is the second smoothing window for P_LWLD. Same criteria, but does not need to be the same number.
calculate.regularity <- function(imm_df, window, SDthreshold, window2){
  
  # Set buffer of size window/2
  # buffer unit = rows, window unit = mins, dt = secs
  dt <- as.numeric(difftime(imm_df$DateTime[2], imm_df$DateTime[1], units = "secs"))
  buffer <- floor(window*60/(dt*2))
  
  # Classify periods of all dry or all wet so we don't waste time
  # putting a rolling window through those periods.
  
  # Identify wet and dry periods
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
  
  # Put columns for several metrics into dataframe
  durations <- imm_df %>%
    dplyr::select(!c(transition, period)) %>%
    mutate(P_wet = ifelse(imm_factor == "all dry", 1,
                          ifelse(imm_factor == "all wet", 0, NA)),
           mean_wet = ifelse(imm_factor == "all dry", 0,
                             ifelse(imm_factor == "all wet", window*60, NA)),
           mean_dry = ifelse(imm_factor == "all dry", window*60,
                             ifelse(imm_factor == "all wet", 0, NA)),
           N_landings = 0,
           sd_wet = NA,
           sd_dry = NA,
           P_LWLD = 0)
  
  # Get vector of "intermediate" rows to iterate through
  inter_rows0 <- which(durations$imm_factor == "intermediate")
  inter_rows <- inter_rows0[inter_rows0 >= buffer & inter_rows0 < (max(inter_rows0) - buffer)]
  
  # Iterate through time series with rolling window, skipping known all dry/all wet periods
  for(i in inter_rows){
    
    # Find lower and upper datetimes of the current window
    lower <- durations$DateTime[i - buffer + 1] # - minutes(window/2) + seconds(dt)
    upper <- durations$DateTime[i + buffer] # + minutes(window/2)
    
    # For now, add one extra row on either side to help calculate landings and takeoffs accurately
    if(i != inter_rows[1]){
      start_row <- which(imm_df$DateTime == lower) - 1
    } else {
      start_row <- which(imm_df$DateTime == lower)
    }
    if(i != inter_rows[length(inter_rows)]){
      end_row <- which(imm_df$DateTime == upper) + 1
    } else {
      end_row <- which(imm_df$DateTime == upper)
    }
    
    # Get portion of buffered time series within the window
    df_subset <- imm_df[start_row:end_row, ] # Using imm_df here because fewer columns, same rows
    
    # Define landings and takeoffs (start and end of wet periods)
    df_subset <- mutate(df_subset, 
                        landing = ifelse(df_subset$leg_wet_C == "wet" & 
                                           (lag(df_subset$leg_wet_C) == "dry"), 
                                         1, 0), 
                        takeoff = ifelse(df_subset$leg_wet_C == "wet" & 
                                           (lead(df_subset$leg_wet_C) == "dry"), 
                                         1, 0))
    
    # Now that landings/takeoffs done, take out first and last row
    start_row <- which(df_subset$DateTime == lower)
    end_row <- which(df_subset$DateTime == upper)
    df_subset <- df_subset[start_row:end_row, ]
    
    # If any NAs left for takeoffs and landings, switch to 0
    df_subset$takeoff <- ifelse(is.na(df_subset$takeoff), 0, df_subset$takeoff)
    df_subset$landing <- ifelse(is.na(df_subset$landing), 0, df_subset$landing)
    
    # Calculate proportion of wet (0-1)
    P_wet <- sum(df_subset$leg_wet_C == "wet")/nrow(df_subset)
    durations$P_wet[i] <- P_wet
    
    # Calculate number of landings
    durations$N_landings[i] <- sum(df_subset$landing)
    
    # Initiate vectors for durations of wet and dry
    wet_vec <- c()
    dry_vec <- c()
    dur_wet <- 0
    dur_dry <- 0
    
    # Calculate durations (in seconds) and add to vectors
    for(j in 1:nrow(df_subset)){
      
      # If end of wet period
      if (df_subset$takeoff[j] == 1){
        dur_wet <- dur_wet + dt
        wet_vec <- c(wet_vec, dur_wet)
        dur_wet <- 0
        
      # If start or middle of wet period  
      } else if (df_subset$leg_wet_C[j] == "wet"){
        dur_wet <- dur_wet + dt
        
      # If end of dry period
      } else if (df_subset$leg_wet_C[j] == "dry" & df_subset$landing[j+1] %in% c(1, NA)){
        dur_dry <- dur_dry + dt
        dry_vec <- c(dry_vec, dur_dry)
        dur_dry <- 0
        
      # If start or middle of dry period
      } else if (df_subset$leg_wet_C[j] == "dry"){
        dur_dry <- dur_dry + dt
      }
    }
    
    # Calculate mean and SD of wet and dry durations, add to durations df
    durations$mean_wet[i] <- ifelse(!is.null(wet_vec), mean(wet_vec), 0)
    durations$sd_wet[i] <- ifelse(!is.null(wet_vec), sd(wet_vec), NA)
    durations$mean_dry[i] <- ifelse(!is.null(dry_vec), mean(dry_vec), 0)
    durations$sd_dry[i] <- ifelse(!is.null(dry_vec), sd(dry_vec), NA)
  }
  
  # Calculate proportion of LWLD (low wet low dry SD) in another sliding window
  for(i in inter_rows){
    startDT <- durations$DateTime[i] - minutes(window2/2) + seconds(dt)
    start_i <- min(which(abs(durations$DateTime - startDT) == 
                           min(abs(durations$DateTime - startDT))))
    endDT <- durations$DateTime[i] + minutes(window2/2)
    end_i <- min(which(abs(durations$DateTime - endDT) == 
                         min(abs(durations$DateTime - endDT))))
    LWLD <- durations[start_i:end_i, ] %>%
      filter(sd_wet < (SDthreshold*60) & sd_dry < (SDthreshold*60)) # Convert SD Threshold to seconds
    durations$P_LWLD[i] <- nrow(LWLD)/nrow(durations[start_i:end_i, ])
  }
  
  # Return durations df
  #message(paste0(durations$Bird[1], " - finished calculate.regularity."))
  return(durations)
}








############### 3. FUNCTION: CALCULATE REGULARITY AND OTHER METRICS FOR UNCORRECTED IMMERSION #############################

# Function to calculate P_LWLD_gls (regularity for uncorrected immersion) and other metrics in sliding window
# Key differences from above: Uses leg_wet data not leg_wet_C, mutates differently named columns
# window is to calculate most metrics. Should be in minutes, and should be cleanly divisible by imm_df interval
# SDthreshold is what defines a "low" SD of wet or dry duration. Should be in minutes.
# window2 is the second smoothing window for P_LWLD. Same criteria, but does not need to be the same number.
calculate.regularity.gls <- function(imm_df, window, SDthreshold, window2){
  
  # Set leg_wet_C to leg_wet value, just within this function
  if(is.na(imm_df$leg_wet[1])){
    message(paste0(imm_df$Bird[1], " did not have a GLS, returned NA."))
    return(NA)
  } else {
    imm_df$leg_wet_C <- imm_df$leg_wet
  }
  
  # Set buffer of size window/2
  # buffer unit = rows, window unit = mins, dt = secs
  dt <- as.numeric(difftime(imm_df$DateTime[2], imm_df$DateTime[1], units = "secs"))
  buffer <- floor(window*60/(dt*2))
  
  # Classify periods of all dry or all wet so we don't waste time
  # putting a rolling window through those periods.
  
  # Identify wet and dry periods
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
  
  # Put columns for several metrics into dataframe
  durations <- imm_df %>%
    mutate(P_wet_gls = ifelse(imm_factor == "all dry", 1,
                          ifelse(imm_factor == "all wet", 0, NA)),
           mean_wet_gls = ifelse(imm_factor == "all dry", 0,
                             ifelse(imm_factor == "all wet", window*60, NA)),
           mean_dry_gls = ifelse(imm_factor == "all dry", window*60,
                             ifelse(imm_factor == "all wet", 0, NA)),
           N_landings_gls = 0,
           sd_wet_gls = NA,
           sd_dry_gls = NA,
           P_LWLD_gls = 0)
  
  # Get vector of "intermediate" rows to iterate through
  inter_rows0 <- which(durations$imm_factor == "intermediate")
  inter_rows <- inter_rows0[inter_rows0 >= buffer & inter_rows0 < (max(inter_rows0) - buffer)]
  
  # Iterate through time series with rolling window, skipping known all dry/all wet periods
  for(i in inter_rows){
    
    # Find lower and upper datetimes of the current window
    lower <- durations$DateTime[i - buffer + 1] # - minutes(window/2) + seconds(dt)
    upper <- durations$DateTime[i + buffer] # + minutes(window/2)
    
    # For now, add one extra row on either side to help calculate landings and takeoffs accurately
    if(i != inter_rows[1]){
      start_row <- which(imm_df$DateTime == lower) - 1
    } else {
      start_row <- which(imm_df$DateTime == lower)
    }
    if(i != inter_rows[length(inter_rows)]){
      end_row <- which(imm_df$DateTime == upper) + 1
    } else {
      end_row <- which(imm_df$DateTime == upper)
    }
    
    # Get portion of buffered time series within the window
    df_subset <- imm_df[start_row:end_row, ] # Using imm_df here because fewer columns, same rows
    
    # Define landings and takeoffs (start and end of wet periods)
    df_subset <- mutate(df_subset, 
                        landing = ifelse(df_subset$leg_wet_C == "wet" & 
                                           (lag(df_subset$leg_wet_C) == "dry"), 
                                         1, 0), 
                        takeoff = ifelse(df_subset$leg_wet_C == "wet" & 
                                           (lead(df_subset$leg_wet_C) == "dry"), 
                                         1, 0))
    
    # Now that landings/takeoffs done, take out first and last row
    start_row <- which(df_subset$DateTime == lower)
    end_row <- which(df_subset$DateTime == upper)
    df_subset <- df_subset[start_row:end_row, ]
    
    # If any NAs left for takeoffs and landings, switch to 0
    df_subset$takeoff <- ifelse(is.na(df_subset$takeoff), 0, df_subset$takeoff)
    df_subset$landing <- ifelse(is.na(df_subset$landing), 0, df_subset$landing)
    
    # Calculate proportion of wet (0-1)
    P_wet <- sum(df_subset$leg_wet_C == "wet")/nrow(df_subset)
    durations$P_wet_gls[i] <- P_wet
    
    # Calculate number of landings
    durations$N_landings_gls[i] <- sum(df_subset$landing)
    
    # Initiate vectors for durations of wet and dry
    wet_vec <- c()
    dry_vec <- c()
    dur_wet <- 0
    dur_dry <- 0
    
    # Calculate durations (in seconds) and add to vectors
    for(j in 1:nrow(df_subset)){
      
      # If end of wet period
      if (df_subset$takeoff[j] == 1){
        dur_wet <- dur_wet + dt
        wet_vec <- c(wet_vec, dur_wet)
        dur_wet <- 0
        
        # If start or middle of wet period  
      } else if (df_subset$leg_wet_C[j] == "wet"){
        dur_wet <- dur_wet + dt
        
        # If end of dry period
      } else if (df_subset$leg_wet_C[j] == "dry" & df_subset$landing[j+1] %in% c(1, NA)){
        dur_dry <- dur_dry + dt
        dry_vec <- c(dry_vec, dur_dry)
        dur_dry <- 0
        
        # If start or middle of dry period
      } else if (df_subset$leg_wet_C[j] == "dry"){
        dur_dry <- dur_dry + dt
      }
    }
    
    # Calculate mean and SD of wet and dry durations, add to durations df
    durations$mean_wet_gls[i] <- ifelse(!is.null(wet_vec), mean(wet_vec), 0)
    durations$sd_wet_gls[i] <- ifelse(!is.null(wet_vec), sd(wet_vec), NA)
    durations$mean_dry_gls[i] <- ifelse(!is.null(dry_vec), mean(dry_vec), 0)
    durations$sd_dry_gls[i] <- ifelse(!is.null(dry_vec), sd(dry_vec), NA)
  }
  
  # Calculate proportion of LWLD (low wet low dry SD) in another sliding window
  for(i in inter_rows){
    startDT <- durations$DateTime[i] - minutes(window2/2) + seconds(dt)
    start_i <- min(which(abs(durations$DateTime - startDT) == 
                           min(abs(durations$DateTime - startDT))))
    endDT <- durations$DateTime[i] + minutes(window2/2)
    end_i <- min(which(abs(durations$DateTime - endDT) == 
                         min(abs(durations$DateTime - endDT))))
    LWLD <- durations[start_i:end_i, ] %>%
      filter(sd_wet_gls < (SDthreshold*60) & sd_dry_gls < (SDthreshold*60)) # Convert SD Threshold to seconds
    durations$P_LWLD_gls[i] <- nrow(LWLD)/nrow(durations[start_i:end_i, ])
  }
  
  # Return durations df
  message(paste0(durations$Bird[1], " - finished calculate.regularity.gls."))
  return(durations)
}









###################### 4. FUNCTIONS TO TEST REGULARITY ACROSS PARAMETER SPACE ###################


# Function to run a single threshold model to test all P_LWLD threshold from 0-1
# data: dataframe that includes follow and P_LWLD columns
# thresholds: Can be vector or single number 0-1
# false_pos_threshold: Max proportion of false pos to observed pos
# peak_expansion: Logical, should peak expansion of predicted values be done before metrics calculated?
run.threshold.model <- function(data, thresholds, 
                                false_pos_threshold, peak_expansion = FALSE){
  
  # If there is no best threshold, return NAs, else carry on
  if(is.na(thresholds[1])){
    return(data.frame(threshold = NA, 
                      sensitivity = NA, 
                      specificity = NA,
                      false_pos = NA,
                      obs_pos = NA,
                      false_pos_percent = NA,
                      met_criteria = NA))
  }
  
  # Extract observed and predicted probabilities
  obs <- data$follow
  P_LWLD <- data$P_LWLD
  
  # Optimize the threshold based on the following criteria:
  # (1) Num false positive bouts must be < 10% total observed positive bouts
  # (2) Maximise sensitivity
  # This will help account for the highly inflated specificity of these models
  SS_list <- lapply(thresholds, function(threshold){
    
    # Get vector predicted response variable
    pred <- ifelse(P_LWLD >= threshold, "F", "N")
    
    # Peak expansion:
    # Expand pred to include beginning and end of regularity "peaks" above threshold
    if(peak_expansion == TRUE){
      
      for(i in 2:(length(pred) - 1)){
        
        # If tip of spike is above threshold
        if(pred[i] == "F"){
          
          # If previous predicted was "N", expand backwards
          if(pred[i-1] == "N"){
            
            for(j in 1:i){
              
              # If regularity > 0 and still going down the peak
              if(P_LWLD[i-j] > 0 & P_LWLD[i-j] < P_LWLD[i-j+1]){
                
                # Switch to 1
                pred[i-j] <- "F"
                
              } else {break}
            }
          }
          
          # If next predicted was "N", expand forwards
          if(pred[i+1] == "N"){
            
            for(k in 1:(length(pred) - i)){
              
              # If regularity > 0 and still going down the peak
              if(P_LWLD[i+k] > 0 & P_LWLD[i+k] < P_LWLD[i+k-1]){
                
                # Switch to 1
                pred[i+k] <- "F"
                
              } else {break}
            }
          }
        }
      }
    }

    
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
    
    # Calculate false positive percent (want this to be above threshold)
    false_pos_percent = false_pos/obs_pos
    
    # Prevent dividing by 0
    if(obs_pos == 0){
      return(c(threshold = threshold, 
               sensitivity = NA, 
               specificity = NA,
               false_pos = NA,
               obs_pos = NA,
               false_pos_percent = NA,
               met_criteria = NA))
    } else {
      
      # If false positives sufficiently low, meets criteria
      if(false_pos_percent < false_pos_threshold){
        return(c(threshold = threshold, 
                 sensitivity = sensitivity, 
                 specificity = specificity,
                 false_pos = false_pos,
                 obs_pos = obs_pos,
                 false_pos_percent = false_pos_percent,
                 met_criteria = TRUE))
        
      # Else if too many false positives, does not meet criteria
      } else {
        return(c(threshold = threshold, 
                 sensitivity = sensitivity, 
                 specificity = specificity,
                 false_pos = false_pos,
                 obs_pos = obs_pos,
                 false_pos_percent = false_pos_percent,
                 met_criteria = FALSE))
      }
    }
  })
  
  # Combine into single df and return
  SS <- as.data.frame(do.call(rbind, SS_list))
  #message("Finished run.threshold.model.")
  return(SS)
}




# Function to implement the above threshold model using k-fold cross-validation
# Tests how well different variations of P_LWLD predict vessel following
# df: A dataframe from calculate.regularity function above
# k: Number of folds
# Returns a list: 
# [[1]] Training data model performance across all P_LWLD thresholds
# [[2]] Best threshold
# [[3]] Test data model performance for best threshold
test.regularity <- function(df, k, false_pos_threshold, peak_expansion = FALSE){
  
  # Establish vector of P_LWLD thresholds to examine
  thresholds <- seq(0, 1, by = 0.01)
  
  # Reclassify "follow" dataset so that follow = "F" and non-follow = "N"
  df <- df %>%
    mutate(follow = ifelse(!is.na(follow), "F", "N")) %>%
    mutate(follow = relevel(as.factor(follow), ref = "N"))
  
  # Create train and test indices for k-fold cross-validation
  # By default, ensures all folds have same number of follow = "F"
  folds <- createFolds(df$follow, k = k, list = TRUE, returnTrain = TRUE)
  
  # Manual k-fold cross-validation
  results <- lapply(folds, function(train_i){
    
    # Create train and test datasets for this fold
    train_data <- df[train_i, ]
    test_data <- df[-train_i, ]
    
    # Train and test model
    train_results <- run.threshold.model(train_data, 
                                         thresholds = thresholds,
                                         false_pos_threshold = false_pos_threshold,
                                         peak_expansion = peak_expansion)
    good_results <- train_results[train_results$met_criteria == 1, ]
    best_threshold <- ifelse(nrow(good_results > 0), 
                             good_results$threshold[which.max(good_results$sensitivity)], NA)
    test_results <- run.threshold.model(test_data, 
                                        thresholds = best_threshold,
                                        false_pos_threshold = false_pos_threshold,
                                        peak_expansion = peak_expansion)
    return(list(train_results, best_threshold, test_results))
  })
  
  # Training dataframe from aggregating metrics across all folds
  agg_train_results <- data.frame(
    threshold = thresholds,
    sensitivity = NA,
    specificity = NA,
    false_pos = NA,
    obs_pos = NA,
    false_pos_percent = NA,
    met_criteria = NA
  )
  
  # Fill in training dataset columns by getting means across folds
  for(i in 2:7){
    agg_train_list <- lapply(results, function(x){as.numeric(x[[1]][, i])})
    agg_train_df <- do.call(cbind, agg_train_list)
    means <- apply(agg_train_df, 1, mean)
    agg_train_results[ ,i] <- means
  }
  
  # Get average best threshold
  avg_best_threshold <- mean(sapply(results, function(x){as.numeric(x[[2]])}))
  
  # Test dataframe from aggregating metrics across all folds
  agg_test_results <- data.frame(
    threshold = avg_best_threshold,
    sensitivity = mean(sapply(results, function(x){as.numeric(x[[3]]$sensitivity)})),
    specificity = mean(sapply(results, function(x){as.numeric(x[[3]]$specificity)})),
    false_pos = mean(sapply(results, function(x){as.numeric(x[[3]]$false_pos)})),
    obs_pos = mean(sapply(results, function(x){as.numeric(x[[3]]$obs_pos)})),
    false_pos_percent = mean(sapply(results, function(x){as.numeric(x[[3]]$false_pos_percent)})),
    met_criteria = mean(sapply(results, function(x){as.numeric(x[[3]]$met_criteria)}))
  )
  
  # Return list of aggregated results
  #message("test.regularity done.")
  return(list(agg_train_results, avg_best_threshold, agg_test_results))
}





################## 5. OVERARCHING FUNCTION TO RUN IT ALL TOGETHER ########################################



# Overarching function to run it all together
# i: Parameter combo from results_summary
test.regularity.all <- function(i, imm_list, 
                                results_summary,
                                k, false_pos_threshold){
  
  # Get parameters for this iteration
  window <- results_summary$window[i]
  SDthreshold <- results_summary$SDthreshold[i]
  window2 <- results_summary$window2[i]
  peak_expansion <- results_summary$peak_expansion[i]
  
  # Run the main functions
  imm_metrics_list <- lapply(imm_list, calculate.regularity,
                             window = window, SDthreshold = SDthreshold, window2 = window2)
  imm_metrics_df <- do.call(rbind, imm_metrics_list)
  test_results_list <- test.regularity(imm_metrics_df,
                                       k = k, false_pos_threshold = false_pos_threshold, 
                                       peak_expansion = peak_expansion)
  
  # Return list of results for this parameter set
  message(paste0("Parameter combination ", i, " done."))
  return(test_results_list)
}






######################## 6. RUN THE FUNCTIONS IN PARALLEL: WIDE PARAMETER SPACE ####################################

# Create results dataframe
results_summary1 <- data.frame(
  param_combo = 1:48,
  window = c(rep(2.5, 6), 
              rep(5, 6),
              rep(7.5, 6),
              rep(10, 6), 
              rep(12.5, 6),
              rep(15, 6),
              rep(20, 6), 
              rep(25, 6)),
  SDthreshold = rep(c(0.5, 1, 1.5, 2, 2.5, 3), 8),
  window2 = rep(10, 48),
  peak_expansion = rep(FALSE, 48),
  sensitivity = NA,
  false_pos_percent = NA
)

# Prep for parallel processing
numCores <- detectCores() - 2
cl <- makeCluster(numCores)
clusterEvalQ(cl, {
  library(tidyverse)
  library(lubridate)
  library(pbapply)
  library(parallel)
  library(caret)
  theme_set(theme_bw())
})
clusterExport(cl = cl, c("imm_list", "int_summary", "results_summary1",
                         "calculate.regularity", "test.regularity", "run.threshold.model",
                         "test.regularity.all"))

# Run the overarching function in parallel
all_results_list1 <- pblapply(results_summary1$param_combo, # This is "i"
                              test.regularity.all, # This is the function
                              results_summary = results_summary1,
                              imm_list = imm_list, 
                              k = 2, # FOR FULL ANALYSIS, k = 5
                              false_pos_threshold = 0.1, # Extra arguments
                              cl = cl)
stopCluster(cl)
#saveRDS(all_results_list1, file = "Output Dataset Files/all_results_list1.RData")

# Summarise results
results_summary1 <- results_summary1 %>%
  mutate(threshold = NA, met_criteria = NA) %>%
  relocate(threshold, .before = sensitivity)
for(i in results_summary1$param_combo){
  results_summary1$threshold[i] <- all_results_list1[[i]][[2]]
  results_summary1$sensitivity[i] <- all_results_list1[[i]][[3]]$sensitivity[1]
  results_summary1$false_pos_percent[i] <- all_results_list1[[i]][[3]]$false_pos_percent[1]
  results_summary1$met_criteria[i] <- all_results_list1[[i]][[3]]$met_criteria[1]
}
results_summary1
#saveRDS(results_summary1, file = "Output Dataset Files/immreg_results_summary1.RData")

# Heat map of parameter space
ggplot(results_summary1, aes(x = window, y = SDthreshold, fill = sensitivity)) + 
  geom_tile() +
  labs(x = "1st rolling window (a)", 
       y = "SD threshold (z)",
       title = "a. Wide parameter space")

# Which parameter combination looks best?
best1 <- results_summary1[which.max(results_summary1$sensitivity),]
best1







######################## 7. RUN THE FUNCTIONS IN PARALLEL: NARROW PARAMETER SPACE ####################################

# Create results dataframe
results_summary2 <- data.frame(
  param_combo = 1:16,
  window = c(rep(11, 4), 
             rep(12, 4),
             rep(13, 4),
             rep(14, 4)),
  SDthreshold = rep(c(0.8, 1, 1.2, 1.4), 4),
  window2 = rep(10, 16),
  peak_expansion = rep(FALSE, 16),
  sensitivity = NA,
  false_pos_percent = NA
)

# Prep for parallel processing
numCores <- detectCores()
cl <- makeCluster(numCores)
clusterEvalQ(cl, {
  library(tidyverse)
  library(lubridate)
  library(pbapply)
  library(parallel)
  library(caret)
  theme_set(theme_bw())
})
clusterExport(cl = cl, c("imm_list", "int_summary", "results_summary2",
                         "calculate.regularity", "test.regularity", "run.threshold.model",
                         "test.regularity.all"))

# Run the overarching function in parallel
all_results_list2 <- pblapply(results_summary2$param_combo, # This is "i"
                              test.regularity.all, # This is the function
                              results_summary = results_summary2,
                              imm_list = imm_list, 
                              k = 2, # FOR FULL ANALYSIS, k = 5
                              false_pos_threshold = 0.1, # Extra arguments
                              cl = cl)
stopCluster(cl)
#saveRDS(all_results_list2, file = "Output Dataset Files/all_results_list2.RData")

# Summarise results
results_summary2 <- results_summary2 %>%
  mutate(threshold = NA, met_criteria = NA) %>%
  relocate(threshold, .before = sensitivity)
for(i in results_summary2$param_combo){
  results_summary2$threshold[i] <- all_results_list2[[i]][[2]]
  results_summary2$sensitivity[i] <- all_results_list2[[i]][[3]]$sensitivity[1]
  results_summary2$false_pos_percent[i] <- all_results_list2[[i]][[3]]$false_pos_percent[1]
  results_summary2$met_criteria[i] <- all_results_list2[[i]][[3]]$met_criteria[1]
}

# Paste in some of the rows from original results_summary
results_summary2_special <- results_summary2 %>%
  filter((window == 10 | window == 15), (SDthreshold >= 0.5 & SDthreshold <= 1.5))
results_summary2 <- rbind(results_summary2, results_summary2_special) %>%
  arrange(window, SDthreshold)
results_summary2
#saveRDS(results_summary2, file = "Output Dataset Files/immreg_results_summary2.RData")

# Heat map of parameter space
ggplot(results_summary2, aes(x = window, y = SDthreshold, fill = sensitivity)) + 
  geom_tile() +
  labs(x = "1st rolling window (a)", 
       y = "SD threshold (z)",
       title = "b. Narrow parameter space")


# Which parameter combination looks best?
best2 <- results_summary2[which.max(results_summary2$sensitivity),]
best2








########################## 8. TEST OUT VARIATIONS ON WINDOW2 AND PEAK EXPANSION ######################

# Create results dataframe
results_summary3 <- data.frame(
  param_combo = 1:8,
  window = rep(best2$window[1], 8),
  SDthreshold = rep(best2$SDthreshold[1], 8),
  window2 = rep(c(8, 10, 12, 14), 2),
  peak_expansion = c(rep(TRUE, 4), rep(FALSE, 4)),
  sensitivity = NA,
  false_pos_percent = NA
)

# Prep for parallel processing
numCores <- detectCores()
cl <- makeCluster(numCores)
clusterEvalQ(cl, {
  library(tidyverse)
  library(lubridate)
  library(pbapply)
  library(parallel)
  library(caret)
  theme_set(theme_bw())
})
clusterExport(cl = cl, c("imm_list", "int_summary", "results_summary3",
                         "calculate.regularity", "test.regularity", "run.threshold.model",
                         "test.regularity.all"))

# Run the overarching function in parallel
all_results_list3 <- pblapply(results_summary3$param_combo, # This is "i"
                              test.regularity.all, # This is the function
                              results_summary = results_summary3,
                              imm_list = imm_list, 
                              k = 2, # FOR FULL ANALYSIS, k = 5
                              false_pos_threshold = 0.1, # Extra arguments
                              cl = cl)
stopCluster(cl)
saveRDS(all_results_list3, file = "Output Dataset Files/all_results_list3.RData")

# Summarise results
results_summary3 <- results_summary3 %>%
  mutate(threshold = NA, met_criteria = NA) %>%
  relocate(threshold, .before = sensitivity)
for(i in results_summary3$param_combo){
  results_summary3$threshold[i] <- all_results_list3[[i]][[2]]
  results_summary3$sensitivity[i] <- all_results_list3[[i]][[3]]$sensitivity[1]
  results_summary3$false_pos_percent[i] <- all_results_list3[[i]][[3]]$false_pos_percent[1]
  results_summary3$met_criteria[i] <- all_results_list3[[i]][[3]]$met_criteria[1]
}
results_summary3
#saveRDS(results_summary3, file = "Output Dataset Files/immreg_results_summary3.RData")

# Heat map of parameter space
ggplot(results_summary3, aes(x = window2, y = peak_expansion, fill = sensitivity)) + 
  geom_tile() +
  labs(x = "2nd rolling window (b)", 
       y = "Peak Expansion",
       title = "c. Testing 2nd rolling window")
  


# Which parameter combination looks best?
results_summary3[which.max(results_summary3$sensitivity),]


# Which parameter combination looks best?
best3 <- results_summary3[which.max(results_summary3$sensitivity),] %>%
  mutate(param_combo = 1)
best3











############################ 9. PRODUCE FINAL IMMERSION DATAFRAME WITH CORRECTED IMMERSION #######################

# Prep for parallel processing
numCores <- detectCores()
cl <- makeCluster(numCores)
clusterEvalQ(cl, {
  library(tidyverse)
  library(lubridate)
  library(pbapply)
  library(parallel)
  library(caret)
  theme_set(theme_bw())
})
clusterExport(cl = cl, c("imm_list", "int_summary",
                         "calculate.regularity", "test.regularity", "run.threshold.model",
                         "test.regularity.all"))

# Run the overarching function in parallel
# Produce immersion metrics list
imm_metrics_list1 <- pblapply(imm_list, calculate.regularity, 
                              window = 12, SDthreshold = 1, window2 = 12, # These metrics came from full dataset
                              cl = cl)
stopCluster(cl)
saveRDS(imm_metrics_list1, file = "Output Dataset Files/imm_metrics_list1.RData")

# Bind into single dataframe
imm_metrics_df1 <- do.call(rbind, imm_metrics_list1)














############################ 10. PRODUCE FINAL IMMERSION DATAFRAME WITH UNCORRECTED IMMERSION #######################

# Narrow list to birds with GLS
w_gls <- which(sapply(imm_list, function(x){!is.na(x$leg_wet[1])}))
imm_list_gls <- imm_list[w_gls]

# Prep for parallel processing
numCores <- detectCores()
cl <- makeCluster(numCores)
clusterEvalQ(cl, {
  library(tidyverse)
  library(lubridate)
  library(pbapply)
  library(parallel)
  library(caret)
  theme_set(theme_bw())
})
clusterExport(cl = cl, c("imm_list_gls", "int_summary",
                         "calculate.regularity.gls"))

# Run the overarching function in parallel
# Produce immersion metrics list
imm_metrics_list1_gls <- pblapply(imm_list_gls, calculate.regularity.gls, 
                              window = 12, SDthreshold = 1, window2 = 12, # These metrics came from full dataset
                              cl = cl)
stopCluster(cl)
saveRDS(imm_metrics_list1_gls, file = "Output Dataset Files/imm_metrics_list1_gls.RData")

# Bind into single dataframe
imm_metrics_gls1 <- do.call(rbind, imm_metrics_list1_gls)


