###################################################
############## GLS READ IN ########################
###################################################

# Reads in GLS-immersion data
# Joins it to GPS tracks

#################### TABLE OF CONTENTS ####################

# 1. Packages & WD
# 2. Read in GLS immersion data
# 3. Join immersion data to GPS data



############## 1. Packages & WD ###############################

# Set WD
setwd("___________")

# Packages
library(tidyverse)
library(lubridate)



############# 2. READ IN GLS IMMERSION DATA #####################

# Specify where immersion files live
wet_files <- list.files(path = "Immersion CSVs")

# Split the filename by the '_' character so as we add metadata to the GLS trace
splits <- strsplit(wet_files, "_")

# Add metadata to parameter list
param_gls <- list()
param_gls$Bird <- sapply(splits, function(x){
  colour <- gsub("\\d.*", "", str_to_lower(x[5])) # Extract colour from name
  number <- sub("\\D+", "", str_to_lower(x[5])) # Extract number from name
  new_name <- ifelse(nchar(number) > 0, 
                     paste(toupper(number), colour, sep = " "), # Flip, recombine
                     colour) # Except for the unringed bird (displacement trial)
  return(new_name)
})
param_gls$colony <- sapply(splits, "[", 1)
param_gls$GLSloc <- sapply(splits, function(x){sub("GLSIMM", "", x[2])})
param_gls$GLSid <- sapply(splits, "[", 3)
param_gls$nest <- sapply(splits, "[", 4)
param_gls$DateOn <- sapply(splits, function(x){dmy(sub("on", "", x[6]))})
param_gls$DateOff <- sapply(splits, function(x){dmy(sub("off", "", x[8]))})

# Loop through wet_files to read in immersion data
# (Ignore warning message that says header and col.names are of different lengths)
wet_data <- lapply(1:length(wet_files), function(i){
  
  # Read in .deg data like it's a CSV
  # (Suppress warning message that says header and col.names are of different lengths)
  suppressWarnings({
    z <- read.csv(paste("Immersion CSVs/", wet_files[i], sep = ""), 
                  col.names = c("DateTime", "duration", "state"), 
                  sep = "\t", fill = TRUE)
  })
  
  # Remove metadata rows, convert data types, add metadata
  z <- z[-(1:which(z$DateTime == "DD/MM/YYYY HH:MM:SS")),] %>%
    mutate(DateTime = dmy_hms(DateTime), 
           duration = as.numeric(duration), 
           state = as.factor(state),
           colony = param_gls$colony[i],
           GLSloc = param_gls$GLSloc[i],
           GLSid = param_gls$GLSid[i],
           nest = param_gls$nest[i],
           Bird = param_gls$Bird[i],
           DateOn = param_gls$DateOn[i],
           DateOff = param_gls$DateOff[i])
  
  # Return the dataframe
  return(z)
  
})



# Save the immersion data to an RDS file
#saveRDS(wet_data, "Output Dataset Files/Wet_Data.RData")

# NOTE: 57L Yellow files were mislabelled as 27L Yellow. I changed all to 57L Yellow.






################## 3. JOIN IMMERSION DATA TO GPS DATA #####################

# Read in GPS tracks
# Read in the data (all tracks in 1 CSV)
BBA_Tracks <- read.csv("Output Dataset Files/BBA_All_Tracks.csv", header = TRUE) %>%
  mutate(DateTime = ymd_hms(DateTime))

# Join immersion data to BBA_Tracks
BBA_wet_list <- lapply(levels(as.factor(BBA_Tracks$Bird)), function(x){
  
  # Full df, 1 bird
  df1 <- BBA_Tracks[BBA_Tracks$Bird == x, ]
  
  # Wet df, 1 bird, LEG ONLY
  which_bird <- which(sapply(1:length(wet_data), function(y){
    (wet_data[[y]]$Bird[1] == x) & (wet_data[[y]]$GLSloc[1] == "leg")
  }))
  if(length(which_bird) > 0){
    wet_df <- wet_data[[which_bird]][, c("DateTime", "state")]
    
    # Join dataframes, order by DateTime, fill up the states, then remove excess rows
    df2 <- full_join(df1, wet_df, by = "DateTime") %>%
      rename(leg_wet = state) %>%
      arrange(DateTime) %>%
      fill(leg_wet, .direction = "updown") %>%
      filter(!is.na(Bird))
    
  } else {
    df2 <- df1 %>% mutate(leg_wet = NA)
  }
  
  # Wet df, 1 bird, BACK ONLY
  which_bird <- which(sapply(1:length(wet_data), function(y){
    (wet_data[[y]]$Bird[1] == x) & (wet_data[[y]]$GLSloc[1] == "back")
  }))
  if(length(which_bird) > 0){
    wet_df <- wet_data[[which_bird]][, c("DateTime", "state")]
    
    # Join dataframes, order by DateTime, fill up the states, then remove excess rows
    df3 <- full_join(df2, wet_df, by = "DateTime") %>%
      rename(back_wet = state) %>%
      arrange(DateTime) %>%
      fill(back_wet, .direction = "updown") %>%
      filter(!is.na(Bird))
    
  } else {
    df3 <- df2 %>% mutate(back_wet = NA)
  }
  
  return(df3)
})

# Combine into one dataframe
BBA_wet <- do.call(rbind, BBA_wet_list)

# Write dataframe to CSV
write.csv(BBA_wet, "Output Dataset Files/BBA_All_Tracks_wImmersion.csv")







# 
# 
# ############## CONVERT IMMERSION DATA TO REGULAR TIME SERIES ###########
# 
# # If necessary, read in immersion data created above
# #wet_data <- readRDS("Output Dataset Files/Wet_Data.RData")
# 
# # Narrow down just to leg GLS (remove back GLS)
# sapply(wet_data, function(x){x$GLSloc[1]})
# leg_data <- wet_data
# 
# # Ensure resolution is 6 seconds
# sapply(leg_data, function(x){min(x$duration)})
# 
# # For each bird...
# imm_reg <- lapply(leg_data, function(x){
#   
#   # Find start and end DateTimes
#   start <- x$DateTime[1] - seconds(x$duration[1])
#   end <- x$DateTime[nrow(x)]
#   
#   # Create data.frame of DateTimes every 6 seconds, join immersion data.frame
#   DTs <- data.frame(DateTime = seq(from = start, to = end, by = "6 secs"))
#   new_df <- left_join(DTs, x, by = "DateTime")
#   
#   # Fill in wet/dry values, remove unnecessary columns
#   new_df <- new_df[, c("Bird", "DateTime", "state")] %>%
#     rename(leg_wet = state) %>%
#     fill(Bird, .direction = "updown") %>%
#     fill(leg_wet, .direction = "updown")
#   return(new_df)
# })
# 
# # Save to RDS file
# saveRDS(imm_reg, "Output Dataset Files/Wet_Data_Reg.RData")
# 
# 



