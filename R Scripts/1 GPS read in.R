##############################################################
########## READ IN FK BBA DATA AND VISUALISE #################
##############################################################

# Reads in black-browed albatross GPS tracks 


#################### TABLE OF CONTENTS ####################

# 1. Packages and WD
# 2. Read in the tracks




########################## 1. Packages and WD ####################

# Packages
library(tidyverse)
library(data.table)
library(lubridate)


# Set WD
setwd("___________")




####################### 2. Read in the tracks ######################

# Read in Summary of Deployments
summary <- read.csv("SUMMARY OF DEPLOYMENTS CSV.csv", header = TRUE)

# Create file list of BBA CSVs
files <- list.files(path = "Indiv GPS Tracks - CSV", pattern = "*.csv")

# Create function to add metadata to a track
add.metadata <- function(x){
  track <- read.csv(paste("Indiv GPS Tracks - CSV", x, sep = "/"))
  
  # Extract bird name and find its index in "summary"
  name <- str_to_lower(strsplit(x, "[_]")[[1]][3]) # Extract bird name
  colour <- gsub("\\d.*", "", name) # Extract colour from name
  number <- sub("\\D+", "", name) # Extract number from name
  new_name <- ifelse(nchar(number) > 0, 
                     paste(toupper(number), colour, sep = " "), # Flip, recombine
                     colour) # Except for the unringed bird (displacement trial)
  row <- which(summary$BIRD == new_name) # Find row number of bird
  
  # Add metadata columns
  track <- track %>% 
    mutate(Bird = summary[row, "BIRD"],
           GPS_Schedule = summary[row, "GPS.schedule"], 
           TDR = summary[row, "TDR"], 
           GLS_Back = summary[row, "GLS.Back"], 
           GLS_Leg = summary[row, "GLS.leg"], 
           GPS_Tracks = summary[row, "GPS.tracks"]) %>%
    relocate(Bird, GPS_Schedule)
  
  # Return dataframe
  return(track)
}

# Add metadata to all files and read into list
tracks_list <- lapply(files, add.metadata)

# Row bind into a single dataframe
BBA_Tracks <- do.call(rbind, tracks_list)

# Detect and correct dates with different formats
ymd <- ymd(BBA_Tracks$Date)
dmy <- dmy(BBA_Tracks$Date)
dmy[is.na(dmy)] <- ymd[is.na(dmy)]
BBA_Tracks$Date <- dmy

# Add a DateTime column
BBA_Tracks <- BBA_Tracks %>%
  mutate(DateTime = ymd_hms(paste(BBA_Tracks$Date, BBA_Tracks$Time, sep = ""))) %>%
  relocate(DateTime, .after = Time) %>%
  arrange(GPS_Schedule, Bird, DateTime)

# Write to CSV and View
write.csv(BBA_Tracks, "Output Dataset Files/BBA_All_Tracks.csv")
#View(BBA_Tracks)

