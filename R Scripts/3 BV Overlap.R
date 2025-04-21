################################################################
############ BIRD-VESSEL OVERLAP SCRIPT ########################
################################################################

# Takes in bird trajectories (GPS) and Global Fishing Watch daily effort layers
# Produces number of hours of potential overlap per day, per bird

#################### TABLE OF CONTENTS ####################

# 1. Libraries, WD, Global Objects, and Data
# 2. Downsample data
# 3. POINT-BASED OVERLAP FUNCTIONS
# 4. Loop through BBA Tracks
# 5. Point-based overlap results
# 6. Inspect and visualise data



################### 1. Libraries, WD, Global Objects, and Data #########################

# Packages
library(tidyverse)
library(lubridate)
library(data.table)
library(move)
library(sf)
library(gfwr)
library(pbapply)

# Set WD
setwd("___________")

# Save GFW API Token (can get for free on GFW website)
GFW_TOKEN <- "___________"

# Read in the data (all tracks in 1 CSV)
BBA_Tracks <- read.csv("Output Dataset Files/BBA_All_Tracks.csv", header = TRUE)

# Set CRS of input data: WGS 1984
WGS1984 <- st_crs(4326)[[2]]




######################## 2. Downsample data ###############################
# Goal: Get to a temporal resolution that roughly matches 1deg x 1deg (~1 km x 1 km) when travelling at speed
# For black-browed albatrosses this was 10 mins

# Convert to move object
BBA_mv <-  move(x = BBA_Tracks$Longitude,
                y = BBA_Tracks$Latitude,
                time = as.POSIXct(BBA_Tracks$DateTime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"),
                proj = WGS1984,
                animal = as.factor(BBA_Tracks$Bird),
                data = BBA_Tracks,
                removeDuplicatedTimestamps = TRUE)


# Calculate speed (m/s), put into BBA_Tracks dataframe (speed distribution code below)
speed_list <- speed(BBA_mv)
for(i in 1:(length(speed_list))){
  speed_list[[i]] <- c(0, speed_list[[i]])
}
BBA_Tracks <- BBA_Tracks %>% mutate(speed_ms = unlist(speed_list))


# Interpolate tracks to 1 fix per 10 mins, return as a dataframe
for(i in 1:n.indiv(BBA_mv)){
  BBA_mv_interp10 <- interpolateTime(BBA_mv[[i]], 
                                  time = as.difftime(10, units = "mins"), 
                                  spaceMethod = "greatcircle")
  if(i == 1){
    BBA_interp10 <- as.data.frame(BBA_mv_interp10)
  } else {
    BBA_interp10 <- rbind(BBA_interp10, as.data.frame(BBA_mv_interp10))
  }
}

# Clean up the dataframe: subset to relevant columns, specify interpolated positions
BBA_interp10 <- BBA_interp10[ , c("Bird", "GPS_Schedule", "timestamps", "coords.x2", "coords.x1", "TDR", "GLS_Back", "GLS_Leg", "GPS_Tracks", "sensor")]
BBA_interp10$sensor <- ifelse(BBA_interp10$sensor == "interpolateTime", "Interpolated", "Real")

# Rename columns, format columns, calculate Birdhour (number of hours bird spent at each point)
BBA_interp10 <- BBA_interp10 %>%
  rename(Longitude = coords.x1, Latitude = coords.x2, DateTime = timestamps, Interpolated = sensor) %>%
  mutate(DateTime = ymd_hms(DateTime), Bird = as.factor(Bird)) %>%
  mutate(birdhour = 10/60) # After interpolateTime, all time lags should be equal to 10 mins

# Create track list
track_list <- split(BBA_interp10, f = BBA_interp10$Bird)





####################### 3. POINT-BASED OVERLAP FUNCTIONS ############################

# FUNCTION to get aeqd projection, centered on center of track(s)
# "track" should be a single dataframe of 1 or more bird tracks with at least the following columns:
# Bird, DateTime, Latitude, Longitude
get.aeqd.proj <- function(track){
  
  # Set CRS of input data: WGS 1984
  WGS1984 <- st_crs(4326)[[2]]
  
  # Convert to move object
  B_mv <-  move(x = track$Longitude,
                y = track$Latitude,
                time = as.POSIXct(track$DateTime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"),
                proj = WGS1984,
                animal = as.factor(track$Bird),
                data = track,
                removeDuplicatedTimestamps = TRUE)
  
  # Get aeqd projection from original MoveStack object, centered on center of tracks
  B_mv <- spTransform(B_mv, center = TRUE)
  aeqd <- st_crs(B_mv@proj4string)
  return(aeqd)
}


# FUNCTION to get point-based overlap
  # Returns hours of overlap within distance threshold of each bird point (1 column per threshold)
  # Resolution of fisheries data: 1 deg x 1 deg x 1 day
# "track" must be a dataframe of 1 Bird track (ideally regular) with at least the following columns:
  # Bird, DateTime, Latitude, Longitude, birdhour
# "proj" is the projection for distance measurements, recommended to be aeqd from function above
  # If not provided, then calculate aeqd within the function itself (more efficient to provide it up front)
# "dist.threshold" is the distance in meters from each bird point that vessel effort is calculated
  # Can be a vector. Should ideally multiple of 1000.
# "GFW_TOKEN" is your unique GFW token to access fisheries data from Global Fishing Watch

# START OF FUNCTION
get.point.based.overlap <- function(track, proj = NA, dist.threshold, GFW_TOKEN){
  
  # Set CRS of input data: WGS 1984
  WGS1984 <- st_crs(4326)[[2]]
  
  # If no proj provided, calculate aeqd
  if(is.na(proj)){
    proj <- get.aeqd.proj(track)
  }
  
  # Create point layer of track, project into Azimuthal Equi-Distance Projection so that map units are meters
  track_points <- st_as_sf(track, coords = c("Longitude", "Latitude"), crs = WGS1984)
  track_points <- st_transform(track_points, crs = proj)
  
  # Create polygon of buffer around track (radius of max distance threshold + 1km), transform back to WGS1984
  track_buffer <- st_union(st_buffer(track_points, dist = max(dist.threshold) + 1000, endCapStyle = "ROUND"))
  track_buffer <- st_transform(track_buffer, crs = WGS1984) %>% st_sf()
  
  # Create vector of days tracked, find min and max
  start_date <- as.character(date(min(track$DateTime)))
  end_date <- as.character(date(max(track$DateTime)) + days(1))
  
  # Download GFW CSV of apparent fishing effort during those days, within shapefile
  effort <- get_raster(spatial_resolution = 'HIGH',
                       temporal_resolution = 'DAILY',
                       group_by = 'VESSEL_ID',
                       start_date = start_date,
                       end_date = end_date,
                       region = track_buffer,
                       region_source = "USER_SHAPEFILE",
                       key = GFW_TOKEN)
  
  # Abort function if no fishing nearby
  if(nrow(effort) > 0){
    
    # Convert effort df to points layer (saving time by not converting to raster)
    effort_points <- st_as_sf(effort, coords = c("Lon", "Lat"), crs = WGS1984)
    effort_points <- st_transform(effort_points, crs = proj)
    
    # Overlap hours within distance threshold(s), for each day
    overlap_list <- lapply(min(date(track_points$DateTime)):max(date(track_points$DateTime)), FUN = function(track_date){
      
      # Narrow GFW data to just the effort on that day
      effort_on_date <- effort_points[effort_points$`Time Range` == track_date, ]
      
      # Spatial join effort within distance threshold, sum effort hours
      # Overlap within dist.threshold of each location = daily effort in hours / 24, multiply by birdhour
      tracks_on_date <- lapply(dist.threshold, function(d){
        track_on_date <- track_points[date(track_points$DateTime) == track_date, ] %>%
          st_join(effort_on_date, join = st_is_within_distance, dist = d) %>%
          group_by(DateTime) %>%
          summarise(effort = sum(`Apparent Fishing Hours`, na.rm = TRUE),
                    overlap = (sum(`Apparent Fishing Hours`, na.rm = TRUE)/24) * mean(birdhour)) %>%
          st_drop_geometry() %>%
          as.data.frame()
        return(track_on_date)
      })
      
      # Join data together for each distance threshold
      new_col_names <- c()
      joined_df <- data.frame(DateTime = tracks_on_date[[1]]$DateTime)
      for(i in 1:length(dist.threshold)){
        
        # Join together
        joined_df <- left_join(joined_df, tracks_on_date[[i]], by = "DateTime")
        
        # Rename columns
        effort_name <- paste0("effort_", format(dist.threshold[i]/1000, trim = TRUE, scientific = FALSE), "km")
        overlap_name <- paste0("overlap_", format(dist.threshold[i]/1000, trim = TRUE, scientific = FALSE), "km")
        colnames(joined_df)[ncol(joined_df)-1] <- effort_name
        colnames(joined_df)[ncol(joined_df)] <- overlap_name
        
      }
      
      # Return joined dataframe
      return(joined_df)
      
    })
    
    # Bind the dataframes together, join to track dataframe
    overlap_df <- do.call(rbind, overlap_list) %>%
      right_join(track, by = "DateTime") %>%
      relocate(Bird, DateTime, Latitude, Longitude)
    
    # Else, reproduce track dataframe with overlap columns set to 0
  } else {
    
    # Create dataframe with new column names
    new_col_names <- c()
    joined_df <- data.frame(DateTime = track$DateTime)
    for(i in 1:length(dist.threshold)){
      
      # Join together
      joined_df <- joined_df %>% mutate(effort = 0, overlap = 0)
      
      # Rename columns
      effort_name <- paste0("effort_", format(dist.threshold[i]/1000, trim = TRUE, scientific = FALSE), "km")
      overlap_name <- paste0("overlap_", format(dist.threshold[i]/1000, trim = TRUE, scientific = FALSE), "km")
      colnames(joined_df)[ncol(joined_df)-1] <- effort_name
      colnames(joined_df)[ncol(joined_df)] <- overlap_name
      
    }
    
    # Create output dataframe
    overlap_df <- joined_df %>%
      right_join(track, by = "DateTime") %>%
      relocate(Bird, DateTime, Latitude, Longitude)
    
  }
  
  # Return dataframe
  message(paste0(track$Bird[1], " done."))
  return(overlap_df)
}



# FUNCTION to get total overlap for each bird
# "track_list" is a list of tracks with at least 1 overlap column (eg "overlap_5")
# ... should be 1 or more column names to group by (not as character)
get.total.overlap <- function(track_list, ...){
  
  # Combine into single dataframe
  tracks_df <- do.call(rbind, track_list)
  
  # Get all overlap columns (equal to number of distance thresholds)
  overlap_cols <- colnames(tracks_df)[which(grepl("overlap_", colnames(tracks_df)))]
  
  # Produce summary dataframe
  output_df <- tracks_df %>% group_by(...) %>% summarise(tmp = 0)
  for(name in overlap_cols){
    
    # Add a general "overlap" column that equals the focus column
    col1 <- which(colnames(tracks_df) == name)
    tracks_df$overlap <- tracks_df[,col1]
    
    # Sum overlap for groups
    tmp_df <- tracks_df %>%
      group_by(...) %>%
      summarise(tmp = sum(overlap))
    
    # Add to output dataframe
    if(name == overlap_cols[1]){
      output_df <- tmp_df
    } else {
      output_df <- cbind(output_df, tmp_df[, ncol(tmp_df)])
    }
    
    # Correct column name
    colnames(output_df)[which(colnames(output_df) == "tmp")] <- name
  }
  
  # Return output df, arranged by overlap
  output_df <- output_df %>% arrange(desc(!!sym(overlap_cols[1])))
  return(as.data.frame(output_df))
  
}










################### 4. Loop through BBA Tracks ##############################
# These cannot run in parallel because GFW API refuses to process simultaneous requests
# Should return a list of dataframes with effort and overlap values

# Calculate vessel overlap (on same day) within 5km and 30km
final_list <- pblapply(track_list, get.point.based.overlap, dist.threshold = c(5000, 30000), GFW_TOKEN = GFW_TOKEN)

# Save list as RDS file
#saveRDS(final_list, "Output Dataset Files/final_overlap_list.RData")







################ 5. Point-based overlap results ##################

# If needed, load final_list RDS
# final_list <- readRDS("Output Dataset Files/final_overlap_list.RData")

# Get total overlap for each bird
get.total.overlap(final_list, Bird)

# Get total overlap for each date
# We used this list to make a query to GFW for raw AIS data, through BirdLife International
# We sent the list of dates alongside the tracks, to query nearby vessels (final_list)
get.total.overlap(final_list, date(DateTime))




################ 6. Inspect and visualise data ###################

# Plot move object (no background)
plot(BBA_mv, xlab="Longitude", ylab="Latitude", type="l", pch=16, lwd=0.5)
points(BBA_mv, pch=20, cex=0.5)

# Show distribution of speeds
BBA_Tracks %>% filter(speed_ms > 0) %>%
  summarise(min = min(speed_ms), max = max(speed_ms), mean = mean(speed_ms), median = median(speed_ms))
BBA_Tracks %>%
  filter(speed_ms > 0) %>%
  ggplot(aes(x = speed_ms, fill = Bird)) + 
  geom_histogram(binwidth = 1) +
  facet_grid(rows = vars(Bird))


# Function to plot a time series of overlap
# "track_list" is a list of tracks with at least 1 overlap column (eg "overlap_5")
# "overlap_column" is the name of the overlap column to look at over time (character)
# "colour" is how to colour the lines (eg "BirdID")
plot.overlap.over.time <- function(track_list, overlap_column, colour = "BirdID"){
  
  # Combine into single dataframe
  tracks_df <- do.call(rbind, track_list)
  
  # Plot overlap over time
  plot <- ggplot(tracks_df, aes(x = DateTime, y = .data[[overlap_column]], colour = .data[[colour]])) + 
    geom_line() + 
    labs(title = paste(min(tracks_df$DateTime), " to ", max(tracks_df$DateTime), sep = "")) +
    theme_bw()
  return(plot)
}



# Function to plot a map of overlap
# "track_list" is a list of tracks with at least 1 overlap column (eg "overlap_5")
# Plus "Latitude", "Longitude", "DateTime", and a colour column like "BirdID"
# "overlap_column" is the name of the overlap column to look at over time (character)
# "colour" is how to colour the lines (eg "BirdID")
# "continent" can be a vector of continents to pass to ne_countries
# "xlim" and "ylim" are optional vectors of Longitude and Latitude limits respectively.
# By default, will use max and min Lat and Long of the whole dataset
# "polygon" should be a simple feature containing land polygons (eg from rnaturalearth::ne_countries)
plot.overlap.map <- function(track_list, overlap_column, colour = "BirdID", adjust.for.IDL = FALSE,
                             xlim = NA, ylim = NA, polygon){
  
  # Combine into single dataframe
  tracks_df <- do.call(rbind, track_list)
  
  # Adjust longitude to account for IDL if necessary
  if(adjust.for.IDL == TRUE){
    tracks_df <- tracks_df %>%
      mutate(Longitude_adj = ifelse(Longitude < 0, Longitude + 360, Longitude))
    polygon <- st_shift_longitude(polygon)
  } else {
    tracks_df <- tracks_df %>% mutate(Longitude_adj = Longitude)
  }
  
  
  # Plot overlap over space
  plot <- ggplot(data = polygon) + 
    geom_sf(color = "black", fill = "gray83") +
    geom_path(data = tracks_df, aes(x = Longitude_adj, y = Latitude, colour = .data[[colour]]), alpha = 0.4) +
    geom_point(data = tracks_df[tracks_df[, overlap_column] > 0, ], 
               aes(x = Longitude_adj, y = Latitude, colour = .data[[colour]], size = .data[[overlap_column]])) +
    labs(x = "Longitude", y = "Latitude", 
         title = paste(min(tracks_df$DateTime), " to ", max(tracks_df$DateTime), sep = "")) +
    theme_bw()
  if(is.na(xlim) & is.na(ylim)){
    plot <- plot +
      coord_sf(xlim = c(min(tracks_df$Longitude_adj), max(tracks_df$Longitude_adj)), 
               ylim = c(min(tracks_df$Latitude), max(tracks_df$Latitude)), 
               expand = TRUE)
  } else {
    plot <- plot +
      coord_sf(xlim = xlim[1], 
               ylim = ylim[1], 
               expand = FALSE)
  }
  return(plot)
}

# Load land polygon
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
sa <- ne_countries(scale = "large", continent = "South America")

# Plot time series and map
plot.overlap.over.time(final_list, overlap_column = "overlap_5km", colour = "Bird")
plot.overlap.map(final_list, overlap_column = "overlap_5km", colour = "Bird",
                 polygon = sa)



