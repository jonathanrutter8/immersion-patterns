#####################################################################
################## BBA BIRD VESSEL INTERACTION ######################
#####################################################################

# Identifies bird-vessel interactions (periods of close proximity)
# Produces maps and animations
# Note that vessel AIS tracks from GFW must be downloaded manually from the GFW website, 
# or queried from the organisation. 
# We have provided 2 simulated tracks for demonstration purposes.

################# CONTENTS #########################

# 1. Packages, WD, Global Objects, and Data
# 2. Function to detect bird-vessel interactions
# 3. Run the function in parallel
# 4. Summarise interactions
# 5. Map interactions
# 6. Animate interactions
# 7. Test and visualise data further
# 8. Explore maps further





################### 1. Packages, WD, Global Objects, and Data #########################

# Packages
library(tidyverse)
library(lubridate)
library(move)
library(pbapply)
library(parallel)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(sf)

# Set WD
wd <- "___________"
setwd(wd)

# Read in the bird data (all tracks in 1 CSV)
BBA_Tracks <- read.csv("Output Dataset Files/BBA_All_Tracks_wImmersion.csv", header = TRUE) %>% 
  mutate(DateTime = ymd_hms(DateTime))

# Read in the vessel data (here, just a small redacted set of downsampled points for demonstration)
AIS_Tracks <- read.csv("Output Dataset Files/BBA_gfw_vessels_redacted.csv", header = TRUE) %>%
  mutate(Bird = as.factor(Bird), vessel_class = as.factor(vessel_class), fishing = as.factor(fishing), SubTraj = as.factor(SubTraj),
         #flag = as.factor(flag), ssvid = as.factor(ssvid)
         ) %>%
  mutate(DateTime = ymd_hms(DateTime, tz = "UTC"))

# Set CRS of input data: WGS 1984
WGS1984 <- st_crs(4326)[[2]]

# Prep for parallel processing
numCores <- 2 #detectCores()
cl <- makeCluster(numCores)
clusterExport(cl = cl, c("BBA_Tracks", "AIS_Tracks", "WGS1984", "wd"))
clusterEvalQ(cl, {
  
  # Load packages
  library(tidyverse)
  library(lubridate)
  library(move)
  library(sf)
  
  # Set WD
  setwd(wd)
})





################# 2. FUNCTION TO DETECT BIRD-VESSEL INTERACTION ###############
# Input bird ID, uses full bird and vessel datasets (here BBA_Tracks, AIS_Tracks)
# Output list of bird tracks with number of vessels encountered, distance to nearest

get.interactions <- function(bird){
  
  # Subset bird (B) data to bird
  B <- BBA_Tracks %>% filter(Bird == bird)
  
  # Subset vessel (V) data to vessels near the bird
  V <- AIS_Tracks %>% filter(Bird == bird)
  
  # If there is at least 1 vessel trajectory in vicinity of bird
  if(nrow(V) > 1){
    
    # Convert vessel data to move object
    V_mv <- move(x = V$Longitude,
                 y = V$Latitude,
                 time = as.POSIXct(V$DateTime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"),
                 proj = WGS1984,
                 animal = as.factor(V$SubTraj),
                 data = V,
                 removeDuplicatedTimestamps = TRUE)
    
    # Get aeqd projection from MoveStack object, centered on center of tracks
    V_mv_aeqd <- spTransform(V_mv, center = TRUE)
    aeqd <- V_mv_aeqd@proj4string
    
    # Get vector of DateTimes
    B_DateTimes <- B$DateTime
    
    # For each vessel, calculate distance from bird and get vessel metadata
    nearby_vessels <- lapply(1:n.indiv(V_mv), function(i){
      
      # If only 1 indiv
      if(n.indiv(V_mv) == 1){
        V_mv_i <- V_mv
      } else {
        V_mv_i <- V_mv[[i]]
      }
      
      # If that iteration's vessel has a trajectory in a move object
      if(class(V_mv_i) == "Move"){
        
        # Get vector of DateTimes
        V_DateTimes_all <- V_mv_i@timestamps
        start <- min(V_DateTimes_all)
        end <- max(V_DateTimes_all)
        V_DateTimes <- B_DateTimes[(B_DateTimes >= start) & (B_DateTimes <= end)]
        
        # If there are at least 2 B DateTimes that fall within V DateTime range
        if(length(V_DateTimes) > 1){
          
          # Interpolate vessel data to DateTimes from bird track
          V_mv_interp <- interpolateTime(V_mv_i, 
                                         time = V_DateTimes, 
                                         spaceMethod = "greatcircle")
          
          # Convert back to dataframe, subset and rename columns
          V_interp <- as.data.frame(V_mv_interp) %>%
            dplyr::select(timestamps, coords.x2, coords.x1) %>%
            rename(Longitude = coords.x1, Latitude = coords.x2, DateTime = timestamps) %>%
            mutate(DateTime = ymd_hms(DateTime))
          
          # Create point layer of vessel track, transform to aeqd projection
          V_points <- st_as_sf(V_interp, 
                               coords = c("Longitude", "Latitude"), crs = WGS1984)
          V_points <- st_transform(V_points, crs = aeqd)
          
          # Create point layer of bird track near vessel
          B_points <- st_as_sf(B[B$DateTime %in% V_DateTimes, c("DateTime", "Latitude", "Longitude")], 
                               coords = c("Longitude", "Latitude"), crs = WGS1984)
          B_points <- st_transform(B_points, crs = aeqd)
          
          # Get distance (in m) between B and V points at each DateTime
          distances <- as.numeric(st_distance(B_points, V_points, by_element = TRUE))
          
          # Get a df of fishing activity for interpolated DateTimes
          if("fishing" %in% colnames(V_mv_i@data)){
            fishing_df1 <- data.frame(DateTime = V_DateTimes_all, fishing = V_mv_i@data$fishing, df = 1)
          } else if("fishing" %in% colnames(V_mv_i@idData)){
            fishing_df1 <- data.frame(DateTime = V_DateTimes_all, fishing = V_mv_i@idData$fishing, df = 1)
          }
          fishing_df2 <- data.frame(DateTime = V_DateTimes, fishing = NA, df = 2)
          fishing_df3 <- rbind(fishing_df1, fishing_df2) %>%
            arrange(DateTime) %>%
            fill(fishing, .direction = "downup")
          fishing_df3 <- fishing_df3[fishing_df3$DateTime %in% V_DateTimes, -3]
          fishing_df3 <- fishing_df3[!duplicated(fishing_df3$DateTime), ]
          
          # Put distances into dataframe and add vessel metadata
          distances_df <- data.frame(DateTime = V_DateTimes,
                                     distance = distances,
                                     dist_class = as.character(rep(NA, length(V_DateTimes))),
                                     SubTraj = as.character(rep(V_mv_i@idData$SubTraj, length(V_DateTimes))),
                                     #ssvid = as.character(rep(V_mv_i@idData$ssvid, length(V_DateTimes))),
                                     #vesselID = as.character(rep(V_mv_i@idData$vesselID, length(V_DateTimes))),
                                     #flag = as.character(rep(V_mv_i@idData$flag, length(V_DateTimes))),
                                     vessel_class = as.character(rep(V_mv_i@idData$vessel_class, length(V_DateTimes))),
                                     class_fishing = as.character(rep(ifelse(V_mv_i@idData$vessel_class %in% c("fishing", "fixed_gear", "set_longlines", "squid_jigger", "trawlers"), 
                                                                             "Fishing", "Non-Fishing"), length(V_DateTimes))),
                                     fishing = as.character(fishing_df3$fishing))
          distances_df$dist_class <- ifelse(distances_df$distance <= 1000, "1 km", 
                                            ifelse(distances_df$distance <= 5000, "5 km", 
                                                   ifelse(distances_df$distance <= 30000, "30 km", "No Interaction")))
          
          # Return that dataframe
          return(distances_df)
          
        # Else if there are no DateTimes in range, return NULL  
        } else {
          return(NULL)
        }
        
      # Else if there is no trajectory for that iteration's vessel, return NULL
      } else {
        return(NULL)
      }
      
        
      
      
    })
    
    # Add empty columns to Bird Dataframe
    B <- B %>% mutate(SubTraj = as.character(NA), distance = as.numeric(NA), dist_class = "No Interaction", 
                      #ssvid = as.character(NA), vesselID = as.character(NA), flag = as.character(NA), 
                      vessel_class = as.character(NA), class_fishing = as.character(NA), fishing = as.character(NA))
    
    # For each DateTime, fill in info of nearest vessel
    for(j in B$DateTime){
      
      # Get vector of distances to each vessel at that DateTime
      dists <- sapply(1:n.indiv(V_mv), function(k){
        if(!is.null(nearby_vessels[[k]])){
          if(j %in% nearby_vessels[[k]]$DateTime){
            return(nearby_vessels[[k]][nearby_vessels[[k]]$DateTime == j, "distance"])
          } else {
            return(as.numeric(NA))
          }
        } else {
          return(as.numeric(NA))
        }
      })
      
      # Figure out which SubTraj is nearest
      if(sum(!is.na(dists)) > 0){
        SubTraj_index <- which.min(dists)
        DT_index <- which(nearby_vessels[[SubTraj_index]]$DateTime == j)
      } else { # If no nearby SubTrajs, skip to next DateTime
        next
      }
      
      # Add values to Bird Dataframe
      B$SubTraj[B$DateTime == j] <- nearby_vessels[[SubTraj_index]]$SubTraj[DT_index]
      B$distance[B$DateTime == j] <- dists[SubTraj_index]
      B$dist_class[B$DateTime == j] <- nearby_vessels[[SubTraj_index]]$dist_class[DT_index]
      #B$ssvid[B$DateTime == j] <- nearby_vessels[[SubTraj_index]]$ssvid[DT_index]
      #B$vesselID[B$DateTime == j] <- nearby_vessels[[SubTraj_index]]$vesselID[DT_index]
      #B$flag[B$DateTime == j] <- nearby_vessels[[SubTraj_index]]$flag[DT_index]
      B$vessel_class[B$DateTime == j] <- nearby_vessels[[SubTraj_index]]$vessel_class[DT_index]
      B$class_fishing[B$DateTime == j] <- nearby_vessels[[SubTraj_index]]$class_fishing[DT_index]
      B$fishing[B$DateTime == j] <- nearby_vessels[[SubTraj_index]]$fishing[DT_index]
    }
    
    # Add columns for numeric GPS resolution and Interaction Index
    B <- B %>% mutate(resolution = with(., case_when(
      (GPS_Schedule == "1 sec") ~ 1,
      (GPS_Schedule == "5 sec") ~ 5,
      (GPS_Schedule == "10 sec") ~ 10,
      (GPS_Schedule == "1 min") ~ 60,
      (GPS_Schedule == "7 min") ~ 420
    ))) %>%
      mutate(interaction = as.numeric(NA))
    
    # Start interaction index and duration index
    int_index <- 0
    dur_index <- 0
    
    # How many fixes equals 5 minutes?
    min_fixes <- 300/B$resolution[1]
    
    # Look for interactions of <1 km that last for > 5 mins, number them
    for(i in 1:nrow(B)){
      
      # If it's not an interaction
      if(is.na(B[i, "dist_class"]) | B[i, "dist_class"] != "1 km"){
        
        # If the interaction just ended, but it's not long enough
        if(dur_index < min_fixes & dur_index > 0){
          dur_index <- 0
          
          # Else if the interaction just ended, and it's long enough
        } else if(dur_index >= min_fixes){
          int_index <- int_index + 1
          B[int_start:(i-1), "interaction"] <- int_index
          dur_index <- 0
          
          # Else, it's not part of an interaction
        } else {
          next
        }
      }
      
      # If beginning of an interaction
      if(dur_index == 0 & B[i, "dist_class"] == "1 km"){
        int_start <- i
        dur_index <- dur_index + 1
        
        # Else if last row and there is an interaction ongoing
      } else if(i == nrow(B) & B[i, "dist_class"] == "1 km"){
        # If interaction is long enough
        if(dur_index >= min_fixes){
          dur_index <- dur_index + 1
          int_index <- int_index + 1
          B[int_start:i, "interaction"] <- int_index
        }
        
        # Else if middle of an interaction
      } else if(dur_index > 0 & B[i, "dist_class"] == "1 km" & B[i, "SubTraj"] == B[i-1, "SubTraj"]){
        dur_index <- dur_index + 1
        
        # Backup else (this should never be needed)
      } else {
        dur_index <- 0
      }
    }
  
  # Else there are no vessels nearby    
  } else {
    
    # Add columns but don't fill them up
    B <- B %>% 
      mutate(SubTraj = as.character(NA), distance = as.numeric(NA), dist_class = "No Interaction", 
             #ssvid = as.character(NA), vesselID = as.character(NA), flag = as.character(NA), 
             vessel_class = as.character(NA), class_fishing = as.character(NA), fishing = as.character(NA)) %>%
      mutate(resolution = with(., case_when(
        (GPS_Schedule == "1 sec") ~ 1,
        (GPS_Schedule == "5 sec") ~ 5,
        (GPS_Schedule == "10 sec") ~ 10,
        (GPS_Schedule == "1 min") ~ 60,
        (GPS_Schedule == "7 min") ~ 420
      ))) %>%
      mutate(interaction = as.numeric(NA))
    
  }
  
  # Return bird dataframe vessel and interaction info
  return(B)
}








#################### 3. Run the function in parallel #########################

# Get vector of bird names
birds <- levels(as.factor(BBA_Tracks$Bird))

# Loop through bird tracks (arguments?)
BBA_int_list <- pblapply(birds, get.interactions, cl = cl)

# Close the cluster
stopCluster(cl)

# Write to file
saveRDS(BBA_int_list, "Output Dataset Files/BBA_Track_List_Imm_Int.RData")

# # Save CSVs to folder (useful to explore in ArcGIS Pro)
# for(i in 1:length(BBA_int_list)){
#   name <- gsub(" ", "", BBA_int_list[[i]]$Bird[1])
#   write.csv(BBA_int_list[[i]], file = paste("Track_CSVs_Imm_Int/", name, "_Imm_Int.csv", sep = ""))
# }
# 
# # Save combined CSV to folder
# all_tracks <- do.call(rbind, BBA_int_list)
# write.csv(all_tracks, file = "Track_CSVs_Imm_Int/AllTracks_Imm_Int.csv")






##################### 4. Summarise interactions #####################

# If necessary, read in data
BBA_int_list <- readRDS("Output Dataset Files/BBA_Track_List_Imm_Int.RData")

# Summary of interactions (duration in mins)
int_summary_list <- lapply(BBA_int_list, function(x){
  
  # If there are interactions
  if(sum(!is.na(x$interaction)) > 0){
    
    # Create vectors for each column of summary dataframe
    int_index <- levels(as.factor(x$interaction))
    bird_name <- rep(x$Bird[1], length(int_index))
    #int_dur <- sapply(int_index, function(y){
    #  interacting <- x$interaction[!is.na(x$interaction) & x$interaction == y]
    #  return(length(interacting) * x$resolution[1] / 60)
    #}) # This way of calculating duration sometimes comes up with different numbers
    start <- sapply(int_index, function(y){
      as_datetime(min(x$DateTime[!is.na(x$interaction) & x$interaction == y]), tz = "UTC")
    })
    start <- as_datetime(start)
    end <- sapply(int_index, function(y){
      max(x$DateTime[!is.na(x$interaction) & x$interaction == y])
    })
    end <- as_datetime(end)
    int_dur <- abs(difftime(end, start, units = "mins"))
    resolution <- rep(x$resolution[1], length(int_index))
    SubTraj <- sapply(int_index, function(y){
      x$SubTraj[!is.na(x$interaction) & x$interaction == y][1]
    })
    # ssvid <- sapply(int_index, function(y){
    #   x$ssvid[!is.na(x$interaction) & x$interaction == y][1]
    # })
    # vesselID <- sapply(int_index, function(y){
    #   x$vesselID[!is.na(x$interaction) & x$interaction == y][1]
    # })
    # flag <- sapply(int_index, function(y){
    #   x$flag[!is.na(x$interaction) & x$interaction == y][1]
    # })
    vessel_class <- sapply(int_index, function(y){
      x$vessel_class[!is.na(x$interaction) & x$interaction == y][1]
    })
    class_fishing <- sapply(int_index, function(y){
      x$class_fishing[!is.na(x$interaction) & x$interaction == y][1]
    })
    fishing <- sapply(int_index, function(y){
      vec <- x$fishing[!is.na(x$interaction) & x$interaction == y]
      if(max(vec, na.rm = TRUE) == 1){
        return(1)
      } else {
        return(0)
      }
    })
    
    # Combine vectors into dataframe and return
    return(data.frame(bird_name, int_index, 
                      start, end, int_dur, resolution, SubTraj, #ssvid, vesselID, flag,
                      vessel_class, class_fishing, fishing))
    
  # Else if no interactions, return NULL  
  } else {
    return(NULL)
  }
})
int_summary <- do.call(rbind, int_summary_list)
int_summary

# Save as R object and CSV
saveRDS(int_summary, "Output Dataset Files/BBA_Interactions_Summary.RData")
#write.csv(int_summary, "Output Dataset Files/BBA_Interactions_Summary.csv")







############################ 5. Map Interactions ###################

# Define function to get min and max vessel DateTimes that correspond to bird DateTimes
# Looking for nearest V DateTimes outside min and max B DateTimes
# b is bird track, v is vessel track
get.vdt <- function(b, v){
  
  # Get min V DateTime
  minB <- min(b$DateTime)
  nearest_min <- v$DateTime[which.min(abs(v$DateTime - minB))]
  if(nearest_min > minB){
    min_index <- which.min(abs(v$DateTime - minB)) - 1
    min_index <- ifelse(min_index < 1, 1, min_index)
    minV <- v$DateTime[min_index]
  } else {
    minV <- nearest_min
  }
  
  # Get max V DateTime
  maxB <- max(b$DateTime)
  nearest_max <- v$DateTime[which.min(abs(v$DateTime - maxB))]
  if(nearest_max < maxB){
    max_index <- which.min(abs(v$DateTime - maxB)) + 1
    max_index <- ifelse(max_index > nrow(v), nrow(v), max_index)
    maxV <- v$DateTime[max_index]
  } else {
    maxV <- nearest_max
  }
  
  # Return vector of min and max
  return(c(minV, maxV))
}

# Get world polygons, if haven't already done
world <- ne_countries(scale = "large", returnclass = "sf", continent = "South America")

# Get first and last DateTime of interactions, find nearest vessel DateTimes, plot together
maps <- lapply(BBA_int_list, function(x){
  
  # If there are interactions
  if(sum(!is.na(x$interaction)) > 0){
    
    maps1bird <- lapply(1:max(x$interaction, na.rm = TRUE), function(y){
      
      # Get bird track just for that interaction, with 2 mins buffer on either side (unless at end of track)
      min_index <- max(0, min(which(x$interaction == y)) - ceiling(120/x$resolution[1]))
      max_index <- min(nrow(x), max(which(x$interaction == y) + ceiling(120/x$resolution[1])))
      Btrack <- x[min_index:max_index, ]
      
      # Get first and last DateTime
      min_Bdt <- Btrack$DateTime[1]
      max_Bdt <- Btrack$DateTime[nrow(Btrack)]
      
      # Get vessel track associated with that interaction
      Vtrack <- AIS_Tracks %>% filter(Bird == Btrack$Bird[1], SubTraj == Btrack$SubTraj[1])
      min_Vdt <- get.vdt(Btrack, Vtrack)[1]
      max_Vdt <- get.vdt(Btrack, Vtrack)[2]
      Vtrack <- Vtrack %>% filter(DateTime >= min_Vdt, DateTime <= max_Vdt)
      
      # Get boundary coordinates
      minLong <- min(c(Btrack$Longitude, Vtrack$Longitude))
      maxLong <- max(c(Btrack$Longitude, Vtrack$Longitude))
      rangeLong <- abs(maxLong - minLong)
      minLat <- min(c(Btrack$Latitude, Vtrack$Latitude))
      maxLat <- max(c(Btrack$Latitude, Vtrack$Latitude))
      rangeLat <- abs(maxLat - minLat)
      
      # Plot if there is GLS leg data
      if(sum(!is.na(Btrack$leg_wet)) > 0){
        
        ggplot(data = world) + 
          geom_sf(color = "black", fill = "gray83") +
          geom_point(data = Vtrack, aes(x = Longitude, y = Latitude), size = 0.7) +
          geom_path(data = Vtrack, aes(x = Longitude, y = Latitude, group = SubTraj)) +
          geom_point(data = Btrack, aes(x = Longitude, y = Latitude, colour = leg_wet), size = 0.7) +
          geom_path(data = Btrack, aes(x = Longitude, y = Latitude, colour = leg_wet, group = SubTraj)) +
          coord_sf(xlim = c(minLong, maxLong), ylim = c(minLat, maxLat), expand = TRUE, clip = "off") +
          labs(x = "Longitude", y = "Latitude", colour = "GLS - Leg") +
          annotate(
            geom = "text", x = maxLong + 0.1*rangeLong, y = minLat - 0.12*rangeLat, 
            label = paste("Bird: ", Btrack$Bird[1], "\n",
                          "Resolution: ", Btrack$GPS_Schedule[1], "\n",
                          "Interaction Num: ", max(Btrack$interaction, na.rm = TRUE), "\n",
                          "SubTraj: ", Btrack$SubTraj[1], "\n", 
                          # "MMSI: ", Btrack$ssvid[1], "\n",
                          # "Flag: ", Btrack$flag[1], "\n", 
                          "Class: ", Btrack$vessel_class[1], "\n",
                          "Active fishing: ", ifelse(length(Btrack$fishing[Btrack$fishing == 1]) > 0, "Yes", "No"), "\n",
                          sep = ""), 
            hjust = 0, vjust = 0, size = 3) +
          geom_text(aes(label = with_tz(min_Bdt, "America/Argentina/Buenos_Aires"),
                        x = Btrack$Longitude[1], y = Btrack$Latitude[1]),
                    position = position_dodge(0.05*rangeLong),
                    size = 3) +
          geom_text(aes(label = with_tz(max_Bdt, "America/Argentina/Buenos_Aires"),
                        x = Btrack$Longitude[nrow(Btrack)], y = Btrack$Latitude[nrow(Btrack)]),
                    size = 3,
                    position = position_dodge(0.05*rangeLong)) +
          theme_bw()
        
      # Plot if there is no GLS leg data  
      } else {
        
        ggplot(data = world) + 
          geom_sf(color = "black", fill = "gray83") +
          geom_point(data = Vtrack, aes(x = Longitude, y = Latitude), size = 0.7) +
          geom_path(data = Vtrack, aes(x = Longitude, y = Latitude, group = SubTraj)) +
          geom_point(data = Btrack, aes(x = Longitude, y = Latitude, colour = Bird), size = 0.7) +
          geom_path(data = Btrack, aes(x = Longitude, y = Latitude, colour = Bird, group = SubTraj)) +
          coord_sf(xlim = c(minLong, maxLong), ylim = c(minLat, maxLat), expand = TRUE, clip = "off") +
          labs(x = "Longitude", y = "Latitude") +
          annotate(
            geom = "text", x = maxLong + 0.1*rangeLong, y = minLat - 0.12*rangeLat, 
            label = paste("Bird: ", Btrack$Bird[1], "\n",
                          "Resolution: ", Btrack$GPS_Schedule[1], "\n",
                          "Interaction Num: ", max(Btrack$interaction, na.rm = TRUE), "\n",
                          "SubTraj: ", Btrack$SubTraj[1], "\n", 
                          # "MMSI: ", Btrack$ssvid[1], "\n",
                          # "Flag: ", Btrack$flag[1], "\n", 
                          "Class: ", Btrack$vessel_class[1], "\n",
                          "Active fishing: ", ifelse(length(Btrack$fishing[Btrack$fishing == 1]) > 0, "Yes", "No"), "\n",
                          sep = ""), 
            hjust = 0, vjust = 0, size = 3) +
          geom_text(aes(label = with_tz(min_Bdt, "America/Argentina/Buenos_Aires"),
                        x = Btrack$Longitude[1], y = Btrack$Latitude[1]),
                    position = position_dodge(0.05*rangeLong),
                    size = 3) +
          geom_text(aes(label = with_tz(max_Bdt, "America/Argentina/Buenos_Aires"),
                        x = Btrack$Longitude[nrow(Btrack)], y = Btrack$Latitude[nrow(Btrack)]),
                    size = 3,
                    position = position_dodge(0.05*rangeLong)) +
          theme_bw()
        
      }
    })
    
  # Else if no interactions, return NULL  
  } else {
    return(NULL)
  }
})


# Save map images to folder
for(i in 1:length(maps)){
  
  if(!is.null(maps[[i]])){
    
    for(j in 1:length(maps[[i]])){
      jpeg(filename = paste("BV Interactions Maps/BVInt_", gsub(" ", "", BBA_int_list[[i]]$Bird[1]), "_Int", j, ".jpg", sep = ""), 
           width = 20, height = 15, units = "cm", res = 300)
      print(maps[[i]][[j]])
      dev.off()
    }
    
  } else {
    next
  }
}











####################### 6. Animate interactions ###############################

# Load packages
library(gganimate)
library(transformr)

# Get world polygons
world <- ne_countries(scale = "large", returnclass = "sf", continent = "South America")

# Prep for parallel processing
numCores <- 3 #detectCores()
cl <- makeCluster(numCores)
clusterExport(cl = cl, c("BBA_Tracks", "AIS_Tracks", "WGS1984", "world", "wd"))
clusterEvalQ(cl, {
  
  # Load packages
  library(tidyverse)
  library(lubridate)
  library(move)
  library(sf)
  library(gganimate)
  library(transformr)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(rnaturalearthhires)
  library(pbapply)
  
  # Set WD
  setwd(wd)
  
})

# Get first and last DateTime of interactions, find nearest vessel DateTimes, plot together
pbsapply(BBA_int_list, cl = cl, function(x){
  
  # If there are interactions
  if(sum(!is.na(x$interaction)) > 0){
    
    anims1bird <- sapply(1:max(x$interaction, na.rm = TRUE), function(y){
      
      # Get bird track just for that interaction, with 2 mins buffer on either side (unless at end of track)
      min_index <- max(0, min(which(x$interaction == y)) - ceiling(120/x$resolution[1]))
      max_index <- min(nrow(x), max(which(x$interaction == y) + ceiling(120/x$resolution[1])))
      Btrack <- x[min_index:max_index, ]
      
      # Get first and last DateTime
      min_Bdt <- Btrack$DateTime[1]
      max_Bdt <- Btrack$DateTime[nrow(Btrack)]
      
      # Get vessel track associated with that interaction
      Vtrack <- AIS_Tracks %>% filter(Bird == Btrack$Bird[1], SubTraj == Btrack$SubTraj[1])
      
      # Get interpolated Vtrack points to smooth animation
      V_mv <- move(x = Vtrack$Longitude,
                   y = Vtrack$Latitude,
                   time = as.POSIXct(Vtrack$DateTime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"),
                   proj = WGS1984,
                   animal = as.factor(Vtrack$SubTraj),
                   data = Vtrack,
                   removeDuplicatedTimestamps = TRUE)
      V_mv_interp <- interpolateTime(V_mv, 
                                     time = Btrack$DateTime, 
                                     spaceMethod = "greatcircle")
      Vtrack2 <- as.data.frame(V_mv_interp) %>%
        dplyr::select(coords.x1, coords.x2, timestamps, Bird, 
                      #flag, ssvid, 
                      vessel_class, fishing, 
                      #vesselID, 
                      SubTraj) %>%
        rename(Longitude = coords.x1, Latitude = coords.x2, DateTime = timestamps) %>%
        mutate(DateTime = ymd_hms(DateTime))
      min_Vdt <- Vtrack2$DateTime[which.min(abs(Vtrack2$DateTime - min_Bdt))]
      max_Vdt <- Vtrack2$DateTime[which.min(abs(Vtrack2$DateTime - max_Bdt))]
      
      # Keep original Vtrack to show AIS fixes?
      Vtrack <- Vtrack %>% filter(DateTime >= min_Vdt, DateTime <= max_Vdt)
      
      # Get boundary coordinates
      minLong <- min(c(Btrack$Longitude, Vtrack$Longitude))
      maxLong <- max(c(Btrack$Longitude, Vtrack$Longitude))
      rangeLong <- abs(maxLong - minLong)
      minLat <- min(c(Btrack$Latitude, Vtrack$Latitude))
      maxLat <- max(c(Btrack$Latitude, Vtrack$Latitude))
      rangeLat <- abs(maxLat - minLat)
      
      # Plot if there is GLS leg data
      if(sum(!is.na(Btrack$leg_wet)) > 0){
        
        anim1 <- ggplot(data = world) + 
          geom_sf(color = "black", fill = "gray83") +
          geom_point(data = Vtrack2, aes(x = Longitude, y = Latitude, group = DateTime), size = 2) +
          #geom_point(data = Vtrack, aes(x = Longitude, y = Latitude, group = DateTime), size = 2, alpha = 0.5) + 
          geom_point(data = Btrack, aes(x = Longitude, y = Latitude, colour = leg_wet, group = DateTime), size = 3) +
          coord_sf(xlim = c(minLong, maxLong), ylim = c(minLat, maxLat), expand = TRUE, clip = "off") +
          labs(x = "Longitude", y = "Latitude", colour = "GLS - Leg", title = "DateTime (UTC): {frame_time}") +
          annotate(
            geom = "text", x = maxLong + 0.08*rangeLong, y = minLat - 0.12*rangeLat, 
            label = paste("Bird: ", Btrack$Bird[1], "\n",
                          "Resolution: ", Btrack$GPS_Schedule[1], "\n",
                          "Interaction Num: ", max(Btrack$interaction, na.rm = TRUE), "\n",
                          "SubTraj: ", Btrack$SubTraj[1], "\n", 
                          # "MMSI: ", Btrack$ssvid[1], "\n",
                          # "Flag: ", Btrack$flag[1], "\n", 
                          "Class: ", Btrack$vessel_class[1], "\n",
                          "Active fishing: ", ifelse(length(Btrack$fishing[Btrack$fishing == 1]) > 0, "Yes", "No"), "\n",
                          sep = ""), 
            hjust = 0, vjust = 0, size = 3) +
          theme_bw() +
          transition_time(DateTime) +
          ease_aes('linear') +
          shadow_trail(distance = 0.01, alpha = 0.05) #+
          #shadow_wake(wake_length = 1, alpha = 0.05, size = FALSE)
        
        #animate(anim1, nframes = 200, duration = 20, renderer = gifski_renderer())
        
        anim_save(filename = paste("BV Interactions Animations/BVint_", gsub(" ", "", x$Bird[1]), "_Int", y, ".gif", sep = ""), 
                  anim1, nframes = 200, duration = 20, renderer = gifski_renderer(),
                  width = 20, height = 15, units = "cm", res = 250)
        
        return(paste("Finished anim for ", x$Bird[1], " Int ", y, ".\n", sep = ""))
        
      # Plot if there is no GLS leg data  
      } else {
        
        anim1 <- ggplot(data = world) + 
          geom_sf(color = "black", fill = "gray83") +
          geom_point(data = Vtrack2, aes(x = Longitude, y = Latitude, group = DateTime), size = 2) +
          #geom_point(data = Vtrack, aes(x = Longitude, y = Latitude, group = DateTime), size = 2, alpha = 0.5) + 
          geom_point(data = Btrack, aes(x = Longitude, y = Latitude, colour = Bird, group = DateTime), size = 3) +
          coord_sf(xlim = c(minLong, maxLong), ylim = c(minLat, maxLat), expand = TRUE, clip = "off") +
          labs(x = "Longitude", y = "Latitude", title = "DateTime (UTC): {frame_time}") +
          annotate(
            geom = "text", x = maxLong + 0.08*rangeLong, y = minLat - 0.12*rangeLat, 
            label = paste("Bird: ", Btrack$Bird[1], "\n",
                          "Resolution: ", Btrack$GPS_Schedule[1], "\n",
                          "Interaction Num: ", max(Btrack$interaction, na.rm = TRUE), "\n",
                          "SubTraj: ", Btrack$SubTraj[1], "\n", 
                          # "MMSI: ", Btrack$ssvid[1], "\n",
                          # "Flag: ", Btrack$flag[1], "\n", 
                          "Class: ", Btrack$vessel_class[1], "\n",
                          "Active fishing: ", ifelse(length(Btrack$fishing[Btrack$fishing == 1]) > 0, "Yes", "No"), "\n",
                          sep = ""), 
            hjust = 0, vjust = 0, size = 3) +
          theme_bw() +
          transition_time(DateTime) +
          ease_aes('linear') +
          shadow_trail(distance = 0.01, alpha = 0.05) #+
        #shadow_wake(wake_length = 1, alpha = 0.05, size = FALSE)
        
        #animate(anim1, nframes = 200, duration = 20, renderer = gifski_renderer())
        
        anim_save(filename = paste("BV Interactions Animations/BVint_", gsub(" ", "", x$Bird[1]), "_Int", y, ".gif", sep = ""), 
                  anim1, nframes = 200, duration = 20, renderer = gifski_renderer(),
                  width = 20, height = 15, units = "cm", res = 250)
        
        return(paste("Finished anim for ", x$Bird[1], " Int ", y, ".\n", sep = ""))
      }
    })
    
    # Else if no interactions, return NULL  
  } else {
    return(NULL)
  }
})
stopCluster(cl)










####################### 7. Test and visualise data further ##################

# Inspect interaction summary
int_summary %>%
  group_by(bird_name) %>%
  summarise(duration = sum(int_dur), num_interations = n()) %>%
  arrange(desc(duration))

# Look at distance to vessels over time (function outputs)
ggplot(BBA_int_list[[1]], aes(DateTime, distance, colour = fishing)) +
  geom_point()
ggplot(BBA_int_list[[2]], aes(DateTime, distance, colour = fishing)) +
  geom_point()

# How many discrete interactions for each bird?
n_interactions <- sapply(BBA_int_list, function(x){
  max(x$interaction, na.rm = TRUE)
})
bird_names <- sapply(BBA_int_list, function(x){x$Bird[1]})
data.frame(bird_names, n_interactions)

# Duration of interactions in minutes for 75F red
dur_interactions <- sapply(1:6, function(x){
  vec75F <- BBA_int_list[[2]]$interaction[!is.na(BBA_int_list[[2]]$interaction) & BBA_int_list[[2]]$interaction == x]
  length(vec75F)
})
data.frame(int_index = 1:6, dur_interactions)







####################### 8. Explore maps further ########################

library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(sf)

#if you want lower res for quicker plotting, you can use scale = "small" or "medium"
world <- ne_countries(scale = "large", returnclass = "sf", continent = "South America")

# Show zoomed out view of all vessels
# This shows lots of vessels turn their AIS off for long periods of time
# Also, it looks like some vessels have multiple ssvids
vessels <- AIS_Tracks
ggplot(data = world) + 
  geom_sf(color = "black", fill = "gray83") +
  geom_path(data = vessels, aes(x = Longitude, y = Latitude, group = SubTraj), show.legend = FALSE) +
  geom_path(data = vessels, aes(x = Longitude, y = Latitude, colour = SubTraj), show.legend = FALSE) +
  coord_sf(xlim = c(min(vessels$Longitude), max(vessels$Longitude)), ylim = c(min(vessels$Latitude), max(vessels$Latitude)), expand = TRUE) +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw()

# 75F red: Show zoomed in view (loopy bit)
vessels <- AIS_Tracks %>% filter(DateTime > ymd_hm("2017-12-26 03:50"), DateTime < ymd_hm("2017-12-26 05:30"))
track <- BBA_Tracks %>% filter(Bird == "75F red", DateTime > ymd_hm("2017-12-26 03:50"), DateTime < ymd_hm("2017-12-26 05:30"))
ggplot(data = world) + 
  geom_sf(color = "black", fill = "gray83") +
  geom_path(data = vessels, aes(x = Longitude, y = Latitude, colour = vessel_class)) +
  geom_point(data = vessels, aes(x = Longitude, y = Latitude)) +
  geom_path(data = track, aes(x = Longitude, y = Latitude, colour = Bird)) +
  coord_sf(xlim = c(-64.22, -64.15), ylim = c(-53.86, -53.82), expand = FALSE) +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw()

# 75F red: plot areas with vessel less than 1km away
track <- BBA_int_list[[2]] %>% filter(distance < 200)
vessels <- AIS_Tracks %>% filter(Bird == "75F red")
ggplot(data = world) + 
  geom_sf(color = "black", fill = "gray83") +
  #geom_point(data = vessels, aes(x = Longitude, y = Latitude, colour = vessel_class)) +
  geom_path(data = vessels, aes(x = Longitude, y = Latitude, group = SubTraj)) +
  geom_path(data = track, aes(x = Longitude, y = Latitude, colour = SubTraj)) +
  coord_sf(xlim = c(min(track$Longitude), max(track$Longitude)), ylim = c(min(track$Latitude), max(track$Latitude)), expand = TRUE) +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw()


# Figure out middle coordinate of each interaction period
int_coords_list <- lapply(BBA_int_list, function(df){
  df1 <- df[!is.na(df$interaction), 
            c("Bird", "DateTime", "Latitude", "Longitude", "leg_wet", "interaction")]
  df3 <- data.frame(Bird = character(), interaction = double(), Latitude = double(), Longitude = double())
  if(length(df1$interaction) > 0){
    for(i in 1:max(df1$interaction)){
      df2 <- df1[df1$interaction == i, ]
      row <- ceiling(nrow(df2)/2)
      df3 <- add_row(df3, 
                     Bird = df2[row, "Bird"],
                     interaction = df2[row, "interaction"],
                     Latitude = df2[row, "Latitude"],
                     Longitude = df2[row, "Longitude"])
    }
  }
  return(df3)
})
int_coords <- do.call(rbind, int_coords_list)
int_summary <- read.csv("Visual Interactions Analysis/BBA Interactions Summary.csv") %>%
  mutate(Real_start = dmy_hm(Real_start), Real_end = dmy_hm(Real_end)) %>%
  dplyr::select(Bird, Int_Num, Int_Class) %>%
  rename(interaction = Int_Num)
int_coords <- left_join(int_coords, int_summary, by = join_by(Bird, interaction))
int_coords <- int_coords[int_coords$Int_Class != "Non-Interaction", ]

# Figure for manuscript: Plot all GPS tracks of albatrosses, with vessel interactions included
ggplot(data = world) + 
  geom_sf(color = "black", fill = "darkseagreen3") +
  geom_path(data = BBA_Tracks, aes(x = Longitude, y = Latitude, group = Bird), 
            show.legend = FALSE, linewidth = 0.2, alpha = 0.6) +
  geom_point(data = int_coords, aes(x = Longitude, y = Latitude, colour = Int_Class),
             size = 1.8) +
  coord_sf(xlim = c(min(BBA_Tracks$Longitude), max(BBA_Tracks$Longitude)), 
           ylim = c(min(BBA_Tracks$Latitude), max(BBA_Tracks$Latitude)), expand = TRUE) +
  labs(x = "Longitude", y = "Latitude", colour = "Type of\nInteraction") +
  theme_bw() +
  theme(text = element_text(size = 14))


# Figure for manuscript: Plot 75F red with vessel
world <- ne_countries(scale = "large", returnclass = "sf", continent = "South America")
vessels <- AIS_Tracks %>% filter(DateTime > ymd_hm("2017-12-26 03:50"), DateTime < ymd_hm("2017-12-26 05:30"))
track <- BBA_Tracks %>% filter(Bird == "75F red", DateTime > ymd_hm("2017-12-26 03:50"), DateTime < ymd_hm("2017-12-26 05:30"))
ggplot(data = world) + 
  geom_sf(color = "black", fill = "gray83") +
  geom_path(data = vessels, aes(x = Longitude, y = Latitude)) +
  geom_point(data = vessels, aes(x = Longitude, y = Latitude)) +
  geom_path(data = track, aes(x = Longitude, y = Latitude, colour = leg_wet, group = Bird)) +
  geom_point(data = track, aes(x = Longitude, y = Latitude, colour = leg_wet)) +
  coord_sf(xlim = c(-64.21, -64.16), ylim = c(-53.85, -53.825), expand = FALSE) +
  labs(x = "Longitude", y = "Latitude", colour = "Immersion\nState") +
  theme_bw() +
  theme(text = element_text(size = 14))

# # Figure for poster: 75F red with vessel, zoomed out a bit and high res
# png(file = "High res figures/Poster 75Fred.png", width = 2800, height = 2100, units = "px")
# ggplot(data = world) + 
#   geom_sf(color = "black", fill = "gray83") +
#   coord_sf(xlim = c(-64.22, -64.15), ylim = c(-53.855, -53.82), expand = FALSE) +
#   geom_path(data = vessels, aes(x = Longitude, y = Latitude), size = 3.5) +
#   geom_point(data = vessels, aes(x = Longitude, y = Latitude), size = 11.3) +
#   geom_path(data = track, aes(x = Longitude, y = Latitude, colour = leg_wet, group = Bird), size = 3.5) +
#   geom_point(data = track, aes(x = Longitude, y = Latitude, colour = leg_wet), size = 11.3) +
#   labs(x = "Longitude", y = "Latitude", colour = "Immersion\nState") +
#   theme_bw() +
#   theme(text = element_text(size = 40),
#         panel.grid = element_blank())
# dev.off()
