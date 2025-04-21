######################################################################
##################### GLS IMMERSION CORRECTION  ######################
######################################################################

############### GENERAL DESCRIPTION #######################

# Remove points on colony, split tracks with gaps over 10 mins
# Use 2-state Hidden Markov Models to predict wet/dry states for all high res tracks
# Merge these into immersion dataset as leg_wet_HMM column
# Compare immersion leg_wet to leg_wet_HMM and leg_wet_TDR
# Make a final corrected immersion column leg_wet_C based on rules


############### CONTENTS ##################################

# 1. Packages, WD, and polygons
# 2. Read in data
# 3. Initial inspection
# 4. Initial data prep
# 5. Inspect modified tracks
# 6. Prep for HMMs
# 7. HMM 5 sec: 2 vs 3 state, regularised tracks
# 8. HMM 10 sec: 2 vs 3 state, regularised tracks
# 9. Compare models and decode behavioural states
# 10. Add original immersion data
# 10a. Add TDR data
# 11. Compare immersion from GLS, TDR, and HMM 
# 12. Correct 52B brown (time shift GLS)
# 13. Functions to correct immersion as necessary
# 14. Run the correction functions
# 15. Final results to quantify GLS error



############## 1. PACKAGES, WD, AND POLYGONS ############################

# Load packages
library(tidyverse)
library(lubridate)
library(move)
library(momentuHMM)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(sf)
library(pbapply)
library(parallel)
theme_set(theme_bw())

# Set WD
wd <- "___________"
setwd(wd)

# Set CRS of input data: WGS 1984
WGS1984 <- st_crs(4326)[[2]]

# Get Falklands and South America polygons
falklands <- ne_countries(scale = "large", returnclass = "sf", country = "Falkland Islands")
s_america <- ne_countries(scale = "large", returnclass = "sf", continent = "South America")

# Buffer Falklands polygon by 1km
aeqd <- st_crs("+proj=aeqd +lat_0=-51.7 +lon_0=-61.3 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
falklands.proj <- st_transform(falklands, crs = aeqd)
falklands.buffer <- st_buffer(falklands.proj, dist = 1000)
falklands.buffer <- st_transform(falklands.buffer, crs = WGS1984)

# I don't want to hit <Return> to see next plot every time
par(ask = F)



################ 2. READ IN DATA ############################

# Read in GPS data
tracks0 <- readRDS("Output Dataset Files/BBA_Track_List_Imm_Int.RData")

# Narrow down the track
tracks <- pblapply(tracks0, function(x){
  
  # Narrow down to relevant columns
  x2 <- x[, c("Bird", "DateTime", "Latitude", "Longitude", "resolution")]
  
  # Find and cut out bits of track when bird is near colony, with 1km buffer
  track_sf <- st_as_sf(x2, coords = c("Longitude", "Latitude"), crs = WGS1984)
  new_track <- x2[which(!st_intersects(track_sf, falklands.buffer, sparse = FALSE, prepared = FALSE)), ]
  return(new_track)
})

# Read in raw immersion data
imm_data <- readRDS("Output Dataset Files/Wet_Data_Reg.Rdata")
birds_w_imm <- sapply(imm_data, function(x){x$Bird[1]})

# Read in visual interactions analysis data
int_summary <- read.csv("Visual Interactions Analysis/BBA Interactions Summary.csv") %>%
  mutate(Real_start = dmy_hm(Real_start), Real_end = dmy_hm(Real_end))

# Read in final datasets if necessary
# tracks.reg.list5 <- readRDS("Output Dataset Files/BBA_tracks_reg_list5.RData")
# tracks.comb3 <- do.call(rbind, tracks.reg.list5)
# percents <- readRDS("Output Dataset Files/imm_corr_percents.RData")




################## 3. INITIAL INSPECTION #####################################
# Note: This section is not necessary for rest of script to run

# # Visualise each track individually
# maps <- pbsapply(tracks, function(track){
#   
#   # Create map of track
#   map <- ggplot(data = s_america) + 
#     geom_sf(color = "black", fill = "gray83") +
#     geom_path(data = track, aes(x = Longitude, y = Latitude)) +
#     geom_point(data = track, aes(x = Longitude, y = Latitude), colour = "green", size = 0.2) +
#     coord_sf(xlim = c(min(track$Longitude), max(track$Longitude)), 
#              ylim = c(min(track$Latitude), max(track$Latitude)), 
#              expand = TRUE) +
#     labs(x = "Longitude", y = "Latitude",
#          title = paste0(track$Bird[1], ", resolution = ", track$resolution, " secs"))
#   
#   # Save map of track
#   jpeg(filename = paste("BBA Track Maps/Track_", gsub(" ", "", track$Bird[1]), ".jpg", sep = ""), 
#        width = 20, height = 15, units = "cm", res = 300, type = "cairo")
#   print(map)
#   dev.off()
#   
#   # Return confirmation
#   message(paste0(track$Bird[1], " done."))
#   return(paste0(track$Bird[1], " done."))
# })
# 
# # Find gaps of >60 sec for each track (dt = time interval)
# gaps.list <- pblapply(tracks, function(x){
#   x <- x %>% mutate(dt = as.numeric(difftime(DateTime, lag(DateTime)), units = "secs"))
#   x <- x %>% mutate(gap = ifelse(dt > 60 & resolution < 60, 1, 0))
#   gaps.df <- x[x$gap == 1, c("Bird", "DateTime", "dt")]
#   gaps.df <- gaps.df %>% filter(!is.na(dt))
#   if(nrow(gaps.df) == 0){
#     add_row(gaps.df, Bird = x$Bird[1], DateTime = NA, dt = NA)
#   }
#   return(gaps.df)
# })
# 
# # Dataframe of gap counts and maxes
# max_gaps <- data.frame(index = 1:length(gaps.list), 
#                        bird = sapply(gaps.list, function(x){x$Bird[1]}), 
#                        n.gaps = sapply(gaps.list, nrow),
#                        max.gap = sapply(gaps.list, function(x){max(x$dt)})) %>%
#   arrange(desc(max.gap)) %>%
#   filter(!is.na(bird))
# max_gaps # 12 birds have at least 1 gap over 1 min (not including low-res tracks)
# 
# # Dataframe of all gaps
# all_gaps <- do.call(rbind, gaps.list)
# sum(all_gaps$dt > 600) # 13 gaps total over 10 mins
# sum(all_gaps$dt > 300) # 19 gaps total over 5 mins
# sum(all_gaps$dt > 120) # 45 gaps total over 2 mins
# sum(all_gaps$dt > 60) # 94 gaps total over 1 mins
# 
# # For tracks with gaps, what is the proportion of points with dt>1min?
# data.frame(bird = max_gaps$bird, 
#            prop_low_res = sapply(1:nrow(max_gaps), function(i){
#              index <- which(sapply(tracks, function(y){y$Bird[1]}) == max_gaps$bird[i])
#              return(max_gaps$n.gaps[i] / nrow(tracks[[index]]))
#            }))




################ 4. INITIAL DATA PREP ##########################

# If this section already run, skip to PREP FOR HMMs

# Function to split tracks with missing data into sub-tracks
# max.gap = gap threshold in sec over which we will split the track
# max.res = resolution over which we will ignore gaps, as not doing HMMs on low res tracks
# remove.obs = Remove any sub-tracks with less than how many obs (< ~6 not useful for HMMs)
split.track <- function(track, max.gap, max.res, remove.obs){
  
  # Create dataframe of gaps, provided resolution is sufficiently high
  track <- track %>% mutate(dt = as.numeric(difftime(DateTime, lag(DateTime)), units = "secs"))
  track <- track %>% 
    mutate(gap = ifelse(dt > max.gap & resolution <= max.res, 1, 0)) %>%
    mutate(ID = Bird, gap = ifelse(is.na(gap), 0, gap)) %>%
    relocate(ID)
  
  # If this track has any gaps over max.gap, split the track
  if(sum(track$gap, na.rm = TRUE) > 0){
    
    # Number the sub-tracks
    track <- track %>% mutate(gap_num = cumsum(gap))
    
    # Initialise list of sub-tracks
    st_list <- vector(mode = "list", length = length(unique(track$gap_num)))
    
    # Split the track into sub-tracks
    for(i in unique(track$gap_num)){
      st_list[[i+1]] <- track[track$gap_num == i, ]
    }
    
    # Prep to label any sub-tracks with too few observations (<6 not useful for HMM)
    too_short <- which(sapply(st_list, nrow) < remove.obs)
    n_obs_removed <- 0
    
    # If any of the sub-tracks are too short
    if(length(too_short) > 0){
      normal <- which(sapply(st_list, nrow) >= remove.obs)
      n_obs_removed <- sum(sapply(st_list[too_short], nrow))
      
      # Label short sub-track IDs with "_noHMM"
      st_list[too_short] <- lapply(st_list[too_short], function(x){
        x$ID <- rep(paste0(x$Bird[1], "_noHMM"), times = nrow(x))
        return(x)
      })
      
      # Label the rest of them (the normal ones) with "_st#"
      st_list[normal] <- lapply(1:length(st_list[normal]), function(i){
        st_list[[normal[i]]]$ID <- rep(paste0(st_list[[normal[i]]]$Bird[1], "_st", i),
                                     times = nrow(st_list[[normal[i]]]))
        return(st_list[[normal[i]]])
      })
      
    # Else if sub-tracks are all reasonable length for HMM
    } else {
      
      # Label bird ID to indicate sub-track on remaining segments
      st_list <- lapply(1:length(st_list), function(i){
        st_list[[i]]$ID <- rep(paste0(st_list[[i]]$Bird[1], "_st", i), 
                               times = nrow(st_list[[i]]))
        return(st_list[[i]])
      })
    }
    
    # Remove some columns and return as a list
    st_list <- lapply(st_list, function(x){x[,-c(8,9)]})
    message(paste0(track$Bird[1], " done, ", length(st_list), " sub-tracks created, ", 
                   n_obs_removed, " observations marked for removal."))
    return(st_list)
    
  # Else return track within a single-element list
  } else {
    st_list <- list(track[ , -8])
    message(paste0(track$Bird[1], " done, no splits."))
    return(st_list)
  }
}


# Function to regularise bird tracks (interpolate or downsample)
# Returns regular track with just xyt columns
regularise.track <- function(track){
  
  # Get the resolution for that bird
  res <- track$resolution[1]
  
  # Convert to move object
  mv <- move(x = track$Longitude,
             y = track$Latitude,
             time = as.POSIXct(track$DateTime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"),
             proj = WGS1984,
             animal = track$Bird,
             data = track,
             removeDuplicatedTimestamps = TRUE)
  
  # Regularise track (DateTime intervals must be completely consistent)
  # If res == 1, then downsample to res == 5 to match scale of wet/dry states better
  if(res == 1){
    mv <- interpolateTime(mv, 
                          time = as.difftime(5, units = "secs"), 
                          spaceMethod = "greatcircle")
    
    # Else regularise track as normal
  } else {
    mv <- interpolateTime(mv, 
                          time = as.difftime(res, units = "secs"), 
                          spaceMethod = "greatcircle")
  }
  
  
  # Convert back to data.frame
  track.reg <- as.data.frame(mv)
  track.reg <- track.reg[, c("ID", "Bird", "timestamps", "coords.x2", "coords.x1", "resolution", "dt")]
  colnames(track.reg) <- c("ID", "Bird", "DateTime", "Latitude", "Longitude", "resolution", "dt")
  
  # If there are sub-tracks in this track
  if(track$ID[1] != track$Bird[1]){
    
    # If there are observations in short segments that should be removed
    IDs <- gsub(".*_", "", track$ID)
    if(sum(IDs == "noHMM") > 0){
      
      # Get start datetimes of non-HMM segments
      start <- ifelse(IDs == "noHMM" & (lag(IDs) != "noHMM" | is.na(lag(IDs))),
                    track$DateTime, NA)
      start_ind <- which(!is.na(start)) - 1
      if(start_ind[1] == 0){start_ind[1] <- 1}
      start_time <- as_datetime(track$DateTime[start_ind])
      
      # Get end datetimes of non-HMM segments
      end <- ifelse(IDs == "noHMM" & (lead(IDs) != "noHMM" | is.na(lead(IDs))),
                      track$DateTime, NA)
      end_ind <- which(!is.na(end)) + 1
      if(end_ind[length(end_ind)] > nrow(track)){end_ind[length(end_ind)] <- 1}
      end_time <- as_datetime(track$DateTime[end_ind])
      
      # Fill in ID column for non-HMM sub-tracks (data gaps)
      for(i in 1:length(start_time)){
        track.reg$ID <- ifelse(track.reg$DateTime > start_time[i] &
                                 track.reg$DateTime < end_time[i],
                               paste0(track.reg$Bird, "_noHMM"),
                               track.reg$ID)
      }
    }
  }
  
  # Fill in ID column for rest of sub-tracks
  track.reg <- track.reg %>%
    fill(ID, .direction = "down")
  
  # Return regularised track
  message(paste0(track.reg$Bird[1], " done."))
  return(track.reg)
}

# Split tracks if gaps are over 2 mins, then expand out the sub-tracks into the main list
tracks.split1 <- pblapply(tracks, split.track, max.gap = 120, max.res = 10, remove.obs = 6)

# Expand out the sub-tracks into the main list and combine sub-tracks for each bird
tracks.split1 <- unlist(tracks.split1, recursive = FALSE)
tracks.split.df <- do.call(rbind, tracks.split1)
tracks.split2 <- split(tracks.split.df, f = tracks.split.df$Bird)

# Regularise just res == 1 tracks (downsample to 5 second res)
# res1 <- which(sapply(tracks.split1, function(x){x$resolution[1]}) == 1)
# tracks.split2 <- tracks.split1
# tracks.split2[res1] <- pblapply(tracks.split2[res1], regularise.track)

# Prep for parallel processing
# numCores <- detectCores() - 4
# cl <- makeCluster(numCores)
# clusterEvalQ(cl, {
#   library(tidyverse)
#   library(lubridate)
#   library(move)
#   library(pbapply)
# })
# clusterExport(cl = cl, c("wd", "WGS1984", "tracks.split2"))

# Regularise all sub-tracks
#tracks.split2 <- tracks.split1
setwd(wd)
tracks.reg.split1 <- pblapply(tracks.split2, regularise.track) #, cl = cl)
#stopCluster()

# Save tracks just in case
#saveRDS(tracks.reg.split1, "Output Dataset Files/BBA_tracks_reg_split_list.RData")





################## 5. INSPECT MODIFIED TRACKS ###########################
# Note: This section is not necessary for the rest of the script to run

# Inspect how different the regular tracks are from the originals
data.frame(subtrack = sapply(tracks.reg.split1, function(x){x$ID[1]}), 
           res = sapply(tracks.reg.split1, function(x){x$resolution[1]}),
           confirm_res = sapply(tracks.reg.split1, 
                                function(x){as.numeric(difftime(x$DateTime[2], 
                                                                x$DateTime[1],
                                                                units = "secs"))}),
           diff = sapply(1:length(tracks.reg.split1), 
                         function(x){nrow(tracks.reg.split1[[x]]) - nrow(tracks.split1[[x]])}))
# Large differences in nrow: 00G blue, 14T blue, 33T blue, 38F red, 39T blue, 41A red, 
# 47F red, 59T blue, 60G blue, 60L yellow, 61T blue, 64 brown, 66F red, 
# 69G blue, 75L yellow, 78B brown, 81E white, 87K white, unringed
# Confirmed that 1 sec tracks have been downsampled to 5 sec.

# Inspect by comparing maps
# maps_reg <- pbsapply(tracks2, function(track){
#   
#   # Create map of track
#   map <- ggplot(data = s_america) + 
#     geom_sf(color = "black", fill = "gray83") +
#     geom_path(data = track, aes(x = Longitude, y = Latitude)) +
#     geom_point(data = track, aes(x = Longitude, y = Latitude), colour = "green", size = 0.2) +
#     coord_sf(xlim = c(min(track$Longitude), max(track$Longitude)), 
#              ylim = c(min(track$Latitude), max(track$Latitude)), 
#              expand = TRUE) +
#     labs(x = "Longitude", y = "Latitude",
#          title = paste0(track$ID[1], ", resolution = ", 
#                         abs(as.numeric(difftime(track$DateTime[1], track$DateTime[2], units = "secs"))), 
#                         " secs"))
#   
#   # Save map of track
#   jpeg(filename = paste("BBA Track Maps/Track_", gsub(" ", "", track$ID[1]), "_reg.jpg", sep = ""), 
#        width = 20, height = 15, units = "cm", res = 300, type = "cairo")
#   print(map)
#   dev.off()
#   
#   # Return confirmation
#   message(paste0(track$ID[1], " done."))
#   return(paste0(track$ID[1], " done."))
# })
# Most look okay
# 61T blue (60 secs) had a lot of points added
# 65T blue (10 secs) has a weird gap where bird appears too fast
# 69G blue (10 secs) has several big gaps that will hopefully count as dry flight
# 72E white (420 secs) has several big gaps
# 81E white (10 secs) has 2 massive gaps that will hopefully count as dry flight
# 85K white (420 secs) has less than 10 points, GPS is almost useless
# 92F red (60 secs) has several big gaps
# CONCLUSION: High res interpolated tracks should hopefully be ok for HMMs






#################### 6. PREP FOR HMMS ##############################################

# If necessary, load the split/regularised tracks
#tracks.reg.split1 <- readRDS("Output Dataset Files/BBA_tracks_reg_split_list.RData")

# Bind the dfs together
tracks.comb.all <- do.call(rbind, tracks.reg.split1)

# Remove the bits of track that are within data gaps
tracks.comb <- tracks.comb.all %>%
  mutate(ID_temp = gsub(".*_", "", ID)) %>%
  filter(ID_temp != "noHMM") %>%
  dplyr::select(!ID_temp)

# Inspect to make sure this worked
t1 <- as.data.frame(tracks.comb.all %>% group_by(Bird) %>% summarise(n = n()))
t2 <- as.data.frame(tracks.comb %>% group_by(Bird) %>% summarise(n = n()))
t3 <- data.frame(t2, comb.all = t1[,2]) %>% mutate(diff = n - comb.all)

# Add columns for leg_wet as defined by HMM
tracks.comb <- tracks.comb %>% 
  mutate(leg_wet_HMM2 = NA, leg_wet_HMM3.1 = NA, leg_wet_HMM3.2 = NA)

# Create vector of sub-tracks
birds <- unique(tracks.comb$ID)

# Which sub-tracks have 1 sec res?
birds1 <- unique(tracks.comb.all[tracks.comb.all$resolution == 1, "ID"])
tracks1 <- tracks.comb[tracks.comb$ID %in% birds1, ]

# Which sub-tracks have 5 sec res?
birds5 <- unique(tracks.comb.all[tracks.comb.all$resolution == 5, "ID"])
tracks5 <- tracks.comb[tracks.comb$ID %in% birds5, ]

# Which sub-tracks have 5 sec res for HMM (both 1 and 5 sec res)
birds5a <- c(birds1, birds5)
tracks5a <- rbind(tracks1, tracks5)

# Which sub-tracks have 10 sec res?
birds10 <- unique(tracks.comb.all[tracks.comb.all$resolution == 10, "ID"])
tracks10 <- tracks.comb[tracks.comb$ID %in% birds10, ]

# Which birds have low res (LR)?
birdsLR <- birds[!(birds %in% c(birds1, birds5, birds10))]

# Reference dataframe of which birds have which resolution
data.frame(bird = sapply(tracks, function(x){x$Bird[1]}),
           res = sapply(tracks, function(x){x$resolution[1]}))



# ################# HMM 1 sec #########################
# 
# # Create momentuHMMData object from data.frame
# BBA_data1 <- prepData(data = tracks1,
#                        type = "LL",
#                        coordNames = c("Longitude", "Latitude"))
# 
# # Label state names
# stateNames <- c("dry", "wet")
# 
# # Choose reasonable Par0 values for step length (gamma dist)
# # (Remember this has been downsampled to 5 sec res for the HMM)
# hist(BBA_data1$step, xlim = c(0, .2), breaks = 100)
# sum(BBA_data1$step == 0, na.rm = TRUE)/nrow(BBA_data1) # Proportion 0s = .28
# step_mean0 <- c(.06,.001)
# step_sd0 <- c(.06,.001)
# step_zm0 <- c(0,.28)
# step_Par0 <- c(step_mean0, step_sd0, step_zm0)
# 
# # Choose reasonable Par0 values for turning angle (wrapped Cauchy dist)
# hist(BBA_data1$angle)
# angle_Par0 <- c(.4,.9) # Concentration values
# 
# # Put initial params into list (from models above)
# Par0_hmm1 <- list(step = step_Par0,
#                    angle = angle_Par0)
# 
# # Distributions for observation processes
# dist = list(step = "gamma",
#             angle = "wrpcauchy")
# 
# # Fit model
# hmm1 <- fitHMM(data = BBA_data1,
#                 nbStates = 2,
#                 dist = dist,
#                 Par0 = Par0_hmm1,
#                 stateNames = stateNames)
# hmm1
# 
# # Add leg_wet_HMM to tracks1
# viterbi1 <- viterbi(hmm1)
# viterbi1 <- ifelse(viterbi1 == 1, "dry", "wet")
# tracks1$leg_wet_HMM <- viterbi1






################# 7. HMM 5 sec: 2 vs 3 state, regularised tracks #####################
# (This includes downsampled 1sec GPS data, having unsuccessfully tried 1sec in commented-out section above)

# Create momentuHMMData object from data.frame
BBA_data5 <- prepData(data = tracks5a,
                      type = "LL",
                      coordNames = c("Longitude", "Latitude"))

# Label state names
stateNames2 <- c("dry", "wet")
stateNames3 <- c("dryfast", "dryslow", "wet")

# Choose reasonable Par0 values for step length (gamma dist) - 2 state
hist(BBA_data5$step, xlim = c(0, .15), breaks = 200, xlab = "Step Length (km)", main = NULL)
sum(BBA_data5$step == 0, na.rm = TRUE)/nrow(BBA_data5) # Proportion 0s = .24
step_mean0 <- c(.06,.003)
step_sd0 <- c(.06,.003)
step_zm0 <- c(0,.24)
step_Par0_2 <- c(step_mean0, step_sd0, step_zm0)

# Choose reasonable Par0 values for step length (gamma dist) - 3 state
step_mean0 <- c(.06,.008,.001)
step_sd0 <- c(.06,.008,.001)
step_zm0 <- c(0,0,.24)
step_Par0_3 <- c(step_mean0, step_sd0, step_zm0)

# Choose reasonable Par0 values for turning angle (wrapped Cauchy dist) - 2 state
hist(BBA_data5$angle, xlab = "Turning Angle (degrees)", breaks = 200, main = NULL)
angle_Par0_2 <- c(.2,.9) # Concentration values

# Choose reasonable Par0 values for turning angle (wrapped Cauchy dist) - 3 state
angle_Par0_3 <- c(.2,.3,.9) # Concentration values

# Put initial params into list (from models above)
Par0_hmm5_2 <- list(step = step_Par0_2,
                     angle = angle_Par0_2)
Par0_hmm5_3 <- list(step = step_Par0_3,
                     angle = angle_Par0_3)

# Distributions for observation processes
dist = list(step = "gamma",
            angle = "wrpcauchy")

# Fit 2-state model
hmm5_2 <- fitHMM(data = BBA_data5,
                  nbStates = 2,
                  dist = dist,
                  Par0 = Par0_hmm5_2,
                  stateNames = stateNames2)

# Fit 3-state model
hmm5_3 <- fitHMM(data = BBA_data5,
                  nbStates = 3,
                  dist = dist,
                  Par0 = Par0_hmm5_3,
                  stateNames = stateNames3)




################# 8. HMM 10 sec: 2 vs 3 state, regularised tracks ######################
# 
# # Create momentuHMMData object from data.frame
# BBA_data10 <- prepData(data = tracks10,
#                        type = "LL",
#                        coordNames = c("Longitude", "Latitude"))
# 
# # Label state names
# stateNames2 <- c("dry", "wet")
# stateNames3 <- c("dryfast", "dryslow", "wet")
# 
# # Choose reasonable Par0 values for step length (gamma dist) - 2 state
# hist(BBA_data10$step, xlim = c(0, .3), breaks = 800, xlab = "Step Length (km)", main = NULL)
# sum(BBA_data10$step == 0, na.rm = TRUE)/nrow(BBA_data10) # Proportion 0s = .12
# step_mean0 <- c(.11,.005)
# step_sd0 <- c(.11,.005)
# step_zm0 <- c(0,.12)
# step_Par0_2 <- c(step_mean0, step_sd0, step_zm0)
# 
# # Choose reasonable Par0 values for step length (gamma dist) - 3 state
# step_mean0 <- c(.11,.02,.005)
# step_sd0 <- c(.11,.02,.005)
# step_zm0 <- c(0,0,.12)
# step_Par0_3 <- c(step_mean0, step_sd0, step_zm0)
# 
# # Choose reasonable Par0 values for turning angle (wrapped Cauchy dist) - 2 state
# hist(BBA_data10$angle)
# angle_Par0_2 <- c(.2,.9) # Concentration values
# 
# # Choose reasonable Par0 values for turning angle (wrapped Cauchy dist) - 3 state
# angle_Par0_3 <- c(.2,.3,.9) # Concentration values
# 
# # Put initial params into list (from models above)
# Par0_hmm10_2 <- list(step = step_Par0_2,
#                      angle = angle_Par0_2)
# Par0_hmm10_3 <- list(step = step_Par0_3,
#                      angle = angle_Par0_3)
# 
# # Distributions for observation processes
# dist = list(step = "gamma",
#             angle = "wrpcauchy")
# 
# # Fit 2-state model
# hmm10_2 <- fitHMM(data = BBA_data10,
#                 nbStates = 2,
#                 dist = dist,
#                 Par0 = Par0_hmm10_2,
#                 stateNames = stateNames2)
# 
# # Fit 3-state model
# hmm10_3 <- fitHMM(data = BBA_data10,
#                   nbStates = 3,
#                   dist = dist,
#                   Par0 = Par0_hmm10_3,
#                   stateNames = stateNames3)










################### 9. Compare models and decode behavioural states #################

# # If necessary, load models
# hmm5_2 <- readRDS("Output Dataset Files/HMMs/hmm5_2.RData")
# hmm5_3 <- readRDS("Output Dataset Files/HMMs/hmm5_3.RData")
# hmm10_2 <- readRDS("Output Dataset Files/HMMs/hmm10_2.RData")
# hmm10_3 <- readRDS("Output Dataset Files/HMMs/hmm10_3.RData")

# Save models
# saveRDS(hmm5_2, "Output Dataset Files/HMMs/hmm5_2.RData")
# saveRDS(hmm5_3, "Output Dataset Files/HMMs/hmm5_3.RData")
# saveRDS(hmm10_2, "Output Dataset Files/HMMs/hmm10_2.RData")
# saveRDS(hmm10_3, "Output Dataset Files/HMMs/hmm10_3.RData")

# Compare models
hmm5_2
hmm5_3
# hmm10_2
# hmm10_3
AIC(hmm5_2, hmm5_3)
#AIC(hmm10_2, hmm10_3)

# Look at plots
# par(ask = F)
# plot(hmm5_2)
# plot(hmm5_3)
# plot(hmm10_2)
# plot(hmm10_3)

# Decode behavioural states for all four models
viterbi5_2 <- viterbi(hmm5_2)
viterbi5_2 <- ifelse(viterbi5_2 == 1, "dry", "wet")
tracks5a$leg_wet_HMM2 <- viterbi5_2

viterbi5_3 <- viterbi(hmm5_3)
viterbi5_3.1 <- ifelse(viterbi5_3 < 3, "dry", "wet")
viterbi5_3.2 <- ifelse(viterbi5_3 < 2, "dry", "wet")
tracks5a$leg_wet_HMM3.1 <- viterbi5_3.1
tracks5a$leg_wet_HMM3.2 <- viterbi5_3.2

# viterbi10_2 <- viterbi(hmm10_2)
# viterbi10_2 <- ifelse(viterbi10_2 == 1, "dry", "wet")
# tracks10$leg_wet_HMM2 <- viterbi10_2
# 
# viterbi10_3 <- viterbi(hmm10_3)
# viterbi10_3.1 <- ifelse(viterbi10_3 < 3, "dry", "wet")
# viterbi10_3.2 <- ifelse(viterbi10_3 < 2, "dry", "wet")
# tracks10$leg_wet_HMM3.1 <- viterbi10_3.1
# tracks10$leg_wet_HMM3.2 <- viterbi10_3.2







############## 10. ADD ORIGINAL IMMERSION DATA ##################

# Bind the tracks together into 1 big df
# (This is smaller than tracks.comb.all because doesn't include noHMM gaps)
tracks.comb1 <- rbind(tracks5a, tracks10, 
                     tracks.comb[tracks.comb$ID %in% birdsLR, ])

# Join to larger df (same except with the data gaps filled in)
tracks.comb2 <- left_join(tracks.comb.all,
                          tracks.comb1[ , c("ID", "Bird", "DateTime", 
                                            "leg_wet_HMM2", "leg_wet_HMM3.1", "leg_wet_HMM3.2")],
                          by = c("ID", "Bird", "DateTime"))

# Add leg_wet and leg_wet_C column
tracks.comb2 <- tracks.comb2 %>% 
  mutate(leg_wet = NA, leg_wet_C = NA) %>%
  relocate(leg_wet, .before = leg_wet_HMM2)

# Make tracks into lists (split by Bird, not ID)
tracks.reg.list1 <- split(tracks.comb2, f = tracks.comb2$Bird)

# Dataframe of resolutions
resolution_df <- data.frame(Bird = sapply(tracks, function(x){x$Bird[1]}), 
                            resolution = sapply(tracks, function(x){x$resolution[1]}))

# Function to add the GLS immersion data to leg_wet
add.imm <- function(track){
  
  # Which bird is it?
  bird <- track$Bird[1]
  
  # What resolution is it?
  res <- track$resolution[1]
  
  # If there is immersion data for that bird...
  if(bird %in% birds_w_imm){
    
    # Get immersion dataframe
    imm_index <- which(birds_w_imm == bird)
    imm_df <- imm_data[[imm_index]]
    
    # Full join track df to imm_df
    df1 <- full_join(track[ , -which(colnames(track) == "leg_wet")], 
                     imm_df[ , -which(colnames(imm_df) == "Bird")], 
                     by = "DateTime")
    df1 <- df1 %>% 
      arrange(DateTime) %>%
      relocate(leg_wet, .before = leg_wet_HMM2)
    
    # Remove rows at beginning and end outside start/end of track
    # (Mainly for testing HMM leg_wet, and removing bits on land)
    if(min(which(!is.na(df1$resolution))) > 1){
      rows <- (min(which(!is.na(df1$resolution))) - 1):nrow(df1)
      df1 <- df1[rows, ]
    }
    if(max(which(!is.na(df1$resolution))) < nrow(df1)){
      rows <- 1:(max(which(!is.na(df1$resolution))) + 1)
      df1 <- df1[rows, ]
    }
    
    # Loop through rows to fill in leg_wet with nearest immersion value
    for(i in 1:nrow(df1)){
      
      # If there is already a wet/dry value, go to next row
      if(!is.na(df1$leg_wet[i])){
        next
        
      # Else if there is an <NA>
      } else {
        
        # Figure out previous DateTime with a wet/dry value
        for(j in 1:nrow(df1)){
          if(i - j > 0){
            if(!is.na(df1$leg_wet[i-j])){
              prevDT <- df1$DateTime[i-j]
              break
            } else {next}
          } else {break}
        }
        
        # Figure out next DateTime with a wet/dry value
        for(k in 1:nrow(df1)){
          if(i + k <= nrow(df1)){
            if(!is.na(df1$leg_wet[i+k])){
              nextDT <- df1$DateTime[i+k]
              break
            } else {next}
          } else {break}
        }
      }
      
      # Fill leg_wet with value from closest row in time
      if(as.numeric(abs(df1$DateTime[i] - prevDT)) <= 
         as.numeric(abs(df1$DateTime[i] - nextDT))){
        df1$leg_wet[i] <- df1$leg_wet[df1$DateTime == prevDT]
      } else {
        df1$leg_wet[i] <- df1$leg_wet[df1$DateTime == nextDT]
      }
    }
    
    # Remove rows originally from imm_df, return updated track
    originals <- which(df1$DateTime %in% track$DateTime)
    message(paste0(track$Bird[1], "done."))
    return(df1[originals, ])
    
    # Else return original track
  } else {
    message(paste0(track$Bird[1], "done."))
    return(track)
  }
}

# # Prep for parallel processing
# numCores <- detectCores() - 1
# cl <- makeCluster(numCores)
# clusterEvalQ(cl, {
#   library(tidyverse)
#   library(lubridate)
#   library(pbapply)
# })
# clusterExport(cl = cl, c("add.imm", "tracks.reg.list1", 
#                          "birds_w_imm", "imm_data", "resolution_df"))





# Run function to add immersion data
# tracks.reg.list4 includes mixed HMM leg_wet column
tracks.reg.list5 <- pblapply(tracks.reg.list1, add.imm)
#stopCluster(cl)
tracks.comb3 <- do.call(rbind, tracks.reg.list5)
saveRDS(tracks.reg.list5, "Output Dataset Files/BBA_tracks_reg_list5.RData")

# Load same data file here if needed
# tracks.reg.list5 <- readRDS("Output Dataset Files/BBA_tracks_reg_list5.RData")
# tracks.comb3 <- do.call(rbind, tracks.reg.list5)





################ 10a. ADD TDR IMMERSION DATA ########################
# TDR immersion error was so high we ignored it for main analysis

# # TDR immersion data already calculated, so add to list manually
# tracks.reg.list3 <- readRDS("Output Dataset Files/BBA_tracks_reg_list3.RData")
# tracks.reg.list4 <- lapply(tracks.reg.list4, function(x){
#   bird <- x$Bird[1]
#   index <- which(sapply(tracks.reg.list3, function(y){y$Bird[1]}) == bird)
#   tdr_vec <- tracks.reg.list3[[index]]$leg_wet_TDR
#   return(mutate(x, leg_wet_TDR = tdr_vec))
# })
# 
# # Save and remove tracks.reg.list3 object to clear space
# saveRDS(tracks.reg.list4, "Output Dataset Files/BBA_tracks_reg_list4.RData")
# rm(tracks.reg.list3)





# # Read in TDR immersion data
# tdr_data2 <- readRDS("Output Dataset Files/TDRimm_Data.RData")
# 
# # Dataframe of resolutions
# resolution_df <- data.frame(Bird = sapply(tracks, function(x){x$Bird[1]}), 
#                             resolution = sapply(tracks, function(x){x$resolution[1]}))
# 
# # Which birds have TDR immersion?
# birds_w_tdrimm <- sapply(tdr_data2, function(x){x$Bird[1]})
# 
# # Function to add TDR immersion data to leg_wet
# add.tdr.imm <- function(track){
#   
#   # Which bird is it?
#   bird <- track$Bird[1]
#   
#   # What resolution is it? Add to df
#   res <- resolution_df$resolution[resolution_df$Bird == bird]
#   track <- track %>% mutate(resolution = rep(res, nrow(track)))
#   
#   # Add TDR immersion column
#   track <- track %>%
#     mutate(leg_wet_TDR = NA) %>%
#     relocate(leg_wet_TDR, .before = leg_wet_C)
#   
#   # If there is TDR immersion for that bird...
#   if(bird %in% birds_w_tdrimm){
#     
#     # Get TDR dataframe
#     tdr_index <- which(birds_w_tdrimm == bird)
#     tdr_df <- tdr_data2[[tdr_index]]
#     
#     # Set all TDR values to dry to start
#     track$leg_wet_TDR <- "dry"
#     
#     # Loop through rows to fill in leg_wet_TDR with immersion values
#     # (Assuming here that TDR immersion is higher res than GPS and GLS immersion)
#     for(i in 1:nrow(track)){
#       
#       # If datetime is within wet zone, set leg_wet_TDR to wet
#       DT <- track$DateTime[i]
#       if(sum(DT %within% interval(tdr_df$Wet, tdr_df$Dry)) >= 1){
#         track$leg_wet_TDR[i] <- "wet"
#       }
#     }
#     
#     # Remove rows at beginning and end outside start/end of track
#     # (Mainly for testing HMM leg_wet, and removing bits on land)
#     if(min(which(!is.na(track$resolution))) > 1){
#       rows <- (min(which(!is.na(track$resolution))) - 1):nrow(track)
#       track <- track[rows, ]
#     }
#     if(max(which(!is.na(track$resolution))) < nrow(track)){
#       rows <- 1:(max(which(!is.na(track$resolution))) + 1)
#       track <- track[rows, ]
#     }
#     
#     # Return updated track
#     return(track)
#     
#     # Else return original track    
#   } else {
#     return(track)
#   }
# }
# 
# 
# # Prep for parallel processing
# numCores <- detectCores() - 2
# cl <- makeCluster(numCores)
# clusterEvalQ(cl, {
#   library(tidyverse)
#   library(lubridate)
#   library(pbapply)
# })
# clusterExport(cl = cl, c("add.tdr.imm", "tracks.reg.list2", 
#                          "birds_w_tdrimm", "tdr_data2", "resolution_df"))
# 
# # Run function to add immersion data
# tracks.reg.list3 <- pblapply(tracks.reg.list2, add.tdr.imm, cl = cl)
# stopCluster(cl)
# saveRDS(tracks.reg.list3, "Output Dataset Files/BBA_tracks_reg_list3.RData")







############# 11. COMPARE IMMERSION FROM GLS, TDR, and HMM #########

# Read in spreadsheet with coordinates of focus areas
# These were produced by looking at the track maps produced in section 3
# and specifying coordinates of areas that looked like they'd have both wet and dry immersion
coords_df <- read.csv("Visual Interactions Analysis/GLS vs HMM Immersion Inspection Coords.csv", 
                      header = TRUE)

# Vector of resolutions
resolution_vec <- sapply(tracks, function(x){x$resolution[1]})

# Vector of birds
birds <- sapply(tracks, function(x){x$Bird[1]})

# Function to visualise each track's GLS immersion
map.gls.immersion <- function(track){
  
  # Bird and resolution
  bird <- track$Bird[1]
  resolution <- resolution_vec[which(birds == bird)]
  
  # If track has GLS immersion
  if(!is.na(track$leg_wet[1])){
    
    # Create map of track
    map <- ggplot(data = s_america) + 
      geom_sf(color = "black", fill = "gray83") +
      geom_path(data = track, aes(x = Longitude, y = Latitude)) +
      geom_point(data = track, aes(x = Longitude, y = Latitude, colour = leg_wet), size = 0.4) +
      coord_sf(xlim = c(coords_df$xmin[coords_df$Bird == bird], 
                        coords_df$xmax[coords_df$Bird == bird]), 
               ylim = c(coords_df$ymin[coords_df$Bird == bird], 
                        coords_df$ymax[coords_df$Bird == bird]), 
               expand = TRUE) +
      labs(x = "Longitude", y = "Latitude",
           title = paste0(bird, ": GLS immersion (res ", resolution, ")"),
           colour = "State")
    
    # Save map of track
    jpeg(filename = paste("BBA Immersion Correction Maps/Immersion_", gsub(" ", "", track$Bird[1]), "_GLS.jpg", sep = ""), 
         width = 20, height = 15, units = "cm", res = 300, type = "cairo")
    print(map)
    dev.off()
  }
  
  # Return confirmation
  message(paste0(bird, " done."))
  return(paste0(bird, " done."))
}


# Function to visualise each track's HMM immersion
map.hmm.immersion <- function(track){
  
  # Bird and resolution
  bird <- track$Bird[1]
  resolution <- resolution_vec[which(birds == bird)]
  
  # If track has 2-state HMM immersion
  if(!is.na(track$leg_wet_HMM2[1])){
    
    # Create map of track
    map <- ggplot(data = s_america) + 
      geom_sf(color = "black", fill = "gray83") +
      geom_path(data = track, aes(x = Longitude, y = Latitude)) +
      geom_point(data = track, aes(x = Longitude, y = Latitude, colour = leg_wet_HMM2), size = 0.4) +
      coord_sf(xlim = c(coords_df$xmin[coords_df$Bird == bird], 
                        coords_df$xmax[coords_df$Bird == bird]), 
               ylim = c(coords_df$ymin[coords_df$Bird == bird], 
                        coords_df$ymax[coords_df$Bird == bird]), 
               expand = TRUE) +
      labs(x = "Longitude", y = "Latitude",
           title = paste0(bird, ": 2-state HMM immersion (res ", resolution, ")"),
           colour = "State")
    
    # Save map of track
    jpeg(filename = paste("BBA Immersion Correction Maps/Immersion_", gsub(" ", "", track$Bird[1]), "_HMM2.jpg", sep = ""), 
         width = 20, height = 15, units = "cm", res = 300, type = "cairo")
    print(map)
    dev.off()
  }
  
  # If track has 3-state HMM immersion
  if(!is.na(track$leg_wet_HMM3.1[1])){
    
    # Create map of track (2 dry states)
    map <- ggplot(data = s_america) + 
      geom_sf(color = "black", fill = "gray83") +
      geom_path(data = track, aes(x = Longitude, y = Latitude)) +
      geom_point(data = track, aes(x = Longitude, y = Latitude, colour = leg_wet_HMM3.1), size = 0.4) +
      coord_sf(xlim = c(coords_df$xmin[coords_df$Bird == bird], 
                        coords_df$xmax[coords_df$Bird == bird]), 
               ylim = c(coords_df$ymin[coords_df$Bird == bird], 
                        coords_df$ymax[coords_df$Bird == bird]), 
               expand = TRUE) +
      labs(x = "Longitude", y = "Latitude",
           title = paste0(bird, ": 3-state HMM immersion, 2 dry states (res ", resolution, ")"),
           colour = "State")
    
    # Save map of track
    jpeg(filename = paste("BBA Immersion Correction Maps/Immersion_", gsub(" ", "", track$Bird[1]), "_HMM3.jpg", sep = ""), 
         width = 20, height = 15, units = "cm", res = 300, type = "cairo")
    print(map)
    dev.off()
    
    # Create map of track (2 wet states)
    map <- ggplot(data = s_america) + 
      geom_sf(color = "black", fill = "gray83") +
      geom_path(data = track, aes(x = Longitude, y = Latitude)) +
      geom_point(data = track, aes(x = Longitude, y = Latitude, colour = leg_wet_HMM3.2), size = 0.4) +
      coord_sf(xlim = c(coords_df$xmin[coords_df$Bird == bird], 
                        coords_df$xmax[coords_df$Bird == bird]), 
               ylim = c(coords_df$ymin[coords_df$Bird == bird], 
                        coords_df$ymax[coords_df$Bird == bird]), 
               expand = TRUE) +
      labs(x = "Longitude", y = "Latitude",
           title = paste0(bird, ": 3-state HMM immersion, 2 wet states (res ", resolution, ")"),
           colour = "State")
    
    # Save map of track
    jpeg(filename = paste("BBA Immersion Correction Maps/Immersion_", gsub(" ", "", track$Bird[1]), "_HMM3.2.jpg", sep = ""), 
         width = 20, height = 15, units = "cm", res = 300, type = "cairo")
    print(map)
    dev.off()
  }
  
  # Return confirmation
  message(paste0(bird, " done."))
  return(paste0(bird, " done."))
}


# # Function to visualise each track's TDR immersion
# map.tdr.immersion <- function(track){
#   
#   # Bird and resolution
#   bird <- track$Bird[1]
#   resolution <- resolution_vec[which(birds == bird)]
#   
#   # If track has TDR immersion
#   if(!is.na(track$leg_wet_TDR[1])){
#     
#     # Create map of track
#     map <- ggplot(data = s_america) + 
#       geom_sf(color = "black", fill = "gray83") +
#       geom_path(data = track, aes(x = Longitude, y = Latitude)) +
#       geom_point(data = track, aes(x = Longitude, y = Latitude, colour = leg_wet_TDR), size = 0.4) +
#       coord_sf(xlim = c(coords_df$xmin[coords_df$Bird == bird], 
#                         coords_df$xmax[coords_df$Bird == bird]), 
#                ylim = c(coords_df$ymin[coords_df$Bird == bird], 
#                         coords_df$ymax[coords_df$Bird == bird]), 
#                expand = TRUE) +
#       labs(x = "Longitude", y = "Latitude",
#            title = paste0(bird, ": TDR immersion (res ", resolution, ")"),
#            colour = "State")
#     
#     # Save map of track
#     jpeg(filename = paste("BBA Immersion Correction Maps/Immersion_", gsub(" ", "", track$Bird[1]), "_TDR.jpg", sep = ""), 
#          width = 20, height = 15, units = "cm", res = 300, type = "cairo")
#     print(map)
#     dev.off()
#   }
#   
#   # Return confirmation
#   message(paste0(bird, " done."))
#   return(paste0(bird, " done."))
# }

# # Function to visualise each track's MHMM1 immersion (continuous random effects included)
# map.hmm.mixed.immersion <- function(track){
#   
#   # Bird and resolution
#   bird <- track$Bird[1]
#   resolution <- resolution_vec[which(birds == bird)]
#   
#   # If track has HMM immersion
#   if(!is.na(track$leg_wet_MHMM1[1])){
#     
#     # Create map of track
#     map <- ggplot(data = s_america) + 
#       geom_sf(color = "black", fill = "gray83") +
#       geom_path(data = track, aes(x = Longitude, y = Latitude)) +
#       geom_point(data = track, aes(x = Longitude, y = Latitude, colour = leg_wet_HMM), size = 0.4) +
#       coord_sf(xlim = c(coords_df$xmin[coords_df$Bird == bird], 
#                         coords_df$xmax[coords_df$Bird == bird]), 
#                ylim = c(coords_df$ymin[coords_df$Bird == bird], 
#                         coords_df$ymax[coords_df$Bird == bird]), 
#                expand = TRUE) +
#       labs(x = "Longitude", y = "Latitude",
#            title = paste0(bird, ": MHMM immersion - cont random effect (res ", resolution, ")"),
#            colour = "State")
#     
#     # Save map of track
#     jpeg(filename = paste("BBA Immersion Correction Maps 1-9-2023/Immersion_", gsub(" ", "", track$Bird[1]), "_MHMM1.jpg", sep = ""), 
#          width = 20, height = 15, units = "cm", res = 300, type = "cairo")
#     print(map)
#     dev.off()
#   }
#   
#   # Return confirmation
#   message(paste0(bird, " done."))
#   return(paste0(bird, " done."))
# }



# Run functions to visualise GLS vs HMM tracks
pbsapply(tracks.reg.list5, map.gls.immersion)
pbsapply(tracks.reg.list5, map.hmm.immersion)
#pbsapply(tracks.reg.list3, map.tdr.immersion)
#pbsapply(tracks.reg.list4, map.hmm.mixed.immersion)

# # Compare HMM with MHMM1 results: tracks5a
# sum(tracks5a$leg_wet_HMM != tracks5a$leg_wet_MHMM1, na.rm = TRUE)/nrow(tracks5a) #.01% mismatch
# sum(tracks5a$leg_wet_HMM != tracks5a$leg_wet_MHMM1, na.rm = TRUE) # 115 time steps mismatch
# sum(tracks5a[tracks5a$Bird == "74F red", "leg_wet_HMM"] != tracks5a[tracks5a$Bird == "74F red", "leg_wet_MHMM1"], na.rm = TRUE) # No mismatch 74F red
# 
# # Compare HMM with MHMM1 results: tracks10
# sum(tracks10$leg_wet_HMM != tracks10$leg_wet_MHMM1, na.rm = TRUE)/nrow(tracks10) # 2.2% mismatch
# sum(tracks10$leg_wet_HMM != tracks10$leg_wet_MHMM1, na.rm = TRUE) # 3612 time steps mismatch

# Compare HMM2 with HMM3 results: 5-sec res
tracks5a.1 <- tracks.comb3 %>% filter(resolution == 5)
sum(tracks5a.1$leg_wet_HMM2 != tracks5a.1$leg_wet_HMM3.1, na.rm = TRUE)/
  nrow(tracks5a.1) #13.0% mismatch between wet/dry and wet/dry/dry
sum(tracks5a.1$leg_wet_HMM2 != tracks5a.1$leg_wet_HMM3.2, na.rm = TRUE)/
  nrow(tracks5a.1) #0.3% mismatch between wet/dry and wet/wet/dry

# Compare HMM2 with HMM3 results: 10-sec res
# tracks10.1 <- tracks.comb3 %>% filter(resolution == 10)
# sum(tracks10.1$leg_wet_HMM2 != tracks10.1$leg_wet_HMM3.1, na.rm = TRUE)/
#   nrow(tracks10.1) #1.8% mismatch between wet/dry and wet/dry/dry
# sum(tracks10.1$leg_wet_HMM2 != tracks10.1$leg_wet_HMM3.2, na.rm = TRUE)/
#   nrow(tracks10.1) #23.8% mismatch between wet/dry and wet/wet/dry

# # Look more closely at problematic spot on 74F red
# red74F <- tracks.comb3 %>% filter(Bird == "74F red")
# ggplot(data = s_america) + 
#   geom_sf(color = "black", fill = "gray83") +
#   geom_path(data = red74F, aes(x = Longitude, y = Latitude)) +
#   geom_point(data = red74F, aes(x = Longitude, y = Latitude, colour = leg_wet), size = 0.5) +
#   coord_sf(xlim = c(-64.52, -64.45), 
#            ylim = c(-53.9, -53.85), 
#            expand = TRUE) +
#   labs(x = "Longitude", y = "Latitude",
#        title = "74F red GLS",
#        colour = "State")
# ggplot(data = s_america) + 
#   geom_sf(color = "black", fill = "gray83") +
#   geom_path(data = red74F, aes(x = Longitude, y = Latitude)) +
#   geom_point(data = red74F, aes(x = Longitude, y = Latitude, colour = leg_wet_HMM3.2), size = 0.5) +
#   coord_sf(xlim = c(-64.52, -64.45), 
#            ylim = c(-53.9, -53.85), 
#            expand = TRUE) +
#   labs(x = "Longitude", y = "Latitude",
#        title = "74F red HMM",
#        colour = "State")


# Determine % of leg_wet not matched by leg_wet_HMM for each bird
percents_list <- pblapply(tracks.reg.list5, function(track){
  
  # Remove bits of track with noHMM (they were interpolated within large gaps)
  track <- track %>% filter(!grepl("_noHMM", ID))
  
  # Coerce leg_wet to character
  track$leg_wet <- as.character(track$leg_wet)
  
  # If the bird has immersion values from both GLS and HMM2
  if(!is.na(track$leg_wet[1]) & !is.na(track$leg_wet_HMM2[1])){
    
    # Add new columns for immersion check
    track <- track %>% mutate(imm_check2 = NA, imm_check3 = NA)
    
    # Fill imm_check columns for 2-state HMMs
    track$imm_check2 <- ifelse(track$leg_wet == track$leg_wet_HMM2, "match",
                              ifelse(track$leg_wet == "dry" & track$leg_wet_HMM2 == "wet", "GLS dry HMM wet", 
                                     ifelse(track$leg_wet == "wet" & track$leg_wet_HMM2 == "dry", "GLS wet HMM dry", NA)))
    
    # If track is res 10, 3-state HMM should be wet/dry/dry (HMM3.1).
    if(track$resolution[1] == 10){
      track$imm_check3 <- ifelse(track$leg_wet == track$leg_wet_HMM3.1, "match",
                                 ifelse(track$leg_wet == "dry" & track$leg_wet_HMM3.1 == "wet", "GLS dry HMM wet", 
                                        ifelse(track$leg_wet == "wet" & track$leg_wet_HMM3.1 == "dry", "GLS wet HMM dry", NA)))
    
    # Else if track is res 1 or 5, 3-state HMM should be wet/wet/dry (HMM3.2)
    } else {
      track$imm_check3 <- ifelse(track$leg_wet == track$leg_wet_HMM3.2, "match",
                                 ifelse(track$leg_wet == "dry" & track$leg_wet_HMM3.2 == "wet", "GLS dry HMM wet", 
                                        ifelse(track$leg_wet == "wet" & track$leg_wet_HMM3.2 == "dry", "GLS wet HMM dry", NA)))
    }
    
    # Calculate percents of each imm_check category
    percents <- data.frame(Bird = track$Bird[1],
                           resolution = track$resolution[1],
                           match2 = sum(track$imm_check2 == "match")/nrow(track),
                           GLSdry_HMM2wet = sum(track$imm_check2 == "GLS dry HMM wet")/sum(track$leg_wet == "dry"),
                           GLSwet_HMM2dry = sum(track$imm_check2 == "GLS wet HMM dry")/sum(track$leg_wet == "wet"),
                           match3 = sum(track$imm_check3 == "match")/nrow(track),
                           GLSdry_HMM3wet = sum(track$imm_check3 == "GLS dry HMM wet")/sum(track$leg_wet == "dry"),
                           GLSwet_HMM3dry = sum(track$imm_check3 == "GLS wet HMM dry")/sum(track$leg_wet == "wet"))
    percents <- percents %>% 
      mutate(incorrect2 = 1 - match2, incorrect3 = 1 - match3)
    
    # Return percents df
    return(percents)
    
    # Else return dataframe empty
  } else {
    
    percents <- data.frame(Bird = track$Bird[1], resolution = track$resolution[1], 
                           match2 = NA, GLSdry_HMM2wet = NA, GLSwet_HMM2dry = NA, 
                           match3 = NA, GLSdry_HMM3wet = NA, GLSwet_HMM3dry = NA, 
                           incorrect2 = NA, incorrect3 = NA)
    return(percents)
  }
})
percents <- do.call(rbind, percents_list) %>%
  filter(resolution < 60)

# If resolutions wrong, correct these
percents$resolution <- sapply(tracks, function(x){x$resolution[1]})
percents

# Save percents df
saveRDS(percents, file = "Output Dataset Files/imm_corr_percents.Rdata")

# # If necessary, read in same df
# percents <- readRDS("Output Dataset Files/imm_corr_percents.Rdata")





# Plot discrepancy between GLS, HMM2, and HMM3 immersion

# 2- vs 3- state HMM: Which is closer to GLS?
ggplot(percents, aes(x = Bird, y = match2 - match3, fill = as.factor(resolution))) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major.x = element_blank()) +
  labs(y = "Difference between 2-state and 3-state HMM\nproportion of matches with GLS",
       fill = "Fix interval (s)")

# GLS dry HMM wet: may suggest slow flight speed?
ggplot(percents, aes(x = Bird, y = GLSdry_HMM2wet, fill = as.factor(resolution))) + # HMM2 only
  geom_col() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major.x = element_blank())
ggplot(percents, aes(x = Bird, y = GLSdry_HMM3wet, fill = as.factor(resolution))) + # HMM3 only
  geom_col() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major.x = element_blank()) +
  labs(y = "Proportion of GLS dry timesteps with HMM wet",
       fill = "Fix interval (s)")
ggplot(percents, aes(x = Bird, 
                     y = ifelse(resolution < 10, GLSdry_HMM2wet, GLSdry_HMM3wet), # Uses correct HMM
                     fill = as.factor(resolution))) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major.x = element_blank()) +
  labs(y = "Proportion of GLS dry timesteps\nwith HMM wet",
       fill = "Fix interval (s)")
# With HMM3, error increases for green (res5), decreases for blue (res10)
# Conclusion: HMM2 better for res 5, HMM3 better for res 10
# 53T blue GLS (res 5) is completely messed up
# Apart from 53T blue, error may be low enough to suggest GLS dry is correct

# GLS wet HMM dry: may suggest drop of water on device while bird flying
ggplot(percents, aes(x = Bird, y = GLSwet_HMM2dry, fill = as.factor(resolution))) + # HMM2 only
  geom_col() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major.x = element_blank())
ggplot(percents, aes(x = Bird, y = GLSwet_HMM3dry, fill = as.factor(resolution))) + # HMM3 only
  geom_col() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major.x = element_blank()) +
  labs(y = "Proportion of GLS wet timesteps with HMM dry",
       fill = "Fix interval (s)")
ggplot(percents, aes(x = Bird, 
                     y = ifelse(resolution < 10, GLSwet_HMM2dry, GLSwet_HMM3dry), # Uses correct HMM
                     fill = as.factor(resolution))) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major.x = element_blank()) +
  labs(y = "Proportion of GLS wet timesteps\nwith HMM dry",
       fill = "Fix interval (s)")
# Much higher errors, likely because GLS stays wet for too long
# With HMM3, error increases for blue (res10), decreases for green (res5)
# With HMM3, res10 tracks have slightly less wet, res5 have slightly more wet than HMM2
# So, if GLS has too much wet, a more correct res10 HMM3 with less wet
# could result in higher GLSwetHMMdry error.


# Incorrect total
ggplot(percents, aes(x = Bird, y = incorrect2, fill = as.factor(resolution))) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major.x = element_blank())
ggplot(percents, aes(x = Bird, y = incorrect3, fill = as.factor(resolution))) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major.x = element_blank())

# Proportion of birds that may need correction
sum((percents$GLSwet_HMM2dry > 0.05 & percents$resolution == 5) | 
      (percents$GLSwet_HMM3dry > 0.05 & percents$resolution == 10),
    na.rm= TRUE)/sum(!is.na(percents$GLSwet_HMM2dry))
# 60.5% of birds (23) with GLS and high res have >5% GLSwet with HMM2dry
sum((percents$GLSwet_HMM2dry > 0.02 & percents$resolution == 5) | 
      (percents$GLSwet_HMM3dry > 0.02 & percents$resolution == 10),
    na.rm= TRUE)/sum(!is.na(percents$GLSwet_HMM2dry))
# 71.1% of birds (27) with GLS and high res have >2% GLSwet_HMM2dry
sum((percents$GLSdry_HMM2wet > 0.05 & percents$resolution == 5) | 
      (percents$GLSdry_HMM3wet > 0.05 & percents$resolution == 10) , na.rm= TRUE)
# Only 1 bird (53T blue) has GLS dry HMM wet error over 5%
sum((percents$GLSdry_HMM2wet > 0.02 & percents$resolution == 5) | 
      (percents$GLSdry_HMM3wet > 0.02 & percents$resolution == 10) , na.rm= TRUE)
# 9 birds have GLS dry HMM wet error over 2%


# Need to inspect some of these in ArcGIS to fully understand errors

# # Inspect 14T Blue closer (points on nest may not have been cut out sufficiently)
# blue49T <- tracks.comb3 %>% filter(Bird == "49T blue")
# blue49T_DW <- blue49T %>% filter(leg_wet == "dry" & leg_wet_HMM2 == "wet")
# nrow(blue49T_DW)/nrow(blue49T)
# ggplot(data = falklands.buffer) + 
#   geom_sf(color = "black", fill = "gray83") +
#   geom_path(data = blue49T, aes(x = Longitude, y = Latitude)) +
#   geom_point(data = blue49T_DW, aes(x = Longitude, y = Latitude, colour = "red"), size = 0.4) +
#   coord_sf(xlim = c(min(blue49T$Longitude),
#                     max(blue49T$Longitude)),
#            ylim = c(min(blue49T$Latitude),
#                     max(blue49T$Latitude)),
#            expand = TRUE) +
#   # coord_sf(xlim = c(-61.4, -61.2), ylim = c(-51.8, -51.6)) +
#   labs(x = "Longitude", y = "Latitude",
#        title = "49T Blue GLS dry HMM wet")
# ggplot(blue49T_DW, aes(x = DateTime, y = !is.na(leg_wet))) + geom_point()

# # Inspect 47F red closer (same reason as 14T blue)
# red47F <- tracks.reg.list2[[which(sapply(tracks.reg.list2, function(x){x$Bird[1]}) == "47F red")]]
# red47F_DW <- red47F %>% filter(leg_wet == "dry" & leg_wet_HMM == "wet")
# nrow(red47F_DW)/nrow(red47F)
# ggplot(data = falklands.buffer) + 
#   geom_sf(color = "black", fill = "gray83") +
#   geom_path(data = red47F, aes(x = Longitude, y = Latitude)) +
#   geom_point(data = red47F_DW, aes(x = Longitude, y = Latitude, colour = "red"), size = 0.4) +
#   coord_sf(xlim = c(-61.4, -61.2), ylim = c(-51.8, -51.6)) +
#   labs(x = "Longitude", y = "Latitude",
#        title = "47F Red GLS dry HMM wet")

# # Inspect whether multiple tracks have incorrect HMMs because of gaps (> 1 min) or nest
# tracks.reg5 <- do.call(rbind, tracks.reg.list5) #do.call(rbind, tracks.reg.list2)
# gaps.list <- pblapply(tracks, function(x){
#   x <- x %>% mutate(dt = as.numeric(difftime(DateTime, lag(DateTime)), units = "secs"))
#   x <- x %>% mutate(gap = ifelse(dt > 60 & resolution < 60, 1, 0))
#   gaps.df <- x[x$gap == 1, c("Bird", "DateTime", "dt")]
#   gaps.df <- gaps.df %>% filter(!is.na(dt))
#   if(nrow(gaps.df) == 0){
#     add_row(gaps.df, Bird = x$Bird[1], DateTime = NA, dt = NA)
#   }
#   return(gaps.df)
# })
# # Dataframe of gap counts and maxes
# max_gaps <- data.frame(index = 1:length(gaps.list), 
#                        bird = sapply(gaps.list, function(x){x$Bird[1]}), 
#                        count = sapply(gaps.list, nrow),
#                        max.gap = sapply(gaps.list, function(x){max(x$dt, na.rm = TRUE)})) %>%
#   arrange(desc(max.gap))
# 
# # Dataframe of all gaps
# all_gaps <- do.call(rbind, gaps.list)
# sum(all_gaps$dt > 600)
# 
# # Figure out which GPS positions are within gaps
# tracks.reg.list5.gaps <- lapply(tracks.reg.list5, function(x){
#   bird <- x$Bird[1]
#   
#   if(!(bird %in% sapply(gaps.list, function(y){y$Bird[1]}))){
#     return(mutate(x, gap = 0))
#     
#   } else {
#     gaps.df <- gaps.list[[which(sapply(gaps.list, function(y){y$Bird[1]}) == bird)]]
#     gaps.df <- gaps.df %>% 
#       filter(DateTime > min(x$DateTime), DateTime < max(x$DateTime)) %>%
#       mutate(Latitude = NA, Longitude = NA, resolution = NA, leg_wet = NA, 
#              leg_wet_HMM2 = NA, leg_wet_HMM3.1 = NA, leg_wet_C = NA,
#              leg_wet_HMM3.2 = NA, ID = NA) %>%
#       relocate(dt, .after = resolution) %>%
#       relocate(ID, .before = Bird)
#     x <- mutate(x, dt = NA)
#     x.gaps <- rbind(x, gaps.df)
#     x.gaps <- x.gaps %>% arrange(DateTime) %>% mutate(gap = NA)
#     
#     # Classify gaps
#     for(row in 1:nrow(x.gaps)){
#       dt1 <- x.gaps$dt[row]
#       if(!is.na(dt1)){
#         start <- x.gaps$DateTime[row] - seconds(dt1)
#         startrow <- min(which(x.gaps$DateTime >= start))
#         x.gaps$gap[startrow:row] <- 1
#       } else {
#         x.gaps$gap[row] <- 0
#       }
#     }
#     
#     # Remove extra rows and columns
#     x.gaps <- x.gaps[is.na(x.gaps$dt),]
#     #x.gaps <- x.gaps[, -9]
#     return(x.gaps)
#   }
#   
# })
# 
# tracks.comb3.gaps <- do.call(rbind, tracks.reg.list5.gaps)

# # Inspect birds with biggest GLS dry HMM wet error
# ggplot(tracks.comb3.gaps %>% filter(Bird %in% c("14T blue", "47F red", "52B brown", "53T blue", "69G blue", "81E white")), 
#        aes(x = DateTime, 
#            y = ifelse((leg_wet == "dry" & leg_wet_HMM2 == "wet"), 1, NA))) + 
#   facet_wrap(~Bird, scales = "free", ncol = 1) +
#   geom_vline(xintercept = tracks.comb3.gaps$DateTime[tracks.comb3.gaps$gap == 1], col = "lightgreen", alpha = 0.1) +
#   geom_point() +
#   labs(title = "GLS dry HMM wet")
# 
# # Inspect other birds with smaller error
# ggplot(tracks.comb3.gaps %>% filter(Bird %in% c("65T blue", "33T blue", "00G blue", "22V red", "66F red", "74F red")), 
#        aes(x = DateTime, 
#            y = ifelse((leg_wet == "dry" & leg_wet_HMM2 == "wet"), 1, NA))) + 
#   facet_wrap(~Bird, scales = "free", ncol = 1) +
#   geom_vline(xintercept = tracks.comb3.gaps$DateTime[tracks.comb3.gaps$gap == 1], col = "lightgreen", alpha = 0.1) +
#   geom_point() +
#   labs(title = "GLS dry HMM wet")






# ################ 12. CORRECT 52B BROWN (TIME SHIFT GLS) ##################
# 
# # Check if GLS timing was shifted, to explain high GLS dry HMM wet birds
# pbsapply(tracks.reg.list5, function(track){
#   
#   if(!is.na(track$leg_wet_HMM2[1])){
#     # Make time series
#     timeseries <- ggplot(track, 
#                          aes(x = DateTime, 
#                              y = ifelse((leg_wet_HMM2 == "wet"), 1, NA))) + 
#       geom_vline(xintercept = track$DateTime[track$leg_wet == "wet"], 
#                  col = "lightblue", alpha = 0.1) +
#       geom_point(size = 1) +
#       labs(title = paste0(track$Bird[1], ": GLS wet (blue) vs HMM wet (black dots)"), 
#            x = "DateTime", y = "")
#     
#     # Save map of track
#     jpeg(filename = paste("BBA Immersion Correction Time Series/Immersion_", gsub(" ", "", track$Bird[1]), "_TDR.jpg", sep = ""), 
#          width = 20, height = 8, units = "cm", res = 300, type = "cairo")
#     print(timeseries)
#     dev.off()
#   }
# })
# # 52B brown: Looks like GLS was ~1 hr behind
# # 53T blue: Timing matches up, but gets dodgy halfway through track
# 
# # Isolate 52B brown track
# brown52B <- tracks.comb3 %>% 
#   filter(Bird == "52B brown") %>%
#   mutate(leg_wet2 = NA)
# 
# # Create df of all possible shift times between 35 and 60 mins
# # (Convert to seconds, then divide by 5 sec fix rate)
# 35*60/5; 60*60/5
# shift_df <- data.frame(shift = 600:800,
#                        GLSdryHMMwet = NA,
#                        GLSwetHMMdry = NA)
# 
# # Loop through those shift times to fill df
# for(shift in 600:800){
#   shifted <- as.character(brown52B$leg_wet[1:(nrow(brown52B) - shift)])
#   shifted <- c(rep("dry", times = shift), shifted)
#   brown52B$leg_wet2 <- shifted
#   GLSdryHMMwet <- sum(ifelse(brown52B$leg_wet2 == "dry" & brown52B$leg_wet_HMM2 == "wet",1,0))/
#     nrow(brown52B)
#   GLSwetHMMdry <- sum(ifelse(brown52B$leg_wet2 == "wet" & brown52B$leg_wet_HMM2 == "dry",1,0))/
#     nrow(brown52B)
#   shift_df[shift_df$shift == shift, "GLSdryHMMwet"] <- GLSdryHMMwet
#   shift_df[shift_df$shift == shift, "GLSwetHMMdry"] <- GLSwetHMMdry
# }
# 
# # Plot errors depending on shift
# ggplot(shift_df, aes(x = shift)) +
#   geom_line(aes(y = GLSdryHMMwet), col = "darkred") +
#   geom_line(aes(y = GLSwetHMMdry), col = "darkblue")
# 
# # Lowest error suggests ideal backwards shift: 718-720 rows or ~ 1 hr
# 610*5/60; 660*5/60
# shift_df$shift[which.min(shift_df$GLSdryHMMwet)]*5/60
# 
# # Adjust shift around 650 rows to line up with leg_wet_HMM2
# shift <- 720
# shifted <- as.character(brown52B$leg_wet[1:(nrow(brown52B) - shift)])
# shifted <- c(rep("dry", times = shift), shifted)
# brown52B$leg_wet2 <- shifted
# min(which(brown52B$leg_wet_HMM2 == "wet"))
# brown52B[1030:1050, c("DateTime", "leg_wet", "leg_wet2", "leg_wet_HMM2")]
# 
# # Correct Brown 52B and put back into main list and dataframe
# brown52B <- brown52B %>%
#   mutate(leg_wet_orig = leg_wet) %>%
#   mutate(leg_wet = as.factor(leg_wet2)) %>%
#   dplyr::select(!leg_wet2)
# saveRDS(brown52B, "Output Dataset Files/brown52B_imm_shifted.Rdata")
# brown52B <- brown52B %>% dplyr::select(-leg_wet_orig)
# tracks.reg.list5[[
#   which(sapply(tracks.reg.list5, function(x){x$Bird[1] == "52B brown"}))]] <-
#   brown52B
# tracks.comb3 <- do.call(rbind, tracks.reg.list5)
# saveRDS(tracks.reg.list5, "Output Dataset Files/BBA_tracks_reg_list5.RData")
# 
# # Check in plot
# ggplot(brown52B, 
#        aes(x = DateTime, 
#            y = ifelse((leg_wet_HMM2 == "wet"), 1, NA))) + 
#   geom_vline(xintercept = brown52B$DateTime[brown52B$leg_wet == "wet"], 
#              col = "lightblue", alpha = 0.1) +
#   geom_point(size = 1) +
#   labs(title = paste0(brown52B$Bird[1], ": GLS wet (blue) vs HMM wet (black dots)"), 
#        x = "DateTime", y = "")





################# 13. FUNCTIONS TO CORRECT IMMERSION AS NECESSARY ########################

# Current rules (as of 16.7.2024):
# (1) If bird has no GLS and low res GPS, leave blank (until ACC data analysed)
# (2) If bird has GLS but low res GPS, go with GLS
# (3) For all other birds (high res GPS), go with HMM*
# (4) Then, for all those birds EXCEPT GLS_wrong birds, correct with GLS dry if available
# *If bird has res 5, go with HMM2
# *If bird has res 10, go with HMM3.1 (2 dry states)
# *If bird has gaps (HMM = NA), go with GLS

# Create corrected immersion imm_corr dataframe
# Add columns leg_wet_HMM and leg_wet_C

# Create function to correct immersion and add interaction
# GLS_wrong = vector of bird names, for those with very high GLS dry HMM wet error
correct.imm.reg <- function(track, GLS_wrong){
  
  # What bird and resolution
  bird <- track$Bird[1]
  res <- track$resolution[1]
  
  # Coerce leg_wet to character
  track$leg_wet <- as.character(track$leg_wet)
  
  # Get corresponding immersion dataframe, if bird had a GLS
  imm_index <- which(sapply(imm_data, function(x){x$Bird[1]}) == bird)
  if(length(imm_index) > 0){
    imm_df <- imm_data[[imm_index]] %>%
      mutate(leg_wet = as.character(leg_wet),
             ID = Bird, leg_wet_HMM2 = NA, leg_wet_HMM3.1 = NA, leg_wet_HMM3.2 = NA,
             leg_wet_C = NA, follow = NA, stop = NA) %>%
      relocate(ID)
    
  # Else if no GLS, create a new df to match other imm_df
  # This one will have whatever resolution the GPS had
  } else {
    imm_df <- data.frame(ID = track$Bird,
                         Bird = track$Bird,
                         DateTime = track$DateTime,
                         leg_wet = as.character(track$leg_wet),
                         leg_wet_HMM2 = track$leg_wet_HMM2,
                         leg_wet_HMM3.1 = track$leg_wet_HMM3.1,
                         leg_wet_HMM3.2 = track$leg_wet_HMM3.2,
                         leg_wet_C = NA, 
                         follow = NA, stop = NA)
  }
  
  # (1) If no GLS and low res GPS, leave leg_wet_C blank (until ACC data analysed)
  if(is.na(track$leg_wet[1]) & res > 10){
    imm_df <- imm_df
  
  # (2) If GLS but low res GPS, go with GLS
  } else if(!is.na(track$leg_wet[1]) & res > 10){
    imm_df$leg_wet_C <- imm_df$leg_wet

  # (3) For all other birds (high res GPS), go with HMM
  } else {
    
    # Full join track df to imm df
    df1 <- full_join(track[ , -which(colnames(track) %in% 
                              c("dt", "Bird", "leg_wet", "Latitude", "Longitude"))], 
                     imm_df[ , -which(colnames(imm_df) %in% 
                              c("ID", "leg_wet_HMM2", "leg_wet_HMM3.1", "leg_wet_HMM3.2", "leg_wet_C"))], 
                     by = "DateTime")
    df1 <- df1 %>% 
      arrange(DateTime) %>%
      fill(ID, .direction = "downup") %>%
      relocate(leg_wet, .before = leg_wet_HMM2) %>%
      relocate(leg_wet_C, .after = leg_wet_HMM3.2) %>%
      relocate(Bird, .before = DateTime)
    
    # Remove rows at beginning and end outside start/end of track
    # (Mainly for testing HMM leg_wet, and removing bits on land)
    if(min(which(!is.na(df1$resolution))) > 1){
      rows <- (min(which(!is.na(df1$resolution))) - 1):nrow(df1)
      df1 <- df1[rows, ]
    }
    if(max(which(!is.na(df1$resolution))) < nrow(df1)){
      rows <- 1:(max(which(!is.na(df1$resolution))) + 1)
      df1 <- df1[rows, ]
    }
    
    # Loop through rows to fill in leg_wet_HMM columns with nearest wet/dry value
    for(i in 1:nrow(df1)){
      
      # If there is already a wet/dry value, go to next row
      if(!is.na(df1$leg_wet_HMM2[i])){
        next
        
      # Else if there is an <NA>
      } else {
        
        # If in a noHMM sub-track, leave as NA and go to next row
        if(grepl("noHMM", df1$ID[i])){
          next
        }
        
        # Figure out previous DateTime with a wet/dry value
        for(j in 1:nrow(df1)){
          if(i - j > 0){
            if(!is.na(df1$leg_wet_HMM2[i-j])){
              prevDT <- df1$DateTime[i-j]
              break
            } else {next}
          } else {break}
        }
        
        # Figure out next DateTime with a wet/dry value
        for(k in 1:nrow(df1)){
          if(i + k <= nrow(df1)){
            if(!is.na(df1$leg_wet_HMM2[i+k])){
              nextDT <- df1$DateTime[i+k]
              break
            } else {next}
          } else {break}
        }
      }
      
      # Fill leg_wet_HMM columns with value from closest row in time
      if(exists("prevDT") & exists("nextDT")){
        if(as.numeric(abs(df1$DateTime[i] - prevDT)) <= 
           as.numeric(abs(df1$DateTime[i] - nextDT))){
          df1$leg_wet_HMM2[i] <- df1$leg_wet_HMM2[df1$DateTime == prevDT]
          df1$leg_wet_HMM3.1[i] <- df1$leg_wet_HMM3.1[df1$DateTime == prevDT]
          df1$leg_wet_HMM3.2[i] <- df1$leg_wet_HMM3.2[df1$DateTime == prevDT]
        } else {
          df1$leg_wet_HMM2[i] <- df1$leg_wet_HMM2[df1$DateTime == nextDT]
          df1$leg_wet_HMM3.1[i] <- df1$leg_wet_HMM3.1[df1$DateTime == nextDT]
          df1$leg_wet_HMM3.2[i] <- df1$leg_wet_HMM3.2[df1$DateTime == nextDT]
        }
      
      # If at beginning of df
      } else if (!exists("prevDT") & exists("nextDT")) {
        df1$leg_wet_HMM2[i] <- df1$leg_wet_HMM2[df1$DateTime == nextDT]
        df1$leg_wet_HMM3.1[i] <- df1$leg_wet_HMM3.1[df1$DateTime == nextDT]
        df1$leg_wet_HMM3.2[i] <- df1$leg_wet_HMM3.2[df1$DateTime == nextDT]
        
      # # If at end of df  
      # } else if(exists("prevDT") & !exists("nextDT")){
      #   df1$leg_wet_HMM2[i] <- df1$leg_wet_HMM2[df1$DateTime == prevDT]
      #   df1$leg_wet_HMM3.1[i] <- df1$leg_wet_HMM3.1[df1$DateTime == prevDT]
      #   df1$leg_wet_HMM3.2[i] <- df1$leg_wet_HMM3.2[df1$DateTime == prevDT]
        
      # If all else fails, set leg_wet_HMM columns to NA (should not happen)
      } else {
        df1$leg_wet_HMM2[i] <- NA
        df1$leg_wet_HMM3.1[i] <- NA
        df1$leg_wet_HMM3.2[i] <- NA
      }
    }
    
    # Remove rows not originally from imm_df
    originals <- which(df1$DateTime %in% imm_df$DateTime)
    imm_df <- df1[originals, -which(colnames(df1) == "resolution")]
    
    # Update leg_wet_C:
    # Use 3-state HMM (2 dry states) for res 10
    # Use 2-state HMM for res 1 and 5
    if(res == 10){
      imm_df$leg_wet_C <- imm_df$leg_wet_HMM3.1
    } else if (res < 10 & !is.na(imm_df$leg_wet_HMM2[1])) {
      imm_df$leg_wet_C <- imm_df$leg_wet_HMM2
    } else {
      imm_df$leg_wet_C <- imm_df$leg_wet
    }
    
    # For noHMM sub-tracks, fill in leg_wet_C with GLS data (leg_wet)
    imm_df$leg_wet_C <- ifelse(is.na(imm_df$leg_wet_C) & !is.na(imm_df$leg_wet[1]), 
                               imm_df$leg_wet, imm_df$leg_wet_C)
  }
  
  # (4) For all those birds EXCEPT GLS_wrong birds, correct with GLS dry if available
  if(!(bird %in% GLS_wrong) & !is.na(imm_df$leg_wet[1])){
    imm_df$leg_wet_C <- ifelse(imm_df$leg_wet == "dry", "dry", imm_df$leg_wet_C)
  }
  
  
  ### Add in corrected interaction
  
  # Get interaction times for this bird
  int_bird <- int_summary[int_summary$Bird == bird, 1:7]
  
  # If the bird had interaction
  # (Must have had either immersion or GPS <= 1 min res, excludes accelerometer birds for now)
  if(nrow(int_bird) > 0 & res <= 60){
    
    for(i in 1:nrow(int_bird)){
      
      # If bird is following vessel, add Int_Num to track and immersion dataframes
      if(int_bird$Int_Class[i] == "Follow"){
        imm_df_Fvec <- ifelse(imm_df$DateTime >= int_bird$Real_start[i] &
                                imm_df$DateTime <= int_bird$Real_end[i],
                              int_bird$Int_Num[i], NA)
        #imm_df$follow <- ifelse(is.na(imm_df$follow), imm_df_Fvec, imm_df$follow)
        imm_df$follow <- imm_df_Fvec
        
      # Else if bird stops near vessel, add Int_Num to track and immersion dataframes
      } else if(int_bird$Int_Class[i] == "Stop"){
        imm_df_Svec <- ifelse(imm_df$DateTime >= int_bird$Real_start[i] &
                                imm_df$DateTime <= int_bird$Real_end[i],
                              int_bird$Int_Num[i], NA)
        #imm_df$stop <- ifelse(is.na(imm_df$stop), imm_df_Svec, imm_df$stop)
        imm_df$stop <- imm_df_Svec
      }
    }
  }
  
  # Return updated dataframe
  message(paste0(bird, " done."))
  return(imm_df)
}




# Create function to correct immersion and add interaction for tracks
# GLS_wrong = vector of bird names, for those with very high GLS dry HMM wet error
correct.imm.tracks <- function(track, GLS_wrong){
  
  # What bird and resolution
  bird <- track$Bird[1]
  res <- track$resolution[1]
  
  # Coerce leg_wet to character
  track$leg_wet <- as.character(track$leg_wet)
  
  # (1) If no GLS and low res GPS, leave blank (until ACC data analysed)
  if(is.na(track$leg_wet[1]) & res > 10){
    track <- track
  
  # (2) If GLS but low res GPS, go with GLS
  } else if(!is.na(track$leg_wet[1]) & res > 10){
    track$leg_wet_C <- track$leg_wet
  
  # (3) For all other birds (high res GPS), go with HMM*
  } else {
    
    # Update leg_wet_C:
    # Use 3-state HMM (2 dry states) for res 10
    # Use 2-state HMM for res 1 and 5
    if(res == 10){
      track$leg_wet_C <- track$leg_wet_HMM3.1
    } else if (res < 10 & !is.na(track$leg_wet_HMM2[1])) {
      track$leg_wet_C <- track$leg_wet_HMM2
    } else {
      track$leg_wet_C <- track$leg_wet
    }
    
    # For noHMM sub-tracks, fill in leg_wet_C with GLS data (leg_wet)
    track$leg_wet_C <- ifelse(is.na(track$leg_wet_C) & !is.na(track$leg_wet[1]), 
                              track$leg_wet, track$leg_wet_C)
    
  }
  
  # (4) For all those birds EXCEPT GLS_wrong birds, correct with GLS dry if available
  if(!(bird %in% GLS_wrong) & !is.na(track$leg_wet[1])){
    track$leg_wet_C <- ifelse(track$leg_wet == "dry", "dry", track$leg_wet_C)
  }
  

  
  ### Add in corrected interaction
  
  # Add follow and stop columns
  track <- track %>% mutate(follow = NA, stop = NA)
  
  # Get interaction times for this bird
  int_bird <- int_summary[int_summary$Bird == bird, 1:7]
  
  # If the bird had interaction
  # (Must have had either immersion or GPS <= 1 min res, excludes accelerometer birds for now)
  if(nrow(int_bird) > 0 & track$resolution[1] <= 60){
    
    for(i in 1:nrow(int_bird)){
      
      # If bird is following vessel, add Int_Num to track and immersion dataframes
      if(int_bird$Int_Class[i] == "Follow"){
        track_Fvec <- ifelse(track$DateTime >= int_bird$Real_start[i] &
                               track$DateTime <= int_bird$Real_end[i],
                             int_bird$Int_Num[i], NA)
        #track$follow <- ifelse(is.na(track$follow), track_Fvec, track$follow)
        track$follow <- track_Fvec
        
        # Else if bird stops near vessel, add Int_Num to track and immersion dataframes
      } else if(int_bird$Int_Class[i] == "Stop"){
        track_Svec <- ifelse(track$DateTime >= int_bird$Real_start[i] &
                               track$DateTime <= int_bird$Real_end[i],
                             int_bird$Int_Num[i], NA)
        #track$stop <- ifelse(is.na(track$stop), track_Svec, track$stop)
        track$stop <- track_Svec
      }
    }
  }
  
  # Return updated dataframe
  message(paste0(bird, " done."))
  return(track)
}









############### 14. RUN THE CORRECTION FUNCTIONS ##########################

# Prep for parallel processing
numCores <- 2 #detectCores()
cl <- makeCluster(numCores)
clusterEvalQ(cl, {
  library(tidyverse)
  library(lubridate)
  library(pbapply)
})
clusterExport(cl = cl, c("wd", "WGS1984", "tracks.reg.list5", "imm_data", "int_summary"))

# Correct immersion
corrected_imm <- pblapply(tracks.reg.list5, correct.imm.reg, GLS_wrong = "53T blue", cl = cl)

# Correct tracks
tracks.reg.list6 <- pblapply(tracks.reg.list5, correct.imm.tracks, GLS_wrong = "53T blue", cl = cl)

# Stop the cluster
stopCluster(cl)

# Double check this worked
imm.comb1 <- do.call(rbind, corrected_imm)
tracks.comb4 <- do.call(rbind, tracks.reg.list6)
as.data.frame(imm.comb1 %>%
  group_by(Bird) %>%
  summarise(n = n(),
            NAlegwet = sum(is.na(leg_wet)),
            NAlegwetHMM = sum(is.na(leg_wet_HMM2)),
            NAlegwetC = sum(is.na(leg_wet_C)),
            follow = sum(follow, na.rm = TRUE)))
as.data.frame(tracks.comb4 %>%
                group_by(Bird) %>%
                summarise(n = n(),
                          res = resolution[1],
                          NAlegwet = sum(is.na(leg_wet)),
                          NAlegwetHMM = sum(is.na(leg_wet_HMM2)),
                          NAlegwetC = sum(is.na(leg_wet_C)),
                          follow = sum(follow, na.rm = TRUE)))

# Save files
saveRDS(corrected_imm, "Output Dataset Files/Wet_Data_Reg_Corr3.RData")
saveRDS(tracks.reg.list6, "Output Dataset Files/BBA_Track_List_Imm_Int_Corr3.RData")








############### 15. FINAL RESULTS TO QUANTIFY GLS ERROR #################################

# # If necessary, read in percents df and corrected_imm df
# percents <- readRDS("Output Dataset Files/imm_corr_percents.Rdata")
# corrected_imm <- readRDS("Output Dataset Files/Wet_Data_Reg_Corr3.RData")
# tracks.reg.list6 <- readRDS("Output Dataset Files/BBA_Track_List_Imm_Int_Corr3.RData")
# tracks.comb4 <- do.call(rbind, tracks.reg.list6)

# Add columns to percents df, to compare leg_wet (GLS) with leg_wet_C (final corrected)
GLSwet_Cdry_vec <- sapply(percents$Bird, function(x){
  i <- which(sapply(tracks.reg.list6, function(y){y$Bird[1] == x}))
  imm_df <- tracks.reg.list6[[i]]
  imm_check <- ifelse(imm_df$leg_wet == "wet" & imm_df$leg_wet_C == "dry", TRUE, FALSE)
  prop <- sum(imm_check)/nrow(imm_df) # Proportion of whole track
  return(prop)
})
GLSdry_Cwet_vec <- sapply(percents$Bird, function(x){
  i <- which(sapply(tracks.reg.list6, function(y){y$Bird[1] == x}))
  imm_df <- tracks.reg.list6[[i]]
  imm_check <- ifelse(imm_df$leg_wet == "dry" & imm_df$leg_wet_C == "wet", TRUE, FALSE)
  prop <- sum(imm_check)/nrow(imm_df) # Proportion of whole track
  return(prop)
})
GLSwet_Cdry_p_vec <- sapply(percents$Bird, function(x){
  i <- which(sapply(tracks.reg.list6, function(y){y$Bird[1] == x}))
  imm_df <- tracks.reg.list6[[i]]
  imm_check <- ifelse(imm_df$leg_wet == "wet" & imm_df$leg_wet_C == "dry", TRUE, FALSE)
  prop <- sum(imm_check)/nrow(imm_df[imm_df$leg_wet == "wet",]) # Proportion of GLS wet
  return(prop)
})
percents <- percents %>%
  mutate(GLSwet_Cdry = GLSwet_Cdry_vec,
         GLSdry_Cwet = GLSdry_Cwet_vec,
         GLSwet_Cdry_p = GLSwet_Cdry_p_vec)

# Plot GLS error as compared to "correct" immersion
ggplot(percents, aes(x = Bird, 
                     y = GLSwet_Cdry,
                     fill = as.factor(resolution))) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major.x = element_blank()) +
  labs(y = "Proportion of timesteps\nwith GLS wet when actually dry",
       fill = "Fix interval (s)")
ggplot(percents, aes(x = Bird, 
                     y = GLSwet_Cdry_p,
                     fill = as.factor(resolution))) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major.x = element_blank()) +
  labs(y = "Proportion of GLS wet timesteps\nthat are actually dry",
       fill = "Fix interval (s)")

# Descriptive statistics

# Calculate amount of time after take-off that GLS tends to stay wet
wet_time_list <- pblapply(tracks.reg.list6, function(df){ # Can do with tracks.reg.list6 or corrected_imm
  
  # Make sure df has required data
  if(sum(is.na(df$leg_wet_C)) > 0 | sum(is.na(df$leg_wet)) > 0){
    message(paste0(df$Bird[1], " had NAs in leg_wet or leg_wet_C."))
    return(NULL)
  }
  
  # Identify takeoffs
  df <- df %>% 
    mutate(takeoff = ifelse(df$leg_wet_C == "dry" & lag(df$leg_wet_C) == "wet", 1, 0),
           landing = ifelse(df$leg_wet_C == "wet" & lag(df$leg_wet_C) == "dry", 1, 0)) %>%
    fill(takeoff, landing, .direction = "updown")
  
  # Calculate amount of time until GLS goes dry again
  wet_time <- c()
  takeoffs <- which(df$takeoff == 1)
  landings <- which(df$landing == 1)
  for(i in takeoffs){
    if(i == max(takeoffs)){
      nearest_landing <- nrow(df)
    } else {
      nearest_landing <- min(landings[landings > i])
    }
    flight <- df[i:nearest_landing,]
    if(sum(flight$leg_wet == "dry") == 0){
      t_until <- abs(as.numeric(difftime(max(flight$DateTime), min(flight$DateTime), units = "secs")))
    } else {
      t_until <- abs(as.numeric(difftime(min(flight$DateTime[flight$leg_wet == "dry"]), 
                                         min(flight$DateTime), units = "secs")))
    }
    wet_time <- c(wet_time, t_until)
  }
  
  # Return dataframe of durations until GLS goes dry after takeoff
  message(paste0(df$Bird[1], " done."))
  new_df <- data.frame(Bird = df$Bird[1], wet_time = wet_time)
  return(new_df)
})

# Visualise amount of time after takeoff that GLS stays wet, assuming leg_wet_C is correct
# (This could be an underestimate as once next takeoff happens the clock restarts)
wet_time_comb <- do.call(rbind, wet_time_list)
wet_time_bird <- data.frame(Bird = sapply(corrected_imm, function(x){x$Bird[1]}),
                            mean_wet_time = sapply(wet_time_list, function(x){mean(x$wet_time)}),
                            sd_wet_time = sapply(wet_time_list, function(x){sd(x$wet_time)}),
                            max_wet_time = sapply(wet_time_list, function(x){max(x$wet_time)})) %>%
  filter(!is.na(mean_wet_time), Bird %in% percents$Bird) %>%
  mutate(resolution = percents$resolution[percents$Bird %in% Bird])

# Quantify GLS error
sum(percents$GLSwet_Cdry_p > 0.05,
    na.rm= TRUE)/sum(!is.na(percents$GLSwet_Cdry_p))
# 63.2% of birds (24) with GLS and high res have >5% GLSwet with corrected dry
sum(percents$GLSwet_Cdry > 0.05,
    na.rm = TRUE)/sum(!is.na(percents$GLSwet_Cdry))
# 39.5% of birds (24) with GLS and high res have >5% GLSwet with corrected dry
tracks.comb4 %>%
  filter(!is.na(leg_wet)) %>%
  group_by(leg_wet) %>%
  summarise(Cdry = sum(leg_wet_C == "dry")/n())
# 12.6% of GLS wet timesteps were in fact dry
wet_time_comb %>%
  summarise(mean = mean(wet_time), median = median(wet_time), sd = sd(wet_time), max = max(wet_time))

# Plots
ggplot(percents, aes(y = GLSwet_Cdry_p)) +
  geom_boxplot() +
  labs(x = "Proportion per bird", y = "") +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
ggplot(wet_time_comb[wet_time_comb$wet_time < 150, ], aes(x = wet_time)) +
  geom_histogram(bins = 30) +
  labs(x = "Time GLS wet following takeoff (s)", y = "Frequency")
ggplot(wet_time_bird, aes(x = Bird, 
                     y = mean_wet_time,
                     fill = as.factor(resolution))) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major.x = element_blank()) +
  labs(y = "Mean time GLS wet following takeoff",
       fill = "Fix interval (s)")
ggplot(wet_time_comb, aes(x = Bird, y = wet_time)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major.x = element_blank()) +
  labs(y = "Time GLS wet following takeoff (s)")
