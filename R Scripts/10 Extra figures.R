###########################################################
########### EXTRA MANUSCRIPT FIGURES ######################
###########################################################

# Produces the figures for the Supporting Information of this study.

############### CONTENTS #######################

# 1. Packages and WD
# 2. Read in imm_metrics data
# 3. Extra maps
# 4. Time series of regularity
# 5. Time series of random forest probabilities
# 6. False positive vessel detections
# 7. False negative vessel detections




################# 1. Packages and WD ##########################

# Packages
library(tidyverse)
library(lubridate)
library(pbapply)
library(caret)
library(randomForest)
library(ggpubr)
library(gganimate)
library(move)
library(sf)
theme_set(theme_bw())

# Set WD
setwd("C:/Users/jonat/OneDrive/Desktop/University of Oxford DPhil Biology/DPhil Chapters/Ch 1 BBA-vessel interactions/Data and Scripts/")



#################### 2. READ IN IMM_METRICS DATA #########################

# Read in data
imm_metrics_list1 <- readRDS("Output Dataset Files/imm_metrics_list1.RData")
imm_metrics_df1 <- do.call(rbind, imm_metrics_list1)

# Filter to only intermediate immersion timestamps
imm_metrics_df2 <- imm_metrics_df1 %>%
  filter(imm_factor == "intermediate") %>%
  mutate(follow = ifelse(!is.na(follow), "F", "N")) %>%
  mutate(follow = relevel(as.factor(follow), ref = "N"))

# Narrow dataset to rows without NA
imm_metrics_df2.1 <- imm_metrics_df2[!is.na(imm_metrics_df2$mean_wet),]







########################### 3. EXTRA MAPS ######################################

# Mapping packages
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(sp)
library(sf)

# Read in the bird data
BBA_int_list <- readRDS("Output Dataset Files/BBA_Track_List_Imm_Int.RData")

# Read in the vessel data
AIS_Tracks <- read.csv("Output Dataset Files/BBA_gfw_vessels_redacted.csv", header = TRUE) %>%
  mutate(Bird = as.factor(Bird), vessel_class = as.factor(vessel_class), fishing = as.factor(fishing), SubTraj = as.factor(SubTraj)) %>%
  mutate(DateTime = ymd_hms(DateTime, tz = "UTC"))

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


# Get first and last DateTime of interactions, find nearest vessel DateTimes, plot together
# x is a bird track dataframe with interactions numbered in the column "interaction"
# polygon is a land polygon
plot.bv.int.map <- function(x, polygon){
  
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
        
        ggplot(data = polygon) + 
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
                          #"MMSI: ", Btrack$ssvid[1], "\n",
                          #"Flag: ", Btrack$flag[1], "\n", 
                          "Class: ", Btrack$vessel_class[1], "\n",
                          "Active fishing: ", ifelse(length(Btrack$fishing[Btrack$fishing == 1]) > 0, "Yes", "No"), "\n",
                          sep = ""), 
            hjust = 0, vjust = 0, size = 3.5) +
          geom_text(aes(label = with_tz(min_Bdt, "America/Argentina/Buenos_Aires"),
                        x = Btrack$Longitude[1], y = Btrack$Latitude[1]),
                    position = position_dodge(0.05*rangeLong),
                    size = 3.5) +
          geom_text(aes(label = with_tz(max_Bdt, "America/Argentina/Buenos_Aires"),
                        x = Btrack$Longitude[nrow(Btrack)], y = Btrack$Latitude[nrow(Btrack)]),
                    size = 3.5,
                    position = position_dodge(0.05*rangeLong)) +
          scale_x_continuous(breaks = pretty(Btrack$Longitude, n = 4)) +
          scale_y_continuous(breaks = pretty(Btrack$Latitude, n = 4)) +
          theme_bw() +
          theme(text = element_text(size = 12))
        
        # Plot if there is no GLS leg data  
      } else {
        
        ggplot(data = polygon) + 
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
                          #"MMSI: ", Btrack$ssvid[1], "\n",
                          #"Flag: ", Btrack$flag[1], "\n", 
                          "Class: ", Btrack$vessel_class[1], "\n",
                          "Active fishing: ", ifelse(length(Btrack$fishing[Btrack$fishing == 1]) > 0, "Yes", "No"), "\n",
                          sep = ""), 
            hjust = 0, vjust = 0, size = 3.5) +
          geom_text(aes(label = with_tz(min_Bdt, "America/Argentina/Buenos_Aires"),
                        x = Btrack$Longitude[1], y = Btrack$Latitude[1]),
                    position = position_dodge(0.05*rangeLong),
                    size = 3.5) +
          geom_text(aes(label = with_tz(max_Bdt, "America/Argentina/Buenos_Aires"),
                        x = Btrack$Longitude[nrow(Btrack)], y = Btrack$Latitude[nrow(Btrack)]),
                    size = 3.5,
                    position = position_dodge(0.05*rangeLong)) +
          scale_x_continuous(breaks = pretty(Btrack$Longitude, n = 4)) +
          scale_y_continuous(breaks = pretty(Btrack$Latitude, n = 4)) +
          theme_bw() +
          theme(text = element_text(size = 12))
        
      }
    })
    
    # Else if no interactions, return NULL  
  } else {
    return(NULL)
  }
}

# Get South America polygons, if haven't already done
s_america <- ne_countries(scale = "large", returnclass = "sf", continent = "South America")


# Figure S1: 17B brown Int7
brown17b <- which(sapply(BBA_int_list, function(y){y$Bird[1]}) == "17B brown")
plot.bv.int.map(BBA_int_list[[brown17b]], polygon = s_america)[[7]]

# To assist with Figure S6: 74F red
red74f <- which(sapply(BBA_int_list, function(y){y$Bird[1]}) == "74F red")
maps <- plot.bv.int.map(BBA_int_list[[red74f]], polygon = s_america)






##################### 4. TIME SERIES OF REGULARITY ######################################

# Plot sensitivity and specificity example
df <- readRDS("Output Dataset Files/all_results_list3.RData")[[7]][[1]]
df %>% ggplot(aes(x = threshold)) +
  geom_line(aes(y = sensitivity), col = "blue") +
  geom_line(aes(y = false_pos_percent), col = "red") +
  ylim(0, 1) +
  labs(x = "Regularity threshold", y = "Sensitivity (blue) | FP Ratio (red)")

# Bird indices
birds <- sapply(imm_metrics_list1, function(x){x$Bird[1]})

# 75F red
F75 <- which(sapply(imm_metrics_list1, function(x){x$Bird[1]}) == "75F red")
ggplot() +
  geom_line(data = imm_metrics_list1[[F75]][6000:19000,], 
            aes(x = DateTime, y = P_LWLD)) +
  geom_rect(data = int_summary[int_summary$Bird == imm_metrics_list1[[F75]]$Bird[1], ],
            aes(xmin = Real_start,
                xmax = Real_end,
                ymin = -Inf,
                ymax = Inf),
            alpha = .2, fill = "green") +
  labs(title = imm_metrics_list1[[F75]]$Bird[1], y = "Regularity of Immersion") +
  theme_bw() +
  theme(text = element_text(size = 13))

# # 74F red
# ggplot() +
#   geom_line(data = imm_metrics_list1[[34]], 
#             aes(x = DateTime, y = P_LWLD)) +
#   geom_rect(data = int_summary[int_summary$Bird == imm_metrics_list1[[34]]$Bird[1], ],
#             aes(xmin = Real_start,
#                 xmax = Real_end,
#                 ymin = -Inf,
#                 ymax = Inf),
#             alpha = .2, fill = "green") +
#   labs(title = imm_metrics_list1[[34]]$Bird[1], y = "Regularity of Immersion") +
#   theme_bw() +
#   theme(text = element_text(size = 13))

# View multiple birds: vessel following and P_LWLD
ggplot() +
  geom_vline(data = imm_metrics_df1[imm_metrics_df1$Bird %in% birds[1:9] & !is.na(imm_metrics_df1$follow), ], 
             aes(xintercept = DateTime), col = "lightgreen", alpha = 0.3) +
  geom_line(data = imm_metrics_df1[imm_metrics_df1$Bird %in% birds[1:9], ],
            aes(x = DateTime, y = P_LWLD)) +
  facet_wrap(~Bird, scales = "free_x", ncol = 1) +
  labs(y = "Regularity of immersion")
# ggplot() +
#   geom_vline(data = imm_metrics_df1[imm_metrics_df1$Bird %in% birds[10:18] & !is.na(imm_metrics_df1$follow), ], 
#              aes(xintercept = DateTime), col = "lightgreen", alpha = 0.3) +
#   geom_line(data = imm_metrics_df1[imm_metrics_df1$Bird %in% birds[10:18], ],
#             aes(x = DateTime, y = P_LWLD)) +
#   facet_wrap(~Bird, scales = "free_x", ncol = 1) +
#   labs(y = "Regularity of immersion")
# ggplot() +
#   geom_vline(data = imm_metrics_df1[imm_metrics_df1$Bird %in% birds[19:27] & !is.na(imm_metrics_df1$follow), ], 
#              aes(xintercept = DateTime), col = "lightgreen", alpha = 0.3) +
#   geom_line(data = imm_metrics_df1[imm_metrics_df1$Bird %in% birds[19:27], ],
#             aes(x = DateTime, y = P_LWLD)) +
#   facet_wrap(~Bird, scales = "free_x", ncol = 1) +
#   labs(y = "Regularity of immersion")
# ggplot() +
#   geom_vline(data = imm_metrics_df1[imm_metrics_df1$Bird %in% birds[28:36] & !is.na(imm_metrics_df1$follow), ], 
#              aes(xintercept = DateTime), col = "lightgreen", alpha = 0.3) +
#   geom_line(data = imm_metrics_df1[imm_metrics_df1$Bird %in% birds[28:36], ],
#             aes(x = DateTime, y = P_LWLD)) +
#   facet_wrap(~Bird, scales = "free_x", ncol = 1) +
#   labs(y = "Regularity of immersion")
# ggplot() +
#   geom_vline(data = imm_metrics_df1[imm_metrics_df1$Bird %in% birds[37:45] & !is.na(imm_metrics_df1$follow), ], 
#              aes(xintercept = DateTime), col = "lightgreen", alpha = 0.3) +
#   geom_line(data = imm_metrics_df1[imm_metrics_df1$Bird %in% birds[37:45], ],
#             aes(x = DateTime, y = P_LWLD)) +
#   facet_wrap(~Bird, scales = "free_x", ncol = 1) +
#   labs(y = "Regularity of immersion")
# 
# # View select birds
# select_birds <- birds[c(8,9,11,16,31,33,35,37,41)]
# ggplot() +
#   geom_vline(data = imm_metrics_df1[imm_metrics_df1$Bird %in% select_birds & !is.na(imm_metrics_df1$follow), ], 
#              aes(xintercept = DateTime), col = "lightgreen", alpha = 0.3) +
#   geom_line(data = imm_metrics_df1[imm_metrics_df1$Bird %in% select_birds, ],
#             aes(x = DateTime, y = P_LWLD)) +
#   facet_wrap(~Bird, scales = "free_x", ncol = 1) +
#   labs(y = "Regularity of immersion")





################ 5. TIME SERIES OF RANDOM FOREST PROBABILITIES ############################

# Read in rf models (time-based, all predictors)
rf2 <- readRDS("Output Dataset Files/rf2.RData") # Corrected immersion
#rf4 <- readRDS("Output Dataset Files/rf4.RData") # Uncorrected immersion

# Put rf probabilities back into original data
rf2_pred <- rf2$pred %>% arrange(rowIndex)
imm_metrics_rf2 <- cbind(imm_metrics_df2.1, rf2_pred)

# Bird indices
birds <- unique(imm_metrics_rf2$Bird)

# View multiple birds: vessel following and random forest prediction
ggplot() +
  geom_vline(data = imm_metrics_rf2[imm_metrics_rf2$Bird %in% birds[1:9] & imm_metrics_rf2$follow == "F", ], 
             aes(xintercept = DateTime, col = pred), alpha = 0.2) +
  geom_line(data = imm_metrics_rf2[imm_metrics_rf2$Bird %in% birds[1:9], ],
            aes(x = DateTime, y = F)) +
  facet_wrap(~Bird, scales = "free_x", ncol = 1) +
  scale_color_manual(
    values = c("N" = "orange", "F" = "turquoise2"), # Define colours for N and F
    labels = c("N" = "Not predicted", "F" = "Predicted")) +
  labs(y = "Probability of following (random forest)", col = "Confirmed\nfollowing")
# ggplot() +
#   geom_vline(data = imm_metrics_rf2[imm_metrics_rf2$Bird %in% birds[10:18] & imm_metrics_rf2$follow == "F", ], 
#              aes(xintercept = DateTime, col = pred), alpha = 0.2) +
#   geom_line(data = imm_metrics_rf2[imm_metrics_rf2$Bird %in% birds[10:18], ],
#             aes(x = DateTime, y = F)) +
#   facet_wrap(~Bird, scales = "free_x", ncol = 1) +
#   scale_color_manual(
#     values = c("N" = "orange", "F" = "turquoise2"), # Define colours for N and F
#     labels = c("N" = "Not predicted", "F" = "Predicted")) +
#   labs(y = "Probability of following (random forest)", col = "Confirmed\nfollowing")
# ggplot() +
#   geom_vline(data = imm_metrics_rf2[imm_metrics_rf2$Bird %in% birds[19:27] & imm_metrics_rf2$follow == "F", ], 
#              aes(xintercept = DateTime, col = pred), alpha = 0.2) +
#   geom_line(data = imm_metrics_rf2[imm_metrics_rf2$Bird %in% birds[19:27], ],
#             aes(x = DateTime, y = F)) +
#   facet_wrap(~Bird, scales = "free_x", ncol = 1) +
#   scale_color_manual(
#     values = c("N" = "orange", "F" = "turquoise2"), # Define colours for N and F
#     labels = c("N" = "Not predicted", "F" = "Predicted")) +
#   labs(y = "Probability of following (random forest)", col = "Confirmed\nfollowing")
# ggplot() +
#   geom_vline(data = imm_metrics_rf2[imm_metrics_rf2$Bird %in% birds[28:36] & imm_metrics_rf2$follow == "F", ], 
#              aes(xintercept = DateTime, col = pred), alpha = 0.2) +
#   geom_line(data = imm_metrics_rf2[imm_metrics_rf2$Bird %in% birds[28:36], ],
#             aes(x = DateTime, y = F)) +
#   facet_wrap(~Bird, scales = "free_x", ncol = 1) +
#   scale_color_manual(
#     values = c("N" = "orange", "F" = "turquoise2"), # Define colours for N and F
#     labels = c("N" = "Not predicted", "F" = "Predicted")) +
#   labs(y = "Probability of following (random forest)", col = "Confirmed\nfollowing")
# ggplot() +
#   geom_vline(data = imm_metrics_rf2[imm_metrics_rf2$Bird %in% birds[37:45] & imm_metrics_rf2$follow == "F", ], 
#              aes(xintercept = DateTime, col = pred), alpha = 0.2) +
#   geom_line(data = imm_metrics_rf2[imm_metrics_rf2$Bird %in% birds[37:45], ],
#             aes(x = DateTime, y = F)) +
#   facet_wrap(~Bird, scales = "free_x", ncol = 1) +
#   scale_color_manual(
#     values = c("N" = "orange", "F" = "turquoise2"), # Define colours for N and F
#     labels = c("N" = "Not predicted", "F" = "Predicted")) +
#   labs(y = "Probability of following (random forest)", col = "Confirmed\nfollowing")
# 
# # View select birds
# select_birds <- birds[c(8,9,11,16,31,33,35,37,41)]
# ggplot() +
#   geom_vline(data = imm_metrics_rf2[imm_metrics_rf2$Bird %in% select_birds & imm_metrics_rf2$follow == "F", ], 
#              aes(xintercept = DateTime, col = pred), alpha = 0.2) +
#   geom_line(data = imm_metrics_rf2[imm_metrics_rf2$Bird %in% select_birds, ],
#             aes(x = DateTime, y = F)) +
#   facet_wrap(~Bird, scales = "free_x", ncol = 1) +
#   theme(legend.position = "none") +
#   labs(y = "Probability of following (random forest)", col = "Confirmed\nfollowing")









################# 6. FALSE POSTIVE VESSEL DETECTIONS ###################
# Including highly regular natural foraging

# Read in rf2 from previous section, plus this:
tracks.reg.list6 <- readRDS("Output Dataset Files/BBA_Track_List_Imm_Int_Corr3.RData")

# Identify the false positive periods (code as 2)
# Also include periods of highly regular natural foraging (code as 1)
imm_metrics_FP <- imm_metrics_rf2 %>%
  mutate(FP = ifelse(obs == "N" & P_LWLD > 0.9, 1, 0),
         FP = ifelse(obs == "N" & pred == "F", 2, FP), # False positives
         FP_buffer = 0)
imm_metrics_FPlist <- split(imm_metrics_FP, f = imm_metrics_FP$Bird)

# Buffer to 8 mins before and after each FP prediction
buffer <- (60*8)/
  as.numeric(difftime(imm_metrics_rf2$DateTime[2], imm_metrics_rf2$DateTime[1], units = "secs"))

# Buffer for each bird, then recombine into 1 dataframe
imm_metrics_FPlist <- lapply(imm_metrics_FPlist, function(x){
  for(i in 1:nrow(x)){
    # Fill in vector for new column (1 if within 8 mins of FP, otherwise 0)
    if(x$FP[i] > 0){
      start <- ifelse(i-buffer > 0, i-buffer, 1)
      end <- ifelse(i+buffer < nrow(x), i+buffer, nrow(x))
      x$FP_buffer[start:end] <- 1
    } else {
      next
    }
  }
  return(x)
})
imm_metrics_FP <- do.call(rbind, imm_metrics_FPlist)

# Create list of false positive "bouts"
imm_metrics_FP2 <- imm_metrics_FP %>%
  group_by(Bird) %>%
  mutate(tmp1 = ifelse(FP_buffer == 1 & (lag(FP_buffer) != 1 | is.na(lag(FP_buffer))), 
                       1, 0)) %>%
  ungroup() %>% as.data.frame() %>%
  filter(FP_buffer == 1) %>%
  mutate(FP_bout = cumsum(tmp1))
imm_metrics_FPlist2 <- lapply(split(imm_metrics_FP2, f = imm_metrics_FP2$FP_bout), 
                              function(x){x %>% dplyr::select(-tmp1)})

# Join GPS data onto those bouts
imm_metrics_FPlist2 <- pblapply(imm_metrics_FPlist2, function(df){
  
  # Which bird
  bird <- unique(df$Bird)
  if(length(bird) > 1){
    message(paste0("2 birds in one dataframe: ", paste(unique(df$Bird), collapse = ", ")))
    return(NULL)
  }
  
  # Get segment of GPS track from same bird and time period
  gps <- tracks.reg.list6[[which(sapply(tracks.reg.list6, function(x){x$Bird[1]}) == bird)]] %>%
    filter(DateTime >= min(df$DateTime) - minutes(1), 
           DateTime <= max(df$DateTime) + minutes(1))
  
  # Interpolate to exact same timesteps as immersion data
  gps_mv <- move(x = gps$Longitude,
                 y = gps$Latitude,
                 time = as.POSIXct(gps$DateTime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"),
                 proj = st_crs(4326)[[2]],
                 animal = gps$Bird,
                 removeDuplicatedTimestamps = TRUE)
  gps_interp <- interpolateTime(gps_mv,
                            time = df$DateTime,
                            spaceMethod = "greatcircle") %>%
    as.data.frame() %>%
    rename(Longitude = coords.x1, Latitude = coords.x2) %>%
    mutate(resolution = gps$resolution[1]) %>%
    dplyr::select(Latitude, Longitude, resolution)
  
  # Put Lat and Long columns into dataframe and return
  return(cbind(df, gps_interp))
  
})

# Maps
pblapply(imm_metrics_FPlist2, function(df){
  
  # Get boundary coordinates
  minLong <- min(df$Longitude)
  maxLong <- max(df$Longitude)
  rangeLong <- abs(maxLong - minLong)
  minLat <- min(df$Latitude)
  maxLat <- max(df$Latitude)
  rangeLat <- abs(maxLat - minLat)
  
  # Static map
  map1 <- ggplot(df, aes(x = Longitude, y = Latitude, col = leg_wet_C, group = 1)) +
    geom_path() +
    geom_point() +
    coord_sf(xlim = c(minLong, maxLong), 
             ylim = c(minLat, maxLat), 
             expand = TRUE, clip = "off") +
    labs(x = "Longitude", y = "Latitude", colour = "Corrected\nImmersion", 
         title = paste0(min(df$DateTime), " to ", max(df$DateTime))) +
    annotate(
      geom = "text", #x = maxLong + 0.08*rangeLong, y = minLat - 0.12*rangeLat, 
      x = minLong, y = minLat,
      label = paste("Bird: ", df$Bird[1], "\n",
                    "Resolution: ", df$resolution[1], "\n",
                    "False pos: ", ifelse(2 %in% df$FP, "Yes", "No"), "\n",
                    "Highly reg: ", ifelse(1 %in% df$FP, "Yes", "No"), "\n",
                    "Includes vessel: ", ifelse("F" %in% df$obs, "Yes", "No"), "\n",
                    "Max prob follow = ", max(df$F), "\n",
                    "Max regularity = ", round(max(df$P_LWLD), digits = 2), "\n",
                    sep = ""), 
      hjust = 0, vjust = 0, size = 3)
  ggsave(paste0("BBA False Positive Maps/FPbout", df$FP_bout[1], "_", gsub(" ", "", df$Bird[1]), ".jpg"),
         plot = map1, device = "jpeg")
  
  # Animation
  anim1 <- ggplot(df, aes(x = Longitude, y = Latitude, col = leg_wet_C, group = 1)) +
    geom_path() +
    geom_point() +
    coord_sf(xlim = c(minLong, maxLong),
             ylim = c(minLat, maxLat),
             expand = TRUE, clip = "off") +
    labs(x = "Longitude", y = "Latitude", colour = "Corrected\nImmersion", title = "DateTime (UTC): {frame_time}") +
    annotate(
      geom = "text", #x = maxLong + 0.08*rangeLong, y = minLat - 0.12*rangeLat,
      x = minLong, y = minLat,
      label = paste("Bird: ", df$Bird[1], "\n",
                    "Resolution: ", df$resolution[1], "\n",
                    "False pos: ", ifelse(2 %in% df$FP, "Yes", "No"), "\n",
                    "Highly reg: ", ifelse(1 %in% df$FP, "Yes", "No"), "\n",
                    "Includes vessel: ", ifelse("F" %in% df$obs, "Yes", "No"), "\n",
                    "Max prob follow = ", max(df$F), "\n",
                    "Max regularity = ", round(max(df$P_LWLD), digits = 2), "\n",
                    sep = ""),
      hjust = 0, vjust = 0, size = 3) +
    transition_time(DateTime) +
    ease_aes('linear') +
    shadow_trail(distance = 0.01, alpha = 0.05)
  anim_save(filename = paste0("BBA False Positive Animations/FPbout", df$FP_bout[1], "_", gsub(" ", "", df$Bird[1]), ".gif"),
            anim1, nframes = 200, duration = 20, renderer = gifski_renderer(),
            width = 20, height = 15, units = "cm", res = 250)
  
  # Message progress and return nothing
  message(paste("Finished anim for FP bout ", df$FP_bout[1], " (", df$Bird[1], ").\n", sep = ""))
  return(NULL)
  
})








################# 7. FALSE NEGATIVE VESSEL DETECTIONS ####################
# Including vessel following with low regularity

# Read in rf2 from above, plus this:
tracks.reg.list6 <- readRDS("Output Dataset Files/BBA_Track_List_Imm_Int_Corr3.RData")

# Identify the false negative periods (code as 2)
# Also include periods of not regular vessel following (code as 1)
imm_metrics_FN <- imm_metrics_rf2 %>%
  mutate(FN = ifelse(obs == "F" & P_LWLD < 0.9, 1, 0),
         FN = ifelse(obs == "F" & pred == "N", 2, FN), # False negatives
         buffer = 0) %>%
  dplyr::select(-follow) %>%
  left_join(imm_metrics_df1[,c("Bird", "DateTime","follow")], by = join_by(Bird, DateTime)) %>%
  rename(int_num = follow)

# Identify the confirmed vessel interactions that are false negative or not regular
int_FN_df <- imm_metrics_FN %>%
  filter(!is.na(int_num)) %>%
  group_by(Bird, int_num) %>%
  summarise(
    duration = difftime(max(DateTime), min(DateTime), units = "mins"),
    N_landings = sum(leg_wet_C == "wet" & (lag(leg_wet_C == "dry") | is.na(lag(leg_wet_C)))),
    prop_predicted = sum(pred == "F")/n(),
    prop_FN = sum(FN == 2)/n(),
    mean_reg = mean(P_LWLD),
    max_reg = max(P_LWLD)
  ) %>% 
  arrange(desc(prop_FN)) %>%
  mutate_if(is.numeric, ~ round(.,3)) %>%
  mutate_if(is.difftime, ~ round(.,0)) %>%
  as.data.frame()
int_FN_df
#View(int_FN_df)

# Buffer to 2 mins before and after each following period
buffer <- (60*2)/
  as.numeric(difftime(imm_metrics_FN$DateTime[2], imm_metrics_FN$DateTime[1], units = "secs"))

# Buffer for each bird, then recombine into 1 dataframe
imm_metrics_FNlist <- split(imm_metrics_FN, f = imm_metrics_FN$Bird)
imm_metrics_FNlist <- lapply(imm_metrics_FNlist, function(x){
  for(i in 1:nrow(x)){
    # Fill in vector for new column (1 if within 5 mins of FP, otherwise 0)
    if(!is.na(x$int_num[i])){
      start <- ifelse(i-buffer > 0, i-buffer, 1)
      end <- ifelse(i+buffer < nrow(x), i+buffer, nrow(x))
      x$buffer[start:end] <- 1
    } else {
      next
    }
  }
  return(x)
})
imm_metrics_FN <- do.call(rbind, imm_metrics_FNlist)

# Narrow down to segments with false negatives
FNcut <- imm_metrics_FN %>%
  mutate(tmp = ifelse(buffer == 1 & (lag(buffer) == 0 | is.na(lag(buffer))), 1, 0)) %>%
  filter(Bird %in% c("17B brown", "69G blue"), buffer == 1) %>%
  #mutate(tmp = ifelse(as.numeric(difftime(DateTime, lag(DateTime))) > 6 | is.na(lag(DateTime)), 1, 0)) %>%
  mutate(split = cumsum(tmp)) %>%
  dplyr::select(-tmp)
FNcut_list <- split(FNcut, f = FNcut$split)

# Join GPS data onto immersion data
FNcut_list <- pblapply(FNcut_list, function(df){
  
  # Which bird
  bird <- unique(df$Bird)
  if(length(bird) > 1){
    message(paste0("2 birds in one dataframe: ", paste(unique(df$Bird), collapse = ", ")))
    return(NULL)
  }
  
  # Get segment of GPS track from same bird and time period
  gps <- tracks.reg.list6[[which(sapply(tracks.reg.list6, function(x){x$Bird[1]}) == bird)]] %>%
    filter(DateTime >= min(df$DateTime) - minutes(1), 
           DateTime <= max(df$DateTime) + minutes(1))
  
  # Interpolate to exact same timesteps as immersion data
  gps_mv <- move(x = gps$Longitude,
                 y = gps$Latitude,
                 time = as.POSIXct(gps$DateTime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"),
                 proj = st_crs(4326)[[2]],
                 animal = gps$Bird,
                 removeDuplicatedTimestamps = TRUE)
  gps_interp <- interpolateTime(gps_mv,
                                time = df$DateTime,
                                spaceMethod = "greatcircle") %>%
    as.data.frame() %>%
    rename(Longitude = coords.x1, Latitude = coords.x2) %>%
    mutate(resolution = gps$resolution[1]) %>%
    dplyr::select(Latitude, Longitude, resolution)
  
  # Put Lat and Long columns into dataframe and return
  return(cbind(df, gps_interp))
  
})

# Read in the vessel data
AIS_Tracks <- read.csv("Output Dataset Files/BBA_gfw_vessels_redacted.csv", header = TRUE) %>%
  mutate(Bird = as.factor(Bird), vessel_class = as.factor(vessel_class), fishing = as.factor(fishing), SubTraj = as.factor(SubTraj)) %>%
  mutate(DateTime = ymd_hms(DateTime, tz = "UTC"))

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



# Function to produce BV interaction map with corrected immersion
# Takes in list of tracks 
BV.int.map <- function(Btrack){
  
  # Get first and last DateTime
  min_Bdt <- Btrack$DateTime[1]
  max_Bdt <- Btrack$DateTime[nrow(Btrack)]
  
  # Get vessel track associated with that interaction
  Vtrack <- AIS_Tracks %>% filter(Bird == Btrack$Bird[1])
  min_Vdt <- get.vdt(Btrack, Vtrack)[1]
  max_Vdt <- get.vdt(Btrack, Vtrack)[2]
  Vtrack <- Vtrack %>% filter(DateTime >= min_Vdt, DateTime <= max_Vdt)
  
  # Get boundary coordinates
  minLong <- min(Btrack$Longitude)
  maxLong <- max(Btrack$Longitude)
  rangeLong <- abs(maxLong - minLong)
  minLat <- min(Btrack$Latitude)
  maxLat <- max(Btrack$Latitude)
  rangeLat <- abs(maxLat - minLat)
  
  # Plot map
  map <- ggplot() + 
    geom_point(data = Vtrack, aes(x = Longitude, y = Latitude), size = 0.7) +
    geom_path(data = Vtrack, aes(x = Longitude, y = Latitude, group = SubTraj)) +
    geom_point(data = Btrack, aes(x = Longitude, y = Latitude, colour = leg_wet_C), size = 0.7) +
    geom_path(data = Btrack, aes(x = Longitude, y = Latitude, colour = leg_wet_C, group = 1)) +
    coord_sf(xlim = c(minLong, maxLong), ylim = c(minLat, maxLat), expand = TRUE) +
    labs(x = "Longitude", y = "Latitude", colour = "Corrected\nImmersion") +
    annotate(
      geom = "text", x = maxLong - 0.15*rangeLong, y = maxLat - 0.4*rangeLat,
      label = paste("Bird: ", Btrack$Bird[1], "\n",
                    "Resolution: ", Btrack$resolution[1], "\n",
                    "Int num: ", max(Btrack$int_num, na.rm = TRUE), "\n",
                    "Max regularity = ", round(max(Btrack$P_LWLD), digits = 2), "\n",
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
  
  # Return map
  return(map)
}

# Make some maps
maps <- pblapply(FNcut_list, BV.int.map)
maps[c(1,4,5,6)]


