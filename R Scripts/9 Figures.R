######################################################################
############ FIGURES - IMMERSION PREDICTS VESSEL FOLLOWING ###########
######################################################################

# Produces the figures for this study.


############### CONTENTS #######################

# 1. Load packages and data
# 2. Plot important maps
# 3. Descriptive statistics
# 4. Immersion metrics figures
# 5. Plot high-res immersion metrics figure for poster
# 6. Random forest figures and tables





################## 1. LOAD PACKAGES AND DATA ##################################

# Packages
library(tidyverse)
library(lubridate)
library(pbapply)
library(parallel)
library(caret)
library(randomForest)
library(ggpubr)
theme_set(theme_bw())

# Set seed
set.seed(123)

# Set WD
setwd("C:/Users/jonat/OneDrive/Desktop/University of Oxford DPhil Biology/DPhil Chapters/Ch 1 BBA-vessel interactions/Data and Scripts/")

# Read in immersion metrics (corrected)
imm_metrics_list1 <- readRDS("Output Dataset Files/imm_metrics_list1.RData")
imm_metrics_df1 <- do.call(rbind, imm_metrics_list1) # All timesteps
imm_metrics_df2 <- imm_metrics_df1 %>%
  filter(imm_factor == "intermediate") %>% # Only foraging timesteps
  mutate(follow = ifelse(!is.na(follow), "F", "N")) %>%
  mutate(follow = relevel(as.factor(follow), ref = "N"))

# Read in immersion metrics (uncorrected - GLS)
imm_metrics_list1_gls <- readRDS("Output Dataset Files/imm_metrics_list1_gls.RData")
imm_metrics_gls1 <- do.call(rbind, imm_metrics_list1_gls) %>%
  filter(imm_factor == "intermediate") %>%
  mutate(follow = ifelse(!is.na(follow), "F", "N")) %>%
  mutate(follow = relevel(as.factor(follow), ref = "N"))

# Read in bouts metrics
imm_metrics_bouts1 <- readRDS("Output Dataset Files/imm_metrics_bouts1.RData")
imm_metrics_bouts1_gls <- readRDS("Output Dataset Files/imm_metrics_bouts1_gls.RData")









################ 2. PLOT IMPORTANT MAPS #######################################

# Mapping packages
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(sp)
library(sf)

# Read in the bird data
BBA_int_list <- readRDS("Output Dataset Files/BBA_Track_List_Imm_Int.RData")
BBA_Tracks <- read.csv("Output Dataset Files/BBA_All_Tracks_wImmersion.csv", header = TRUE) %>%
  mutate(DateTime = ymd_hms(DateTime))

# Read in the vessel data
AIS_Tracks <- read.csv("Output Dataset Files/BBA_gfw_vessels_redacted.csv", header = TRUE) %>%
  mutate(Bird = as.factor(Bird), vessel_class = as.factor(vessel_class), fishing = as.factor(fishing), SubTraj = as.factor(SubTraj)) %>%
  mutate(DateTime = ymd_hms(DateTime, tz = "UTC"))

# Read in interactions summmary
int_summary <- read.csv("Visual Interactions Analysis/BBA Interactions Summary.csv") %>%
  mutate(Real_start = dmy_hm(Real_start), Real_end = dmy_hm(Real_end)) %>%
  dplyr::select(Bird, Int_Num, Int_Class) %>%
  rename(interaction = Int_Num)

# Set CRS of input data: WGS 1984
WGS1984 <- st_crs(4326)[[2]]

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
int_coords <- left_join(int_coords, int_summary, by = join_by(Bird, interaction))
int_coords <- int_coords[int_coords$Int_Class != "Non-Interaction", ]

# Land polygons
samerica <- ne_countries(scale = "large", returnclass = "sf", continent = "South America")

# Figure for manuscript: Plot all GPS tracks of albatrosses, with vessel interactions included
ggplot(data = samerica) + 
  geom_sf(color = "black", fill = "grey83") +
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
vessels <- AIS_Tracks %>% filter(DateTime > ymd_hm("2017-12-26 03:50"), DateTime < ymd_hm("2017-12-26 05:30"))
track <- BBA_Tracks %>% filter(Bird == "75F red", DateTime > ymd_hm("2017-12-26 03:50"), DateTime < ymd_hm("2017-12-26 05:30"))
ggplot(data = samerica) + 
  geom_sf(color = "black", fill = "gray83") +
  geom_path(data = vessels, aes(x = Longitude, y = Latitude)) +
  geom_point(data = vessels, aes(x = Longitude, y = Latitude)) +
  geom_path(data = track, aes(x = Longitude, y = Latitude, colour = leg_wet, group = Bird)) +
  geom_point(data = track, aes(x = Longitude, y = Latitude, colour = leg_wet)) +
  coord_sf(xlim = c(-64.21, -64.16), ylim = c(-53.85, -53.825), expand = FALSE) +
  labs(x = "Longitude", y = "Latitude", colour = "Immersion\nState") +
  theme_bw() +
  theme(text = element_text(size = 14))

# Figure for poster: 75F red with vessel, zoomed out a bit and high res
#png(file = "High res figures/Poster 75Fred.png", width = 2800, height = 2100, units = "px")
ggplot(data = samerica) + 
  geom_sf(color = "black", fill = "gray83") +
  coord_sf(xlim = c(-64.22, -64.15), ylim = c(-53.855, -53.82), expand = FALSE) +
  geom_path(data = vessels, aes(x = Longitude, y = Latitude), linewidth = 3.5) +
  geom_point(data = vessels, aes(x = Longitude, y = Latitude), size = 11.3) +
  geom_path(data = track, aes(x = Longitude, y = Latitude, colour = leg_wet, group = Bird), size = 3.5) +
  geom_point(data = track, aes(x = Longitude, y = Latitude, colour = leg_wet), size = 11.3) +
  labs(x = "Longitude", y = "Latitude", colour = "Immersion\nState") +
  theme_bw() +
  theme(text = element_text(size = 40),
        panel.grid = element_blank())
#dev.off()











##################### 3. DESCRIPTIVE STATISTICS ################################

# Read in the bird data
BBA_int_list <- readRDS("Output Dataset Files/BBA_Track_List_Imm_Int.RData")
BBA_Tracks <- read.csv("Output Dataset Files/BBA_All_Tracks_wImmersion.csv", header = TRUE) %>%
  mutate(DateTime = ymd_hms(DateTime))

# Read in the vessel data
AIS_Tracks <- read.csv("Output Dataset Files/BBA_gfw_vessels_redacted.csv", header = TRUE) %>%
  mutate(Bird = as.factor(Bird), vessel_class = as.factor(vessel_class), fishing = as.factor(fishing), SubTraj = as.factor(SubTraj)) %>%
  mutate(DateTime = ymd_hms(DateTime, tz = "UTC"))

# Read in interactions summmary
int_summary <- read.csv("Visual Interactions Analysis/BBA Interactions Summary.csv") %>%
  mutate(Real_start = dmy_hm(Real_start), Real_end = dmy_hm(Real_end)) %>%
  rename(interaction = Int_Num)

# Descriptive stats on foraging trips
imm_metrics_df1 %>%
  group_by(Bird) %>%
  summarise(track_length = difftime(max(DateTime), min(DateTime), units = "days")) %>%
  as.data.frame() %>%
  summarise(mean = mean(track_length),
            min = min(track_length),
            max = max(track_length))

# From visual analysis:
# Birds that stayed in Falklands EEZ: 07E white, 09G blue, 40T blue, unringed
# Birds where battery died before return to the colony: 20L yellow, 22V red, 26F red, 38T blue, 47F red, 52B brown, 59T blue, 60G blue, 77T blue, 78B brown, 86G blue

# Descriptive stats on interactions
A <- int_summary[, c(1,2,3,5,6,7)] %>% 
  filter(Res <= 10 | Bird %in% c("92F red", "38F red", "00L yellow", "47A red"), 
         Int_Class != "Non-Interaction") 
A %>% group_by(Int_Class) %>% 
  summarise(N = n(), prop = n()/nrow(A))
A %>% group_by(Bird) %>% summarise(N = n())
A %>% filter(Int_Class == "Follow") %>% group_by(Bird) %>% summarise(N = n())
# 4:45-22:05 daytime
A_local <- A %>%
  mutate(Local_start = with_tz(Real_start, tzone = "America/Argentina/Buenos_Aires"),
         Local_end = with_tz(Real_end, tzone = "America/Argentina/Buenos_Aires")) %>%
  mutate(night = ifelse(hm(format(as.POSIXct(Local_start), format = "%H:%M")) < hm("04:45") |
                          hm(format(as.POSIXct(Local_start), format = "%H:%M")) > hm("22:05") |
                          hm(format(as.POSIXct(Local_end), format = "%H:%M")) < hm("04:45") |
                          hm(format(as.POSIXct(Local_end), format = "%H:%M")) > hm("22:05"),
                        "dark", "")) %>%
  mutate(elapsed = Local_end - Local_start)
A_local
mean(as.numeric(A_local$elapsed, units = "mins"), na.rm = TRUE)
sd(as.numeric(A_local$elapsed, units = "mins"), na.rm = TRUE)
min(as.numeric(A_local$elapsed, units = "mins"), na.rm = TRUE)
max(as.numeric(A_local$elapsed, units = "mins"), na.rm = TRUE)

# Following events duration
A_follow <- A_local %>% filter(Int_Class == "Follow")
mean(as.numeric(A_follow$elapsed, units = "mins"), na.rm = TRUE)
sd(as.numeric(A_follow$elapsed, units = "mins"), na.rm = TRUE)
min(as.numeric(A_follow$elapsed, units = "mins"), na.rm = TRUE)
max(as.numeric(A_follow$elapsed, units = "mins"), na.rm = TRUE)

# Get duration of wet and dry periods in minutes (just within foraging periods)
imm_metrics_df3 <- imm_metrics_df2 %>% mutate(duration = 6) # This is in seconds
for(i in levels(as.factor(imm_metrics_df3$Bird))){
  df_i <- imm_metrics_df3 %>% filter(Bird == i)
  
  if(length(df_i) > 0){
    for(j in levels(as.factor(df_i$follow))){
      df_j <- df_i %>% filter(follow == j)
      
      if(nrow(df_j) > 1){
        for(k in 2:nrow(df_j)){
          
          if(df_j$leg_wet_C[k] == df_j$leg_wet_C[k-1]){
            df_j$duration[k] <- df_j$duration[k-1] + 6
          }
        }
      }
      imm_metrics_df3[imm_metrics_df3$Bird == i & imm_metrics_df3$follow == j, ] <- df_j
    }
  }
}
imm_metrics_df4 <- imm_metrics_df3 %>%
  mutate(duration1 = ifelse(duration >= lead(duration), duration, NA)) %>%
  filter(!is.na(duration1)) %>%
  filter(follow == "F")
imm_metrics_df4 %>%
  group_by(leg_wet_C) %>%
  summarise(min_dur = min(duration1)/60, # Convert to minutes
            max_dur = max(duration1)/60, 
            mean_dur = mean(duration1)/60,
            sd_dur = sd(duration1)/60)
summary(imm_metrics_df4$duration1[imm_metrics_df4$leg_wet_C == "wet"])/60
summary(imm_metrics_df4$duration1[imm_metrics_df4$leg_wet_C == "dry"])/60

# Following bouts metrics
imm_metrics_bouts1 %>%
  filter(follow == "F") %>%
  summarise(mean_wet = mean(mean_wet)/60,
            mean_dry = mean(mean_dry)/60)

# Number of landings
summary(imm_metrics_df2$N_landings)
boxplot(imm_metrics_df2$N_landings)
N_landings_df <- imm_metrics_df2 %>%
  group_by(Bird, follow) %>%
  summarise(N_landings = mean(N_landings))
  #summarise(N_landings = sum(leg_wet_C == "wet" & lag(leg_wet_C) == "dry", na.rm = TRUE))
print(N_landings_df, n=23)
summary(N_landings_df$N_landings[N_landings_df$follow == "F"])
sd(N_landings_df$N_landings[N_landings_df$follow == "F"])

# What proportion of all time steps had nonzero PLWLD?
sum(imm_metrics_df1$P_LWLD > 0)/nrow(imm_metrics_df1)











####################### 4. IMMERSION METRICS FIGURES #############################

# Read in interactions summmary
int_summary <- read.csv("Visual Interactions Analysis/BBA Interactions Summary.csv") %>%
  mutate(Real_start = dmy_hm(Real_start), Real_end = dmy_hm(Real_end)) %>%
  dplyr::select(Bird, Int_Num, Real_start, Real_end, Int_Class) %>%
  rename(interaction = Int_Num)

# Read in immersion metrics data
imm_metrics_list1 <- readRDS("Output Dataset Files/imm_metrics_list1.RData")
imm_metrics_df1 <- do.call(rbind, imm_metrics_list1)
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

# Get dataframes of bouts
imm_metrics_bouts1 <- readRDS("Output Dataset Files/imm_metrics_bouts1.RData")
imm_metrics_bouts1_gls <- readRDS("Output Dataset Files/imm_metrics_bouts1_gls.RData")

# Plot time series of regularity for 75F red
imm_metrics_gls1 %>%
  filter(Bird == "75F red",
         DateTime > ymd_hm("2017-12-26 01:00"),
         DateTime < ymd_hm("2017-12-26 23:00")) %>%
  ggplot() +
  geom_line(aes(x = DateTime, y = P_LWLD_gls)) +
  geom_rect(data = int_summary[int_summary$Bird == "75F red", ],
            aes(xmin = Real_start,
                xmax = Real_end,
                ymin = -Inf,
                ymax = Inf),
            alpha = .2, fill = "#253582FF") +
  labs(y = "Regularity of Immersion") +
  theme_bw() +
  theme(text = element_text(size = 13))

# Plot histograms of predictors (time-based), split by following or not
p_mw <- ggplot(imm_metrics_df2, aes(x = mean_wet/60, fill = follow, colour = follow)) +
  geom_density(alpha = 0.3) +
  labs(title = "Mean duration of wet periods (min)",
       colour = "Foraging behaviour: ",
       fill = "Foraging behaviour: ") +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.title = element_blank(),
        text = element_text(size = 12)) +
  scale_fill_manual(values = c("#F9A242FF", "#253582FF"),
                    aesthetics = c("colour", "fill"),
                    labels = c("Natural foraging", "Vessel following"))
p_md <- ggplot(imm_metrics_df2, aes(x = mean_dry/60, fill = follow, colour = follow)) +
  geom_density(alpha = 0.3) +
  labs(title = "Mean duration of dry periods (min)",
       colour = "Foraging behaviour: ",
       fill = "Foraging behaviour: ") +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.title = element_blank(),
        text = element_text(size = 12)) +
  scale_fill_manual(values = c("#F9A242FF", "#253582FF"),
                    aesthetics = c("colour", "fill"),
                    labels = c("Natural foraging", "Vessel following"))
p_nl <- ggplot(imm_metrics_df2, aes(x = N_landings, fill = follow, colour = follow)) +
  geom_density(alpha = 0.3, adjust = 8) +
  labs(title = "Number of landings",
       colour = "Foraging behaviour: ",
       fill = "Foraging behaviour: ") +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.title = element_blank(),
        text = element_text(size = 12)) +
  scale_fill_manual(values = c("#F9A242FF", "#253582FF"),
                    aesthetics = c("colour", "fill"),
                    labels = c("Natural foraging", "Vessel following"))
p_pw <- ggplot(imm_metrics_df2, aes(x = P_wet, fill = follow, colour = follow)) +
  geom_density(alpha = 0.3) +
  labs(title = "Proportion of time spent wet",
       colour = "Foraging behaviour: ",
       fill = "Foraging behaviour: ") +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.title = element_blank(),
        text = element_text(size = 12)) +
  scale_fill_manual(values = c("#F9A242FF", "#253582FF"),
                    aesthetics = c("colour", "fill"),
                    labels = c("Natural foraging", "Vessel following"))
p_reg <- ggplot(imm_metrics_df2, aes(x = P_LWLD, fill = follow, colour = follow)) +
  geom_density(alpha = 0.3) +
  labs(title = "Regularity",
       colour = "Foraging behaviour: ",
       fill = "Foraging behaviour: ") +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.title = element_blank(),
        text = element_text(size = 12)) +
  scale_fill_manual(values = c("#F9A242FF", "#253582FF"),
                    aesthetics = c("colour", "fill"),
                    labels = c("Natural foraging", "Vessel following"))
hists_time <- ggarrange(p_mw, p_md, p_pw, p_nl, p_reg, ncol = 1, common.legend = TRUE, legend = "bottom") %>%
  annotate_figure(left = text_grob("Density", rot = 90, size = 18),
                  top = text_grob("Time-based metrics", size = 24))


# Plot histograms of predictors (bout-based), split by following or not
pb_mw <- ggplot(imm_metrics_bouts1, aes(x = mean_wet/60, fill = follow, colour = follow)) +
  geom_density(alpha = 0.3) +
  labs(title = "Mean duration of wet periods (min)",
       colour = "Foraging behaviour: ",
       fill = "Foraging behaviour: ") +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.title = element_blank(),
        text = element_text(size = 12)) +
  scale_fill_manual(values = c("#F9A242FF", "#253582FF"),
                    aesthetics = c("colour", "fill"),
                    labels = c("Natural foraging", "Vessel following"))
pb_md <- ggplot(imm_metrics_bouts1, aes(x = mean_dry/60, fill = follow, colour = follow)) +
  geom_density(alpha = 0.3) +
  labs(title = "Mean duration of dry periods (min)",
       colour = "Foraging behaviour: ",
       fill = "Foraging behaviour: ") +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.title = element_blank(),
        text = element_text(size = 12)) +
  scale_fill_manual(values = c("#F9A242FF", "#253582FF"),
                    aesthetics = c("colour", "fill"),
                    labels = c("Natural foraging", "Vessel following"))
pb_nl <- ggplot(imm_metrics_bouts1, aes(x = N_landings, fill = follow, colour = follow)) +
  geom_density(alpha = 0.3, adjust = 4) +
  labs(title = "Number of landings",
       colour = "Foraging behaviour: ",
       fill = "Foraging behaviour: ") +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.title = element_blank(),
        text = element_text(size = 12)) +
  scale_fill_manual(values = c("#F9A242FF", "#253582FF"),
                    aesthetics = c("colour", "fill"),
                    labels = c("Natural foraging", "Vessel following"))
pb_pw <- ggplot(imm_metrics_bouts1, aes(x = P_wet, fill = follow, colour = follow)) +
  geom_density(alpha = 0.3) +
  labs(title = "Proportion of time spent wet",
       colour = "Foraging behaviour: ",
       fill = "Foraging behaviour: ") +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.title = element_blank(),
        text = element_text(size = 12)) +
  scale_fill_manual(values = c("#F9A242FF", "#253582FF"),
                    aesthetics = c("colour", "fill"),
                    labels = c("Natural foraging", "Vessel following"))
pb_reg <- ggplot(imm_metrics_bouts1, aes(x = P_LWLD, fill = follow, colour = follow)) +
  geom_density(alpha = 0.3) +
  labs(title = "Max regularity",
       colour = "Foraging behaviour: ",
       fill = "Foraging behaviour: ") +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.title = element_blank(),
        text = element_text(size = 12)) +
  scale_fill_manual(values = c("#F9A242FF", "#253582FF"),
                    aesthetics = c("colour", "fill"),
                    labels = c("Natural foraging", "Vessel following"))
hists_bout <- ggarrange(pb_mw, pb_md, pb_pw, pb_nl, pb_reg, ncol = 1, common.legend = TRUE, legend = "bottom") %>%
  annotate_figure(top = text_grob("Bout-based metrics", size = 24))

# Plot both time-based and bout-based density historgrams all together
ggarrange(hists_time, hists_bout, nrow = 1, ncol = 2, common.legend = TRUE, legend = "bottom")





##################### 5. Plot high-res immersion metrics figure for poster ##############################

p_mw <- ggplot(imm_metrics_df2, aes(x = mean_wet/60, fill = follow, colour = follow)) +
  geom_density(alpha = 0.3) +
  labs(colour = "Foraging behaviour: ",
       fill = "Foraging behaviour: ") +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        text = element_text(size = 16)) +
  scale_fill_manual(values = c("#F9A242FF", "#253582FF"),
                    aesthetics = c("colour", "fill"),
                    labels = c("Not following", "Following"))
p_md <- ggplot(imm_metrics_df2, aes(x = mean_dry/60, fill = follow, colour = follow)) +
  geom_density(alpha = 0.3) +
  labs(colour = "Foraging behaviour: ",
       fill = "Foraging behaviour: ") +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        text = element_text(size = 16)) +
  scale_fill_manual(values = c("#F9A242FF", "#253582FF"),
                    aesthetics = c("colour", "fill"),
                    labels = c("Not following", "Following"))
p_pw <- ggplot(imm_metrics_df2, aes(x = P_wet, fill = follow, colour = follow)) +
  geom_density(alpha = 0.3) +
  labs(colour = "Foraging behaviour: ",
       fill = "Foraging behaviour: ") +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        text = element_text(size = 16)) +
  scale_fill_manual(values = c("#F9A242FF", "#253582FF"),
                    aesthetics = c("colour", "fill"),
                    labels = c("Not following", "Following"))
p_reg <- ggplot(imm_metrics_df2, aes(x = P_LWLD, fill = follow, colour = follow)) +
  geom_density(alpha = 0.3) +
  labs(colour = "Foraging behaviour: ",
       fill = "Foraging behaviour: ") +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        text = element_text(size = 16)) +
  scale_fill_manual(values = c("#F9A242FF", "#253582FF"),
                    aesthetics = c("colour", "fill"),
                    labels = c("Not following", "Following"))


hists_time_HR <- ggarrange(p_mw, p_md, p_pw, p_reg, ncol = 1, common.legend = TRUE, legend = "bottom")

# ggsave("High res figures/Poster density plots.png", plot = hists_time_HR,
#        width = 2400, height = 2800, units = "px")






######################### 6. RANDOM FOREST FIGURES AND TABLES #####################################

# Set FP threshold
fp_threshold <- 0.1

# Custom summary function (from 7 Random forest time.R)
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


# Read in random forest models
rf1 <- readRDS("Output Dataset Files/rf1.RData")
rf2 <- readRDS("Output Dataset Files/rf2.RData")
rf3 <- readRDS("Output Dataset Files/rf3.RData")
rf4 <- readRDS("Output Dataset Files/rf4.RData")
rf5 <- readRDS("Output Dataset Files/rf5.RData")
rf6 <- readRDS("Output Dataset Files/rf6.RData")
rf7 <- readRDS("Output Dataset Files/rf7.RData")
rf8 <- readRDS("Output Dataset Files/rf8.RData")
rf_list <- list(rf1, rf2, rf3, rf4, rf5, rf6, rf7, rf8)


# Confusion matrix for rf2
cm2 <- confusionMatrix(data = rf2$pred$pred, 
                       reference = rf2$pred$obs,
                       positive = "F", mode = "everything")

# Function to get N from each RF model
get.rf.n <- function(rf){nrow(rf$pred)}
rf_summaries_n <- pblapply(rf_list, get.rf.n)

# Function to get dataframe of results from each fold of RF model
get.rf.summary <- function(rf){
  summaries <- lapply(unique(rf$pred$Resample), function(resample){
    custom_summary(rf$pred[rf$pred$Resample == resample,])
  })
  summaries_df <- as.data.frame(do.call(rbind, summaries))
  return(summaries_df)
}
rf_summaries_custom <- pblapply(rf_list, get.rf.summary)

# Function to get more standard default summary metrics (accuracy, kappa) from each fold
get.rf.summary.default <- function(rf){
  summaries <- lapply(unique(rf$pred$Resample), function(resample){
    defaultSummary(rf$pred[rf$pred$Resample == resample,])
  })
  summaries_df <- as.data.frame(do.call(rbind, summaries))
  return(summaries_df)
}
rf_summaries_default <- pblapply(rf_list, get.rf.summary.default)

# Function to get F1 statistic from each fold
get.rf.summary.f1 <- function(rf){
  summaries <- lapply(unique(rf$pred$Resample), function(resample){
    prSummary(rf$pred[rf$pred$Resample == resample,], lev = c("F", "N"))
  })
  summaries_df <- as.data.frame(do.call(rbind, summaries))
  return(summaries_df)
}
rf_summaries_f1 <- pblapply(rf_list, get.rf.summary.f1)


# Combine summaries into a table for reporting
rf_summaries <- lapply(1:8, function(x){cbind(rf_summaries_n[[x]],
                                              rf_summaries_custom[[x]], 
                                              rf_summaries_default[[x]], 
                                              rf_summaries_f1[[x]])})
rf_mean <- round(do.call(rbind, lapply(1:8, function(x){apply(rf_summaries[[x]], 2, mean)})), 3)
rf_sd <- round(do.call(rbind, lapply(1:8, function(x){apply(rf_summaries[[x]], 2, sd)})), 3)
rf_min <- round(do.call(rbind, lapply(1:8, function(x){apply(rf_summaries[[x]], 2, min)})), 3)
rf_max <- round(do.call(rbind, lapply(1:8, function(x){apply(rf_summaries[[x]], 2, max)})), 3)
rf_results <- lapply(1:8, function(x){
  c(n = rf_mean[x,1],
    sensitivity = paste0(rf_mean[x, "sensitivity"], " (", rf_min[x, "sensitivity"], ", ", rf_max[x, "sensitivity"], ")"),
    fpp = paste0(rf_mean[x, "false_pos_percent"], " (", rf_min[x, "false_pos_percent"], ", ", rf_max[x, "false_pos_percent"], ")"),
    specificity = paste0(rf_mean[x, "specificity"], " (", rf_min[x, "specificity"], ", ", rf_max[x, "specificity"], ")"),
    accuracy = paste0(rf_mean[x, "Accuracy"], " (", rf_min[x, "Accuracy"], ", ", rf_max[x, "Accuracy"], ")"),
    kappa = paste0(rf_mean[x, "Kappa"], " (", rf_min[x, "Kappa"], ", ", rf_max[x, "Kappa"], ")"),
    F = paste0(rf_mean[x, "F"], " (", rf_min[x, "F"], ", ", rf_max[x, "F"], ")"))
})
rf_results_df <- do.call(rbind, rf_results) %>%
  as.data.frame() %>%
  mutate(model = 1:8,
         data = c(rep("Time", 4), rep("Bout", 4)),
         immersion = rep(c("Corrected", "Corrected", "Uncorrected", "Uncorrected"), 2),
         predictors = rep(c("Regularity", "All"), 4)) %>%
  rename(N = `n.rf_summaries_n[[x]]`) %>%
  relocate(model, data, immersion, predictors)
rf_results_df
#write.csv(rf_results_df, "Output Dataset Files/rf_results_df.csv")


# Variable importance plots
varImp2 <- varImp(rf2)$importance
varImp4 <- varImp(rf4)$importance
varImp6 <- varImp(rf6)$importance
varImp8 <- varImp(rf8)$importance
varImp_all <- list(varImp2, varImp4, varImp6, varImp8)
vip2 <- ggplot(varImp2, aes(x = Overall, y = reorder(rownames(varImp2), Overall))) + 
  geom_col() +
  scale_y_discrete(labels = c("Number of landings", "Proportion wet", "Mean dry", "Mean wet", "Regularity")) +
  labs(x = "Importance (scaled mean decrease in accuracy)", y = "Predictor",
       title = "2. Time-based RF, Corrected Immersion") +
  theme(panel.grid.major.y = element_blank(),
        axis.title = element_blank())
vip4 <- ggplot(varImp4, aes(x = Overall, y = reorder(gsub("_gls", "", rownames(varImp4)), Overall))) + 
  geom_col() +
  scale_y_discrete(labels = c("Number of landings", "Proportion wet", "Mean dry", "Mean wet", "Regularity")) +
  labs(x = "Importance (scaled mean decrease in accuracy)", y = "Predictor",
       title = "4. Time-based RF, Uncorrected Immersion") +
  theme(panel.grid.major.y = element_blank(),
        axis.title = element_blank())
vip6 <- ggplot(varImp6, aes(x = Overall, y = reorder(rownames(varImp6), c(5,2,4,3,1)))) + 
  geom_col() +
  scale_y_discrete(labels = c("Number of landings", "Proportion wet", "Mean dry", "Mean wet", "Regularity")) +
  labs(x = "Importance (scaled mean decrease in accuracy)", y = "Predictor",
       title = "6. Bout-based RF, Corrected Immersion") +
  theme(panel.grid.major.y = element_blank(),
        axis.title = element_blank())
vip8 <- ggplot(varImp8, aes(x = Overall, y = reorder(rownames(varImp8), c(5,2,4,3,1)))) + 
  geom_col() +
  scale_y_discrete(labels = c("Number of landings", "Proportion wet", "Mean dry", "Mean wet", "Regularity")) +
  labs(x = "Importance (scaled mean decrease in accuracy)", y = "Predictor",
       title = "8. Bout-based RF, Uncorrected Immersion") +
  theme(panel.grid.major.y = element_blank(),
        axis.title = element_blank())
ggarrange(vip2, vip4, vip6, vip8, ncol = 1) %>%
  annotate_figure(left = text_grob("Predictors", rot = 90, size = 18), 
                  bottom = text_grob("Importance (scaled mean decrease in accuracy)"))





# Plot RF prediction probabilities compared to vessel following for all birds
birds <- sapply(imm_metrics_list1, function(x){x$Bird[1]})
imm_metrics_df2.1 <- imm_metrics_df2[!is.na(imm_metrics_df2$mean_wet),]
imm_metrics_df2.1 <- imm_metrics_df2.1[!is.na(imm_metrics_df2.1$mean_wet),] %>%
  mutate(rf2_probs = rf2$pred$F[order(rf2$pred$rowIndex)])
ggplot() +
  geom_vline(data = imm_metrics_df1[imm_metrics_df1$Bird %in% birds[1:2] & !is.na(imm_metrics_df1$follow), ],
             aes(xintercept = DateTime), col = "lightgreen", alpha = 0.3) +
  geom_line(data = imm_metrics_df2.1[imm_metrics_df2.1$Bird == birds[1:2], ],
            aes(x = DateTime, y = rf2_probs)) +
  facet_wrap(~ Bird, scales = "free_x", ncol = 1) +
  labs(y = "RF Probability")
# ggplot() +
#   geom_vline(data = imm_metrics_df1[imm_metrics_df1$Bird == birds[1:9] & !is.na(imm_metrics_df1$follow), ],
#              aes(xintercept = DateTime), col = "lightgreen", alpha = 0.3) +
#   geom_line(data = imm_metrics_df2.1[imm_metrics_df2.1$Bird == birds[1:9], ],
#             aes(x = DateTime, y = rf2_probs)) +
#   facet_wrap(~Bird, scales = "free_x", ncol = 1) +
#   labs(y = "RF Probability")
# ggplot() +
#   geom_vline(data = imm_metrics_df1[imm_metrics_df1$Bird == birds[10:18] & !is.na(imm_metrics_df1$follow), ],
#              aes(xintercept = DateTime), col = "lightgreen", alpha = 0.3) +
#   geom_line(data = imm_metrics_df2.1[imm_metrics_df2.1$Bird == birds[10:18], ],
#             aes(x = DateTime, y = rf2_probs)) +
#   facet_wrap(~Bird, scales = "free_x", ncol = 1) +
#   labs(y = "RF Probability")
# ggplot() +
#   geom_vline(data = imm_metrics_df1[imm_metrics_df1$Bird == birds[19:27] & !is.na(imm_metrics_df1$follow), ],
#              aes(xintercept = DateTime), col = "lightgreen", alpha = 0.3) +
#   geom_line(data = imm_metrics_df2.1[imm_metrics_df2.1$Bird == birds[19:27], ],
#             aes(x = DateTime, y = rf2_probs)) +
#   facet_wrap(~Bird, scales = "free_x", ncol = 1) +
#   labs(y = "RF Probability")
# ggplot() +
#   geom_vline(data = imm_metrics_df1[imm_metrics_df1$Bird == birds[28:36] & !is.na(imm_metrics_df1$follow), ],
#              aes(xintercept = DateTime), col = "lightgreen", alpha = 0.3) +
#   geom_line(data = imm_metrics_df2.1[imm_metrics_df2.1$Bird == birds[28:36], ],
#             aes(x = DateTime, y = rf2_probs)) +
#   facet_wrap(~Bird, scales = "free_x", ncol = 1) +
#   labs(y = "RF Probability")
# ggplot() +
#   geom_vline(data = imm_metrics_df1[imm_metrics_df1$Bird == birds[37:45] & !is.na(imm_metrics_df1$follow), ],
#              aes(xintercept = DateTime), col = "lightgreen", alpha = 0.3) +
#   geom_line(data = imm_metrics_df2.1[imm_metrics_df2.1$Bird == birds[37:45], ],
#             aes(x = DateTime, y = rf2_probs)) +
#   facet_wrap(~Bird, scales = "free_x", ncol = 1) +
#   ylim(0, 1) +
#   labs(y = "RF Probability")
# 

