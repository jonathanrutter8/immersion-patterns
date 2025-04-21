# immersion-patterns
R scripts and sample data from Rutter et al. 2025: Immersion patterns alone can predict vessel following by albatrosses

# Contents
- Folders labelled as containing "animations" or "maps" contain output animated and static maps from the FULL dataset, not just the demonstration dataset. "BV Interactions Animations" is the most important for understanding our study's findings.
- "R scripts" contains 10 R scripts for this study's analysis, described below.
- "Immersion CSVs" and "Indiv GPS Tracks - CSV" contain a subset of the raw data (from 2 birds) used for this study, for demonstration of the scripts. The full dataset can be obtained via Movebank, the BirdLife Seabird Tracking Database, and Global Fishing Watch.
- "Output Dataset Files" contains intermediate files produced by the R scripts, using the demonstration dataset. The content of this folder allows each R script to be run independently. For analysis of other datasets, scripts should be run in order.
- "Visual Interactions Analysis" contains CSVs produced through manual annotation of the FULL dataset, not just the demonstration dataset. For example, "BBA Interactions Summary.csv" contains start and end times of the interactions based on animations.
- "SUMMARY OF DEPLOYMENTS CSV.csv" contains the metadata for the FULL dataset.

# R Scripts
All data processing and analysis for this study were performed using R 4.4.1 (R Core Team, 2024).
1.	GPS read in: Reads in black-browed albatross GPS tracks (2 full tracks provided).
2.	GLS read in: Reads in GLS-immersion data (2 full .deg files provided) and joins it to the GPS tracks.
3.	BV overlap: Takes in seabird GPS tracks to find daily overlap with fishing effort using the Global Fishing Watch (GFW) API. We used this script to facilitate a data query to Global Fishing Watch for individual vessel AIS tracks.
4.	BBA BV interaction: Overlaps bird GPS tracks with vessel AIS tracks to identify interactions within a distance threshold. Produces maps and animations of interactions. Note that vessel AIS tracks from GFW must be downloaded manually from the GFW website, or queried from the organisation. We have provided 2 simulated tracks for demonstration purposes.
5.	GLS immersion correction: Inspects geolocator (GLS) and TDR-immersion data using maps. Uses 2- and 3-state Hidden Markov Models (HMMs) to derive immersion data from GPS tracks. Produces a correct immersion dataset by integrating GPS and GLS-immersion data. Produces several diagnostic figures to inspect GLS-immersion error.
6.	Immersion metrics: Calculates several rolling metrics, including immersion regularity, from GLS-immersion data. Tests the predictive performance of the regularity metric across parameter space using a binary threshold classifier model. Implements custom k-fold cross-validation to assess model performance (for full dataset, k = 5; for demonstration dataset, k = 2).
7.	Random forest time: Fits random forest models to predict vessel following from immersion metrics on a timestep-by-timestep basis. Implements custom k-fold cross-validation to assess model performance (for full dataset, k = 5; for demonstration dataset, k = 2). Demonstration dataset will not produce same results.
8.	Random forest bouts: Equivalent to “Random forest time”, but works on a bout-by-bout basis.
9.	Figures: Produces the figures for this study.
10.	Extra figures: Produces the figures for the Supporting Information of this study.
