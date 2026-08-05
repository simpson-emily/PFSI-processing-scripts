#### PFSI-processing-scripts

# Data is supplementary to: Simpson et al. 2026 Fire-Moisture Interactions Shape Recovery Trajectories Across Temperate Vegetation Communities

## Dataset overview

This dataset contains the data and scripts from a decade-long assessment of vegetation change following fire in the Greater Blue Mountains World Heritage Area, NSW, Australia. The dataset includes soil moisture and post-fire stability index values from 760 sampling sites stratified across five fire severity classes and seven vegetation formations. The R script and JavaScript used for processing and analysis are also provided.

## Files included

### 1_PFSI_Calculation.js

JavaScript for processing satellite images and calculating post-fire stability index values across sites and seasons.

### 2_Site_Sampling.js

JavaScript for site creation (stratification across fire severity and vegetation classes) and site sampling (both post-fire stability index and soil moisture values). 

### 3_Model_code.R

R script to generate final Bayesian model. 

### Final_Dataset.csv

Site-level measurements of vegetation change (PFSI) and soil moisture (SM) and site-level descriptors (vegetation community, fire severity). Output of 1_PFSI_Calculation and 2_Site_Sampling. Used to run 3_Model code.R.

#### Variable descriptions:
* sampleId: Unique sample id
* year: Sampling year
* season: Sampling season
* fireSeverity: Numeric fire severity class (0 - 5).
* fireSevN: Qualitative fire severity class (0 = Unburnt, 2 = Low: burnt understorey & unburnt canopy, 3 = Moderate: partial canopy scorch, 4 = High: complete canopy scorch & partial canopy consumption, 5 = Extreme: full canopy consumption)
* vegType: Numeric vegetation code (5, 6, 8, 9, 11, 34, 1516)
* vegTypeN: Qualitative vegetation class (5 = Forested Wetland, 6 = Freshwater Wetland, 8 = Grassy Woodland, 9 = Heathland, 11 = Rainforest, 34 = Dry Sclerophyll Forest, 1516 = Wet Sclerophyll Forest)
* PFSI = Post-fire Stability Index values
* SM = Soil Moisture values 
* SMLag1 = Soil moisture values lagged by one year
* SMLag1Seas = Soil moisture values lagged by one season

## Code / software

JavaScript was run in Google Earth Engine. Model creation was done in R version 4.5.0 (R Code Team 2025).

