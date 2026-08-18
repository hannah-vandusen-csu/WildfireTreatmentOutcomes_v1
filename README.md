# Automated Workflow for Wildfire-Treatment Outcomes: Version 1

This repository introduces a standardized, scalable workflow to quantify wildfire-treatment outcomes across diverse landscapes. This workflow supports the study "Automated workflow demonstrates mixed treatment impacts on burn severity across ecoregion scales" (Van Dusen et al. 2026). These scripts evaluate federal fuel treatments across 50 wildfires in southern Oregon and Northern California (2012–2021). Specifically, the workflow aims to predict the probability of high severity fire in treated and untreated areas across two distinct landscapes and environmental conditions: the productive, western Klamath mountains (west) and the more arid southeastern Cascades (east). More detailed information and methodology from the study can be found at Van Dusen et al. 2024 (doi: ). 

Broadly, this study and associated repository proposes a semi-automated workflow to understand the wildfire-treatment interactions. Although these scripts are applied to understand the association between fuel treatments and subsequent high severity wildfire, the authors developed them with the hope that the scripts can be applied to study any wildfire-treatment interaction with a measurable outcome. For example, treatment's effect on infrastructure loss, containment success or another ecological wildfire outcome. This is the first published version of the wildfire-treatment outcomes workflow (WTO-V1). Please contact the cooresponding author (h.vandusen@colostate.edu) if you are curious about further developments and want to apply these script to other outcomes or fires. 


> **Note:** Raw and processed data are **not** tracked in this Github repository. All data can be found in the dryad repository: ADDD!

---

## Workflow Overview

The analysis follows three sequential stages. Set `user_input/user_input_general.csv` and `user_input/treatment_selection.csv` before running any scripts.

![Analysis Workflow](https://github.com/user-attachments/assets/12e3a1fe-a7f5-4dc0-902e-4ee043e24282)

---

## Stage 1 — Data Preparation

Select fires, compile treatment histories, map predictor variables, and extract sample points.

| Script | Description |
|--------|-------------|
| `01_clean_treatments.R` | Download and clean FACTS/NFPORS treatment data; clip to fire perimeters; assign treatment categories |
| `02_climateNA.R` | Extract 30-year climate normals using ClimateNA |
| `03_dob_parks.R` | Compute day-of-burning (DOB) surface using Parks et al. (2013) method *(we cannot redistribute this script -- Please, contact Sean Parks for latest version)* |
| `04_upload_dob_gee_assets.ipynb` | Upload DOB raster assets to Google Earth Engine |
| `05_rasterize_geoprocess.R` | Rasterize treatment polygons; geoprocess upwind predictor layers |
| `06_get_previous_fire_severity.ipynb` | Retrieve pre-fire burn severity from GEE |
| `07_previous_fire_severity_expand.R` | Expand previous fire severity rasters to analysis extent |
| `08_create_treatment_and_fire_masks.R` | Build treatment presence/absence and previous-fire masks |
| `09_create_initial_severity_points.R` | Generate candidate sample-point grid from burn severity layer |
| `10_get_fire_severity.ipynb` | Extract dNBR / RdNBR burn-severity values from GEE for each sample point |
| `11_GEE_data_extraction.ipynb` | Extract fire weather, topography, climate, and spatial predictor variables from GEE |
| `12_R_data_extraction.R` | Supplement GEE extraction with R-based predictor layers |

---

## Stage 2 — Fire-level Matching

Determine the dominant response drivers within each individual fire, then match treated and untreated sample points on those drivers to create balanced comparison pairs.

| Script | Description |
|--------|-------------|
| `13_fire_clustering.R` | Cluster fires by climate similarity to ensure adequate post-matching sample sizes |
| `14_cluster_level_models.R` | Fit fire-level random forest models; rank and weight predictor variables by importance |
| `15_paired_matching.R` | Match treated to untreated points by minimizing weighted multivariate Euclidean distance |

---

## Stage 3 — Landscape-level Modeling

Pool matched pairs across fires, reduce correlated predictors, and fit multi-fire random forest models to evaluate treatment effects on each outcome.

| Script | Description |
|--------|-------------|
| `16_landscape_df_manipulation.R` | Assemble landscape-level dataframes from matched fire pairs|
| `17_landscape_variable_reduction.R` | Remove multicollinear predictors; apply recursive feature elimination; balance binary response variable |
| `18_final_landscape_level_models.R` | Fit landscape-level random forest classification models; evaluate treatment effects through model-based inference across treatment age and biophysical gradients |

---

## User Inputs

Configured in `user_input/` before running scripts:

| File | Purpose |
|------|---------|
| `user_input_general.csv` | Script and data directory paths, coordinate reference system *The user must change this to point to data location on their local machine*|
| `treatment_selection.csv` | Outcome-specific parameters for treatments: fire perimeter buffer, minimum treatment size, and treatment age window for burn severity, containment, and infrastructure outcomes |

---

## Function Scripts

Helper functions called by the numbered scripts are stored in `function_scripts/`, organized by stage:

```
function_scripts/
├── treatments/           # Treatment cleaning and classification helpers
├── rasterize_geoprocess/ # Rasterization and upwind geoprocessing helpers
├── gee/                  # GEE extraction helpers
├── R_extraction/         # R-based data extraction helpers
└── variable_reduction/   # Variable selection and model helpers
```

---

## Required Data (Box)

All data is publically available. Data specifically used for this project is located in the data a data responsitory. 
```

---

## Dependencies

Scripts require **R** (≥ 4.2) and **Python** (for GEE notebooks). Key R packages:

`terra` · `sf` · `tidyverse` · `tidymodels` · `ranger` · `vip` · `DALEXtra` · `data.table` · `doParallel` · `future.apply` · `ggplot2` · `patchwork` · `lubridate` · `tigris`

GEE notebooks require a Google Earth Engine account and the `earthengine-api` Python package.

---

## Contact

**Hannah Van Dusen** — h.vandusen@colostate.edu, hannah.vandusen@usda.gov  
Rocky Mountain Research Station, USDA Forest Service
