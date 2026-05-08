####################################.
###
### Script for Treatment Cleaning and Fire Selection
### Uses Raw FACTS and NFPORS data
### Hannah Van Dusen --- 12/01/2023
### hannah.vandusen@usda.gov
### Adapted for FACTS and NFPORS by BB (elizabeth.black@usda.gov) 3/2024
### re-organized/edited by Mark Kreider 2/19/2026
###
####################################.

# Script set up 
# load packages -----
# Necessary packages
pkgs <- c("terra",
          "lubridate",
          "tidyterra",
          "tidyverse",
          "data.table",
          "doParallel",
          "future.apply",
          "crayon",
          "tigris", 
          "rstudioapi")

# Install any that are missing
missing <- pkgs[!(pkgs %in% installed.packages()[,"Package"])]
if (length(missing) > 0) {install.packages(missing)}

# Load all packages
lapply(pkgs, library, character.only = TRUE)

# read in user inputs -----
# set working directory
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# Read in user input
user_input_csv <- read.csv("user_input/user_input_general.csv", stringsAsFactors = FALSE)

# loop through user input csv and assign the values provided to the named objects
for (i in seq_len(nrow(user_input_csv))){
  val <- user_input_csv$value[i]
  if(user_input_csv$treat_as_R_expression[i] == T){val <- eval(parse(text = val))}
  assign(user_input_csv$name[i], val, envir = .GlobalEnv)
  rm(val)
}

# load required functions from other scripts
function_scripts <- list.files(file.path(script_directory, "function_scripts/treatments"), full.names = TRUE) 
sapply(function_scripts, source)

# input data directory defined by above csv file, stored in object "global_data_directory"
# data directory for processed data defined by above csv file, stored in object "processed_data_directory"
# CRS defined by above csv file, stored in object "CRS"

# Step 0: Fire Selection -------------------------------------------------------

#### Actions Taken in this Script:
#### 1) Run Fire selection script: Follow these prompts
####
#### A) Ask for project AOI - 2 options: A) write 2-character state codes (ex. WA, OR)
#### or B) Wildfire Crisis Strategy Landscape Name(s)
####
#### B) Ask for time period of interest minimum and maximum year of interest (ex. 2021 - 2021, only 2021 fires)
####
#### C) Ask if there is pre-selected csv file with Fire names listed in a column titles "Names" ex.(fires.list.csv)
#### fire csv list must match MTBS ID, or Script will select all available fires with treatment data


# Set Source Script and Run Function for Fire Selection
fire_selection(global_data_directory, processed_data_directory, CRS)

# Step 1: Automated Treatment download-----------------------------
# This script extracts treatments that occurred w/in fire perimeters and excludes fires with no treatment overlap
treatment_download(global_data_directory, processed_data_directory)

# Step 2: Script For Automated Treatment cleaning
# This script gets rid of non-fuels activities-------------------------------
treatment_cleaning(global_data_directory,  processed_data_directory)


# Step 3: Script For Spatial join of csv and shapefile for treatment
treatment_spatialjoin(global_data_directory, processed_data_directory , the_prj = CRS)

# Run Step 4 script for Rasterize Treatment Vectors
# This rasterizes classed treatments and condenses instances of multiple treatments on a pixel----------------------
treatment_to_raster(global_data_directory, processed_data_directory , the_prj = CRS)

# Run Step 5 script for Rasterize Treatment Vectors
# This creates a separate raster for each outcome - based on outcome specific thresholds
treatments_by_outcome(global_data_directory, 
                     processed_data_directory , 
                     script_directory,
                     the_prj = CRS, 
                     wto_outcome = 'burn_severity')


print('Treatment Cleaning Complete')




