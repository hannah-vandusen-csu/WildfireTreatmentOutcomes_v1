#####################################.
##  Script for variable reduction  ##.
#####################################.

# load packages -----

# List of required packages
packages <- c(
  "sf", "tidyverse", "spdep", "spmoran", "terra", "tidymodels", "vip", "mapview",
  "data.table", "DALEXtra", "cowplot", "patchwork", "ggplot2", "parallel",
  "paletteer", "future", "doFuture", "ranger", "rstudioapi"
)

# Install any packages that aren't already installed
installed <- packages %in% rownames(installed.packages())
if (any(!installed)) {
  install.packages(packages[!installed])
}

# Load all packages
lapply(packages, library, character.only = TRUE)

# read in user input and set directories -------

# set working directory based on this script (needs to be within the broader folder)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

#read in user input
user_input_csv <- read.csv("user_input/user_input_general.csv", stringsAsFactors = FALSE)

# loop through user input csv and assign the values provided to the named objects
for (i in seq_len(nrow(user_input_csv))){
  val <- user_input_csv$value[i]
  if(user_input_csv$treat_as_R_expression[i] == T){val <- eval(parse(text = val))}
  assign(user_input_csv$name[i], val, envir = .GlobalEnv)
  rm(val)
}

# Load in wrapper functions 
source(paste0(script_directory, 'function_scripts/variable_reduction/var_reduction_wrappers.R'))
source(paste0(script_directory, 'function_scripts/variable_reduction/rf_rfe_helper_functions.R'))


# User Input -------------------------------------------------------------------

# Assign Main Directory & Data Filepath-----------------------------------------
# This does not need to be changed to follow the directory format of this project

# Directory to workout of 
local <- file.path(processed_data_directory, 'processed_outputs')

# Sub-directory and dataframe to conduct variable reduction 
# Must be within the local directory 
data_subdirectory <-  "clustering_outputs"
filename <- "all_samples_w_clusterid.csv"

# Variable crosswalk 
variable_crosswalk <-  paste0(global_data_directory, "/variable_descriptions/full_variable_list.csv")

# Assign Variable Reduction Inputs ---------------------------------------------
# Workflow step refers to the fire-level matching process ('Matching') or 
# the landscape-level final variable reduction ('Predictor')
# 'Matching' or 'Predictor' (column in 'variable_crosswalk' dataframe)
analysis_step <- "Matching" 

# Outcome set in the user input:
outcome <-  wto_outcome 

# Use a column in the dataframe to group the data: 
# In matching, we conduct variable reduction by cluster 
group_by <- 'cluster_id'

# Variables to retain in the variable reduction process:
# a) Retain the top variable from each selected variable group
retained_groups = c('Weather', 'Short-term Climate', 'Long-term Climate','Topography')
# b) specific variable to retain
retained_vars = c('class_cod')

# Should spatial layers be used:
# In matching, spatial layers are used (TRUE)
spatial_layers <- TRUE

# Should recursive feature elimination (rfe) be used for variable reduction:
# In matching, no rfe is used (FALSE)
rfe <-  FALSE

# What is the spatial distance of sampling 
# For burn severity, 270. 
sampling_dist <- 270

# Should the data be subsampled until there is no significant spatial autocorrelation 
# For burn severity set to FALSE
sp_subsampling  <-  FALSE

# Should the hyperparameters for the random forest models be tuned
# For matching (FALSE)
hyper_tune = FALSE


# Assign Other Directories & Filepaths------------------------------------------
# This does not need to change to follow the directory format of the project

# Figure output directory
FIGURE_DIR <- file.path(local, outcome, 'matching_model_outputs/figures')

# CSV output directory
CSV_DIR <- file.path(local, outcome, 'matching_model_outputs/csv')

# Model output directory
MODEL_DIR <- file.path(local, outcome, "matching_model_outputs/models")

# End of User Input ------------------------------------------------------------

# Set seed for reproducibility
set.seed(1234)

# Create output directories if they do not exist
if (!dir.exists(FIGURE_DIR)) {
  dir.create(FIGURE_DIR, recursive = TRUE)
}

if (!dir.exists(CSV_DIR)) {
  dir.create(CSV_DIR, recursive = TRUE)
}

if (!dir.exists(MODEL_DIR)) {
  dir.create(MODEL_DIR, recursive = TRUE)
}


# Set the response variable name
# This is set in earlier scripts, no need to change
if(response_type == 'binary' & outcome == 'burn_severity'){
  response_var <-  "response_bin"
  output_suffix <- paste0('_', response_var)
}else{
  response_var <- "response"
  output_suffix <- ''
}



# Load in data 
data_list <- load_fire_data()
df_raw <- data_list$df_raw
var_cw <- data_list$var_filtered


# Check if grouping column exists
if (!group_by %in% colnames(df_raw)) {
  stop(sprintf("Grouping column '%s' not found in data", group_by))
}else{
  # Get unique group values
  group_ids <- unique(df_raw[[group_by]])
  
  # Process each group
  final <- lapply(group_ids, function(id) {
    variable_reduction(
      group_id = id,
      group_col = group_by,
      df_raw = df_raw,
      var_filtered = data_list$var_filtered,
      resp_type = response_type,
      retained_vars = retained_vars,
      retained_groups = retained_groups,
      wto_outcome = outcome,
      response = response_var,
      spatial = spatial_layers,
      rfe = rfe, 
      hyper_tune = hyper_tune,
      sp_subsampling = sp_subsampling, 
      sampling_dist = sampling_dist
    )
  })
  
  names(final) <- group_ids
  
}
  
 

