####################################################.
##                                                ##.
## Extract data from potential sample points in R ##. 
##                                                ##.
####################################################.

# load packages -----

# Necessary packages
pkgs <- c("terra",
          "tidyverse",
          "rstudioapi",
          "sf",
          "stringr",
          "lubridate",
          "DescTools",
          "fuzzyjoin",
          "data.table",
          "tidyterra"
)


# Install any that are missing
missing <- pkgs[!(pkgs %in% installed.packages()[,"Package"])]
if (length(missing) > 0) {install.packages(missing)}

# Load all packages
lapply(pkgs, library, character.only = TRUE)

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

# set local directory 
local_directory <- processed_data_directory

# set inputs
# Outcome : ('burn_severity', 'containment', 'infrastructure' )
wto_outcome = 'burn_severity'
# Outcome type : ('burn_severity', 'containment', 'infrastructure', 'infrastructure_buff' )
column = 'burn_severity'

## Sampling points csvs
# Location of sampling points & file suffix (should be the same as gee)
# each point should have a column "point_id" that is a unique identifier
point_directory = paste0(local_directory, '/initial_sampling_points/', wto_outcome) 
filename_suffix =  '_initial_pts' 
unique_identifer = 'point_id'

## GEE data path where data is extracted by a per fire basis
# The default format fo the file path should be:
gee_data <- paste0('/processed_outputs/' , wto_outcome, '/gee_point_extracts')

# 5) suffix of gee files: the default is "_df". 
# suffix will also be assigned to final dataframe: default is "final_df"
# If creating multiple final dataframes there might by other suffixes
gee_suffix <- 'initial_pts_df'

# 7) Final Output directory for data frame
# Default '/final_dataframes/ + WTO outcome + /simple_extract/
# additionally sub-directories could be necessary
ouput_fp <- paste0('/processed_outputs/' , wto_outcome, '/extracted_dataframes/point_level_extracts')


# 8) List of fires saved in the fire_list folder
fire_list_fp <- paste0(wto_outcome, '_fire_list')

# 9) Variable crosswalks
variable_crosswalk <- '/variable_crosswalks_simple.csv'

# 10) Buffer
buffer <- FALSE # TRUE or FALSE
buffer_size <- 0

# load functions -----
source(file.path(script_directory, 'function_scripts/R_extraction/simple_extraction_function.R')) 


# run function to extract data for points
simple_extraction_local(local_directory = local_directory,
                         wto_outcome = wto_outcome,
                         column = column,
                         buffer,
                         buffer_size,
                         point_directory = point_directory,
                         filename_suffix = filename_suffix,
                         gee_data = gee_data,
                         gee_suffix = gee_suffix,
                         ouput_fp = ouput_fp,
                         fire_list_fp = fire_list_fp, 
                         var_df = variable_crosswalk)


