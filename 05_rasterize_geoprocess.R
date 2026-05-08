####################################.
###
### Rasterize data for upwind analysis
###
### Written by Betsy Black April 2024 (elizabeth.black@usda.gov)
### Re-organized into single script by Mark Kreider February 2026
###
####################################.

rm(list=ls())

# load packages -----

# Necessary packages
pkgs <- c("terra",
          "tidyterra",
          "tidyverse",
          "sf",
          "rstudioapi", 
          "data.table",
          "tigris")


# Install any that are missing
missing <- pkgs[!(pkgs %in% installed.packages()[,"Package"])]
if (length(missing) > 0) {install.packages(missing)}

# Load all packages
lapply(pkgs, library, character.only = TRUE)

sf_use_s2(FALSE)

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

# load functions -----
# load required functions from other scripts
function_scripts <- list.files(file.path(script_directory, "function_scripts/rasterize_geoprocess"), full.names = TRUE) 
sapply(function_scripts, source)

# Rasterize data -------

# each of these functions runs a component of the rasterizing data workflow

# Download small fires and clip to fire perimeters 
smfire_download(data_directory = global_data_directory,
                output_directory = processed_data_directory,
                the.prj = CRS)

# extract rasters for small fires 
smfire_download_rast(data_directory = global_data_directory,
                     output_directory = processed_data_directory,
                     the.prj = CRS)

# create blank rasters for fires without previous small fires
smfire_download_dummy_rast(data_directory = global_data_directory,
                           output_directory = processed_data_directory,
                           the.prj = CRS)

# create blank rasters for treatments for fires without them
create_dummy_treatment_rasters(data_directory = global_data_directory,
                               output_directory = processed_data_directory,
                               the.prj = CRS)

# extract land ownership polygons for each fire
create_land_ownership_shapefiles(data_directory = global_data_directory,
                                 output_directory = processed_data_directory,
                                 the.prj = CRS)

# create land ownership rasters
create_land_ownership_rasters(output_directory = processed_data_directory)


# get roads data
get_roads_data(global_data_directory = global_data_directory,
               strategy_landscapes_of_interest = strategy_landscapes_of_interest,
               output_directory = processed_data_directory,
               the.prj = CRS)

# rasterize roads data
rasterize_roads_by_fire(output_directory = global_data_directory,
                        the.prj = CRS)

# make all areas outside of perimeter 24 hours later burn time 
# (this is necessary to ensure upwind fans face into the fire)
add_plus_one_day(processed_data_directory = processed_data_directory)
