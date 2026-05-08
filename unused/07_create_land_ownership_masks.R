####################################
###
### Land Ownership Data Download
### Exports a shapefile for each Fire Perimeter
###
### Written by Betsy Black May 2024 (elizabeth.black@usda.gov)
###
####################################

rm(list=ls())

# load packages -----

# Necessary packages
pkgs <- c("terra",
          "tidyterra",
          "tidyverse",
          "rstudioapi",
          "sf",
          "stringr",
          "raster"
)


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



#Load land ownership gdb####
fp <- paste0(global_data_directory, "/land_ownership/SMA_WM.gdb")  
own <- vect(fp, layer="SurfaceManagementAgency")
own <- project(own, CRS)

# Import CSV of names and Ignition Dates
fires <- read.csv(paste0(processed_data_directory,"/fire_lists/selectedfires.csv"))

# Create list of fire names from CSV
list_f <- fires$MTBS_ID
list_f = gsub(' ','_', list_f)

# create necessary directories for desired files paths
if(dir.exists(paste0(processed_data_directory,'/masks')) == FALSE){
  dir.create(paste0(processed_data_directory,'/masks'))
} 
# create necessary directories for desired files paths
if(dir.exists(paste0(processed_data_directory,'/masks/land_ownership')) == FALSE){
  dir.create(paste0(processed_data_directory,'/masks/land_ownership'))
} 

# Function to clip FACTS + NFPORS treatment layers to fire perimeters, 
# output shapefile and csv of clipped treatment data
ClipLO2Fire <- function(this_fire){
  
  #START A FIRE####
  print(paste0("Processing ",this_fire))
  
  # Load fire perimeter
  fire_perim <- vect(paste0(processed_data_directory, '/fire_perimeters/merged/',this_fire,'_Extent.shp'))
  fire_perim <- project(fire_perim, crs(own))
  
  #Add buffer to fire perim
  buffered_extent<-buffer(fire_perim,1000) #in meters bc epsg 26910 projection
  
  #Clip ownership to a 1000m buffer around fire extent
  own.clipped <- terra::crop(own, buffered_extent)

  ####
  
  # create file path inside folder for shapefiles
  own.fp <- (paste0(processed_data_directory,'/masks/land_ownership/'))
  own.nm <- paste0(this_fire, "_land_own.shp")

  # Write Shape file
  writeVector(own.clipped, paste0(own.fp, own.nm), overwrite = T)
  print(paste0('land ownership for ',this_fire,' clipped'))
  
}

# Iterate over each fire and wui piece
for (i in seq_along(list_f)) {
  ClipLO2Fire(list_f[i])
}
