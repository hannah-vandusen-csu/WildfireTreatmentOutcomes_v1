# WE CANNOT REDISTRIBUTE THIS SCRIPT
# PLEASE GAIN ACCESS TO LATEST CODE DIRECTLY FROM AUTHOR (IN CASE ANY MODIFICATIONS WERE MADE)

# Author: Sean Parks (sean.parks@usda.gov)
# Rationale for this procedure, validation of results, and comparison to 
# other interpolation methods can be found in (please cite this paper if you use this code): 
# Parks. 2014. Mapping day-of-burning with coarse-resolution satellite fire-detection data. International Journal of Wildland Fire. 23:215-223.

# Parameters used for the day of burning calculation (DOB) are coded below
# This workflow uses hotspot data from VIIRS and MODIS,
# the user can then add the DOB calculations once they have gain permission from the orginal author





rm(list=ls(all=T))

pkgs <- c(
  "terra",
  "plyr",
  "dplyr",
  "tidyr",
  "stringr",
  "data.table",
  "rstudioapi", 
  "lubridate",
  "sf",
  "igraph",
  "timeDate", 
  "FNN"
)

# Install missing CRAN packages
missing_pkgs <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(missing_pkgs) > 0) {
  install.packages(missing_pkgs)
}
# Load all packages
lapply(pkgs, library, character.only = TRUE)

# directories for data and scripts are defined in "user_input_general.csv"
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# read in user input
user_input_csv <- read.csv("user_input/user_input_general.csv", stringsAsFactors = FALSE)

# loop through user input csv and assign the values provided to the named objects
for (i in seq_len(nrow(user_input_csv))){
  val <- user_input_csv$value[i]
  if(user_input_csv$treat_as_R_expression[i] == T){val <- eval(parse(text = val))}
  assign(user_input_csv$name[i], val, envir = .GlobalEnv)
  rm(val)
}

# Navigate to whatever folder you're working
#5/29- weird bug right now where choose.dir() wont run until choose.files() is called
processed_data_directory <- file.path(choose.dir())
processed_data_directory <- gsub('\\\\', '/', processed_data_directory) #for windows
setwd(processed_data_directory)

########### You need 4 items to run this script, in their respective folders  #################
#these are all folders that need to be located within the processed_data_directory folder, change the folder name as needed
sub_folders <- c(
  firms = 'fire_detections/firms_fire_detections', #FIRMS data from https://firms.modaps.eosdis.nasa.gov/download
  perim = 'merged_perims', #fire perimeters of an area of interest, or all from https://mtbs.gov/direct-download
  fire_list = 'selectedfires', #List of fires of interest with an mtbs ID column
  fire_id_column = 'MTBS_ID', #name of the column that corresponds to the name/ID of the fire in the .csv of fires
  )

###
# Import CSV of names and Ignition Dates
###

fires <- read.csv(file.path(processed_data_directory, 'fire_lists', 'selectedfires.csv'))
fires$year <- lubridate::year(fires$Ignition)
fires$yday <- lubridate::yday(fires$Ignition)
#fires <- fires |> filter(yday > 100) # Filter out fires that occurred before April 1st

years <- unique(fires$year) %>%
  sort()# Run all years 2001 - 20221


##### Set your processing parameters here #########
pixel_size = 30       #resolution, in meters
day_frac = 1          #fraction of day to divide FIRMS data, 0.5 = half day
round_type = floor    #options: round, floor, ceiling see ?round_any() for more info
shift_hour = 4        #number of hours to shift the hotspot timestamp, moves 3AM VIIRS pass to day prior
time_zone = 8 #7 for MST, need to figure out way to automate this (DST)
regions = 5 #the size of regions in hectares (and smaller to remove)
buffer = 1000 #size to buffer the fire perim by, 0 for no buffer
the_prj = CRS
selected_years = 'all years' #need to work on this
nifc_perims = F
date_format = 'yymmdd' # or 'julian'
output_directory = processed_data_directory


if (day_frac < 1) { #changes the raster type to float values (0.00..) if opting for sub daily
  data_type = 'FLT4S'
} else {
  data_type = 'FLT8S'
}


###
# Load hotspot data for all years and all of US
###

## Get MODIS fire detections
{
hotspotsModis = c(); hotspotsViirsSpp = c(); hotspotsViirsNoaa20 = c(); hotspotsViirsNoaa21 = c()
hotspotsModis <- st_read(list.files(file.path(global_data_directory,'fire_detections/firms_fire_detections', 'modis'), pattern = '.shp', full.names = T))
hotspotsModis <- subset(hotspotsModis, select = -c(TYPE))

## Get VIIRS fire detections
hotspotsViirsSpp <- st_read(paste0(global_data_directory, '/', sub_folders[['firms']], '/viirs'), 'fire_archive_SV-C2_431932')
hotspotsViirsSpp <- subset(hotspotsViirsSpp, select = -c(TYPE))
hotspotsViirsNoaa20 <- st_read(paste0(global_data_directory, '/', sub_folders[['firms']], '/viirs'), 'fire_nrt_J1V-C2_458502')
hotspotsViirsNoaa21 <- st_read(paste0(global_data_directory, '/', sub_folders[['firms']], '/viirs'), 'fire_nrt_J2V-C2_458503')

# Combine/clean/transfrom into one folder, 'hotspots'
hotspots <- rbind(hotspotsModis, hotspotsViirsSpp, hotspotsViirsNoaa20, hotspotsViirsNoaa21)
rm(hotspotsModis, hotspotsViirsSpp, hotspotsViirsNoaa20, hotspotsViirsNoaa21)
hotspots <- hotspots[, c('LATITUDE', 'LONGITUDE', 'ACQ_DATE', 'ACQ_TIME', 'SATELLITE', 'INSTRUMENT', 'FRP')]
hotspots <- st_transform(hotspots, crs = crs(the_prj))
hotspots$years <- lubridate::year(hotspots$ACQ_DATE) # Set column for year
}

##
## Start of function.... 
##

##
## Please copy and paste DOB functions below 
##

##
## End of DOB function here.... 
##
################### Aggregate all DOB rasters by year #######################

# List all dob fires 
setwd(paste0(processed_data_directory, '/data/interim_data/date_of_burn/dob_tif'))
if (nifc_perims == T) {
  dob_files <- list.files(pattern = paste0('}_', day_frac, '_dob.tif$'))
} else {
  dob_files <-  list.files(pattern = paste0('^[^}]*_', day_frac, '_dob.tif$'))
}



# instead of julian day make date in format yyMMDD
#LEO: my suggestion would be to re-write all code to retain the yy/mm attribute in DOB per fire as julian date
add_yymmdd <- function(x){
  
  # download fire ID
  fire <- gsub(paste0('_', day_frac, '_dob.tif'),'', x)
  
  # Fire burn year
  year <- fires[fires[[sub_folders[['fire_id_column']]]] %in% fire,]$year
  
  dob <- terra::rast(paste0(output_directory,'/date_of_burn/dob_tif/', x))
  
  # Add year
  origin_year <- as.Date(paste0(year, '-01-01'))
  
  
  # Function to convert Julian Days to YYMMDD
  julian_day_to_yymmdd <- function(juliandays){
    #hack for leap year
    if(year %in% seq(2000,2032, 4)){
      print("accounting for leap year")
      date_obj <- origin_year + as.integer(juliandays)
    }else(date_obj <- origin_year + as.integer(juliandays) - 1)
    
    ##date_obj <- as.Date(juliandays, origin = origin_year)
    formatted_date <- format(date_obj, '%Y%m%d')
    return(formatted_date)
  }
  
  # as DF
  dob_df <- as.data.frame(dob, xy = T)
  
  # Run function to add year
  DOB <- julian_day_to_yymmdd(dob_df$DOB) |> as.numeric()
  
  # Range of days 
  # days <- na.omit(DOB)%% 100
  # range_days <- range(days)
  # print(range_days)

  #add coords and DOB
  dob_df$DOB <- DOB

  # Create a Raster
  modeled.dob.year <- rast(dob_df, type = 'xyz', crs=(the_prj))
  #writeRaster(modeled.dob.year, "Test.tif", overwrite = T, datatype = 'INT4U')
  
  return(modeled.dob.year)
}

# Function that creates a raster stack per year with formated date
stack_per_year <- function(y){
  print(paste0("year: ", y))
  IDs <- fires[fires$year == y,][[sub_folders[['fire_id_column']]]]
  
  dob_files_year <- dob_files[gsub(paste0('_', day_frac, '_dob.tif'),'', dob_files) %in% IDs]
  
  if (length(dob_files_year) > 0) {
    # Apply Function to all rasters
    if (date_format == "julian") {
      print("    converting to yyyymmdd")
      dob_stack <- lapply(dob_files_year, add_yymmdd)
    } else {
      dob_stack <- lapply(dob_files_year, rast)
    }
  
    #large_rast <- do.call(terra::mosaic, args = c(dob_stack, fun="modal"))
    # Combine into one large raster for the year
    if(length(dob_stack) > 1){
      # large_rast <- do.call(merge, dob_stack)
      large_sprc <- sprc(dob_stack)
      # CAUTION having trouble with merge rounding dates to even values, may be a RAM issue (restart R with clean env)
      large_rast <- terra::merge(large_sprc)
      #large_rast <- terra::mosaic(large_sprc, fun="min")
    } else {
      large_rast <- dob_stack[[1]]
    }
    
  
    #large_rast <- app(large_rast, function(x){x*10^6})
    if(!dir.exists(paste0(output_directory,'/date_of_burn/dob_by_years'))){
      dir.create(paste0(output_directory,'/date_of_burn/dob_by_years'))
    }
  
    # Write raster for year
    # Must save a 4-byte integer to preserve data values
  
    writeRaster(large_rast, paste0(output_directory,'/date_of_burn/dob_by_years/', sub_folders[[2]], '_', day_frac, '_dob_', y,'.tif'),  datatype = data_type, overwrite = T)

  }
  else{
    print(paste0('no fires for ', y))
  }
}

lapply(years, stack_per_year)



#
