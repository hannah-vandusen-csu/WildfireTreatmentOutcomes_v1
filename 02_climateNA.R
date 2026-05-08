######################################################.
#
# Climate NA Climate Normal Extraction
#
# Author: Betsy Black elizabeth.black@usda.gov June 2024
# Edited By: Sarah Hettema July 2025
# re-structured and edited by: Mark Kreider, February 2026
#
######################################################.



# load packages -----

# Vector of CRAN packages
pkgs <- c(
  "terra",
  "dplyr",
  "tidyr",
  "stringr",
  "data.table",
  "rstudioapi"
)

# Install missing CRAN packages
missing_pkgs <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(missing_pkgs) > 0) {
  install.packages(missing_pkgs)
}

# Install climatenaR from GitHub (if needed)
if (!("climatenaR" %in% installed.packages()[, "Package"])) {
  if (!("devtools" %in% installed.packages()[, "Package"])) {
    install.packages("devtools")
  }
  library(devtools)
  devtools::install_github("burnett-m/climatenaR", build_vignettes = TRUE)
}

# Load all packages
lapply(c("climatenaR", pkgs), library, character.only = TRUE)

# read in user inputs and set directories -----

# directories for data and scripts are defined in "user_input_general.csv"
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

# set directory where climateNA application and input data lives
climateNA_directory <- file.path(global_data_directory, "climate/ClimateNA_v750/")

# Create necessary directories
if(dir.exists(paste0(processed_data_directory,'/climate')) == FALSE){
  dir.create(paste0(processed_data_directory,'/climate'))}
if(dir.exists(paste0(processed_data_directory,'/climate/mid_processing')) == FALSE){
  dir.create(paste0(processed_data_directory,'/climate/mid_processing'))}
if(dir.exists(paste0(processed_data_directory,'/climate/climate_normals')) == FALSE){
  dir.create(paste0(processed_data_directory,'/climate/climate_normals'))}

# Crop DEM to strategy landscape and buffered fires ----

#Bring in global DEM
dem_path <- paste0(global_data_directory,'/climate/dem90_hf/dem90_hf.tif') 
dem_global <- rast(dem_path)

#Select Strategy Landscapes
Strategy_landscapes <- vect(paste0(global_data_directory,'/strategy_landscapes/AllPriorityLandscapes.gdb'))
unique(Strategy_landscapes$NAME)

# Select names to include the shape of the landscape of interest
partial_AOI <- Strategy_landscapes[Strategy_landscapes$NAME %in% strategy_landscapes_of_interest,]

# Define path to fire perimeter shapefiles
fire_extents_path <- file.path(processed_data_directory, "/fire_perimeters/merged/")  

# Extract MTBS_IDs from file names
list_f <- list.files(fire_extents_path, pattern = "\\.shp$", full.names = FALSE) %>%
  str_remove("_Extent.shp") %>% 
  unique()

# Read and buffer each fire perimeter, combine with strategy landscapes shapes
fires <- vect()
for (fire_id in list_f) {
  message("Adding ", fire_id)
  shp_path <- file.path(fire_extents_path, paste0(fire_id, "_Extent.shp"))
  extent <- vect(shp_path)
  extent <- terra::buffer(extent, 1000)
  fires <- rbind(fires, extent)
}
partial_AOI<-project(partial_AOI,crs(fires))
AOI <- partial_AOI %>% 
  rbind(fires) %>% 
  aggregate() %>%
  fillHoles() %>%
  project(crs(dem_global))


#Crop DEM to AOI
dem<-crop(dem_global,AOI,mask=TRUE)

# #Save DEM
dem_path <- paste0(processed_data_directory,'/climate/mid_processing/dem_AOI.tif')
writeRaster(dem, dem_path, overwrite=T)

# Extract coordinates and elevation to csv files --------

#Output filepath
out_fp <- paste0(processed_data_directory,'/climate/mid_processing')

#Extract lat/long and elevation CSV
dem<-project(dem,"EPSG:4326") #DEM must be in "EPSG:4326"

# convert raster to coordinates (with elevation data associated)
dem_df <- raster::as.data.frame(dem, xy = T, na.rm=  T) %>%
  mutate(x = x * -1) %>% # Positive numbers required for longtitude
  mutate(ID1 = rownames(.),
         ID2 = rownames(.) # A second ID column must be made to produce the climate variables
         )

# get dataframe ready to export as csv
dem_csv <- dem_df[c(4,5,2,1,3)]
names(dem_csv) <- c('ID1','ID2','lat','long','el')
data.table::fwrite(dem_csv, file = paste0(out_fp,'/dem_AOI.csv'),
                   row.names=FALSE)


# Use lat/long/elevation to extract climate data -------

#lat long elevation (lle) csv path
out_fp <- paste0(processed_data_directory,'/climate/mid_processing')
lle_path <- paste0(out_fp,'/dem_AOI.csv')
lle <- read.csv(lle_path) #Inspect output

#Path to climatena package
exe <- file.path(climateNA_directory,"ClimateNA_v7.50.exe")

#SELECT TEMPORAL PARAMETERS

# DATE RANGE: 30, 10, or single year period across which to get mean climate attributes
### if year range, must be formatted for whole decade chunks: YYY1_YYY0
### if single year: YYYY
date_range<-'1991_2020'

# TIME FRAME: Climate attribute sampling windows: annual (Y), seasonal (S), or monthly (M)
time_frame<-'Y'

#Specify output directory
out_path <- paste0(out_fp,'/extracted_',time_frame)
if(dir.exists(paste0(out_path)) == FALSE){
  dir.create(paste0(out_path))}

#Set WD to where climate NA app is
setwd(climateNA_directory)

#EXTRACT CLIMATE DATA - THIS STEP TAKES A LONG TIME
histClimateNA(lle_path,date_range,time_frame,exe,out_path)



# Export rasters of climate normals for individual fires ------ 


#Must export rasters to WD so make a folder there
if(dir.exists(paste0(climateNA_directory,'/raster_exports')) == FALSE){
  dir.create(paste0(climateNA_directory,'/raster_exports'))}

#SELECT TEMPORAL PARAMETERS TO EXPORT

# DATE RANGE: 30, 10, or single year period across which to get mean climate attributes
### if year range, must be formatted for whole decade chunks: YYY1_YYY0
### if single year: YYYY
date_range<-'1991_2020'

# TIME FRAME: Climate attribute sampling windows: annual (Y), seasonal (S), or monthly (M)
time_frame<-'Y'

#Read in extracted point data
out_fp <- paste0(processed_data_directory,'/climate/mid_processing')
climate_path <- paste0(out_fp,'/extracted_',time_frame,'/dem_AOI_',date_range,'.csv')
climate<-read.csv(climate_path) #Inspect

#SELECT METRICS TO EXPORT

#List available metrics (for input climate CSV)
names(climate)

#Variable names and definitions for each range (annual, seasonal, or monthly) can be found on https://climatena.ca/Help2

#Metric overview:

##Annual-scale climate variables##

# MAT	mean annual temperature (°C)
# MWMT	mean warmest month temperature (°C)
# MCMT	mean coldest month temperature (°C)
# TD	temperature difference between MWMT and MCMT, or continentality (°C)
# MAP	mean annual precipitation (mm)
# MSP May to September precipitation (mm)
# AHM	annual heat-moisture index (MAT+10)/(MAP/1000))
# SHM	summer heat-moisture index ((MWMT)/(MSP/1000))  	 
# DD_0  degree-days below 0°C, chilling degree-days
# DD5  degree-days above 5°C, growing degree-days
# DD_18  degree-days below 18°C, heating degree-days
# DD18	degree-days above 18°C, cooling degree-days
# NFFD	the number of frost-free days
# FFP	frost-free period
# bFFP	the day of the year on which FFP begins
# eFFP	the day of the year on which FFP ends
# PAS	precipitation as snow (mm). For individual years, it covers the period between August in the previous year and July in the current year
# EMT	extreme minimum temperature over 30 years (°C)
# EXT	extreme maximum temperature over 30 years (°C)
# Eref	Hargreaves reference evaporation (mm)
# CMD	Hargreaves climatic moisture deficit (mm)
# MAR	mean annual solar radiation (MJ m‐2 d‐1)
# RH	mean annual relative humidity (%)
# CMI	Hogg’s climate moisture index (mm)
# DD1040	degree-days above 10°C and below 40°C


##Seasonal-scale climate variables##

# Tave_wt 	winter mean temperature (°C)
# Tave_sp	spring mean temperature (°C)
# Tave_sm	summer mean temperature (°C)
# Tave_at	autumn mean temperature (°C)
# Tmax_wt	winter mean maximum temperature (°C)
# Tmax_sp	spring mean maximum temperature (°C)
# Tmax_sm	summer mean maximum temperature (°C)
# Tmax_at  	autumn mean maximum temperature (°C)
# Tmin_wt	winter mean minimum temperature (°C)
# Tmin_sp	spring mean minimum temperature (°C)
# Tmin_sm	summer mean minimum temperature (°C)
# Tmin_at	autumn mean minimum temperature (°C)
# PPT_wt	winter precipitation (mm)
# PPT_sp	spring precipitation (mm)
# PPT_sm	summer precipitation (mm)
# PPT_at	autumn precipitation (mm)
# RAD_wt	winter solar radiation (MJ m-2 d-1)
# RAD_sp	spring solar radiation (MJ m-2 d-1)
# RAD_sm	summer solar radiation (MJ m-2 d-1)
# RAD_at	autumn solar radiation (MJ m-2 d-1)
# DD_0_wt ~ DD_0_at	winter ~ Autumn degree-days below 0°C
# DD5_wt ~ DD5_at	winter ~ Autumn degree-days above 5°C
# DD_18_wt ~ DD_18_at	winter ~ Autumn degree-days below 18°C
# DD18_wt ~ DD18_at	winter ~ Autumn degree-days above 18°C
# NFFD_wt ~ NFFD_at	winter ~ Autumn number of frost-free days
# PAS_wt ~ PAS_at	winter ~ Autumn precipitation as snow (mm)
# Eref_wt ~ Eref_at	winter ~ Autumn Hargreaves reference evaporation (mm)
# CMD_wt ~ CMD_at	winter ~ Autumn Hargreaves climatic moisture deficit (mm)
# RH_wt ~ RH_at	winter ~ Autumn relative humidity (%)
# CMI_wt ~ CMI_at	winter ~ Autumn Hogg’s climate moisture index (mm)


##Monthly-scale climate variables##

# Tave01 – Tave12 January - December mean temperatures (°C)
# TMX01 – TMX12 January - December maximum mean temperatures (°C)
# TMN01 – TMN12 January - December minimum mean temperatures (°C)
# PPT01 – PPT12 January - December precipitation (mm)
# RAD01 – RAD12 January - December solar radiation (MJ m‐2 d‐1)
# DD_0_01 – DD_0_12 January - December degree-days below 0°C
# DD5_01 – DD5_12 January - December degree-days above 5°C
# DD_18_01 – DD_18_12 January - December degree-days below 18°C
# DD18_01 – DD18_12 January - December degree-days above 18°C
# NFFD01 – NFFD12 January - December number of frost-free days
# PAS01 – PAS12 January – December precipitation as snow (mm)
# Eref01 – Eref12 January – December Hargreaves reference evaporation (mm)
# CMD01 – CMD12 January – December Hargreaves climatic moisture deficit (mm)
# RH01 – RH12 January – December relative humidity (%)
# CMI01 – CMI12 January – December Hogg’s climate moisture index (mm)

#Create list of selected variables to export rasters of
var_list <- list("MAT","MWMT","MCMT","MAP","MSP","AHM","SHM","PAS","Eref","CMD","RH","CMI")

#Export whole-DEM rasters of selected climate metrics
CreateNormalsRasters <- function(var){
  print(paste0("Working on ",var))
  
  #Specify output path
  vartif_fp <-paste0('/raster_exports/',var,'_',date_range,'_',time_frame,'.tif')
  
  #Transform xy point csvs into rasters for specified climate metric
  CSVtoTIFF(
    Longitude=climate$Longitude,
    Latitude=climate$Latitude,
    Value=climate[[var]],
    crsref = "+proj=longlat +datum=WGS84",
    vartif_fp)
  
}

for(var in var_list){
  CreateNormalsRasters(var)
}


# Import csv of fire names
fires <- read.csv(paste0(processed_data_directory, '/fire_lists/selectedfires.csv'))

# Create list of fire names from csv
list_f <- fires$MTBS_ID

#Bring in fire perimeters
fire_perims <- lapply(list_f, function(fire_id) {
  print(fire_id)
  extent_fp <- file.path(processed_data_directory, "/fire_perimeters/merged", paste0(fire_id, "_Extent.shp")) 
  extent <- vect(extent_fp)
  extent <- project(extent, crs(CRS))
  buffer_extent <- terra::buffer(extent, width = 1000)
  return(buffer_extent)
})
names(fire_perims) <- list_f  # Name the list by fire ID


Normals <- function(var){
  print(paste0("Working on ",var))
  
  #Path to current var raster
  vartif_fp <-paste0(climateNA_directory,'/raster_exports/',var,'_',date_range,'_',time_frame,'.tif')
  
  #Pull raster in
  climate_raster <- rast(paste0(vartif_fp))
  
  #Turn NA values (-9999) into actual NAs
  climate_raster[climate_raster == -9999] <- NA
  
  #Reproject raster to fire perimeters' projections
  climate_raster=project(climate_raster,crs(CRS))
  
  #Loop through fire list and clip out each fire perim
  for (fire_id in names(fire_perims)) {
    message("  Clipping to fire: ", fire_id)
    
    cropped <- crop(climate_raster, fire_perims[[fire_id]], mask = TRUE)
    
    out_fp <- file.path(processed_data_directory, "climate/climate_normals", paste0(fire_id, "_", var, ".tif"))
    writeRaster(cropped, out_fp, overwrite = TRUE)
  }
  
}

for(var in var_list){
  Normals(var)
}

##Seasonal climate variables##
# Tave_wt 	winter mean temperature (°C)
# Tave_sp	spring mean temperature (°C)
# Tave_sm	summer mean temperature (°C)
# Tave_at	autumn mean temperature (°C)
# Tmax_wt	winter mean maximum temperature (°C)
# Tmax_sp	spring mean maximum temperature (°C)
# Tmax_sm	summer mean maximum temperature (°C)
# Tmax_at  	autumn mean maximum temperature (°C)
# Tmin_wt	winter mean minimum temperature (°C)
# Tmin_sp	spring mean minimum temperature (°C)
# Tmin_sm	summer mean minimum temperature (°C)
# Tmin_at	autumn mean minimum temperature (°C)
# PPT_wt	winter precipitation (mm)
# PPT_sp	spring precipitation (mm)
# PPT_sm	summer precipitation (mm)
# PPT_at	autumn precipitation (mm)
# RAD_wt	winter solar radiation (MJ m-2 d-1)
# RAD_sp	spring solar radiation (MJ m-2 d-1)
# RAD_sm	summer solar radiation (MJ m-2 d-1)
# RAD_at	autumn solar radiation (MJ m-2 d-1)
# DD_0_wt ~ DD_0_at	winter ~ Autumn degree-days below 0°C
# DD5_wt ~ DD5_at	winter ~ Autumn degree-days above 5°C
# DD_18_wt ~ DD_18_at	winter ~ Autumn degree-days below 18°C
# DD18_wt ~ DD18_at	winter ~ Autumn degree-days above 18°C
# NFFD_wt ~ NFFD_at	winter ~ Autumn number of frost-free days
# PAS_wt ~ PAS_at	winter ~ Autumn precipitation as snow (mm)
# Eref_wt ~ Eref_at	winter ~ Autumn Hargreaves reference evaporation (mm)
# CMD_wt ~ CMD_at	winter ~ Autumn Hargreaves climatic moisture deficit (mm)
# RH_wt ~ RH_at	winter ~ Autumn relative humidity (%)
# CMI_wt ~ CMI_at	winter ~ Autumn Hogg’s climate moisture index (mm)

var_list <- list("Tave_sp","Tave_sm","Tave_at","Tave_wt","PPT_sp","PPT_sm","PPT_at","PPT_wt")

# DATE RANGE: 30, 10, or single year period across which to get mean climate attributes
### if year range, must be formatted for whole decade chunks: YYY1_YYY0
### if single year: YYYY
date_range<-'1991_2020'

# TIME FRAME: Climate attribute sampling windows: annual (Y), seasonal (S), or monthly (M)
time_frame<-'S'

#Read in extracted point data
out_fp <- paste0(processed_data_directory,'/climate/mid_processing')
climate_path <- paste0(out_fp,'/extracted_',time_frame,'/dem_AOI_',date_range,'.csv')
climate<-read.csv(climate_path) #Inspect

#SELECT METRICS TO EXPORT

#List available metrics (for input climate CSV)
names(climate)

for(var in var_list){
  CreateNormalsRasters(var)
}

for(var in var_list){
  Normals(var)
}
# calculate evapotranspiration --------

#ET is not a supplied metric, but ET = Eref - CMD
#Optional extract ET

#Loop through fire list and clip out each fire perim
for(fire in list_f){
  print(fire)
  #Load Eref
  Eref<- rast(paste0(processed_data_directory, "/climate/climate_normals/", fire, "_Eref.tif"))
  
  #Load CMD
  CMD<- rast(paste0(processed_data_directory, "/climate/climate_normals/", fire, "_CMD.tif"))
  
  #Check that geom matches
  if (!compareGeom(Eref,CMD)){
    stop("Eref and CMD rasters do not have the same extent and resolution.")
  }
  
  #Calculate ET
  ET <- Eref - CMD
  
  #Load fire perimeter and make sure crs matches
  extent<- vect(paste0(processed_data_directory, "/fire_perimeters/merged/", fire, "_Extent.shp"))
  extent=project(extent,CRS)
  
  #Also grab 1km buffer
  buffer_distance <- 1000
  buffered_extent <- terra::buffer(extent, width = buffer_distance)
  
  #Clip out fire perim
  cropped <- crop(ET,buffered_extent,mask=TRUE)
  
  #Write out rasters
  this_path <- paste0(processed_data_directory,"/climate/climate_normals/",fire,"_ET.tif")
  writeRaster(cropped,this_path,overwrite=TRUE)
  
}




