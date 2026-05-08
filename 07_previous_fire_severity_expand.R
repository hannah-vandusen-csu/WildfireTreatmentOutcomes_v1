#######################################################.
###                                                 ###.
### Script to Expand previous fire severity rasters ###.
###                                                 ###.
#######################################################.

rm(list=ls())

# load packages -----

# Necessary packages
pkgs <- c("terra",
          "tidyterra",
          "tidyverse",
          "rstudioapi",
          "sf"
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

# Set directories
data_directory <- global_data_directory
output_directory <- processed_data_directory
the.prj <- CRS


# expand small fire severity -------

#Call in list of fires
fires <- read.csv(paste0(output_directory,"/fire_lists/selectedfires.csv"))

# Create list of fire names from CSV
list_f <- fires$MTBS_ID

# create directory to store expanded fire severity rasters
if(dir.exists(paste0(processed_data_directory,'/smfires/raster/past_sev_corrected/')) == FALSE){
  dir.create(paste0(processed_data_directory,'/smfires/raster/past_sev_corrected/'))
  }

ByFire <- function(this_fire){
  
  print(paste0("Processing ",this_fire))
  
  #Bring in fire perim
  perim <- vect(paste0(output_directory,"/fire_perimeters/merged/",this_fire,"_Extent.shp")) %>%
    project(crs(the.prj))

  #buffer
  per_buf <- terra::buffer(perim, 1000)
  
  # get classified past severity raster path
  raster_path <- paste0(output_directory, "/smfires/raster/", this_fire, "_classified_past_sev.tif")
  
  #Create Buffered, fixed raster
  if (file.exists(raster_path)){
    print(paste0("File exists for ",this_fire))
    
    #If yes previous fires and we exported prev fire severity from GEE, fix and expand raster:
    pf <- rast(raster_path) %>%
      project(crs(per_buf))

    # Create a new raster with per_buf extent, same resolution, and fill with 0
    template <- rast(ext=ext(per_buf), 
                     resolution=res(pf), 
                     crs=crs(per_buf))
    template[] <- 0
    
    # Align pf to template (expand pf to template extent, filling with NA)
    pf <- extend(pf, template)
    
    # Fill NA (from new cells) with 0
    pf[is.na(pf)] <- 0
    
  }else{
    print(paste0("File does not exist for ",this_fire))
    #If no previous fires and we did not export prev fire severity from GEE make dummy raster:
    dummy_res <- 26.36652
    
    # Create dummy raster matching per_buf's extent and CRS, one band, name "CBI_bc"
    pf <- rast(ext=ext(per_buf), 
                       resolution=dummy_res, 
                       crs=crs(per_buf), 
                       nlyrs=1)
    names(pf) <- "CBI_bc"
    pf[] <- 0
  
  }
  
  writeRaster(pf,paste0(output_directory,"/smfires/raster/past_sev_corrected/",this_fire,"_classified_past_sev.tif"),overwrite=T)
  print(paste0("Saved raster for ",this_fire))
}

DoThatThang<-lapply(1:length(list_f), function(i) {
  ByFire(list_f[i])
  
})
