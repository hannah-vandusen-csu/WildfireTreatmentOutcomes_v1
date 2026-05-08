##########################################################.
###                                                    ###.
### Script to create treatment and previous fire masks ###.
###                                                    ###.
##########################################################.

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



# create treatment masks ---------

# define filepaths
input_dir_treatments <- file.path(processed_data_directory, "treatments/shapefile")
output_dir_treatments <- file.path(processed_data_directory, "masks/treatment")

# inputs for mask creation
res <- 30            # raster resolution in meters
buffer_m <- 30      # buffer distance in meters for treatments

# get fire perimeters
shp_files <- list.files(input_dir_treatments, pattern = "\\_cleaned.shp$", full.names = TRUE)

# create directory for treatment masks
dir.create(output_dir_treatments,recursive = T ,showWarnings = T)

# process each fire and make treatment mask
for (file in shp_files) {
  
  # Extract MTBS ID (everything before first underscore)
  mtbs_id <- str_extract(basename(file), "^[^_]+")
  
  # Read shapefile
  v <- vect(file)
  
  # Check if CRS is geographic (degrees)
  if (is.lonlat(v)) {
    message("Reprojecting ", mtbs_id, " to projected CRS for buffering...")
    
    # Automatically choose UTM zone based on centroid
    centroid <- geom(v, "centroid")[1, c("x", "y")]
    lon <- centroid[1]; lat <- centroid[2]
    utm_zone <- floor((lon + 180) / 6) + 1
    utm_crs <- paste0("EPSG:", 32600 + utm_zone)   # WGS84 UTM North
    
    v <- project(v, utm_crs)
  }
  
  # buffer 30 meters
  v_buff <- terra::buffer(v, width = buffer_m)
  
  # Create raster template for this buffered shape
  r_template <- rast(v_buff, resolution = res)
  
  # Rasterize buffered polygon (1 inside, NA outside)
  r_mask <- rasterize(v_buff, r_template, field = 1, touches = TRUE)
  
  # Convert NA to 0 for binary mask
  r_mask[is.na(r_mask)] <- 0
  
  # Define output path
  output_file <- file.path(output_dir_treatments, paste0(mtbs_id, "_treatment.tif"))
  
  # Write raster
  writeRaster(r_mask, output_file, overwrite = TRUE)
  
  message("Created treatment mask for: ", mtbs_id)
}

message("All treatment masks created and saved to: ", output_dir_treatments)


# create previous fire masks ----------

# define filepaths for creating masks of previous fires
input_dir_previousfire <- file.path(processed_data_directory, "smfires/raster")
output_dir_previousfire <- file.path(processed_data_directory, "masks/wildfire/")

# create a directory in which to save previous fire masks
dir.create(output_dir_previousfire, showWarnings = F)

res <- 30 # raster resolution in meters

# get filepaths of previous fires
files <- list.files(input_dir_previousfire, pattern = "\\_sf_cnt.tif$", full.names = TRUE)

# create previous fire masks for all fires
for (file in files) {
  
  # Extract MTBS ID (everything before first underscore)
  mtbs_id <- str_extract(basename(file), "^[^_]+")
  
  # Read shapefile
  r <- rast(file)
  
  # create new raster and write values
  out_r <- r
  
  # write output
  out_name <- paste0(mtbs_id, "_wildfire.tif")
  out_path <- file.path(output_dir_previousfire, out_name)
  
  # Create a binary clas
  out_r[out_r > 0] <- 1
  out_r[is.nan(out_r)] <- 0
  
  writeRaster(out_r, out_path, overwrite = TRUE)
  
  message(paste("Written:", out_path))
}

