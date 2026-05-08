##############################################.
##                                          ##.
## Create sampling grid from burn severity  ##. 
##                                          ##.
##############################################.



# load packages -----

# Necessary packages
pkgs <- c("terra",
          "tidyverse",
          "rstudioapi",
          "sf",
          "stringr"
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


# set other inputs --------

# Set outcome 
wto_outcome <- 'burn_severity'

# Read in fire list
fire_list <- read.csv(paste0(processed_data_directory, "/fire_lists/selectedfires.csv"))
fire_list$year <- lubridate::year(as.Date(fire_list$Ignition, format = '%m/%d/%Y'))

# Filter to 2012
fire_list <- fire_list %>% filter(year >= 2012)

# Burn severity fp
burn_fp <- '/fire_severity/RBR'

# resolution of burn severity (30 m is default)
sev_res <- 30 

# desired sampling resolution
resol <- 270 

# read in RAVG table
RAVG_basal_area_lookup_table <- file.path(global_data_directory, "RAVG_lookup_table/SWModel_Lookup_BA100_EA.txt") %>%
  read.table(header = T) %>%
  as.data.frame()
names(RAVG_basal_area_lookup_table) <- names(RAVG_basal_area_lookup_table) %>% tolower()

# get high severity threshold (based on user input "basal_area_high_severity_percent_threshold" in "user_inputs_general.csv")
high_severity_threshold <- RAVG_basal_area_lookup_table %>% 
  filter(meanba == basal_area_high_severity_percent_threshold) %>%
  reframe(threshold = min(.data[[tolower(burn_severity_metric)]], na.rm = TRUE)) %>%
  pull(threshold)


# Read in burn severity files 
fils <- list.files(paste0(processed_data_directory, burn_fp),full.names = T, pattern = paste0(burn_severity_metric, ".tif"))

# Keep only files ending with ".tif"
fils <- fils[grepl("\\.tif$", fils)]

# Filter files to only those that match the fire list
fils <- fils[grepl(paste0(fire_list$MTBS_ID, collapse="|"), fils)]

# Read in rasters
rasts <- lapply(fils, rast)

# Loop through each raster and create points at the desired resolution
for(i in 1:length(rasts)) {
  # subset to the current raster
  this_rast <- rasts[[i]]
  
  # Aggregate to desired resolution
  agg_factor <- resol / sev_res
  this_rast <- aggregate(this_rast, agg_factor, fun = mean)
  
  # Create points from raster
  these_pnts <- st_as_sf(as.points(this_rast))
  
  # Extract the fire ID from the file name
  id <- basename(fils[[i]])
  id <- str_extract(id, "[^_]+")
  
  # load fire perimeter
  if(wto_outcome == "burn_severity"){
    these_pnts_v <- these_pnts %>% vect() %>% project(CRS)
    fire_perim <- vect(paste0(processed_data_directory, "/fire_perimeters/merged/", id, "_Extent.gpkg"))
    these_pnts_v <- these_pnts_v[fire_perim, ] # keep only points that fall within perimeter 
    these_pnts <- these_pnts_v %>% st_as_sf() %>% st_transform("WGS84")
  }
  
  # Add MTBS_ID and point_id columns to the points
  these_pnts$MTBS_ID <- id
  these_pnts$point_id <- 1:nrow(these_pnts)
  
  # convert to data frame and add coordinates
  gcoords <- as.data.frame(st_coordinates(these_pnts))
  names(gcoords) <- c("long", "lat")
  these_pnts <- cbind(these_pnts, gcoords)
  
  # add UTM coordinates
  ucoords <- st_transform(these_pnts, CRS)
  ucoords <- as.data.frame(st_coordinates(ucoords))
  names(ucoords) <- c("X", "Y")
  these_pnts <- cbind(these_pnts, ucoords)
  
  # extract the burn severity values for each point
  these_pnts$response <- terra::extract(this_rast, these_pnts)[,2]
  
  # if points are less than burn severity threshold, set response to 0, else 1
  these_pnts$response_bin <- ifelse(these_pnts$response < high_severity_threshold, 0, 1)
  
  # select and rename columns
  st_geometry(these_pnts) <- NULL
  these_pnts <- dplyr::select(these_pnts, c("MTBS_ID", "point_id", "long", "lat", "X", "Y", "response", "response_bin"))


  # If the number of points is greater than 15000, 
  # save data frames into smaller chunks of ~8000 pts 
  if (nrow(these_pnts) > 14000) {
    
    # Calculate the number of chunks needed
    num_chunks <- ceiling(nrow(these_pnts) / 8000)
    
    # Split the points into chunks
    chunk_size <- ceiling(nrow(these_pnts) / num_chunks)
    pnts_chunks <- split(these_pnts, ceiling(seq_len(nrow(these_pnts)) / chunk_size))
    
    # Create directory
    if(!dir.exists(paste0(processed_data_directory, "/initial_sampling_points/",wto_outcome, '/'))){
      dir.create(paste0(processed_data_directory, "/initial_sampling_points/",wto_outcome, '/'), recursive = TRUE)
    }
    
    # Save each chunk as a separate shapefile
    for (i in seq_along(pnts_chunks)) {
      output_fp <- paste0(processed_data_directory, "/initial_sampling_points/",wto_outcome, '/', id,"_", i, "_initial_pts.csv")
      write.csv(pnts_chunks[[i]], output_fp, row.names = F)
      print(paste0(i,'/', num_chunks, ' ', id))
    }
  } else {
    
    # Create directory
    if(!dir.exists(paste0(processed_data_directory, "/initial_sampling_points/",wto_outcome, '/'))){
      dir.create(paste0(processed_data_directory, "/initial_sampling_points/",wto_outcome, '/'), recursive = TRUE)
    }
    
    # Save the centroids as a single shapefile
    output_fp <- paste0(processed_data_directory, "/initial_sampling_points/",wto_outcome, '/', id, "_initial_pts.csv")
    write.csv(these_pnts, output_fp, row.names = F)
    print(paste0('no chunking ', id))
  }

}

