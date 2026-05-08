####################################
###
### Step 1: Treatment Download
### The script acquires forest activity data from FACTS and NFPORS geodatabases for each Fire Perimeter
### Exports a shapefile and csv of activity attributes for each Fire Perimeter
###
### Written by Deborah Nemens May 2023 (dnemens@uw.edu)
### Edited by Hannah Van Dusen January 2024 (hannah.vandusen@usda.gov)
### Adapted for FACTS data by Betsy Black February 2024 (elizabeth.black@usda.gov)
###
####################################


# Set working Directory to Current R project location

# Get directory where scripts are saved
# current_directory <- dirname(rstudioapi::getActiveDocumentContext()$path)

# Navigate to parent and project directory
# parent_directory <- dirname(dirname(dirname(current_directory)))
# setwd(parent_directory)
# data_directory <- parent_directory



# Function to Load global datasets and run each fire 
treatment_download <- function(data_directory,
                               output_directory) {


  # Import CSV of names and Ignition Dates
  fires <- read.csv(paste0(output_directory,"/fire_lists/selectedfires_raw.csv"))
  
  # Create list of fire names from CSV
  list_f <- fires$MTBS_ID
  
  # List of filepaths for each fire's perimeter
  fire_fp <- list.files(path = (paste0( output_directory, '/fire_perimeters/merged')), all.files = T, full.names = T,
                        recursive = T, pattern = "Extent.shp")
  fire_fp <- fire_fp[gsub('_Extent.shp', '', gsub(paste0(output_directory, '/fire_perimeters/merged/'), '', fire_fp)) %in% list_f]
  
  # Pull in perimeter shapefiles for fires of interest
  extents <- lapply(fire_fp, FUN = vect)
  study_ext <- ext(vect(extents))

  # Read in FACTS geodatabases
  # FACTS databases are saved by region, so we must figure which region to read in since
  # FACTS is such a large database
  
  # Forest service regions
  usfs_regions <- terra::vect(paste0(data_directory, '/data/global_data/geopolitical_units/S_USA.AdministrativeRegion/S_USA.AdministrativeRegion.shp'))
  usfs_regions <- terra::project(usfs_regions, crs(vect(extents)))
  usfs_regions <- terra::crop(usfs_regions, study_ext)
  
  usfs_Rcode <- usfs_regions$REGION
  
  facts <- lapply(usfs_Rcode, function(x){
    vect(paste0(data_directory, "/data/global_data/treatments/facts/Actv_CommonAttribute_PL_Region", x ,".gdb"), proxy = T)})
  
  # Read in NFPORS treatments
  nfpors <- file.path(paste0(data_directory, "/data/global_data/treatments/nfpors/fuelperims")) |>
    vect(proxy = T)
  
  # Function to clip FACTS + NFPORS treatment layers to fire perimeters,
  # output shapefile and csv of clipped treatment data
  ClipTreatment2Fire <- function(extent) {
    
    # add a 2000 m buffer to the fire perimeter for infrastructure needs
    region <- vect(intersect(usfs_regions, extent))$REGION
    extent <- project(extent, crs(facts[[1]]))
    extent <- terra::buffer(extent, 1000)
    
    
    # Clip Treatment Layer to Fire perimeter
    facts_vect <- lapply(facts, function(x){
      facts_region <- str_extract(sources(x), "(?<=Region)\\d{2}(?=\\.gdb)")
      if(facts_region %in% region){
        facts1 <- query(x, vars =  c("SUID", "GIS_ACRES", "ACTIVITY", "DATE_COMPLETED", "UKCN"), extent = ext(extent))
        return(facts1)
      }else{
        return(NULL)
      }
    })
    
    facts_vect <- facts_vect[!lengths(facts_vect) == 0]
    # no treatments are found return NULL
    if(length(facts_vect) == 0){
      print(paste0('No Treatments found for ', extent$MTBS_ID, ' fire.'))
      return(NULL)
    }
    treats_f <- vect(facts_vect)
    
    # Make sure all time stamps are 10 characters
    treats_f$DATE_COMPLETED <- (str_sub(treats_f$DATE_COMPLETED, 1, 10))
    
    if (nrow(treats_f) > 0) {
      treats_f <- treats_f %>%
        mutate(
          cent_lat = as.numeric(NA),
          cent_lon = as.numeric(NA),
          nfporsfid = as.numeric(NA),
          type_name = as.character(NA),
          cat_nm = as.character(NA)
        )
    }
    
    # Read in NFPORS
    extent <- project(extent, crs(nfpors))
    nfpors <- query(nfpors, vars =  c("objectid", "gis_acres", "cent_lat", "cent_lon", "act_comp_d", "nfporsfid", "type_name", "cat_nm"), extent = ext(extent))
    treats_n <- project(nfpors, crs(treats_f))
    
    # rename colummns to match FACTS
    if (nrow(treats_n) > 0) {
      treats_n <- treats_n %>%
        rename(
          UKCN = objectid,
          DATE_COMPLETED = act_comp_d,
          GIS_ACRES = gis_acres) %>%
        mutate(
          cent_lat = as.numeric(NA),
          cent_lon = as.numeric(NA),
          nfporsfid = as.numeric(NA),
          type_name = as.character(NA),
          cat_nm = as.character(NA)
        )
    }
    
    # Combine facts and nfpors treatment vectors
    treats_clipped <- rbind(treats_f, treats_n)
    
    # Naming conventions
    names(treats_clipped) <- sub("_", "", sub("_", "", names(treats_clipped))) # remove spacing
    
    # create a csv of treatment attributes
    this_csv <- as.data.frame(treats_clipped)
    
    
    if(dim(this_csv)[1] > 0) {
      
      # Naming conventions
      this_nm <- extent$MTBS_ID[1]
      this_fp <- paste0(output_directory, '/treatments')
      
      # Create Treatments Folder
      if(!dir.exists(paste0(output_directory, '/treatments'))){
        dir.create(paste0(output_directory, '/treatments'))
      }
      
      
      # create file path inside fire folder
      this_shp_fp <- (paste0(this_fp, '/shapefile/'))
      this_shp_nm <- paste0(this_nm, "_treatments_raw.shp")
      
      # create necessary directories for desired files paths
      if(dir.exists(paste0(this_fp, '/shapefile')) == FALSE) {
        dir.create(paste0(this_fp, '/shapefile'))
      }
      
      # output shapefile of treatments for each fire
      # Write Shape file
      writeVector(treats_clipped, paste0(this_shp_fp, this_shp_nm), overwrite = T)
      
      
      
      # create necessary directories for desired files paths\
      if(dir.exists(paste0(this_fp, '/csv')) == FALSE) {
        dir.create(paste0(this_fp, '/csv'))
      }
      
      #creates filepath for csv inside fire folder
      this_csv_path <- (paste0(this_fp, '/csv/'))
      this_csv_nm <- (paste0(this_nm, "_treat_attributes_raw.csv"))
      
      #output csv of treatment attributes
      write.csv(this_csv, paste0(this_csv_path, this_csv_nm), row.names = FALSE)
      print(paste0('Downloading Treatments of ', this_nm, ' fire complete'))
      
      return(this_nm)
      
    } else {
      
      print(paste0('No Treatments found for ', extent$MTBS_ID, ' fire.'))
      
      return(NULL)
    }
    
  }

  # Apply Function to Extents
  ClippedTreatments <- lapply(extents, ClipTreatment2Fire)

  # Return a New Selected Fire Shape file with fires with Successfully Downloaded Treatments
  # Select fire extents with viable clipped treatments
  completed_fires <- unlist(ClippedTreatments[lengths(ClippedTreatments) != 0])

  # Subset Selected Fires to those Successfully Downloaded
  fires$step1_retention <- fires$MTBS_ID %in% completed_fires

  # Save Selected Fires as CSV
  write.csv(fires, paste0(output_directory, '/fire_lists/selectedfires_raw.csv'), row.names=FALSE)

}






