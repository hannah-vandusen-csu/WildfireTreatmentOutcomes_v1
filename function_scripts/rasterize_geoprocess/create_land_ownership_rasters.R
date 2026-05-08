# function to create land ownership rasters
create_land_ownership_rasters <- function(output_directory){
  
  
  # Import CSV of names and Ignition Dates
  fires <- read.csv(paste0(output_directory,"/fire_lists/selectedfires.csv"))
  
  # Create list of fire names from CSV
  list_f <- fires$MTBS_ID
  list_f = gsub(' ','_', list_f)
  
  
  # create necessary directories for desired files paths
  if(dir.exists(paste0(output_directory,'/masks/land_ownership/rasters')) == FALSE){
    dir.create(paste0(output_directory,'/masks/land_ownership/rasters'))
  } 
  
  #Specify agency acronyms we have TX data for
  public_lands <- c("USFS","FWS","BLM", "NPS","USBR") #Excluded BIA because there is no data for BIA in our landscape
  
  # Create Public Land (1) or not (0) rasters
  ClipLO2Fire <- function(this_fire){
    
    #START A FIRE####
    print(paste0("Processing ",this_fire))
    
    # Pull in shapefile
    ownership<-vect(paste0(output_directory,'/masks/land_ownership/shapes/',this_fire,'_land_own.shp'))
    
    #Assign 1 if public land 0 if not
    ownership$pub_lnd <- ifelse(ownership$ADMIN_AGE0 %in% public_lands, 1, 0)
    
    
    # Create an empty raster with the desired extent and resolution
    r_template <- rast(ownership, resolution = 30)  # Adjust resolution as needed
    
    # Rasterize: Assign 1 if ADMIN_AGE0 is in public_lands, 0 otherwise
    pub_lnd <- rasterize(ownership, r_template, field = "pub_lnd", 
                         fun = "max")
    
    # Set NA values to 0 (optional, to make all non-matching areas explicitly 0)
    pub_lnd[is.na(pub_lnd)] <- 0
    
    writeRaster(pub_lnd,paste0(output_directory,'/masks/land_ownership/rasters/',this_fire,'_pub_lnd.tif'), overwrite = T)
    
    print(paste0('public land raster for ',this_fire,' created'))
    
  }
  
  # Iterate over each fire and wui piece
  for (i in seq_along(list_f)) {
    ClipLO2Fire(list_f[i])
  }
}