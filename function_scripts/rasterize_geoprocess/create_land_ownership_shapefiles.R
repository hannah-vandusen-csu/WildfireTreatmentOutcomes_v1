# function to extract land ownership polygons for each fire
create_land_ownership_shapefiles <- function(data_directory,
                                             output_directory,
                                             the.prj){
  
  
  #Load land ownership gdb####
  fp <- paste0(data_directory, "/land_ownership/SMA_WM.gdb")  
  own <- vect(fp, layer="SurfaceManagementAgency") %>%
    project(crs(the.prj))
  
  # Import CSV of names and Ignition Dates
  fires <- read.csv(paste0(output_directory,"/fire_lists/selectedfires.csv"))
  
  # Create list of fire names from CSV
  list_f <- fires$MTBS_ID
  list_f = gsub(' ','_', list_f)
  
  # create necessary directories for desired files paths
  if(dir.exists(paste0(output_directory,'/masks')) == FALSE){
    dir.create(paste0(output_directory,'/masks'))
  } 
  # create necessary directories for desired files paths
  if(dir.exists(paste0(output_directory,'/masks/land_ownership')) == FALSE){
    dir.create(paste0(output_directory,'/masks/land_ownership'))
  } 
  
  # create necessary directories for desired files paths
  if(dir.exists(paste0(output_directory,'/masks/land_ownership/shapes')) == FALSE){
    dir.create(paste0(output_directory,'/masks/land_ownership/shapes'))
  } 
  
  # Function to clip FACTS + NFPORS treatment layers to fire perimeters, 
  # output shapefile and csv of clipped treatment data
  ClipLO2Fire <- function(this_fire){
    
    #START A FIRE####
    print(paste0("Processing ",this_fire))
    
    # Load fire perimeter
    fire_perim <- vect(paste0(output_directory, '/fire_perimeters/merged/',this_fire,'_Extent.shp')) %>%
      project(crs(the.prj))
    
    #Add buffer to fire perim
    buffered_extent<-buffer(fire_perim,1000) #in meters bc epsg 26912 projection
    
    #Clip ownership to a 1000m buffer around fire extent
    own.clipped <- terra::crop(own, buffered_extent)
    
    ####
    
    # create file path inside folder for shapefiles
    own.fp <- (paste0(output_directory,'/masks/land_ownership/shapes/'))
    own.nm <- paste0(this_fire, "_land_own.gpkg")
    
    # Write Shape file
    writeVector(own.clipped, paste0(own.fp, own.nm), overwrite = T)
    print(paste0('land ownership for ',this_fire,' clipped'))
    
  }
  
  # Iterate over each fire and wui piece
  for (i in seq_along(list_f)) {
    ClipLO2Fire(list_f[i])
  }
}
