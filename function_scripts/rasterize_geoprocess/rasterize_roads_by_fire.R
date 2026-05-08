# function to rasterize roads
rasterize_roads_by_fire <- function(output_directory,
                                    the.prj){
  # load in roads data again
  rds_raw <- vect(file.path(output_directory, "roads/roads_AOI.gpkg"))
  
  
  # Ensure output path exists
  dir <- paste0(output_directory, "/roads/roads_by_fire/")
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  dir <- paste0(output_directory, "/roads/roads_by_fire/shapes/")
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  dir <- paste0(output_directory, "/roads/roads_by_fire/rasters/")
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  
  
  # load list of fires
  fires <- read.csv(paste0(output_directory,"/fire_lists/selectedfires.csv"))
  list_f <- fires$MTBS_ID
  
  
  rasterize_roads <- function(this_fire){
    
    print(paste0("Processing ",this_fire))
    
    #Bring in fire perim
    perim<-vect(paste0(output_directory,"/fire_perimeters/merged/",this_fire,"_Extent.gpkg")) %>%
      project(crs(the.prj))
    
    #buffer
    per_buf<-buffer(perim,1000)
    
    #clip rds to buffer
    this_rds<-crop(rds_raw,per_buf)
    
    #Save containment shapes associated with this fire (using sf, as terra can't write if empty)
    st_write(st_as_sf(this_rds),paste0(output_directory,"/roads/roads_by_fire/shapes/",this_fire,"_rds.gpkg"),append = F)
    
    #Buffer rds
    this_rds<-buffer(this_rds,5)
    
    #add column for raster
    this_rds$rds<-1
    
    # Create an empty raster
    r_template <- rast(per_buf, resolution = 5)
    
    # Rasterize
    rds <- rasterize(this_rds, r_template, field = "rds",
                     fun = "max")
    
    # Set NA values to 0 (optional, to make all non-matching areas explicitly 0)
    rds[is.na(rds)] <- 0
    
    rds<-mask(rds,per_buf)
    
    writeRaster(rds,paste0(output_directory,"/roads/roads_by_fire/rasters/",this_fire,"_rds.tif"),overwrite=T)
    
    rm(perim,per_buf,this_rds,r_template,rds)
    
    gc()
    
    
    
  }
  
  lapply(1:length(list_f), function(i) {
    rasterize_roads(list_f[i])
    
  })
}