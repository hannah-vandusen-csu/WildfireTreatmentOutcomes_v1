# function to create blank rasters for fires without previous small fires

smfire_download_dummy_rast <- function(
  data_directory,
  output_directory,
  the.prj
){
  
  # Import CSV of selected fire names and Ignition Dates
  fires <- read.csv(paste0(output_directory,"/fire_lists/selectedfires.csv"))
  
  #load which fires have small wildfires
  all_fp <- list.files(path = (paste0(output_directory,"/smfires/csv")), all.files = T, full.names = F, recursive = T)
  all_fp <-sub("_smfires.csv","",all_fp)
  all_fp=as.data.frame(all_fp)
  all_fp=all_fp%>%
    rename(MTBS_ID=all_fp)
  
  # subset to only fires that don't have small wildfires
  fires <- fires[!fires$MTBS_ID %in% all_fp$MTBS_ID,]
  
  # Create list of fire names from CSV
  list_f <- fires$MTBS_ID
  list_f = gsub(' ','_', list_f)
  
  # Function to clip FACTS + NFPORS treatment layers to fire perimeters, 
  # output shapefile and csv of clipped treatment data
  ClipSmFires2Fire <- function(this_fire){
    
    #Pull in fire perimeter
    perim<-vect(paste0(output_directory,"/fire_perimeters/merged/",this_fire,"_Extent.shp")) %>%
      project(the.prj)
    
    #Create blank raster
    this_rast <- rast(extent = ext(perim), res = 5, crs = crs(the.prj))
    
    sf_cnt<-this_rast
    sf_cnt$sf_cnt<-0
    sf_cnt<-mask(sf_cnt,perim)
    
    
    sf_max_yr<-this_rast
    sf_max_yr$sf_max_yr<-0
    sf_max_yr<-mask(sf_max_yr,perim)
    sf_max_yr[sf_max_yr$sf_max_yr==0]<-NA
    
    
    # write raster for individual fires 
    writeRaster(sf_max_yr, paste0(output_directory,'/smfires/raster/',this_fire,"_sf_max_yr.tif"), overwrite = T)
    writeRaster(sf_cnt, paste0(output_directory,'/smfires/raster/',this_fire,"_sf_cnt.tif"), overwrite = T)
  }
  
  
  ClippedFires<-lapply(1:length(list_f), function(i) {
    ClipSmFires2Fire(list_f[i])
    
  })
  
  
}