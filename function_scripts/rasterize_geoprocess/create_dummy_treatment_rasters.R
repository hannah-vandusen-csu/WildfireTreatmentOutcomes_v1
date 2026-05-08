# function to create blank rasters for treatments for fires without them
create_dummy_treatment_rasters <- function(
    data_directory,
    output_directory,
    the.prj){
  
  #Call in list of fires
  fires <- read.csv(paste0(output_directory,"/fire_lists/selectedfires.csv"))
  
  #Subset to only include fires that don't already have tx rasters
  all_fp_MTBS <- list.files(path = (paste0(output_directory,"/treatments/csv")), all.files = T, full.names = F, recursive = T, pattern = "_treatments_cleaned.csv") %>%
    sub("_treatments_cleaned.csv","",.)
  # find what fires don't already have treatment rasters
  fires <- fires[!fires$MTBS_ID %in% all_fp_MTBS]
  
  # Create list of fire names from CSV
  list_f <- fires$MTBS_ID
  
  if(length(list_f) > 0){
    
    ByFire <- function(this_fire){
      
      print(paste0("Processing ",this_fire))
      
      #Bring in fire perim
      perim<-vect(paste0(output_directory,"/fire_perimeters/merged/",this_fire,"_Extent.shp"))
      perim<-project(perim,crs(the.prj))
      
      #buffer
      per_buf<-terra::buffer(perim,1000)
      
      #Specify all layers as 0 (later replace w NAs where necessary)
      per_buf$class.cod<-0
      per_buf$tr_max_year<-0
      per_buf$uniqueID<-0
      per_buf$txYrs_cnt<-0
      
      # Create an empty raster
      r_template <- rast(per_buf, resolution = 5) 
      
      # Rasterize: Assign 1 if ADMIN_AGE0 is in public_lands, 0 otherwise
      class.cod <- rasterize(per_buf, r_template, field = "class.cod", 
                             fun = "max")
      
      class.cod<-mask(class.cod,per_buf)
      class.cod[class.cod$class.cod==0]<-NA
      
      
      # Rasterize: Assign 1 if ADMIN_AGE0 is in public_lands, 0 otherwise
      tr_max_year <- rasterize(per_buf, r_template, field = "tr_max_year", 
                               fun = "max")
      
      tr_max_year<-mask(tr_max_year,per_buf)
      tr_max_year[tr_max_year$tr_max_year==0]<-NA
      
      
      # Rasterize: Assign 1 if ADMIN_AGE0 is in public_lands, 0 otherwise
      uniqueID <- rasterize(per_buf, r_template, field = "uniqueID", 
                            fun = "max")
      
      uniqueID<-mask(uniqueID,per_buf)
      uniqueID[uniqueID$uniqueID==0]<-NA
      
      
      # Rasterize: Assign 1 if ADMIN_AGE0 is in public_lands, 0 otherwise
      txYrs_cnt <- rasterize(per_buf, r_template, field = "txYrs_cnt", 
                             fun = "max")
      
      txYrs_cnt<-mask(txYrs_cnt,per_buf)
      
      txs_classed<-c(class.cod,tr_max_year)
      txs_classed<-c(txs_classed,uniqueID)
      txs_classed<-c(txs_classed,txYrs_cnt)
      
      writeRaster(txs_classed,paste0(output_directory,"/treatments/rasters/",this_fire,"_treatments_classed.tif"),overwrite=T)
      
    }
    
    DoThatThang<-lapply(1:length(list_f), function(i) {
      ByFire(list_f[i])
      
    })
  }
  
}