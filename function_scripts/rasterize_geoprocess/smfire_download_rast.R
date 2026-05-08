# function to extract rasters for previous small fires 
smfire_download_rast <- function(
    data_directory,
    output_directory,
    the.prj
    
){
  
  # Import CSV of selected fire names and Ignition Dates
  fires <- read.csv(paste0(output_directory,"/fire_lists/selectedfires.csv"))
  
  
  #Subset to only include fires that have small wildfires
  all_fp <- list.files(path = (paste0(output_directory,"/smfires/csv")), all.files = T, full.names = F, recursive = T)
  all_fp <-sub("_smfires.csv","",all_fp)
  all_fp=as.data.frame(all_fp)
  all_fp=all_fp%>%
    rename(MTBS_ID=all_fp)
  fires <- fires[fires$MTBS_ID %in% all_fp$MTBS_ID,]
  
  # Create list of fire names from CSV
  list_f <- fires$MTBS_ID
  list_f = gsub(' ','_', list_f)
  
  # Function to clip FACTS + NFPORS treatment layers to fire perimeters, 
  # output shapefile and csv of clipped treatment data
  ClipSmFires2Fire <- function(this_fire){
    
    #Pull in vector of fires
    smfires.clipped<-vect(paste0(output_directory,'/smfires/shapefile/',this_fire,'_smfires.shp'))
    
    #Name and file path
    this_fp <- paste0(output_directory,'/smfires')
    
    #Create a raster too####
    # create necessary directories for desired files paths
    if(dir.exists(paste0(this_fp,'/raster')) == FALSE){
      dir.create(paste0(this_fp,'/raster'))} 
    
    # create file path inside fire folder
    this_rast.fp <- file.path(this_fp,'raster',(paste0(this_fire,"_smfires.tif")))
    
    #Assign ID to polygons
    smfires.clipped$ID <- 1:nrow(smfires.clipped)
    
    #Create blank raster
    this_rast <- rast(extent = ext(smfires.clipped), res = 5, crs = crs(smfires.clipped))
    
    #Extract attributes we want from vector
    smfires.clipped <- smfires.clipped %>%
      tidyterra::select(FIREYEAR, ID)
    
    #Rasterize
    rasterize_extract <- function(smf_i) {
      this_r <- rasterize(smf_i, this_rast, field = "ID", fun = max)
      this_r <- crop(this_r, smf_i)
      xy <- as.data.frame(xyFromCell(this_r, 1:ncell(this_r)))
      extracted <- terra::extract(this_r, xy, ID = FALSE)
      data.table(x = xy[,1], y = xy[,2], extracted)
    }
    
    raster_coords <- rbindlist(lapply(1:nrow(smfires.clipped), function(p) {
      this_smf_i <- smfires.clipped[p, ]
      rasterize_extract(this_smf_i)
    }), use.names = TRUE, fill = TRUE)
    
    raster_coords <- na.omit(raster_coords)
    
    df_new <- raster_coords %>%
      group_by(x, y) %>%
      summarize(Combined = paste(unique(ID), collapse = "+"))
    
    overlapped_pixels <- df_new %>%
      filter(str_detect(Combined, "\\+"))
    
    # continue with data manipulation if there is overlapping pixels
    if (nrow(overlapped_pixels) > 0) {
      # Max number of overlapping pixels
      max_splits <- max(str_count(overlapped_pixels$Combined, '\\+')) + 1
      
      # Columns Names for UCKN
      col_names <- paste0("ID", 1:max_splits)
      
      
      # Seperate strings into individual UCKNS
      overlapped_pixels_2 <- tidyr::separate(overlapped_pixels, Combined, into = col_names, sep ='\\+', extra = 'merge')
      
      # rbind with non_overlapping pixels
      nonoverlapped_pixels <- df_new %>%
        filter(!str_detect(Combined, "\\+"))
      
      # Colnames
      colnames(nonoverlapped_pixels)[3] <- "ID1"
      
      # Combined
      ID_pixels <- bind_rows(nonoverlapped_pixels, overlapped_pixels_2)
      
      #Grab csv info
      this_csv<-as.data.frame(smfires.clipped)%>%
        select(FIREYEAR)
      this_csv$ID<-as.character(1:nrow(this_csv))
      
      # Join with pixel values 
      for (n in seq_along(col_names)){
        csv_clean_n <- this_csv
        colnames(csv_clean_n) <- paste0(colnames(csv_clean_n), n)
        ID_pixels <- left_join(ID_pixels, csv_clean_n, by = col_names[n])
      }
      
      # Get Max Year 
      ID_pixels$sf_max_yr <- do.call(pmax, c(ID_pixels[grep("FIREYEAR", names(ID_pixels))], na.rm = TRUE))
      
      # Count Overlaps
      # Get all unique ID Columns 
      SFCnt_columns <- grep("^ID", names(ID_pixels), value=T)
      # Convert these columns to character type
      ID_pixels[SFCnt_columns] <- lapply(ID_pixels[SFCnt_columns], as.character)
      # Instead of NA put in blanks 
      ID_pixels[,SFCnt_columns][is.na(ID_pixels[,SFCnt_columns])] <- ''
      #Get a count of overlapping fires
      ID_pixels$sf_cnt <- apply(ID_pixels[,SFCnt_columns], 1, function(x) length(unique(unlist(strsplit(as.character(x), '\\s*\\+\\s*')))))
      
      # select desired columns
      ID_pixels_xyz <- ID_pixels[,c('x','y','sf_max_yr', 'sf_cnt')]
      
      #Grab csv info
      this_csv<-as.data.frame(smfires.clipped)%>%
        select(FIREYEAR)
      this_csv$ID<-as.character(1:nrow(this_csv))
      
      
      # Join with pixel values
      for (n in seq_along(col_names)){
        csv_clean_n <- this_csv
        colnames(csv_clean_n) <- paste0(colnames(csv_clean_n), n)
        ID_pixels <- left_join(ID_pixels, csv_clean_n, by = col_names[n])
      }
      
      # Create Raster this is probably the slowest step!
      sf.comb <- rast(ID_pixels_xyz, type = "xyz", digits=0, crs=(the.prj))
      
      # write raster for individual fires 
      writeRaster(sf.comb[['sf_max_yr']], file.path(this_fp,'raster',(paste0(this_fire,"_sf_max_yr.tif"))), overwrite = T)
      writeRaster(sf.comb[['sf_cnt']], file.path(this_fp,'raster',(paste0(this_fire,"_sf_cnt.tif"))), overwrite = T)
      
      # Print complete statement
      print(paste0('Overlapping small fires found. Raster complete for fire ', this_fire))
    }else{
      
      # Edit attributes to match raster format
      this_csv=as.data.frame(smfires.clipped)
      this_csv$sf_cnt=1
      this_csv$ID <- as.character(this_csv$ID)
      
      counts<-this_csv%>%
        select(ID,FIREYEAR,sf_cnt)
      
      #merge csv columns to shp to add counts in
      smfires.clipped <- merge(smfires.clipped, counts)
      
      
      # Rasterize ID, Final Class, and max year from treatment shape
      smfires.clipped$sf_max_yr <- values(smfires.clipped[,'FIREYEAR'])
      smfires.clipped$sf_cnt <- values(smfires.clipped[,'sf_cnt'])
      
      rast_max <- rasterize(smfires.clipped, this_rast, field = 'sf_max_yr', fun = max)
      rast_cnt <- rasterize(smfires.clipped, this_rast, field = 'sf_cnt', fun = max)
      
      # write raster for individual fires 
      writeRaster(rast_max, file.path(this_fp,'raster',(paste0(this_fire,"_sf_max_yr.tif"))), overwrite = T)
      writeRaster(rast_cnt, file.path(this_fp,'raster',(paste0(this_fire,"_sf_cnt.tif"))), overwrite = T)
      
      # Print complete statement
      print(paste0('No overlapping fires found. Raster complete for fire ', this_fire))
      
    }
    
    gc()
    tmpFiles(current=TRUE, orphan=TRUE, old=TRUE, remove=TRUE)
    
  }
  
  
  ClippedFires<-lapply(1:length(list_f), function(i) {
    ClipSmFires2Fire(list_f[i])
    
  })
  
  
}