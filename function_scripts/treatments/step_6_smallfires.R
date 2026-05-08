####################################
###
### Step 2.5: Previous Small Wildfire Download
### The script acquires previous small wildfire data from USGS fire geodatabase for each Fire Perimeter 
### Exports a shapefile and csv of wildfire attributes for each Fire Perimeter
###
### Written by Betsy Black April 2024 (elizabeth.black@usda.gov)
###
####################################

smfire_download <- function(
    
  data_directory,
  output_directory,
  the_prj
  
){
  # Import CSV of names and Ignition Dates

  # fires <- read.csv(paste0(output_directory,"/fire_lists/selectedfires_list_v5.csv"))
  fires <- read.csv(paste0(output_directory,"/fire_lists/selectedfires_raw.csv"))

  
  # Create list of fire names from CSV
  list_f <- fires$MTBS_ID
  list_f = gsub(' ','_', list_f)
  
  # Create a list of ignition dates from CSV
  ign_f <- as.Date(fires$Ignition)
  
  #Load in other wildfires
  smfires_shp <- file.path(paste0(data_directory,"/fire_perimeters/usgs_small_fires/InterAgencyFirePerimeterHistory_All_Years_View2"))
  smfires <- vect(smfires_shp)
  smfires <- project(smfires, the_prj)
  
  # Function to clip FACTS + NFPORS treatment layers to fire perimeters, 
  # output shapefile and csv of clipped treatment data
  ClipSmFires2Fire <- function(this_fire, this_ign){
    
    #Get fire year
    this_year = as.numeric(year(this_ign))
    
    # Load fire extent
    extent<- vect(paste0(output_directory, "/fire_perimeters/merged/", this_fire, "_Extent.shp"))
    extent <- project(extent, the_prj)
    
    #Crop small wildfires spatvector to a 1000m buffer around fire extent
    buffer_distance <- 1000
    buffered_extent <- buffer(extent, width = buffer_distance)
    smfires.clipped <- terra::crop(smfires, buffered_extent)
    
    if (nrow(smfires.clipped) >0){
      
      ####
      smfires.clipped$FIRE_YEAR_<-NULL
      names(smfires.clipped) <- sub("_","", sub("_", "", names(smfires.clipped))) # remove spacing
      
      ####
      
      #Remove fires outside temporal ROI and prescribed fires
      lower_bound <- this_year - 21 # 20 years previous
      
      smfires.clipped <- smfires.clipped%>%
        filter(FEATURECA == "Prescribed Fire"|FEATURECA == "Wildfire Final Perimeter"|FEATURECA == "Wildfire")%>% #Checked the 11 instances in which there were "prescribed fires" downloaded in the KT and none overlapped with RX fires in cleaned treatment data
        filter(FIREYEAR < this_year)%>% # must ocurr before fire
        filter(FIREYEAR >= lower_bound)
      
      if (nrow(smfires.clipped) >0){
        
        #Filter out duplicate fires
        smfires.clipped <- smfires.clipped%>%
          mutate(INCIDENT =str_replace_all(tolower(INCIDENT),"\\s+",""))
        
        smfires.clipped=unique(smfires.clipped)
        smfires.clipped <- smfires.clipped%>%
          mutate(FEATURECA=factor(FEATURECA,levels=c("Wildfire Final Perimeter","Wildfire","Prescribed Fire")))%>%
          arrange(INCIDENT, FIREYEAR, FEATURECA)%>%
          group_by(INCIDENT, FIREYEAR)%>%
          slice(1)%>%
          ungroup()
        
        #Clean up columns
        smfires.clipped<- smfires.clipped%>%
          tidyterra::select(INCIDENT,FIREYEAR,FEATURECA,SOURCE,AGENCY,UNITID,IRWINID,COMMENTS,GISACRES,ShapeAre,ShapeLen)
        
        #Extract CSV of small fires
        this_csv <- as.data.frame(smfires.clipped)
        
        #Name and file path
        this_fp <- paste0(output_directory,'/smfires')
        
        # Create necessary directories
        if(!dir.exists(paste0(output_directory,'/smfires'))){
          dir.create(paste0(output_directory,'/smfires'))}
        if(dir.exists(paste0(this_fp,'/shapefile')) == FALSE){
          dir.create(paste0(this_fp,'/shapefile'))} 
        if(dir.exists(paste0(this_fp,'/csv')) == FALSE){
          dir.create(paste0(this_fp,'/csv'))}
        
        # Write shapefile
        this_shp.fp <- file.path(this_fp,'shapefile',paste0(this_fire,"_smfires.shp"))
        writeVector(smfires.clipped, this_shp.fp, overwrite = T)
        
        # Write csv
        this_csv.fp <- file.path(this_fp,'csv',paste0(this_fire,"_smfires.csv"))
        if (file.exists(this_csv.fp)) { file.remove(this_csv.fp) }
        write.csv(this_csv, this_csv.fp)
        
        print(paste0('Downloading Small fire shapes for ',this_fire, ' fire complete'))
        
        return(TRUE)
      }else{
        
        print(paste0('No Small fires found for ',this_fire))
        
        return(FALSE)
      }
      
    }else{
      
      print(paste0('No Small fires found for ',this_fire))
      
      return(FALSE)
    }
    
  }
  
  ClippedFires<-lapply(1:length(list_f), function(i) {
    ClipSmFires2Fire(list_f[i],ign_f[i])
    
  })
  
}







#RASTER EXTRACTION####

smfire_download_rast <- function(

  data_directory,
  output_directory ,
  the_prj

){

  # Import CSV of selected fire names and Ignition Dates
  # fires <- read.csv(paste0(output_directory,"/fire_lists/selectedfires_list_v5.csv"))
  fires <- read.csv(paste0(output_directory,"/fire_lists/selectedfires.csv"))


  #Subset to only include fires that have small wildfires
  all_fp <- list.files(path = (paste0(output_directory,"/smfires/csv")), all.files = T, full.names = F, recursive = T)
  all_fp <-sub("_smfires.csv","",all_fp)
  all_fp=as.data.frame(all_fp)
  all_fp=all_fp%>%
    rename(MTBS_ID=all_fp)
  fires <- fires[fires$MTBS_ID %in% all_fp$MTBS_ID,]
  #fires=left_join(all_fp,fires)

  # Create list of fire names from CSV
  list_f <- fires$MTBS_ID
  list_f = gsub(' ','_', list_f)

  # Create a list of ignition dates from CSV
  # ign_f <- as.Date(fires$Ignition)

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
        dplyr::select(FIREYEAR)
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
      #ID_pixels_xyz <-

      #Grab csv info
      this_csv<-as.data.frame(smfires.clipped)%>%
        dplyr::select(FIREYEAR)
      this_csv$ID<-as.character(1:nrow(this_csv))


      # Join with pixel values --??
      for (n in seq_along(col_names)){
        csv_clean_n <- this_csv
        colnames(csv_clean_n) <- paste0(colnames(csv_clean_n), n)
        ID_pixels <- left_join(ID_pixels, csv_clean_n, by = col_names[n])
      }

      # Create Raster this is probably the slowest step!
      sf.comb <- rast(ID_pixels_xyz, type = "xyz", digits=0, crs=(the_prj))

      # write raster for individual fires
      writeRaster(sf.comb[['sf_max_yr']], file.path(this_fp,'raster',(paste0(this_fire,"_sf_max_yr.tif"))), overwrite = T)
      writeRaster(sf.comb[['sf_cnt']], file.path(this_fp,'raster',(paste0(this_fire,"_sf_cnt.tif"))), overwrite = T)

      # Write Raster
      # writeRaster(sf.comb, this_rast.fp, overwrite = T)

      # Print complete statement
      print(paste0('Overlapping small fires found. Raster complete for fire ', this_fire))
    }else{

      # Edit attributes to match raster format
      this_csv=as.data.frame(smfires.clipped)
      this_csv$sf_cnt=1
      # this_csv$ID<-1:nrow(this_csv)
      this_csv$ID <- as.character(this_csv$ID)

      counts<-this_csv%>%
        dplyr::select(ID,FIREYEAR,sf_cnt)

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

      # sf.comb <- c(rast_max, rast_cnt)
      # crs(sf.comb) <- the_prj
      # sf.comb <- crop(sf.comb, smfires.clipped)

      # Write Raster
      # writeRaster(sf.comb, this_rast.fp, overwrite = T)

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
# smfire_download_rast(data_directory)
# 
# 
# 
# 
# 
# 
# #DUMMY RASTER EXTRACTION####
# 
# # library(terra)
# # tmpFiles(current=TRUE, orphan=FALSE, old=TRUE, remove=TRUE)
# rm(list=ls())
# 
# # Load Libraries 
# # require(sf)
# # sf::sf_use_s2(FALSE)#this switches off spherical geom, removes errors from facts layers not loading
# library(tidyverse)
# library(terra)
# library(tidyterra)
# library(data.table)
# # library(parallel)
# 
# # Get data directories
# data_directory <- "C:/Users/ElizabethBlack/Box/External Wildfire Treatment Outcomes"
# setwd(data_directory)
# output_directory = paste0(data_directory, '/data/interim_data')
# 
# smfire_download_rast <- function(
#     
#   data_directory,
#   output_directory = paste0(parent_directory, '/data/interim_data'),
#   the_prj = crs('epsg:26910')
#   
# ){
#   
#   # Import CSV of selected fire names and Ignition Dates
#   # fires <- read.csv(paste0(output_directory,"/fire_lists/selectedfires_list_v5.csv"))
#   fires <- read.csv(paste0(output_directory,"/fire_lists/fire_list_all_projects.csv"))
#   
#   #Subset to only include fires that have small wildfires
#   all_fp <- list.files(path = (paste0(data_directory,"/data/interim_data/smfires/csv")), all.files = T, full.names = F, recursive = T)
#   #all_fp <-sub("C:/Users/ElizabethBlack/OneDrive - USDA/WildfireTreatmentOutcomes/data/interim_data/smfires/csv/","",all_fp)
#   all_fp <-sub("_smfires.csv","",all_fp)
#   all_fp=as.data.frame(all_fp)
#   all_fp=all_fp%>%
#     rename(MTBS_ID=all_fp)
#   fires <- fires[!fires$MTBS_ID %in% all_fp$MTBS_ID,]
#   #fires=left_join(all_fp,fires)
#   
#   # Create list of fire names from CSV
#   list_f <- fires$MTBS_ID
#   list_f = gsub(' ','_', list_f)
#   
#   # Create a list of ignition dates from CSV
#   # ign_f <- as.Date(fires$Ignition)
#   
#   # Function to clip FACTS + NFPORS treatment layers to fire perimeters, 
#   # output shapefile and csv of clipped treatment data
#   ClipSmFires2Fire <- function(this_fire){
#     
#     # #Test
#     # this_fire<-list_f[[1]]
#     
#     #Pull in fire perimeter
#     perim<-vect(paste0(output_directory,"/fire_perimeters/merged/",this_fire,"_Extent.shp"))
#     perim<-project(perim,the_prj)
#     
#     #Consider adding 1km buffer?
#     
#     #Create blank raster
#     this_rast <- rast(extent = ext(perim), res = 5, crs = crs(the_prj))
#     
#     sf_cnt<-this_rast
#     sf_cnt$sf_cnt<-0
#     sf_cnt<-mask(sf_cnt,perim)
#     
#     
#     sf_max_yr<-this_rast
#     sf_max_yr$sf_max_yr<-0
#     sf_max_yr<-mask(sf_max_yr,perim)
#     sf_max_yr[sf_max_yr$sf_max_yr==0]<-NA
#     
#     
#     # write raster for individual fires 
#     writeRaster(sf_max_yr, paste0(output_directory,'/smfires/raster/',this_fire,"_sf_max_yr.tif"), overwrite = T)
#     writeRaster(sf_cnt, paste0(output_directory,'/smfires/raster/',this_fire,"_sf_cnt.tif"), overwrite = T)
#   }
#   
#   
#   ClippedFires<-lapply(1:length(list_f), function(i) {
#     ClipSmFires2Fire(list_f[i])
#     
#   })
#   
#   
# }
# smfire_download_rast(data_directory)
# 
# 
# 
# 
# 
# 
# 

