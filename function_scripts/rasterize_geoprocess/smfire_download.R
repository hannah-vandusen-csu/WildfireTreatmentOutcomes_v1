# function to download small fires and clip to fire perimeters
smfire_download <- function(
    
  data_directory,
  output_directory,
  the.prj
  
){
  
  # Import CSV of fire names and Ignition Dates
  fires <- read.csv(paste0(output_directory,"/fire_lists/selectedfires.csv"))
  
  # Create list of fire names from CSV
  list_f <- fires$MTBS_ID
  
  # Create a list of ignition dates from CSV
  ign_f <- as.Date(fires$Ignition)
  
  #Load in other wildfires
  smfires_shp <- file.path(paste0(data_directory,"/fire_perimeters/usgs_small_fires/InterAgencyFirePerimeterHistory_All_Years_View.shp"))
  smfires <- vect(smfires_shp) %>%
    project(the.prj)
  
  # Function to clip FACTS + NFPORS treatment layers to fire perimeters, 
  # output shapefile and csv of clipped treatment data
  ClipSmFires2Fire <- function(this_fire, this_ign){
    
    #Get fire year
    this_year = as.numeric(year(this_ign))
    
    # Load fire extent
    extent <- vect(paste0(output_directory, "/fire_perimeters/merged/", this_fire, "_Extent.shp")) %>%
      project(the.prj)
    
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
        filter(FEATURECA == "Prescribed Fire"|FEATURECA == "Wildfire Final Perimeter"|FEATURECA == "Wildfire") %>% 
        filter(FIREYEAR < this_year)%>% # must occur before fire
        filter(FIREYEAR >= lower_bound)
      
      if (nrow(smfires.clipped) >0){
        
        #Filter out duplicate fires
        smfires.clipped <- smfires.clipped%>%
          mutate(INCIDENT =str_replace_all(tolower(INCIDENT),"\\s+","")) %>%
          unique() %>%
          mutate(FEATURECA=factor(FEATURECA,levels=c("Wildfire Final Perimeter","Wildfire","Prescribed Fire")))%>%
          arrange(INCIDENT, FIREYEAR, FEATURECA)%>%
          group_by(INCIDENT, FIREYEAR)%>%
          slice(1)%>%
          ungroup() %>%
          select(INCIDENT,FIREYEAR,FEATURECA,SOURCE,AGENCY,UNITID,IRWINID,COMMENTS,GISACRES,ShapeAre,ShapeLen)
        
        #Extract CSV of small fires
        this_csv <- as.data.frame(smfires.clipped)
        
        #Name and file path
        this_fp <- paste0(output_directory,'/smfires')
        
        # Create necessary directories
        if(!dir.exists(this_fp)){
          dir.create(this_fp)}
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