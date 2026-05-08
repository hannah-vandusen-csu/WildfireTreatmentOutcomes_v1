# function to get roads data 
get_roads_data <- function(global_data_directory,
                           strategy_landscapes_of_interest,
                           output_directory,
                           the.prj){
  # load state boundary polygon
  states <- states()
  
  #Select Strategy Landscapes
  Strategy_landscapes <- vect(paste0(global_data_directory,'/strategy_landscapes/AllPriorityLandscapes.gdb'))
  
  # Select names to include the shape of the landscape of interest
  partial_AOI <- Strategy_landscapes[Strategy_landscapes$NAME %in% strategy_landscapes_of_interest,] %>%
    project(the.prj)
  
  # get the state in question
  state_name_vec <- partial_AOI$STATE %>% unique()
  
  # make a union of AOI polygon
  partial_AOI_union <- partial_AOI %>% terra::aggregate()
  
  # load counties within the state(s)
  counties_in_state <- state_name_vec %>% lapply(function(x){
    counties_state <- counties(state = x) %>%
      vect() %>%
      project(the.prj)
    return(counties_state)
  }) %>% bind_spat_rows()
  
  # find all counties that intersect with the area of interest
  counties_AOI <- counties_in_state %>% 
    terra::intersect(partial_AOI) %>%
    as.data.frame() %>%
    dplyr::select(STATEFP, COUNTYFP) %>%
    distinct()
  
  # download roads for each county and clip to AOI
  roads_AOI <- 1:nrow(counties_AOI) %>% lapply(function(x){
    print(paste("processing county", x, "of", nrow(counties_AOI)))
    roads_subset <- roads(counties_AOI$STATEFP[x], counties_AOI$COUNTYFP[x])
    roads_clip <- roads_subset %>% st_intersection(partial_AOI_union %>% 
                                                     st_as_sf() %>%
                                                     st_transform(crs(roads_subset)))
    roads_clip_vect <- roads_clip %>% vect() %>% project(the.prj)
    gc()
    return(roads_clip_vect)
  }) %>% bind_spat_rows()
  
  # save roads layer
  dir <- paste0(output_directory, "/roads/")
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  writeVector(roads_AOI, file.path(output_directory, "roads", "roads_AOI.gpkg"))
}
