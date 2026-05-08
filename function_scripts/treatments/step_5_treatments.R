####################################
###
### Filtering Treatment Rasters For Each out come
###
####################################


treatments_by_outcome <- function(data_directory, output_directory, script_directory, the_prj, wto_outcome) {
  
  # Import CSV of names and Ignition Dates
  fires <- file.path(output_directory, "fire_lists/selectedfires.csv") %>%
    read.csv()
  fires$year <- year(fires$Ignition)
  
  # Create list of fire names from CSV
  list_f <- fires$MTBS_ID
  
  # Read in outcome specifications
  outcome_params <- file.path(script_directory, 'user_input/treatment_selection.csv') %>%
    read.csv()
  
  # Function to clip FACTS + NFPORS treatment layers to fire perimeters,
  treatments4outcomes <- function(fireID, wto_outcome) {
  
  # Read in Treatment Raster
  treatments <- rast(paste0(output_directory, '/treatments/rasters/', fireID, '_treatments_classed.tif'))
 
  # Read in Fire perimeter
  extent <- vect(paste0(output_directory, '/fire_perimeters/merged/', fireID , '_Extent.shp'))
  extent <- project(extent, the_prj)
  
  # 1. Clip treatments set distance for outcome
  buffer_distance <- outcome_params |> filter(parameter == 'perim_buffer') |> dplyr::select(all_of(wto_outcome))  |> pull()
  buffered_extent <- terra::buffer(extent, width = buffer_distance)
  
  if (relate(treatments, buffered_extent, "intersects")[,1]) {
    treatments_cropped <- mask(crop(treatments, buffered_extent), buffered_extent)
    rm(treatments)
  } else {
    print(paste0('No treatments found for ', fireID))
    return(NULL)
  }
  
  # treatments_cropped <- mask(terra::crop(treatments, buffered_extent), buffered_extent)
  
  # 2. Filter treatments to maximum and minimum age of treatments
  if(global(treatments_cropped$uniqueID, fun = "notNA")$notNA > 0){
    
    treatments_cropped_age <- treatments_cropped
    rm(treatments_cropped)
    
    fire_year <- fires[fires$MTBS_ID == fireID,]$year
    min_age <- outcome_params |> filter(parameter == 'min_treatment_age') |> dplyr::select(c(wto_outcome))  |> pull()
    max_age <- outcome_params |> filter(parameter == 'max_treatment_age') |> dplyr::select(c(wto_outcome))  |> pull()
    treatments_cropped_age[treatments_cropped_age$tr_max_year > (fire_year - min_age)] <- NA
    treatments_cropped_age[treatments_cropped_age$tr_max_year <= (fire_year - max_age)] <- NA
    
    # 3. Filter treatments by treatment area
    if(global(treatments_cropped_age$uniqueID, fun="notNA")$notNA > 0){
      
      # Filter for minimum size 
      min_size_m <- outcome_params |> filter(parameter == 'min_treatment_size ') |> dplyr::select(c(wto_outcome))  |> pull()
      
      ## --- existing steps (comb + patches) ---
      treatments_size <- treatments_cropped_age
      rm(treatments_cropped_age)
      # comb <- treatments_size$class_cod * 10000 + treatments_size$tr_max_year
      # #comb <- as.factor(comb) # make sure layer is categorical
      # 
      # patch_area <- landscapemetrics::lsm_p_area(comb) 
      # # 1) Build combined categorical raster
      comb <- treatments_size$class_cod * 10000 + treatments_size$tr_max_year
      
      # # 2) Patch raster per class, saved to disk (memory-friendly)
      # if(!dir.exists(paste0(output_directory, '/treatments/rasters/temp'))){
      #   dir.create(paste0(output_directory, '/treatments/rasters/temp'))
      # }
      # 
      # setwd(paste0(output_directory, '/treatments/rasters/temp'))
      # options(to_disk = TRUE)
      # patch_list <- get_patches(comb, directions = 8, to_disk = TRUE, return_raster = TRUE)[[1]][[1]]
      # 
      # # 3) Patch areas (ha) by patch id
      pa <- lsm_p_area(comb, directions = 8)
      # pa has columns: class, id, metric, value (ha)
      
      # 4) Keep only patches >= minimum size threshold
      keep_classes <- pa$class[pa$value >= (min_size_m/10000)]
      
      # 5) Build mask from patch IDs
      mask_ok <- comb
      if(length(keep_classes) > 0){
        mask_ok[!mask_ok %in% keep_classes] <- NA
        
        # Now check if all treatments have been masked out
        if(global(mask_ok, "notNA")$notNA > 0){
          mask_ok[!is.na(mask_ok)] <- 1
        }else{
          print(paste0('No treatments found for ', fireID))
          return(NULL)
        }
      }
      
      treatments_size <- terra::mask(treatments_size, mask_ok)
      
      # units hectares
      # comb <- treatments_size$class_cod * 10000 + treatments_size$tr_max_year
      # comb <- lapp(treatments_size[[c("class_cod","tr_max_year")]],
      #              fun = function(x, y) as.integer(x) * 10000L + as.integer(y))
      # comb <- as.integer(comb)          # keep as numeric/int
      # levels(comb) <- NULL              # drop categories if any
      # 
      # p <- terra::patches(comb, directions = 8)
      # p
      # p <- terra::patches(comb, directions = 8)
      # 
      # # patch areas
      # patch_area <- as.data.frame(freq(p))[,c(2,3)]
      # names(patch_area) <- c("patch_id", "area_cells")
      # 
      # cell_area <- prod(res(p)) # For large extents (larger than a single fire perimeter) maybe use cellSize to incorporate differences in pixel areas
      # patch_area$area_m2 <- patch_area$area_cells * cell_area
      
      # --- filter patches based on threshold ---
      # small_ids <- patch_area$class[patch_area$value < (min_size_m/10000)]
      # 
      # # create mask: NA for small patches, 1 for others
      # mask_ok <- comb
      # 
      # if(length(small_ids) > 0){
      #   mask_ok[mask_ok %in% small_ids] <- NA
      #   
      #   # Now check if all treatments have been masked out
      #   if(global(mask_ok, "notNA")$notNA > 0){
      #     mask_ok[is.na(mask_ok)] <- 1
      #   }else{
      #     print(paste0('No treatments found for ', fireID))
      #     return(NULL)
      #   }
      # }
      # 
      # # apply to all layers of original raster
      # treatments_size <- terra::mask(treatments_size, mask_ok)

            # 4. Write Raster
      if(global(treatments_size$uniqueID, fun="notNA")$notNA > 0){
        print(paste0('Writing rasters for ', fireID))
        writeRaster(treatments_size, paste0(output_directory, '/treatments/rasters/', fireID, '_treatments_', wto_outcome,'.tif'), overwrite = TRUE)
        return(fireID)
       
      }else{
        print(paste0('No treatments found for ', fireID))
        return(NULL)
      }
    }else{
      print(paste0('No treatments found for ', fireID))
      return(NULL) 
    }
  }else{
    print(paste0('No treatments found for ', fireID))
    return(NULL) 
  }
  
}


  # Apply Function to Fires
  fires_complete <- lapply(list_f, treatments4outcomes, wto_outcome)
  
  
  # Return a New Selected Fire Shape file with fires with Successfully Downloaded Treatments
  CleanFires <- unlist(fires_complete)

  # Subset Selected Fires to those Successfully Downloaded
  fires[[wto_outcome]] <- fires$MTBS_ID %in% CleanFires 
    
  fires <-  dplyr::select(fires, -contains("_retention"))
  
  # Save Selected Fires to raw CSV
  write.csv(fires, paste0(output_directory, '/fire_lists/', wto_outcome, '_fire_list.csv'), row.names=FALSE)


}


