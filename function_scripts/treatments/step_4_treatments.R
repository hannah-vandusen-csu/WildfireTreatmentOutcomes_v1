####################################
###
### Step 4: Treatment Rasterization
### This step rasterizes the classified treatment polygons 
### Combines overlapping treatment polygons for final raster class 
### Outputs:
### Raster of treatments per fire
### with layers: Final Classs Factorize (final_class_factor), 
### Most recent treatment year (max_tr_yr), 
### Unique ID or Combined Unique ID for each treament or overlapping treatment
###
### Written by Hannah Van Dusen April 2024 (hannah.vandusen@usda.gov)
###
####################################


treatment_to_raster <- function(
    data_directory,
    output_directory,
    the_prj = CRS
){

  # Create Folder for Treatment Rasters 
  if(!dir.exists(paste0(output_directory, '/treatments/rasters'))){
    dir.create(paste0(output_directory, '/treatments/rasters'))
  }
  
  # Import csv of fire names and ignition dates
  fires_raw <- read.csv(paste0(output_directory, '/fire_lists/selectedfires_raw.csv'))
  
  # Filter fires from step 3
  fires <- fires_raw |> dplyr::filter(step3_retention == 'TRUE')
  
  # Create list of fire names from csv, formatting can be deleted once we have edited this fire. list
  list_f <- fires$MTBS_ID

  # Set which treatment ID to group by
  ID <- 'ID'

# Function to rasterize treatments
rasterize_treaments <- function(fire){
  
  # Fire name, ID, and year
  this_fire <- fire
  this_fire_ID <- fires[fires$MTBS_ID == this_fire,]$numID
  fire_year <- lubridate::year(fires[fires$MTBS_ID == this_fire,]$Ignition)
  
  # Load Treatment Shapes
  this_shp <- vect(paste0(output_directory, "/treatments/shapefile/", this_fire, "_treatments_cleaned.shp"))
  if(file.exists(paste0(output_directory, '/treatments/rasters/', this_fire, "_treatments_classed.tif"))){
    print(paste0('Raster already exists for ', this_fire))
    return(paste0('Raster already exists for ', this_fire))
  }else{
    # Project shape
    this_shpVect <- project(this_shp, the_prj)
    
    # Create blank raster at 5m resolution over the shapefile extent
    this_rast <- rast(extent = ext(this_shpVect), res = 5, crs = the_prj)
    
    # ---- Step 1: Count overlaps per cell ----
    count_rast <- rasterize(this_shpVect, this_rast, field = 1, fun = "sum")
    
    has_overlap <- global(count_rast > 1, "max", na.rm = TRUE)$max > 0
    
    if (has_overlap) {
      
      # ---- Step 2: Rasterize ALL cells using terra::extract on a point grid ----
      # But ONLY for cells that are touched by at least one polygon
      # This avoids the enormous xyFromCell(1:ncell) approach
      
      # Get cells touched by any polygon (non-NA in count raster)
      valid_cells <- cells(count_rast)
      valid_xy <- xyFromCell(this_rast, valid_cells)
      
      # Extract polygon attributes at each valid cell
      # relate = "intersects" ensures all overlapping polygons are returned per point
      valid_pts <- vect(valid_xy, crs = the_prj)
      
      # terra::extract returns one row PER polygon-point intersection
      # So a pixel overlapped by 11 shapes -> 11 rows with the same point id
      extracted <- terra::extract(this_shpVect, valid_pts)
      
      # Attach coordinates via the point ID
      extracted$x <- valid_xy[extracted$id.y, 1]
      extracted$y <- valid_xy[extracted$id.y, 2]
      
      df <- as.data.table(extracted)
      df[[ID]] <- as.character(df[[ID]])
      rm(extracted, valid_pts, valid_xy)
      
      # ---- Step 3: Separate overlap vs non-overlap pixels ----
      pixel_counts <- df[, .N, by = .(x, y)]
      
      # Non-overlapping pixels (single polygon)
      single_xy <- pixel_counts[N == 1, .(x, y)]
      df_single <- df[single_xy, on = .(x, y)]
      
      # Overlapping pixels (multiple polygons)
      multi_xy <- pixel_counts[N > 1, .(x, y)]
      df_multi <- df[multi_xy, on = .(x, y)]
      
      # ---- Step 4: Combine overlapping polygon IDs per pixel ----
      ID_col = 'ID'
      df_combined <- df_multi[, .(
        Combined = paste(sort(unique(.SD[[ID_col]])), collapse = "+")
      ), by = .(x, y), .SDcols = ID_col]
      
      # Non-overlapping: just the single ID
      df_single_combined <- df_single[, .(
        Combined = as.character(.SD[[ID_col]])
      ), by = .(x, y), .SDcols = ID_col]
      
      rm(df, df_multi,df_single, single_xy, pixel_counts, multi_xy)
      # # ---- Step 5: Existing overlap resolution logic ----
      # overlapped_pixels <- df_combined
      
      if (nrow(df_combined) > 0) {
        max_splits <- max(str_count(df_combined$Combined, "\\+")) + 1
        # max_splits could be up to 11+
        
        col_names <- paste0(ID, 1:max_splits)
        
        
        overlapped_pixels <- copy(df_combined)
        overlapped_pixels[, (col_names) := tstrsplit(Combined, "\\+", fixed = FALSE, fill = NA)]
        
        
        nonoverlapped_pixels <- df_single_combined
        colnames(nonoverlapped_pixels)[3] <- paste0(ID, "1")
        
        ID_pixels <- bind_rows(nonoverlapped_pixels, overlapped_pixels)
        rm(df_combined, df_single_combined, nonoverlapped_pixels, overlapped_pixels)
        
        # Read in cleaned and categorized CSV for treatments 
        this_csv <- read.csv(paste0(output_directory, "/treatments/csv/", this_fire, "_treatments_cleaned.csv"))
        this_csv <- this_csv |>
          dplyr::select(year_max, class_all, ID, year_all) |>
          mutate(ID = as.character(ID), 
                 year_all = as.character(year_all))
      
        # Join with pixel values 
        for (n in 1:length(col_names)){
          csv_clean_n <- this_csv
          colnames(csv_clean_n) <- paste0(colnames(csv_clean_n), n)
          ID_pixels <- left_join(ID_pixels, csv_clean_n, by = col_names[n])
        }
    
        # Get all Final Class Columns 
        Final_Class_columns <- grep("^class_all", names(ID_pixels), value=T)
        
        # Combine final Class
        ID_pixels$FinalClassComb <- apply(
          ID_pixels[, ..Final_Class_columns], 1,
          function(x) paste0(sort(unique(unlist(strsplit(x, '\\s*\\+\\s*')))), collapse = ' + ')
        )   
       
        # Get Max Year 
        year_cols <- grep("^year_max", names(ID_pixels), value = TRUE)
        ID_pixels[, tr_max_year := do.call(pmax, c(.SD, na.rm = TRUE)), .SDcols = year_cols]
        
        # Combine All Years 
        TxYrs_columns <- grep("^year_all", names(ID_pixels), value=T)
        
        # Set NA's as ''
        for (col in TxYrs_columns) {
          set(ID_pixels, which(is.na(ID_pixels[[col]])), col, "")
        }
        
        # Count number of unique treatment years 
        ID_pixels$txYrs_cnt <- apply(
          ID_pixels[, ..TxYrs_columns], 1,
          function(x) length(unique(unlist(strsplit(x, '\\s*\\+\\s*'))))
        )
      
      # Add Pixel ID
       ID_pixels[, ID_combined := fifelse(is.na(Combined), ID1, Combined)]
      
      # Make unique id for strings of pixel IDs
      shpID <- unique(ID_pixels$ID_combined)
      uniqueID <- dim(this_shpVect)[1] + 1:length(shpID) + this_fire_ID
      ids_df <- as.data.frame(cbind(uniqueID,shpID))
      ID_pixels <- merge(ID_pixels, ids_df, by.x = "ID_combined", by.y = "shpID", all.x = T)
      
      # Read in classification csv
      classes <- read.csv(paste0(data_directory, '/treatments/treatment_categories/class_synthesis_tcremoval.csv'), header=T)
      
      # Assign true 'final class': class.cod
      ID_pixels <- merge(ID_pixels, classes, by.x = "FinalClassComb", by.y = "class_all", all.x = T)
      
      # Select desired columns
      ID_pixels_xyz <- ID_pixels[,c('x','y','class_cod', 'tr_max_year', 'uniqueID','txYrs_cnt')]
      
      # Make sure uniqueID is an integer
      ID_pixels_xyz[, uniqueID := as.integer(uniqueID)]
      
      # Define the output raster template from extent/res BEFORE loading all data
      # (requires you know the CRS, resolution, and extent ahead of time)
      r_template <- rast(
        xmin = min(ID_pixels_xyz$x),
        xmax = max(ID_pixels_xyz$x),
        ymin = min(ID_pixels_xyz$y),
        ymax = max(ID_pixels_xyz$y),
        resolution = 5,          
        crs = the_prj
      )
      
      # Write directly to disk — never fully materializes in RAM
      setwd(paste0(output_directory, '/treatments/rasters'))
      trt_comb <- rasterize(
        x    = ID_pixels_xyz[, c("x", "y")],
        y    = r_template,
        values = ID_pixels_xyz[, c("uniqueID", "class_cod", "tr_max_year", 'txYrs_cnt')],
        fun  = "last",
        filename = paste0(this_fire, '_treatments_classed.tif'),
        overwrite = TRUE,
        wopt = list(gdal = c("COMPRESS=LZW", "TILED=YES"))  # smaller file
      )
      
      # Save Raster Key as csv
      unique_treatments <- unique(ID_pixels[,-c(1,2)])
      write.csv(unique_treatments, paste0(this_fire, '_raster_key.csv'),row.names=FALSE)
      
      # Print complete statement
      print(paste0('Overlaps found. Raster complete for fire ', this_fire))
      }else{
    print(paste0('Overlaps found, but not valid overlaps. Check fire ', this_fire))
  }
    }else{
   
      # Read in cleaned and categorized CSV for treatments 
      this_csv <- read.csv(paste0(output_directory, "/treatments/csv/", this_fire, "_treatments_cleaned.csv"))
      this_csv$ID <- as.character(this_csv$ID)
  
      this_csv <- this_csv %>%
        mutate(txYrs_cnt = str_count(year_all, "\\+") + 1)
      
      counts<-this_csv%>%
        dplyr::select(ID,txYrs_cnt)
  
      # Merge csv columns to shp to add counts in
      this_shpVect <- merge(this_shpVect, counts, by = 'ID')
      

      # Rasterize ID, Final Class, and max year from treatment shape
      this_shpVect$tr_max_year <- values(this_shpVect[,'year_max'])
      this_shpVect$ID_combined <- values(this_shpVect[,ID])
      this_shpVect$uniqueID <- values(this_shpVect[,ID])
      
      rast_class <- rasterize(this_shpVect, this_rast, field = 'class_cod', fun = max)
      rast_year <- rasterize(this_shpVect, this_rast, field = 'tr_max_year', fun = max)
      rast_ID <- rasterize(this_shpVect, this_rast, field = 'uniqueID', fun = max)
      rast_cnt <- rasterize(this_shpVect, this_rast, field = 'txYrs_cnt', fun = max)
      
      trt_comb <- c(rast_class, rast_year, rast_ID, rast_cnt)
      crs(trt_comb) <- the_prj
      trt_comb <- crop(trt_comb, this_shpVect)
      
  
    # Write Raster
    setwd(paste0(output_directory, '/treatments/rasters'))
    writeRaster(trt_comb, paste0(this_fire, '_treatments_classed.tif'), overwrite = T)
    
    # Match raster key format
    this_csv = this_csv %>%
      rename(FinalClassComb = class_all,
             tr_max_year = year_max,
             ID_combined = ID)
    this_csv$uniqueID = this_csv$ID_combined
    
    # Select columns from csv
    this_csv <- this_csv%>%
      dplyr::select(year_all,FinalClassComb,tr_max_year,txYrs_cnt,ID_combined,uniqueID,class_fin,class_cod)
    
    # Save Raster Key as csv
    write.csv(this_csv, paste0(this_fire, '_raster_key.csv'),row.names=FALSE)
    
    # Print complete statement
    print(paste0('No overlaps found. Raster complete for fire ', this_fire))

  }
  
}
}
  

  # lapply
  lapply(list_f, rasterize_treaments)

}

