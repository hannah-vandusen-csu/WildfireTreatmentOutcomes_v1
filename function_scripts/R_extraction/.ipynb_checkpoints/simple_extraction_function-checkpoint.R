# R functions for data extraction at point or shape
# Hannah Van Dusen 2/10/2025
# Adapted to use bilinear interpolation - 11/20/2025
# updated by Mark Kreider March 2026


calculate_daily_area_spread <- function(fireID) { #calculate total fire spread of particular day
  
  dob_raw <- rast(paste0(local_directory,'/date_of_burn/dob_tif/', fireID, '_1_dob.tif'))
  growth <- expanse(dob_raw$DOB, unit = "ha", byValue = T) #calculates area of dob raster in hectares
  
  growth <- growth |> mutate(
    total_area = sum(area),
    cumulative = cumsum(area),
    percent_burned = area / total_area,
    relative_cumulation = area / cumulative) |> 
    select(!c(layer))
  return(growth)
}
simple_extraction_local <- function(local_directory,
                                     wto_outcome,
                                     column,
                                     buffer,
                                     buffer_size,
                                     point_directory,
                                     filename_suffix,
                                     gee_data,
                                     gee_suffix,
                                     ouput_fp,
                                     fire_list_fp,
                                     var_df){
  
  # fire list CSV
  fires <- read.csv(paste0(local_directory, '/fire_lists/', fire_list_fp, '.csv')) %>%
    tibble()
  fires$year <- lubridate::year( as.Date(fires$Ignition, format = "%m/%d/%Y"))
  if(any( is.na(fires$year))){
    fires$year <- lubridate::year(as.Date(
        fires$Ignition))  # this makes this flexible to excel-caused discrepancies in ignition date format
  }

  # Possible fires
  files <- list.files(point_directory)
  potential_fires <- unique(str_extract(files[!grepl("\\.txt$", files)], "^[^_]+"))
  fires <- fires[fires$MTBS_ID %in% potential_fires,]
  
  # variable cross walk sheet
  var_df <- read.csv(paste0(global_data_directory, '/variable_descriptions/variable_crosswalks_simple.csv'))
  var_outcome_df <- var_df[var_df[[column]] == TRUE,]
  var_r_df <-var_outcome_df %>%
    filter(extraction_source == 'r')
  
  # variable full cross walk sheet for final variable selection
  var_df_full <- read.csv(paste0(global_data_directory, '/variable_descriptions/full_variable_list.csv'))
  var_outcome_full <- var_df_full[var_df_full[[column]] == TRUE,]
  
  # read in roads
  if (var_df[var_df$Abbreviation == 'roads', ][column][1, ] == 'TRUE'){
    roads <- vect(paste0(local_directory, "/roads/roads_AOI.gpkg"))
  }
  
  if(var_df[var_df$Abbreviation == 'BPS',][column][1,] == 'TRUE'){
    BPS <- rast(paste0(dirname(local_directory), "/global_data/vegetation/bps/LF2016_BPS_200_CONUS/Tif/LC16_BPS_200.tif"))
    activeCat(BPS) <- "BPS_CODE"
    BPS_csv <- read.csv(paste0(dirname(local_directory), "/global_data/vegetation/bps/LF2016_BPS_200_CONUS/CSV_Data/LF16_BPS_200.csv")) |>
      select(c('BPS_CODE', 'BPS_NAME'))
  }
  
  
  # Extract data on a per fire basis
  f_list <- fires$MTBS_ID
  #fireID <- f_list[1]
  local_extraction_ls <- f_list %>% lapply(function(fireID){
    
    # Get fire year
    fire_year <- fires[fires$MTBS_ID == fireID,]$year[1]
    
    # read in sampling points
    if(length(list.files(paste0(point_directory), pattern=fireID)) > 0){
      
      files <- list.files(paste0(point_directory), pattern=fireID, full.names = TRUE)
      files <- files[str_detect(files, paste0(filename_suffix, "\\.csv$"))]
      sample_vect <- lapply(files, function(x){
        df <-read.csv(x)
        vect(df, geom=c("X", "Y"), crs=CRS)
      }) %>% bind_spat_rows()
      
      if(buffer){
        sample_vect <-buffer(sample_vect, buffer_size)
      }
      
      # 1) Extract treatment data from the treatment rasters
      if('Treatments' %in% var_r_df$Group & var_r_df[var_r_df$Group == 'Treatments', ][column][1, ] == 'TRUE'){
        
        # Extract treatment data
        treatments_fp <- paste0(local_directory, '/', var_r_df[var_r_df$Group == 'Treatments',]$file_path[1],
                                '/', fireID, '_', var_r_df[var_r_df$Group == 'Treatments',]$Abbreviation[1],'_', column,'.tif')
        
        if(wto_outcome == "containment"){
          sample_vect_buffer <- buffer(sample_vect, 25)  # buffer containment points by 25 meters
          treatments_df <- terra::extract(rast(treatments_fp), sample_vect_buffer, fun=function(x)
              Mode(x, na.rm = T)[1], ID = F, bind = F)
          
        }else{
          treatments_df <- terra::extract(rast(treatments_fp), sample_vect, method = "bilinear", ID=F, bind=F)
        }
        
        treatments_df <- treatments_df %>%
          as.data.frame() %>%
          tibble() %>%
          mutate_all(~ifelse( is.nan(.), NA,.)) %>%
          mutate(
            class_cod = ifelse(is.na(class_cod), 0, class_cod),  # all untreated areas are 0
            txYrs_cnt = ifelse( is.na(txYrs_cnt), 0, txYrs_cnt),  # 0 treatments in untreated areas
            tr_max_year = ifelse( is.na(tr_max_year), 1900, tr_max_year),  # 1900 for untreated areas
            uniqueID = ifelse( is.na(uniqueID), 999, uniqueID),  # 999 for unique ID of untreated areas
            tr_age = ifelse(class_cod != 0, fire_year - tr_max_year, 21) %>% as.integer(),
            treated = ifelse( class_cod == 0, 0, 1)) %>%
          bind_cols(values(sample_vect))
      } 
      
      # 2) Extract other simple extract rasters grouping by crs
      var_r_df_simple <- var_r_df %>%
        filter(Group != 'Treatments') %>%
        filter(Abbreviation != 'roads' & Abbreviation != 'frs' ) %>%
        mutate(pattern = paste0(file_path, '/.+_', Abbreviation, '\\.tif'))
      
      # raster tif file paths
      tif_files <- list.files(path = paste(sep = '/', local_directory, na.omit(unique(var_r_df_simple$file_path))),
                              recursive = T, full.names = T,
                              pattern = paste0('.*\\.tif$'))
      
      # Reading in Raster files for fire
      fire_tif_files <- str_subset(tif_files, fireID)
      fire_tif_files <- as.data.frame(unique(fire_tif_files)); colnames(fire_tif_files) <- 'fire_tif_files'
      
      rasters_df <- regex_inner_join(x = fire_tif_files, y = var_r_df_simple, by = c(fire_tif_files = "pattern"))  # INNER JOIN so some things are getting dropped, not sure if they actually should be or not
      rasters <- mapply(FUN = rast, unique(rasters_df$fire_tif_files), SIMPLIFY = F, USE.NAMES = F); names(rasters) <- unique(rasters_df$Abbreviation)
      rasters <- mapply(FUN = function(raster, layer_name){names(raster) <- layer_name;
        return (raster)}, rasters, layer_name = unique(rasters_df$Abbreviation))
      rasters <- rasters[order(names(rasters))]
      
      # grouping rasters by CRS
      crs_layers <- vapply(
        rasters,
        function(x) terra::crs(x, proj = TRUE),
        character(1)
      )
      
      unique_crs <- unique(crs_layers)
      
      crs_stacks <- setNames(
        vector("list", length(unique_crs)),
        unique_crs
      )
      
      
      for(j in 1:length(unique(crs_layers))){  # one loop for each CRS
        
        crs_stacks[[unique(crs_layers)[j]]] <- rasters[names(crs_layers)[which(crs_layers == unique(crs_layers)[j])]] %>% sprc()  # adding sublist of rasters matching the current CRS
        names(crs_stacks[[unique(crs_layers)[j]]]) <- names(crs_layers)[which(crs_layers == unique(crs_layers)[j])]
        
        
        crs_stacks[[unique(crs_layers)[j]]] <- terra::extract(crs_stacks[[unique(crs_layers)[j]]], sample_vect, method = "bilinear", na.rm = T, ID = F, bind = F) %>%
          bind_cols() %>%
          tibble() %>%
          bind_cols(values(sample_vect)) %>%
          mutate_all(~ifelse( is.nan(.), NA,.))
        
      }  # end CRS loop
      
      # join r data with treatment raster data
      if (!length(crs_stacks) == 0){
        preds <-tibble(!!unique_identifer := sample_vect %>% pull(!!unique_identifer))
        for (df in crs_stacks) {
          out <- df %>% 
            dplyr::select(!any_of(c(names(sample_vect)[! names(sample_vect) %in%unique_identifer])))
          preds = full_join(preds, out)
        }
      } else {
        preds <-tibble(!!unique_identifer := sample_vect %>% pull(!!unique_identifer))
      }
      
      if (exists("treatments_df")){
        full_df <- left_join(treatments_df, preds, by = unique_identifer) %>%
          mutate(sf_max_yr = ifelse( is.na(sf_max_yr), 0, sf_max_yr))
        
      } else {
        full_df <-preds
        full_df <- full_df %>% left_join( as.data.frame(sample_vect), by = unique_identifer)
      }
      
      # Calculate distance to roads
      if (var_df[var_df$Abbreviation == 'roads', ][column][1, ] == 'TRUE'){
        
        # Clip to fire perimeter
        f_extent <- vect(paste0(local_directory, '/fire_perimeters/merged/', fireID, '_Extent.shp'))
        f_extent <- project(f_extent, crs(sample_vect))
        
        roads_f <- crop(roads, terra::buffer(f_extent, 1000))
        
        if (dim(roads_f)[1] == 0){
          full_df$nearest_road <- 999999
        } else {
          dist <- nearest(sample_vect, roads_f)
          full_df$nearest_road <- dist$distance
        }
      }
      
      
      # Read in BPS
      if(var_df[var_df$Abbreviation == 'BPS',][column][1,] == 'TRUE'){
        
        BPS_df <-  as.data.frame(terra::extract(BPS, sample_vect, fun = median, na.rm = T, ID = F, bind = T)) |>
          select(!c(names(sample_vect)[! names(sample_vect) %in% unique_identifer ]))  |>
          mutate(BPS_CODE = as.numeric(as.character(BPS_CODE))) |>
          left_join(unique(BPS_csv), by = 'BPS_CODE') |> select(!c('BPS_CODE'))
        
        bps_davis <- read.csv(paste0(global_data_directory, '/vegetation/bps/orginal_davis_bps.csv')) |> 
          select(c('BpS_name', 'Cluster_name', 'Coarse_cluster'))
        
        BPS_df <- left_join(BPS_df, bps_davis, by = c('BPS_NAME' = 'BpS_name')) |>
          mutate(Cluster_name = ifelse(is.na(Cluster_name), BPS_NAME, Cluster_name),
                 Coarse_cluster = ifelse(is.na(Coarse_cluster), BPS_NAME, Coarse_cluster))
        
        full_df <- left_join(full_df, BPS_df, by = unique_identifer)
      }
      
      
      # Join in GEE data List all files
      gee_files <- list.files(paste0(local_directory, gee_data), pattern=fireID, full.names = TRUE)
      
      # Function to check if a file is empty
      is_empty_file <- function(file_path)
      {
        file_info <- file.info(file_path)
        return (file_info$size == 0)
      }
      empty_files <- sapply(gee_files, is_empty_file)
      gee_files <- gee_files[!empty_files]
      
      if (length(gee_files) == 0){
        
        print(paste0('No GEE sample points found for ', fireID))
        return (NULL)
        
      } else {
        gee_df <- lapply(gee_files, read.csv) %>% bind_rows() %>%
          dplyr::select(!c('.geo', 'system.index', 'ee_date', 'day', 'month')) %>%
          # mutate(SCF = ifelse(is.na(SCF), 0, SCF),
          #        SDD = ifelse(is.na(SDD), 0, SDD)) |>
          dplyr::select(-any_of(c(names(sample_vect)[! names(sample_vect) %in%unique_identifer], 'BPS')))
        
        if ("SCF" %in%names(gee_df)){
          gee_df <-gee_df %>%
            mutate(SCF = ifelse( is.na(SCF), 0, SCF))
        }
        
        if ("SDD" %in%names(gee_df)){
          gee_df <-gee_df %>%
            mutate(SDD = ifelse( is.na(SDD), 0, SDD))
        }
        
        full_df <- inner_join(full_df, gee_df, by=unique_identifer)
        
        # Calculate dob growth calculations
        dob_growth <- calculate_daily_area_spread(fireID)
        full_df <- left_join(full_df, dob_growth, by=c('min' = 'value')) %>%
          rename('DOB' = 'min')
        
        
        full_df$MTBS_ID <- fireID
        full_df$fire_year <- fire_year
        
        # Filter out any points that 'previously' burn at high severity
        if (any(colnames(full_df) %in% 'classified_past_sev')){
          full_df <- full_df %>% filter(classified_past_sev != 2)
        }
        
        # Subset to only variables in full variable list
        full_df_clean <- full_df %>% dplyr::select(MTBS_ID, response, response_bin, lat, long, X, Y,
                                                   any_of(na.omit(var_outcome_full$Abbreviation)))
        
        print(paste0('done with ', fireID))
        if (!dir.exists(paste0(local_directory, ouput_fp, '/', 'individual_fires'))){
          dir.create(paste0(local_directory, ouput_fp, '/', 'individual_fires'), recursive = TRUE)
        }
        fwrite(full_df_clean, paste0(local_directory, ouput_fp, '/individual_fires/', fireID, '_simple_extract.csv'))
        
        return (full_df_clean)
      }
      
    } else {
      print(paste0('No sample points found for ', fireID))
      return (NULL)
    }
    })
  

  
  local_extraction_ls <- local_extraction_ls[lengths(local_extraction_ls) != 0]
  all_fires <- rbindlist(local_extraction_ls, fill = T)
  if (!dir.exists(paste0(local_directory, ouput_fp))){
    dir.create(paste0(local_directory, ouput_fp), recursive=FALSE)
  }
  fwrite(all_fires, paste0(local_directory, ouput_fp, '/final_simple_extract.csv'))
  print('final simple extract done!')
}

