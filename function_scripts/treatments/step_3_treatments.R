####################################
###
### Step 3: Join FACTS and NFPORS Treatment Attributes to Shapefiles
### The script joins treatment csv attributes to the cleaned shapefile
### Exports a shapefile with fuel treatment attributes for each Fire Perimeter
###
### Written by Deborah Nemens May 2023 (dnemens@uw.edu)
### Edited by Hannah Van Dusen January 2024 (hannah.vandusen@usda.gov)
### Adapted for FACTS data by Betsy Black March 2024 (elizabeth.black@usda.gov)
###
####################################



# this script joins the csv of cleaned treatments produced by the step_2_treatments script 
# to the existing shapefile of treatment attributes for each fire produced by the step_1_treatments script 
# script to produce a shapefile of filtered polygons  
# Written by Deborah Nemens May 2023 (dnemens@uw.edu)
# Edited by Hannah Van Dusen January 2023 (hannah.vandusen@usda.gov)


sf::sf_use_s2(FALSE)#this switches off spherical geom, removes errors from facts layers not loading


treatment_spatialjoin <- function(
    data_directory,
    output_directory, 
    the_prj){

  # Import csv of fire names and ignition dates
  fires_raw <- read.csv(paste0(output_directory, '/fire_lists/selectedfires_raw.csv'))
  
  #Assign unique ID to every fire
  fires_raw$numID <- seq_len(nrow(fires_raw))
  options(scipen=10)
  set.seed(1234)
  fires_raw$numID <- (fires_raw$numID)*100000
  
  # Filter fires from step 2
  fires <- fires_raw |> dplyr::filter(step2_retention == 'TRUE')
  
  # Create list of fire names from csv, formatting can be deleted once we have edited this fire. list
  list_f <- fires$MTBS_ID
  
  # Read in classification CSV
  classification <- read.csv(paste0(data_directory, '/treatments/treatment_categories/class_synthesis_tcremoval.csv'), header=T)
  
  # Function to loop through each fire, join raw FACTS and NFPORS treatment shp to new "cleaned" treatment data csv, 
  # and output a "clean" shp with attributes
  Csv2ShapeJoin <- function(fire){
    
    # Select this fire
    this_fire = fire
    
    # List and Read csv of cleaned treatment attributes
    this_csv <- read.csv(paste0(output_directory, "/treatments/csv/", this_fire, "_treatments_step2.csv")) |>
      dplyr::select(UKCN, nfporsfid, year, trtcatego, trttype, sankey, class)
    
    
    # List and read shp of raw FACTS + NFPORS treatments
    this_shp_raw <- vect(paste0(output_directory, "/treatments/shapefile/", this_fire, "_treatments_raw.shp")) 
    
    this_shp <- tidyterra::inner_join(this_shp_raw,  this_csv, by = c('UKCN', 'nfporsfid'))
    
    # Disaggregate multipolgyons to single polygons
    this_shp <- terra::disagg(this_shp)
    
    # Set crs to NAD83 (reproject before calculating area)
    this_shp <- project(this_shp, the_prj)
    
    # Load fire perimeter and make sure crs matches
    extent <- vect(paste0(output_directory, "/fire_perimeters/merged/", this_fire, "_Extent.shp"))
    extent <- project(extent, the_prj)
    extent <- terra::buffer(extent, 1000) # add a 1000 m buffer for infrastructure
    
    # Re-clip disaggregated polygons to fire perimeter
    this_shp <- this_shp[extent, ]
    
    # Since we are cropping to 2000m from fire perimeter some fires might fall off
    if(dim(this_shp)[1] == 0){
      return(NULL)
    }else{
      rm(extent)
      
      # Get centroids and round
      
      shp_centroids <- this_shp |> centroids() |> crds() |> 
        as.data.frame() |> 
        dplyr::rename(cent_lat = x, cent_long = y) |>
        dplyr::mutate(cent_lat = round(cent_lat,digits = 0),
               cent_long = round(cent_long,digits = 0))
      
      
      # Bind centroids lats and longs to polygons in spatvector
      this_shp <- cbind(this_shp, shp_centroids)
      
      # Calculate polygon area
      this_shp$area_m <- this_shp |> terra::expanse(unit="m") |>
        round(digits = 0)
      
      # Calculate polygon perimeter
      this_shp$perim_m <- this_shp |> terra::perim() |>
        round(digits = 0)
      
      # Extract new version of CSV
      this_csv <- this_shp |> as.data.frame() |>
        dplyr::mutate(ACTIVITY = ifelse(ACTIVITY == "Piling of Fuels, Hand or Machine", "Piling of Fuels", ACTIVITY)) |>
        dplyr::select(any_of(c('ACTIVITY', 'trttype', 'class', 'year', 'cent_lat', 'cent_long','area_m', 'perim_m')))
      
      #  Keep record of treatment types performed on each polygon in each year
      types <- this_csv |>
        group_by(cent_lat, cent_long, area_m, perim_m, year) |>
        summarise(activitiesComb = toString(sort(unique(ACTIVITY))),
                  typesComb = toString(sort(unique(trttype))),
                  classesComb = toString(sort(unique(class))), .groups='drop') |>
        dplyr::mutate(activitiesComb = str_replace_all(activitiesComb, ",", " +"),
                      typesComb = str_replace_all(typesComb, ",", " +"),
                      classesComb = str_replace_all(classesComb, ",", " +")) |>
        dplyr::mutate(yr_activity = str_c(year, activitiesComb, sep = ": "), 
                      yr_type = str_c(year, typesComb, sep = ": "),
                      yr_class = str_c(year, classesComb, sep = ": "))
      
      
      # Combine year and activity
      types1 <- types |>
        group_by(cent_lat, cent_long, area_m, perim_m) |>
        summarise(activ_yr = toString(sort(unique(yr_activity))),
                  type_yr = toString(sort(unique(yr_type))),
                  class_yr = toString(sort(unique(yr_class))), .groups='drop')
      
      #And tidy list of all treatment types on each unique polygon
      types2 <- this_csv |>
        group_by(cent_lat, cent_long, area_m, perim_m) |>
        summarise(activ_all = toString(sort(unique(ACTIVITY))),
                  type_all = toString(sort(unique(trttype))),
                  class_all = toString(sort(unique(class))),
                  year_all = toString(sort(unique(year))), 
                  year_max = max(year, na.rm = T), .groups='drop') |>
        mutate(activ_all = str_replace_all(activ_all, ",", " +"),
               type_all = str_replace_all(type_all, ",", " +"),
               class_all = str_replace_all(class_all, ",", " +"),
               year_all = str_replace_all(year_all, ",", " +"))|>
        mutate(type_all = str_replace_all(type_all, ",", " +"))
      
      types <- full_join(types1, types2, by = c('cent_lat', 'cent_long', 'area_m', 'perim_m'))
      
      # Remove duplicates of unique polygons
      this_csv_clean <- this_csv |>
        distinct(cent_lat, cent_long, area_m, perim_m, .keep_all= TRUE)
      this_csv_clean <- left_join(this_csv_clean, types, by = c('cent_lat', 'cent_long', 'area_m', 'perim_m'))
      
      # For NFPORS rows where there are no activities, turn activ.all and activ.yr into NAs (instead of random meaningless text strings or blanks)
      fix <- function(column) {
        column[!grepl("[a-zA-Z]", column)] <- NA
        return(column)
      }
      
      this_csv_clean$activ_all <- fix(this_csv_clean$activ_all)
      this_csv_clean$activ_yr <- fix(this_csv_clean$activ_yr)
      
      # Merge with classification csv to assign "final class"
      this_csv_clean <- left_join(this_csv_clean, classification, by = 'class_all')
      
      # Assign unique ID to every unique polygon for each fire
      this_csv_clean <- this_csv_clean |>
        mutate(ID = seq_len(nrow(this_csv_clean)) + fires[fires$MTBS_ID==this_fire,]$numID[1]) |>
        dplyr::select(cent_lat:ID) #get rid of non-polygon specific columns
      
      # Name for output file
      this_csv_nm <- (paste0(output_directory, "/treatments/csv/", this_fire, "_treatments_cleaned.csv"))
      
      # write CSV
      write.csv(this_csv_clean, this_csv_nm, row.names=FALSE)
      
      # Specify the columns you want to keep
      columns_to_keep <- c("cent_lat","cent_long","area_m","perim_m")
      
      # Extract the subset of columns
      this_shp_clean <- this_shp[, columns_to_keep]
      
      # Join clean csv
      this_shp_clean <- tidyterra::inner_join(this_shp_clean, this_csv_clean, by= c("cent_lat","cent_long","area_m","perim_m"))
      
      # Remove duplicates
      this_shp_clean <- unique(this_shp_clean)
      
      # Fix invalid geometries
      # use sf package for spherical geometry issues
      if(sum(!is.valid(this_shp_clean)) > 0){
        this_shp_clean <- st_as_sf(this_shp_clean) |>
          sf::st_make_valid() |>
          sf::st_collection_extract(type = "POLYGON") |> # take only polygons from fixed geometries
          vect() 
      }
      
      this_shp_clean <- terra:: aggregate(this_shp_clean, by = c("cent_lat","cent_long","area_m","perim_m", "ID"), fun = 'modal', dissolve  =T)
      
      # Aggregate takes the "mode" but in all cases rows are identical 
      names(this_shp_clean) <- gsub("^modal_", "", names(this_shp_clean))
      
      # Fix invalid geometries
      if(sum(!is.valid(this_shp_clean)) > 0){
        this_shp_clean <- st_as_sf(this_shp_clean) |>
          sf::st_make_valid() |>
          sf::st_collection_extract(type = "POLYGON") |> # take only polygons from fixed geometries
          vect() 
      }
      
      # Name for output file
      this_shp_nm <- paste0(output_directory, "/treatments/shapefile/", this_fire, "_treatments_cleaned.shp")
      
      # Remove file if it exits -- weird overwrite 
      if(file.exists(this_shp_nm)){
        file.remove(this_shp_nm)
      }
      
      # Output shapefile of clean treatments for each fire
      writeVector(this_shp_clean, this_shp_nm,insert = TRUE, overwrite = TRUE)
      
      print(paste0('Joining Treatments ', this_fire))
      return(this_fire)
    }
    
  }
  
  fires_found <- lapply(list_f, Csv2ShapeJoin)
  
  # Return a New Selected Fire Shape file with fires with Successfully Downloaded Treatments
  CleanFires <- unlist(fires_found)
  
  # Subset Selected Fires to those Successfully Downloaded
  fires_raw$step3_retention <- fires_raw$MTBS_ID %in% CleanFires
  
  # Save Selected Fires to raw CSV
  write.csv(fires_raw, paste0(output_directory, '/fire_lists/selectedfires_raw.csv'),row.names=FALSE)
  
  # Save Final selected fires as CSV
  # Subset Selected Fires to those Successfully Downloaded
  fires_final <- fires_raw |> filter(step3_retention == TRUE)
  write.csv(fires_final, paste0(output_directory, '/fire_lists/selectedfires.csv'),row.names=FALSE)
  
}

