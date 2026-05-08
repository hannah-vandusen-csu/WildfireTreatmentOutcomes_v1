####################################
###
###
### Fire Selection
### Hannah Van Dusen --- 12/01/2023
### hannah.vandusen@usda.gov
###
###
####################################


fire_selection <- function(data_directory, output_directory , CRS) {

  # Receive User Input
  if(!exists("userInput")) userInput <- getUserInput()
  
  ####
  ####
  # 00. Select fires with treatments and Save Selected Fires as a csv and a shape file -----------------------------
  ###
  ###

  # Load in MTBS shapefile
  mtbs <- file.path(paste0(data_directory, "/data/global_data/fire_perimeters/mtbs/mtbs_perims_DD.shp")) |>
    vect() |>
    select(Event_ID, irwinID, Incid_Name, Incid_Type, Ig_Date) |>
    mutate(Ig_Date = as.Date(Ig_Date)) |>
    mutate(Year = year(Ig_Date)) |>
    # Narrow down shapefile to user-selected year range
    filter(between(Year, userInput[['Year_range']][1], userInput[['Year_range']][2])) |>
    project(CRS)
  
  # Remove Prescribed fires in MTBS data set
  mtbs <- filter(mtbs,
                 (is.na(Incid_Type) | Incid_Type != "Prescribed Fire"),
                 (is.na(Incid_Name) | str_detect(Incid_Name, "RX", negate = T)))
  
  

  # Cropped Fire Perimeters to AOI: depending on users input
  Strategy_landscapes <- data_directory |>
    paste0('/data/global_data/strategy_landscapes/AllPriorityLandscapes.gdb') |>
    vect()

  if(sum(userInput[['AOI']] %in% Strategy_landscapes$NAME) > 1) { # if the user selected one or more strategy landscapes

    Strategy_landscapes <- Strategy_landscapes |>
      filter(NAME %in% userInput[['AOI']]) |>
      project(CRS)

    mtbs <- mtbs[Strategy_landscapes, ] # Crop to selected strategy landscape

  } else { # assumes the user selected a state that's in the dataset, and that they did NOT also select a strategy landscape

    mtbs$State_ID <- sapply(mtbs$Event_ID, function(x) { substr(x, 1, 2) } )
    mtbs <- mtbs[mtbs$State_ID %in% userInput[['AOI']], ]

  }

  # Read in NIFC perimeters for more accurate perimeters 
  # NIFC Perimeters currently there is only perimeters for years 2018 - 2023
  nifc_all <- vect(paste0(data_directory, r'(/data/global_data/fire_perimeters/nifc/perims_2017_2023.shp)')) |>
    project(CRS)
  
  # Irwin ID's for MTBS perimeters
  irwin_mtbs <- na.omit(mtbs$irwinID)
  nifc_all <- dplyr::filter(nifc_all,
                     poly_IRWIN %in% paste0("{", irwin_mtbs, "}")) |>
    terra::aggregate(by = 'poly_IRWIN', dissolve = TRUE)
  
  # Combine perimeters datasets 
  selectedfires <- rbind(mtbs, nifc_all) |>
    mutate(irwinID = if_else(is.na(irwinID), str_replace_all(poly_IRWIN, "[{}]", ""), irwinID)) |> 
    select(c('Event_ID', 'irwinID', 'Incid_Name', 'Ig_Date'))
  
  # separate fires with an irwin ID and without and irwinID 
  selectedfires_w_irwin <- selectedfires[!is.na(selectedfires$irwinID), ]
  selectedfires_wo_irwin <- selectedfires[is.na(selectedfires$irwinID), ]
  # dissolve perimeters with same irwin ID & add MTBS information back in 
  selectedfires_w_irwin <- selectedfires_w_irwin |> select('irwinID') |> terra::aggregate(by = 'irwinID', count = FALSE) |> 
    left_join(as.data.frame(selectedfires_w_irwin[!duplicated(selectedfires_w_irwin$irwinID),]), by= 'irwinID')
 
  # Selected fires with correct headers
  selectedfires_merged <- rbind(selectedfires_w_irwin, selectedfires_wo_irwin) |>
    mutate(irwinID = paste0("{", irwinID, "}")) |> 
    select(c('Event_ID', 'Incid_Name', 'Ig_Date', 'irwinID')) |>
    rename(c('MTBS_ID' = 'Event_ID', 'Name' = 'Incid_Name', 'Ignition' = 'Ig_Date'), 'poly_IRWINID' = 'irwinID')
    
  
  # Extract selectedfires attribute table and rename columns with Pilot Project Syntax
  SelectedFiresCSV <- as.data.frame(selectedfires_merged)

  # Save Selected Fire List as a single CSV for downstream use
  dir.create(paste0(output_directory, '/fire_lists'), recursive = T)
  
  SelectedFiresCSV <- SelectedFiresCSV[order(SelectedFiresCSV$Ignition),]
  write.csv(SelectedFiresCSV,
            paste0(output_directory, '/fire_lists/selectedfires_raw.csv'),
            row.names = F)

  # Create 3 folders to save fire perimeters
  # 1) mtbs perimeter
  # 2) nifc perimeter
  # 3) merged perimeter 
  dir.create(paste0(output_directory, '/fire_perimeters/mtbs'), recursive = T)
  dir.create(paste0(output_directory, '/fire_perimeters/nifc'), recursive = T)
  dir.create(paste0(output_directory, '/fire_perimeters/merged'), recursive = T)

  save_fire_perimeters <- function(x){
    # mtbs 
    mtbs |> filter(Event_ID ==  x[['MTBS_ID']]) |> 
      writeVector(paste0(output_directory, '/fire_perimeters/mtbs/', x[['MTBS_ID']], '_Extent.shp'),
                  overwrite = T)
    # if a nifc perimeter exists save a nifc perimeter
    if(x[['poly_IRWINID']] %in% nifc_all$poly_IRWIN){
      nifc_all |> filter(poly_IRWIN ==  x[['poly_IRWINID']]) |> 
        writeVector(paste0(output_directory, '/fire_perimeters/nifc/', x[['MTBS_ID']], '_Extent.shp'),
                    overwrite = T)
    }
    
    # merged
    selectedfires_merged |> filter(MTBS_ID ==  x[['MTBS_ID']]) |> 
      writeVector(paste0(output_directory, '/fire_perimeters/merged/', x[['MTBS_ID']], '_Extent.shp'),
                  overwrite = T)
  }
  
  apply(SelectedFiresCSV, 1, FUN = save_fire_perimeters , simplify = F)

  print(paste0('Number of Potential Fires Selected before Treatment Cleaning: ', nrow(SelectedFiresCSV)))

}

