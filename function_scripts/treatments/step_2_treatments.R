####################################
###
### Step 2: FACTS + NFPORS Treatment Cleaning
### The script cleans raw treatment data for each fire deleting treatment replicates
### Exports a shapefile and csv of fuel treatment attributes for each Fire Perimeter
###
### Written by Deborah Nemens May 2023 (dnemens@uw.edu)
### Edited by Hannah Van Dusen January 2024 (hannah.vandusen@usda.gov)
### Adapted for raw FACTS and NFPORS data by Betsy Black March 2024 (elizabeth.black@usda.gov)
###
###
####################################



# this script automates the process of cleaning FACTS and NFPORS treatment data for fires of interest
# Treatment_raw files for each fire were created in R using the "step_1_treatments.R" script
# Written by Deborah Nemens May 2023 dnemens@uw.edu
# Edited by Hannah Van Dusen January 2023 hannah.vandusen@usda.gov

# library(stringr)
# library(dplyr)
# library(sf)
# library(terra)
# library(lubridate)
# # Set working Directory to Current R project location 
# 
# 
# # Get directory where scripts are saved
# current_directory <- dirname(rstudioapi::getActiveDocumentContext()$path)
# 
# # Navigate to parent and project directory 
# parent_directory <- dirname(dirname(dirname(current_directory)))
# data_directory<-parent_directory
# setwd(parent_directory)

# Function Creation

treatment_cleaning <- function(
    data_directory,
    output_directory){
  
  # Import CSV of names and Ignition Dates
  fires <- read.csv(paste0(output_directory,"/fire_lists/selectedfires_raw.csv")) |> 
    filter(step1_retention == 'TRUE') |>
    mutate(Ignition = as.Date(Ignition))
  
  # Read in csv of activity types that = fuel treatments
  facts_cats <- read.csv(paste0(data_directory, '/data/global_data/treatments/treatment_categories/facts.classes.csv'), header=T)
  nfpors_cats <- read.csv(paste0(data_directory, '/data/global_data/treatments/treatment_categories/nfpors.classes.csv'), header=T)
  
  # Functions to clean FACTS and NFPORS df
  clean_facts <- function(facts_df){
    # Only FACTS
    this_csv_clean <- full_join(facts_df, facts_cats, by = join_by(ACTIVITY)) |>
      filter(!(is.na(SUID))) |>  # Merge cats with this this_csv fire activities
      distinct(year, SUID, trttype, ACTIVITY, .keep_all= TRUE) |> # keep only unique treatments
      filter(!((is.na(class))|class=="other")) # Filter out activities that were not fuel treatments
    
    
    
    # Remove facts treatments that were actually wildfires or "Preparation" and wildfire/preparation
    this_csv_clean <- this_csv_clean |>
      filter(!ACTIVITY %in% c("Wildfire - Natural Ignition",
                              "Wildfire - Human Ignition", 
                              "Wildfire - Fuels Benefit",
                              "Planned Treatment Burned in Wildfire",
                              "Wildland Fire Use")) |> 
      filter(!trtcatego %in% c("Preparation for Treatment",
                               "Wildfire Non-Treatment"))
    return(this_csv_clean)
  }
  clean_nfpors <- function(nfpors_df){
    # Only NFPORS
    this_csv_clean <- full_join(nfpors_df, nfpors_cats) |>
      filter(!(is.na(GISACRES))) |> #remove rows with no tx
      distinct(year, GISACRES, centlat, centlon, catnm, trttype, .keep_all= TRUE) |> 
      filter(!((is.na(class))|class=="other")) |> # Filter out activities that were not fuel treatments
      filter(!trtcatego %in% c("Preparation for Treatment","Wildfire Non-Treatment"))
    return(this_csv_clean)
  }
  
  cleaned_csvs <- function(x){
    
    this_fire <- x[['MTBS_ID']]
    this_ign <-  x[['Ignition']] 
    
    tryCatch({
      
      # Read a csv of treatment attributes
      this_csv_all <- read.csv(paste0(output_directory, "/treatments/csv/", this_fire, '_treat_attributes_raw.csv'))
      
      this_csv_all <- this_csv_all |>
        filter(!is.na(DATECOMPLETED)) |>                   # Remove treatments that are not yet completed
        mutate(DATECOMPLETED = as.Date(DATECOMPLETED)) |>  
        filter(DATECOMPLETED < this_ign) |>                # Filter out treatments that occurred after the ignition date
        mutate(year = year(DATECOMPLETED)) |>              
        filter(year > (year(this_ign) - 21))                # Only include treatments that have occured in the last 20 years
      
      
      # Split up facts and nfpors data for processing
      facts <- filter(this_csv_all,!is.na(SUID))
      nfpors <- filter(this_csv_all,is.na(SUID))
      
      
      # Facts and nifpors have different attributes
      if(dim(nfpors)[1] == 0){
        
        # Only FACTS
        clean_df <- clean_facts(facts)
        
      }else{ 
        if(dim(facts)[1]==0){
          
          # Only NFPORS
          clean_df <- clean_nfpors(nfpors)
          
        }else{
          # Both FACTS and NFPORS
          clean_facts_df <- clean_facts(facts)
          clean_nfpors_df <- clean_nfpors(nfpors)
          clean_df <- bind_rows(clean_facts_df, clean_nfpors_df)
        }
      } 
      
      # Output cleaned list of treatments as CSV
      
      this_csv_path <- paste0(output_directory, "/treatments/csv/")
      
      #name for output file
      this_csv_nm <- paste0(this_fire, "_treatments_step2.csv")
      
      #output csv of treatment attributes to individual treatment folder
      write.csv(clean_df, paste0(this_csv_path, this_csv_nm),row.names=FALSE)
      
      # If no clean treatments are found delete folder for data storage
      if(dim(clean_df)[1] == 0){
        
        # Add to list 
        print(paste0('No cleaned treatements found in ',this_fire,' fire.'))
        return(NULL)
        
        
      }else{
        
        
        print(paste0('Cleaning Treatments of ',this_fire, ' complete.'))
        # Add to list 
        return(this_fire)
      }
      
      
    }, error = function(e){
      
      # Error Message
      print(paste0("Error occured for ", this_fire, ":", e$message))
      
      # Add to list 
      return(NULL)
    }
    
    ) 
  }
  CleanTreatmentsFound <- apply(fires, 1, cleaned_csvs)
 
  # Return a New Selected Fire Shape file with fires with Successfully Downloaded Treatments
  CleanFires <- unlist(CleanTreatmentsFound)
  
  # Subset Selected Fires to those Successfully Downloaded
  fires$step2_retention <- fires$MTBS_ID %in% CleanFires
  
  # Save Selected Fires as CSV
  write.csv(fires, paste0(output_directory, '/fire_lists/selectedfires_raw.csv'),row.names=FALSE)
  
  
}



