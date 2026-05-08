# Script to create landscape-level dataframes

# This script is up to this user desire to split the current selected fires 
# and respective paired treated and untreated points in respective landscapes. 
# This gives the user multiple dataframes to run landscape-level models 

# In the Klamath study, we divided the selected fires into the "east" and "west"
# by spliting the data by the east and west eco-regions. This was informed by discussion 
# with local managers

# For data exploration, we also split data into treated and untreated
# and days with extreme fire growth vs normal fire grouwth 


# This script  double checks that there is enough sample points in each 
# treatment class for each reach to retain in the modeling process. 


# load packages -----

# List of required packages
packages <- c("sf", "terra", "tidyverse", "ggplot2", "rdist", "leaflet")

# Install any packages that aren't already installed
installed <- packages %in% rownames(installed.packages())
if (any(!installed)) {
  install.packages(packages[!installed])
}

# Load all packages
lapply(packages, library, character.only = TRUE)

# read in user input and set directories -------

# set working directory based on this script (needs to be within the broader folder)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

#read in user input
user_input_csv <- read.csv("user_input/user_input_general.csv", stringsAsFactors = FALSE)

# loop through user input csv and assign the values provided to the named objects
for (i in seq_len(nrow(user_input_csv))){
  val <- user_input_csv$value[i]
  if(user_input_csv$treat_as_R_expression[i] == T){val <- eval(parse(text = val))}
  assign(user_input_csv$name[i], val, envir = .GlobalEnv)
  rm(val)
}


# Set dataframe and directories ------------------------------------------------
# Load in dataframe
matched_filepath <- file.path(processed_data_directory, 
                              "processed_outputs", wto_outcome, 'matching_outputs/csv')

# could be the same as input
output_filepath <- matched_filepath

# Read in dataframe
df <- read.csv(file.path(matched_filepath, 'final_matched_df.csv'))


# Split dataframes into east and west landscapes -------------------------------
# Now eperate east and west extracted dataframes
# Create spatial df
df_spatial <- st_as_sf(df, coords = c("X", "Y"), crs = 26910)

# Ecoregions LIII
ecos <- st_read(paste0(global_data_directory, '/vegetation/ecoregions/us_eco_l3/us_eco_l3.shp'))
ecos <- st_transform(ecos, st_crs(df_spatial))
ecos <- dplyr::select(ecos, "US_L3CODE")

# Intersect
df_spatial <- st_intersection(df_spatial, ecos)
#mapview(df_spatial, zcol = "epa_3")

df_spatial$epa_3 <- df_spatial$US_L3CODE
df_spatial <- dplyr::select(df_spatial, !US_L3CODE)

df_spatial <- cbind(df_spatial, st_coordinates(df_spatial))
st_geometry(df_spatial) <- NULL
head(df_spatial)

df_region <- df_spatial |> mutate(region = ifelse(epa_3 == 4 | epa_3 == 9, 'East', 'West'))


# Treatment classes to drop ----------------------------------------------------
# Allows users to visually look at each treatment class 
# For burn severity we dropped a treatment class if there was < 20 unique 
# treatments for a given treatment class

# Number of treatments for a give class to retain
num_tr_threshold <- 20 

tr_summary <- df_region |> filter(class_cod != 0) |>
  filter(!duplicated(uniqueID)) |> # unique treatments
  group_by(class_cod, region) |>
  summarise(total_unique_tr = n())
print(tr_summary)


drop_treatments <- unique(tr_summary[tr_summary$total_unique_tr <= num_tr_threshold,]$class_cod)


matched_ids_tr_drops <- df_region |> filter(class_cod %in% drop_treatments) |>
  pull(matched_id)

df_region <- df_region |> filter(!matched_id %in% matched_ids_tr_drops)

# Filter out treatments that are 100x larger than the median treatment value for each class in each region
# These are most likely invalid / inaccurage treatment shapes
tr_df <- df_region |> filter(class_cod != 0)
tr_unique <- tr_df[!duplicated(tr_df$uniqueID),]

X <- 100
tr_extreme_rows_med <- tr_unique %>%
  mutate(tr_area = tr_area / 10000) %>%
  group_by(class_cod, region) %>%
  mutate(
    med = median(tr_area, na.rm = TRUE)) %>%
  dplyr::select(class_cod, region, uniqueID, tr_area, med)

tr_extreme_rows_med <- tr_extreme_rows_med[!is.na(tr_extreme_rows_med$med) &
                                             tr_extreme_rows_med$med > 0 &
                                             tr_extreme_rows_med$tr_area / tr_extreme_rows_med$med > X,]

matched_ids_size_drops <- df_region |> filter(uniqueID %in% tr_extreme_rows_med$uniqueID) |>
  pull(matched_id)

df_region <- df_region |> filter(!matched_id %in% matched_ids_tr_drops)

# Split dataframe by region ----------------------------------------------------

df_east <- df_region |> filter(region == 'East')
df_west <- df_region |> filter(region == 'West')

# write csv -- clean full, east, and west
write.csv(df_region, file.path(output_filepath, "all_pnts_region.csv"))
write.csv(df_east, file.path(output_filepath, "all_pnts_east.csv"))
write.csv(df_west, file.path(output_filepath, "all_pnts_west.csv"))


df_east <- read.csv(file.path(output_filepath, "all_pnts_east.csv"))
df_west <- read.csv(file.path(output_filepath, "all_pnts_west.csv"))



# # Treatment only dataframes 
# df_east_tr <- df_east |> filter(class.cod != 0)
# df_west_tr <- df_west |> filter(class.cod != 0)
# write.csv(df_east_tr, file.path(working_directory, "interim_outputs/burn_severity/extracted_dataframes/point_and_upwind_extracts/all_pnts_east_only_treatments.csv"))
# write.csv(df_west_tr, file.path(working_directory, "interim_outputs/burn_severity/extracted_dataframes/point_and_upwind_extracts/all_pnts_west_only_treatments.csv"))
# 
# 
# # Extreme growth dataframes 
# df_east_exgrowth <- df_east |> filter(area > 4900)
# df_west_exgrowth <- df_west |> filter(area > 4900)
# write.csv(df_east_exgrowth, file.path(working_directory, "interim_outputs/burn_severity/extracted_dataframes/point_and_upwind_extracts/all_pnts_east_ex_growth.csv"))
# write.csv(df_west_exgrowth, file.path(working_directory, "interim_outputs/burn_severity/extracted_dataframes/point_and_upwind_extracts/all_pnts_west_ex_growth.csv"))
# 
# 
# # Treatment only dataframes 
# df_east_smgrowth <- df_east |> filter(area <= 4900)
# df_west_smgrowth <- df_west |> filter(area <= 4900)
# write.csv(df_east_smgrowth, file.path(working_directory, "interim_outputs/burn_severity/extracted_dataframes/point_and_upwind_extracts/all_pnts_east_sm_growth.csv"))
# write.csv(df_west_smgrowth, file.path(working_directory, "interim_outputs/burn_severity/extracted_dataframes/point_and_upwind_extracts/all_pnts_west_sm_growth.csv"))