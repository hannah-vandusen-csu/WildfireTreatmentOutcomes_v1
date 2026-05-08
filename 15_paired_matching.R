#####################################################################.
##  For each treated point, find a matched (paired) control point  ##.
#####################################################################.

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

# End of User Input ------------------------------------------------------------

# set local directory
local_directory <- file.path(processed_data_directory, "processed_outputs")

# Final extracted data to match with
df <- read.csv(paste0(local_directory, "/", wto_outcome, "/clustering_outputs/all_samples_w_clusterid.csv")) %>%
  mutate(point_id_old = point_id) %>%
  mutate(point_id = 1:nrow(.))

# Weights dataframes
var_dfs <- list.files(paste0(local_directory, "/", wto_outcome, "/matching_model_outputs/csv/weights/"), full.names=T, pattern = ".csv$")

# Output folder
out_fold <- paste0(local_directory, "/", wto_outcome, "/matching_outputs/")

# Exact matching variables, besides MTBS_ID
exact_match <- c("Coarse_cluster")

# Masks for controls
trt_fils <- list.files(file.path(processed_data_directory, "masks/treatment"),
                       full.names = T, pattern = "treatment.tif$")

pfire_fils <- list.files(file.path(processed_data_directory, "masks/wildfire"),
                         full.names = T, pattern = "wildfire.tif$")


# ------------------FUNCTIONS---------------------#

findMatch <- function(trt_row, conts, this_df, exact_match, other_vars) {
  
  #print(trt_row$point_id)
  
  tmat <- as.matrix(dplyr::select(trt_row, !other_vars))
  
  these_conts <- filter(conts, !!sym(exact_match) == trt_row[,exact_match])
  if (nrow(these_conts)>0) {
    cmat <- as.matrix(dplyr::select(these_conts, !other_vars))
    
    these_dists <- as.vector(cdist(tmat, cmat))
    
    selected_control <- these_conts[which.min(these_dists),]
    
    mtbs <- substr(trt_row$MTBS_ID, 3, 11)
    this_trt_id <- paste0(mtbs, "_", trt_row$cluster_id, "_", trt_row$point_id)
      
    nw_trt_row <- this_df %>% 
      filter(MTBS_ID == trt_row$MTBS_ID & cluster_id == trt_row$cluster_id & point_id == trt_row$point_id) %>%
      mutate(matched_id = this_trt_id)
    
    nw_cont_row <- this_df %>% 
      filter(MTBS_ID == selected_control$MTBS_ID & cluster_id == selected_control$cluster_id & point_id == selected_control$point_id) %>%
      mutate(matched_id = this_trt_id)
    
    nw_rows <- rbind(nw_trt_row, nw_cont_row)
    
    return(nw_rows)
  }
}

pairedMatching <- function(id, df, var_dfs, exact_match="Coarse_cluster", wto_outcome=wto_outcome) {
  print(id)
  if (wto_outcome=="infrastructure") {
    this_df <- dplyr::filter(df, cluster_id == id)
  } else if (wto_outcome=="containment" | wto_outcome=="burn_severity") {
    this_df <- dplyr::filter(df, MTBS_ID == id)
  }
  
  # masking possible controls
  if (wto_outcome == "burn_severity" | wto_outcome == "containment") {
    # convert to points and pull out controls
    these_pnts <- st_as_sf(this_df, coords=c("X","Y"), crs=CRS)
    
    cpnts <- filter(these_pnts, treated == 0)
    
    if(dim(cpnts)[1] == 0){
      print(paste0('No control points found for ', id))
      return(NULL)
    }
    
    # read in masks
    trt_mask <- rast(trt_fils[grepl(id, trt_fils)])
    fire_mask <- rast(pfire_fils[grepl(id, pfire_fils)])
    
    # extract masks to points
    cpnts <- st_transform(cpnts, crs(trt_mask))
    cpnts$trt_mask <- terra::extract(trt_mask, cpnts)[,2]
    cpnts$trt_mask <- ifelse(cpnts$trt_mask == 1, 1, NA)
    
    cpnts <- st_transform(cpnts, crs(fire_mask))
    cpnts$fire_mask <- terra::extract(fire_mask, cpnts)[,2]
    cpnts$fire_mask <- ifelse(cpnts$fire_mask == 1, 1, NA)
    
    # filter out masked points
    cpnts <- filter(cpnts, is.na(trt_mask) & is.na(fire_mask))
    
    # combine back with original dataframe
    cpnts <- st_transform(cpnts, CRS)
    cpnts <- cbind(cpnts, st_coordinates(cpnts))
    st_geometry(cpnts) <- NULL
    cpnts <- dplyr::select(cpnts, colnames(this_df))
    this_df <- rbind(filter(this_df, treated==1), cpnts)
  }

  # get cluster id
  this_clust_id <- unique(this_df$cluster_id)
  
  # define variables
  other_vars <- c("treated", "point_id","MTBS_ID", "cluster_id",  exact_match)
  
  # get var df associated with this fire's cluster
  var_df_nms <- as.numeric(str_extract(basename(var_dfs), "^[^_]+"))
  this_var_df <- read.csv(var_dfs[which(var_df_nms == this_clust_id)])
  this_var_df <- arrange(this_var_df, desc(proportional_weights))
  
  # select only matching variables
  this_df_matching <- dplyr::select(this_df, all_of(c(other_vars, this_var_df$Variable)))
  
  # Filter out variables that don't have any variation
  no_variation <- this_df_matching |> dplyr::select(!all_of(other_vars))|> dplyr::select(where(~ n_distinct(.) == 1)) |> colnames()
  
  # rearrange weights df by order of columns
  these_vars <- colnames(dplyr::select(this_df_matching, -any_of(c(other_vars, no_variation))))
  these_wts <- this_var_df[match(these_vars, this_var_df$Variable),]$proportional_weights
  
  # scale matching vars
  scaled_df <- this_df_matching %>% dplyr::select(-any_of(c(other_vars,no_variation))) %>%
    mutate(across(everything(), .fns = scale))
  # apply weights to matching vars
  scaled_df <- sweep(scaled_df, 2, these_wts, FUN="*")
  scaled_df <- cbind(dplyr::select(this_df_matching, other_vars), scaled_df)
  
  trts <- scaled_df %>% filter(treated == 1)
  conts <- scaled_df %>% filter(treated == 0)
  
  if(dim(trts)[1] == 0){
    print(paste0('No treated points found for ', id))
    return(NULL)
  }
  if(dim(conts)[1] == 0){
    print(paste0('No untreated points found for ', id))
    return(NULL)
  }
  
  trts <- split(trts, 1:nrow(trts))
  
  # find matches
  results <- lapply(trts, findMatch,
                    conts=conts,
                    this_df=this_df,
                    exact_match=exact_match,
                    other_vars=other_vars)
  results <- do.call(rbind, results)
  if(is.null(results)){
    return(NULL)
  }
  # export comparison figure
  plt_df_allcs <- this_df %>% filter(treated==0) %>% dplyr::select(all_of(these_vars)) %>% mutate("Type" = "All Controls")
  plt_df_trts <- this_df %>% filter(treated==1) %>% dplyr::select(all_of(these_vars)) %>% mutate("Type" = "Treated")
  plt_df_conts <- results %>% filter(treated==0) %>% dplyr::select(all_of(these_vars)) %>% mutate("Type" = "Matched Controls")
  
  plt_df <- rbind(plt_df_allcs, plt_df_trts, plt_df_conts)
  plt_df$Type <- factor(plt_df$Type, levels = c("All Controls", "Treated", "Matched Controls"))
  
  plt_df <- plt_df %>% pivot_longer(!Type)
  
  vars_ordered <- this_var_df[match(these_vars, this_var_df$Variable),]$Variable
  plt_df$name <- factor(plt_df$name, levels = vars_ordered)
  
  
  plt <- ggplot(plt_df) +
    geom_violin(aes(x = Type, y = value, fill = Type)) +
    facet_wrap(~name, scales = "free_y", nrow = 3) +
    theme_bw() +
    scale_fill_manual(values=alpha(c("green","orange","blue"))) +
    theme(legend.position = "bottom")

    ggsave(paste0(out_fold, "/figures/", id, "_comparison_plots.png"),
           plt, width = 6, height = 8)
  return(results)
  
}


#---------------APPLY FUNCTIONS---------------------#



# Create a treated column if necessary 
if(!'treated' %in% colnames(df)){
  df <- df |> mutate(treated = ifelse(class_cod == 0 , 0, 1))
}

# Filter out extremely large treatments x10 the next maximum treatment
# These are most likely fake treatments 
# before matching 
df_untreated <- df |> filter(treated == 0)
df_treated <- df |> filter(treated == 1)

tr_unique <- df_treated[!duplicated(df_treated$uniqueID),]
tr_summary <- tr_unique %>%
  mutate(tr_area = tr_area/10000) |> # convert to ha for easier viewing 
  group_by(class_cod) %>%
  summarise(
    tr_area_max = max(tr_area, na.rm = TRUE),
    tr_area_max_no_max = ifelse(
      n() > 1,
      ifelse(
        sum(tr_area == max(tr_area, na.rm = TRUE), na.rm = TRUE) > 1,
        max(tr_area, na.rm = TRUE), # If max appears more than once, removing one still leaves the same max
        max(tr_area[tr_area != max(tr_area, na.rm = TRUE)], na.rm = TRUE)
      ),
      NA_real_
    ), 
    max_unqiueID = ifelse(
      n() > 1,
      unique(uniqueID[tr_area == max(tr_area, na.rm = TRUE)]),
      NA_real_
    )
  )

# find treatments that go over the 10x greater then next largest treatment 
tr_extremes <- tr_summary |> filter(tr_area_max / tr_area_max_no_max > 10 ) 


# Filter out the extreme treatments and their matched points
tr_extremes <- tr_summary |> filter(tr_area_max / tr_area_max_no_max > 10 ) |> pull(max_unqiueID)
if(length(tr_extremes) > 0){
  df <- df |> filter(!uniqueID %in% tr_extremes)
}

if (wto_outcome=="infrastructure") {
  ids <- unique(df$cluster_id)
} else if (wto_outcome=="containment" | wto_outcome=="burn_severity") {
  ids <- unique(df$MTBS_ID)
}

matched_dfs <- lapply(ids, pairedMatching,
                      df=df,
                      var_dfs=var_dfs,
                      exact_match="Coarse_cluster",
                      wto_outcome=wto_outcome)
matched_dfs <- matched_dfs[lengths(matched_dfs)!=0 ] 
final_df <- do.call(rbind, matched_dfs)

if (!dir.exists(paste0(out_fold, "/csv"))) {
  dir.create(paste0(out_fold, "/csv"), recursive = TRUE)
}
write.csv(final_df, paste0(out_fold, "/csv/final_matched_df.csv"))

