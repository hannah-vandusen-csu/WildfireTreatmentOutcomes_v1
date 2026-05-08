##########################################################################
# this function clusters fires based on similarity of climate variables to
# ensure sufficient post-balancing sample sizes for modeling
##########################################################################

# load packages -----

# Necessary packages
pkgs <- c("ggplot2",
          "tidyverse",
          "rstudioapi",
          "data.table",
          "purrr",
          "tidyterra",
          "terra",
          "raster"
)




# Install any that are missing
missing <- pkgs[!(pkgs %in% installed.packages()[,"Package"])]
if (length(missing) > 0) {install.packages(missing)}

# Load all packages
lapply(pkgs, library, character.only = TRUE)

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


# outcome - should be one of burn_severity, containment, or infrastructure
project <- "burn_severity" 

# Define response type 
response_type <- "binary" 

# response
# Set the response to the correct column name dynamically based on response_type 
#"response" for containment and infrastructure, "response_bin" or "response" for burn severity
# Burn severity can be a binary or continuous response
response <- ifelse(
  project == "burn_severity",
  ifelse(response_type == "binary", "response_bin", "response"), # Choose based on response type
  ifelse(project %in% c("containment", "infrastructure"), "response", NA) # Default to "response"
)

# csv output from point-level extraction
df <- fread(paste0(processed_data_directory, '/processed_outputs/', project, 
                   '/extracted_dataframes/point_level_extracts/final_simple_extract.csv'))


# treatment column
treat_col <- 'treated'

# output folder - you shouldn't need to change this
outfold <- paste0(processed_data_directory, '/processed_outputs')

# function to cluster data --------

df %>% as.data.frame() %>%
  group_by(MTBS_ID) %>%
  reframe(n = n())

#---------------FUNCTION-------------------#
cluster_data <- function(df, id_col= "MTBS_ID", response, mets = c("CMD", "Eref"), min0=75, min1=75, trt_per_fire = TRUE) {
  
  # # Make sure there is at least one treatments in each fire
  # if(trt_per_fire){
  #   df <- df %>%
  #     group_by(!!sym(id_col)) %>%
  #     filter(any(get(treat_col) != 0)) 
  # }
  
  # compute evapotranspiration
  df$et <- df$Eref - df$CMD
  
  mets[mets=="Eref"] <- "et"
  
  df <- df[complete.cases(dplyr::select(df, all_of(mets))), ]#remove rows with NAs (all columns)
  
  fire_data <- df %>%
    dplyr::select(all_of(c(id_col, mets))) %>%
    group_by(!!sym(id_col)) %>%
    summarise_at(vars(contains(mets)), ~median(as.numeric(.x))) %>%
    ungroup()
  names(fire_data) <- c(id_col, mets)
  
  cntdf <- df %>% group_by(!!sym(id_col)) %>% 
    summarise(cnt0 = sum(!!sym(response) == 0),
              cnt1 = sum(!!sym(response) == 1),
              tot = n())
  names(cntdf) <- c(id_col, c("cnt0", "cnt1", "tot"))
  
  fire_data <- left_join(fire_data, cntdf, by = id_col)
  
  fire_data$cluster_id <- 1:nrow(fire_data)
  
  fire_data <- fire_data %>% mutate_at(vars(contains(mets)), ~scale(.x))
  
  check_criteria <- function(id, clust_df, m0=min0, m1=min1) {
    cluster_data <- filter(clust_df, cluster_id == id)
    count0_tot <- sum(cluster_data$cnt0) < m0
    count1_tot <- sum(cluster_data$cnt1) < m1
    return(any(count0_tot, count1_tot))
  }
  
  ids <- unique(fire_data$cluster_id)
  clusters_to_merge <- sapply(ids, check_criteria, clust_df = fire_data)
  clusters_to_merge <- ids[clusters_to_merge]
  
  while(length(clusters_to_merge) > 1) {
    cluster_data_to_merge <- filter(fire_data, cluster_id %in% clusters_to_merge)
    
    cluster_centroids <- cluster_data_to_merge %>%
      dplyr::select(all_of(c("cluster_id", mets))) %>%
      group_by(get("cluster_id")) %>%
      summarise_at(vars(contains(mets)), ~mean(as.numeric(.x))) %>%
      ungroup()
    names(cluster_centroids) <- c("cluster_id", mets)
    
    dmat <- as.matrix(dist(dplyr::select(cluster_centroids, all_of(mets))))
    diag(dmat) <- Inf
    
    closest_clusters <- which.min(dmat)
    cluster_pair <- arrayInd(closest_clusters, dim(dmat))
    
    cluster_id_1 <- cluster_centroids$cluster_id[cluster_pair[1]]
    cluster_id_2 <- cluster_centroids$cluster_id[cluster_pair[2]]
    
    fire_data[fire_data$cluster_id == cluster_id_2, "cluster_id"] <- cluster_id_1
    
    ids <- unique(fire_data$cluster_id)
    clusters_to_merge <- sapply(ids, check_criteria, clust_df = fire_data)
    clusters_to_merge <- ids[clusters_to_merge]
    
    print(length(clusters_to_merge))
  }
  
  # merge final cluster with closest if criteria not met
  ids <- unique(fire_data$cluster_id)
  final_clust_to_merge <- sapply(ids, check_criteria, clust_df = fire_data)
  final_clust_to_merge <- ids[final_clust_to_merge]
  
  if(length(final_clust_to_merge) > 0) {
    cluster_centroids <- fire_data %>%
      dplyr::select(all_of(c("cluster_id", mets))) %>%
      group_by(get("cluster_id")) %>%
      summarise_at(vars(contains(mets)), ~mean(as.numeric(.x))) %>%
      ungroup()
    names(cluster_centroids) <- c("cluster_id", mets)
    
    cluster_centroids <- rbind(filter(cluster_centroids, cluster_id == final_clust_to_merge),
                               filter(cluster_centroids, cluster_id != final_clust_to_merge))
    dmat <- as.matrix(dist(dplyr::select(cluster_centroids, all_of(mets))))
    diag(dmat) <- Inf
    nw_cluster_id <- cluster_centroids$cluster_id[which.min(as.vector(dmat[,1]))]
    
    fire_data[fire_data$cluster_id == final_clust_to_merge, "cluster_id"] <- nw_cluster_id
  }
  
  fire_data <- dplyr::select(fire_data, all_of(c(id_col, "cluster_id")))
  
  new_df <- df
  new_df <- left_join(new_df, fire_data, by = id_col)
  
  outpath <- file.path(outfold, project, "clustering_outputs/all_samples_w_clusterid.csv")
  
  if(!dir.exists(file.path(outfold, project, "clustering_outputs"))){
    dir.create(file.path(outfold, project, "clustering_outputs"))
  }
  write.csv(new_df, outpath, row.names = F)
  
  # plotting
  
  plt_df <- new_df %>% group_by(MTBS_ID) %>% summarise(cnt=n(),
                                                       CMD = mean(CMD),
                                                       et= mean(et),
                                                       cluster_id = first(cluster_id))
  
  plt <- ggplot() +
    geom_point(data = plt_df, aes(x = CMD, y = et, color = as.character(cluster_id), size = cnt)) +
    theme_bw() +
    labs(x = "climatic moisture deficit", y = "actual evapotranspiration") +
    guides(color=guide_legend(ncol=4, title = "Cluster ID"),
           size=guide_legend(title = "Total Sample Size"))
  
  ggsave(paste0(outfold, "\\", project, "\\clustering_outputs\\clustering_plot.png"), plt, width = 10, height = 8)
  
}

#---------------APPLY-------------------#

cluster_data(df=df, id_col="MTBS_ID",response=response,trt_per_fire = TRUE)

