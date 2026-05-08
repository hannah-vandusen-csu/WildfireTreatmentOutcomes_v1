##
##
## Function that conduction variable reduction 
##
##

#' Process a single group of data
#'
#' @param group_id The value of the group identifier (or NULL if no grouping)
#' @param group_col The name of the grouping column (or NULL if no grouping)
#' @param df_raw The raw data frame
#' @param response Response variable name
#' @param var_filtered Variable matching information
#' @param resp_type Response type ('continuous' or 'binary')
#' @param retained_vars Variables to force back into final model
#' @param retained_groups Top variable of each group will be retained
#' @param variable_crosswalk Path to variable crosswalk file
#' @param wto_outcome Outcome type
#' @param spatial Boolean: spatial layers added
#' @param rfe Boolean: RFE is conducted
#' @param hyper_tune Boolean: tuing hyperparameters of random forest
#' @param sp_subsampling Boolean: spatial subsampling is done
#' @param sampling_dist numeric: sampling distance
#'
#' @return Analysis results for the group
variable_reduction <- function(group_id, group_col, df_raw, var_filtered, resp_type, response, 
                               retained_vars, retained_groups, variable_crosswalk, wto_outcome, 
                               spatial, rfe, hyper_tune, sp_subsampling, sampling_dist) {
  set.seed(1234)
  
  # Filter data for this group if grouping is used
  if (!is.null(group_id) && !is.null(group_col)) {
    df <- df_raw |> dplyr::filter(!!sym(group_col) == group_id)
  } else {
    df <- df_raw
    group_id = ''
  }

  # Balance data if a binary variable and create response as factor
  if(resp_type == 'binary'){
    df <- balance_binary(df, analysis_step, wto_outcome) |> as.data.frame()
  }
  
  # Prep variable names 
  preds <- var_filtered$Abbreviation
  cvars <- var_filtered[var_filtered$Continuous == 'yes',]$Abbreviation
  fvars <- var_filtered[var_filtered$Continuous == 'no',]$Abbreviation
  
  # add EVs if spatial
  if (spatial) {
    eigs <- create_eigen(df, num=30)
    df <- cbind(df, eigs)
    cvars <- c(cvars, colnames(eigs))
  } 
  
  # Run functions 
  morans_i = 1
  while(morans_i > 0.1){
    # Step 1) run correlation 
    prepped_data <- prepare_data_for_model(data = df, 
                                           numeric_cols = cvars, 
                                           factor_cols = fvars, 
                                           exclude_cor_vars = if(spatial){colnames(eigs)}else{c()}, 
                                           run_cor = TRUE, 
                                           var_type = resp_type,
                                           coords = c("X", "Y"), 
                                           response = "response", 
                                           plot_cor = TRUE)
    
    removed_var_df <- data.frame('Variable' = prepped_data$removed_vars, 'removal_step'= 'corr')
    
    # 2) Recursive variable reduction
    if(rfe){
      
      result <- slim_model_recursive(data_train = prepped_data$data_train,
                                     data_test = prepped_data$data_test,
                                     data_folds = prepped_data$data_folds,
                                     predictor_vars = prepped_data$predictor_vars,
                                     response  = "response",
                                     resp_type = resp_type)
      
      
      kept_vars = result$retained_vars
      removed_vars = result$removed_vars
      removed_vars_rfe = removed_vars[!removed_vars %in% removed_var_df$Variable]
      ranked_vars = result$ranked_vars
      removed_var_df <- bind_rows(removed_var_df, data.frame('Variable' = removed_vars_rfe, 'removal_step'= 'rfe'))
    }else{
      
      kept_vars = prepped_data$predictor_vars
      removed_vars = prepped_data$removed_vars
      ranked_vars = c()
    }
    
    # Create Initial Final fit
    final_fit <- final_fit_weighted_vars(kept_vars = kept_vars,
                                         removed_vars = removed_vars,
                                         rfe = rfe,
                                         ranked_vars = ranked_vars,
                                         data_train = prepped_data$data_train,
                                         data_test = prepped_data$data_test,
                                         data_folds = prepped_data$data_folds,
                                         hyper_tune = hyper_tune,
                                         variable_crosswalk = var_filtered,
                                         eig_nms = if(spatial){colnames(eigs)}else{NULL},
                                         retained_groups = retained_groups,
                                         retained_vars = retained_vars,
                                         resp_type = resp_type,
                                         cluster_id = group_id, 
                                         wto_outcome = wto_outcome,
                                         output_suffix = output_suffix, 
                                         spatial = spatial,
                                         sampling_dist = sampling_dist, 
                                         sp_subsampling = sp_subsampling)
    if(!sp_subsampling){
      morans_i = 0 # if sub-sampling is false stop while statement
    }else{
      if(is.null(final_fit$min_distance)){
        morans_i = 0 # no SA is found stop while statement
        write.csv(removed_var_df, paste0(CSV_DIR, '/', cluster_id, 'removed_vars', output_suffix,'.csv'), row.names = F)
      }else{
        
        # If there is still a distance where moran's I is > 0.1 then subsample data 
        df <-  interative_subsample(df, wto_outcome = wto_outcome, min_dist = final_fit$min_distance)
        print(paste0('data has been sub sampled to ', final_fit$min_distance))
        sampling_dist <- final_fit$min_distance
        morans_i = 1 # continue with a high morans I
      
        }
    }
    
  }
}
  
