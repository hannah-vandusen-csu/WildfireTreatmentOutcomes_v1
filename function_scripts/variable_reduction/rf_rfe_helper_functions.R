##
##
## Random Forest Functions 
## Used for RFE and matching
##
##

# functions
balance_binary <- function(df, analysis_step, wto_outcome){
  
  # Just balance based on Response
  if(analysis_step == 'Matching'){
    
    cntdf <- df %>% group_by(response) %>% summarise(cnt = n())
    minn <- min(cntdf$cnt)
    set.seed(1234)
    df <- df %>% 
      group_by(response) %>% 
      slice_sample(n = minn) %>% 
      ungroup()
    
    # make sure response is still a factor
    df[['response']] <- factor(df[['response']])
    return(df)
    
    
    # Balance based on paired points and response
  } else if(analysis_step == 'Predictor'){
    
    # If there is no untreated or treated points just balance regular
    if(length(unique(df$treated)) == 1){
      cntdf <- df %>% group_by(response) %>% summarise(cnt = n())
      minn <- min(cntdf$cnt)
      set.seed(1234)
      df <- df %>% 
        group_by(response) %>% 
        slice_sample(n = minn) %>% 
        ungroup()
      
      # make sure response is still a factor
      df[['response']] <- factor(df[['response']])
      return(df)
    }else{
      
  
    
    # if (wto_outcome=="infrastructure") {
    #   df <- df |> mutate(pair_id = paste0(cluster_id, '_', treatment_id))
    # } else if (wto_outcome=="containment" | wto_outcome=="burn_severity") {
    #   df <- df |> mutate(pair_id =  paste0(MTBS_ID, '_', treatment_id))
    # }
    
    set.seed(1234)
    valid_pairs <- df %>%
      group_by(matched_id) %>%
      summarise(treated_count = sum(treated == 1), 
                untreated_count = sum(treated == 0), 
                .groups = "drop") %>%
      filter(treated_count > 0 & untreated_count > 0) %>%
      pull(matched_id)
    
    df_valid <- df %>% filter(matched_id %in% valid_pairs)
    
    # Step 3: Assign each pair to response groups
    pair_response_counts <- df_valid %>%
      group_by(matched_id) %>%
      summarise(response_0 = sum(response == 0), 
                response_1 = sum(response == 1), .groups = "drop")
   
    # We can keep all pairs with already balanced counts 
    balanced_pairs <- pair_response_counts |> filter(response_0 == 1) |> pull(matched_id)
    
    # Sample points that are not balanced
    pairs_0 <- pair_response_counts |> filter(response_0 == 2) |> pull(matched_id)
    pairs_1 <- pair_response_counts |> filter(response_1 == 2) |> pull(matched_id)
    
    minn <- min(length(pairs_0), length(pairs_1))
    pairs_0_samp <- sample(pairs_0, minn)
    pairs_1_samp <- sample(pairs_1, minn)
    
    # Combine ids
    balanced_ids <- c(balanced_pairs, pairs_0_samp, pairs_1_samp)

    df_balanced <- df_valid %>%
      filter(matched_id %in% balanced_ids)
    
    df_balanced$response <- factor(df_balanced$response)
    return(df_balanced)
  }}
}

create_eigen <- function(df, num){
  # Create spatial layers: we have deem 30 spatial layers in sufficient
  if(dim(df)[1] > 3000){
    final_eigs <- meigen_f(as.matrix(dplyr::select(df, c("X","Y"))), enum = num)
  }else{
    final_eigs <- meigen(as.matrix(dplyr::select(df, c("X","Y"))), enum = num)
  }
  final_eigs <- as.data.frame(final_eigs$sf)
  eig_nms <- paste0("eig", rep(1:ncol(final_eigs)))
  names(final_eigs) <- eig_nms
  return(final_eigs)
}

simple_rf <- function(df, VAR_TYPE) {
  set.seed(1234)
  # Define the model mode and performance metric depending on a continuous or categorical variable
  if(VAR_TYPE == 'binary'){
    MODE <- 'classification'
    METRIC <- "accuracy"}
  if(VAR_TYPE == 'continuous'){
    MODE <- 'regression'
    METRIC <- "rsq"}
  
  # Split data
  colnames(df) <- c("response", "predictor")
  if (VAR_TYPE == "binary") {
    data_split <- initial_split(df, prop = 0.8, strata="response")
    }
  if (VAR_TYPE == "continuous") {
    data_split <- initial_split(df, prop = 0.8)
  }
  
  data_train <- training(data_split) 
  data_test <- testing(data_split)
  
  # Create model 
  model <-rand_forest(trees = 100) %>% 
    set_engine("ranger", importance = "permutation", num.threads = parallel::detectCores()/2-2) %>%
    set_mode(MODE)
  model_recipie <- recipe(response ~ predictor, data = data_train)
  model_fit <- workflow() %>% add_model(model) %>% add_recipe(model_recipie) %>% fit(data_train)
  
  # Test predictions 
  test_predictions <- data_test %>%
    dplyr::select(all_of('response')) %>%
    cbind(model_fit %>% predict(new_data = data_test %>% dplyr::select(-all_of('response')), type = ifelse(VAR_TYPE == 'binary', 'class', "numeric")))
  
  met <- get(METRIC)(test_predictions, truth = 'response', estimate = ifelse(MODE == 'classification', ".pred_class", ".pred"), event_level = "second")[3] %>% as.numeric()
  return(met) #returns the performance metric as a number
}

reduce_cor_vars <- function(data, vars, response, var_type, plot){
  set.seed(1234)
  
  numeric_data <- data[, vars]
  
  # remove variables with no variable
  numeric_distinct <- numeric_data |> dplyr::select(where(~ n_distinct(.) > 1))
  sd_vars <- colnames(numeric_data)[!colnames(numeric_data) %in% colnames(numeric_distinct)]
  
  if(dim(numeric_data)[2] > 0){
    fullcor <- cor(numeric_distinct)
    fullcor[is.na(fullcor)] <- 0
    predcor <- abs(fullcor)
    pred_nms <- colnames(predcor)
    
    for (i in 1:nrow(predcor)) {
      predcor[i,i] = 0
    }
    
    metrics <- c()
    progress_bar <- txtProgressBar(min = 0, max = length(vars), style = 3, char = "-")
    for (i in 1:ncol(numeric_distinct)) { #run an RF for each numerical predictor
      setTxtProgressBar(progress_bar, value = i)
      df <- cbind(data[, response, drop = F], numeric_distinct[,i, drop = F])
      this_mod <- simple_rf(df, VAR_TYPE = var_type)
      metrics <- c(metrics, this_mod)
    }
    names(metrics) <- colnames(numeric_distinct)
    close(progress_bar)
    
    removed_vars = c()
    
    while (max(predcor)>=0.6) { #for two correlated variables, select the var that can predict engagement better
      w <- which(predcor==max(predcor))[1]-1
      col <- (w %% ncol(predcor))+1
      row <- (w %/% ncol(predcor))+1
      colmetric <- metrics[names(metrics)==colnames(predcor)[col]]
      rowmetric <- metrics[names(metrics)==colnames(predcor)[row]]
      toremove <- ifelse(colmetric < rowmetric, col, row)
      removed_vars <- c(removed_vars,colnames(predcor)[toremove])
      predcor <- predcor[-toremove,-toremove]
    }
    
    # plot correlation matrix
    if(plot == TRUE){
      
      # Convert correlation matrix to long format using tidyr and filter for lower triangle
      melted_cor <- as.data.frame(fullcor) %>%
        mutate(Var1 = rownames(fullcor)) %>%
        pivot_longer(cols = -Var1, 
                     names_to = "Var2", 
                     values_to = "value") %>%
        # Convert Var1 and Var2 to same type for comparison
        mutate(
          Var1 = factor(Var1, levels = colnames(fullcor)),
          Var2 = factor(Var2, levels = colnames(fullcor))
        ) %>%
        # Keep only lower triangle (exclude diagonal)
        filter(as.numeric(Var1) > as.numeric(Var2))
      
      # Create the correlation plot with custom formatting
      ggplot(data = melted_cor, aes(x = Var1, y = Var2, fill = value)) +
        geom_tile() +
        scale_fill_gradient2(
          low = "#4575B4", 
          mid = "white",
          high = "#D73027",
          midpoint = 0,
          limits = c(-1,1),
          name = "Correlation"
        ) +
        theme_minimal() +
        theme(
          axis.text.x = element_text(
            angle = 45, 
            hjust = 1,
            size = ifelse(levels(factor(melted_cor$Var1)) %in% unlist(unique(removed_vars)), 10, 6),
            face = ifelse(levels(factor(melted_cor$Var1)) %in% unlist(unique(removed_vars)), "bold", "plain"),
            color = ifelse(levels(factor(melted_cor$Var1)) %in% unlist(unique(removed_vars)), "black", "grey50")
          ),
          axis.text.y = element_text(
            size = ifelse(levels(factor(melted_cor$Var2)) %in% unlist(unique(removed_vars)), 10, 6),
            face = ifelse(levels(factor(melted_cor$Var2)) %in% unlist(unique(removed_vars)), "bold", "plain"),
            color = ifelse(levels(factor(melted_cor$Var2)) %in% unlist(unique(removed_vars)), "black", "grey50")
          ),
          plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
        ) +
        coord_fixed() +
        labs(
          title = "Correlation Matrix Heatmap",
          subtitle = "Strong correlations (|r| > 0.6) shown in black",
          x = "",
          y = ""
        ) +
        # Add two separate geom_text layers for different formatting
        geom_text(data = subset(melted_cor, abs(value) <= 0.6),
                  aes(label = sprintf("%.2f", value)), 
                  color = "grey70", 
                  size = 1.5) +
        geom_text(data = subset(melted_cor, abs(value) > 0.6),
                  aes(label = sprintf("%.2f", value)), 
                  color = "black", 
                  size = 3, 
                  fontface = "bold")
      
      # Save the plot with appropriate dimensions
      ggsave(file.path(FIGURE_DIR,"temp_correlation_matrix.png"), width = 12, height = 10, dpi = 300)
    }
    # lists of removed and kept variables
    removed_vars <- c(unlist(unique(removed_vars)), sd_vars)
    kept_vars <- pred_nms[!pred_nms %in% removed_vars]
  }else{
    removed_vars = sd_vars
    kept_vars = c()
  }
  
  
  
  return(list(removed_vars = removed_vars,
              kept_vars = kept_vars))
}


prepare_data_for_model <- function(data, numeric_cols, factor_cols, exclude_cor_vars, run_cor, var_type, coords, response, plot_cor) {
  set.seed(1234)
  # Reduce correlated variables
  if (run_cor) {
    if (length(exclude_cor_vars) > 0 ){
      numeric_cols <- numeric_cols[!numeric_cols %in% exclude_cor_vars]
    }
    vars_reduced <- reduce_cor_vars(data = data, vars = numeric_cols, response = response, var_type = var_type, plot = plot_cor)
    removed_vars <- vars_reduced[[1]]
    kept_vars <- vars_reduced[[2]]
    
    # Combine predictors and remove unwanted ones
    if (length(exclude_cor_vars) > 0 ){
      predictors <- c(kept_vars, factor_cols, exclude_cor_vars)
    }else{
      predictors <- c(kept_vars, factor_cols)
    }
  } else {
    predictors <- c(numeric_cols, factor_cols)
    removed_vars <- NULL
    kept_vars <- NULL
  }
  
  # convert factor columns to factors
  data <-  data %>% 
    dplyr::mutate(across(all_of(factor_cols), factor))
  
  # Create a formula
  terms <- reformulate(termlabels = predictors, response = response)
  
  # Split data (ensuring training uses only the selected predictors)
  if (var_type == "binary") {
    data_split <- initial_split(data, prop = 0.8, strata="response")
  }
  if (var_type == "continuous") {
    
    data_split <- initial_split(data, prop = 0.8, strata="response_bin")
  }
  data_train <- training(data_split) #%>%
    #select(all_of(c(predictors, coords, response)))
  data_test <- testing(data_split)
  data_folds <- vfold_cv(data_train, v = 5)
  
  list(
    predictor_vars = predictors,
    terms = terms,
    data_train = data_train,
    data_test = data_test,
    data_folds = data_folds,
    removed_vars = removed_vars
  )
}

safe_tune_grid <- function(wf, folds, grid, metric) {
  tune_grid(
    wf,
    resamples = folds,
    grid = grid,
    metrics = metric_set(!!sym(metric)),
    control = control_grid(parallel_over = "resamples")
  )
}

run_the_model <- function(TRAIN, TEST, DATA_FOLDS, TERMS, PREDICTORS, VAR_TYPE, SAMPLE, TUNE = FALSE, LEVELS = 3, RESPONSE = "response"){
  set.seed(1234)
  # Define the model mode and performance metric
  if(VAR_TYPE == 'binary'){
    MODE <- 'classification'
    METRIC <- "accuracy"}
  if(VAR_TYPE == 'continuous'){
    MODE <- 'regression'
    METRIC <- "rsq"}

  
  formula <- recipe(TERMS, data = TRAIN) # Create the recipe for pre-processing
  
  if (SAMPLE == "under") { #undersample by treatment?
    formula <- formula %>%
      step_downsample(tx_class, skip = T)
  }
  
  # Grid search the hyper-parameters for tunning
  if (TUNE) {
    
    # Define the model with possible hyperparameter tuning
    model <- rand_forest(mtry = tune(), min_n = tune(), trees = 1000) %>% 
      set_engine("ranger", importance = "permutation") %>%
      set_mode(MODE)
    
    # Define the workflow
    workflow <- workflow() %>%
      add_model(model) %>%
      add_recipe(formula)
    
    rf_grid <- grid_regular(mtry(range = c(1, length(PREDICTORS))), 
                            min_n(range = c(5,50)), 
                            levels = LEVELS)
    # Perform parallel processing
    #system.time({
      
      # Parallelization isn't working...
      #cors <- parallel::detectCores()/2 # Detect the number of cores
      #cl <- parallel::makeCluster(cors)
      #doParallel::stopImplicitCluster()
      #doParallel::registerDoParallel(cl)


      # Tune the model with grid search
      # tune_fit <- workflow %>%
      #   tune_grid(resamples = DATA_FOLDS,
      #             grid = rf_grid,
      #             metrics = metric_set(!!sym(METRIC)),
      #             control = control_grid(parallel_over = "resamples"))
      tune_fit <- safe_tune_grid(workflow, DATA_FOLDS, rf_grid, METRIC)
      
      best_fit <- tune_fit %>%  # Select the best model based on the specified metric
        select_best(metric = METRIC)
      final_fit <- workflow %>% # Finalize the workflow with the best hyperparameters
        finalize_workflow(best_fit) %>%
        fit(data = TRAIN)
      #parallel::stopCluster(cl) # Stop the cluster
    #})
  } else { 
    
    # Add cores
    num_cores <- (detectCores()-2)/2
    model <- rand_forest() %>% # Just default 
      set_engine("ranger", importance = "permutation", num.threads = num_cores) %>%
      set_mode(MODE)
    
    workflow <- workflow() %>%
      add_model(model) %>%
      add_recipe(formula)
    
    final_fit <- workflow %>%
      #finalize_workflow(default_hyp) %>%
      fit(data = TRAIN)
    
  }
  
  vi <- final_fit %>% extract_fit_parsnip %>% vip(num_features = length(PREDICTORS))
  
  test_predictions <- TEST %>%
    dplyr::select(all_of(RESPONSE)) %>%
    cbind(final_fit %>% predict(new_data = TEST %>% dplyr::select(-all_of(RESPONSE)), ifelse(VAR_TYPE == 'binary', 'class', "numeric")))
  
  met <-  get(METRIC)(test_predictions, truth = !!sym(RESPONSE), estimate = ifelse(MODE == 'classification', ".pred_class", ".pred"), event_level = "second")[3] %>% as.numeric()
  
  #gc() #free up RAM
  return(list(final_fit = final_fit, tune_fit = if(TUNE){tune_fit}else{NA}, vi = vi, met = met))
}


slim_model_recursive <- function(data_train, data_test, data_folds, predictor_vars, response, resp_type, best_performance = NULL, removed_vars = c(), ranked_vars = NULL) {
  set.seed(1234)
  # Update the model formula with the current slim set of predictors
  terms <- reformulate(termlabels = predictor_vars, response = response)
  
  set.seed(1234)
  # Run random forest model without tuning hyper-parameters
  model_list <- run_the_model(TRAIN = data_train,
                              TEST = data_test,
                              DATA_FOLDS = data_folds,
                              TERMS = terms,
                              PREDICTORS = predictor_vars,
                              VAR_TYPE = resp_type, 
                              SAMPLE = 'none')
  
  set.seed(1234)
  tune_fit <- model_list$tune_fit
  set.seed(1234)
  final_fit <- model_list$final_fit
  set.seed(1234)
  final_vi <- model_list$vi
  set.seed(1234)
  final_vi_df <- final_vi$data %>% arrange(Importance) %>% filter(Importance >= 0) # removes negative importance variables
  set.seed(1234)
  performance <- model_list$met
  
  # Initialize best_performance if  NULL
  if (is.null(best_performance)) {
    set.seed(1234)
    best_performance <- performance
  }
  
  
  if (performance/best_performance > 0.95) {
    cat("    Slimming model run:", length(removed_vars))
    
    # Recalculate variable importance with the slimmer model
    ranked_vars <- final_vi_df
    removed_vars <- c(removed_vars, ranked_vars$Variable[1])
    set.seed(1234)
    predictor_vars <- ranked_vars %>% slice(-1) %>% pull(Variable)
    
    cat(", removing variable: ", tail(removed_vars, 1), "\n", sep = "")
    
    # Update the best performance if the current performance is better
    if (performance > best_performance) {
      best_performance <- performance
    }
    
    # Recursively call the function with updated parameters
    set.seed(1234)
    slim_model_recursive(data_train, data_test, data_folds, predictor_vars, response, resp_type, best_performance, removed_vars, ranked_vars)
  } else {
    cat(", model performance sub-optimal retaining last model variable weights\n")
    
    # Return the final model and removed variables
    return(list(final_fit = final_fit, retained_vars = predictor_vars, removed_vars = removed_vars, final_vi = final_vi, ranked_vars = ranked_vars))
  }
}

final_fit_weighted_vars <- function(kept_vars,
                                    removed_vars,
                                    rfe,
                                    ranked_vars, 
                                    data_train, 
                                    data_test, 
                                    data_folds,
                                    hyper_tune, 
                                    variable_crosswalk , 
                                    eig_nms, 
                                    retained_groups,
                                    retained_vars,
                                    resp_type, 
                                    cluster_id,
                                    wto_outcome,
                                    output_suffix, 
                                    spatial,
                                    sampling_dist,
                                    sp_subsampling = NULL){
  set.seed(1234)
  
  # Kept variables
  pred_vars <- kept_vars
  
  # Read in variable crosswalk and join in with importance table
  cwlks <- variable_crosswalk
  cwlks <- dplyr::select(cwlks, c("Group", "Abbreviation", "Figure_Name"))
  
  if(rfe){
    
    # Variable importance ranking (removed the first since that has been dropped)
    vidf <- ranked_vars[-1,]
    
    vidf <- left_join(vidf, cwlks, join_by("Variable" == "Abbreviation"))
    
    # Determine which Variable Groups are still in the model and which aren't 
    groups_df <- vidf |> filter(!is.na(Group))
    top_groups <- unique(groups_df$Group)
    needed_groups <- retained_groups[!retained_groups %in% top_groups]
    
    # If all variables in Group(s) were removed add the last variable in that Group back in
    # *note: if all variables in a group were removed in correlation then that group will not be included
    if(length(needed_groups) > 0){
      
      removed_vars <-  data.frame("Variable" = removed_vars, "Importance" = rev(1:length(removed_vars)))
      add_back_vars <- left_join(removed_vars, cwlks, join_by("Variable" == "Abbreviation")) |>
        filter(!is.na(Group)) |>
        filter(Group %in% needed_groups) |>
        group_by(Group) |>
        slice(which.min(Importance))
      
      pred_vars <- c(pred_vars, add_back_vars$Variable, retained_vars[!retained_vars %in% pred_vars])
    }else{
      
      pred_vars <- c(pred_vars, retained_vars[!retained_vars %in% pred_vars])}
    
    }else{
      
      # Add forced variables that aren't included in predictor variables
      pred_vars <- c(pred_vars, retained_vars[!retained_vars %in% pred_vars])
    }
  
  # Refit the model 
  
  # Run random forest model without tuning hyper-parameters
  terms <- reformulate(termlabels = pred_vars, response = 'response')
  set.seed(1234)
  final_mod <- run_the_model(TRAIN = data_train,
                             TEST = data_test,
                             DATA_FOLDS = data_folds,
                             TERMS = terms,
                             PREDICTORS = pred_vars,
                             VAR_TYPE = resp_type, 
                             SAMPLE = 'none', 
                             TUNE = hyper_tune,
                             LEVELS = 3)
  
  
  # Output generation
  # Define the model mode and performance metric
  if(resp_type == 'binary'){
    MODE <- 'classification'
    METRIC <- "accuracy"}
  if(resp_type == 'continuous'){
    MODE <- 'regression'
    METRIC <- "rsq"}
  
  # 1) Model Performance 
  test_met <- final_mod$met
  
  if(hyper_tune){
    args <- select_best(final_mod$tune_fit, metric = METRIC)
    mets <- c()
    for (i in 1:5) {
      this_df <-  final_mod$tune_fit$.metrics[[i]]
      this_df <- filter(this_df, mtry == args$mtry & min_n == args$min_n)
      mets <- c(mets, filter(this_df, .metric == METRIC)$.estimate)
    }
    
    mn_met <- mean(mets)
  }else{
    mn_met <- NA
  }
  
  
  mets <- data.frame('test_met' = test_met, 'mean_met' = mn_met, 'cluster_id' = cluster_id)
  
  # 2) Variable Importance  
  vidf <- final_mod$vi$data %>% arrange(Importance)
  vidf <- left_join(vidf, cwlks, join_by("Variable" == "Abbreviation"))
  
  if(spatial){
    vidf <- vidf %>%
      mutate(
        Group = if_else(str_detect(Variable, "^eig"), "Spatial", Group),
        Figure_Name = if_else(str_detect(Variable, "^eig"),
                              str_to_title(Variable),
                              Figure_Name))
  }
 
 
  # 2a) Variable Importance by weights for non-spatial variable
  vi_topvars <- vidf |>  filter(Importance > 0) %>% filter(Group %in% retained_groups) |> 
    mutate(proportional_weights = Importance/sum(Importance),
           cluster_id = cluster_id) 
  
  # set color palette
  #cols <- c('#2E8289FF','#91D5DEFF', '#5C4033', '#EAAE37FF', '#565F41FF', '#C94C24FF', '#6B4C9AFF','#009E73', '#F0E442', '#B6B4E9') #,
  cols <-  paletteer_d("MoMAColors::Warhol", n = length( unique(vidf$Group)))
  group_cols <- data.frame(group = unique(vidf$Group), color = cols)
  
  
  # Arrange vidf by Importance and set factor levels for Figure_Name
  vidf <- vidf %>%
    arrange(Importance) %>%
    filter(Importance > 0) %>% 
    mutate(Figure_Name = factor(Figure_Name, levels = unique(Figure_Name)))
  
  # Set factor levels for Group
  vidf$Group <- factor(vidf$Group, levels = unique(vidf$Group))
  # Variable Importance plots and pdp 
  (final_imps <- ggplot(vidf) +
      geom_col(aes(x = Figure_Name, y = Importance, fill = Group), width = 0.7) +
      coord_flip() +
      scale_fill_manual("Predictor Category", values = setNames(group_cols$color, group_cols$group))+ 
      theme_cowplot(16) +
      theme(axis.line.x  = element_blank(),
            axis.title.y  = element_blank(),
            axis.text.x=element_blank(),
            axis.line.y  = element_blank(),
            panel.background = element_rect(fill = 'transparent'),
            plot.background = element_rect(fill = 'transparent'),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = 'right',
            legend.title = element_text(hjust = 0.5),
            legend.justification = "left",
            axis.line = element_line(colour = "black", size=1.2),
            panel.border = element_rect(colour = "black", fill=NA, size=1.2),
            strip.background = element_rect(colour = "black", fill = "transparent", size = 1.2), 
            strip.text.x = element_text(colour = "black", face = "bold")) +
      guides(fill = guide_legend(ncol = 1), color = guide_legend(ncol = 1)))
  
  # Create PDP for top predictors 
  
  # Create an explainer using DALEXtra
  explainer <- explain_tidymodels(
    final_mod$final_fit,
    data = data_train %>% dplyr::select(-(!!sym('response'))),
    y = if(resp_type  == 'binary'){as.integer(pull(data_train['response']))}else{as.numeric(pull(data_train['response']))}, 
    label = "Random Forest"
  )
  
  # Generate a Partial Dependence Plot
  top_vars <- vi_topvars |> group_by(Group)|> slice_max(order_by = Importance, n = 1) |> dplyr::select(Variable)|> pull()
  pdp_rf <- model_profile(
    explainer,
    variables = top_vars,
    N = 500,
    type = "partial"
  )
  
  
  pdp_rf$agr_profiles$`_label_` = "partial dependence"
  
  vi_topvars_plot <- vi_topvars %>%
    mutate(Figure_Name = str_wrap(Figure_Name, width = 20))
  
  pdp_plot <- as_tibble(pdp_rf$agr_profiles) %>%
    mutate(`_label_` = stringr::str_remove(`_label_`, "random forest_")) %>%
    left_join(vi_topvars_plot, by = c('_vname_' = 'Variable')) |>
    ggplot(aes(`_x_`, `_yhat_`)) +
    geom_line(linewidth = 1.5, alpha = 0.8) +
    coord_cartesian(ylim = c(0,1)) +
    facet_wrap(~ Figure_Name, nrow = 3, scales = "free_x") +
    labs(y = "Predicted Probability", x = "", color = '', fill = '') +
    theme_cowplot(16) +
    theme(axis.line = element_line(colour = "black", size=1.2),
          strip.background = element_rect(colour = "black", fill = "grey", size = 1.2), 
          strip.text.x = element_text(colour = "black", face = "bold")) +
    guides(fill = guide_legend(nrow = 3), color = guide_legend(nrow = 3))+
    guides(fill = guide_legend(nrow = 3), color = guide_legend(nrow = 3))
  
  combined_plot <- final_imps + pdp_plot
  
  if(!dir.exists(paste0(FIGURE_DIR, '/var_importance'))){
    dir.create(paste0(FIGURE_DIR, '/var_importance')) 
  }
  
  cluster_id <- ifelse(cluster_id == '', '', paste0(cluster_id, "_"))
  ggsave2(file.path(FIGURE_DIR, 'var_importance', paste0(cluster_id, 'vi_pdp_matching', output_suffix, '.png')), width = 16,
          height = 8)
  
  # Moran's I plot to assess for spatial auto-correlation 
  resids_df <- data_train
  
  p_hat <- final_mod$final_fit %>% 
    predict(new_data = resids_df %>% dplyr::select(-all_of('response')), ifelse(resp_type == 'binary', 'prob', "numeric")) |>
    dplyr::select(ifelse(resp_type == 'binary',".pred_1", '.pred')) |> pull()
  
  resids <- as.numeric(as.character(resids_df$response)) - p_hat
  resids_df$resids <- resids
  
  sdf <- st_as_sf(resids_df, coords = c("X", "Y"), crs = 32610)
  #mapview(sdf[,"resids"])
  
  max_sample <- ifelse(sampling_dist < 200, 1000, ifelse(sampling_dist < 400, 2000, 3000))
  neighbs <- seq(sampling_dist , max_sample, sampling_dist)
  mis <- c()
  ps <- c()
  
  for (k in neighbs) {
    part_nb = dnearneigh(sdf, d1=0, d2=k)
    attributes(part_nb)$region_id<-as.character(1:length(part_nb))
    part_dlist <- nbdists(part_nb, sdf)
    part_wt2 = nb2listw(part_nb, glist=part_dlist, style="W", zero.policy=TRUE)
    
    m_out <- moran.test(sdf$resids, part_wt2, zero.policy = T)
    mis <- c(mis, unname(m_out$estimate[1]))
    ps <- c(ps, unname(m_out$p.value[1]))
  }
  
  ps_fact <- ifelse(ps < 0.05, "1", "0")
  
  final_morans <- data.frame("lag" = neighbs, "mi" = mis, "pval" = factor(ps_fact, levels = c(0,1)), "pval_num" = ps)
  
  ggplot(final_morans) +
    geom_line(aes(x = lag, y = mi)) +
    geom_point(aes(x = lag, y = mi, color = pval)) +
    scale_color_manual(values = c("blue", "red"), labels = c("0"=">=0.05", "1" = "<0.05")) +
    labs(x = "lag distances (m)", y = "Moran's I") +
    theme_bw()
  
  if(!dir.exists(paste0(FIGURE_DIR,'/morans'))){
    dir.create(paste0(FIGURE_DIR, '/morans')) 
  }
  ggsave2(file.path(FIGURE_DIR, '/morans/', paste0(cluster_id, 'moransI_matching', output_suffix, '.png')), width = 9,
          height = 7)
  
  if(!dir.exists(paste0(CSV_DIR, '/morans'))){
    dir.create(paste0(CSV_DIR, '/morans')) 
  }
  
  num_EVs <- vidf |> filter(Group == 'Spatial') |> nrow()
  final_morans |> mutate(cluster_id = cluster_id, final_EVs = num_EVs) |> 
    write.csv(paste0(CSV_DIR, '/morans/', cluster_id, 'moransI', output_suffix,'.csv'), row.names = F)
  
  if(!dir.exists(paste0(CSV_DIR, '/mets'))){
    dir.create(paste0(CSV_DIR, '/mets')) 
  }
  write.csv(mets, paste0(CSV_DIR, '/mets/', cluster_id, 'mets', output_suffix,'.csv'), row.names = F)
  
  # If Matching step write out the variable weights
  if(analysis_step == 'Matching'){
    if(!dir.exists(paste0(CSV_DIR, '/weights'))){
      dir.create(paste0(CSV_DIR, '/weights')) }
    write.csv(vi_topvars, paste0(CSV_DIR, '/weights/', cluster_id, 'weights', output_suffix,'.csv'), row.names = F)
  }
  
  # If predictor step write out the final dataframe
  if(analysis_step == 'Predictor'){
    
    # If no spatial sub sampling save final df and write model
    if(!sp_subsampling){
      
      df_whole <- bind_rows(data_train, data_test)
      
      preds_whole <- final_mod$final_fit  %>%
        predict(new_data = df_whole %>% dplyr::select(-all_of('response')), type = ifelse(resp_type == 'binary', 'prob', "numeric"))
      
      if(resp_type == 'binary'){
        df_whole$preds <- pull(preds_whole[,2])
      }else{
        df_whole$preds <- pull(preds_whole[,1])
      }
      
      final_df <-dplyr:: select(df_whole, any_of(c("response", "X", "Y", "lat", "long", "MTBS_ID", "point_id", "treatment_id",'matched_id' , 'preds', pred_vars)))
      write.csv(final_df, paste0(CSV_DIR, '/', cluster_id, 'final_df_', wto_outcome, output_suffix,'.csv'), row.names = F)
      
      if(!is.null(MODEL_DIR)){
        # Save the finalized workflow
        saveRDS(final_mod$final_fit, paste0(MODEL_DIR, "/final_workflow.rds"))}
      
      return(list(final_df = final_df))
    }else if(sp_subsampling){
      min_distance <- final_morans |> filter(mi >= 0.1 & pval_num < 0.05) |> dplyr::select(lag) |> pull() |> min() 
      
      # If the current distance is the smallest distance where SA matters, sub-sample to the next step
      if(min_distance == sampling_dist){
        next_dist <- final_morans[which(final_morans$lag == min_distance) + 1 , ]$lag
        min_distance <- next_dist
      }
      # If there is spatial autocorrelation return sampling distance
      if(min_distance != Inf){
        return(list(min_distance = min_distance)) 
      }else{
        df_whole <- bind_rows(data_train, data_test)
        final_df <- dplyr::select(df_whole, any_of(c(get(response), "X", "Y", "lat", "long", "MTBS_ID", "point_id", "treatment_id", "matched_id", pred_vars)))
        write.csv(final_df, paste0(CSV_DIR, '/', cluster_id, 'final_df_', wto_outcome, output_suffix,'.csv'), row.names = F)
        
        if(!is.null(MODEL_DIR)){
          # Save the finalized workflow
          saveRDS(final_mod$final_fit, paste0(MODEL_DIR, "/final_workflow.rds"))}
      }
    }
      
    }else{
    return(list(mets, vi_topvars))
  }
}

sub_sample_by_group <- function(f, min_dist){
  
  if (wto_outcome=="infrastructure") {
    df <- df |> filter(cluster_id == f)
  } else if (wto_outcome=="containment" | wto_outcome=="burn_severity") {
    df <- df |> filter(MTBS_ID == f)
  }
  
  # Convert to sf object
  points_sf <- st_as_sf(df, coords = c("X", "Y"), crs = st_crs(26910))
  
  # Create empty vector to store kept point_ids
  kept_ids <- c()
  
  # Get unique point_ids
  unique_pairs <- unique(points_sf$treatment_id)
  
  # Process one pair at a time
  for (ids in unique_pairs){
    current_pair_id <- ids
    current_pair <- points_sf |> filter(treatment_id == current_pair_id)
    
    # If this is the first pair, keep it
    if(length(kept_ids) == 0) {
      kept_ids <- c(kept_ids, current_pair_id)
    }else {
      # Get all previously kept points
      kept_points <- points_sf |> filter(treatment_id %in% kept_ids)
      
      # Calculate distances between current pair and all kept points
      distances <- drop_units(st_distance(current_pair, kept_points))
      # Convert self distance to Inf
      distances[distances == 0] <- Inf
      
      # If all distances are greater than min_distance, keep the pair
      if(all(distances > min_dist)) {
        kept_ids <- c(kept_ids, current_pair_id)
      }
    }
    
  }
  return(df |> filter(treatment_id %in% kept_ids))
}

interative_subsample <- function(df, wto_outcome, min_dist){
  if (wto_outcome=="infrastructure") {
    group_id <- unique(df$cluster_id)
  } else if (wto_outcome=="containment" | wto_outcome=="burn_severity") {
    group_id <-  unique(df$MTBS_ID)
  }
  df_list <- lapply(group_id, sub_sample_by_group, min_dist = min_dist)
  return(bind_rows(df_list))
}

add_treatment_id <- function(df, wto_outcome){
  if (wto_outcome=="infrastructure") {
    group_id <- unique(df$cluster_id)
  } else if (wto_outcome=="containment" | wto_outcome=="burn_severity") {
    group_id <-  unique(df$MTBS_ID)
  }
  df_list <- lapply(group_id, function(id){
    if (wto_outcome=="infrastructure") {
      df <- df |> filter(cluster_id == id)
    } else if (wto_outcome=="containment" | wto_outcome=="burn_severity") {
      df <- df |> filter(MTBS_ID == id)
    }
    
    
    df_treated <- df |> filter(treated == 1)
    df_untreated <- df |> filter(treated == 0)
    
    lapply(seq_along(1:dim(df_treated)[1]), function(row){
      trt_row <- df_treated[row,]
      
      paired_row <- df |>
        filter(treated == 0) |># untreated
        filter(control_id == trt_row$control_id) |> # Find the control for this treatment
        first() |> # select the first if duplicates
        bind_rows(trt_row) |>
        mutate(treatment_id = trt_row$point_id)
      return(paired_row)
    }) |> rbindlist()
  })
  df <- rbindlist(df_list)
  return(df)
}

sub_sample_by_group <- function(f){
  
  if (wto_outcome=="infrastructure") {
    df1 <- df |> filter(cluster_id == f)
  } else if (wto_outcome=="containment" | wto_outcome=="burn_severity") {
    df1 <- df |> filter(MTBS_ID == f)
  }
  
  # Convert to sf object
  points_sf <- st_as_sf(df1, coords = c("X", "Y"), crs = st_crs(26910))
  
  # Create empty vector to store kept point_ids
  kept_ids <- c()
  
  # Get unique point_ids
  unique_pairs <- unique(points_sf$treatment_id)
  
  # Process one pair at a time
  for (ids in unique_pairs){
    current_pair_id <- ids
    current_pair <- points_sf |> filter(treatment_id == current_pair_id)
    
    # If this is the first pair, keep it
    if(length(kept_ids) == 0) {
      kept_ids <- c(kept_ids, current_pair_id)
    }else {
      # Get all previously kept points
      kept_points <- points_sf |> filter(treatment_id %in% kept_ids)
      
      # Calculate distances between current pair and all kept points
      distances <- drop_units(st_distance(current_pair, kept_points))
      # Convert self distance to Inf
      distances[distances == 0] <- Inf
      
      # If all distances are greater than min_distance, keep the pair
      if(all(distances > min_dist)) {
        kept_ids <- c(kept_ids, current_pair_id)
      }
    }
    
  }
  return(df1 |> filter(treatment_id %in% kept_ids))
}

interative_subsample <- function(df, wto_outcome, min_dist){
  if (wto_outcome=="infrastructure") {
    group_id <- unique(df$cluster_id)
  } else if (wto_outcome=="containment" | wto_outcome=="burn_severity") {
    group_id <-  unique(df$MTBS_ID)
  }
  df_list <- lapply(group_id, sub_sample_by_group)
  return(bind_rows(df_list))
}


# Config based load fire data 
load_fire_data <- function() {
  
  # Construct data path
  data_path <- file.path(
    local,
    outcome,
    data_subdirectory,
    filename
  )
  
  # Load raw data
  df_raw <- fread(data_path)
  
  # Rename response column if needed
  if(response_var != 'response'){
    df_raw <- df_raw |> dplyr::select(-any_of('response'))
    df_raw <- dplyr::rename(df_raw, 'response' = !!sym(response_var))
  }
  
  # Add these lines before the return statement
  cat("\nFirst few rows of the loaded data frame:\n")
  print(head(df_raw))
  
  
  # Process data based on grouping
  if (!is.null(variable_crosswalk)) {
    # Check if crosswalk file exists
    if (!file.exists(variable_crosswalk)) {
      stop(sprintf("variable crosswalk files does not exist"))
    }
    
    var_df <- read.csv(variable_crosswalk)
    var_proj <- var_df[var_df$Abbreviation %in%  colnames(df_raw),]
    
    # Check if grouping column exists
    if (!analysis_step %in% colnames(var_proj)) {
      stop(sprintf("Column '%s' not found in crosswalk", analysis_step))
    }
    
    # Filter to desired column
    var_filtered <- var_proj  |> filter(!!sym(analysis_step) == 'yes'& !!sym(outcome) == TRUE)
    
    cat(paste0("\nVariables included ", analysis_step, " step:\n"))
    print(var_filtered$Abbreviation)
    
    return(list(
      df_raw = df_raw,
      var_filtered = var_filtered
    ))
  }else{
    return(list(
      df_raw = df_raw
    ))
  }
  
}
