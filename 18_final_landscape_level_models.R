
# List of required packages
packages <- c("sf", 
              "terra",
              'ggpubr', 
              'paletteer', 
              'patchwork', 
              'ggridges', 
              'future.apply', 
              'future',
              'ingredients',
              'gridExtra',
              'grid', 
              'parallel', 
              'cowplot',
              'tidymodels', 
              'vip', 
              'DALEXtra', 
              'RColorBrewer', 
              "tidyverse", 
              "ggplot2", 
              "rdist", 
              "leaflet",
              'R.utils',
              'rsample',
              'doParallel')


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




# User Inputs -----------------------------------------------------------------
region <- 'east' # set to blank if it is the entire dataset
fig_suffix <- region

# Name of response
response = "response"

# Treatment (or any categorical you want to split the data by)
tr_var <- if (wto_outcome == "infrastructure") "class_cod.uw" else "class_cod"
age_var <- "tr_age"  

# number of bootstraps, probably should be 1000 for final runs
n_boots <- 1000

# Define the treatment age constraints 
untreated_label <- 0 
untreated_age <- 21
treated_age_range <- 2:20

#----------------- Load in input and output directories ------------------------
# Final dataframe directory 
df_dir <- file.path(processed_data_directory, 
                    'processed_outputs', wto_outcome, 
                    'final_dataframes', region, 'csv')

# Model directory
model_dir <- file.path(processed_data_directory, 'processed_outputs',
                       wto_outcome, "final_dataframes", region, "models")

# Figure directory 
fig_dir <- file.path(processed_data_directory, 'processed_outputs',
                     wto_outcome, "final_dataframes", region, "figures")

# Variable Crosswalk csv 
variable_crosswalk <- file.path(global_data_directory, 'variable_descriptions/full_variable_list.csv') 

# Color and naming schema for treatment classes
class_code <- file.path(global_data_directory, 'treatments/treatment_categories/treatment_code_naming_test.csv')

# Set output directories 

df_dir_output <- df_dir

# Figure output directory
FIGURE_DIR <- fig_dir

# CSV output directory
CSV_DIR <- df_dir

# Model output directory
MODEL_DIR <- model_dir


# Read in necessary data ---------------------------------------------------------

# Read in dataframe 
df <- read.csv(file.path(df_dir, paste0('final_df_', wto_outcome, 
                                        ifelse(response_type == 'binary', '_response_bin', '_response' ), '.csv')))

# Read in Model fit
mod <- readRDS(file.path(MODEL_DIR, 'final_workflow.rds'))

# Read variable crosswalk 
variable_crosswalk <- read.csv(variable_crosswalk) |> dplyr::filter(Predictor == "yes", burn_severity == T)

# Read in naming for treatments 
class_code <- read.csv(class_code) |> dplyr::filter(class_cod %in% df$class_cod)
class_code$class_cod <- as.factor(class_code$class_cod)


# predictors to drop
# in case you want to drop some variable for this model run that was forced into the dataframe
pred_drops <- c("treated", "p_uwt_all_2000", 'MTBS_ID')

# Create custom theme function
theme_hannah <- function(base_size = 16) {
  theme_cowplot(base_size) +
    theme(
      axis.line.x = element_blank(),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
      axis.line.y = element_blank(),
      panel.background = element_rect(fill = 'transparent'),
      plot.background = element_rect(fill = 'transparent'),
      panel.grid.major = element_blank(),
      #panel.grid.major = element_line(color = "grey90", size = 0.5),                          
      panel.grid.minor = element_blank(),
      legend.position = 'right',
      legend.title = element_text(hjust = 0.5),
      legend.justification = "center",
      axis.line = element_line(colour = "black", size = 1.2),
      panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
      strip.background = element_rect(colour = "black", fill = "transparent", size = 1.2),
      strip.text.x = element_text(colour = "black", face = "bold"),
      aspect.ratio = 6/5
    )
}
#-----------------PROCESSING-------------------#

# A) Dataframe manipulation ----------------------------------------------------

preds <- colnames(df)
preds <- preds[preds %in% variable_crosswalk$Abbreviation]
preds <- preds[!preds %in% pred_drops]

fvars <- dplyr::filter(variable_crosswalk, Continuous == "no")$Abbreviation

fvars <- fvars[fvars %in% preds]

cvars <- dplyr::filter(variable_crosswalk, Continuous == "yes")$Abbreviation
cvars <- cvars[cvars %in% preds]

mod_df <- df %>% dplyr::select(all_of(c("response", "X", "Y", preds))) %>%
  mutate(across(c(fvars, response), as.factor))

# Extract hyper-parameters for model fitting in variable reduction
system.time({final_fit <- mod |> fit(mod_df)})
trees <- as.numeric(quo_name(final_fit$fit$fit$spec$args$trees))
mtry <- as.numeric(quo_name(final_fit$fit$fit$spec$args$mtry))
min_n <- as.numeric(quo_name(final_fit$fit$fit$spec$args$min_n))

hyps <- tibble(trees=trees, mtry=mtry, min_n=min_n)

# B) Bootstrap Model fitting to create Confidence Intervals---------------------

# Create bootstraps

# 1) paste all categorical variables to define strata
if (length(fvars) > 0) {
  mod_df$strata <- do.call(paste, c(mod_df[c(fvars, response)], sep = "_"))
} else {
  stop("fvars must contain at least one variable name")
}

# 2) create bootstraps
boots <- bootstraps(mod_df, times = n_boots,  strata = "strata")

# fit each bootstrap to the workflow 
wf_final <- mod %>% finalize_workflow(hyps)

options(future.globals.maxSize = 2 * 1024^3)  
plan(cluster, workers = parallel::detectCores() - 2)


safe_boot_fit <- function(split, idx, timeout = 300) {
  tryCatch({
    res <- withTimeout({
      wf_final %>% fit(analysis(split))
    }, timeout = timeout, onTimeout = "error")
    
    list(fit = res, index = idx, error = NULL)
    
  }, error = function(e) {
    list(fit = NULL, index = idx, error = conditionMessage(e))
  })
}

out_dir <- file.path(df_dir_output, paste0("bootstraps", n_boots), "boot_fits")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Attach index safely
splits <- Map(
  function(split, i) {
    attr(split, "boot_id") <- i
    split
  },
  boots$splits,
  seq_along(boots$splits)
)

bootstrap_results <- future_lapply(
  splits,
  function(split) {
    
    idx <- attr(split, "boot_id")
    
    res <- safe_boot_fit(split, idx)
    
    saveRDS(
      res,
      file = sprintf("%s/fit_%04d.rds", out_dir, idx)
    )
    
    list(index = idx, error = res$error)
  },
  future.seed = TRUE
)


# C) Variable importance and Partial dependence plots with CI's-----------------

# Base explainer using the first fit
fit1 <- readRDS(
  sprintf("%sbootstraps%d/boot_fits/fit_%04d.rds",
          paste0(df_dir_output, '/'), n_boots, 1)
)$fit

base_explainer <- explain_tidymodels(
  model = fit1,
  data  = mod_df %>% dplyr::select(-all_of(response)),
  y     = if (response_type == "binary") {
    as.integer(mod_df[[response]])
  } else {
    mod_df[[response]]
  },
  verbose = FALSE
)


base_data <- base_explainer$data
all_classes <- unique(base_data[[tr_var]])

process_bootstrap <- function(i, base_explainer) {
  tryCatch(
    {
      # --- Load fit ---
      this_fit <- readRDS(
        sprintf("%sbootstraps%d/boot_fits/fit_%04d.rds",
                paste0(df_dir_output, '/'), n_boots, i)
      )$fit
      
      # --- Variable importance ---
      vidf <- vip::vip(
        this_fit,
        num_features = length(preds)
      )$data |>
        dplyr::arrange(tolower(Variable)) |>
        dplyr::select(Importance)
      
      names(vidf) <- paste0("boot", i)
      
      # --- Reuse explainer ---
      this_explainer <- base_explainer
      this_explainer$model <- this_fit
      this_data <- base_data

      
      N_FAST <- 100
      
     
      
      # --- Treatment PDP (manual ICE-style, vectorized over classes) ---
      tr_pdp <- dplyr::bind_rows(lapply(all_classes, function(cls) {
        
        df_subset <- if (nrow(this_data) > N_FAST) {
          this_data |> dplyr::slice_sample(n = N_FAST)
        } else {
          this_data
        }
        
        # --- Set treatment class to the one we're profiling ---
        df_subset[[tr_var]] <- factor(cls, levels = levels(this_data[[tr_var]]))
        
        # --- Assign realistic age given treatment class ---
        if (cls == untreated_label) {
          df_subset[[age_var]] <- untreated_age
        } 
        if (cls != untreated_label) {
          tr_age_sample <- this_data |>
            dplyr::filter(
              .data[[age_var]] >= min(treated_age_range),
              .data[[age_var]] <= max(treated_age_range)
            ) |>
            dplyr::pull(!!sym(age_var)) |>
            sample(size = nrow(df_subset))
          
          df_subset[[age_var]] <- tr_age_sample
        }
        
        # --- Mean prediction across the background ---
        avg_pred <- mean(
          dplyr::pull(predict(this_fit, df_subset, 'prob')[, 2]),
          na.rm = TRUE
        )
        
        tibble(var = tr_var, var_val = cls, pred = avg_pred)
      }))
      
      # --- Continuous PDPs ---
      profile <- model_profile(
        this_explainer,
        variables = cvars,
        N = N_FAST,
        groups = tr_var
      )
      
      pdp_df <- as_tibble(profile$agr_profiles)
      names(pdp_df) <- c("var", "lab", "var_val", "group", "pred", "id")
      
      pdp_df <- pdp_df |>
        dplyr::select(var, var_val, group, pred) |>
        dplyr::arrange(var, group, var_val)
      
      # --- Success return ---
      list(
        ok     = TRUE,
        boot   = i,
        vidf   = vidf,
        tr_pdp = tr_pdp,
        pdp_df = pdp_df,
        error  = NULL
      )
    },
    
    error = function(e) {
      message(sprintf("❌ Bootstrap %d failed: %s", i, e$message))
      list(
        ok     = FALSE,
        boot   = i,
        vidf   = NULL,
        tr_pdp = NULL,
        pdp_df = NULL,
        error  = e$message
      )
    }
  )
}


# Start parallel processing using mclapply (won't work without mac/linux)
system.time({
  results <-mclapply(seq_len(n_boots), function(i){
    process_bootstrap(i, base_explainer)})
})


vis <- data.frame("var" = sort(preds))
for (i in seq_along(results)) {
  if(is.null(results[[i]]$vidf)){
    next
  }
  vis <- cbind(vis, results[[i]]$vidf)
}

tr_pdps <- lapply(results, function(x) x$tr_pdp)
pdp_dfs <- lapply(results, function(x) x$pdp_df)

# Save bootsraps as outputs: 
dir.create(file.path(df_dir_output, paste0('bootstraps', n_boots )))

# Variable Importance
write_csv(vis, paste0(df_dir_output, '/bootstraps', n_boots ,'/vi_bootsraps.csv'))
vis <- read_csv(paste0(df_dir_output, '/bootstraps', n_boots ,'/vi_bootsraps.csv'))

# tr pdps 
tr_df_bind <- bind_rows(tr_pdps) 
write_csv(tr_df_bind, paste0(df_dir_output, '/bootstraps', n_boots ,'/tr_mean_bootsraps.csv')) # this isn't work
tr_df_bind <- read_csv(paste0(df_dir_output, '/bootstraps', n_boots ,'/tr_mean_bootsraps.csv')) 

# pdps 
pdp_df_bind <- bind_rows(pdp_dfs) 
write_csv(pdp_df_bind, paste0(df_dir_output, '/bootstraps', n_boots ,'/pdp_bootsraps.csv')) # this isn't work
pdp_df_bind <- read_csv(paste0(df_dir_output, '/bootstraps', n_boots ,'/pdp_bootsraps.csv')) # this isn't work

# 1) Variable Importance --------------------------------------------------
vi_cis <- vis %>% 
  pivot_longer(!var) %>%
  group_by(var) %>%
  dplyr::summarise(lwr = quantile(value, 0.0275),
            median = median(value),
            upr = quantile(value, 0.975))

vi_cis <- left_join(vi_cis, variable_crosswalk, join_by("var"=="Abbreviation"))
vi_cis <- arrange(vi_cis, median)
vi_cis$Figure_Name <- factor(vi_cis$Figure_Name, levels = vi_cis$Figure_Name)

cols <-  paletteer_d("MoMAColors::Klein", n = length( unique(vi_cis$Group)))
group_cols <- data.frame(group = c("Fire Stats",
                                   "Short-term Climate", 
                                   "Topography", 
                                   "Weather", 
                                   "Long-term Climate",
                                   "Treatments"), 
                         color = c('#E45A5A',
                                   '#678096',
                                   '#803233',
                                   '#F2D56F',
                                   '#2BBBA4',
                                   '#80792B'))
(final_imps <- ggplot(vi_cis) +
    geom_col(aes(x = Figure_Name, y = median*100, fill = Group), width = 0.8) +
    geom_point(aes(x = Figure_Name, y = median*100)) +
    geom_errorbar(aes(x = Figure_Name, ymin=lwr*100, ymax=upr*100), width = 0, linewidth = 0.7) +
    scale_fill_manual("Predictor Category", values = setNames(group_cols$color, group_cols$group))+ 
    theme_bw() +
    #coord_flip() +
    labs(y = "% Change in Accuracy", x = '') +
    theme_hannah(16))
ggsave(file.path(FIGURE_DIR,paste0("vi_", fig_suffix,".png")), final_imps,width = 8, height =9, dpi = 300)
write_csv(vi_cis, paste0(df_dir_output, '/bootstraps', n_boots ,'/vi_cis.csv'))

# 2) Partial Dependence of treatment class -------------------------------------
tr_pdp_ci <- tr_df_bind |>
  group_by(var_val) %>%
  dplyr::summarise(lwr = quantile(pred, 0.025, na.rm = T),
            median = mean(pred, na.rm = T),
            upr = quantile(pred, 0.975, na.rm = T))
tr_pdp_ci$var_val <- as.factor(tr_pdp_ci$var_val)

class_code$class_cod <- as.factor(class_code$class_cod)
tr_pdp_ci <- left_join(tr_pdp_ci, class_code, join_by("var_val" == "class_cod"))

tr_pdp_ci <- tr_pdp_ci |> arrange(intensity_rank) |> mutate(class_name = factor(class_name, levels = class_name))

(mean_preds <- tr_pdp_ci |> 
    ggplot() +
    geom_hline(yintercept = 0.5, linewidth = 1, color = 'grey', linetype='dashed') +
    geom_point(aes(x = class_name, y = median, color = class_name),
               size = 4, shape = 16, position=position_dodge(width=0.5)) +
    geom_errorbar(aes(x = class_name, ymin = lwr, ymax = upr, color = class_name),
                  linewidth = 1.5, width=0, position=position_dodge(width=0.5)) +
    labs(y = "Predicted Probability \nof High Severity", x = "", color = '', fill = '') +
    coord_cartesian(ylim = c(0,1)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0.2, 0.7)) +
    scale_color_manual(values = c(tr_pdp_ci$treatment_color_hvd), 
                       labels = c(tr_pdp_ci$class_name)) +
    theme_hannah(16) +
    coord_flip())
ggsave(file.path(FIGURE_DIR,paste0("class_code_mean_prob_pdp",fig_suffix,".png")), mean_preds,width = 10.5, height =5, dpi = 300)
write_csv(tr_pdp_ci, paste0(df_dir_output, '/bootstraps', n_boots ,'/tr_pdps_cis.csv'))

# 3) Partial Dependence Plots --------------------------------------------------

pdp_cis <- pdp_df_bind |>  
  # mutate(pred = 1 - pred) |> # first calculate the right prediction 
  group_by(var, var_val, group) %>%
  dplyr::summarise(lwr = quantile(pred, 0.025, na.rm = T),
            median = median(pred, na.rm = T),
            upr = quantile(pred, 0.975, na.rm = T))

pdp_cis <- left_join(pdp_cis, variable_crosswalk, join_by("var"=="Abbreviation"))# get treatment names
class_code$class_cod <- as.factor(class_code$class_cod)
pdp_cis$group <- as.factor(pdp_cis$group)
pdp_cis <- left_join(pdp_cis, class_code, join_by("group" == "class_cod"))
write_csv(pdp_cis, paste0(df_dir_output, '/bootstraps', n_boots ,'/pdp_cis'))
# create plots for all top predictors

all_plts <- list()
for (i in 1:(length(cvars))) {
  this_df <- dplyr::filter(pdp_cis, var == cvars[i])
  
  if(cvars[i] != 'tr_age'){
    minval <- quantile(this_df$var_val, 0.025)
    maxval <- quantile(this_df$var_val, 0.975)
    this_df <- filter(this_df, var_val > minval & var_val < maxval)
  }
  
  classes <- unique(this_df$group)[-1]
  for (j in 1:length(classes)) {
    if(cvars[i] != 'tr_age'){
      
      
      this_df_plot <-  this_df |> filter(group == classes[j]| group == 0 )
      ymin <- ifelse(min(this_df_plot$lwr)< 0.2,0, 0.2)
      ymax <- ifelse( max(this_df_plot$upr)> 0.8, 1, 0.8)
      (this_pdp <- this_df_plot |> 
        ggplot() +
        geom_line(aes(x=var_val , y=median, color = class_name), linewidth = 1, alpha = 0.8) +
        geom_ribbon(aes(x=var_val, ymin=lwr, ymax=upr, fill = class_name), alpha = 0.2) +
        labs(x = paste0(first(this_df$Figure_Name), ''),y ="Predicted Probability \nof High Severity",color = NULL, fill = NULL) +
        scale_color_manual(values = setNames(unique(this_df_plot$treatment_color_hvd), 
                                             unique(this_df_plot$class_name))) +
        scale_fill_manual(values = setNames(unique(this_df_plot$treatment_color_hvd), 
                                            unique(this_df_plot$class_name)))+
        scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(ymin, ymax)) +
        theme_hannah(15))
      print(this_pdp)
      
      this_full_df <- dplyr::select(mod_df, c(cvars[i], 'class_cod')) |> filter(class_cod == as.integer(classes[j])| class_cod == 0)
      names(this_full_df)[1] <- "x_val"
      this_hist <- ggplot(this_full_df) +
        geom_density(aes(x = x_val, fill = class_cod), color = NA, alpha = 0.3) +
        scale_fill_manual(values = c(unique(this_df_plot$treatment_color_hvd)), 
                          labels = c(unique(this_df_plot$class_name))) +
        theme(axis.text = element_blank(),
              axis.title = element_blank(),
              axis.ticks = element_blank(),
              panel.background = element_rect(fill="white"),
              plot.background = element_rect(fill="white", color=NA),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              legend.position = "none",
              plot.margin = unit(c(0,.2,0,.55), "cm"))
      
      fplt <- as_ggplot(grid.arrange(this_hist, this_pdp, nrow = 2, heights = c(1,5)))
      ggsave(file.path(FIGURE_DIR, paste0(cvars[i], "_",classes[j],"_", fig_suffix,".png")),fplt , dpi = 300, width = 8, height =  6)
      all_plts[[i]] <- fplt
    }else{
      
      
      this_con <- dplyr::filter(this_df, group == 0 & var_val == 21) 
      
      this_df_plot <-  this_df |> filter(group == classes[j] & var_val != 21)
      
      ymin <- ifelse(min(this_df_plot$lwr)< 0.2,0, 0.2)
      ymax <- ifelse( max(this_df_plot$upr)> 0.8, 1, 0.8)
      
      xmin <- min(this_df_plot$var_val)
      xmax <- max(this_df_plot$var_val)
      
      (this_pdp <- this_df_plot |> 
          ggplot() +
          geom_line(aes(x=var_val , y=median, color = class_name), linewidth = 1, alpha = 0.8) +
          geom_ribbon(aes(x=var_val, ymin=lwr, ymax=upr, fill = class_name), alpha = 0.2) +
          annotate("rect",
                   xmin = -Inf, xmax = Inf,
                   ymin = this_con$lwr, ymax = this_con$upr,
                   fill = this_con$treatment_color_hvd, alpha = 0.1) +
          geom_hline(yintercept = this_con$median,
                     color = this_con$treatment_color_hvd, linetype = "dashed", linewidth = 1.5) +
          labs(x = paste0(first(this_df$Figure_Name), ''),y ="Predicted Probability \nof High Severity",color = NULL, fill = NULL) +
          scale_color_manual(values = setNames(unique(this_df_plot$treatment_color_hvd), 
                                               unique(this_df_plot$class_name))) +
          scale_fill_manual(values = setNames(unique(this_df_plot$treatment_color_hvd), 
                                              unique(this_df_plot$class_name)))+
          scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(ymin, ymax)) +
          scale_x_continuous(limits = c(xmin, xmax), expand = c(0, 0)) +
          theme_hannah(15))
      print(this_pdp)
      
      this_full_df <- dplyr::select(mod_df, c(cvars[i], 'class_cod')) |> filter(class_cod == as.integer(classes[j]))
      names(this_full_df)[1] <- "x_val"
      this_hist <- ggplot(this_full_df) +
        geom_density(aes(x = x_val, fill = class_cod), color = NA, alpha = 0.3) +
        scale_fill_manual(values = c(unique(this_df_plot$treatment_color_hvd)), 
                          labels = c(unique(this_df_plot$class_name))) +
        theme(axis.text = element_blank(),
              axis.title = element_blank(),
              axis.ticks = element_blank(),
              panel.background = element_rect(fill="white"),
              plot.background = element_rect(fill="white", color=NA),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              legend.position = "none",
              plot.margin = unit(c(0,.2,0,.55), "cm"))
      
      fplt <- as_ggplot(grid.arrange(this_hist, this_pdp, nrow = 2, heights = c(1,5)))
      ggsave(file.path(FIGURE_DIR, paste0(cvars[i], "_",classes[j],"_", fig_suffix,".png")), fplt , dpi = 300, width = 8, height =  6)
      all_plts[[i]] <- fplt
    }
  }
  
  
}#


rcvars <- rev(vi_cis$var[vi_cis$var %in% cvars])

all_plts <- list()
for (i in 1:(length(rcvars))) {
  this_df <- dplyr::filter(pdp_cis, var == rcvars[i])
  
  if(rcvars[i] != 'tr_age'){
    
    minval <- quantile(this_df$var_val, 0.025)
    maxval <- quantile(this_df$var_val, 0.975)
    this_df <- filter(this_df, var_val > minval & var_val < maxval)
    
    classes <- unique(this_df$group)[-1]
    
    this_df_plot <-  this_df 
    ymin <- 0
    ymax <- 0.8
    # ymin <- ifelse(min(this_df_plot$lwr)< 0.2,0, 0.2)
    # ymax <- ifelse( max(this_df_plot$upr)> 0.8, 1, 0.8)
    
    levs <- unique(this_df_plot$class_name)
    levs <- setdiff(levs, c("Untreated", "Rearrangement + Fire")) # Just some order
    new_levels <- c(levs, "Rearrangement + Fire", "Untreated")
    this_df_plot <- this_df_plot %>%
      mutate(class_name = factor(class_name, levels = new_levels))
    
    xlims <- range(this_df_plot$var_val)
    
    (this_pdp <- this_df_plot |> 
        ggplot() +
        geom_ribbon(aes(x=var_val, ymin=lwr, ymax=upr, fill = class_name), alpha = 0.2) +
        geom_line(aes(x=var_val , y=median, color = class_name), linewidth = 1.5, alpha = 1) +
        labs(x = paste0(str_wrap(first(this_df$Figure_Name), width = 25), ''),y ="Predicted Probability \nof High Severity",color = NULL, fill = NULL) +
        scale_color_manual(values = setNames(unique(this_df_plot$treatment_color_hvd), 
                                             unique(this_df_plot$class_name))) +
        scale_fill_manual(values = setNames(unique(this_df_plot$treatment_color_hvd), 
                                            unique(this_df_plot$class_name)))+
        scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
        scale_x_continuous(limits = xlims, expand = c(0, 0)) +
        coord_cartesian(ylim = c(ymin, ymax)) +
        theme_hannah(15))
    print(this_pdp)
    
    # Get the legend of the pdp  
    pdp_leg <- get_legend(this_pdp)
    
    this_full_df <- dplyr::select(mod_df, c(rcvars[i], 'class_cod')) 
    this_full_df <- left_join(this_full_df, class_code)
    names(this_full_df)[1] <- "x_val"
    this_hist <- ggplot(this_full_df) +
      geom_density(aes(x = x_val, fill = class_name), color = NA, alpha = 0.3) +
      scale_fill_manual(values = setNames(unique(this_df_plot$treatment_color_hvd), 
                                          unique(this_df_plot$class_name))) + 
      scale_x_continuous(limits = xlims) +
      theme(axis.text = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            panel.background = element_rect(fill="white"),
            plot.background = element_rect(fill="white", color=NA),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = "none",
            plot.margin = unit(c(0,.2,0,.55), "cm"))
    this_hist <- ggplot(this_full_df) +
      geom_density_ridges(aes(x = x_val, y = class_name, fill = class_name), 
                          color = NA, alpha = 0.5, jittered_points = TRUE) +
      scale_fill_manual(values = setNames(unique(this_df_plot$treatment_color_hvd), 
                                          unique(this_df_plot$class_name))) + 
      scale_x_continuous(limits = xlims) +
      theme(axis.text = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            panel.background = element_rect(fill="white"),
            plot.background = element_rect(fill="white", color=NA),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = "none",
            plot.margin = unit(c(0,.2,0,.55), "cm"))
    
    combined_plot <- this_hist / this_pdp + 
      plot_layout(heights = c(1, 5), guides = "collect")
    
    # Save with appropriate width/height (increase width if legend is wide)
    # ggsave(file.path(FIGURE_DIR, paste0(cvars[i], "_", fig_suffix,"2.png")),
    #        combined_plot,
    #        dpi = 300,
    #        width = 8, height = 6)
    # 
    all_plts[[i]] <- combined_plot
  }else{
    
    # calculate the control mean value 
    this_con <- dplyr::filter(this_df, group == 0 & var_val == 21) 
    
    minval <- quantile(this_df$var_val, 0.025)
    maxval <- quantile(this_df$var_val, 0.975)
    this_df <- filter(this_df, var_val > minval & var_val < maxval)
    
    classes <- unique(this_df$group)[-1]
    
    this_df_plot <-  this_df 
    ymin <- 0
    ymax <- 0.8
    # ymin <- ifelse(min(this_df_plot$lwr)< 0.2,0, 0.2)
    # ymax <- ifelse( max(this_df_plot$upr)> 0.8, 1, 0.8)
    
    this_df_plot <- dplyr::filter(this_df_plot, group != 0)

    
    levs <- unique(this_df_plot$class_name)
    levs <- setdiff(levs, c("Untreated", "Rearrangement + Fire")) # Just some order
    new_levels <- c(levs, "Rearrangement + Fire", "Untreated")
    this_df_plot <- this_df_plot %>%
      mutate(class_name = factor(class_name, levels = new_levels))
    
    xlims <- range(this_df_plot$var_val)
    
    (this_pdp <- this_df_plot |> 
        ggplot() +
        geom_ribbon(aes(x=var_val, ymin=lwr, ymax=upr, fill = class_name), alpha = 0.2) +
        geom_line(aes(x=var_val , y=median, color = class_name), linewidth = 1.5, alpha = 0.9) +
        annotate("rect",
                 xmin = -Inf, xmax = Inf,
                 ymin = this_con$lwr, ymax = this_con$upr,
                 fill = this_con$treatment_color_hvd, alpha = 0.1) +
        geom_hline(yintercept = this_con$median,
                   color = this_con$treatment_color_hvd, linetype = "dashed", linewidth = 1.5) +
        labs(x = paste0(first(this_df$Figure_Name), ''),y ="Predicted Probability \nof High Severity",color = NULL, fill = NULL) +
        scale_color_manual(values = setNames(unique(this_df_plot$treatment_color_hvd), 
                                             unique(this_df_plot$class_name))) +
        scale_fill_manual(values = setNames(unique(this_df_plot$treatment_color_hvd), 
                                            unique(this_df_plot$class_name)))+
        scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
        scale_x_continuous(limits = xlims, expand = c(0, 0)) +
        coord_cartesian(ylim = c(ymin, ymax)) +
        theme_hannah(15))
    print(this_pdp)
    
    # Get the legend of the pdp  
    pdp_leg <- get_legend(this_pdp)
    
    # Filter data
    this_full_df <- dplyr::select(mod_df, c(rcvars[i], 'class_cod')) |> filter(class_cod != 0)
    this_full_df <- left_join(this_full_df, class_code)
    names(this_full_df)[1] <- "x_val"

    this_hist <- ggplot(this_full_df) +
      geom_density_ridges(aes(x = x_val, y = class_name, fill = class_name), 
                          color = NA, alpha = 0.5, jittered_points = TRUE) +
      scale_fill_manual(values = setNames(unique(this_df_plot$treatment_color_hvd), 
                                          unique(this_df_plot$class_name))) + 
      scale_x_continuous(limits = xlims) +
      theme(axis.text = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            panel.background = element_rect(fill="white"),
            plot.background = element_rect(fill="white", color=NA),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = "none",
            plot.margin = unit(c(0,.2,0,.55), "cm"))
    
    library(patchwork)
    combined_plot <- (this_hist / this_pdp) + 
      plot_layout(heights = c(3, 5), guides = "collect")
    
    # Save with appropriate width/height (increase width if legend is wide)
    # ggsave(file.path(FIGURE_DIR, paste0(cvars[i], "_", fig_suffix,".png")),
    #        combined_plot,
    #        dpi = 300,
    #        width = 8, height = 6)
    all_plts[[i]] <- combined_plot
  }
  
}#



# Extract legend from the first plot (or any plot)
legend_plot <- get_legend(
  all_plts[[1]] +
    theme(legend.position = "right") # ensure the legend is present
)
legend_plot <- as_ggplot(legend_plot)

# Combine the plots and add the legend as the last plot
total_plots <- length(all_plts) + 1
ncol <- 3


all_plts_adj <- lapply(seq_along(all_plts), function(i) {
  p <- all_plts[[i]]
  # Leave legend on only the first plot
  p <- p + theme(legend.position = "none")

  #remove y-axis for non-leftmost columns
  if ((i - 1) %% ncol != 0) {
    p <- p + theme(axis.title.y = element_blank(),
                   axis.text.y = element_blank(),
                   axis.ticks.y = element_blank())
  }
  p
})

# Compose the plot list including the legend
patch_list <- c(all_plts_adj, list(legend_plot))

# Make the layout fill by row, legend will be in the last available grid cell
combined_patchwork <- wrap_plots(patch_list, ncol = ncol) +
  plot_layout(guides = "collect")

# Save all 
ggsave(file.path(FIGURE_DIR, "all_combined_pdps.png"),
       combined_patchwork, width = 13.5, height = 18.5, dpi = 300)



ncol <- 2
which_vars <- rcvars %in%  rcvars[rcvars != "tr_age"][1:4] 
all_plts_adj <- lapply(seq_along(all_plts[which_vars]), function(i) {
  p <- all_plts[which_vars][[i]]
  # Leave legend on only the first plot
  p <- p + theme(legend.position = "none")
  
  #remove y-axis for non-leftmost columns
  if ((i - 1) %% ncol != 0) {
    p <- p + theme(axis.title.y = element_blank(),
                   axis.text.y = element_blank(),
                   axis.ticks.y = element_blank())
  }
  p
})

legend_plot <- get_legend(
  all_plts_adj[[1]] +
    theme(
      legend.position = "right",
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 18)
    )
) |> as_ggplot()

# Make the layout fill by row, legend will be in the last available grid cell
combined_patchwork <- wrap_plots(all_plts_adj, nrow = 2, ncol = 2) +
  plot_layout(guides = "collect")

final_plot <- combined_patchwork | legend_plot +
  plot_layout(widths = c(3, 1))

# Save all 
ggsave(file.path(FIGURE_DIR, "top4_pdps.png"),
       final_plot, width = 13, height = 9.5, dpi = 300)


# Save treatments 

ggsave(file.path(FIGURE_DIR, "tr_age_pdp.png"),
       all_plts[rcvars == "tr_age"][[1]] , width = 8, height = 7, dpi = 300)

# # Final carts 
# data_selected <- mod_df |> select(all_of(c('response', preds))) #|>
#  # mutate(tmmn = tmmn - 273.15)
# 
# # if adding in extreme growth instead of ha's per day
# extreme_growth = TRUE
# if (extreme_growth == TRUE) {
#   data_selected <- data_selected %>%
#     dplyr::mutate(extreme_growth = ifelse(area > 5900, 1, 0)) |>
#     select(-c(area))
# }
# 
# 
# cart_model <- rpart(
#   formula = response ~ .,  # replace with your response variable
#   data = data_selected,
#   method = "class",       
#   control = rpart.control(
#     cp = 0.01, 
#     maxdepth = 7
#   )
# )
# 
# printcp(cart_model)
# plotcp(cart_model)
# 
# 
# 
# # Create a mapping for pretty names
# variable_crosswalk_model <- variable_crosswalk |> filter(Abbreviation %in% colnames(data_selected))
# pretty_names <- setNames(variable_crosswalk_model$Figure_Name, variable_crosswalk_model$Abbreviation)
# 
# replace_var_names <- function(x) {
#   x$frame$var <- ifelse(x$frame$var %in% names(pretty_names),
#                         pretty_names[x$frame$var],
#                         x$frame$var)
#   x
# }
# named_model <- replace_var_names(cart_model)
# # Visualize the tree with custom hex colors
# mod_cart <- rpart.plot(
#   named_model,
#   type   = 5,          # Show decision splits and labels
#   extra  = 106,        # Add probabilities and observation percentages
#   faclen = 0,          # Don't truncate factor names
#   tweak  = 1.5,        # Increase font size
#   box.palette = c("#FFFF00", "#8B0000"),  # Yellow for the first class and dark red for the second
#   main  = ""
# )
# 
# cart_model_py <- as.party(named_model)
# 
# 
# # Plot the classification tree
# ggparty(cart_model_py) +
#   geom_edge() +
#   geom_edge_label() +
#   geom_node_splitvar() +
#   geom_node_plot(
#     gglist = list(
#       geom_bar(aes(x = "", fill = response), 
#                position = position_fill()),  # Proportional response bar
#       xlab(""),                           # Remove x-axis label
#       theme(axis.text.y = element_blank(), # Remove y-axis labels
#             axis.ticks.y = element_blank(), # Remove y-axis ticks
#             panel.background = element_blank(), # Remove background
#             panel.grid.major = element_blank(), # Remove major grid lines
#             panel.grid.minor = element_blank(), # Remove minor grid lines
#             axis.line = element_blank())       # Remove axis lines
#     ),
#     shared_axis_labels = TRUE,                 # Shared axis labels across nodes
#     legend_separator = TRUE                    # Add line between tree and legend
#   ) 
# 
# # Print the custom tree
# print(custom_tree_plot)
# 
# data_selected |>
#   filter(tr_age == 1, srad > 270, rmin < 21) |>
#   ggplot(aes(x=response, fill = response)) +
#   geom_bar(stat = "count", position = "stack")
# 
# data_selected |>
#   filter(tr_age < 13, extreme_growth == 0, MCMT >= -0.16, SCF = ) |>
#   ggplot(aes(x=SCF, fill = response)) +
#   geom_histogram(position = "stack", alpha = 0.6, bins = 30)
# data_selected |>
#   filter(extreme_growth == 1, srad > 270) |>
#   ggplot(aes(x = rmin, fill = response)) +
#   geom_histogram(position = "stack", alpha = 0.6, bins = 30)
# 
# 
# # Predict
# train_pred <- predict(tree_depth_r, train_data, type = "class")
# test_pred  <- predict(tree_depth_r, test_data, type = "class")
# 
# # Accuracy
# mean(train_pred == train_data$response) # Training accuracy
# mean(test_pred == test_data$response)   # Test accuracy
# 
# ggsave2(file.path(FIGURE_DIR, "data_cart.png"),width = 15, height = 19, dpi = 300)
# 
# #-----------------
# 
# # Load required libraries
# library(tidymodels)
# 
# # Split data into training and testing sets
# set.seed(123)  # For reproducibility
# split <- initial_split(data_selected, prop = 0.8)  # 80% train, 20% test
# train_data <- training(split)
# test_data <- testing(split)
# 
# # Create a recipe for preprocessing
# tree_recipe <- recipe(response ~., data = train_data) %>%
#   step_normalize(all_numeric_predictors()) %>%  # Normalize numeric predictors
#   step_dummy(all_nominal_predictors())         # Dummy encode categorical predictors
# 
# # Specify the model
# tree_spec <- decision_tree(
#   mode = "classification",         # Change to "regression" if response is numeric
#   tree_depth = tune(),             # Tuning parameter
#   min_n = tune(),                  # Minimum observations in a node (minsplit)
#   cost_complexity = tune()         # Tuning `cp` parameter
# ) %>%
#   set_engine("rpart")              # Use the rpart engine
# 
# # Combine the recipe and model into a workflow
# tree_workflow <- workflow() %>%
#   add_model(tree_spec) %>%
#   add_recipe(tree_recipe)
# 
# # Create cross-validation folds
# cv_folds <- vfold_cv(train_data, v = 10)  # 10-fold cross-validation
# 
# # Define the hyperparameter grid for tuning
# # Define the parameter objects
# tree_depth_param <- tree_depth(range = c(4, 12))       # Tree depth range
# min_n_param <- min_n(range = c(50, 200))               # Minimum observations per node
# cost_complexity_param <- cost_complexity(c(0.001, 0.05), NULL)  # Fixed values for `cp`
# 
# tree_grid <- grid_regular(
#   tree_depth_param,    
#   min_n_param,
#   cost_complexity_param,
#   levels = 3                       
# )
# 
# # Tune the model using grid search
# set.seed(123)
# tree_tuned <- tune_grid(
#   tree_workflow,
#   resamples = cv_folds,
#   grid = tree_grid,
#   metrics = metric_set(roc_auc, accuracy)  # Choose relevant metrics
# )
# 
# # Fit with best parameters
# best_params <- select_best(tree_tuned, metric = "accuracy")
# final_tree_workflow <- finalize_workflow(tree_workflow, best_params)
# final_fit <- fit(final_tree_workflow, data = train_data)
# 
# # Test dat model accuracy 
# tree_results <- predict(final_fit, new_data = test_data) %>%
#   bind_cols(test_data %>% select(response)) %>%
#   metrics(truth = response, estimate =.pred_class)
# 
# print(tree_results)
# final_rpart <- extract_fit_engine(final_fit)  # Extract the rpart model
# rpart.plot(final_rpart, type = 2, extra = 1)
# 
# install.packages('partykit')
# library(partykit)



# calculate the control mean value 
this_con <- dplyr::filter(this_df, group == 0 & var_val == 21) 

minval <- quantile(this_df$var_val, 0.025)
maxval <- quantile(this_df$var_val, 0.975)
this_df <- filter(this_df, var_val > minval & var_val < maxval)

classes <- unique(this_df$group)[-1]

this_df_plot <-  this_df 
ymin <- 0
ymax <- 0.8
# ymin <- ifelse(min(this_df_plot$lwr)< 0.2,0, 0.2)
# ymax <- ifelse( max(this_df_plot$upr)> 0.8, 1, 0.8)

this_df_plot <- dplyr::filter(this_df_plot, group != 0)

tr_mean_plot <- this_df_plot |> 
  summarise(lwr_mean = mean(lwr),
            tr_mean = mean(median),
            upr_mean = mean(upr))

levs <- unique(this_df_plot$class_name)
levs <- setdiff(levs, c("Untreated", "Rearrangement + Fire")) # Just some order
new_levels <- c(levs, "Rearrangement + Fire", "Untreated")
this_df_plot <- this_df_plot %>%
  mutate(class_name = factor(class_name, levels = new_levels))

(this_pdp <- this_df_plot |> 
    ggplot() +
    annotate(
      "rect",
      xmin = -Inf, xmax = Inf,
      ymin = this_con$lwr, ymax = this_con$upr,
      fill = this_con$treatment_color_hvd, alpha = 0.1
    ) +
    geom_hline(
      yintercept = this_con$median,
      color = this_con$treatment_color_hvd,
      linetype = "solid",
      linewidth = 1.5
    ) +
    geom_ribbon(
      data = tr_mean_plot,
      aes(x = var_val, ymin = lwr_mean, ymax = upr_mean),
      alpha = 0.2
    ) +
    geom_line(
      data = tr_mean_plot,
      aes(x = var_val, y = tr_mean),
      linewidth = 1.5
    ) +
    geom_line(aes(x = var_val, y = median, color = class_name), linewidth = 1.5) +
    labs(x = paste0(first(this_df$Figure_Name), ''),y ="Predicted Probability \nof High Severity",color = NULL, fill = NULL) +
    scale_color_manual(values = setNames(unique(this_df_plot$treatment_color_hvd), 
                                         unique(this_df_plot$class_name))) +
    scale_fill_manual(values = setNames(unique(this_df_plot$treatment_color_hvd), 
                                        unique(this_df_plot$class_name)))+
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    coord_cartesian(ylim = c(ymin, ymax)) +
    theme_hannah(15))
print(this_pdp)


library(dplyr)
library(DALEX)




bg_data <- this_explainer$data

treatment_pdp_list <- list()

for (cls in all_classes) {
  
  if (nrow(bg_data) < 5) next
  
  # --- Sample down if desired (keeps it fast inside the bootstrap) ---
  N_FAST = 100
  if (nrow(bg_data) > N_FAST) {
    bg_subset <- bg_data %>% slice_sample(n = N_FAST)
  }
  
  # --- Set treatment class to the one we're profiling ---
  bg_subset[[tr_var]] <- cls
  bg_subset[[tr_var]] <- factor(bg_subset[[tr_var]], levels = levels(bg_data[[tr_var]]))
  
  
  # --- Assign realistic treatment ---
  if (cls == untreated_label) {
    bg_subset[[age_var]] <- untreated_age
  } else {
    bg_subset[[age_var]] <- bg_data |> filter(.data[[age_var]] >= 0, .data[[age_var]] <= 20) |> pull(!!sym(age_var)) |> sample(size = dim(bg_subset)[[1]])
  }
  
  tmp <- bg_subset
  
  # --- Mean prediction across the (valid) background ---
  avg_pred <- mean(pull(predict(this_fit, tmp, 'prob')[,2]), na.rm = T)
  # or: mean(this_explainer$predict_function(this_fit, bg_subset))
  
  treatment_pdp_list[[cls]] <- tibble(
    var     = tr_var,
    var_val = cls,
    pred    = avg_pred
  )
}

treatment_class_pdp <- bind_rows(treatment_pdp_list)

# --- Quick plot ---
library(ggplot2)

ggplot(treatment_class_pdp,
       aes(x = reorder(var_val, pred), y = pred)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(
    x = "Treatment Class",
    y = "Predicted P(High Severity)",
    title = "Conditional PDP: Treatment Class\n(averaged over realistic treatment ages only)"
  ) +
  theme_minimal()

# C) Variable importance and Partial dependence plots with CI's-----------------

# Base explainer using the first fit
fit1 <- readRDS(
  sprintf("%sbootstraps%d/boot_fits/fit_%04d.rds",
          paste0(df_dir_output, '/'), n_boots, 1)
)$fit

base_explainer <- explain_tidymodels(
  model = fit1,
  data  = mod_df %>% dplyr::select(-all_of(response)),
  y     = if (response_type == "binary") {
    as.integer(mod_df[[response]])
  } else {
    mod_df[[response]]
  },
  verbose = FALSE
)

base_data <- base_explainer$data


# Define the treatment age onstraints 
untreated_label <- 0 
untreated_age <- 21
treated_age_range <- 2:20

# Determine correct variable names
tr_var <- if (wto_outcome == "infrastructure") "class_cod.uw" else "class_cod"
age_var <- "tr_age"  # adjust if your column name differs

# Identify treated classes
all_classes <- unique(base_data[[tr_var]])
untreated_label <- 0  # adjust to match your factor level


process_bootstrap <- function(i, base_explainer) {
  tryCatch(
    {
      # --- Load fit ---
      this_fit <- readRDS(
        sprintf("%sbootstraps%d/boot_fits/fit_%04d.rds",
                paste0(df_dir_output, '/'), n_boots, i)
      )$fit
      
      
      # --- Reuse explainer ---
      this_explainer <- base_explainer
      this_explainer$model <- this_fit
      this_data <- base_data
      
      N_FAST <- 500
      
      tr_var <- if (wto_outcome == "infrastructure") "class_cod.uw" else "class_cod"
      
      # --- Treatment PDP (manual ICE-style, vectorized over classes) ---
      tr_pdp <- dplyr::bind_rows(lapply(all_classes, function(cls) {
        
        df_subset <- if (nrow(this_data) > N_FAST) {
          this_data |> dplyr::slice_sample(n = N_FAST)
        } else {
          this_data
        }
        
        # --- Set treatment class to the one we're profiling ---
        df_subset[[tr_var]] <- factor(cls, levels = levels(this_data[[tr_var]]))
        
        # --- Assign realistic age given treatment class ---
        if (cls == untreated_label) {
          df_subset[[age_var]] <- untreated_age
        } 
        if (cls != untreated_label) {
          tr_age_sample <- this_data |>
            dplyr::filter(
              .data[[age_var]] >= min(treated_age_range),
              .data[[age_var]] <= max(treated_age_range)
            ) |>
            dplyr::pull(!!sym(age_var)) |>
            sample(size = nrow(df_subset))
          
          df_subset[[age_var]] <- tr_age_sample
        }
        
        # --- Mean prediction across the background ---
        avg_pred <- mean(
          dplyr::pull(predict(this_fit, df_subset, 'prob')[, 2]),
          na.rm = TRUE
        )
        
        tibble(var = tr_var, var_val = cls, pred = avg_pred)
      }))
      
     
      # --- Success return ---
      list(
        ok     = TRUE,
        boot   = i,
        tr_pdp = tr_pdp,
        error  = NULL
      )
    },
    
    error = function(e) {
      message(sprintf("❌ Bootstrap %d failed: %s", i, e$message))
      list(
        ok     = FALSE,
        boot   = i,
        tr_pdp = NULL,
        error  = e$message
      )
    }
  )
}
system.time({results2 <- lapply(
  seq_len(n_boots), process_bootstrap, base_explainer = base_explainer)})

tr_pdps <- lapply(results2, function(x) x$tr_pdp)

tr_df_bind <-bind_rows(tr_pdps)

tr_pdp_ci <- tr_df_bind |>
  group_by(var_val) %>%
  dplyr::summarise(lwr = quantile(pred, 0.025, na.rm = T),
                   median = mean(pred, na.rm = T),
                   upr = quantile(pred, 0.975, na.rm = T))
tr_pdp_ci$var_val <- as.factor(tr_pdp_ci$var_val)

class_code$class_cod <- as.factor(class_code$class_cod)
tr_pdp_ci <- left_join(tr_pdp_ci, class_code, join_by("var_val" == "class_cod"))

tr_pdp_ci <- tr_pdp_ci |> arrange(intensity_rank) |> mutate(class_name = factor(class_name, levels = class_name))

(mean_preds <- tr_pdp_ci |> 
    ggplot() +
    geom_hline(yintercept = 0.5, linewidth = 1, color = 'grey', linetype='dashed') +
    geom_point(aes(x = class_name, y = median, color = class_name),
               size = 4, shape = 16, position=position_dodge(width=0.5)) +
    geom_errorbar(aes(x = class_name, ymin = lwr, ymax = upr, color = class_name),
                  linewidth = 1.5, width=0, position=position_dodge(width=0.5)) +
    labs(y = "Predicted Probability \nof High Severity", x = "", color = '', fill = '') +
    coord_cartesian(ylim = c(0,1)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0.2, 0.7)) +
    scale_color_manual(values = c(tr_pdp_ci$treatment_color_hvd), 
                       labels = c(tr_pdp_ci$class_name)) +
    theme_hannah(16) +
    coord_flip())
