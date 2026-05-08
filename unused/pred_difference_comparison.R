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
              'rsample')


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
region <- 'west' # set to blank if it is the entire dataset
fig_suffix <- region

# Name of response
response = "response"

# Treatment (or any categorical you want to split the data by)
trt_var <- 'class_cod'

# number of bootstraps, probably should be 1000 for final runs
n_boots <- 1000

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

#class_code <- read.csv(r"(~/Library/CloudStorage/Box-Box/External Wildfire Treatment wto_outcomes/data/global_data/treatments/treatment_categories/treatment_code_naming.csv)")

#fire_list <- read.csv(r"(~/Library/CloudStorage/Box-Box/External Wildfire Treatment wto_outcomes/data/interim_data/fire_lists/burn_severity_fire_list.csv)")

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


library(tidymodels)
library(iml)
library(dplyr)
library(ggplot2)


train_df = mod_df

# ── 1. Define a prediction wrapper ──────────────────────────────────────────
# iml needs a function that returns a numeric vector or data frame of probs.
# Adjust `.pred_<classname>` to match your actual positive class label.

pred_fun <- function(model, newdata) {
  predict(model, new_data = newdata, type = "prob")[,2]
}

# ── 2. Create the iml Predictor object ──────────────────────────────────────
# `data`  = your feature data (no outcome column)
# `y`     = your outcome vector

predictor <- Predictor$new(
  model            = fit1,
  data             = train_df %>% select(-!!sym(response)),   # replace `outcome` with your target column name
  y                = train_df[[response]],
  predict.function = pred_fun
)

# ── 3. Order Categorical Variables ──────────────────────────────────────
# ── Option A: Order by predicted mean probability ───────────────────────────
# Most theoretically sound for ALE

pred_probs <- predict(fit1, new_data = train_df, type = "prob")

# Adjust `.pred_<your_class>` to whichever class probability you're modeling
mean_pred_by_class <- train_df %>%
  mutate(pred_prob = pull(pred_probs[,2])) %>%  
  group_by(class_cod) %>%
  summarise(mean_pred = mean(pred_prob)) %>%
  arrange(mean_pred)

ordered_levels_pred <- mean_pred_by_class$class_cod
cat("Order by predicted mean:\n"); print(ordered_levels_pred)

# ── Option B: Order by user set "intesity  ───────────────────────────
# Most theoretically sound for ALE

ordered_levels_pred <- class_code$class_cod
cat("Order by predicted mean:\n"); print(ordered_levels_pred)

# ── Option C: Order by numberic value  ───────────────────────────
# Not good, just to test

ordered_levels_pred <- sort(as.character(class_code$class_cod))
cat("Order by predicted mean:\n"); print(ordered_levels_pred)

# ── 3. Calculate ALE Plots ──────────────────────────────────────
# Apply your chosen ordering as factor levels BEFORE passing to iml
# Here we use Option A (predicted mean) — recommended for ALE

train_df_ordered <- train_df %>%
  mutate(class_cod = factor(class_cod, levels = ordered_levels_pred))

# Rebuild Predictor with the re-leveled data
predictor_ordered <- Predictor$new(
  model            = fit1,
  data             = train_df_ordered %>% select(-!!sym(response)),
  y                = train_df_ordered[[response]],
  predict.function = pred_fun
)

# Compute ALE for class_cod
ale_class_cod <- FeatureEffect$new(
  predictor = predictor_ordered,
  feature   = "class_cod",
  method    = "ale"
)
#ale_class_cod$plot()

# iml can still return results in a different order, so pin it in the plot
ale_data <- ale_class_cod$results %>%
  mutate(class_cod = factor(as.character(class_cod), levels = ordered_levels_pred))
ale_data <- left_join(ale_data, class_code)
ale_data <- ale_data |> arrange(intensity_rank) |> mutate(class_name = factor(class_name, levels = class_name))

ale_data |> 
ggplot() +
  geom_hline(yintercept = 0, linewidth = 1, color = 'grey', linetype='dashed') +
  geom_point(aes(x = class_name, y = mean, color = class_name),
             size = 4, shape = 16, position=position_dodge(width=0.5)) +
  geom_errorbar(aes(x = class_name, ymin = lwr, ymax = upr, color = class_name),
                linewidth = 1.5, width=0, position=position_dodge(width=0.5)) +
  labs(y = "Mean Difference in Predicted Probability \nof High Severity (P(treated) - P(untreated))", x = "", color = '', fill = '') +
  #coord_cartesian(ylim = c(0,1)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(-0.3,0.25)) +
  scale_color_manual(values = c(pred_dif_ci$treatment_color_hvd), 
                     labels = c(pred_dif_ci$class_name)) +
  theme_hannah(16) +
  coord_flip())

post_model_preds <- function(i){
  tryCatch(
      {
    # --- Load fit ---
    readRDS(
      sprintf("%sbootstraps%d/boot_fits/fit_%04d.rds",
              paste0(df_dir_output, '/'), n_boots, i)
    )$fit |>
      predict(mod_df, type = 'prob') |>
        bind_cols(mod_df) |>
        mutate(matched_id = df$matched_id) |>
          group_by(matched_id) |> 
          dplyr::summarise(pred_diff = .pred_1[class_cod != 0] - .pred_1[class_cod == 0],
              class_cod = class_cod[class_cod != 0]) |>
      group_by(class_cod) |>
      dplyr::summarise(mean_pred_diff = mean(pred_diff))
      },

    error = function(e) {
      message(sprintf("Bootstrap %d failed: %s", i, e$message))
      
      return(NULL)
    }
  )
}
  

results <- lapply(seq_len(n_boots)[1:100], post_model_preds)

pred_dif_ci <- bind_rows(results) |>
  group_by(class_cod) %>%
    dplyr::summarise(lwr = quantile(mean_pred_diff, 0.0275),
                   mean = mean(mean_pred_diff),
                   upr = quantile(mean_pred_diff, 0.975)) |>
  mutate(class_cod = factor(as.character(class_cod), levels = class_code$class_cod))
pred_dif_ci <- left_join(pred_dif_ci, class_code)

pred_dif_ci <- pred_dif_ci |> arrange(intensity_rank) |> mutate(class_name = factor(class_name, levels = class_name))
pred_dif_ci <- pred_dif_ci |> mutate(region = region)

# code if I manually changed region to bind both regions together
pred_dif_ci_regions <- bind_rows(pred_dif_ci_old, pred_dif_ci)
(pred_dif <- pred_dif_ci_regions |> 
    ggplot() +
    geom_hline(yintercept = 0, linewidth = 1, color = 'grey', linetype='dashed') +
    geom_point(aes(x = class_name, y = mean, color = class_name, shape = region),
               size = 4, position=position_dodge(width=0.5)) +
    geom_errorbar(aes(x = class_name, ymin = lwr, ymax = upr, color = class_name),
                  linewidth = 1.5, width=0, position=position_dodge(width=0.5)) +
    labs(y = "Mean Difference in Predicted Probability \nof High Severity (P(treated) - P(untreated))", x = "", color = '', fill = '') +
    #coord_cartesian(ylim = c(0,1)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(-0.7,0.05)) +
    scale_color_manual(values = c(pred_dif_ci$treatment_color_hvd), 
                       labels = c(pred_dif_ci$class_name)) +
    theme_hannah(16) +
    coord_flip())

(pred_dif <- pred_dif_ci_regions |> 
    ggplot() +
    geom_hline(yintercept = 0, linewidth = 1, color = 'grey', linetype = 'dashed') +
    # ✅ errorbar drawn FIRST so points render on top
    # geom_errorbar(aes(x = class_name, ymin = lwr, ymax = upr,
    #                   color = class_name, group = region),
    #               linewidth = 0.5, width = 0, position = position_dodge(width = 0.5)) +
    # ✅ points drawn AFTER so they sit on top of the bars
    geom_point(aes(x = class_name, y = mean, color = class_name, shape = region),
               size = 4,position = position_dodge(width = 0.5)) +
    labs(y = "Mean Difference in Predicted Probability \nof High Severity (P(treated) - P(untreated))",
         x = "", color = '', fill = '', shape = '') +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(-0.7, 0.05)) +
    scale_color_manual(values = c(pred_dif_ci$treatment_color_hvd),
                       labels = c(pred_dif_ci$class_name)) +
    # ✅ solid circle (16) for first region, open triangle (2) for second
    scale_shape_manual(values = c(16, 2)) +
    theme_hannah(16) +
    coord_flip())






