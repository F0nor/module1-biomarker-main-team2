library(tidyverse)
library(rsample)
library(modelr)
library(yardstick)
library(randomForest)
library(infer)

# Navigate to the folder on your computer containing this repo in the "Files" window
# (it should say something like "Home > Desktop > (whatever folder you have the repo in) > module1-biomarker-team2 > data")
#. Then, click "More" -> "Set as Working Directory", so you can run the lines below:

load("../data/biomarker-clean.RData") # With trimming
load("../data/biomarker-clean_no_trim.RData") # Without trimming

source("../scripts/Q3-TestPanelEval-Functions.R")
source("../scripts/inclass-analysis.R")

# evaluate errors on test set
class_metrics <- metric_set(sensitivity, 
                            specificity, 
                            accuracy,
                            roc_auc)



original_metrics <- testing(biomarker_split) %>%
  add_predictions(fit, type = 'response') %>% 
  mutate(estimate = factor(pred >0.5, levels = c(FALSE, TRUE)),
         truth = factor(class, levels = c(FALSE, TRUE))) %>%
  class_metrics(estimate = estimate,
                truth = truth, pred,
                event_level = 'second')

# For training
set.seed(101422)

## Without Trimming (Hannah)

# Training partition method
partition2 <- biomarker_clean_out %>%
  initial_split(prop = 0.8)

biomarker_train <- training(partition2)
biomarker_test <- testing(partition2)



Q3_results_no_trim <- map(list_ns, function(n){
  # get the panel of n proteins 
  p <- get_train_panels(n)
  
  # iterate through the panels by type
  map(names(p), function(type){
    metrics <- evaluate(p[[type]])
    proteins <- p[[type]]
    # if empty then return NULL
    if (is.null(metrics)) return (NULL)
    
    # How many proteins captured by fuzzy or hard intersect
    length <- length(p[[type]])
    
    # Add descriptive columns to differentiate metrics
    metrics |>
      mutate(type = type, n = n, length = length, panel = list(proteins)) 
  }) |> list_rbind() 
}) |> list_rbind() # combine each list by rows into one data frame


results <- map(list_ns, function(n){
  # get the panel of n proteins 
  p <- get_full_panels(n)
  
  # iterate through the panels by type
  map(names(p), function(type){
    metrics <- evaluate(p[[type]])
    proteins <- p[[type]]
    # if empty then return NULL
    if (is.null(metrics)) return (NULL)
    
    # How many proteins captured by fuzzy or hard intersect
    length <- length(p[[type]])
    
    # Add descriptive columns to differentiate metrics
    metrics |>
      mutate(type = type, n = n, length = length, panel = list(proteins)) 
  }) |> list_rbind() 
}) |> list_rbind() # combine each list by rows into one data frame


original_accuracy <- original_metrics[original_metrics$.metric == "accuracy",]$.estimate

# combine into one data frame

all_results_no_trim <- bind_rows(
  results |> mutate(source = "full data"),
  Q3_results_no_trim |> mutate(source = "split data")
)

# results ranked
all_results_no_trim |>
  filter(`.estimate` >= original_accuracy & `.metric` == "accuracy") |>
  arrange(desc(.estimate), length)

# split data is more honest

all_results_no_trim |>
  filter(source == "split data" & `.estimate` & .metric %in% c("accuracy", "roc_auc")) |>
  select(type, n, length, .metric, .estimate) |>
  pivot_wider(names_from = .metric, values_from = .estimate) |>
  arrange(desc(accuracy), desc(roc_auc), length)

split.simp_betterAC_panel_proteins <- all_results_no_trim |>
  filter(source == "split data" & `.estimate` >= original_accuracy & `.metric` == "accuracy") |>
  arrange(desc(.estimate), length) |>
  head(1) |>
  pull(panel) |> pluck(1)

split.simp_betterAC_panel <- all_results_no_trim |>
  filter(source == "split data" & `.estimate` >= original_accuracy & `.metric` %in% c("accuracy") ) |>
  arrange(desc(.estimate), length) |>
  head(1)


full.simp_betterAC_panel_proteins <- all_results_no_trim |>
  filter(source == "full data" & `.estimate` >= original_accuracy & `.metric` == "accuracy") |>
  arrange(desc(.estimate), length) |>
  head(1) |>
  pull(panel) |> pluck(1)

full.simp_betterAC_panel <- all_results_no_trim |>
  filter(source == "full data" & `.estimate` >= original_accuracy & `.metric` %in% c("accuracy") ) |>
  arrange(desc(.estimate), length) |>
  head(1)

# Save data, might want to save other dfs

save(list = "all_results_no_trim", file = "../data/all_results_no_trim.RData" )
