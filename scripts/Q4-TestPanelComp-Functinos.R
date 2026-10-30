library(tidyverse)
library(rsample)
library(modelr)
library(yardstick)
library(randomForest)
library(infer)
load("../data/biomarker-clean.RData")
source("Q3-TestPanelEval-Functions.R")
source("inclass-analysis.R")

# For training
set.seed(101422)

partition <- biomarker_clean %>%
  initial_split(prop = 0.8)

biomarker_train <- training(partition)
biomarker_test <- testing(partition)



get_full_panels <- function(n){
  ttest_proteins <- ttests_out |> slice_min(p.adj, n=n) |> pull(protein)
  rf_proteins <- rf_out$importance |> as_tibble() |>
    mutate(protein = rownames(rf_out$importance)) |> slice_max(MeanDecreaseGini, n = n) |> pull(protein)
  
  fuzzy <- fuzzy_panel(n, ttests_out, rf_out)
  hard <- intersect(ttest_proteins, rf_proteins)
  list(hard = hard, fuzzy = fuzzy)
}




list_ns <- c(10, 15, 20, 25, 30, 35)

Q3_results <- map(list_ns, function(n){
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

all_results <- bind_rows(
  results |> mutate(source = "full data"),
  Q3_results |> mutate(source = "split data")
)

# results ranked
all_results |>
  filter(`.estimate` >= original_accuracy & `.metric` == "accuracy") |>
  arrange(desc(.estimate), length)

# split data is more honest

all_results |>
  filter(source == "split data" & `.estimate` & .metric %in% c("accuracy", "roc_auc")) |>
  select(type, n, length, .metric, .estimate) |>
  pivot_wider(names_from = .metric, values_from = .estimate) |>
  arrange(desc(accuracy), desc(roc_auc), length)

split.simp_betterAC_panel_proteins <- all_results |>
  filter(source == "split data" & `.estimate` >= original_accuracy & `.metric` == "accuracy") |>
  arrange(desc(.estimate), length) |>
  head(1) |>
  pull(panel) |> pluck(1)

split.simp_betterAC_panel <- all_results |>
  filter(source == "split data" & `.estimate` >= original_accuracy & `.metric` %in% c("accuracy") ) |>
  arrange(desc(.estimate), length) |>
  head(1)


full.simp_betterAC_panel_proteins <- all_results |>
  filter(source == "full data" & `.estimate` >= original_accuracy & `.metric` == "accuracy") |>
  arrange(desc(.estimate), length) |>
  head(1) |>
  pull(panel) |> pluck(1)

full.simp_betterAC_panel <- all_results |>
  filter(source == "full data" & `.estimate` >= original_accuracy & `.metric` %in% c("accuracy") ) |>
  arrange(desc(.estimate), length) |>
  head(1)

# Save data, might want to save other dfs

save(list = "all_results", file = "../data/all_results.RData" )

