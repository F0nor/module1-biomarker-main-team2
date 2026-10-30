library(tidyverse)
library(rsample)
library(modelr)
library(yardstick)
library(randomForest)
library(infer)

# Navigate to the folder on your computer containing this repo in the "Files" window
# (it should say something like "Home > Desktop > (whatever folder you have the repo in) > module1-biomarker-team2 > data")
#. Then, click "More" -> "Set as Working Directory", so you can run the lines below:

load("../data/biomarker-clean.RData")
set.seed(101422)

partition <- biomarker_clean %>%
  initial_split(prop = 0.8)

biomarker_train <- training(partition)
biomarker_test <- testing(partition)

# function to compute tests
test_fn <- function(.df){
  t_test(.df, 
         formula = level ~ group,
         order = c('ASD', 'TD'),
         alternative = 'two-sided',
         var.equal = F)
}


train.ttests_out <- biomarker_train %>%
  # drop ADOS score
  select(-ados) %>%
  # arrange in long format
  pivot_longer(-group, 
               names_to = 'protein', 
               values_to = 'level') %>%
  # nest by protein
  nest(data = c(level, group)) %>% 
  # compute t tests
  mutate(ttest = map(data, test_fn)) %>%
  unnest(ttest) %>%
  # sort by p-value
  arrange(p_value) %>%
  # multiple testing correction
  mutate(m = n(),
         hm = log(m) + 1/(2*m) - digamma(1),
         rank = row_number(),
         p.adj = m*hm*p_value/rank)

## RANDOM FOREST
##################

# store predictors and response separately
train.predictors <- biomarker_train %>%
  select(-c(group, ados))

train.response <- biomarker_train %>% pull(group) %>% factor()

# fit RF
set.seed(101422)
train.rf_out <- randomForest(x = train.predictors, 
                       y = train.response, 
                       ntree = 1000, 
                       importance = T)


fuzzy_panel <- function(n, ttests, rf, thresh = 0.25){
  mem <- function(rank) pmax(0, 1 - (rank - 1) / (2 * n))
  
  t_rank <- ttests |>
    arrange(p.adj) |>
    transmute(protein, m_t = mem(row_number()))
  
  r_rank <- rf$importance |>
    as_tibble(rownames = "protein") |>
    arrange(desc(MeanDecreaseGini)) |>
    transmute(protein, m_r = mem(row_number()))
  
  inner_join(t_rank, r_rank, by = "protein") |>
    mutate(m = pmin(m_t, m_r)) |>
    filter(m >= thresh) |>
    pull(protein)
}


get_train_panels <- function(n){
  ttest_proteins <- train.ttests_out |> slice_min(p.adj, n=n) |> pull(protein)
  rf_proteins <- train.rf_out$importance |> as_tibble() |>
    mutate(protein = rownames(train.rf_out$importance)) |> slice_max(MeanDecreaseGini, n = n) |> pull(protein)
  
  
  fuzzy <- fuzzy_panel(n, train.ttests_out, train.rf_out)
  hard <- intersect(ttest_proteins, rf_proteins)
  list(hard = hard, fuzzy = fuzzy)
}



evaluate <- function(panel){
  
  if(length(panel) == 0) return(NULL) # Just in case panel has no proteins
  
  fit <- glm(class ~., data = biomarker_train |> mutate(class = (group == "ASD")) |> select(class, all_of(panel)), family = 'binomial')
  
  class_metrics <- metric_set(sensitivity, 
                              specificity, 
                              accuracy,
                              roc_auc)
  
  biomarker_test |> 
    add_predictions(fit, type = 'response') |>
    mutate(estimate = factor(pred > 0.5, levels = c(FALSE, TRUE)),
           class = (group == "ASD"),
           truth = factor(class, levels = c(FALSE, TRUE))) |>
    class_metrics(estimate = estimate, truth = truth, pred, event_level = 'second')
}


Q3_results <- map(c(5,10, 15, 20, 25, 30), function(n){
  # get the panel of n proteins 
  p <- get_train_panels(n)
  
  # iterate through the panels by type
  map(names(p), function(type){
    metrics <- evaluate(p[[type]])
    # if empty then return NULL
    if (is.null(metrics)) return (NULL)
    
    # How many proteins captured by fuzzy or hard intersect
    length <- length(p[[type]])
    
    # Add descriptive columns to differentiate metrics
    metrics |>
      mutate(type = type, n = n, length = length)
  }) |> list_rbind() 
}) |> list_rbind() # combine each list by rows into one data frame


save(list = 'Q3_results', 
     file = '../data/Q3_results.RData')
