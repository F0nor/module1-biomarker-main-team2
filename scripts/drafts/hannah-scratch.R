library(tidyverse)
library(rsample)
library(modelr)
library(yardstick)
library(randomForest)
library(infer)

## Combining Trimmed and Untrimmed Results

load("../data/all_results.RData") # With trimming
load("../data/all_results_no_trim.RData") # Without trimming

trimmed_modified <- all_results |> mutate(partition = "Trimmed")
untrimmed_modified <- all_results_no_trim |> mutate(partition = "Untrimmed")

all_results_combined <- bind_rows(trimmed_modified, untrimmed_modified)
  
# Saving the dataset
save(list = "all_results_combined", file = "../data/all_results_combined.RData" )


## Top 10 within each partition

library(gt)   # install.packages(c("gt", "webshot2")) if needed

# Top 10 untrimmed (only split data)
untrimmed_top10 <- all_results_no_trim |>
  filter(source == "split data", .metric == "accuracy") |>
  arrange(desc(.estimate), length) |>
  slice_head(n = 10) |>
  select(type, n, length, accuracy = .estimate) |>
  mutate(accuracy = round(accuracy, 3))

# Seeing all results
all_results_no_trim |>
  filter(.metric == "accuracy") |>
  arrange(desc(.estimate), length)

untrimmed_top10

untrimmed_top10 |> gt() |> gtsave("../results/figures/untrimmed_top10.png")

# Top 10 trimmed
trimmed_top10 <- all_results |>
  filter(source == "split data", .metric == "accuracy") |>
  arrange(desc(.estimate), length) |>
  slice_head(n = 10) |>
  select(type, n, length, accuracy = .estimate) |>
  mutate(accuracy = round(accuracy, 3))

trimmed_top10

trimmed_top10 |> gt() |> gtsave("../results/figures/trimmed_top10.png")



