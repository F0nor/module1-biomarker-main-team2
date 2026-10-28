setwd("C:/Users/Felix/Documents/GitHub/module1-biomarker-team2/scripts")
source("preprocessing.R")
load("../data/biomarker-clean.RData")



biomarker_raw <- read_csv('../data/biomarker-raw.csv', 
                                skip = 2,
                                col_select = -2L,
                                col_names = c('group', 
                                              'empty',
                                              pull(var_names, abbreviation),
                                              'ados'),
                                na = c('-', '')) %>%
  # filter(!is.na(group)) %>%
  # log transform, center and scale, and trim
  # mutate(across(.cols = -c(group, ados), 
  # ~ trim(scale(log10(.x))[, 1], .at = 3))) %>%
  # reorder columns
  select(group, ados, everything())


hist_PreTransform <- biomarker_raw |>
  pivot_longer(
    cols = -c(group, ados),
    names_to = "Proteins",
    values_to = "Values"
  ) |>
  filter(Proteins %in% sample(unique(Proteins), 9)) |>
  ggplot(aes(x = Values)) +
  geom_histogram(stat = "bin", bins = 60) +
  labs(
    title = "Histogram of Sampled Pre-Transformed Protein Values",
    x = "Protein Values",
    y = "Count"
  ) 

ggsave(filename = "../results/figures/Hist-Pretransf.jpg", hist_PreTransform)

biomarker_clean |> 
  pivot_longer(
    cols = -c(group, -ados),
    names_to = "Proteins",
    values_to = "values"
  ) |>
  group_by(Proteins) |>
  summarise(
    min = min(values),
    max = max(values),
    mean = mean(values), 
    median =  median(values)
  ) 

hist_LogTransform <- biomarker_clean |> 
  pivot_longer(
    cols = -c(group, -ados),
    names_to = "Proteins",
    values_to = "values"
  ) |>
  filter(Proteins %in% sample(unique(Proteins), 9)) |>
  ggplot(aes(x = values)) +
  geom_histogram( bins = 60) +
  labs( 
    title = "Histogram of Sampled Log-Transformed Protein Values",
    x = "Protein values",
    y = "Count")

ggsave(filename = "../results/figures/Hist-Logtransf.jpg", hist_LogTransform)
