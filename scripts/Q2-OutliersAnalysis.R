library(tidyverse)
library(ggplot2)

# Navigate to the folder on your computer containing this repo in the "Files" window
# (it would say something like "Home > Desktop > (whatever folder you have the repo in) > module1-biomarker-team2")
#. Then, click "More" -> "Set as Working Directory", so you can run the lines below:

var_names <- read_csv('../data/biomarker-raw.csv', 
                      col_names = F, 
                      n_max = 2, 
                      col_select = -(1:2)) %>%
  t() %>%
  as_tibble() %>%
  rename(name = V1, 
         abbreviation = V2) %>%
  na.omit()

# function for trimming outliers (good idea??)
# trim <- function(x, .at){
#   x[abs(x) > .at] <- sign(x[abs(x) > .at])*.at
#   return(x)
# }

# read in data and remove trimming
biomarker_clean_out <- read_csv('../data/biomarker-raw.csv', 
                                skip = 2,
                                col_select = -2L,
                                col_names = c('group', 
                                              'empty',
                                              pull(var_names, abbreviation),
                                              'ados'),
                                na = c('-', '')) %>%
  filter(!is.na(group)) %>%
  # log transform, center and scale, and trim
  mutate(across(.cols = -c(group, ados),
                ~ scale(log10(.x)))) %>%
  # reorder columns
  select(group, ados, everything())

# Save the data
save(list = 'biomarker_clean_out', 
     file = '../data/biomarker-clean_no_trim.RData')

# Find outliers in the untrimmed data
outliers <- biomarker_clean_out |>
  mutate(subject = row_number()) |>
  pivot_longer( # move proteins and values into one column respectively
    cols = -c(group, ados, subject),
    names_to = "Proteins",
    values_to = "Values"
  ) |>
  filter(abs(Values) > 3) # Same criteria as original (|z| > 3)

# Count the number of outliers each subject has
outliers_cnt <- biomarker_clean_out |>
  mutate(subject = row_number()) |>
  pivot_longer(
    cols = -c(group, ados, subject),
    names_to = "Proteins",
    values_to = "Values"
  ) |>
  group_by(subject, group) |>
  summarise(n_out = sum(abs(Values) > 3)) |>
  arrange(desc(n_out)) # Sort in descending order

# Graph the distribution of num of outliers each subject has by group
OutDist_byGroup <- outliers_cnt |>
  ggplot(aes(x = n_out, fill = group)) +
  geom_histogram(bins = 60) +
  facet_wrap(~ group) +
  scale_fill_manual(values = c("ASD" = "pink", "TD" = "orange")) +
  labs(
    title = "Distribution of Outliers in Each Subject",
    x = "Number of Outliers (|z| > 3)",
    y = "Count of Subjects"
  )

# Save image
ggsave(filename = "../results/figures/OutDist_byGroup.jpg", OutDist_byGroup)


outliers_summary <- outliers_cnt |>
  group_by(group) |>
  summarise(min = min(n_out),
            max = max(n_out),
            mean = mean(n_out),
            median = median(n_out))
# 152 subjects had outliers
# Highest amount of outliers was 157 for TD and 126 for ASD
# mean outliers was 13 for ASD and 17 for TD
# Both with median of 9 outliers

# TD: 1372 outliers
# ASD: 1007 outliers

# Outliers do not seem to be significantly more prevalant in one group than another. TD marginally has more outliers than ASD.

outliers_wilcox <- wilcox.test(n_out ~ group, data = outliers_cnt)