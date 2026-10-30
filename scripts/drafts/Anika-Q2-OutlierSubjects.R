# This script builds on Q2-OutliersAnalysis.R
# adds a baseline for "how many outliers is normal",
# a rule for flagging outlier subjects, and
# a size-adjusted group comparison

library(tidyverse)
library(infer)
library(randomForest)
library(rsample)
library(modelr)
library(yardstick)

load(file.path(root, "data/biomarker-clean_no_trim.RData"))  # biomarker_clean_out (untrimmed)
load(file.path(root, "data/biomarker-clean.RData"))          # biomarker_clean (trimmed, for comparison)

# Works whether your working directory is the project root, scripts/, or scripts/drafts/
root <- c(".", "..", "../..")[file.exists(file.path(c(".", "..", "../.."), "data/biomarker-raw.csv"))][1]

# scale() without [, 1] stored each protein as a 1-col matrix; flatten to plain numeric
# and add a subject ID (row order = row order in biomarker-raw.csv)
untrimmed <- biomarker_clean_out %>%
  mutate(across(-c(group, ados), ~ as.numeric(.x))) %>%
  mutate(subject_id = row_number(), .before = 1)

n_proteins <- ncol(untrimmed) - 3

# outlying values per subject

outlier_long <- untrimmed %>%
  pivot_longer(-c(subject_id, group, ados),
               names_to = "protein",
               values_to = "z") %>%
  mutate(outlier = abs(z) > 3)

subject_counts <- outlier_long %>%
  group_by(subject_id, group, ados) %>%
  summarise(n_out  = sum(outlier),
            n_high = sum(z > 3),     # unusually high protein levels
            n_low  = sum(z < -3),    # unusually low protein levels
            .groups = "drop") %>%
  arrange(desc(n_out))

# if log-levels were exactly normal, 
# how many |z| > 3 per subject would we expect?
expected_per_subject <- n_proteins * 2 * pnorm(-3)
expected_per_subject   # ~3.6

summary(subject_counts$n_out)
# median is 8.5, so even "typical" subjects 
# exceed the normal baseline -> heavy tails overall.
# The question is which subjects are extreme 
# relative to other subjects

# flag outlier subjects
q <- quantile(subject_counts$n_out, c(0.25, 0.75))
cutoff <- q[2] + 1.5 * diff(q)
cutoff   # ~29

subject_counts <- subject_counts %>%
  mutate(outlier_subject = n_out > cutoff)

outlier_subjects <- subject_counts %>% 
  filter(outlier_subject)
outlier_subjects
# 13 subjects flagged. Six of them (rows 154, 108, 9, 121, 52, 77) have 114-157
# outlying values each, a clear break from everyone else (next highest is 77).

# are their outliers mostly high, mostly low, or mixed?
outlier_subjects %>% 
  select(subject_id, group, n_out, n_high, n_low)
# Several are almost entirely one-directional (subject 9: 126 high / 0 low;
# subject 77: 114 / 0; subject 147: 77 / 0). That looks like a whole-sample shift
# (e.g., sample handling or concentration) rather than a few odd proteins.

# Rank plot: each point is a subject
subject_plot <- subject_counts %>%
  mutate(rank = row_number()) %>%
  ggplot(aes(x = rank, y = n_out, color = group)) +
  geom_point() +
  geom_hline(yintercept = cutoff, linetype = "dashed") +
  geom_hline(yintercept = expected_per_subject, linetype = "dotted") +
  geom_text(data = ~ filter(.x, outlier_subject),
            aes(label = subject_id), vjust = -0.7, size = 3, show.legend = FALSE) +
  scale_color_manual(values = c("ASD" = "pink", "TD" = "orange")) +
  labs(title = "Outlying protein values per subject",
       subtitle = "Dashed: IQR cutoff for outlier subjects. Dotted: expected count if data were normal",
       x = "Subject (sorted by number of outliers)",
       y = "Number of values with |z| > 3",
       color = "Group")
subject_plot

ggsave(file.path(root, "results/figures/Q2-OutlierSubjects.jpg"), subject_plot, width = 8, height = 5)

# are outlier subjects more common in one group?

# Group sizes differ (76 ASD vs 78 TD), so compare per-subject rates, not raw totals
group_summary <- subject_counts %>%
  group_by(group) %>%
  summarise(n_subjects = n(),
            total_outliers = sum(n_out),
            mean_per_subject = mean(n_out),
            median_per_subj = median(n_out),
            n_outlier_subj = sum(outlier_subject),
            prop_outlier_subj = mean(outlier_subject))
group_summary

# Proportion of outlier subjects by group
outlier_table <- table(subject_counts$group, subject_counts$outlier_subject)
outlier_table
fisher.test(outlier_table)   
# p is approx. 0.25

# Distribution of counts by group
wilcox.test(n_out ~ group, data = subject_counts)   
# p is approx. 0.82

# Within ASD: do outlier counts track symptom severity (ADOS)?
subject_counts %>%
  filter(group == "ASD") %>%
  with(cor.test(n_out, ados, method = "spearman", exact = FALSE))

# ADOS vs outlier count in ASD: rho ~ -0.04, p ~ 0.72 (no relationship)

# Takeaway: a handful of subjects are clear outliers; they lean TD (9 of 13 vs 4 of 13)
# but the difference is not statistically significant, and outlier count isn't
# related to symptom severity.
