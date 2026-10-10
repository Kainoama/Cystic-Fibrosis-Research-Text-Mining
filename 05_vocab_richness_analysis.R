# Part 5: Vocabulary Richness Analysis

library(tidyverse)
library(ggplot2)

# load in necessary files from Part 4
unigram_counts <- readRDS("unigram_counts.rds")
time_periods <- readRDS("time_periods.rds")

# vocabulary richness analysis
vocab_richness <- unigram_counts %>%
  group_by(time_period) %>%
  summarise(N = sum(n),
            V = n_distinct(unigrams),
            hapax_legomena = sum(n == 1),
            gini = {
              x <- sort(as.numeric(n))
              n_x <- length(x)
              (2 * sum(seq_len(n_x) * x) / (n_x * sum(x))) - (n_x + 1) / n_x
    },
    yules_k = {freq_table <- table(n)
      m <- as.numeric(names(freq_table))
      Vm <- as.numeric(freq_table)
      N_local <- sum(n)
      10000 * (sum(Vm * m^2) - N_local) / N_local^2
    }
  )

vocab_richness

# save for downstream analysis
saveRDS(vocab_richness, "vocab_richness.rds")

# vocab richness table (adds % hapax and orders periods properly)
ordered_periods <- names(time_periods)

vocab_richness_table <- vocab_richness %>%
  mutate(time_period = factor(time_period, levels = ordered_periods)) %>%
  arrange(time_period) %>%
  mutate(pct_hapax = round(100 * hapax_legomena / V, 1),
         gini = round(gini, 3),
         yules_k = round(yules_k, 2)) %>%
  select(time_period, N, V, hapax_legomena, pct_hapax, gini, yules_k)

saveRDS(vocab_richness_table, "vocab_richness_table.rds")

# log odds vs. gini
## median weighted log odds of top 20 terms per period and n-gram type are correlated w Gini

# load in necessary files from part 4
unigram_logodds_filtered <- readRDS("unigram_logodds_filtered.rds")
bigram_logodds_filtered <- readRDS("bigram_logodds_filtered.rds")
trigram_logodds_filtered <- readRDS("trigram_logodds_filtered.rds")

unigram_median <- unigram_logodds_filtered %>%
  group_by(time_period) %>%
  slice_max(log_odds_weighted, n = 20, with_ties = FALSE) %>%
  summarise(median_log_odds = median(log_odds_weighted), .groups = "drop") %>%
  mutate(ngram_type = "unigram")

bigram_median <- bigram_logodds_filtered %>%
  group_by(time_period) %>%
  slice_max(log_odds_weighted, n = 20, with_ties = FALSE) %>%
  summarise(median_log_odds = median(log_odds_weighted), .groups = "drop") %>%
  mutate(ngram_type = "bigram")

trigram_median <- trigram_logodds_filtered %>%
  group_by(time_period) %>%
  slice_max(log_odds_weighted, n = 20, with_ties = FALSE) %>%
  summarise(median_log_odds = median(log_odds_weighted), .groups = "drop") %>%
  mutate(ngram_type = "trigram")

median_logodds_all <- bind_rows(unigram_median, bigram_median, trigram_median) # combine all n-gram types

# join median log odds with gini per period
logodds_gini <- median_logodds_all %>%
  left_join(vocab_richness %>% 
              select(time_period, gini), by = "time_period") %>%
  mutate(time_period = factor(time_period, levels = ordered_periods)) %>%
  arrange(ngram_type, time_period)

logodds_gini

saveRDS(logodds_gini, "logodds_gini.rds") # save

# correlation summary per ngram type
correlation_summary <- logodds_gini %>%
  group_by(ngram_type) %>%
  summarise(r = cor(median_log_odds, gini, method = "spearman"),
            p_value = cor.test(median_log_odds, gini, method = "spearman")$p.value,
            .groups = "drop") %>%
  mutate(ngram_type = factor(ngram_type, levels = c("unigram", "bigram", "trigram")),
         r = round(r, 3),
         p_value = round(p_value, 3),
         significant = p_value < 0.05) %>%
  arrange(ngram_type)

correlation_summary

saveRDS(correlation_summary, "correlation_summary.rds") # save

# gini trend over time graph
ggplot(vocab_richness, aes(x = factor(time_period, levels = ordered_periods), y = gini, group = 1)) +
  geom_line(linewidth = 1, colour = "#2c3e50") +
  geom_point(size = 2, colour = "#2c3e50") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Gini Coefficient (Vocabulary Inequality) by Period", x = "Quinquennial Period", y = "Gini Coefficient")

# vocabulary richness table
library(gt)

vocab_richness_gt <- vocab_richness_table %>%
  gt() %>%
  cols_label(time_period = "Period", N = "N (Tokens)", V = "V (Unique Types)",
             hapax_legomena = "Hapax Legomena", pct_hapax = "% Hapax", gini = "Gini", yules_k = "Yule's K"
             ) %>%
  fmt_number(columns = c(N, V, hapax_legomena), decimals = 0, use_seps = TRUE) %>%
  fmt_number(columns = pct_hapax, decimals = 1) %>%
  fmt_number(columns = gini, decimals = 3) %>%
  fmt_number(columns = yules_k, decimals = 2) %>%
  tab_header(title = "Vocabulary Richness Metrics by Quinquennial Period") %>%
  cols_align(align = "center", columns = everything()) %>%
  cols_align(align = "left", columns = time_period)

vocab_richness_gt

gtsave(vocab_richness_gt, filename = "vocab_richness_table.png", vwidth = 900) # save as png for diss

# vocabulary size and hapax legomena trend, alongside existing gini plot (shows whether vocabulary is
# genuinely diversifying or just scaling with corpus size)
vocab_richness_table %>%
  select(time_period, N, V) %>%
  pivot_longer(cols = c(N, V), names_to = "metric", values_to = "count") %>%
  mutate(metric = recode(metric, N = "Total tokens (N)", V = "Unique types (V)")) %>%
  ggplot(aes(x = time_period, y = count, colour = metric, group = metric)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 14)) +
  labs(title = "Corpus Size (N) vs. Vocabulary Size (V) by Period", x = "Quinquennial Period", y = "Count", colour = NULL)

# yule's k trend graph
ggplot(vocab_richness_table, aes(x = time_period, y = yules_k, group = 1)) +
  geom_line(linewidth = 1, colour = "#8e44ad") +
  geom_point(size = 2, colour = "#8e44ad") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 14)) + 
  labs(title = "Yule's K (Lexical Diversity) by Period", subtitle = "Lower K indicates higher lexical diversity",
       x = "Quinquennial Period", y = "Yule's K")

# MIGHT NEED TO RUN THIS SEPARATE FROM OTHER CODE! For some reason its not producing table if run top to bot
# log odds vs gini correlation results table
library(gt)

# log-odds vs gini correlation table
correlation_gt <- correlation_summary %>%
  gt() %>%
  tab_header(
    title = "Spearman Correlation Between Median Weighted Log-Odds and Gini Coefficient",
    subtitle = "Top 20 terms by period, per n-gram type"
  ) %>%
  cols_label(ngram_type = "N-gram Type", r = "Spearman's rho", p_value = md("*p*-value"), significant = md("Significant (*p* < .05)")
             ) %>%
  fmt_number(columns = c(r, p_value), decimals = 3) %>%
  cols_align(align = "center", columns = everything())

print(correlation_gt)

gtsave(correlation_gt, "logodds_gini_correlation.png") # enter 1 to save it as png