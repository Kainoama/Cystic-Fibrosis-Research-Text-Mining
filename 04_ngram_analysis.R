# Part 4: N-gram Analysis

library(tidyverse)
library(tidytext)
library(tidylo)

# load in necessary files from Part 3
cf_all_raw <- readRDS("cf_all_raw.rds")
all_stopwords <- readRDS("all_stopwords.rds")
cf_all_deduped <- readRDS("cf_all_deduped.rds")
time_periods <- readRDS("time_periods.rds")

# unigram frequency analysis
unigrams <- unnest_tokens(cf_all_raw, unigrams, AB_lower)

unigrams_clean <- unigrams %>%
  anti_join(all_stopwords, by = c("unigrams" = "word")) %>%
  filter(!str_detect(unigrams, "^[0-9]+$")) # remove numeric values

unigram_counts <- unigrams_clean %>%
  count(time_period, unigrams, sort = TRUE)

# bigram frequency analysis
bigrams <- unnest_tokens(cf_all_raw, bigrams, AB_lower, token = "ngrams", n = 2)

bigram_counts <- bigrams %>%
  count(time_period, bigrams, sort = TRUE)

# trigram frequency analysis
trigrams <- unnest_tokens(cf_all_raw, trigrams, AB_lower, token = "ngrams", n = 3)

trigram_counts <- trigrams %>%
  count(time_period, trigrams, sort = TRUE)

# comparing bigrams vs trigrams with vs. without stopwords to see how vocab/frequency changes
# when context is not being preserved
bigrams_no_stopwords <- bigrams %>%
  separate(bigrams, into = c("word1", "word2"), sep = " ", remove = FALSE) %>%
  filter(!word1 %in% all_stopwords$word, !word2 %in% all_stopwords$word) %>%
  select(-word1, -word2)

bigram_counts_no_stopwords <- bigrams_no_stopwords %>%
  count(time_period, bigrams, sort = TRUE)

trigrams_no_stopwords <- trigrams %>%
  separate(trigrams, into = c("word1", "word2", "word3"), sep = " ", remove = FALSE) %>%
  filter(!word1 %in% all_stopwords$word, !word2 %in% all_stopwords$word, !word3 %in% all_stopwords$word) %>%
  select(-word1, -word2, -word3)

trigram_counts_no_stopwords <- trigrams_no_stopwords %>%
  count(time_period, trigrams, sort = TRUE)

# comparing raw occurrence counts and unique term counts, with vs without stopwords
tibble(metric = c("total bigram occurrences", "unique bigrams", "total trigram occurrences", "unique trigrams"),
       with_stopwords = c(nrow(bigrams), n_distinct(bigrams$bigrams), nrow(trigrams), n_distinct(trigrams$trigrams)),
       without_stopwords = c(nrow(bigrams_no_stopwords), n_distinct(bigrams_no_stopwords$bigrams), 
                             nrow(trigrams_no_stopwords), n_distinct(trigrams_no_stopwords$trigrams))
)

# compare content: top 10 bigrams/trigrams overall, with vs without stopwords
bigram_counts %>% group_by(bigrams) %>% summarise(n = sum(n)) %>% slice_max(n, n = 10)
bigram_counts_no_stopwords %>% group_by(bigrams) %>% summarise(n = sum(n)) %>% slice_max(n, n = 10)

trigram_counts %>% group_by(trigrams) %>% summarise(n = sum(n)) %>% slice_max(n, n = 10)
trigram_counts_no_stopwords %>% group_by(trigrams) %>% summarise(n = sum(n)) %>% slice_max(n, n = 10)


# setting a threshold, checking how many bigrams/trigrams remain at different occurence cutoffs
sapply(c(2, 5, 10, 20, 50), function(t) sum(bigram_counts$n >= t))
sapply(c(2, 5, 10, 20, 50), function(t) sum(trigram_counts$n >= t))
sapply(c(2, 5, 10, 20, 50), function(t) sum(bigram_counts$n[bigram_counts$n >= t]))
sapply(c(2, 5, 10, 20, 50), function(t) sum(trigram_counts$n[trigram_counts$n >= t]))

# testing thresholds of bigrams >= 10, trigrams >=5 (per period)
bigram_counts_thresh <- bigram_counts %>%
  filter(n >= 10) %>%
  filter(!is.na(bigrams)) # empty/too short of abstracts create N/A rows; remove

trigram_counts_thresh <- trigram_counts %>%
  filter(n >= 5) %>%
  filter(!is.na(trigrams))

# double check results are meaningful and not noise
bigram_counts_thresh %>% arrange(desc(n)) %>% head(30)
trigram_counts_thresh %>% arrange(desc(n)) %>% head(30)

set.seed(123)
bigram_counts_thresh %>% sample_n(30)
trigram_counts_thresh %>% sample_n(30)

bigram_counts_thresh %>% filter(time_period == "1970-1974") %>% arrange(desc(n)) %>% head(20)
bigram_counts_thresh %>% filter(time_period == "2020-2024") %>% arrange(desc(n)) %>% head(20)

# rate of bigrams/trigrams per period
## simple calculation of uni/bi/trigrams per article's abstract
articles_per_period <- cf_all_deduped %>%
  count(time_period, name = "number_of_articles")

unigram_rate_by_period <- unigram_counts %>%
  group_by(time_period) %>%
  summarise(total_unigrams = sum(n)) %>%
  left_join(articles_per_period, by = "time_period") %>%
  mutate(unigrams_per_article = total_unigrams/number_of_articles)

bigram_rate_by_period <- bigram_counts_thresh %>%
  group_by(time_period) %>%
  summarise(total_bigrams = sum(n)) %>%
  left_join(articles_per_period, by = "time_period") %>%
  mutate(bigrams_per_article = total_bigrams/number_of_articles)

trigram_rate_by_period <- trigram_counts_thresh %>%
  group_by(time_period) %>%
  summarise(total_trigrams = sum(n)) %>%
  left_join(articles_per_period, by = "time_period") %>%
  mutate(trigrams_per_article = total_trigrams/number_of_articles)

# calculating weighted log-odds 
unigram_logodds <- unigram_counts %>%
  bind_log_odds(time_period, unigrams, n)

bigram_logodds <- bigram_counts_thresh %>%
  bind_log_odds(time_period, bigrams, n)

trigram_logodds <- trigram_counts_thresh %>%
  bind_log_odds(time_period, trigrams, n)

# raw log odds are dominated by disease-name/generic/academic filer (these need to
# be removed if they aren't useful for treatment/AE terms!)

# phrases tied to the disease name
cf_phrases <- c("cystic fibrosis", "with cystic", "in cystic", "of cystic", 
                   "cf patients", "with cf", "in cf", "fibrosis cf",
                   "cystic fibrosis cf", "with cystic fibrosis", "in cystic fibrosis",
                   "of cystic fibrosis", "cystic fibrosis patients", "patients with cf")

# single disease descriptors words
cf_words <- c("cf", "cystic", "fibrosis", "patients", "patient", "disease", "study", "clinical", "pwcf", "people")

# any methodology, abstract-related, noisy terms that are not treatment/AE-related
abstract_noise <- c(
  "quality of", "quality evidence", "risk of bias", "trials register",
  "handsearches of relevant", "the cochrane", "conference proceedingsdate of",
  "of the evidence", "very lowquality evidence", "lowquality evidence",
  "very lowcertainty evidence", "lowcertainty evidence", "the certainty of",
  "authors independently", "we searched the", "review to", "in people",
  "study was to", "of this study", "this study was", "an update",
  "of bias", "less than 001", "less than 0001", "less than 005",
  "than 005", "95 ci", "is suggested that", "cochrane cystic fibrosis",
  "the cochrane cystic", "searched the cochrane", "controlled trials",
  "newborn screening nbs", "of conference proceedingsdate",
  "cystic fibrosis pwcf", "elexacaftortezacaftorivacaftor eti",
  "cf pwcf", "with cf pwcf", "of cf", "cf and", "the cf", "the cystic fibrosis", 
  "cystic fibrosis and", "in cf patients", "children with cf", "for cystic fibrosis", 
  "with cf and", "people with cf", "fibrosis cf is", "fibrosis cf patients",
  
  # found more of these on a second pass through the top-10 charts
  "lowquality", "quasirandomised", "quasirandomised controlled trials",
  "conference abstract books", "abstract books", "haemoglobinopathies trials register", 
  "authors independently assessed", "electronic database searches", "recent search", 
  "clinicaltrialsgov", "significant difference", "significant differences", 
  "statistically significant", "randomised controlled trials", "randomised controlled", 
  "quasirandomised controlled","250 words"
)

# more noise found on second pass
academic_noise <- c(
  "significant", "significantly", "data", "quality", "reported", "review",
  "included", "compared", "outcomes", "increased", "participants", "age",
  "identified", "including", "primary", "effect", "studies", "trials",
  "evidence", "potential",
  
  # second pass; more noise found
  "cochrane", "randomised", "searches", "searched", "control", "difference",
  "analysis", "performed", "found", "cohort", "subjects", "test", "role",
  "strategies", "confidence interval", "95 confidence interval", "ci",
  "reviewers independently", "specialised register", "haemoglobinopathies trials"
)

network_noise <- unique(c(abstract_noise, academic_noise)) # combined list of all n-gram levels

# variations of concepts that need to be collapsed/dombined so they're not treated as
# distinct terms
variant_collapse <- c(
  "δf508 mutation" = "f508del mutation", "delta f508 mutation" = "f508del mutation",
  "the delta f508" = "f508del mutation", "mutation delta f508" = "f508del mutation",
  "delta f508" = "f508del mutation", "deltaf508 mutation" = "f508del mutation",
  "deltaf508 cftr" = "f508del cftr", "of deltaf508cftr" = "f508del cftr",
  "δf508 cftr" = "f508del cftr", "the δf508 mutation" = "f508del mutation",
  "the delta" = "f508del mutation", "the δf508" = "f508del mutation",
  "the deltaf508 mutation" = "f508del mutation","cf chromosomes" = "f508del mutation",
  "of the delta" = "f508del mutation","regulator cftr gene" = "cftr gene", "regulator cftr" = "cftr",
  "cystic fibrosis transmembrane" = "cftr", "conductance regulator cftr" = "cftr",
  "transmembrane conductance" = "cftr", "regulator cftr modulators" = "cftr modulators",
  "the cftr" = "cftr", "the cftr gene" = "cftr gene", "fibrosis transmembrane conductance" = "cftr",
  "conductance regulator" = "cftr", "fibrosis transmembrane" = "cftr"
)

unigram_logodds_filtered <- unigram_logodds %>%
  filter(!unigrams %in% cf_words) %>%
  filter(!unigrams %in% all_stopwords$word) %>%
  filter(!str_detect(unigrams, "^[0-9]+$")) %>%
  filter(!unigrams %in% network_noise) %>%
  filter(!str_detect(unigrams, "[^a-z0-9]")) # catches random symbols that survive punctuation removal

bigram_logodds_filtered <- bigram_logodds %>%
  filter(!bigrams %in% cf_phrases) %>%
  mutate(bigrams = recode(bigrams, !!!variant_collapse)) %>%
  filter(str_count(bigrams, "\\S+") == 2) %>% # drops words that are collapsed to 1 word
  separate(bigrams, into = c("word1", "word2"), sep = " ", remove = FALSE, fill = "right") %>%
  filter(!(word1 %in% all_stopwords$word | word2 %in% all_stopwords$word)) %>%
  filter(!(word1 %in% network_noise | word2 %in% network_noise)) %>%
  select(-word1, -word2) %>%
  group_by(time_period, bigrams) %>%
  summarise(n = sum(n), log_odds_weighted = max(log_odds_weighted), .groups = "drop")

trigram_logodds_filtered <- trigram_logodds %>%
  filter(!trigrams %in% cf_phrases) %>%
  mutate(trigrams = recode(trigrams, !!!variant_collapse)) %>%
  filter(str_count(trigrams, "\\S+") == 3) %>%  # drops words that are collapsed to 1 word
  separate(trigrams, into = c("word1", "word2", "word3"), sep = " ", remove = FALSE, fill = "right") %>%
  filter(!(word1 %in% all_stopwords$word | word2 %in% all_stopwords$word | word3 %in% all_stopwords$word)) %>%
  filter(!(word1 %in% network_noise | word2 %in% network_noise | word3 %in% network_noise)) %>%
  select(-word1, -word2, -word3) %>%
  group_by(time_period, trigrams) %>%
  summarise(n = sum(n), log_odds_weighted = max(log_odds_weighted), .groups = "drop")

# remove descriptive/methodology layer before ranking top 10 weighted log odds 
bigram_logodds_filtered <- bigram_logodds_filtered %>%
  filter(!bigrams %in% network_noise)

trigram_logodds_filtered <- trigram_logodds_filtered %>%
  filter(!trigrams %in% network_noise)

# top 10 per period for each n-gram type
unigram_top10 <- unigram_logodds_filtered %>%
  group_by(time_period) %>%
  slice_max(log_odds_weighted, n = 10) %>%
  select(time_period, term = unigrams, n, log_odds_weighted) %>%
  arrange(time_period, desc(log_odds_weighted))

bigram_top10 <- bigram_logodds_filtered %>%
  group_by(time_period) %>%
  slice_max(log_odds_weighted, n = 10) %>%
  select(time_period, term = bigrams, n, log_odds_weighted) %>%
  arrange(time_period, desc(log_odds_weighted))

trigram_top10 <- trigram_logodds_filtered %>%
  group_by(time_period) %>%
  slice_max(log_odds_weighted, n = 10) %>%
  select(time_period, term = trigrams, n, log_odds_weighted) %>%
  arrange(time_period, desc(log_odds_weighted))

ngram_type_table <- bind_rows(
  unigram_top10 %>% mutate(ngram_type = "unigram"),
  bigram_top10 %>% mutate(ngram_type = "bigram"),
  trigram_top10 %>% mutate(ngram_type = "trigram")
)

ngram_type_table # top 10 by log odds per period

## overall freq summary
ngram_summary <- bind_rows(
  unigram_counts %>% 
    group_by(time_period) %>% 
    summarise(total_occurrences = sum(n), unique_terms = n_distinct(unigrams)) %>% 
    mutate(ngram_type = "unigram"),
  bigram_counts %>% 
    group_by(time_period) %>% 
    summarise(total_occurrences = sum(n), unique_terms = n_distinct(bigrams)) %>% 
    mutate(ngram_type = "bigram"),
  trigram_counts %>% 
    group_by(time_period) %>% 
    summarise(total_occurrences = sum(n), unique_terms = n_distinct(trigrams)) %>% 
    mutate(ngram_type = "trigram")
)
ngram_summary # overall totals (not just top 10)

# barcharts for  top 10 uni, bi, tri-grams
library(ggplot2)

periods_to_show <- names(time_periods)

## for each barchart, the intervals are of 10, but the largest value dictates the x-axis range
unigram_top10 %>%
  filter(time_period %in% 
           periods_to_show) %>%
  mutate(term = fct_reorder(term, log_odds_weighted)) %>%
  ggplot(aes(x = term, y = log_odds_weighted, fill = time_period)) +
  geom_col(show.legend = FALSE) +
  facet_wrap(~time_period, scales = "free", ncol = 2) +
  coord_flip() +
  scale_y_continuous(limits = c(0, 110), breaks = seq(0, 110, by = 10)) +
  labs(title = "Top 10 Unigrams by Weighted Log-Odds", x = NULL, y = "Weighted Log-Odds") +
  theme_minimal(base_size = 12)

bigram_top10 %>%
  filter(time_period %in% periods_to_show,
         time_period != "1970-1974") %>% # 1970-1974 dont have enough bi to meet min freq thresh
  mutate(term = fct_reorder(term, log_odds_weighted)) %>%
  ggplot(aes(x = term, y = log_odds_weighted, fill = time_period)) +
  geom_col(show.legend = FALSE) +
  facet_wrap(~time_period, scales = "free", ncol = 2) +
  coord_flip() +
  scale_y_continuous(limits = c(0, 60), breaks = seq(0, 60, by = 10)) +
  labs(title = "Top 10 Bigrams by Weighted Log-Odds", x = NULL, y = "Weighted Log-Odds") +
  theme_minimal(base_size = 12)

trigram_top10 %>%
  filter(time_period %in% periods_to_show,
         time_period != "1970-1974") %>% # 1970-1974 dont have enough tri to meet min freq thresh
  mutate(term = fct_reorder(term, log_odds_weighted)) %>%
  ggplot(aes(x = term, y = log_odds_weighted, fill = time_period)) +
  geom_col(show.legend = FALSE) +
  facet_wrap(~time_period, scales = "free", ncol = 2) +
  coord_flip() +
  scale_y_continuous(limits = c(0, 50), breaks = seq(0, 50, by = 10)) +
  labs(title = "Top 10 Trigrams by Weighted Log-Odds", x = NULL, y = "Weighted Log-Odds") +
  theme_minimal(base_size = 12)

# save so Part 5/6 doesn't need this script re-ran
saveRDS(unigram_counts, "unigram_counts.rds")
saveRDS(unigram_logodds_filtered, "unigram_logodds_filtered.rds")
saveRDS(bigram_logodds_filtered, "bigram_logodds_filtered.rds")
saveRDS(trigram_logodds_filtered, "trigram_logodds_filtered.rds")
saveRDS(bigram_counts, "bigram_counts.rds")
saveRDS(trigram_counts, "trigram_counts.rds")
saveRDS(ngram_type_table, "ngram_type_table.rds")
saveRDS(network_noise, "network_noise.rds")
saveRDS(variant_collapse, "variant_collapse.rds")