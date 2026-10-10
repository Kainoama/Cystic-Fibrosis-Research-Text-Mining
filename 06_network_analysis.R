# Part 6: Network Analysis

library(tidyverse)
library(ggraph)
library(quanteda)
library(igraph)

# load in necessary files from part 5
unigram_counts <- readRDS("unigram_counts.rds")
unigram_logodds_filtered <- readRDS("unigram_logodds_filtered.rds")
bigram_logodds_filtered <- readRDS("bigram_logodds_filtered.rds")
trigram_logodds_filtered <- readRDS("trigram_logodds_filtered.rds")
cf_all_raw <- readRDS("cf_all_raw.rds")
time_periods <- readRDS("time_periods.rds")
network_noise <- readRDS("network_noise.rds")
variant_collapse <- readRDS("variant_collapse.rds")
all_stopwords <- readRDS("all_stopwords.rds")

# build network vocab (treatment and adverse event terms)
# keep these as separate lists so the nodes can later be labelled as either "treatment"
# or "adverse events" in the network graph

treatment_terms <- c(
  # cftr modulators/brands
  "ivacaftor", "lumacaftor", "elexacaftor", "tezacaftor", "vanzacaftor",
  "kaftrio", "orkambi", "trikafta", "symkevi", "symdeko", "kalydeco", "alyftrek",
  "cftr modulator", "cftr modulators",
  # aminoglycosides
  "aminoglycoside", "aminoglycosides", "tobramycin", "amikacin",
  # mucolytics/enzyme therapies
  "dornase alfa", "dornase", "dnase", "pulmozyme", "hypertonic saline",
  "mannitol", "mucoactive", "mucolytic"
)

ae_terms <- c(
  "adverse event", "adverse events", "side effect", "side effects",
  "adverse effect", "adverse effects", "toxicity", "safety",
  "mental health", "depression", "anxiety"
)

base_terms <- c(treatment_terms, ae_terms)
length(base_terms)

# check which ae/treatment terms are available in the uni/bi/trigram data
base_present_uni <- base_terms[base_terms %in% unique(unigram_counts$unigrams)]
base_present_bi  <- base_terms[base_terms %in% unique(bigram_logodds_filtered$bigrams)]
base_present_tri <- base_terms[base_terms %in% unique(trigram_logodds_filtered$trigrams)]

base_present <- unique(c(base_present_uni, base_present_bi, base_present_tri))

length(base_terms)   # total seeded
length(base_present) # how many actually occur in the corpus

# check which of the original terms dont appear
setdiff(base_terms, base_present)
treatment_present <- intersect(treatment_terms, base_present)
ae_present <- intersect(ae_terms, base_present)

# widen candidate pool to top 50 per period (top 10 per period was too narrow to find
# ae/treatment terms specifically since general cf/methodology terms dominated
candidate_bigrams <- bigram_logodds_filtered %>%
  group_by(time_period) %>% 
  slice_max(log_odds_weighted, n = 50) %>%
  ungroup() %>% 
  distinct(bigrams) %>% 
  pull(bigrams)

candidate_trigrams <- trigram_logodds_filtered %>%
  group_by(time_period) %>%
  slice_max(log_odds_weighted, n = 50) %>%
  ungroup() %>% 
  distinct(trigrams) %>% 
  pull(trigrams)

candidate_unigrams <- unigram_logodds_filtered %>%
  group_by(time_period) %>%
  slice_max(log_odds_weighted, n = 50) %>%
  ungroup() %>%
  distinct(unigrams) %>% 
  pull(unigrams)

candidate_terms <- unique(c(candidate_unigrams, candidate_bigrams, candidate_trigrams))
length(candidate_terms)

# create a document term matrix
docs_for_network <- cf_all_raw %>%
  distinct(PMID, .keep_all = TRUE) %>%
  filter(!is.na(AB_lower), AB_lower != "") %>%
  select(PMID, time_period, AB_lower)

nrow(docs_for_network)

dtm_vocab <- unique(c(base_present, candidate_terms))
length(dtm_vocab)

corp <- corpus(docs_for_network, text_field = "AB_lower", docid_field = "PMID")

toks <- tokens(corp) %>% # keep multi-word terms together
  tokens_compound(pattern = phrase(dtm_vocab), concatenator = " ")

dfm_terms_full <- dfm(toks) %>%
  dfm_select(pattern = dtm_vocab, valuetype = "fixed")

dim(dfm_terms_full)

# finds terms that co-occur with the base terms
X_bin <- (as.matrix(dfm_terms_full) > 0) * 1
cooc_full <- t(X_bin) %*% X_bin
diag(cooc_full) <- 0

base_in_dtm <- intersect(base_present, colnames(cooc_full))
stopifnot(length(base_in_dtm) > 0) # double checker

candidate_cooc_with_base <- tibble(
  term = colnames(cooc_full),
  max_cooc_with_base = apply(cooc_full[, base_in_dtm, drop = FALSE], 1, max)
) %>%
  filter(!term %in% base_present)

candidate_cooc_with_base %>% arrange(desc(max_cooc_with_base)) %>% head(30)

# keep terms that co-occur with a base term more than 30 times
expanded_terms <- candidate_cooc_with_base %>%
  filter(max_cooc_with_base > 30) %>%
  pull(term)

length(expanded_terms)
expanded_terms

# combine the OG and expanded terms
network_vocab_raw <- unique(c(base_present, expanded_terms))
length(network_vocab_raw)


# clean the network vocabulary, checking for abbreivations
cf_abbrev_pairs <- list(
  c("cystic fibrosis", "cf"), c("cystic fibrosis", "pwcf"),
  c("cftr modulator", "eti"), c("elexacaftortezacaftorivacaftor", "eti"),
  c("forced expiratory volume", "fev"), c("body mass index", "bmi"),
  c("pseudomonas aeruginosa", "pa"), c("confidence interval", "ci"),
  c("cystic fibrosis transmembrane conductance regulator", "cftr")
)

gloss_artifact_candidates <- network_vocab_raw[
  sapply(network_vocab_raw, function(term) {
    any(sapply(cf_abbrev_pairs, function(pair) {
      str_detect(term, fixed(pair[1])) & str_detect(term, paste0("\\b", pair[2], "\\b"))
    }))
  })
]
gloss_artifact_candidates # confirmed artifacts -> add to network_noise below

# in network, drop any term where any word is a stopword
network_vocab_strict <- network_vocab_raw[
  sapply(network_vocab_raw, function(term) {
    words <- str_split(term, " ")[[1]]
    all(!words %in% all_stopwords$word)
  })
]
length(network_vocab_strict)

# update variant_collapse to merge single and plural forms show as separate graph nodes
variant_collapse <- c(variant_collapse, c(
  "infections" = "infection",
  "mutations" = "mutation",
  "modulators" = "modulator",
  "cftr modulators" = "cftr modulator"
))
saveRDS(variant_collapse, "variant_collapse.rds")

network_vocab_collapsed <- recode(network_vocab_strict, !!!variant_collapse) %>% unique()
length(network_vocab_collapsed)

# remove excluded terms, compare final vocab size with OG, and check that no excluded terms are still here
network_vocab <- setdiff(network_vocab_collapsed, network_noise)
length(network_vocab)

stopifnot(length(intersect(network_vocab, network_noise)) == 0) # should be empty

# start to build network graph
# step 1: make final document-term matrix
toks_final <- tokens(corp) %>%
  tokens_compound(pattern = phrase(network_vocab), concatenator = " ")

dfm_terms <- dfm(toks_final) %>%
  dfm_select(pattern = network_vocab, valuetype = "fixed")

docvars(dfm_terms, "time_period") <- docs_for_network$time_period
dim(dfm_terms)

# step 2: calculate co-occurences for each period
periods <- unique(docvars(dfm_terms, "time_period"))

cooc_by_period <- map_dfr(periods, function(p) {
  dfm_p <- dfm_subset(dfm_terms, time_period == p)
  X_bin <- (as.matrix(dfm_p) > 0) * 1
  
  cooc_matrix <- t(X_bin) %*% X_bin
  diag(cooc_matrix) <- 0
  cooc_matrix[lower.tri(cooc_matrix)] <- 0
  
  term_names <- colnames(X_bin)
  idx <- which(cooc_matrix > 0, arr.ind = TRUE)
  
  if (nrow(idx) == 0) {
    return(tibble(term1 = character(), term2 = character(), n = integer(), time_period = character()))
  }
  
  tibble(term1 = term_names[idx[, 1]], term2 = term_names[idx[, 2]],
         n = cooc_matrix[idx], time_period = p)
})

# step 3: keep edges at >= 30 co-occurences and require each edge to include a base term
sapply(c(2, 3, 5, 10, 20, 30), function(t) sum(cooc_by_period$n >= t))

cooc_thresh <- cooc_by_period %>% # keep stronger connections that include treatment/ae term
  filter(n >= 30) %>%
  filter(term1 %in% base_present | term2 %in% base_present)

nrow(cooc_thresh)
cooc_thresh %>% distinct(term1, term2) %>% nrow()

# step 4: build the graph
cf_network <- graph_from_data_frame(cooc_thresh, directed = FALSE)
E(cf_network)$time_period <- factor(E(cf_network)$time_period, levels = names(time_periods))

V(cf_network)$degree <- degree(cf_network)
V(cf_network)$category <- case_when(
  V(cf_network)$name %in% treatment_present ~ "Treatment",
  V(cf_network)$name %in% ae_present ~ "Adverse Event",
  TRUE ~ "Associated term"
)

# step 5: keep largest connected section of network
comp <- components(cf_network)
cf_network_main <- induced_subgraph(cf_network, which(comp$membership == which.max(comp$csize)))
vcount(cf_network_main); ecount(cf_network_main)


# setup network graph
set.seed(123)

# identify key treatment/mental health edges
# find edges directly linking to cftr modulator treatment to a mental health ae
cftr_modulator_full_list <- c(
  "ivacaftor", "lumacaftor", "elexacaftor", "tezacaftor", "vanzacaftor",
  "kaftrio", "orkambi", "trikafta", "symkevi", "symdeko", "kalydeco", "alyftrek",
  "cftr modulator", "cftr modulators"
)
cftr_modulator_terms <- treatment_present[treatment_present %in% cftr_modulator_full_list]
cftr_modulator_terms # which modulator terms actually appear in the corpus

mental_health_terms <- c("mental health", "depression", "anxiety")

# mark edges that connect a cftr modulator w/a mental health term
E(cf_network_main)$highlight <- (
  (get.edgelist(cf_network_main)[,1] %in% cftr_modulator_terms &
     get.edgelist(cf_network_main)[,2] %in% mental_health_terms) |
    (get.edgelist(cf_network_main)[,2] %in% cftr_modulator_terms &
       get.edgelist(cf_network_main)[,1] %in% mental_health_terms)
)

sum(E(cf_network_main)$highlight) # how many aformentioned edges exist

# edge colours: highlighted relationships are shown in black, while all others
# edges are colored
E(cf_network_main)$edge_colour_group <- ifelse(
  E(cf_network_main)$highlight,
  "CFTR modulator <-> Mental health",
  as.character(E(cf_network_main)$time_period)
)

period_colours <- c(
  "1970-1974" = "#e41a1c", "1975-1979" = "#377eb8", "1980-1984" = "#4daf4a",
  "1985-1989" = "#984ea3", "1990-1994" = "#ff7f00", "1995-1999" = "#ffe100",
  "2000-2004" = "#a65628", "2005-2009" = "#f781bf", "2010-2014" = "#00bcd4",
  "2015-2019" = "#66c2a5", "2020-2024" = "#fc4e2a"
)

edge_colour_values <- c(period_colours, "CFTR modulator <-> Mental health" = "black")

# edge opacity/width increase with co-occurence strength, thus stronger relationships
# should stand out over weaker ones
ggraph(cf_network_main, layout = "stress") +
  geom_edge_link(aes(colour = edge_colour_group, width = n, alpha = n)) +
  scale_edge_width(range = c(0.3, 3)) +
  scale_edge_alpha(range = c(0.08, 0.95), guide = "none") +
  scale_edge_colour_manual(values = edge_colour_values) +
  geom_node_point(aes(size = degree, fill = category), shape = 21, colour = "grey20", stroke = 0.4) +
  geom_node_text(
    aes(label = name),
    repel = TRUE, size = 3, max.overlaps = Inf, force = 5, max.time = 10
  ) +
  scale_size(range = c(1, 10)) +
  scale_fill_manual(values = c(
    "Treatment" = "#2ca02c", "Adverse Event" = "#d62728", "Associated term" = "#a6a6a6"
  )) +
  theme_void() +
  labs(
    title = "From Symptomatic Management to Targeted Therapy: Treatment-Adverse Event\nAssociations in CF Literature (1970-2024)",
    subtitle = "Edge opacity scales with co-occurrence strength; black edges highlight\nCFTR modulator-mental health co-occurrence (Talwalkar et al., 2017).",
    edge_colour = "Period / Highlight", edge_width = "Co-occurrences",
    size = "Degree", fill = "Category"
  )

ggsave("cf_network.png", width = 20, height = 15, dpi = 300)

# save everything just in case
saveRDS(treatment_terms, "treatment_terms.rds")
saveRDS(ae_terms, "ae_terms.rds")
saveRDS(base_terms, "base_terms.rds")
saveRDS(base_present, "base_present.rds")
saveRDS(treatment_present, "treatment_present.rds")
saveRDS(ae_present, "ae_present.rds")
saveRDS(expanded_terms, "expanded_terms.rds")
saveRDS(network_vocab_raw, "network_vocab_raw.rds")
saveRDS(network_vocab, "network_vocab.rds")
saveRDS(cooc_by_period, "cooc_by_period.rds")
saveRDS(cooc_thresh, "cooc_thresh.rds")
saveRDS(cf_network, "cf_network.rds")
saveRDS(cf_network_main, "cf_network_main.rds")

# fixed-term frequency trends and period over period change (i.e. track the 
# treatment/ae terms across all periods)
## load necessary files from part 5
bigram_counts <- readRDS("bigram_counts.rds")
trigram_counts <- readRDS("trigram_counts.rds")

articles_per_period <- readRDS("corpus_composition.rds") %>%
  select(time_period, number_of_articles = n_articles)

ordered_periods <- names(time_periods)

# combine treatment and adverse event terms into 1 fixed list
fixed_terms <- c(treatment_terms, ae_terms)

# combine singular/plural variants (should be treated as same concept_)
fixed_term_collapse <- c(
  "adverse effects" = "adverse effect",
  "adverse events" = "adverse event",
  "side effects" = "side effect",
  "aminoglycosides" = "aminoglycoside",
  "cftr modulators" = "cftr modulator",
  "dornase" = "dornase alfa"
)

fixed_terms_canonical <- unique(recode(fixed_terms, !!!fixed_term_collapse)) # remove and duplicates

# create clean term lists separately for treatment and ae
treatment_terms_canonical <- unique(recode(treatment_terms, !!!fixed_term_collapse))
ae_terms_canonical <- unique(recode(ae_terms, !!!fixed_term_collapse))

# combine counts acorss all n-gram lengths (contains uni/bi/trigrams)
term_counts_raw <- bind_rows(
  unigram_counts %>% filter(unigrams %in% fixed_terms) %>% rename(term = unigrams),
  bigram_counts %>% filter(bigrams %in% fixed_terms) %>% rename(term = bigrams),
  trigram_counts %>% filter(trigrams %in% fixed_terms) %>% rename(term = trigrams)
) %>%
  mutate(term = recode(term, !!!fixed_term_collapse)) %>%
  group_by(time_period, term) %>%
  summarise(n = sum(n), .groups = "drop") %>% # merge singular and plural counts together
  mutate(time_period = factor(time_period, levels = ordered_periods)) %>%
  complete(term = fixed_terms_canonical, time_period = ordered_periods, fill = list(n = 0)) %>% # add rows for terms that are absent (give them a 0)
  mutate(time_period = factor(time_period, levels = ordered_periods)) %>%   # complete() drops the factor levels ordering above (restore it after filling gaps)
  left_join(articles_per_period, by = "time_period") %>%
  mutate(rate_per_article = n / number_of_articles) %>%
  mutate(category = case_when(
    term %in% treatment_terms_canonical ~ "Treatment", # each term will be labelled treatment or ae
    term %in% ae_terms_canonical ~ "Adverse Event"
  ))

# save
saveRDS(term_counts_raw, "term_counts_raw.rds")

# period over period change, calculating change in mean term rate from each period to the preceding period,
# separately for treatments and ae
category_mean_by_period <- term_counts_raw %>%
  group_by(category, time_period) %>%
  summarise(mean_rate_per_article = mean(rate_per_article), .groups = "drop") %>%
  arrange(category, time_period)

period_over_period_change <- category_mean_by_period %>%
  group_by(category) %>%
  mutate(
    change_vs_prior_period = mean_rate_per_article - lag(mean_rate_per_article),
    pct_change_vs_prior_period = round(100 * change_vs_prior_period / lag(mean_rate_per_article), 1)
  ) %>%
  ungroup()

period_over_period_change

# save
saveRDS(period_over_period_change, "period_over_period_change.rds")

## period-over-period change bargraph
period_over_period_change %>%
  filter(!is.na(change_vs_prior_period)) %>%
  ggplot(aes(x = time_period, y = change_vs_prior_period, fill = category)) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = 0, colour = "grey40") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Period-over-Period Change in Mean Term Rate (vs. Preceding Period)",
       x = "Quinquennial Period", y = "Change in Mean Occurrences per Article",
       fill = NULL)


# mental health edge trend, tracking cftr modul and mental health term co-occurence over time
# (using unthresholded data)

# keep only edges linking a cftr modulator to a mental health term
mental_health_edge_trend <- cooc_by_period %>%
  filter(
    (term1 %in% cftr_modulator_terms & term2 %in% mental_health_terms) |
      (term2 %in% cftr_modulator_terms & term1 %in% mental_health_terms)
  ) %>%
  mutate(
    time_period = factor(time_period, levels = ordered_periods),
    edge_label = paste(term1, "<->", term2)
  ) %>%
  arrange(time_period, desc(n))

mental_health_edge_trend

# complete missing periods by adding 0 count rows for periods where edge pair doesnt occur
mental_health_edge_trend_complete <- mental_health_edge_trend %>%
  complete(edge_label, time_period = ordered_periods, fill = list(n = 0)) %>%
  mutate(time_period = factor(time_period, levels = ordered_periods))

# save
saveRDS(mental_health_edge_trend_complete, "mental_health_edge_trend.rds")

## trend line per edge pair graph (only 2010 onward bc prior years are practically at 0)
mental_health_edge_trend_complete %>%
  filter(time_period %in% c("2010-2014", "2015-2019", "2020-2024")) %>%
  filter(edge_label %in% (mental_health_edge_trend_complete %>%
                            group_by(edge_label) %>%
                            summarise(total = sum(n)) %>%
                            filter(total > 0) %>%
                            pull(edge_label))) %>%
  ggplot(aes(x = time_period, y = n, colour = edge_label, group = edge_label)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom") +
  guides(colour = guide_legend(ncol = 2)) +
  labs(title = "CFTR Modulator <-> Mental Health Co-occurrence Over Time",
       subtitle = "Unthresholded counts - 2010 onwards",
       x = "Quinquennial Period", y = "Co-occurrence Count (n)", colour = NULL)