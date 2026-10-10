# Part 2: Corpus Retrieval

library(tidyverse)
library(pubmedR)

# load in necessary files from Part 1
cf_subsets <- readRDS("cf_subsets.rds")
time_periods <- readRDS("time_periods.rds")

# fetching pubmed data, with 121 queries total (11 subsets, 11 periods)
results <- list()

for (x in names(cf_subsets)) {
  for (y in names(time_periods)) {
    pubmed_query <- paste0(cf_subsets[[x]], # insert cf_subsets to test 1 subset before running all
                           ' AND ("',
                           time_periods[[y]][1], # insert specific time_periods to test 1 subset
                           '"[dp] : "',
                           time_periods[[y]][2],
                           '"[dp])'
    )
    # get total record count
    total_count <- pmQueryTotalCount(query = pubmed_query, api_key = NULL) # input API key if you want
    
    # pull actual records
    records <- pmApiRequest(query = pubmed_query, limit = total_count$total) # pubmed has max limit of 9999 records per request, so batch by time periods
    
    # store results
    cf_corpus_results <- pmApi2df(records)
    
    # convert to a df
    results[[x]][[y]] <- cf_corpus_results
  }
}

# add the originating subset/period to each result before combining
for (x in names(cf_subsets)) {
  for (y in names(time_periods)) {
    results[[x]][[y]]<- mutate(results[[x]][[y]], cf_subset = x, time_period = y)
    results[[x]][[y]] <- select(results[[x]][[y]],-any_of("country")) # all values are N/A in "country" column; drop it
  }
}

# deduplication: the 11 subsets aren't mutually exclusive (same PMID shows up more than once across results...)
cf_all_combined <- bind_rows(lapply(names(cf_subsets), function(x) bind_rows(results[[x]])))
cf_all_deduped <- cf_all_combined %>% distinct(PMID, .keep_all = TRUE)

nrow(cf_all_combined) # pre-deduplication count
nrow(cf_all_deduped) # actual unique articles (USE THIS FOR ANALYSIS GOING FORWARD)

# save so Part 3 etc. doesn't need this script re-ran
saveRDS(results, "results.rds")
saveRDS(cf_all_combined, "cf_all_combined.rds")
saveRDS(cf_all_deduped, "cf_all_deduped.rds")

# composition table (i.e. how many articles per period, post dedpulication)
periods_ordered <- names(time_periods)

corpus_comp <- cf_deduped %>%
  count(time_period, name = "n_articles") %>%
  mutate(time_period = factor(time_period, levels = ordered_periods)) %>%
  arrange(time_period) %>%
  mutate(pct_of_total = round(100 * n_articles / sum(n_articles), 1), cumulative_n = cumsum(n_articles)
  )

stopifnot(sum(corpus_composition$n_articles) == nrow(cf_all_deduped)) # should always match, just checking

# totals row for table
comp_total <- corpus_comp %>%
  select(-cumulative_n) %>%
  bind_rows(tibble(time_period = "Total",
                   n_articles = sum(corpus_comp$n_articles),
                   pct_of_total = 100.0
    )
  )

comp_total

# create into table
library(gt)

comp_gt <- comp_total %>%
  gt() %>%
  tab_header(title = "Number of Articles Retrieved per Quinquennial Period",
             subtitle = "After deduplication on PMID") %>%
  cols_label(time_period = "Period",
             n_articles = "Articles (n)",
             pct_of_total = "% of Total") %>%
  fmt_number(columns = n_articles, decimals = 0, use_seps = TRUE) %>%
  fmt_number(columns = pct_of_total, decimals = 1) %>%
  cols_align(align = "center", columns = everything())

comp_gt

gtsave(corpus_composition_gt, "corpus_composition_table.png") # save as png

saveRDS(corpus_comp, "corpus_composition.rds")

# plotting articles per period
ggplot(corpus_composition, aes(x = time_period, y = n_articles)) +
  geom_col(fill = "#2c3e50") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Number of Articles per Quinquennial Period",
       x = "Period", y = "Number of Articles")