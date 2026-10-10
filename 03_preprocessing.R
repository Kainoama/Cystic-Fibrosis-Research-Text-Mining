# Part 3: Pre-processing

library(tidyverse)
library(tidytext)
library(textstem)

# load in necessary files from Part 2
cf_all_deduped <- readRDS("cf_all_deduped.rds")

# extra stopwords to be removed form abstracts
supplementary_stopwords <- data.frame(word = 
                                        c("pmid", "collaborators", "affiliations",
                                          "background", "methods", "results")
)

all_stopwords <- bind_rows(supplementary_stopwords, stop_words)

# pre-processing step 1: lowercasing & punctuation removal
cf_all_clean <- mutate(cf_all_deduped, AB_lower = str_remove_all(tolower(AB), "[[:punct:]]")) # str_remove_all for new lower-cased abstracts (AB_lower)
                                                                                              # [[:punct:]] to remove punctuation
# pre-processing step 2: tokenisation
tokenised_AB <- unnest_tokens(cf_all_clean, tokenisation, AB_lower)

# pre-processing step 3: stopword removal
clean_tokenised_AB <- anti_join(tokenised_AB, all_stopwords, by = c("tokenisation" = "word"))

# pre-processing step 4: lemmatisation
cf_all <- mutate(clean_tokenised_AB, lemmatisation = lemmatize_words(tokenisation))

table(cf_all_deduped$time_period)  # unique article counts per period

# preparation for n-gram analysis
## standardise US spelling to fit UK spelling
spelling_standard <- c(
  "pediatric" = "paediatric", "pediatrics" = "paediatrics", "pediatrician" = "paediatrician", 
  "pediatricians" = "paediatricians", "tumor" = "tumour", "anemia" = "anaemia", 
  "esophageal" = "oesophageal", "hemoglobin" = "haemoglobin", "diarrhea" = "diarrhoea", 
  "leukemia" = "leukaemia", "edema" = "oedema", "fetal" = "foetal", "fetus" = "foetus", 
  "color" = "colour", "favor" = "favour", "labor" = "labour", "behavior" = "behaviour", 
  "defense" = "defence", "center" = "centre", "fiber" = "fibre", "liter" = "litre", 
  "catalog" = "catalogue", "dialog" = "dialogue", "license" = "licence", "aging" = "ageing", 
  "characterize" = "characterise", "characterized" = "characterised",
  "organize" = "organise", "organized" = "organised", "randomize" = "randomise", 
  "randomized" = "randomised", "analyze" = "analyse", "analyzed" = "analysed", 
  "standardize" = "standardise", "standardization" = "standardisation",
  "summarize" = "summarise", "recognize" = "recognise", "categorize" = "categorise",
  "maximize" = "maximise", "minimize" = "minimise", "realize" = "realise",
  "emphasize" = "emphasise", "prioritize" = "prioritise", "utilize" = "utilise",
  "mobilize" = "mobilise", "stabilize" = "stabilise", "sensitize" = "sensitise",
  "customize" = "customise", "specialize" = "specialise", "homogenize" = "homogenise",
  "localize" = "localise", "localization" = "localisation", "immunization" = "immunisation",
  "generalization" = "generalisation", "optimization" = "optimisation",
  "hospitalization" = "hospitalisation", "practice" = "practise", "modeling" = "modelling",
  "signaling" = "signalling", "program" = "programme"
)

# n-gram analysis requires no tokenisation or lemmatisation (to preserve context)
cf_all_raw <- cf_all_deduped %>%
  mutate(AB_lower = str_remove_all(tolower(AB), "[[:punct:]]")) %>%
  mutate(AB_lower = str_replace_all(AB_lower, setNames(spelling_standard, paste0("\\b", names(spelling_standard), "\\b"))))

# save so Part 4 etc. doesn't need this script re-ran
saveRDS(all_stopwords, "all_stopwords.rds")
saveRDS(cf_all, "cf_all.rds")
saveRDS(cf_all_raw, "cf_all_raw.rds")
saveRDS(spelling_standard, "spelling_standard.rds")