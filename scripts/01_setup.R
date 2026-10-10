# Part 1: Setup

# MEMEORY ISSUE! more than default 16gb R memory is required to support lemmatisation of the full corpus size

# increased to 32gb using .Renviron (for mac). This only needs to be done once, not every session!
# 1. file.edit("~/.Renviron")
# 2. add R_MAX_VSIZE = 32gb
# 3. save file and close, and restart R
# 4. check to see if it worked using mem.maxVsize() (should say about 32000,  not 16000)

# all packages needed for scripts 1-6 (and extra)
install.packages("pubmedR")
install.packages("tidyverse")
install.packages("tidytext")
install.packages("igraph")
install.packages("textstem")
install.packages("tidylo")
install.packages("gt")
install.packages("ggraph")
install.packages("quanteda")

library(pubmedR)
library(tidyverse)
library(tidytext)
library(igraph)
library(textstem)
library(tidylo)

# not required (I used one though)...enter your own PubdMed API key here for faster retrieval for records
api_key <- "insert_api_key_here"

# all 11 search subsets, each starting with same base cf/cftr modulator terms and with english/humans filters
cf_subsets <- list(
  parent = '("cystic fibrosis"[tiab] OR "CFTR modul*"[tiab] OR mucoviscidosis[tiab]) AND english[lang] AND humans[MeSH Terms]',
  lung = '("cystic fibrosis"[tiab] OR "CFTR modul*"[tiab] OR mucoviscidosis[tiab]) AND (lung[tiab] OR airway*[tiab] OR pulmonary[tiab] OR respirat*[tiab]) AND english[lang] AND humans[MeSH Terms]',
  liver = '("cystic fibrosis"[tiab] OR "CFTR modul*"[tiab] OR mucoviscidosis[tiab]) AND (liver[tiab] OR hepat*[tiab]) AND english[lang] AND humans[MeSH Terms]',
  pancreas = '("cystic fibrosis"[tiab] OR "CFTR modul*"[tiab] OR mucoviscidosis[tiab]) AND pancrea*[tiab] AND english[lang] AND humans[MeSH Terms]',
  treatment = '(“cystic fibrosis”[tiab] OR "CFTR modul*"[tiab] OR mucoviscidosis[tiab]) AND (treatment[tiab] OR therap*[tiab] OR manag*[tiab] OR care[tiab]) AND english[lang] AND humans[MeSH Terms]',
  life_expectancy = '(“cystic fibrosis”[tiab] OR "CFTR modul*"[tiab] OR mucoviscidosis[tiab]) AND (“life expectancy”[tiab] OR surviv*[tiab] OR mortalit*[tiab]) AND english[lang] AND humans[MeSH Terms]',
  modulators = '("cystic fibrosis"[tiab] OR "CFTR modul*"[tiab] OR mucoviscidosis[tiab]) AND (ivacaftor[tiab] OR lumacaftor[tiab] OR elexacaftor[tiab] OR tezacaftor[tiab] OR kaftrio[tiab] OR orkambi[tiab] OR trikafta[tiab] OR symkevi[tiab] OR symdeko[tiab] OR kalydeco[tiab] OR vanzacaftor[tiab] OR alyftrek[tiab]) AND english[lang] AND humans[MeSH Terms]',
  aminoglycosides = '("cystic fibrosis"[tiab] OR "CFTR modul*"[tiab] OR mucoviscidosis[tiab]) AND (aminoglycoside*[tiab] OR tobramycin[tiab] OR amikacin[tiab]) AND english[lang] AND humans[MeSH Terms]',
  enzyme_dornase = '("cystic fibrosis"[tiab] OR "CFTR modul*"[tiab] OR mucoviscidosis[tiab]) AND (dornase*[tiab] OR DNase[tiab] OR pulmozyme[tiab] OR hypertonic*[tiab] OR mannitol[tiab] OR mucoactiv*[tiab] OR mucolytic[tiab]) AND english[lang] AND humans[MeSH Terms]',
  adverse_event = '("cystic fibrosis"[tiab] OR "CFTR modul*"[tiab] OR mucoviscidosis[tiab]) AND (“adverse event*”[tiab] OR “side effect*”[tiab] OR toxicit*[tiab] OR “adverse effect*”[tiab] OR safety[tiab]) AND english[lang] AND humans[MeSH Terms]',
  mental_health = '("cystic fibrosis"[tiab] OR "CFTR modul*"[tiab] OR mucoviscidosis[tiab]) AND ("mental health"[tiab] OR depression[tiab] OR anxiety[tiab]) AND english[lang] AND humans[MeSH Terms]'
)

# list of quinquennial time periods from 1970-2024
time_periods <- list(
  "1970-1974" = c("1970/01/01", "1974/12/31"),
  "1975-1979" = c("1975/01/01", "1979/12/31"),
  "1980-1984" = c("1980/01/01", "1984/12/31"),
  "1985-1989" = c("1985/01/01", "1989/12/31"),
  "1990-1994" = c("1990/01/01", "1994/12/31"),
  "1995-1999" = c("1995/01/01", "1999/12/31"),
  "2000-2004" = c("2000/01/01", "2004/12/31"),
  "2005-2009" = c("2005/01/01", "2009/12/31"),
  "2010-2014" = c("2010/01/01", "2014/12/31"),
  "2015-2019" = c("2015/01/01", "2019/12/31"),
  "2020-2024" = c("2020/01/01", "2024/12/31")
)

# save so Part 2 etc. doesn't need this script re-ran
saveRDS(cf_subsets, "cf_subsets.rds")
saveRDS(time_periods, "time_periods.rds")
