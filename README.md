# Mapping Treatment Evolution and Adverse Events in Cystic Fibrosis: A Corpus-Based Text-Mining and Network Analysis of PubMed Literature

## Overview

For my MSc dissertation in Health Data Science at University College London, I used R to analyse how treatment and adverse-event terminology in cystic fibrosis (CF) research has changed over time. The project applies text mining and network analysis to PubMed titles and abstracts published between 1970 and 2024.

The aim was to identify changes in the language used to discuss CF treatment, measure how the literature's vocabulary evolved, and explore which treatment and adverse-event terms appeared together. The analysis provides a large-scale view of research terminology; it does not assess clinical effectiveness or establish that a treatment causes an adverse event.

## Key Findings

- Analysed **41,222 unique PubMed articles** across 11 five-year periods, from 1970–1974 through 2020–2024.
- Found a clear shift in terminology from general clinical and biological language toward named CFTR modulator therapies and more specific clinical outcomes. Multiple analyses pointed to a marked rise in modulator-related language beginning around 2010–2014, consistent with ivacaftor's 2012 approval.
- Measured a rise in the Gini coefficient from **0.479 to 0.873**, indicating that word use became more concentrated around frequently used terms as the corpus grew.
- Observed a decline in Yule's K from **37.33 to 15.60**, indicating greater lexical diversity by this measure. Together, the metrics suggest that a more concentrated core vocabulary coexisted with an expanding range of specialised terms.
- Built a treatment/adverse-event co-occurrence network with **73 nodes and 366 edges**, identifying clusters around CFTR modulators, treatment safety, infection, and mental-health terminology.
- Found that CFTR modulator and mental-health terms did not co-occur in the same abstract before 2015–2019. In 2020–2024, the most frequent examined pair co-occurred only **8 times**, despite the individual terms appearing much more often.

These results describe patterns in how topics are discussed in PubMed abstracts. Co-occurrence does not establish causation, direction of effect, or the clinical strength of an association.

## Workflow

### 1. Corpus Retrieval and Deduplication

- Used the `pubmedR` package to retrieve PubMed records programmatically.
- Created 11 overlapping search subsets covering general CF literature, organs, treatment, survival, CFTR modulators, aminoglycosides, dornase/enzyme therapies, adverse events, and mental health.
- Split the study period into 11 quinquennial (five-year) intervals from 1970 through 2024.
- Retrieved 94,597 records across the searches and deduplicated them by PubMed ID, leaving **41,222 unique articles**.

### 2. Text Preprocessing

- Converted abstracts to lowercase and removed punctuation for text analysis.
- Tokenised abstracts, removed standard and supplementary stopwords, and lemmatised tokens for vocabulary analysis.
- Standardised selected US spellings to UK spellings for the n-gram analysis.
- Preserved stopwords when generating bigrams and trigrams where needed to retain phrase context, then filtered common or irrelevant phrases during analysis.

### 3. Term Frequency and Temporal Analysis

- Counted unigrams, bigrams, and trigrams by five-year period.
- Used weighted log-odds to identify terms and phrases distinctive to each period, with filtering to reduce disease-name, systematic-review, and other methodological noise.
- Examined period-over-period changes in the use of selected treatment and adverse-event terms, considering term rates relative to the number of articles rather than relying only on raw counts.

### 4. Vocabulary Richness

Calculated vocabulary metrics for each period:

- **Gini coefficient:** measures how concentrated word frequencies are across the vocabulary.
- **Yule's K:** measures lexical repetitiveness; lower values indicate greater lexical diversity.
- **Hapax legomena:** counts words that occur only once.

I also used Spearman's rank correlation to examine the relationship between period-specific term distinctiveness and vocabulary concentration across unigram, bigram, and trigram analyses.

### 5. Co-occurrence Network Analysis

- Defined a treatment and adverse-event vocabulary using clinically relevant terms and a validated expansion from terms found in the corpus.
- Counted whether terms appeared together in the same abstract and constructed period-specific co-occurrence relationships.
- Used `quanteda` and `igraph`/`ggraph` to build and visualise the network.
- Examined CFTR modulator and mental-health term co-occurrence separately using unfiltered counts, because those pairs did not meet the threshold for inclusion in the main network.

## My Contributions

- Designed and implemented an R workflow to retrieve, clean, and analyse a large corpus of biomedical abstracts.
- Applied text-mining methods to compare terminology across five-year periods.
- Calculated and interpreted vocabulary richness metrics and term-distinctiveness statistics.
- Developed co-occurrence networks to explore how treatment, safety, infection, and mental-health terminology appeared together in the literature.
- Cross-checked the treatment-era terminology shift across multiple analyses and interpreted findings in light of the historical timeline of CFTR modulator therapy.

## Tools and Technologies

- **Language:** R
- **PubMed retrieval:** `pubmedR`
- **Data wrangling:** `tidyverse`, `dplyr`, `stringr`
- **Text mining:** `tidytext`, `textstem`, `tidylo`
- **Network analysis:** `quanteda`, `igraph`, `ggraph`
- **Visualisation and tables:** `ggplot2`, `gt`
- **Methods:** Corpus construction, deduplication, tokenisation, lemmatisation, n-gram analysis, weighted log-odds, vocabulary richness metrics, temporal analysis, co-occurrence networks, Spearman correlation

## Limitations

- The analysis uses PubMed titles and abstracts rather than full-text articles, limiting the context available for interpreting treatment and adverse-event relationships.
- Co-occurrence shows that terms appear in the same abstract; it cannot determine whether a treatment caused, reduced, or was unrelated to an adverse event.
- The network vocabulary was partly curated using prior clinical knowledge. Relevant relationships may be missed if the terminology was not included in the search vocabulary.
- English-language restrictions may underrepresent research published in other languages.
- The corpus includes different article types, including reviews, and some irrelevant phrases or text-processing artefacts remained despite filtering.
- PubMed is one source of biomedical literature and may not capture all relevant research.
- The number of publications increased substantially over the study period. Per-article rates and multiple analytical approaches were used to improve comparisons, but differences in corpus size still require care when interpreting summary statistics.

## Running the Analysis

The R scripts are numbered in the order of the analysis:

1. `01_setup.R` — defines search subsets, time periods, and package setup.
2. `02_corpus_retrieval.R` — retrieves PubMed records, deduplicates articles, and summarises corpus composition.
3. `03_preprocessing.R` — prepares abstracts and creates analysis-ready text.
4. `04_ngram_analysis.R` — calculates n-gram frequencies and weighted log-odds.
5. `05_vocab_richness_analysis.R` — calculates vocabulary metrics and related comparisons.
6. `06_network_analysis.R` — constructs and analyses treatment/adverse-event co-occurrence networks.

Install the packages listed in `01_setup.R` before running the scripts, and run them in order from the same working directory so intermediate `.rds` files can be found. Corpus retrieval requires an internet connection to access PubMed. The setup script includes a placeholder for an optional PubMed API key.

The scripts rely on intermediate data files generated during earlier steps. They are provided as the analysis workflow, but a clean end-to-end run should be tested in the target environment before treating the pipeline as fully reproducible. Retrieving the corpus again may also produce different results if PubMed records or indexing have changed since the dissertation data were collected.

## Author

**Kainoa Martinez**

MSc Health Data Science, University College London
