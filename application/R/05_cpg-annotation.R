# load required packages and scripts
library(here)
library(tidyverse)
library(conflicted)
library(knitr)
library(kableExtra)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
conflict_prefer("filter", "dplyr")
conflict_prefer("slice", "dplyr")

# load the results
full_data <- read_rds(here("results", "full-dataset-biotmle.Rds"))

# Add CpG Annotation ###########################################################
# load and prep annotation data
anno <- IlluminaHumanMethylation450kanno.ilmn12.hg19 %>%
  getAnnotation() %>%
  as_tibble(rownames = "cpg")

# create an annotation table
annotated_results <- full_data %>%
  filter(adj_pval <= 0.05) %>%
  mutate(cpg = ID) %>%
  left_join(anno, by = "cpg") %>%
  select(
    cpg, chr, pos, t, P.Value, adj_pval, UCSC_RefGene_Name, Relation_to_Island
  ) %>%
  arrange(adj_pval)

# save the results
write_csv(
  annotated_results,
  file = here("results", "full-dataset-annotated-results.csv")
)

# save a table of results for the supplement
table_results_anno <- annotated_results %>%
  slice(1:50) %>%
  select(-t, -P.Value) %>%
  mutate(adj_pval = as.character(signif(adj_pval, 3))) %>%
  kable(col.names = c("CpG", "Chromosome", "Position",
                      "Adj. P-value", "Gene Name (UCSC)",
                      "CpG Island Relation"),
        format = "latex",
        booktabs = TRUE,
        caption = "Top 50 CpG sites (ranked by adjusted p-value) tagged by our nonparametric differential methlation analysis.",
        label = "anno_top50",
        digits = 50
  ) %>%
  kable_styling()
write(table_results_anno, here("tables", "cpg_anno_results.tex"))
