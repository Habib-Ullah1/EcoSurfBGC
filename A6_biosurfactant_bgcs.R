#!/usr/bin/env Rscript
library(dplyr)
library(readr)
BASE <- "/data/habib/EcoSurfBGC"
master <- read_csv(file.path(BASE, "results/integrated_master_table_v2.csv"),
                   show_col_types=FALSE) %>%
  mutate(IMG.Genome.ID = as.character(IMG.Genome.ID))
novelty_bgc <- read_csv(file.path(BASE, "results/bgc_mibig_novelty.csv"),
                        show_col_types=FALSE) %>%
  mutate(genome_id = as.character(genome_id))

biosurf_keywords <- c("surfactin","iturin","fengycin","lichenysin",
  "rhamnolipid","serrawettin","arthrofactin","plipastatin",
  "lipopeptide","glycolipid","viscosin","massetolide","putisolvin")

novelty_bgc <- novelty_bgc %>%
  mutate(is_biosurf_related = grepl(
    paste(biosurf_keywords, collapse="|"),
    best_bgc_description, ignore.case=TRUE))

cat("Total BGCs:", nrow(novelty_bgc), "
")
cat("BGCs with biosurfactant hit:", sum(novelty_bgc$is_biosurf_related), "
")
cat("As pct of BGCs with hits:",
    round(sum(novelty_bgc$is_biosurf_related) /
          sum(novelty_bgc$has_mibig_hit==1) * 100, 1), "%
")

cat("
Biosurfactant BGCs by ecosystem:
")
print(novelty_bgc %>% filter(is_biosurf_related) %>%
  left_join(master %>% select(IMG.Genome.ID, Ecosystem_Category),
            by=c("genome_id"="IMG.Genome.ID")) %>%
  count(Ecosystem_Category) %>% arrange(desc(n)), n=20)

cat("
Fraction of BGC types with biosurfactant hits:
")
print(novelty_bgc %>% filter(has_mibig_hit==1) %>%
  group_by(best_bgc_type) %>%
  summarise(n_total=n(), n_biosurf=sum(is_biosurf_related),
    pct=round(sum(is_biosurf_related)/n()*100,1), .groups="drop") %>%
  filter(n_total>=10) %>% arrange(desc(pct)))

write_csv(novelty_bgc,
  file.path(BASE, "results/bgc_mibig_novelty_annotated.csv"))
cat("
Saved annotated file.
")
