#!/usr/bin/env Rscript
library(dplyr)
library(readr)

BASE <- "/data/habib/EcoSurfBGC"
master  <- read_csv(file.path(BASE, "results/integrated_master_table.csv"), show_col_types=FALSE)
density <- read_csv(file.path(BASE, "results/bgc_density.csv"), show_col_types=FALSE) %>%
  mutate(genome_id = as.character(genome_id))

df <- master %>%
  mutate(IMG.Genome.ID = as.character(IMG.Genome.ID)) %>%
  left_join(density, by=c("IMG.Genome.ID"="genome_id")) %>%
  mutate(bgc_category = case_when(
    total_bgcs == 0 ~ "Zero",
    total_bgcs <= 3 ~ "Low (1-3)",
    total_bgcs <= 7 ~ "Medium (4-7)",
    TRUE            ~ "High (8+)"
  ))

zero <- df %>% filter(total_bgcs == 0)
cat("=== Zero-BGC genomes: overview ===\n")
cat("Total zero-BGC genomes:", nrow(zero), "\n")
cat(sprintf("Mean completeness: %.1f%%\n", mean(zero$CheckM2.Completeness, na.rm=TRUE)))
cat(sprintf("Mean contamination: %.1f%%\n", mean(zero$CheckM2.Contamination, na.rm=TRUE)))
cat(sprintf("Mean genome size: %.2f Mb\n", mean(zero$genome_size_mb, na.rm=TRUE)))
cat(sprintf("Mean contigs: %.0f\n", mean(zero$n_contigs, na.rm=TRUE)))
cat(sprintf("Mean N50: %.0f bp\n", mean(zero$n50_bp, na.rm=TRUE)))

cat("\n=== Comparison across BGC categories ===\n")
comparison <- df %>%
  filter(!is.na(genome_size_mb)) %>%
  group_by(bgc_category) %>%
  summarise(
    n                  = n(),
    mean_completeness  = round(mean(CheckM2.Completeness, na.rm=TRUE),1),
    mean_size_mb       = round(mean(genome_size_mb, na.rm=TRUE),2),
    mean_contigs       = round(mean(n_contigs, na.rm=TRUE),0),
    mean_n50_kb        = round(mean(n50_bp, na.rm=TRUE)/1000,1),
    mean_known_genes   = round(mean(Biosurf_Total, na.rm=TRUE),2),
    .groups="drop"
  ) %>% arrange(bgc_category)
print(comparison)

cat("\n=== Zero-BGC by ecosystem ===\n")
print(zero %>% count(Ecosystem_Category) %>% arrange(desc(n)))

cat("\n=== Zero-BGC by phylum ===\n")
print(zero %>% count(C_Phylum) %>% arrange(desc(n)) %>% head(15))

cat("\n=== Are zero-BGC genomes more fragmented? (Wilcoxon) ===\n")
wt_c <- wilcox.test(df$n_contigs[df$total_bgcs==0],
                    df$n_contigs[df$total_bgcs>0], na.action=na.omit)
cat(sprintf("Contig count - W=%.0f, p=%.3e\n", wt_c$statistic, wt_c$p.value))
cat(sprintf("Median contigs zero: %.0f vs non-zero: %.0f\n",
    median(df$n_contigs[df$total_bgcs==0], na.rm=TRUE),
    median(df$n_contigs[df$total_bgcs>0],  na.rm=TRUE)))

cat("\n=== Are zero-BGC genomes smaller? (Wilcoxon) ===\n")
wt_s <- wilcox.test(df$genome_size_mb[df$total_bgcs==0],
                    df$genome_size_mb[df$total_bgcs>0], na.action=na.omit)
cat(sprintf("Genome size - W=%.0f, p=%.3e\n", wt_s$statistic, wt_s$p.value))
cat(sprintf("Median size zero: %.2f Mb vs non-zero: %.2f Mb\n",
    median(df$genome_size_mb[df$total_bgcs==0], na.rm=TRUE),
    median(df$genome_size_mb[df$total_bgcs>0],  na.rm=TRUE)))

cat("\n=== BGC detection by genome size bin ===\n")
size_bgc <- df %>%
  filter(!is.na(genome_size_mb)) %>%
  mutate(size_bin = cut(genome_size_mb,
    breaks=c(0,0.5,1,1.5,2,3,4,5,10,Inf),
    labels=c("<0.5","0.5-1","1-1.5","1.5-2","2-3","3-4","4-5","5-10",">10"))) %>%
  group_by(size_bin) %>%
  summarise(n=n(), pct_zero=round(mean(total_bgcs==0)*100,1),
            mean_bgcs=round(mean(total_bgcs),1), .groups="drop")
print(size_bgc)

write_csv(zero %>% select(IMG.Genome.ID, C_Phylum, C_Genus,
  Ecosystem_Category, total_bgcs, CheckM2.Completeness,
  CheckM2.Contamination, genome_size_mb, n_contigs, n50_bp, Biosurf_Total),
  file.path(BASE, "results/zero_bgc_genomes.csv"))
cat("\nA5 complete.\n")
