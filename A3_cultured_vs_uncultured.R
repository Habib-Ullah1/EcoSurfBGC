#!/usr/bin/env Rscript
library(dplyr)
library(readr)
library(tidyr)

BASE <- "/data/habib/EcoSurfBGC"
master <- read_csv(file.path(BASE, "results/integrated_master_table.csv"), show_col_types=FALSE)
density <- read_csv(file.path(BASE, "results/bgc_density.csv"), show_col_types=FALSE) %>%
  mutate(genome_id = as.character(genome_id))
novelty <- read_csv(file.path(BASE, "results/genome_novelty_summary.csv"), show_col_types=FALSE) %>%
  mutate(genome_id = as.character(genome_id))

df <- master %>%
  mutate(IMG.Genome.ID = as.character(IMG.Genome.ID)) %>%
  left_join(density, by=c("IMG.Genome.ID"="genome_id")) %>%
  left_join(novelty, by=c("IMG.Genome.ID"="genome_id"))

cat("Cultured column unique values:\n")
print(unique(df$Cultured))

df <- df %>%
  mutate(Culture_Status = case_when(
    grepl("Yes|yes", Cultured, ignore.case=TRUE) ~ "Cultured",
    grepl("No|no", Cultured, ignore.case=TRUE) ~ "Uncultured/MAG",
    TRUE ~ "Unknown"
  ))

cat("\n=== Culture Status Distribution ===\n")
print(df %>% count(Culture_Status) %>% arrange(desc(n)))

cat("\n=== BGC counts by culture status ===\n")
culture_bgc <- df %>%
  filter(Culture_Status != "Unknown") %>%
  group_by(Culture_Status) %>%
  summarise(
    n = n(),
    mean_total_bgcs    = round(mean(total_bgcs, na.rm=TRUE), 2),
    median_total_bgcs  = median(total_bgcs, na.rm=TRUE),
    mean_bgc_per_mb    = round(mean(bgc_per_mb, na.rm=TRUE), 2),
    median_bgc_per_mb  = median(bgc_per_mb, na.rm=TRUE),
    mean_genome_size   = round(mean(genome_size_mb, na.rm=TRUE), 2),
    mean_n_contigs     = round(mean(n_contigs, na.rm=TRUE), 0),
    mean_pct_novel     = round(mean(pct_novel, na.rm=TRUE), 1),
    pct_NRPS           = round(mean(NRPS > 0, na.rm=TRUE)*100, 1),
    pct_terpene        = round(mean(terpene > 0, na.rm=TRUE)*100, 1),
    pct_RiPP           = round(mean(`RiPP-like` > 0, na.rm=TRUE)*100, 1),
    .groups="drop"
  )
print(culture_bgc)

cat("\n=== Wilcoxon test: total BGCs cultured vs uncultured ===\n")
cultured_bgc   <- df$total_bgcs[df$Culture_Status == "Cultured"]
uncultured_bgc <- df$total_bgcs[df$Culture_Status == "Uncultured/MAG"]
if(length(cultured_bgc[!is.na(cultured_bgc)]) > 5 &
   length(uncultured_bgc[!is.na(uncultured_bgc)]) > 5) {
  wt <- wilcox.test(cultured_bgc, uncultured_bgc)
  cat(sprintf("Total BGCs - W=%.0f, p=%.3e\n", wt$statistic, wt$p.value))
  cat(sprintf("Median cultured: %.1f  |  Median uncultured: %.1f\n",
      median(cultured_bgc, na.rm=TRUE), median(uncultured_bgc, na.rm=TRUE)))
}

cat("\n=== Wilcoxon test: BGC density cultured vs uncultured ===\n")
cultured_d   <- df$bgc_per_mb[df$Culture_Status == "Cultured"]
uncultured_d <- df$bgc_per_mb[df$Culture_Status == "Uncultured/MAG"]
if(length(cultured_d[!is.na(cultured_d)]) > 5 &
   length(uncultured_d[!is.na(uncultured_d)]) > 5) {
  wt2 <- wilcox.test(cultured_d, uncultured_d)
  cat(sprintf("BGC/Mb - W=%.0f, p=%.3e\n", wt2$statistic, wt2$p.value))
  cat(sprintf("Median cultured: %.2f  |  Median uncultured: %.2f\n",
      median(cultured_d, na.rm=TRUE), median(uncultured_d, na.rm=TRUE)))
}

cat("\n=== Novelty by culture status ===\n")
novelty_culture <- df %>%
  filter(Culture_Status != "Unknown") %>%
  group_by(Culture_Status) %>%
  summarise(mean_pct_novel = round(mean(pct_novel, na.rm=TRUE),1),
            median_pct_novel = median(pct_novel, na.rm=TRUE),
            .groups="drop")
print(novelty_culture)

cat("\n=== Ecosystem breakdown by culture status ===\n")
eco_culture <- df %>%
  filter(Culture_Status != "Unknown") %>%
  group_by(Ecosystem_Category, Culture_Status) %>%
  summarise(n=n(), mean_bgc=round(mean(total_bgcs,na.rm=TRUE),1), .groups="drop") %>%
  pivot_wider(names_from=Culture_Status, values_from=c(n, mean_bgc), values_fill=0) %>%
  arrange(desc(mean_bgc_Cultured))
print(eco_culture, n=25)

write_csv(culture_bgc, file.path(BASE, "results/cultured_vs_uncultured_summary.csv"))
write_csv(df %>% select(IMG.Genome.ID, Culture_Status, genome_size_mb,
  n_contigs, n50_bp, bgc_per_mb, total_bgcs, pct_novel,
  Ecosystem_Category, C_Phylum),
  file.path(BASE, "results/genome_full_stats.csv"))
cat("\nA3 complete. Files saved.\n")
