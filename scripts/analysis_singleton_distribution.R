# ============================================================
# Singleton BGC Phylogenetic Distribution Analysis
# Which phyla harbor the most singleton BGCs per genome?
# EcoSurfBGC — Ullah et al.
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)

set.seed(42)

# ── 1. Load data ─────────────────────────────────────────────
master <- read.csv('/data/habib/EcoSurfBGC/results/integrated_master_table_v2.csv',
                   stringsAsFactors = FALSE)
gcf    <- read.csv('/data/habib/EcoSurfBGC/supplementary/Table_S5_bigscape_gcf_summary.csv',
                   stringsAsFactors = FALSE)

cat("Master table rows:", nrow(master), "\n")
cat("GCF table rows:", nrow(gcf), "\n")
cat("Total singleton GCFs:", sum(gcf$Is_singleton), "\n")
cat("Total non-singleton GCFs:", sum(!gcf$Is_singleton), "\n")
cat("Singleton fraction:", round(mean(gcf$Is_singleton)*100, 1), "%\n\n")

# ── 2. Extract singleton BGC genome IDs ──────────────────────
# Each singleton GCF has exactly one BGC and one genome
singleton_gcfs <- gcf %>% filter(Is_singleton == TRUE)
cat("Singleton GCFs:", nrow(singleton_gcfs), "\n")

# Extract genome IDs from singleton GCFs
singleton_genome_ids <- singleton_gcfs %>%
  select(GCF_id, Genome_IDs, BiGSCAPE_classes) %>%
  mutate(Genome_ID = trimws(Genome_IDs)) %>%
  select(GCF_id, Genome_ID, BiGSCAPE_classes)

cat("Singleton BGCs extracted:", nrow(singleton_genome_ids), "\n\n")

# ── 3. Link to master table via genome ID ────────────────────
# Clean genome IDs in master table
master$IMG.Genome.ID <- as.character(master$IMG.Genome.ID)
singleton_genome_ids$Genome_ID <- as.character(singleton_genome_ids$Genome_ID)

# Join singleton counts per genome to master
singleton_counts <- singleton_genome_ids %>%
  group_by(Genome_ID) %>%
  summarise(
    n_singleton_bgcs = n(),
    singleton_classes = paste(unique(BiGSCAPE_classes), collapse=";"),
    .groups = 'drop'
  )

cat("Genomes with at least one singleton BGC:",
    nrow(singleton_counts), "\n")

# Merge with master
master_with_singletons <- master %>%
  left_join(singleton_counts,
            by = c("IMG.Genome.ID" = "Genome_ID")) %>%
  mutate(
    n_singleton_bgcs = ifelse(is.na(n_singleton_bgcs), 0, n_singleton_bgcs),
    has_singleton    = n_singleton_bgcs > 0,
    pct_singleton    = ifelse(total_bgcs > 0,
                              n_singleton_bgcs / total_bgcs * 100, 0)
  )

cat("Genomes with singleton BGCs:",
    sum(master_with_singletons$has_singleton), "\n")
cat("Mean singleton BGCs per genome:",
    round(mean(master_with_singletons$n_singleton_bgcs), 3), "\n\n")

# ── 4. Phylum-level singleton analysis ───────────────────────
phylum_singleton <- master_with_singletons %>%
  filter(!is.na(GTDB.Phylum) & GTDB.Phylum != "") %>%
  group_by(GTDB.Phylum) %>%
  summarise(
    n_genomes             = n(),
    total_bgcs            = sum(total_bgcs, na.rm = TRUE),
    total_singletons      = sum(n_singleton_bgcs, na.rm = TRUE),
    mean_singletons_per_genome = mean(n_singleton_bgcs, na.rm = TRUE),
    mean_bgc_density      = mean(bgc_per_mb, na.rm = TRUE),
    pct_genomes_with_singleton = mean(has_singleton, na.rm = TRUE) * 100,
    mean_pct_singleton    = mean(pct_singleton, na.rm = TRUE),
    singleton_rate        = ifelse(total_bgcs > 0,
                                   total_singletons / total_bgcs * 100, 0),
    .groups = 'drop'
  ) %>%
  filter(n_genomes >= 5) %>%
  arrange(desc(mean_singletons_per_genome))

cat("============================================\n")
cat("PHYLUM-LEVEL SINGLETON DISTRIBUTION\n")
cat("============================================\n")
print(phylum_singleton, n = 30)

# ── 5. Ecosystem-level singleton analysis ────────────────────
ecosystem_singleton <- master_with_singletons %>%
  filter(!is.na(Ecosystem_Category) & Ecosystem_Category != "") %>%
  group_by(Ecosystem_Category) %>%
  summarise(
    n_genomes             = n(),
    total_bgcs            = sum(total_bgcs, na.rm = TRUE),
    total_singletons      = sum(n_singleton_bgcs, na.rm = TRUE),
    mean_singletons_per_genome = mean(n_singleton_bgcs, na.rm = TRUE),
    singleton_rate        = ifelse(total_bgcs > 0,
                                   total_singletons / total_bgcs * 100, 0),
    .groups = 'drop'
  ) %>%
  arrange(desc(singleton_rate))

cat("\n============================================\n")
cat("ECOSYSTEM-LEVEL SINGLETON DISTRIBUTION\n")
cat("============================================\n")
print(ecosystem_singleton, n = 20)

# ── 6. Key numbers for manuscript ────────────────────────────
cat("\n============================================\n")
cat("KEY NUMBERS FOR MANUSCRIPT\n")
cat("============================================\n")

# Top 3 phyla by singleton BGCs per genome
top3_phyla <- phylum_singleton %>%
  slice_max(mean_singletons_per_genome, n = 3)
cat("Top 3 phyla by mean singleton BGCs per genome:\n")
print(top3_phyla %>% select(GTDB.Phylum, n_genomes,
                              mean_singletons_per_genome,
                              singleton_rate))

# Bottom 3 phyla
bottom3_phyla <- phylum_singleton %>%
  slice_min(mean_singletons_per_genome, n = 3)
cat("\nBottom 3 phyla by mean singleton BGCs per genome:\n")
print(bottom3_phyla %>% select(GTDB.Phylum, n_genomes,
                                mean_singletons_per_genome,
                                singleton_rate))

# Global singleton rate
global_singleton_rate <- sum(gcf$Is_singleton) / nrow(gcf) * 100
cat("\nGlobal singleton rate:", round(global_singleton_rate, 1), "%\n")

# Phyla with above-average singleton rates
above_avg <- phylum_singleton %>%
  filter(singleton_rate > global_singleton_rate) %>%
  arrange(desc(singleton_rate))
cat("\nPhyla with above-average singleton rates:\n")
print(above_avg %>% select(GTDB.Phylum, n_genomes,
                             singleton_rate,
                             mean_singletons_per_genome))

# Top ecosystems by singleton rate
cat("\nTop 5 ecosystems by singleton rate:\n")
print(head(ecosystem_singleton %>%
           select(Ecosystem_Category, n_genomes,
                  singleton_rate, total_singletons), 5))

# ── 7. Save results ──────────────────────────────────────────
write.csv(phylum_singleton,
          '/data/habib/EcoSurfBGC/results/singleton_phylum_distribution.csv',
          row.names = FALSE)

write.csv(ecosystem_singleton,
          '/data/habib/EcoSurfBGC/results/singleton_ecosystem_distribution.csv',
          row.names = FALSE)

# ── 8. Generate figure ───────────────────────────────────────
# Bar chart — phyla ranked by mean singleton BGCs per genome
plot_data <- phylum_singleton %>%
  filter(n_genomes >= 5) %>%
  mutate(GTDB.Phylum = gsub("_", " ", GTDB.Phylum)) %>%
  arrange(mean_singletons_per_genome) %>%
  mutate(GTDB.Phylum = factor(GTDB.Phylum,
                               levels = unique(GTDB.Phylum)))

p <- ggplot(plot_data,
            aes(x = mean_singletons_per_genome,
                y = GTDB.Phylum,
                fill = singleton_rate)) +
  geom_bar(stat = "identity", width = 0.75) +
  geom_text(aes(label = paste0("n=", n_genomes)),
            hjust = -0.1, size = 2.8,
            color = "#2C2C2A") +
  scale_fill_gradient(low  = "#E6F1FB",
                      high = "#042C53",
                      name = "Singleton\nrate (%)") +
  labs(x     = "Mean singleton BGCs per genome",
       y     = NULL,
       title = "Singleton BGC distribution across phyla") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.y        = element_text(size = 9),
    plot.title         = element_text(size = 11, face = "bold"),
    legend.position    = "right"
  ) +
  expand_limits(x = max(plot_data$mean_singletons_per_genome) * 1.15)

ggsave('/data/habib/EcoSurfBGC/figures/FigureS6_singleton_distribution.pdf',
       p, width = 9, height = 8)
ggsave('/data/habib/EcoSurfBGC/figures/FigureS6_singleton_distribution.png',
       p, width = 9, height = 8, dpi = 300)

cat("\nAll results and figures saved.\n")
cat("Analysis complete.\n")

