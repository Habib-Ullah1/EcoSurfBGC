# ============================================================
# Phylogenetic Signal Analysis — Blomberg's K
# BGC density across GTDB phyla using backbone constraint tree
# EcoSurfBGC — Ullah et al.
# ============================================================

library(ape)
library(phytools)
library(vegan)
library(dplyr)
library(ggplot2)

set.seed(42)

# ── 1. Load data ────────────────────────────────────────────
master <- read.csv('/data/habib/EcoSurfBGC/results/integrated_master_table_v2.csv',
                   stringsAsFactors = FALSE)

cat("Total genomes loaded:", nrow(master), "\n")

# Clean GTDB phylum — remove empty or NA
master <- master %>%
  filter(!is.na(GTDB.Phylum) & GTDB.Phylum != "" & GTDB.Phylum != "unclassified")

cat("Genomes with valid GTDB Phylum:", nrow(master), "\n")

# ── 2. Calculate phylum-level mean BGC density ──────────────
phylum_stats <- master %>%
  group_by(GTDB.Phylum) %>%
  summarise(
    n_genomes       = n(),
    mean_bgc_density = mean(bgc_per_mb, na.rm = TRUE),
    sd_bgc_density   = sd(bgc_per_mb, na.rm = TRUE),
    se_bgc_density   = sd(bgc_per_mb, na.rm = TRUE) / sqrt(n()),
    median_bgc_density = median(bgc_per_mb, na.rm = TRUE),
    mean_genome_size = mean(genome_size_mb, na.rm = TRUE),
    mean_total_bgcs  = mean(total_bgcs, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  filter(n_genomes >= 3) %>%   # minimum 3 genomes per phylum for reliability
  arrange(desc(mean_bgc_density))

cat("\nPhyla with >=3 genomes:", nrow(phylum_stats), "\n")
print(phylum_stats, n = 40)

# Save phylum stats
write.csv(phylum_stats,
          '/data/habib/EcoSurfBGC/results/phylum_bgc_density_stats.csv',
          row.names = FALSE)

# ── 3. Build GTDB backbone constraint tree ──────────────────
# GTDB r214 reference topology for bacterial/archaeal phyla
# Based on Parks et al. 2022 — standardized bacterial taxonomy

# Define phyla present in dataset
phyla_in_data <- phylum_stats$GTDB.Phylum
cat("\nPhyla included in phylogenetic analysis:\n")
print(phyla_in_data)

# GTDB r214 backbone topology (Newick format)
# Reference: Parks et al. 2022 Nature Microbiology
# Topology based on concatenated marker gene phylogeny
gtdb_backbone_newick <- "((((((Pseudomonadota,Bdellovibrionota,Myxococcota),
(Desulfobacterota,Campylobacterota,Campylobacterota_A,Aquificota,Chrysiogenota,
Fusobacteriota,Synergistota,Thermotogota,Fibrobacterota)),
(Spirochaetota,Marinisomatota)),
((Actinomycetota,Chloroflexota,Deinococcota,Armatimonadota),
(Acidobacteriota,Verrucomicrobiota,Planctomycetota,Patescibacteria),
(Bacteroidota,Cyanobacteriota))),
((Bacillota,Bacillota_A,Bacillota_B,Bacillota_C,Bacillota_D,
Bacillota_E,Bacillota_F,Bacillota_I),
(Chlamydiota))),
(Halobacteriota,Thermoproteota,Nanoarchaeota,Micrarchaeota,
Methanobacteriota,Methanobacteriota_A,Methanobacteriota_B,
Thermoplasmatota));"

# Parse the backbone tree
backbone_tree <- read.tree(text = gtdb_backbone_newick)
cat("\nBackbone tree tip labels:\n")
print(backbone_tree$tip.label)

# Keep only tips present in our dataset
tips_to_keep <- backbone_tree$tip.label[backbone_tree$tip.label %in% phyla_in_data]
tips_to_drop  <- backbone_tree$tip.label[!backbone_tree$tip.label %in% phyla_in_data]

cat("\nPhyla in tree AND dataset:", length(tips_to_keep), "\n")
cat("Phyla in tree but not dataset:", length(tips_to_drop), "\n")
print(tips_to_drop)

# Prune tree to only phyla in dataset
pruned_tree <- drop.tip(backbone_tree, tips_to_drop)
cat("\nPruned tree tips:", length(pruned_tree$tip.label), "\n")
cat("Tree is rooted:", is.rooted(pruned_tree), "\n")

# Make tree ultrametric (required for Blomberg's K)
# Use equal branch lengths since divergence times unknown at phylum level
pruned_tree$edge.length <- rep(1, nrow(pruned_tree$edge))
ultrametric_tree <- compute.brlen(pruned_tree, method = "Grafen")

cat("Tree is ultrametric:", is.ultrametric(ultrametric_tree), "\n")

# ── 4. Prepare trait vector ──────────────────────────────────
# Match phylum stats to tree tips
trait_df <- phylum_stats %>%
  filter(GTDB.Phylum %in% ultrametric_tree$tip.label)

bgc_density_vector <- setNames(trait_df$mean_bgc_density,
                                trait_df$GTDB.Phylum)

genome_size_vector  <- setNames(trait_df$mean_genome_size,
                                 trait_df$GTDB.Phylum)

cat("\nTrait vector length:", length(bgc_density_vector), "\n")
cat("Tree tips:", length(ultrametric_tree$tip.label), "\n")

# Verify alignment
cat("All tree tips in trait vector:",
    all(ultrametric_tree$tip.label %in% names(bgc_density_vector)), "\n")

# ── 5. Blomberg's K — BGC density ───────────────────────────
cat("\n============================================\n")
cat("BLOMBERG'S K — BGC density\n")
cat("============================================\n")

K_bgc <- phylosig(ultrametric_tree,
                  bgc_density_vector[ultrametric_tree$tip.label],
                  method = "K",
                  test   = TRUE,
                  nsim   = 10000)

cat("K statistic:", round(K_bgc$K, 4), "\n")
cat("P-value:", round(K_bgc$P, 4), "\n")

if (K_bgc$K > 1) {
  cat("Interpretation: STRONG phylogenetic signal\n")
  cat("BGC density is more conserved than expected under Brownian motion\n")
  cat("Secondary metabolic investment is an evolutionarily conserved trait\n")
} else if (K_bgc$K > 0.5 & K_bgc$P < 0.05) {
  cat("Interpretation: MODERATE phylogenetic signal\n")
  cat("BGC density shows partial phylogenetic conservatism\n")
} else if (K_bgc$P < 0.05) {
  cat("Interpretation: WEAK but significant phylogenetic signal\n")
} else {
  cat("Interpretation: NO significant phylogenetic signal\n")
  cat("BGC density is ecologically plastic — not phylogenetically conserved\n")
}

# ── 6. Blomberg's K — Genome size ───────────────────────────
cat("\n============================================\n")
cat("BLOMBERG'S K — Genome size\n")
cat("============================================\n")

K_genome <- phylosig(ultrametric_tree,
                     genome_size_vector[ultrametric_tree$tip.label],
                     method = "K",
                     test   = TRUE,
                     nsim   = 10000)

cat("K statistic:", round(K_genome$K, 4), "\n")
cat("P-value:", round(K_genome$P, 4), "\n")

# ── 7. Lambda (Pagel's lambda) — alternative signal test ────
cat("\n============================================\n")
cat("PAGEL'S LAMBDA — BGC density\n")
cat("============================================\n")

lambda_bgc <- phylosig(ultrametric_tree,
                       bgc_density_vector[ultrametric_tree$tip.label],
                       method = "lambda",
                       test   = TRUE)

cat("Lambda:", round(lambda_bgc$lambda, 4), "\n")
cat("P-value (vs lambda=0):", round(lambda_bgc$P, 4), "\n")
cat("Log-likelihood:", round(lambda_bgc$logL, 4), "\n")

# ── 8. Phylogenetic ANOVA — BGC density by phylum ───────────
cat("\n============================================\n")
cat("PHYLOGENETIC ANOVA — does phylum predict BGC density?\n")
cat("============================================\n")

# Test whether phylum explains BGC density after accounting for phylogeny
# Using phylogenetic generalized least squares (PGLS)
library(nlme)
library(ape)

# Create data frame aligned to tree
pgls_data <- data.frame(
  bgc_density  = bgc_density_vector[ultrametric_tree$tip.label],
  genome_size  = genome_size_vector[ultrametric_tree$tip.label],
  row.names    = ultrametric_tree$tip.label
)

# Fit Brownian motion correlation structure
bm_corr <- corBrownian(phy = ultrametric_tree)

pgls_model <- gls(bgc_density ~ genome_size,
                  data       = pgls_data,
                  correlation = bm_corr)

cat("PGLS: BGC density ~ genome size (phylogenetically corrected)\n")
print(summary(pgls_model))

# ── 9. Save results ─────────────────────────────────────────
results_list <- list(
  K_bgc_density   = data.frame(K = K_bgc$K,   P = K_bgc$P,
                                 interpretation = ifelse(K_bgc$P < 0.05,
                                   "significant", "not significant")),
  K_genome_size   = data.frame(K = K_genome$K, P = K_genome$P,
                                 interpretation = ifelse(K_genome$P < 0.05,
                                   "significant", "not significant")),
  lambda_bgc      = data.frame(lambda = lambda_bgc$lambda,
                                P      = lambda_bgc$P)
)

cat("\n============================================\n")
cat("SUMMARY OF PHYLOGENETIC SIGNAL RESULTS\n")
cat("============================================\n")
cat("BGC density K:", round(K_bgc$K, 4),
    "| P:", round(K_bgc$P, 4), "\n")
cat("Genome size K:", round(K_genome$K, 4),
    "| P:", round(K_genome$P, 4), "\n")
cat("BGC density lambda:", round(lambda_bgc$lambda, 4),
    "| P:", round(lambda_bgc$P, 4), "\n")

# Save pruned tree for figure
write.tree(ultrametric_tree,
           '/data/habib/EcoSurfBGC/results/phylum_pruned_tree.nwk')

# Save results table
results_df <- data.frame(
  Test      = c("Blomberg K - BGC density",
                "Blomberg K - Genome size",
                "Pagel lambda - BGC density"),
  Statistic = c(K_bgc$K, K_genome$K, lambda_bgc$lambda),
  P_value   = c(K_bgc$P, K_genome$P, lambda_bgc$P),
  Significant = c(K_bgc$P < 0.05, K_genome$P < 0.05, lambda_bgc$P < 0.05)
)

write.csv(results_df,
          '/data/habib/EcoSurfBGC/results/phylogenetic_signal_results.csv',
          row.names = FALSE)

cat("\nResults saved to:\n")
cat("  /data/habib/EcoSurfBGC/results/phylogenetic_signal_results.csv\n")
cat("  /data/habib/EcoSurfBGC/results/phylum_bgc_density_stats.csv\n")
cat("  /data/habib/EcoSurfBGC/results/phylum_pruned_tree.nwk\n")
cat("\nAnalysis complete.\n")

