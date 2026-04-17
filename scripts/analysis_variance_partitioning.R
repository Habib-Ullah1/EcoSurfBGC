# ============================================================
# Variance Partitioning Analysis
# Decomposing unique vs shared variance in BGC density
# explained by genome size, phylum, and ecosystem
# EcoSurfBGC — Ullah et al.
# ============================================================

library(vegan)
library(dplyr)
library(ggplot2)

set.seed(42)

# ── 1. Load data ─────────────────────────────────────────────
master <- read.csv('/data/habib/EcoSurfBGC/results/integrated_master_table_v2.csv',
                   stringsAsFactors = FALSE)

cat("Total genomes loaded:", nrow(master), "\n")
cat("Columns available:", paste(names(master), collapse=", "), "\n\n")

# ── 2. Prepare response variable ─────────────────────────────
# BGC density — the primary response variable
bgc_density <- master$bgc_per_mb

cat("BGC density summary:\n")
print(summary(bgc_density))
cat("Any NA:", sum(is.na(bgc_density)), "\n\n")

# ── 3. Prepare predictor matrices ────────────────────────────

# Predictor 1 — Genome size (continuous)
genome_size_matrix <- data.frame(
  genome_size_mb = master$genome_size_mb
)
cat("Genome size matrix dimensions:", dim(genome_size_matrix), "\n")

# Predictor 2 — Phylum (categorical — convert to dummy variables)
phylum_clean <- master$GTDB.Phylum
phylum_clean[is.na(phylum_clean) | phylum_clean == ""] <- "Unknown"
phylum_factor <- as.factor(phylum_clean)
phylum_matrix <- as.data.frame(model.matrix(~ phylum_factor - 1))
cat("Phylum matrix dimensions:", dim(phylum_matrix), "\n")
cat("Number of phyla:", ncol(phylum_matrix), "\n")

# Predictor 3 — Ecosystem (categorical — convert to dummy variables)
ecosystem_clean <- master$Ecosystem_Category
ecosystem_clean[is.na(ecosystem_clean) | ecosystem_clean == ""] <- "Unknown"
ecosystem_factor <- as.factor(ecosystem_clean)
ecosystem_matrix <- as.data.frame(model.matrix(~ ecosystem_factor - 1))
cat("Ecosystem matrix dimensions:", dim(ecosystem_matrix), "\n")
cat("Number of ecosystems:", ncol(ecosystem_matrix), "\n")

# Predictor 4 — Culture status (binary)
culture_matrix <- data.frame(
  culture_status = ifelse(master$Cultured == "Yes" | 
                          master$Cultured == "cultured" |
                          master$Cultured == "Cultured", 1, 0)
)
cat("Culture matrix dimensions:", dim(culture_matrix), "\n")
cat("Cultured genomes:", sum(culture_matrix$culture_status), "\n")
cat("MAG genomes:", sum(culture_matrix$culture_status == 0), "\n\n")

# ── 4. Remove rows with any NA ───────────────────────────────
complete_idx <- complete.cases(bgc_density,
                                genome_size_matrix,
                                phylum_matrix,
                                ecosystem_matrix,
                                culture_matrix)

cat("Complete cases:", sum(complete_idx), "of", length(complete_idx), "\n\n")

bgc_density_clean    <- bgc_density[complete_idx]
genome_size_clean    <- genome_size_matrix[complete_idx, , drop = FALSE]
phylum_clean_mat     <- phylum_matrix[complete_idx, ]
ecosystem_clean_mat  <- ecosystem_matrix[complete_idx, ]
culture_clean_mat    <- culture_matrix[complete_idx, , drop = FALSE]

# ── 5. Three-way variance partitioning ──────────────────────
# Primary partitioning: genome size vs phylum vs ecosystem
# Following Borcard et al. 1992 — the standard vegan approach

cat("============================================\n")
cat("VARIANCE PARTITIONING: Genome size | Phylum | Ecosystem\n")
cat("============================================\n")

vp3 <- varpart(bgc_density_clean,
               genome_size_clean,
               phylum_clean_mat,
               ecosystem_clean_mat)

cat("\nVariance partitioning results:\n")
print(vp3)
cat("\n")

# Extract individual fractions
# [a] = unique to genome size
# [b] = unique to phylum
# [c] = unique to ecosystem
# [d] = shared genome size + phylum
# [e] = shared genome size + ecosystem
# [f] = shared phylum + ecosystem
# [g] = shared all three
# [h] = unexplained

fractions <- vp3$part$indfract
cat("Individual fractions:\n")
print(fractions)

# ── 6. Four-way variance partitioning ───────────────────────
# Including culture status as fourth predictor

cat("\n============================================\n")
cat("VARIANCE PARTITIONING: Genome size | Phylum | Ecosystem | Culture\n")
cat("============================================\n")

vp4 <- varpart(bgc_density_clean,
               genome_size_clean,
               phylum_clean_mat,
               ecosystem_clean_mat,
               culture_clean_mat)

cat("\nFour-way variance partitioning results:\n")
print(vp4)

# ── 7. Test significance of each fraction ───────────────────
cat("\n============================================\n")
cat("SIGNIFICANCE TESTS — RDA permutation tests\n")
cat("============================================\n")

# Test unique contribution of genome size
cat("Unique genome size fraction:\n")
rda_genome <- rda(bgc_density_clean ~ genome_size_mb +
                    Condition(as.matrix(phylum_clean_mat)) +
                    Condition(as.matrix(ecosystem_clean_mat)),
                  data = cbind(genome_size_clean,
                               phylum_clean_mat,
                               ecosystem_clean_mat))
perm_genome <- anova(rda_genome, permutations = 999)
print(perm_genome)

# Test unique contribution of phylum
cat("\nUnique phylum fraction:\n")
rda_phylum <- rda(bgc_density_clean ~ as.matrix(phylum_clean_mat) +
                    Condition(genome_size_mb) +
                    Condition(as.matrix(ecosystem_clean_mat)),
                  data = cbind(genome_size_clean,
                               phylum_clean_mat,
                               ecosystem_clean_mat))
perm_phylum <- anova(rda_phylum, permutations = 999)
print(perm_phylum)

# Test unique contribution of ecosystem
cat("\nUnique ecosystem fraction:\n")
rda_ecosystem <- rda(bgc_density_clean ~ as.matrix(ecosystem_clean_mat) +
                       Condition(genome_size_mb) +
                       Condition(as.matrix(phylum_clean_mat)),
                     data = cbind(genome_size_clean,
                                  phylum_clean_mat,
                                  ecosystem_clean_mat))
perm_ecosystem <- anova(rda_ecosystem, permutations = 999)
print(perm_ecosystem)

# ── 8. Extract key numbers for manuscript ───────────────────
cat("\n============================================\n")
cat("KEY NUMBERS FOR MANUSCRIPT\n")
cat("============================================\n")

# Get adjusted R2 fractions from three-way partitioning
adj_fractions <- vp3$part$indfract

cat("Unique to genome size [a]:",
    round(adj_fractions$Adj.R.squared[1], 4), "\n")
cat("Unique to phylum [b]:",
    round(adj_fractions$Adj.R.squared[2], 4), "\n")
cat("Unique to ecosystem [c]:",
    round(adj_fractions$Adj.R.squared[3], 4), "\n")
cat("Shared genome size + phylum [d]:",
    round(adj_fractions$Adj.R.squared[4], 4), "\n")
cat("Shared genome size + ecosystem [e]:",
    round(adj_fractions$Adj.R.squared[5], 4), "\n")
cat("Shared phylum + ecosystem [f]:",
    round(adj_fractions$Adj.R.squared[6], 4), "\n")
cat("Shared all three [g]:",
    round(adj_fractions$Adj.R.squared[7], 4), "\n")
cat("Unexplained [h]:",
    round(adj_fractions$Adj.R.squared[8], 4), "\n")
cat("Total explained:",
    round(sum(adj_fractions$Adj.R.squared[1:7]), 4), "\n")

# ── 9. Save results ─────────────────────────────────────────
write.csv(as.data.frame(vp3$part$indfract),
          '/data/habib/EcoSurfBGC/results/variance_partitioning_3way.csv',
          row.names = TRUE)

write.csv(as.data.frame(vp4$part$indfract),
          '/data/habib/EcoSurfBGC/results/variance_partitioning_4way.csv',
          row.names = TRUE)

# ── 10. Plot variance partitioning Venn diagram ─────────────
pdf('/data/habib/EcoSurfBGC/figures/variance_partitioning_venn.pdf',
    width = 8, height = 7)
plot(vp3,
     digits    = 2,
     bg        = c("#E1F5EE", "#EEEDFE", "#FAEEDA"),
     col       = c("#0F6E56", "#534AB7", "#BA7517"),
     Xnames    = c("Genome size", "Phylum", "Ecosystem"),
     main      = "Variance partitioning of BGC density\n(adjusted R²)")
dev.off()

png('/data/habib/EcoSurfBGC/figures/variance_partitioning_venn.png',
    width = 2400, height = 2100, res = 300)
plot(vp3,
     digits    = 2,
     bg        = c("#E1F5EE", "#EEEDFE", "#FAEEDA"),
     col       = c("#0F6E56", "#534AB7", "#BA7517"),
     Xnames    = c("Genome size", "Phylum", "Ecosystem"),
     main      = "Variance partitioning of BGC density\n(adjusted R²)")
dev.off()

cat("\nAll results saved.\n")
cat("Analysis complete.\n")

