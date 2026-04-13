suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr)
  library(ggplot2); library(ComplexHeatmap)
  library(circlize); library(RColorBrewer)
  library(forcats); library(stringr); library(scales)
  library(tibble); library(dunn.test)
})

BASE <- "/data/habib/EcoSurfBGC"
SUPP <- file.path(BASE, "supplementary")
dir.create(SUPP, showWarnings=FALSE)

cat("Loading datasets...\n")

master <- read_csv(file.path(BASE,"results/integrated_master_table_v3.csv"),
                   show_col_types=FALSE) %>%
  mutate(
    IMG.Genome.ID = as.character(IMG.Genome.ID),
    Genome_source = ifelse(grepl("^[0-9]+$", IMG.Genome.ID), "JGI", "NCBI")
  )

bgc_detailed <- read_csv(file.path(BASE,"results/bgc_detailed_list.csv"),
                          show_col_types=FALSE) %>%
  mutate(genome_id=as.character(genome_id))

gcf_genome <- read_csv(file.path(BASE,"results/genome_gcf_summary.csv"),
                        show_col_types=FALSE) %>%
  mutate(genome_id=as.character(genome_id))

gcf_bgcs <- read_csv(file.path(BASE,"results/bgc_gcf_assignments.csv"),
                      show_col_types=FALSE) %>%
  mutate(genome_id=as.character(genome_id))

cat("All datasets loaded.\n\n")

# Identify BGC product type columns
# These are numeric columns between total_bgcs and Culture_Status
all_cols <- names(master)
start_idx <- which(all_cols == "total_bgcs") + 1
end_idx   <- which(all_cols == "Culture_Status") - 1
bgc_cols  <- all_cols[start_idx:end_idx]
bgc_cols  <- bgc_cols[sapply(bgc_cols,
               function(x) is.numeric(master[[x]]))]
cat(sprintf("BGC product type columns identified: %d\n\n", length(bgc_cols)))

#=============================================================
# SUPPLEMENTARY TABLE S1
# Complete genome metadata for all 2,057 genomes
#=============================================================
cat("Generating Supplementary Table S1...\n")

S1 <- master %>%
  select(
    Genome_ID                = IMG.Genome.ID,
    Genome_source            = Genome_source,
    Ecosystem                = Ecosystem_Category,
    Phylum                   = C_Phylum,
    Class                    = C_Class,
    Order                    = C_Order,
    Family                   = C_Family,
    Genus                    = C_Genus,
    Culture_Status           = Culture_Status,
    Completeness_pct         = CheckM2.Completeness,
    Contamination_pct        = CheckM2.Contamination,
    Genome_size_Mb           = genome_size_mb,
    N_contigs                = n_contigs,
    N50_bp                   = n50_bp,
    Total_BGCs               = total_bgcs,
    BGC_density_per_Mb       = bgc_per_mb,
    N_GCFs                   = n_gcfs,
    Pct_singleton_BGCs       = pct_singleton,
    Pct_novel_knownCB        = pct_novel,
    Pct_MIBiG_conn_BiGSCAPE = pct_mibig_conn,
    N_BGC_classes            = n_bgc_classes,
    Known_biosurf_genes      = Biosurf_Total
  ) %>%
  arrange(Ecosystem, Phylum, Genome_ID)

write_csv(S1,
  file.path(SUPP, "Table_S1_genome_metadata.csv"), na="NA")
cat(sprintf("  S1 saved: %d genomes x %d columns\n\n",
            nrow(S1), ncol(S1)))

#=============================================================
# SUPPLEMENTARY TABLE S2
# Complete BGC class matrix (2,057 genomes x all product types)
#=============================================================
cat("Generating Supplementary Table S2...\n")

S2 <- master %>%
  select(
    Genome_ID      = IMG.Genome.ID,
    Genome_source  = Genome_source,
    Ecosystem      = Ecosystem_Category,
    Phylum         = C_Phylum,
    Culture_Status = Culture_Status,
    all_of(bgc_cols)
  ) %>%
  arrange(Ecosystem, Phylum)

write_csv(S2,
  file.path(SUPP, "Table_S2_bgc_class_matrix.csv"), na="0")
cat(sprintf("  S2 saved: %d genomes x %d BGC product types\n\n",
            nrow(S2), length(bgc_cols)))

#=============================================================
# SUPPLEMENTARY TABLE S3
# Complete Dunn pairwise comparisons for BGC density
#=============================================================
cat("Generating Supplementary Table S3...\n")

df_dunn <- master %>%
  filter(!is.na(bgc_per_mb), !is.na(Ecosystem_Category))

dunn_result <- dunn.test(
  df_dunn$bgc_per_mb,
  df_dunn$Ecosystem_Category,
  method = "bh",
  altp   = TRUE
)

S3 <- data.frame(
  Ecosystem_1        = str_split_fixed(dunn_result$comparisons," - ",2)[,1],
  Ecosystem_2        = str_split_fixed(dunn_result$comparisons," - ",2)[,2],
  Z_statistic        = round(dunn_result$Z, 4),
  P_raw              = signif(dunn_result$altP, 4),
  P_adj_BH           = signif(dunn_result$altP.adjusted, 4),
  Significant_FDR005 = dunn_result$altP.adjusted < 0.05
) %>%
  arrange(P_adj_BH)

write_csv(S3,
  file.path(SUPP, "Table_S3_dunn_pairwise_complete.csv"))
cat(sprintf("  S3 saved: %d pairs (%d significant FDR<0.05)\n\n",
            nrow(S3), sum(S3$Significant_FDR005)))

#=============================================================
# SUPPLEMENTARY TABLE S4
# Complete BGC class enrichment results
#=============================================================
cat("Generating Supplementary Table S4...\n")

key_classes <- c("NRPS","NRPS-like","T1PKS","T3PKS","terpene",
                 "RiPP-like","betalactone","arylpolyene",
                 "NI-siderophore","ectoine","hserlactone",
                 "lassopeptide","thiopeptide","ranthipeptide","NAPAA")
key_classes <- key_classes[key_classes %in% names(master)]

ecosystems <- unique(master$Ecosystem_Category)
ecosystems <- ecosystems[!is.na(ecosystems)]

S4_rows <- list()
for(eco in ecosystems) {
  in_eco  <- !is.na(master$Ecosystem_Category) &
             master$Ecosystem_Category == eco
  out_eco <- !is.na(master$Ecosystem_Category) &
             master$Ecosystem_Category != eco
  for(cls in key_classes) {
    vals <- master[[cls]]
    vals[is.na(vals)] <- 0
    a  <- sum(vals[in_eco]  > 0)
    b  <- sum(vals[out_eco] > 0)
    cc <- sum(vals[in_eco]  == 0)
    d  <- sum(vals[out_eco] == 0)
    ft <- fisher.test(matrix(c(a,b,cc,d),2,2),
                      alternative="two.sided")
    S4_rows[[length(S4_rows)+1]] <- data.frame(
      Ecosystem        = eco,
      BGC_class        = cls,
      N_in_ecosystem   = sum(in_eco),
      N_with_class_in  = a,
      N_with_class_out = b,
      Pct_in_ecosystem = round(a/sum(in_eco)*100,  1),
      Pct_outside      = round(b/sum(out_eco)*100, 1),
      Odds_ratio       = round(ft$estimate, 3),
      CI_lower         = round(ft$conf.int[1], 3),
      CI_upper         = round(ft$conf.int[2], 3),
      P_raw            = signif(ft$p.value, 4),
      Direction        = ifelse(ft$estimate>1,"Enriched","Depleted")
    )
  }
}

S4 <- bind_rows(S4_rows) %>%
  mutate(
    P_adj_BH           = signif(p.adjust(P_raw, method="BH"), 4),
    Significant_FDR005 = P_adj_BH < 0.05
  ) %>%
  arrange(Ecosystem, P_adj_BH)

write_csv(S4,
  file.path(SUPP, "Table_S4_bgc_enrichment_complete.csv"))
cat(sprintf("  S4 saved: %d combinations (%d significant)\n\n",
            nrow(S4), sum(S4$Significant_FDR005)))

#=============================================================
# SUPPLEMENTARY TABLE S5
# BiG-SCAPE GCF summary
#=============================================================
cat("Generating Supplementary Table S5...\n")

gcf_with_eco <- gcf_bgcs %>%
  left_join(
    master %>% select(IMG.Genome.ID, Ecosystem_Category, C_Phylum),
    by = c("genome_id" = "IMG.Genome.ID")
  )

gcf_eco_summary <- gcf_with_eco %>%
  group_by(GCF_id) %>%
  summarise(
    Ecosystems = paste(unique(na.omit(Ecosystem_Category)),
                       collapse=";"),
    Phyla      = paste(unique(na.omit(C_Phylum)), collapse=";"),
    .groups    = "drop"
  )

S5 <- gcf_bgcs %>%
  group_by(GCF_id) %>%
  summarise(
    GCF_size            = n(),
    Is_singleton        = all(is_singleton),
    MIBiG_connected     = any(mibig_connected),
    N_genomes           = n_distinct(genome_id),
    BiGSCAPE_classes    = paste(unique(`BiG-SCAPE class`),
                                collapse=";"),
    Product_predictions = paste(unique(`Product Prediction`),
                                collapse=";"),
    BGC_IDs             = paste(BGC, collapse=";"),
    Genome_IDs          = paste(unique(genome_id), collapse=";"),
    .groups             = "drop"
  ) %>%
  left_join(gcf_eco_summary, by="GCF_id") %>%
  arrange(desc(GCF_size)) %>%
  select(GCF_id, GCF_size, Is_singleton, MIBiG_connected,
         N_genomes, BiGSCAPE_classes, Product_predictions,
         Ecosystems, Phyla, BGC_IDs, Genome_IDs)

write_csv(S5,
  file.path(SUPP, "Table_S5_bigscape_gcf_summary.csv"))
cat(sprintf("  S5 saved: %d GCFs (%d singletons, %d families)\n\n",
            nrow(S5), sum(S5$Is_singleton), sum(!S5$Is_singleton)))

#=============================================================
# SUPPLEMENTARY TABLE S6
# Zero-BGC genome characterization
#=============================================================
cat("Generating Supplementary Table S6...\n")

S6 <- master %>%
  filter(total_bgcs == 0) %>%
  select(
    Genome_ID           = IMG.Genome.ID,
    Genome_source       = Genome_source,
    Ecosystem           = Ecosystem_Category,
    Phylum              = C_Phylum,
    Class               = C_Class,
    Genus               = C_Genus,
    Culture_Status      = Culture_Status,
    Completeness_pct    = CheckM2.Completeness,
    Contamination_pct   = CheckM2.Contamination,
    Genome_size_Mb      = genome_size_mb,
    N_contigs           = n_contigs,
    N50_bp              = n50_bp,
    Total_BGCs          = total_bgcs,
    BGC_density_per_Mb  = bgc_per_mb,
    Known_biosurf_genes = Biosurf_Total
  ) %>%
  arrange(Ecosystem, Phylum)

write_csv(S6,
  file.path(SUPP, "Table_S6_zero_bgc_genomes.csv"), na="NA")
cat(sprintf("  S6 saved: %d zero-BGC genomes\n\n", nrow(S6)))

#=============================================================
# SUPPLEMENTARY FIGURE S1
# BGC product type discovery rarefaction curves per ecosystem
#=============================================================
cat("Generating Supplementary Figure S1 (rarefaction curves)...\n")

set.seed(42)
N_PERM <- 100

eco_short <- function(x) str_replace_all(x, c(
  "Animal-Associated \\(Non-Mammalian\\)" = "Animal-Assoc.",
  "Microbial/Fermentation"               = "Microbial/Ferm.",
  "Oil-Related Environments"             = "Oil-Related",
  "Extreme Environments"                 = "Extreme",
  "Saline Environments"                  = "Saline",
  "Built-Environment"                    = "Built-Env.",
  "Plant-Associated"                     = "Plant-Assoc.",
  "Food-Products"                        = "Food-Prod."
))

eco_colors <- c(
  "Plant-Assoc."   = "#2E8B57", "Terrestrial"    = "#8B4513",
  "Marine"         = "#1E90FF", "Freshwater"     = "#00BFFF",
  "Clinical"       = "#DC143C", "Wastewater"     = "#9932CC",
  "Air"            = "#87CEEB", "Algae"          = "#3CB371",
  "Industrial/Waste"="#FF8C00", "Built-Env."     = "#708090",
  "Extreme"        = "#FF4500", "Saline"         = "#20B2AA",
  "Geological"     = "#8B6914", "Food-Prod."     = "#DAA520",
  "Microbial/Ferm."= "#9370DB", "Animal-Assoc."  = "#FF69B4",
  "Oil-Related"    = "#4A4A4A"
)

# Build long format of BGC presence per genome
bgc_long <- master %>%
  filter(!is.na(Ecosystem_Category)) %>%
  select(IMG.Genome.ID, Ecosystem_Category, all_of(bgc_cols)) %>%
  pivot_longer(
    cols      = all_of(bgc_cols),
    names_to  = "product",
    values_to = "count"
  ) %>%
  filter(!is.na(count), count > 0)

ecosystems_list <- unique(master$Ecosystem_Category)
ecosystems_list <- ecosystems_list[!is.na(ecosystems_list)]

rarefaction_data <- list()

for(eco in ecosystems_list) {
  eco_genomes <- master %>%
    filter(Ecosystem_Category == eco) %>%
    pull(IMG.Genome.ID)
  if(length(eco_genomes) < 3) next

  max_n <- length(eco_genomes)
  steps <- unique(c(
    seq(1, min(10, max_n), 1),
    seq(10, max_n, max(1, floor(max_n/20))),
    max_n
  ))

  perm_mat <- sapply(1:N_PERM, function(p) {
    shuffled <- sample(eco_genomes)
    sapply(steps, function(n) {
      bgc_long %>%
        filter(IMG.Genome.ID %in% shuffled[1:n]) %>%
        pull(product) %>%
        n_distinct()
    })
  })

  if(is.vector(perm_mat))
    perm_mat <- matrix(perm_mat, nrow=length(steps))

  eco_label <- eco_short(eco)
  rarefaction_data[[eco_label]] <- data.frame(
    Ecosystem        = eco_label,
    N_genomes        = steps,
    Mean_BGC_classes = rowMeans(perm_mat),
    SE               = apply(perm_mat, 1, sd) / sqrt(N_PERM)
  )
  cat(sprintf("  Rarefaction done: %s (%d genomes)\n",
              eco_label, max_n))
}

raref_df <- bind_rows(rarefaction_data)

final_order <- raref_df %>%
  group_by(Ecosystem) %>%
  slice_max(N_genomes, n=1) %>%
  arrange(desc(Mean_BGC_classes)) %>%
  pull(Ecosystem)

raref_df <- raref_df %>%
  mutate(Ecosystem = factor(Ecosystem, levels=final_order))

pS1 <- ggplot(raref_df,
              aes(x=N_genomes, y=Mean_BGC_classes,
                  color=Ecosystem, fill=Ecosystem)) +
  geom_ribbon(aes(ymin=Mean_BGC_classes - SE,
                  ymax=Mean_BGC_classes + SE),
              alpha=0.12, color=NA) +
  geom_line(linewidth=0.85) +
  scale_color_manual(values=eco_colors, name="Ecosystem") +
  scale_fill_manual(values=eco_colors,  guide="none") +
  scale_x_continuous(
    name   = "Number of genomes sampled",
    expand = expansion(mult=c(0, 0.02))
  ) +
  scale_y_continuous(
    name   = "Cumulative BGC product types discovered",
    expand = expansion(mult=c(0, 0.05))
  ) +
  labs(
    title    = "Supplementary Figure S1",
    subtitle = paste("BGC product type discovery rarefaction curves",
                     "per ecosystem (mean +/- SE, 100 permutations)")
  ) +
  theme_classic(base_size=11) +
  theme(
    plot.title   = element_text(face="bold", size=12),
    plot.subtitle= element_text(size=9, color="grey40"),
    axis.title   = element_text(face="bold", size=10),
    axis.text    = element_text(size=9, color="black"),
    legend.title = element_text(face="bold", size=9),
    legend.text  = element_text(size=8),
    legend.key.size = unit(0.5,"cm"),
    panel.grid.major.y = element_line(color="grey92", linewidth=0.3),
    plot.margin  = margin(10,10,10,10)
  ) +
  guides(color=guide_legend(ncol=1, keyheight=0.8))

ggsave(file.path(SUPP,"FigureS1_rarefaction_curves.pdf"),
       pS1, width=10, height=7, dpi=300)
ggsave(file.path(SUPP,"FigureS1_rarefaction_curves.png"),
       pS1, width=10, height=7, dpi=300)
cat("  Figure S1 saved.\n\n")

#=============================================================
# SUPPLEMENTARY FIGURE S2
# Full BGC product type prevalence heatmap across ecosystems
#=============================================================
cat("Generating Supplementary Figure S2 (full heatmap)...\n")

eco_totals <- master %>%
  filter(!is.na(Ecosystem_Category)) %>%
  count(Ecosystem_Category, name="total_genomes")

product_eco_pct <- master %>%
  filter(!is.na(Ecosystem_Category)) %>%
  select(IMG.Genome.ID, Ecosystem_Category, all_of(bgc_cols)) %>%
  pivot_longer(
    cols      = all_of(bgc_cols),
    names_to  = "product",
    values_to = "count"
  ) %>%
  filter(!is.na(count)) %>%
  group_by(Ecosystem_Category, product) %>%
  summarise(n_with=sum(count>0), .groups="drop") %>%
  left_join(eco_totals, by="Ecosystem_Category") %>%
  mutate(pct=round(n_with/total_genomes*100, 1))

# Keep classes present in >= 5% of at least one ecosystem
major_products <- product_eco_pct %>%
  group_by(product) %>%
  summarise(max_pct=max(pct)) %>%
  filter(max_pct >= 5) %>%
  pull(product) %>% sort()

cat(sprintf("  BGC classes with >=5%% in at least one ecosystem: %d\n",
            length(major_products)))

heatmap_mat <- product_eco_pct %>%
  filter(product %in% major_products) %>%
  select(Ecosystem_Category, product, pct) %>%
  pivot_wider(names_from=product, values_from=pct,
              values_fill=0) %>%
  column_to_rownames("Ecosystem_Category") %>%
  as.matrix()

# Order rows by BGC density
eco_order_density <- master %>%
  group_by(Ecosystem_Category) %>%
  summarise(m=mean(bgc_per_mb, na.rm=TRUE), .groups="drop") %>%
  arrange(desc(m)) %>%
  filter(Ecosystem_Category %in% rownames(heatmap_mat)) %>%
  pull(Ecosystem_Category)

heatmap_mat <- heatmap_mat[eco_order_density, ]

col_fun_s2 <- colorRamp2(
  c(0, 10, 30, 60, 100),
  c("#F7FBFF","#BDD7E7","#6BAED6","#2171B5","#08306B")
)

ht_s2 <- Heatmap(
  t(heatmap_mat),
  name                 = "% genomes\nwith BGC class",
  col                  = col_fun_s2,
  cluster_rows         = TRUE,
  cluster_columns      = FALSE,
  show_row_dend        = TRUE,
  column_names_rot     = 45,
  row_names_gp         = gpar(fontsize=7.5),
  column_names_gp      = gpar(fontsize=8.5, fontface="bold"),
  cell_fun             = function(j,i,x,y,width,height,fill) {
    v <- t(heatmap_mat)[i,j]
    if(v >= 15) {
      grid.text(sprintf("%.0f",v), x, y,
        gp=gpar(fontsize=6,
                col=ifelse(v>=50,"white","black")))
    }
  },
  heatmap_legend_param = list(
    title_gp      = gpar(fontsize=9, fontface="bold"),
    labels_gp     = gpar(fontsize=8),
    legend_height = unit(3,"cm")
  ),
  column_title    = "Ecosystem (ordered by BGC density)",
  column_title_gp = gpar(fontsize=10, fontface="bold"),
  row_title       = "BGC product type",
  row_title_gp    = gpar(fontsize=10, fontface="bold"),
  width           = unit(10,"cm"),
  height          = unit(0.38 * length(major_products),"cm")
)

fig_h <- max(10, length(major_products) * 0.30)

pdf(file.path(SUPP,"FigureS2_full_product_type_heatmap.pdf"),
    width=14, height=fig_h)
draw(ht_s2,
  padding=unit(c(5,5,10,5),"mm"),
  column_title=paste("Complete BGC product type prevalence across",
    "17 ecosystems (classes present in >=5% of at least one ecosystem)"),
  column_title_gp=gpar(fontsize=12,fontface="bold"))
dev.off()

png(file.path(SUPP,"FigureS2_full_product_type_heatmap.png"),
    width=14, height=fig_h, units="in", res=300)
draw(ht_s2,
  padding=unit(c(5,5,10,5),"mm"),
  column_title=paste("Complete BGC product type prevalence across",
    "17 ecosystems (classes present in >=5% of at least one ecosystem)"),
  column_title_gp=gpar(fontsize=12,fontface="bold"))
dev.off()
cat("  Figure S2 saved.\n\n")

#=============================================================
# FINAL SUMMARY
#=============================================================
cat("\n=== SUPPLEMENTARY MATERIALS COMPLETE ===\n")
cat("Saved to:", SUPP, "\n\n")

cat("Tables:\n")
for(t in list.files(SUPP, pattern="Table_S.*\\.csv")) {
  sz <- file.size(file.path(SUPP,t))
  cat(sprintf("  %-50s %s bytes\n", t, format(sz, big.mark=",")))
}
cat("\nFigures:\n")
for(f in list.files(SUPP, pattern="Figure.*\\.(pdf|png)")) {
  sz <- file.size(file.path(SUPP,f))
  cat(sprintf("  %-50s %s bytes\n", f, format(sz, big.mark=",")))
}