suppressPackageStartupMessages({
  library(ComplexHeatmap); library(circlize)
  library(dplyr); library(readr); library(tibble)
  library(RColorBrewer); library(grid)
})
BASE <- "/data/habib/EcoSurfBGC"
FIGS <- file.path(BASE, "figures")

master <- read_csv(file.path(BASE, "results/integrated_master_table_v2.csv"),
  show_col_types = FALSE) %>%
  filter(!is.na(C_Phylum), C_Phylum != "", C_Phylum != "NA")

key_classes <- c("NRPS", "NRPS-like", "T1PKS", "T3PKS", "terpene", "RiPP-like",
                 "lanthipeptide-class-i", "lanthipeptide-class-ii",
                 "lassopeptide", "thiopeptide", "ranthipeptide", "betalactone",
                 "arylpolyene", "NI-siderophore", "ectoine", "hserlactone", "NAPAA")
key_classes <- key_classes[key_classes %in% names(master)]

phylum_mat <- master %>%
  group_by(C_Phylum) %>% filter(n() >= 10) %>%
  summarise(n_genomes = n(),
            bgc_density = round(mean(bgc_per_mb, na.rm = TRUE), 2),
            across(all_of(key_classes),
                   ~ round(mean(. > 0, na.rm = TRUE) * 100, 1)),
            .groups = "drop") %>%
  arrange(desc(bgc_density))

mat <- phylum_mat %>% select(all_of(key_classes)) %>% as.matrix()
rownames(mat) <- phylum_mat$C_Phylum

col_labels <- c("NRPS" = "NRPS", "NRPS-like" = "NRPS-like",
  "T1PKS" = "T1PKS", "T3PKS" = "T3PKS", "terpene" = "Terpene",
  "RiPP-like" = "RiPP-like", "lanthipeptide-class-i" = "Lanthi-I",
  "lanthipeptide-class-ii" = "Lanthi-II", "lassopeptide" = "Lasso",
  "thiopeptide" = "Thiopeptide", "ranthipeptide" = "Ranthipeptide",
  "betalactone" = "Betalactone", "arylpolyene" = "Arylpolyene",
  "NI-siderophore" = "Siderophore", "ectoine" = "Ectoine",
  "hserlactone" = "Hserlactone", "NAPAA" = "NAPAA")
colnames(mat) <- col_labels[colnames(mat)]

col_fun <- colorRamp2(c(0, 20, 50, 75, 100),
  c("#F7FBFF", "#BDD7E7", "#6BAED6", "#2171B5", "#08306B"))
density_col <- colorRamp2(c(0, max(phylum_mat$bgc_density)),
  c("#FFF7EC", "#7F2704"))
n_col <- colorRamp2(c(10, max(phylum_mat$n_genomes)),
  c("#F7FCF0", "#00441B"))

ha_right <- rowAnnotation(
  `BGC density\n(BGCs/Mb)` = phylum_mat$bgc_density,
  `N genomes` = phylum_mat$n_genomes,
  col = list(`BGC density\n(BGCs/Mb)` = density_col,
             `N genomes` = n_col),
  annotation_name_gp = gpar(fontsize = 7.5, fontface = "bold"),
  show_legend = TRUE,
  annotation_legend_param = list(
    `BGC density\n(BGCs/Mb)` = list(title_gp = gpar(fontsize = 8, fontface = "bold"),
                                    labels_gp = gpar(fontsize = 7)),
    `N genomes` = list(title_gp = gpar(fontsize = 8, fontface = "bold"),
                       labels_gp = gpar(fontsize = 7))
  ))

col_order <- colnames(mat)
bgc_categories <- c(
  "NRPS" = "NRPS/PKS", "NRPS-like" = "NRPS/PKS",
  "T1PKS" = "NRPS/PKS", "T3PKS" = "NRPS/PKS", "Terpene" = "Terpene",
  "RiPP-like" = "RiPP", "Lanthi-I" = "RiPP", "Lanthi-II" = "RiPP",
  "Lasso" = "RiPP", "Thiopeptide" = "RiPP", "Ranthipeptide" = "RiPP",
  "Betalactone" = "Other", "Arylpolyene" = "Other", "Siderophore" = "Other",
  "Ectoine" = "Other", "Hserlactone" = "Other", "NAPAA" = "Other")
cat_colors <- c("NRPS/PKS" = "#C0392B", "Terpene" = "#27AE60",
                "RiPP" = "#2980B9", "Other" = "#F39C12")

ha_top <- HeatmapAnnotation(
  Category = bgc_categories[col_order],
  col = list(Category = cat_colors),
  annotation_name_gp = gpar(fontsize = 7.5, fontface = "bold"),
  show_legend = TRUE, annotation_height = unit(4, "mm"),
  annotation_legend_param = list(
    Category = list(title_gp = gpar(fontsize = 8, fontface = "bold"),
                    labels_gp = gpar(fontsize = 7))
  ))

ht <- Heatmap(mat,
  name = "% genomes\nwith BGC class",
  col = col_fun,
  top_annotation = ha_top,
  right_annotation = ha_right,
  cluster_rows = TRUE, cluster_columns = TRUE,
  show_row_dend = TRUE, show_column_dend = TRUE,
  row_dend_side = "left", column_names_rot = 45,
  row_names_gp = gpar(fontsize = 8.5, fontface = "italic"),
  column_names_gp = gpar(fontsize = 8),
  cell_fun = function(j, i, x, y, width, height, fill) {
    v <- mat[i, j]
    if (v >= 25) grid.text(sprintf("%.0f", v), x, y,
      gp = gpar(fontsize = 6.5, col = ifelse(v >= 55, "white", "black"),
                fontface = "bold"))
  },
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 8, fontface = "bold"),
    labels_gp = gpar(fontsize = 7),
    legend_height = unit(3, "cm")),
  column_title = "BGC class",
  column_title_gp = gpar(fontsize = 9, fontface = "bold"),
  row_title = "Microbial phylum (n >= 10 genomes)",
  row_title_gp = gpar(fontsize = 9, fontface = "bold"),
  width = unit(10, "cm"),
  height = unit(0.6 * nrow(mat), "cm"))

pdf(file.path(FIGS, "Figure6_phylum_heatmap.pdf"), width = 12, height = 9)
draw(ht, padding = unit(c(5, 5, 10, 5), "mm"),
     heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()
png(file.path(FIGS, "Figure6_phylum_heatmap.png"),
    width = 12, height = 9, units = "in", res = 600)
draw(ht, padding = unit(c(5, 5, 10, 5), "mm"),
     heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()
tiff(file.path(FIGS, "Figure6_phylum_heatmap.tiff"),
     width = 12, height = 9, units = "in", res = 600, compression = "lzw")
draw(ht, padding = unit(c(5, 5, 10, 5), "mm"),
     heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()
cat("Figure 6 saved.\n")
