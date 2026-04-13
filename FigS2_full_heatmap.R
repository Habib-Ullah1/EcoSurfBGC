suppressPackageStartupMessages({
  library(ComplexHeatmap); library(circlize)
  library(dplyr); library(readr); library(tidyr)
  library(tibble); library(grid)
})
BASE <- "/data/habib/EcoSurfBGC"
SUPP <- file.path(BASE, "supplementary")

master <- read_csv(file.path(BASE, "results/integrated_master_table_v2.csv"),
                   show_col_types = FALSE)

all_cols <- names(master)
start_idx <- which(all_cols == "total_bgcs") + 1
end_idx   <- which(all_cols == "Culture_Status") - 1
bgc_cols  <- all_cols[start_idx:end_idx]
bgc_cols  <- bgc_cols[sapply(bgc_cols, function(x) is.numeric(master[[x]]))]

eco_totals <- master %>%
  filter(!is.na(Ecosystem_Category)) %>%
  count(Ecosystem_Category, name = "total_genomes")

product_eco_pct <- master %>%
  filter(!is.na(Ecosystem_Category)) %>%
  select(IMG.Genome.ID, Ecosystem_Category, all_of(bgc_cols)) %>%
  pivot_longer(cols = all_of(bgc_cols), names_to = "product", values_to = "count") %>%
  filter(!is.na(count)) %>%
  group_by(Ecosystem_Category, product) %>%
  summarise(n_with = sum(count > 0), .groups = "drop") %>%
  left_join(eco_totals, by = "Ecosystem_Category") %>%
  mutate(pct = round(n_with / total_genomes * 100, 1))

major_products <- product_eco_pct %>%
  group_by(product) %>%
  summarise(max_pct = max(pct)) %>%
  filter(max_pct >= 5) %>%
  pull(product) %>% sort()

cat(sprintf("BGC classes with >=5%% in at least one ecosystem: %d\n",
            length(major_products)))

heatmap_mat <- product_eco_pct %>%
  filter(product %in% major_products) %>%
  select(Ecosystem_Category, product, pct) %>%
  pivot_wider(names_from = product, values_from = pct, values_fill = 0) %>%
  column_to_rownames("Ecosystem_Category") %>%
  as.matrix()

eco_order_density <- master %>%
  group_by(Ecosystem_Category) %>%
  summarise(m = mean(bgc_per_mb, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(m)) %>%
  filter(Ecosystem_Category %in% rownames(heatmap_mat)) %>%
  pull(Ecosystem_Category)
heatmap_mat <- heatmap_mat[eco_order_density, ]

col_fun_s2 <- colorRamp2(c(0, 10, 30, 60, 100),
  c("#F7FBFF", "#BDD7E7", "#6BAED6", "#2171B5", "#08306B"))

ht_s2 <- Heatmap(
  t(heatmap_mat),
  name = "% genomes\nwith BGC class",
  col = col_fun_s2,
  cluster_rows = TRUE, cluster_columns = FALSE,
  show_row_dend = TRUE, column_names_rot = 45,
  row_names_gp = gpar(fontsize = 7),
  column_names_gp = gpar(fontsize = 8, fontface = "bold"),
  cell_fun = function(j, i, x, y, width, height, fill) {
    v <- t(heatmap_mat)[i, j]
    if (v >= 15) grid.text(sprintf("%.0f", v), x, y,
      gp = gpar(fontsize = 5.5, col = ifelse(v >= 50, "white", "black")))
  },
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 8, fontface = "bold"),
    labels_gp = gpar(fontsize = 7),
    legend_height = unit(3, "cm")),
  column_title = "Ecosystem (ordered by BGC density)",
  column_title_gp = gpar(fontsize = 9, fontface = "bold"),
  row_title = "BGC product type",
  row_title_gp = gpar(fontsize = 9, fontface = "bold"),
  width = unit(10, "cm"),
  height = unit(0.38 * length(major_products), "cm"))

fig_h <- max(10, length(major_products) * 0.30)

pdf(file.path(SUPP, "FigureS2_full_product_type_heatmap.pdf"),
    width = 14, height = fig_h)
draw(ht_s2, padding = unit(c(5, 5, 10, 5), "mm"),
     heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

png(file.path(SUPP, "FigureS2_full_product_type_heatmap.png"),
    width = 14, height = fig_h, units = "in", res = 600)
draw(ht_s2, padding = unit(c(5, 5, 10, 5), "mm"),
     heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

tiff(file.path(SUPP, "FigureS2_full_product_type_heatmap.tiff"),
     width = 14, height = fig_h, units = "in", res = 600,
     compression = "lzw")
draw(ht_s2, padding = unit(c(5, 5, 10, 5), "mm"),
     heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()
cat("Figure S2 saved.\n")
