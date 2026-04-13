suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(readr)
  library(igraph); library(ggraph); library(scales)
  library(RColorBrewer); library(stringr)
})
BASE <- "/data/habib/EcoSurfBGC"
FIGS <- file.path(BASE, "figures")

pairs <- read_csv(file.path(BASE, "results/bgc_cooccurrence_pairs.csv"),
                  show_col_types = FALSE) %>%
  filter(jaccard >= 0.15)
master <- read_csv(file.path(BASE, "results/integrated_master_table_v2.csv"),
                   show_col_types = FALSE)

key_classes <- c("NRPS", "NRPS-like", "T1PKS", "T3PKS", "PKS-like",
                 "terpene", "RiPP-like", "lanthipeptide-class-i",
                 "lanthipeptide-class-ii", "lassopeptide", "thiopeptide",
                 "ranthipeptide", "betalactone", "arylpolyene",
                 "NI-siderophore", "ectoine", "hserlactone", "NAPAA")
key_classes <- key_classes[key_classes %in% names(master)]

node_stats <- sapply(key_classes, function(cls)
  mean(master[[cls]] > 0, na.rm = TRUE) * 100)

node_df <- data.frame(
  name = key_classes,
  prevalence = node_stats[key_classes],
  category = case_when(
    key_classes %in% c("NRPS", "NRPS-like", "lanthipeptide-class-i",
      "lanthipeptide-class-ii", "lassopeptide", "thiopeptide",
      "ranthipeptide", "NAPAA") ~ "RiPP/NRP",
    key_classes %in% c("T1PKS", "T3PKS", "PKS-like") ~ "PKS",
    key_classes %in% c("terpene") ~ "Terpene",
    key_classes %in% c("RiPP-like") ~ "RiPP",
    TRUE ~ "Others"))

g <- graph_from_data_frame(
  d = pairs %>% select(class1, class2, jaccard),
  vertices = node_df, directed = FALSE)

cat_colors <- c("RiPP/NRP" = "#C0392B", "PKS" = "#8E44AD",
                "Terpene" = "#27AE60", "RiPP" = "#2ECC71", "Others" = "#2980B9")

set.seed(42)
p4 <- ggraph(g, layout = "fr") +
  geom_edge_link(aes(width = jaccard, alpha = jaccard), color = "grey50") +
  geom_node_point(aes(size = prevalence, fill = category),
                  alpha = 0.85, shape = 21, color = "white", stroke = 0.5) +
  geom_node_text(aes(label = name), size = 2.5, repel = TRUE,
                 fontface = "bold", max.overlaps = 25,
                 segment.size = 0.2, segment.color = "grey70",
                 box.padding = 0.4, seed = 42) +
  scale_edge_width_continuous(range = c(0.3, 2.5),
    name = "Jaccard index", breaks = c(0.2, 0.3, 0.4)) +
  scale_edge_alpha_continuous(range = c(0.2, 0.8), guide = "none") +
  scale_size_continuous(range = c(3, 13),
    name = "Prevalence (% genomes)", breaks = c(10, 25, 50)) +
  scale_fill_manual(values = cat_colors, name = "BGC category") +
  theme_graph(base_family = "sans", base_size = 9) +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 8),
    legend.text = element_text(size = 7),
    legend.key.size = unit(3.5, "mm"),
    plot.margin = margin(8, 8, 8, 8))

ggsave(file.path(FIGS, "Figure4_cooccurrence_network.pdf"),
       p4, width = 180, height = 160, units = "mm", dpi = 600)
ggsave(file.path(FIGS, "Figure4_cooccurrence_network.png"),
       p4, width = 180, height = 160, units = "mm", dpi = 600)
ggsave(file.path(FIGS, "Figure4_cooccurrence_network.tiff"),
       p4, width = 180, height = 160, units = "mm", dpi = 600,
       compression = "lzw")
cat("Figure 4 saved.\n")
