suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(readr); library(tidyr)
  library(patchwork); library(scales); library(forcats)
  library(stringr); library(ggrepel); library(RColorBrewer)
})
BASE <- "/data/habib/EcoSurfBGC"
FIGS <- file.path(BASE, "figures")
dir.create(FIGS, showWarnings = FALSE)

pub_theme <- theme_classic(base_size = 9) +
  theme(
    text              = element_text(family = "sans", color = "black"),
    axis.title        = element_text(face = "bold", size = 9),
    axis.text         = element_text(size = 8, color = "black"),
    axis.line         = element_line(linewidth = 0.4, color = "black"),
    axis.ticks        = element_line(linewidth = 0.3, color = "black"),
    axis.ticks.length = unit(1.5, "mm"),
    legend.title      = element_text(face = "bold", size = 8),
    legend.text       = element_text(size = 7),
    legend.key.size   = unit(3.5, "mm"),
    legend.background = element_blank(),
    panel.grid.major  = element_line(color = "grey92", linewidth = 0.25),
    panel.grid.minor  = element_blank(),
    plot.margin       = margin(6, 8, 6, 6),
    plot.tag          = element_text(face = "bold", size = 12),
    plot.tag.position = c(0, 1)
  )

eco_colors <- c(
  "Plant-Associated" = "#2E8B57", "Terrestrial" = "#8B4513",
  "Marine" = "#1E90FF", "Freshwater" = "#00BFFF", "Clinical" = "#DC143C",
  "Wastewater" = "#9932CC", "Air" = "#87CEEB", "Algae" = "#3CB371",
  "Industrial/Waste" = "#FF8C00", "Built-Environment" = "#708090",
  "Extreme Environments" = "#FF4500", "Saline Environments" = "#20B2AA",
  "Geological" = "#8B6914", "Food-Products" = "#DAA520",
  "Microbial/Fermentation" = "#9370DB",
  "Animal-Associated (Non-Mammalian)" = "#FF69B4",
  "Oil-Related Environments" = "#4A4A4A")

eco_short <- function(x) str_replace_all(x,
  c("Animal-Associated [(]Non-Mammalian[)]" = "Animal-Assoc.",
    "Microbial/Fermentation" = "Microbial/Ferm.",
    "Oil-Related Environments" = "Oil-Related",
    "Extreme Environments" = "Extreme",
    "Saline Environments" = "Saline",
    "Built-Environment" = "Built-Env.",
    "Plant-Associated" = "Plant-Assoc.",
    "Food-Products" = "Food-Prod."))

eco <- read_csv(file.path(BASE, "results/ecosystem_complete_summary.csv"),
                show_col_types = FALSE) %>%
  mutate(Ecosystem_short = eco_short(Ecosystem_Category),
         Ecosystem_Category = fct_reorder(Ecosystem_Category, mean_BGC_per_Mb))

p1r <- read_csv(file.path(BASE, "results/paper1_vs_paper2_comparison.csv"),
                show_col_types = FALSE) %>%
  select(Ecosystem_Category, P1_rank, P2_rank)
eco <- eco %>% left_join(p1r, by = "Ecosystem_Category")

# Panel A — legend moved to RIGHT SIDE (outside plot)
pA <- ggplot(eco, aes(x = mean_BGC_per_Mb, y = pct_novel_BGCs,
                      size = n, color = Ecosystem_Category,
                      label = Ecosystem_short)) +
  geom_point(alpha = 0.85, stroke = 0.3) +
  geom_text_repel(size = 2.5, fontface = "italic", max.overlaps = 25,
                  show.legend = FALSE, segment.size = 0.2,
                  segment.color = "grey60", box.padding = 0.5,
                  seed = 42, force = 1.2, force_pull = 0.5,
                  min.segment.length = 0.2) +
  scale_size_continuous(name = "N genomes", range = c(2, 12),
                        breaks = c(30, 100, 200, 400)) +
  scale_color_manual(values = eco_colors, guide = "none") +
  scale_x_continuous(name = "Mean BGC density (BGCs per Mb)",
                     expand = expansion(mult = 0.12)) +
  scale_y_continuous(name = "Novel BGCs (% with no MIBiG match)",
                     labels = function(x) paste0(x, "%"),
                     expand = expansion(mult = 0.08)) +
  labs(tag = "A") +
  pub_theme +
  theme(
    # OUTSIDE the plot on the right — never overlaps data
    legend.position = "right")

master <- read_csv(file.path(BASE, "results/integrated_master_table_v2.csv"),
                   show_col_types = FALSE)
eco_order <- levels(fct_reorder(eco$Ecosystem_Category,
                                eco$mean_BGC_per_Mb, .desc = TRUE))

key_classes <- c("NRPS", "NRPS-like", "T1PKS", "T3PKS", "terpene",
                 "RiPP-like", "betalactone", "arylpolyene",
                 "NI-siderophore", "ectoine", "lassopeptide", "hserlactone")

class_category <- c(
  "NRPS" = "NRPS/PKS", "NRPS-like" = "NRPS/PKS",
  "T1PKS" = "NRPS/PKS", "T3PKS" = "NRPS/PKS",
  "terpene" = "Terpene", "RiPP-like" = "RiPP",
  "lassopeptide" = "RiPP", "betalactone" = "Others",
  "arylpolyene" = "Others", "NI-siderophore" = "Others",
  "ectoine" = "Otheris", "hserlactone" = "Others")
cat_colors <- c("NRPS/PKS" = "#C0392B", "Terpene" = "#27AE60",
                "RiPP" = "#2980B9", "Others" = "#F39C12")

dot_data <- master %>%
  group_by(Ecosystem_Category) %>%
  summarise(across(all_of(key_classes),
    ~ round(mean(. > 0, na.rm = TRUE) * 100, 1)), .groups = "drop") %>%
  pivot_longer(-Ecosystem_Category,
               names_to = "BGC_class", values_to = "pct") %>%
  mutate(
    Ecosystem_Category = factor(Ecosystem_Category, levels = eco_order),
    BGC_class = factor(BGC_class, levels = rev(key_classes)),
    category = class_category[as.character(BGC_class)])

pB <- ggplot(dot_data, aes(x = Ecosystem_Category, y = BGC_class,
                           size = pct, color = category, alpha = pct)) +
  geom_point() +
  scale_size_continuous(name = "% genomes with\nBGC class",
    range = c(0.3, 9), breaks = c(10, 25, 50, 75),
    labels = c("10%", "25%", "50%", "75%")) +
  scale_alpha_continuous(range = c(0.3, 1), guide = "none") +
  scale_color_manual(values = cat_colors, name = "BGC category") +
  scale_x_discrete(labels = eco_short) +
  xlab(NULL) + ylab("BGC class") +
  labs(tag = "B") +
  pub_theme +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7.5),
    axis.text.y = element_text(size = 8),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    legend.position = "right")

fig_out <- pA / pB + plot_layout(heights = c(1, 1.1))
ggsave(file.path(FIGS, "Figure1_BGC_landscape.pdf"),
       fig_out, width = 200, height = 210, units = "mm", dpi = 600)
ggsave(file.path(FIGS, "Figure1_BGC_landscape.png"),
       fig_out, width = 200, height = 210, units = "mm", dpi = 600)
ggsave(file.path(FIGS, "Figure1_BGC_landscape.tiff"),
       fig_out, width = 200, height = 210, units = "mm", dpi = 600,
       compression = "lzw")
cat("Figure 1 saved.\n")
