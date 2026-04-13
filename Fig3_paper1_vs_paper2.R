suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(readr); library(tidyr)
  library(patchwork); library(scales); library(ggrepel)
  library(RColorBrewer); library(stringr); library(forcats)
  library(grid); library(gridExtra); library(cowplot)
})

BASE <- "/data/habib/EcoSurfBGC"
FIGS <- file.path(BASE, "figures")

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
    legend.background = element_rect(fill = "white", color = NA),
    legend.box.background = element_blank(),
    panel.grid.major  = element_line(color = "grey92", linewidth = 0.25),
    panel.grid.minor  = element_blank(),
    plot.margin       = margin(6, 8, 6, 6),
    plot.tag          = element_text(face = "bold", size = 12),
    plot.tag.position = c(0, 1),
    plot.subtitle     = element_text(size = 7.5, color = "grey40",
                                     margin = margin(b = 4))
  )

eco_short_fn <- function(x) str_replace_all(x,
  c("Animal-Associated [(]Non-Mammalian[)]" = "Animal-Assoc.",
    "Microbial/Fermentation"               = "Microbial/Ferm.",
    "Oil-Related Environments"              = "Oil-Related",
    "Extreme Environments"                  = "Extreme",
    "Saline Environments"                   = "Saline",
    "Built-Environment"                     = "Built-Env.",
    "Plant-Associated"                      = "Plant-Assoc.",
    "Food-Products"                         = "Food-Prod.",
    "Industrial/Waste"                      = "Ind./Waste"))

col_up   <- "#1B7837"
col_down <- "#C51B7D"
col_nc   <- "grey55"

eco_colors <- c(
  "Plant-Assoc."    = "#2E8B57", "Terrestrial"     = "#8B4513",
  "Marine"          = "#1E90FF", "Freshwater"       = "#00BFFF",
  "Clinical"        = "#DC143C", "Wastewater"       = "#7B3F9E",
  "Air"             = "#87CEEB", "Algae"            = "#3CB371",
  "Ind./Waste"      = "#E67E22", "Built-Env."       = "#5D6D7E",
  "Extreme"         = "#E74C3C", "Saline"           = "#16A085",
  "Geological"      = "#8B6914", "Food-Prod."       = "#D4AC0D",
  "Microbial/Ferm." = "#8E44AD", "Animal-Assoc."    = "#E91E8C",
  "Oil-Related"     = "#2C3E50")

pred_colors <- c(
  "Genome"   = "#3498DB",
  "Taxonomy" = "#E74C3C",
  "Combined" = "#8E44AD",
  "Ecology"  = "#27AE60",
  "Ullah et al. 2026" = "#E67E22")

# ══════════════════════════════════════════════════════════
# PANEL A
# ══════════════════════════════════════════════════════════
rank_data <- read_csv(file.path(BASE, "results/paper1_vs_paper2_comparison.csv"),
                      show_col_types = FALSE) %>%
  mutate(
    rank_change = P1_rank - P2_rank,
    change_dir  = case_when(
      rank_change > 0 ~ "Increased",
      rank_change < 0 ~ "Decreased",
      TRUE            ~ "No change"),
    abs_change  = abs(rank_change),
    eco_short   = eco_short_fn(Ecosystem_Category))

slope_data <- rank_data %>%
  select(eco_short, P1_rank, P2_rank, change_dir, abs_change)

space_labels <- function(df, ycol = "y_label") {
  df <- df %>% arrange(get(ycol))
  min_gap <- 0.58
  for (i in 2:nrow(df)) {
    if (df[[ycol]][i] - df[[ycol]][i-1] < min_gap) {
      df[[ycol]][i] <- df[[ycol]][i-1] + min_gap
    }
  }
  df
}

left_labels <- slope_data %>% arrange(P1_rank) %>% mutate(y_label = P1_rank)
left_labels <- space_labels(left_labels)
right_labels <- slope_data %>% arrange(P2_rank) %>% mutate(y_label = P2_rank)
right_labels <- space_labels(right_labels)

pA <- ggplot(slope_data) +
  geom_segment(aes(x = 1, xend = 3, y = P1_rank, yend = P2_rank,
                   color = change_dir, linewidth = abs_change),
               alpha = 0.70) +
  geom_point(aes(x = 1, y = P1_rank, color = change_dir), size = 2) +
  geom_point(aes(x = 3, y = P2_rank, color = change_dir), size = 2) +
  geom_text(data = left_labels,
    aes(x = 0.93, y = y_label, label = eco_short, color = change_dir),
    hjust = 1, size = 2.2, fontface = "italic", show.legend = FALSE) +
  geom_segment(data = left_labels,
    aes(x = 0.94, xend = 0.99, y = y_label, yend = P1_rank),
    color = "grey75", linewidth = 0.15) +
  geom_text(data = right_labels,
    aes(x = 3.07, y = y_label, label = eco_short, color = change_dir),
    hjust = 0, size = 2.2, fontface = "italic", show.legend = FALSE) +
  geom_segment(data = right_labels,
    aes(x = 3.01, xend = 3.06, y = P2_rank, yend = y_label),
    color = "grey75", linewidth = 0.15) +
  # Split headers onto multiple lines so they stay narrow
  annotate("text", x = 1, y = -0.4,
           label = "bold('Ullah et al.')", parse = TRUE,
           size = 2.7, color = "grey20") +
  annotate("text", x = 1, y = -1.0,
           label = "bold('2026')", parse = TRUE,
           size = 2.7, color = "grey20") +
  annotate("text", x = 1, y = -1.6,
           label = "'(biosurfactant genes)'", parse = TRUE,
           size = 2.1, color = "grey50") +
  annotate("text", x = 3, y = -0.4,
           label = "bold('study')", parse = TRUE,
           size = 2.7, color = "grey20") +
  annotate("text", x = 3, y = -1.0,
           label = "bold('Current')", parse = TRUE,
           size = 2.7, color = "grey20") +
  annotate("text", x = 3, y = -1.6,
           label = "'(BGC density)'", parse = TRUE,
           size = 2.1, color = "grey50") +
  scale_y_reverse(name = "Ecosystem rank", breaks = 1:17,
                  limits = c(18, -2.2), expand = c(0, 0)) +
  scale_x_continuous(limits = c(-0.8, 4.5), breaks = NULL) +
  scale_color_manual(
    values = c("Increased" = col_up, "Decreased" = col_down,
               "No change" = col_nc), guide = "none") +
  scale_linewidth_continuous(range = c(0.3, 2.2), guide = "none") +
  labs(tag = "A") +
  pub_theme +
  theme(
    panel.grid.major.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x  = element_blank(),
    legend.position = "none")

# ══════════════════════════════════════════════════════════
# PANEL B (no legend)
# ══════════════════════════════════════════════════════════
master <- read_csv(file.path(BASE, "results/integrated_master_table_v2.csv"),
                   show_col_types = FALSE)

set.seed(42)
plot_data <- master %>%
  filter(!is.na(bgc_per_mb), !is.na(Biosurf_Total)) %>%
  sample_n(min(nrow(.), 2000)) %>%
  mutate(Eco_short = eco_short_fn(Ecosystem_Category))

rho_val <- cor(master$Biosurf_Total, master$total_bgcs,
               method = "spearman", use = "complete.obs")

pB <- ggplot(plot_data, aes(x = Biosurf_Total, y = bgc_per_mb)) +
  geom_point(aes(color = Eco_short),
             alpha = 0.40, size = 0.9, stroke = 0, shape = 16) +
  geom_smooth(method = "lm", color = "black", linewidth = 0.8,
              se = TRUE, fill = "grey80", alpha = 0.3) +
  scale_color_manual(values = eco_colors, name = "Ecosystem",
                     guide = guide_legend(
                       nrow = 3,
                       override.aes = list(alpha = 1, size = 2),
                       title.position = "left")) +
  scale_x_continuous(name = "Known biosurfactant gene count",
                     breaks = seq(0, 18, 3)) +
  scale_y_continuous(name = "BGC density (BGCs / Mb)") +
  annotate("label", x = 0.3,
           y = max(plot_data$bgc_per_mb, na.rm = TRUE) * 0.96,
           label = paste0("rho == ", round(rho_val, 3),
                          " ~~ italic(p) < 2.2 %*% 10^{-16} ~~ italic(n) == '2,057'"),
           parse = TRUE, hjust = 0, size = 2.6,
           fill = "white", label.size = 0.3, label.r = unit(1, "mm"),
           color = "grey20") +
  labs(tag = "B") +
  pub_theme +
  theme(
    legend.position   = "bottom",
    legend.margin      = margin(t = 2, b = 0),
    legend.key.height  = unit(2.5, "mm"),
    legend.key.width   = unit(3, "mm"),
    legend.text        = element_text(size = 6.5),
    legend.title       = element_text(size = 7.5, face = "bold"),
    legend.spacing.x   = unit(1, "mm"),
    legend.box.spacing = unit(1, "mm"),
    plot.tag = element_text(face = "bold", size = 12))

# ══════════════════════════════════════════════════════════
# PANEL C
# ══════════════════════════════════════════════════════════
var_data <- data.frame(
  Predictor = c("Full model", "Genome size", "Ecosystem + Phylum",
                "Phylum", "Gene count (Ullah et al. 2026)", "Ecosystem only"),
  R2        = c(0.415, 0.288, 0.265, 0.224, 0.221, 0.093),
  Category  = c("Combined", "Genome", "Combined",
                "Taxonomy", "Ullah et al. 2026", "Ecology"),
  stringsAsFactors = FALSE)

var_data$Predictor <- factor(var_data$Predictor,
                             levels = var_data$Predictor[order(var_data$R2)])

pC <- ggplot(var_data, aes(x = Predictor, y = R2, fill = Category)) +
  geom_col(width = 0.6, color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.3f", R2)),
            hjust = -0.1, size = 2.8, fontface = "bold", color = "grey20") +
  scale_fill_manual(values = pred_colors, name = "Predictor type") +
  scale_y_continuous(
    name = expression(bold("Variance explained") ~~ (italic(R)^2)),
    expand = expansion(mult = c(0, 0.18)),
    limits = c(0, 0.55)) +
  xlab(NULL) +
  coord_flip() +
  labs(tag = "C",
       subtitle = "Genome size is the strongest single predictor; phylogeny explains more variance than ecosystem") +
  pub_theme +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey92"),
    legend.position      = c(0.90, 0.30),
    legend.justification = c(1, 0),
    legend.background    = element_rect(fill = alpha("white", 0.9),
                                        color = "grey80", linewidth = 0.3),
    plot.tag = element_text(face = "bold", size = 12))

# ══════════════════════════════════════════════════════════
# EXTRACT LEGENDS INTO DEDICATED FULL-WIDTH ROW
# ══════════════════════════════════════════════════════════
p_rank_leg <- ggplot(data.frame(x = c(1,2), y = c(1,1),
                                d = c("Increased","Decreased")),
                     aes(x, y, color = d)) +
  geom_point(size = 2) + geom_line(linewidth = 1.2) +
  scale_color_manual(name = "Rank change",
    values = c("Increased" = col_up, "Decreased" = col_down),
    guide = guide_legend(override.aes = list(linewidth = 1.2, size = 2))) +
  theme_void() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.title = element_text(face = "bold", size = 8),
        legend.text = element_text(size = 7),
        legend.key.width = unit(6, "mm"))

leg_rank <- cowplot::get_legend(p_rank_leg)
leg_eco  <- cowplot::get_legend(pB)
pB_clean <- pB + theme(legend.position = "none")

legend_row <- wrap_elements(
  full = arrangeGrob(
    grid::textGrob(""),  # tiny left spacer
    leg_rank,
    leg_eco,
    nrow = 1,
    widths = c(0.03, 0.25, 0.72)
  )
)

fig_out <- (pA | pB_clean) / legend_row / pC +
  plot_layout(heights = c(1.4, 0.18, 0.7)) +
  plot_annotation(theme = theme(plot.margin = margin(5, 5, 5, 5)))

ggsave(file.path(FIGS, "Figure3_paper1_vs_paper2_v8.pdf"),
       fig_out, width = 220, height = 220, units = "mm", dpi = 600)
ggsave(file.path(FIGS, "Figure3_paper1_vs_paper2_v8.png"),
       fig_out, width = 220, height = 220, units = "mm", dpi = 600)
ggsave(file.path(FIGS, "Figure3_paper1_vs_paper2_v8.tiff"),
       fig_out, width = 220, height = 220, units = "mm", dpi = 600,
       compression = "lzw")

cat("Figure 3 v8 saved.\n")
