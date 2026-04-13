suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(readr); library(tidyr)
  library(patchwork); library(stringr); library(forcats)
  library(scales); library(ggrepel)
})

BASE <- "/data/habib/EcoSurfBGC"
FIGS <- file.path(BASE, "figures")
SUPP <- file.path(BASE, "supplementary")

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

col_imp <- "#1B7837"
col_dec <- "#C51B7D"
col_nc  <- "grey55"

cat("Loading data...\n")
master <- read_csv(
  file.path(BASE, "results/integrated_master_table_v3.csv"),
  show_col_types = FALSE) %>%
  mutate(IMG.Genome.ID = as.character(IMG.Genome.ID))

full_data   <- master %>% filter(!is.na(Ecosystem_Category))
no_mag_data <- master %>%
  filter(!is.na(Ecosystem_Category),
         Culture_Status == "Cultured")

cat(sprintf("Full dataset: %d | No-MAG: %d\n",
            nrow(full_data), nrow(no_mag_data)))

eco_sum <- function(df, label) {
  df %>%
    group_by(Ecosystem_Category) %>%
    summarise(
      N               = n(),
      mean_BGC_per_Mb = round(mean(bgc_per_mb,  na.rm = TRUE), 3),
      mean_pct_novel  = round(mean(pct_novel,    na.rm = TRUE), 1),
      pct_MAG         = round(mean(
        Culture_Status == "Uncultured/MAG", na.rm = TRUE) * 100, 1),
      .groups = "drop") %>%
    arrange(desc(mean_BGC_per_Mb)) %>%
    mutate(Rank = row_number(), Dataset = label)
}

full_eco   <- eco_sum(full_data,   "Full dataset")
no_mag_eco <- eco_sum(no_mag_data, "Cultured isolates only")

comparison <- full_eco %>%
  select(Ecosystem_Category,
         Full_N           = N,
         Full_Rank        = Rank,
         Full_BGC_density = mean_BGC_per_Mb,
         Full_pct_novel   = mean_pct_novel,
         Full_pct_MAG     = pct_MAG) %>%
  left_join(
    no_mag_eco %>%
      select(Ecosystem_Category,
             NoMAG_N           = N,
             NoMAG_Rank        = Rank,
             NoMAG_BGC_density = mean_BGC_per_Mb,
             NoMAG_pct_novel   = mean_pct_novel),
    by = "Ecosystem_Category") %>%
  mutate(
    Rank_change    = Full_Rank - NoMAG_Rank,
    Density_change = round(NoMAG_BGC_density - Full_BGC_density, 3),
    Novelty_change = round(NoMAG_pct_novel   - Full_pct_novel,   1),
    shift_dir      = case_when(
      Rank_change > 0 ~ "Improved (no-MAG)",
      Rank_change < 0 ~ "Declined (no-MAG)",
      TRUE            ~ "No change"),
    abs_change     = abs(Rank_change),
    Ecosystem_short = eco_short_fn(Ecosystem_Category)
  ) %>%
  arrange(Full_Rank)

# ══════════════════════════════════════════════════════════
# PANEL A — Slope chart (same approach as Figure 3 v8)
# Use x=1 and x=3 to spread columns apart
# Use annotate for headers, no x-axis text
# Legend outside at "right"
# ══════════════════════════════════════════════════════════
space_fn <- function(df, ycol) {
  df <- df %>% arrange(get(ycol))
  for (i in 2:nrow(df))
    if (df[[ycol]][i] - df[[ycol]][i-1] < 0.6)
      df[[ycol]][i] <- df[[ycol]][i-1] + 0.6
  df
}

left_labs  <- space_fn(
  comparison %>% arrange(Full_Rank)  %>% mutate(y_lab = Full_Rank), "y_lab")
right_labs <- space_fn(
  comparison %>% arrange(NoMAG_Rank) %>% mutate(y_lab = NoMAG_Rank), "y_lab")

pA <- ggplot(comparison) +
  geom_segment(
    aes(x = 1, xend = 3, y = Full_Rank, yend = NoMAG_Rank,
        color = shift_dir, linewidth = abs_change + 0.1),
    alpha = 0.75) +
  geom_point(aes(x = 1, y = Full_Rank,  color = shift_dir), size = 2) +
  geom_point(aes(x = 3, y = NoMAG_Rank, color = shift_dir), size = 2) +
  geom_text(data = left_labs,
    aes(x = 0.88, y = y_lab, label = Ecosystem_short, color = shift_dir),
    hjust = 1, size = 2.2, fontface = "italic", show.legend = FALSE) +
  geom_segment(data = left_labs,
    aes(x = 0.89, xend = 0.98, y = y_lab, yend = Full_Rank),
    color = "grey75", linewidth = 0.15) +
  geom_text(data = right_labs,
    aes(x = 3.12, y = y_lab, label = Ecosystem_short, color = shift_dir),
    hjust = 0, size = 2.2, fontface = "italic", show.legend = FALSE) +
  geom_segment(data = right_labs,
    aes(x = 3.02, xend = 3.11, y = NoMAG_Rank, yend = y_lab),
    color = "grey75", linewidth = 0.15) +
  # Column headers as annotations — no x-axis text overlap possible
  annotate("text", x = 1, y = -0.4,
           label = "bold('Full dataset')", parse = TRUE,
           size = 2.5, color = "grey20") +
  annotate("text", x = 1, y = -1.0,
           label = "'(with MAGs)'", parse = TRUE,
           size = 2.1, color = "grey50") +
  annotate("text", x = 3, y = -0.4,
           label = "bold('Cultured only')", parse = TRUE,
           size = 2.5, color = "grey20") +
  annotate("text", x = 3, y = -1.0,
           label = "'(no MAGs)'", parse = TRUE,
           size = 2.1, color = "grey50") +
  scale_y_reverse(name = "BGC density rank", breaks = 1:17,
                  limits = c(18, -1.5), expand = c(0, 0)) +
  scale_x_continuous(limits = c(-0.8, 4.8), breaks = NULL) +
  scale_color_manual(
    name = "Rank shift",
    values = c("Improved (no-MAG)" = col_imp,
               "Declined (no-MAG)" = col_dec,
               "No change"         = col_nc),
    guide = guide_legend(
      override.aes = list(linewidth = 1.2, size = 2))) +
  scale_linewidth_continuous(range = c(0.3, 2.2), guide = "none") +
  labs(tag = "A",
       subtitle = "Rank shifts in BGC density after excluding MAG genomes") +
  pub_theme +
  theme(
    panel.grid.major.x = element_blank(),
    axis.title.x       = element_blank(),
    axis.text.x        = element_blank(),
    axis.ticks.x       = element_blank(),
    axis.line.x        = element_blank(),
    # Legend OUTSIDE — never overlaps "Air" or any other label
    legend.position     = "right")

# ══════════════════════════════════════════════════════════
# PANEL B — Novelty scatter
# ══════════════════════════════════════════════════════════
pB <- ggplot(comparison,
  aes(x = Full_pct_novel, y = NoMAG_pct_novel,
      color = Full_pct_MAG, size = Full_N,
      label = Ecosystem_short)) +
  geom_abline(slope = 1, intercept = 0,
    linetype = "dashed", color = "grey60", linewidth = 0.5) +
  geom_point(alpha = 0.85, shape = 16) +
  geom_text_repel(size = 2.3, show.legend = FALSE,
    segment.size = 0.3, segment.color = "grey60",
    box.padding = 0.5, max.overlaps = 25, seed = 42,
    force = 1.2, force_pull = 0.5) +
  scale_color_gradient2(
    low = "#2166AC", mid = "#F7F7F7", high = "#D6604D",
    midpoint = 5, name = "% MAG in\nfull dataset") +
  scale_size_continuous(name = "N genomes", range = c(2, 8),
    breaks = c(50, 150, 300)) +
  scale_x_continuous(
    name = "% Novel BGCs \u2014 full dataset",
    expand = expansion(mult = 0.08)) +
  scale_y_continuous(
    name = "% Novel BGCs \u2014 cultured isolates only",
    expand = expansion(mult = 0.08)) +
  annotate("text", x = 44, y = 88,
    label = "Above line: novelty increases\nwhen MAGs removed",
    size = 2.1, color = "grey45", hjust = 0, fontface = "italic") +
  annotate("text", x = 44, y = 43,
    label = "Below line: novelty decreases\nwhen MAGs removed",
    size = 2.1, color = "grey45", hjust = 0, fontface = "italic") +
  labs(tag = "B",
       subtitle = "Points above dashed line: MAG removal increases apparent novelty") +
  pub_theme +
  theme(legend.position = "right",
        legend.key.width = unit(4, "mm"))

# ══════════════════════════════════════════════════════════
# PANEL C — BGC density paired dots
# ══════════════════════════════════════════════════════════
density_long <- comparison %>%
  select(Ecosystem_short,
         `Full dataset`  = Full_BGC_density,
         `Cultured only` = NoMAG_BGC_density) %>%
  pivot_longer(-Ecosystem_short,
    names_to = "Dataset", values_to = "BGC_density") %>%
  mutate(Ecosystem_short = fct_reorder(
    Ecosystem_short,
    BGC_density * (Dataset == "Full dataset"),
    .fun = max))

pC <- ggplot(density_long,
  aes(x = Ecosystem_short, y = BGC_density,
      fill = Dataset, color = Dataset)) +
  geom_line(aes(group = Ecosystem_short),
    color = "grey70", linewidth = 0.4) +
  geom_point(size = 2.5, shape = 21, stroke = 0.5) +
  scale_fill_manual(
    values = c("Full dataset"  = "#3498DB",
               "Cultured only" = "#E74C3C"),
    name = "Dataset") +
  scale_color_manual(
    values = c("Full dataset"  = "#3498DB",
               "Cultured only" = "#E74C3C"),
    name = "Dataset", guide = "none") +
  scale_y_continuous(
    name   = "Mean BGC density (BGCs/Mb)",
    expand = expansion(mult = c(0.05, 0.1))) +
  xlab(NULL) +
  labs(tag = "C",
       subtitle = "Blue = full dataset; Red = cultured isolates only; connected by grey line") +
  pub_theme +
  theme(
    axis.text.x = element_text(angle = 40, hjust = 1, size = 7.5),
    # Legend OUTSIDE on the right
    legend.position = "right")

# ══════════════════════════════════════════════════════════
# Layout
# ══════════════════════════════════════════════════════════
fig_out <- (pA + pB + plot_layout(widths = c(1, 1.1))) /
            pC +
            plot_layout(heights = c(1.4, 0.9)) +
            plot_annotation(
              theme = theme(plot.margin = margin(5, 5, 5, 5)))

ggsave(file.path(SUPP, "FigureS3_sensitivity_analysis.pdf"),
       fig_out, width = 220, height = 210,
       units = "mm", dpi = 600)
ggsave(file.path(SUPP, "FigureS3_sensitivity_analysis.png"),
       fig_out, width = 220, height = 210,
       units = "mm", dpi = 600)
ggsave(file.path(SUPP, "FigureS3_sensitivity_analysis.tiff"),
       fig_out, width = 220, height = 210,
       units = "mm", dpi = 600, compression = "lzw")

cat("Figure S3 saved to supplementary/\n")
