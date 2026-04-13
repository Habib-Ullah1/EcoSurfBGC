suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(readr); library(tidyr)
  library(patchwork); library(scales); library(forcats); library(stringr)
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
    legend.background = element_blank(),
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.25),
    panel.grid.minor  = element_blank(),
    plot.margin       = margin(6, 8, 6, 6),
    plot.tag          = element_text(face = "bold", size = 12),
    plot.tag.position = c(0, 1),
    plot.subtitle     = element_text(size = 7.5, color = "grey40",
                                     margin = margin(b = 4))
  )

eco_short <- function(x) str_replace_all(x,
  c("Animal-Associated [(]Non-Mammalian[)]" = "Animal-Assoc.",
    "Microbial/Fermentation" = "Microbial/Ferm.",
    "Oil-Related Environments" = "Oil-Related",
    "Extreme Environments" = "Extreme",
    "Saline Environments" = "Saline",
    "Built-Environment" = "Built-Env.",
    "Plant-Associated" = "Plant-Assoc.",
    "Food-Products" = "Food-Prod."))

master_full <- read_csv(file.path(BASE, "results/integrated_master_table_v2.csv"),
                        show_col_types = FALSE)
master <- master_full %>%
  filter(!is.na(Culture_Status), Culture_Status != "Unknown",
         !is.na(bgc_per_mb), !is.na(Ecosystem_Category))

eco_order <- master %>%
  filter(Culture_Status == "Cultured") %>%
  group_by(Ecosystem_Category) %>%
  summarise(m = median(bgc_per_mb, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(m)) %>% pull(Ecosystem_Category)

master <- master %>%
  mutate(Ecosystem_f = factor(eco_short(Ecosystem_Category),
                              levels = eco_short(eco_order)))

# Panel A — Genome size vs zero BGCs
all_bgcs <- master_full %>%
  mutate(IMG.Genome.ID = as.character(IMG.Genome.ID)) %>%
  select(IMG.Genome.ID, total_bgcs_master = total_bgcs)

density_df <- read_csv(file.path(BASE, "results/bgc_density.csv"),
  show_col_types = FALSE) %>%
  mutate(genome_id = as.character(genome_id)) %>%
  left_join(all_bgcs, by = c("genome_id" = "IMG.Genome.ID")) %>%
  filter(!is.na(total_bgcs_master)) %>%
  mutate(
    size_bin = cut(genome_size_mb,
      breaks = c(0, 0.5, 1, 1.5, 2, 3, 4, 5, 10, Inf),
      labels = c("<0.5", "0.5-1", "1-1.5", "1.5-2", "2-3",
                  "3-4", "4-5", "5-10", ">10")),
    zero_bgc = total_bgcs_master == 0)

size_summary <- density_df %>%
  filter(!is.na(size_bin)) %>%
  group_by(size_bin) %>%
  summarise(n = n(), pct_zero = round(mean(zero_bgc) * 100, 1),
            mean_bgcs = round(mean(total_bgcs_master), 2), .groups = "drop")

max_bgcs <- max(size_summary$mean_bgcs)
scale_factor <- 100 / max_bgcs

pA <- ggplot(size_summary, aes(x = size_bin)) +
  geom_col(aes(y = pct_zero), fill = "#D6604D", alpha = 0.85,
           width = 0.65, color = "white", linewidth = 0.2) +
  geom_line(aes(y = mean_bgcs * scale_factor, group = 1),
            color = "#2166AC", linewidth = 1.2) +
  geom_point(aes(y = mean_bgcs * scale_factor), color = "#2166AC", size = 3) +
  geom_text(aes(y = pmin(pct_zero + 4, 98), label = paste0("n=", n)),
            size = 2.2, color = "grey35", fontface = "italic") +
  scale_y_continuous(
    name = "Genomes with zero BGCs (%)",
    limits = c(0, 103), breaks = seq(0, 100, 20),
    labels = function(x) paste0(x, "%"), expand = c(0, 0),
    sec.axis = sec_axis(transform = ~ . / scale_factor,
      name = "Mean BGC count per genome",
      breaks = pretty(c(0, max_bgcs), n = 6))) +
  xlab("Genome size (Mb)") +
  labs(tag = "A",
       subtitle = "Red bars: % genomes with zero BGCs (left axis)  |  Blue line: mean BGC count (right axis)") +
  pub_theme +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 7.5),
    axis.title.y.left  = element_text(color = "#D6604D", face = "bold"),
    axis.text.y.left   = element_text(color = "#D6604D"),
    axis.title.y.right = element_text(color = "#2166AC", face = "bold",
                                      angle = 270),
    axis.text.y.right  = element_text(color = "#2166AC"))

# Panel B — Cultured vs MAG boxplot
# FIX: Stats go in subtitle (always outside plot). Legend goes to "right" (outside plot).
# This way NOTHING inside the plot area overlaps.
pB <- ggplot(master, aes(x = Ecosystem_f, y = bgc_per_mb,
                          fill = Culture_Status)) +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.3, width = 0.6,
               position = position_dodge(0.75), linewidth = 0.3,
               outlier.color = "grey60") +
  scale_fill_manual(
    values = c("Cultured" = "#2166AC", "Uncultured/MAG" = "#D6604D"),
    name = "Genome origin") +
  scale_y_continuous(name = "BGC density (BGCs/Mb)",
    trans = "log1p", breaks = c(0, 1, 2, 3, 5, 8, 12),
    labels = c("0", "1", "2", "3", "5", "8", "12")) +
  xlab(NULL) +
  labs(tag = "B",
       subtitle = expression(
         "MAGs show lower BGC density (Wilcoxon " * italic(p) *
         " = 1.83e-8; median cultured: 1.15, median MAG: 0.77 BGCs/Mb)")) +
  pub_theme +
  theme(
    axis.text.x = element_text(angle = 40, hjust = 1, size = 7.5),
    # Legend OUTSIDE the plot on the right — can never overlap anything
    legend.position = "right")

fig_out <- pA / pB + plot_layout(heights = c(1, 1.15))
ggsave(file.path(FIGS, "Figure5_methodological.pdf"),
       fig_out, width = 200, height = 220, units = "mm", dpi = 600)
ggsave(file.path(FIGS, "Figure5_methodological.png"),
       fig_out, width = 200, height = 220, units = "mm", dpi = 600)
ggsave(file.path(FIGS, "Figure5_methodological.tiff"),
       fig_out, width = 200, height = 220, units = "mm", dpi = 600,
       compression = "lzw")
cat("Figure 5 saved.\n")
