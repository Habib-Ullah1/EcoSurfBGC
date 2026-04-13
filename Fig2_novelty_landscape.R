suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(readr); library(tidyr)
  library(scales); library(forcats); library(stringr); library(ggrepel)
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
    plot.margin       = margin(6, 10, 6, 6)
  )

eco_short <- function(x) str_replace_all(x,
  c("Animal-Associated [(]Non-Mammalian[)]" = "Animal-Assoc.",
    "Microbial/Fermentation" = "Microbial/Ferm.",
    "Oil-Related Environments" = "Oil-Related",
    "Extreme Environments" = "Extreme Env.",
    "Saline Environments" = "Saline Env.",
    "Built-Environment" = "Built-Env.",
    "Plant-Associated" = "Plant-Assoc."))

novelty_bgc <- read_csv(file.path(BASE, "results/bgc_mibig_novelty_annotated.csv"),
  show_col_types = FALSE) %>% mutate(genome_id = as.character(genome_id))
master <- read_csv(file.path(BASE, "results/integrated_master_table_v2.csv"),
  show_col_types = FALSE) %>% mutate(IMG.Genome.ID = as.character(IMG.Genome.ID))

eco_novelty <- novelty_bgc %>%
  left_join(master %>% select(IMG.Genome.ID, Ecosystem_Category),
            by = c("genome_id" = "IMG.Genome.ID")) %>%
  filter(!is.na(Ecosystem_Category)) %>%
  group_by(Ecosystem_Category) %>%
  summarise(
    novel    = sum(novelty_class == "novel", na.rm = TRUE),
    low_sim  = sum(novelty_class == "low_similarity", na.rm = TRUE),
    mod_sim  = sum(novelty_class == "moderate_similarity", na.rm = TRUE),
    high_sim = sum(novelty_class == "high_similarity", na.rm = TRUE),
    total    = n(), .groups = "drop") %>%
  mutate(across(novel:high_sim, ~ . / total * 100)) %>%
  arrange(desc(novel)) %>%
  mutate(Ecosystem_Category = fct_reorder(Ecosystem_Category, novel))

novelty_long <- eco_novelty %>%
  select(Ecosystem_Category, novel, low_sim, mod_sim, high_sim) %>%
  pivot_longer(-Ecosystem_Category,
               names_to = "class", values_to = "pct") %>%
  mutate(class = factor(class,
    levels = c("high_sim", "mod_sim", "low_sim", "novel"),
    labels = c("High similarity (5+ proteins)",
               "Moderate (3-4 proteins)",
               "Low similarity (1-2 proteins)",
               "Novel (no MIBiG match)")))

novelty_colors <- c(
  "Novel (no MIBiG match)"        = "#D62728",
  "Low similarity (1-2 proteins)" = "#FF7F0E",
  "Moderate (3-4 proteins)"       = "#2CA02C",
  "High similarity (5+ proteins)" = "#1F77B4")

pA <- ggplot(novelty_long,
             aes(x = Ecosystem_Category, y = pct, fill = class)) +
  geom_bar(stat = "identity", position = "stack", width = 0.72,
           color = "white", linewidth = 0.2) +
  scale_fill_manual(values = novelty_colors, name = "MIBiG similarity") +
  scale_x_discrete(labels = eco_short) +
  scale_y_continuous(name = "% of BGCs",
    expand = expansion(mult = c(0, 0.02)),
    labels = function(x) paste0(x, "%")) +
  xlab(NULL) +
  pub_theme +
  theme(
    axis.text.x = element_text(angle = 40, hjust = 1, size = 7.5),
    # OUTSIDE the plot on the right — never overlaps bars
    legend.position = "right") +
  guides(fill = guide_legend(reverse = TRUE))

ggsave(file.path(FIGS, "Figure2A_novelty_stacked_bar.pdf"),
       pA, width = 200, height = 110, units = "mm", dpi = 600)
ggsave(file.path(FIGS, "Figure2A_novelty_stacked_bar.png"),
       pA, width = 200, height = 110, units = "mm", dpi = 600)
ggsave(file.path(FIGS, "Figure2A_novelty_stacked_bar.tiff"),
       pA, width = 200, height = 110, units = "mm", dpi = 600,
       compression = "lzw")
cat("Figure 2A saved.\n")
