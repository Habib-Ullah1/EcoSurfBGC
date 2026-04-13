suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(readr); library(tidyr)
  library(forcats); library(stringr); library(scales)
})
BASE <- "/data/habib/EcoSurfBGC"
SUPP <- file.path(BASE, "supplementary")

set.seed(42)
N_PERM <- 100

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

eco_colors <- c(
  "Plant-Assoc." = "#2E8B57", "Terrestrial" = "#8B4513",
  "Marine" = "#1E90FF", "Freshwater" = "#00BFFF",
  "Clinical" = "#DC143C", "Wastewater" = "#9932CC",
  "Air" = "#87CEEB", "Algae" = "#3CB371",
  "Industrial/Waste" = "#FF8C00", "Built-Env." = "#708090",
  "Extreme" = "#FF4500", "Saline" = "#20B2AA",
  "Geological" = "#8B6914", "Food-Prod." = "#DAA520",
  "Microbial/Ferm." = "#9370DB", "Animal-Assoc." = "#FF69B4",
  "Oil-Related" = "#4A4A4A")

master <- read_csv(file.path(BASE, "results/integrated_master_table_v2.csv"),
                   show_col_types = FALSE)

all_cols <- names(master)
start_idx <- which(all_cols == "total_bgcs") + 1
end_idx   <- which(all_cols == "Culture_Status") - 1
bgc_cols  <- all_cols[start_idx:end_idx]
bgc_cols  <- bgc_cols[sapply(bgc_cols, function(x) is.numeric(master[[x]]))]

bgc_long <- master %>%
  filter(!is.na(Ecosystem_Category)) %>%
  select(IMG.Genome.ID, Ecosystem_Category, all_of(bgc_cols)) %>%
  pivot_longer(cols = all_of(bgc_cols), names_to = "product", values_to = "count") %>%
  filter(!is.na(count), count > 0)

ecosystems_list <- unique(master$Ecosystem_Category)
ecosystems_list <- ecosystems_list[!is.na(ecosystems_list)]

rarefaction_data <- list()
for (eco in ecosystems_list) {
  eco_genomes <- master %>%
    filter(Ecosystem_Category == eco) %>% pull(IMG.Genome.ID)
  if (length(eco_genomes) < 3) next
  max_n <- length(eco_genomes)
  steps <- unique(c(seq(1, min(10, max_n), 1),
                    seq(10, max_n, max(1, floor(max_n / 20))), max_n))
  perm_mat <- sapply(1:N_PERM, function(p) {
    shuffled <- sample(eco_genomes)
    sapply(steps, function(n)
      bgc_long %>% filter(IMG.Genome.ID %in% shuffled[1:n]) %>%
        pull(product) %>% n_distinct())
  })
  if (is.vector(perm_mat))
    perm_mat <- matrix(perm_mat, nrow = length(steps))
  eco_label <- eco_short(eco)
  rarefaction_data[[eco_label]] <- data.frame(
    Ecosystem = eco_label, N_genomes = steps,
    Mean_BGC_classes = rowMeans(perm_mat),
    SE = apply(perm_mat, 1, sd) / sqrt(N_PERM))
}

raref_df <- bind_rows(rarefaction_data)
final_order <- raref_df %>%
  group_by(Ecosystem) %>% slice_max(N_genomes, n = 1) %>%
  arrange(desc(Mean_BGC_classes)) %>% pull(Ecosystem)
raref_df <- raref_df %>%
  mutate(Ecosystem = factor(Ecosystem, levels = final_order))

pS1 <- ggplot(raref_df,
              aes(x = N_genomes, y = Mean_BGC_classes,
                  color = Ecosystem, fill = Ecosystem)) +
  geom_ribbon(aes(ymin = Mean_BGC_classes - SE,
                  ymax = Mean_BGC_classes + SE),
              alpha = 0.12, color = NA) +
  geom_line(linewidth = 0.75) +
  scale_color_manual(values = eco_colors, name = "Ecosystem") +
  scale_fill_manual(values = eco_colors, guide = "none") +
  scale_x_continuous(name = "Number of genomes sampled",
                     expand = expansion(mult = c(0, 0.02))) +
  coord_cartesian(xlim = c(1, NA)) +
  scale_y_continuous(name = "Cumulative BGC product types discovered",
                     expand = expansion(mult = c(0, 0.05))) +
  labs(subtitle = "Mean +/- SE across 100 random permutations") +
  pub_theme +
  theme(legend.position = "right") +
  guides(color = guide_legend(ncol = 1, keyheight = unit(3.5, "mm")))

ggsave(file.path(SUPP, "FigureS1_rarefaction_curves.pdf"),
       pS1, width = 190, height = 120, units = "mm", dpi = 600)
ggsave(file.path(SUPP, "FigureS1_rarefaction_curves.png"),
       pS1, width = 190, height = 120, units = "mm", dpi = 600)
ggsave(file.path(SUPP, "FigureS1_rarefaction_curves.tiff"),
       pS1, width = 190, height = 120, units = "mm", dpi = 600,
       compression = "lzw")
cat("Figure S1 saved.\n")
