#!/usr/bin/env Rscript
# Sensitivity analysis: exclude MAGs and recompute all key metrics
# Compares full dataset vs cultured-isolates-only rankings and statistics

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr)
  library(ggplot2); library(patchwork); library(stringr)
  library(forcats); library(scales)
})

BASE <- "/data/habib/EcoSurfBGC"
FIGS <- file.path(BASE, "figures")
SUPP <- file.path(BASE, "supplementary")

cat("Loading data...\n")
master <- read_csv(file.path(BASE,"results/integrated_master_table_v3.csv"),
                   show_col_types=FALSE) %>%
  mutate(IMG.Genome.ID=as.character(IMG.Genome.ID))

cat(sprintf("Full dataset: %d genomes\n", nrow(master)))
cat(sprintf("MAGs: %d (%.1f%%)\n",
    sum(master$Culture_Status=="Uncultured/MAG",na.rm=TRUE),
    mean(master$Culture_Status=="Uncultured/MAG",na.rm=TRUE)*100))

# ── Two datasets ──────────────────────────────────────────
full_data     <- master %>% filter(!is.na(Ecosystem_Category))
no_mag_data   <- master %>%
  filter(!is.na(Ecosystem_Category),
         Culture_Status == "Cultured")

cat(sprintf("No-MAG dataset: %d genomes\n\n", nrow(no_mag_data)))

# ── Function: compute ecosystem summary ──────────────────
eco_summary <- function(df, label) {
  df %>%
    group_by(Ecosystem_Category) %>%
    summarise(
      N                  = n(),
      mean_BGC_per_Mb    = round(mean(bgc_per_mb,    na.rm=TRUE), 3),
      mean_total_BGCs    = round(mean(total_bgcs,    na.rm=TRUE), 1),
      mean_pct_novel     = round(mean(pct_novel,     na.rm=TRUE), 1),
      mean_n_bgc_classes = round(mean(n_bgc_classes, na.rm=TRUE), 2),
      pct_MAG            = round(mean(Culture_Status=="Uncultured/MAG",
                                      na.rm=TRUE)*100, 1),
      .groups="drop"
    ) %>%
    arrange(desc(mean_BGC_per_Mb)) %>%
    mutate(
      Rank    = row_number(),
      Dataset = label
    )
}

full_eco   <- eco_summary(full_data,   "Full dataset")
no_mag_eco <- eco_summary(no_mag_data, "Cultured isolates only")

# ── Merge for comparison ──────────────────────────────────
comparison <- full_eco %>%
  select(Ecosystem_Category,
         Full_N            = N,
         Full_Rank         = Rank,
         Full_BGC_density  = mean_BGC_per_Mb,
         Full_pct_novel    = mean_pct_novel,
         Full_pct_MAG      = pct_MAG) %>%
  left_join(
    no_mag_eco %>%
      select(Ecosystem_Category,
             NoMAG_N           = N,
             NoMAG_Rank        = Rank,
             NoMAG_BGC_density = mean_BGC_per_Mb,
             NoMAG_pct_novel   = mean_pct_novel),
    by = "Ecosystem_Category"
  ) %>%
  mutate(
    Rank_change         = Full_Rank - NoMAG_Rank,
    Density_change      = round(NoMAG_BGC_density - Full_BGC_density, 3),
    Novelty_change      = round(NoMAG_pct_novel   - Full_pct_novel,   1),
    Rank_shift_category = case_when(
      abs(Rank_change) == 0 ~ "No change",
      abs(Rank_change) <= 2 ~ "Minor (1-2 ranks)",
      abs(Rank_change) <= 4 ~ "Moderate (3-4 ranks)",
      TRUE                  ~ "Major (5+ ranks)"
    )
  ) %>%
  arrange(Full_Rank)

cat("=== SENSITIVITY ANALYSIS: FULL vs NO-MAG DATASET ===\n\n")
cat("Ecosystem rankings and metric changes:\n")
print(comparison %>%
  select(Ecosystem_Category, Full_Rank, NoMAG_Rank,
         Rank_change, Full_BGC_density, NoMAG_BGC_density,
         Density_change, Full_pct_novel, NoMAG_pct_novel,
         Novelty_change, Full_pct_MAG),
  n=17, width=120)

cat("\n=== KEY STATISTICS ===\n")
cat(sprintf("Ecosystems with NO rank change:        %d\n",
    sum(comparison$Rank_change==0)))
cat(sprintf("Ecosystems with minor shift (1-2):     %d\n",
    sum(abs(comparison$Rank_change) %in% 1:2)))
cat(sprintf("Ecosystems with moderate shift (3-4):  %d\n",
    sum(abs(comparison$Rank_change) %in% 3:4)))
cat(sprintf("Ecosystems with major shift (5+):      %d\n",
    sum(abs(comparison$Rank_change) >= 5)))

cat("\nEcosystems affected most by MAG removal:\n")
print(comparison %>%
  filter(abs(Rank_change)>0 | abs(Novelty_change)>5) %>%
  select(Ecosystem_Category, Full_Rank, NoMAG_Rank,
         Rank_change, Density_change, Novelty_change, Full_pct_MAG) %>%
  arrange(desc(abs(Novelty_change))), width=100)

cat("\n=== OVERALL DATASET STATISTICS ===\n")
cat(sprintf("Full dataset    - mean BGC density: %.3f BGCs/Mb | mean novelty: %.1f%%\n",
    mean(full_data$bgc_per_mb,na.rm=TRUE),
    mean(full_data$pct_novel,na.rm=TRUE)))
cat(sprintf("No-MAG dataset  - mean BGC density: %.3f BGCs/Mb | mean novelty: %.1f%%\n",
    mean(no_mag_data$bgc_per_mb,na.rm=TRUE),
    mean(no_mag_data$pct_novel,na.rm=TRUE)))

# ── Kruskal-Wallis on no-MAG dataset ─────────────────────
kw_density <- kruskal.test(bgc_per_mb ~ Ecosystem_Category,
                            data=no_mag_data)
kw_novelty <- kruskal.test(pct_novel  ~ Ecosystem_Category,
                            data=no_mag_data)

cat(sprintf("\nNo-MAG Kruskal-Wallis BGC density: chi2=%.2f, p=%.2e\n",
    kw_density$statistic, kw_density$p.value))
cat(sprintf("No-MAG Kruskal-Wallis novelty:     chi2=%.2f, p=%.2e\n",
    kw_novelty$statistic, kw_novelty$p.value))
cat(sprintf("Full dataset Kruskal-Wallis BGC density p = 7.64e-29 (reference)\n"))

# ── Save comparison table ─────────────────────────────────
write_csv(comparison,
  file.path(SUPP,"Table_S7_sensitivity_analysis_no_MAGs.csv"))
cat("\nSaved: Table_S7_sensitivity_analysis_no_MAGs.csv\n")

# ── FIGURE: Rank comparison full vs no-MAG ───────────────
eco_short <- function(x) str_replace_all(x,
  c("Animal-Associated [(]Non-Mammalian[)]"="Animal-Assoc.",
    "Microbial/Fermentation"="Microbial/Ferm.",
    "Oil-Related Environments"="Oil-Related",
    "Extreme Environments"="Extreme",
    "Saline Environments"="Saline",
    "Built-Environment"="Built-Env.",
    "Plant-Associated"="Plant-Assoc.",
    "Food-Products"="Food-Prod.",
    "Industrial/Waste"="Ind./Waste"))

plot_data <- comparison %>%
  mutate(Ecosystem_short=eco_short(Ecosystem_Category))

# Panel A: rank comparison slope chart
slope_df <- plot_data %>%
  select(Ecosystem_short, Full_Rank, NoMAG_Rank,
         Rank_change, Full_pct_MAG) %>%
  mutate(
    shift_dir=case_when(
      Rank_change > 0 ~ "Improved (no-MAG)",
      Rank_change < 0 ~ "Declined (no-MAG)",
      TRUE            ~ "No change"),
    abs_change=abs(Rank_change))

# Space labels
space_fn <- function(df, ycol) {
  df <- df %>% arrange(get(ycol))
  for(i in 2:nrow(df))
    if(df[[ycol]][i]-df[[ycol]][i-1] < 0.6)
      df[[ycol]][i] <- df[[ycol]][i-1]+0.6
  df
}

left_labs  <- space_fn(slope_df %>%
  arrange(Full_Rank) %>% mutate(y_lab=Full_Rank), "y_lab")
right_labs <- space_fn(slope_df %>%
  arrange(NoMAG_Rank) %>% mutate(y_lab=NoMAG_Rank), "y_lab")

col_imp <- "#1B7837"; col_dec <- "#C51B7D"; col_nc <- "grey55"

pA <- ggplot(slope_df) +
  geom_segment(aes(x=1,xend=2,y=Full_Rank,yend=NoMAG_Rank,
    color=shift_dir,linewidth=abs_change+0.1),alpha=0.75) +
  geom_point(aes(x=1,y=Full_Rank,color=shift_dir),size=2) +
  geom_point(aes(x=2,y=NoMAG_Rank,color=shift_dir),size=2) +
  geom_text(data=left_labs,
    aes(x=0.88,y=y_lab,label=Ecosystem_short,color=shift_dir),
    hjust=1,size=2.2,fontface="italic",show.legend=FALSE) +
  geom_segment(data=left_labs,
    aes(x=0.89,xend=0.98,y=y_lab,yend=Full_Rank),
    color="grey75",linewidth=0.15) +
  geom_text(data=right_labs,
    aes(x=2.12,y=y_lab,label=Ecosystem_short,color=shift_dir),
    hjust=0,size=2.2,fontface="italic",show.legend=FALSE) +
  geom_segment(data=right_labs,
    aes(x=2.02,xend=2.11,y=NoMAG_Rank,yend=y_lab),
    color="grey75",linewidth=0.15) +
  scale_y_reverse(name="BGC density rank",breaks=1:17,
    limits=c(18,0.3),expand=c(0,0)) +
  scale_x_continuous(limits=c(-0.5,3.8),breaks=c(1,2),
    labels=c("Full dataset\n(with MAGs)",
             "Cultured only\n(no MAGs)")) +
  scale_color_manual(name="Rank shift",
    values=c("Improved (no-MAG)"=col_imp,
             "Declined (no-MAG)"=col_dec,
             "No change"=col_nc),
    guide=guide_legend(override.aes=list(linewidth=1.2,size=2))) +
  scale_linewidth_continuous(range=c(0.3,2),guide="none") +
  labs(tag="A",
    subtitle="Rank shifts in BGC density after excluding MAG genomes") +
  theme_classic(base_size=9) +
  theme(
    axis.text.x=element_text(size=8,face="bold"),
    axis.title.x=element_blank(),
    axis.ticks.x=element_blank(),
    panel.grid.major.x=element_blank(),
    panel.grid.major.y=element_line(color="grey92",linewidth=0.25),
    legend.position=c(0.97,0.97),
    legend.justification=c(1,1),
    legend.background=element_rect(fill=alpha("white",0.9),
      color="grey70",linewidth=0.3),
    plot.tag=element_text(face="bold",size=11),
    plot.subtitle=element_text(size=7.5,color="grey40"))

# Panel B: novelty rate comparison scatter
pB <- ggplot(plot_data,
  aes(x=Full_pct_novel, y=NoMAG_pct_novel,
      color=Full_pct_MAG, size=Full_N,
      label=Ecosystem_short)) +
  geom_abline(slope=1,intercept=0,linetype="dashed",
    color="grey60",linewidth=0.6) +
  geom_point(alpha=0.85) +
  ggrepel::geom_text_repel(size=2.3,show.legend=FALSE,
    segment.size=0.3,segment.color="grey60",
    box.padding=0.4,max.overlaps=20) +
  scale_color_gradient2(low="#2166AC",mid="#F7F7F7",
    high="#D6604D",midpoint=5,
    name="% MAG in\nfull dataset") +
  scale_size_continuous(name="N genomes",range=c(2,8)) +
  scale_x_continuous(name="% Novel BGCs — full dataset") +
  scale_y_continuous(name="% Novel BGCs — cultured isolates only") +
  annotate("text",x=44,y=87,label="Above line: novelty increases\nwhen MAGs removed",
    size=2.4,color="grey40",hjust=0,fontface="italic") +
  annotate("text",x=44,y=43,label="Below line: novelty decreases\nwhen MAGs removed",
    size=2.4,color="grey40",hjust=0,fontface="italic") +
  labs(tag="B",
    subtitle="Points above dashed line: MAG removal increases apparent novelty") +
  theme_classic(base_size=9) +
  theme(
    axis.title=element_text(face="bold",size=9),
    axis.text=element_text(size=8),
    legend.title=element_text(face="bold",size=8),
    legend.text=element_text(size=7),
    panel.grid.major=element_line(color="grey92",linewidth=0.25),
    plot.tag=element_text(face="bold",size=11),
    plot.subtitle=element_text(size=7.5,color="grey40"))

# Panel C: BGC density comparison
density_long <- comparison %>%
  mutate(Ecosystem_short=eco_short(Ecosystem_Category)) %>%
  select(Ecosystem_short,
         `Full dataset`=Full_BGC_density,
         `Cultured only`=NoMAG_BGC_density) %>%
  pivot_longer(-Ecosystem_short,
    names_to="Dataset",values_to="BGC_density") %>%
  mutate(Ecosystem_short=fct_reorder(Ecosystem_short,
    BGC_density*(Dataset=="Full dataset"),.fun=max))

pC <- ggplot(density_long,
  aes(x=Ecosystem_short,y=BGC_density,
      fill=Dataset,color=Dataset)) +
  geom_line(aes(group=Ecosystem_short),color="grey70",linewidth=0.4) +
  geom_point(size=2.5,shape=21,stroke=0.5,color="white") +
  scale_fill_manual(values=c("Full dataset"="#3498DB",
                              "Cultured only"="#E74C3C"),
    name="Dataset") +
  scale_color_manual(values=c("Full dataset"="#3498DB",
                               "Cultured only"="#E74C3C"),
    name="Dataset",guide="none") +
  scale_y_continuous(name="Mean BGC density (BGCs/Mb)",
    expand=expansion(mult=c(0.05,0.1))) +
  xlab(NULL) +
  labs(tag="C",
    subtitle="Blue = full dataset; Red = cultured isolates only; connected by grey line") +
  theme_classic(base_size=9) +
  theme(
    axis.text.x=element_text(angle=40,hjust=1,size=7.5),
    axis.title=element_text(face="bold",size=9),
    legend.title=element_text(face="bold",size=8),
    legend.text=element_text(size=7),
    panel.grid.major.y=element_line(color="grey92",linewidth=0.25),
    plot.tag=element_text(face="bold",size=11),
    plot.subtitle=element_text(size=7.5,color="grey40"))

# Combine
fig_sensitivity <- (pA + pB) / pC +
  plot_layout(heights=c(1.3,0.9)) +
  plot_annotation(theme=theme(plot.margin=margin(5,5,5,5)))

ggsave(file.path(FIGS,"FigureS3_sensitivity_analysis.pdf"),
  fig_sensitivity,width=190,height=200,units="mm",dpi=600)
ggsave(file.path(FIGS,"FigureS3_sensitivity_analysis.png"),
  fig_sensitivity,width=190,height=200,units="mm",dpi=600)

cat("Figure S3 saved.\n")
cat("\n=== SENSITIVITY ANALYSIS COMPLETE ===\n")
