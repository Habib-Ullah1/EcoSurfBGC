#!/usr/bin/env Rscript
library(dplyr)
library(readr)
library(tidyr)

BASE <- "/data/habib/EcoSurfBGC"
master <- read_csv(file.path(BASE, "results/integrated_master_table.csv"), show_col_types=FALSE)

key_classes <- c("NRPS","NRPS-like","T1PKS","T3PKS","PKS-like",
                 "terpene","RiPP-like","lanthipeptide-class-i",
                 "lanthipeptide-class-ii","lassopeptide","thiopeptide",
                 "ranthipeptide","betalactone","arylpolyene",
                 "NI-siderophore","ectoine","hserlactone","NAPAA")

key_classes <- key_classes[key_classes %in% names(master)]
cat("BGC classes for co-occurrence:", length(key_classes), "\n")
cat("Classes:", paste(key_classes, collapse=", "), "\n\n")

bgc_bin <- master %>%
  select(IMG.Genome.ID, Ecosystem_Category, C_Phylum,
         Biosurf_Total, all_of(key_classes)) %>%
  mutate(across(all_of(key_classes), ~ifelse(. > 0, 1, 0)))

cat("Genomes:", nrow(bgc_bin), "\n")

# Pairwise Jaccard co-occurrence
cat("\n=== Top 25 co-occurring BGC class pairs (Jaccard) ===\n")
cooc_pairs <- data.frame()
for(i in 1:(length(key_classes)-1)) {
  for(j in (i+1):length(key_classes)) {
    a <- bgc_bin[[key_classes[i]]]
    b <- bgc_bin[[key_classes[j]]]
    intersection <- sum(a==1 & b==1, na.rm=TRUE)
    union_ab     <- sum(a==1 | b==1, na.rm=TRUE)
    jaccard      <- if(union_ab > 0) round(intersection/union_ab, 3) else 0
    cooc_pairs   <- rbind(cooc_pairs, data.frame(
      class1=key_classes[i], class2=key_classes[j],
      jaccard=jaccard, n_cooccur=intersection, n_either=union_ab
    ))
  }
}
cooc_pairs <- cooc_pairs %>% arrange(desc(jaccard))
print(head(cooc_pairs, 25))

# BGC class richness per genome
bgc_bin <- bgc_bin %>%
  mutate(n_bgc_classes = rowSums(across(all_of(key_classes)), na.rm=TRUE))

cat("\n=== BGC class richness distribution ===\n")
print(bgc_bin %>% count(n_bgc_classes) %>%
      mutate(pct=round(n/sum(n)*100,1)) %>% arrange(n_bgc_classes))

cat("\n=== Mean BGC class richness by ecosystem ===\n")
eco_v <- bgc_bin %>%
  group_by(Ecosystem_Category) %>%
  summarise(
    n=n(),
    mean_classes    = round(mean(n_bgc_classes),2),
    median_classes  = median(n_bgc_classes),
    pct_multi       = round(mean(n_bgc_classes >= 3)*100,1),
    pct_NRPS_T1PKS  = round(mean(NRPS==1 & T1PKS==1, na.rm=TRUE)*100,1),
    pct_NRPS_terp   = round(mean(NRPS==1 & terpene==1, na.rm=TRUE)*100,1),
    pct_NRPS_RiPP   = round(mean(NRPS==1 & `RiPP-like`==1, na.rm=TRUE)*100,1),
    .groups="drop"
  ) %>% arrange(desc(mean_classes))
print(eco_v, n=25)

cat("\n=== Correlation P1 known genes vs BGC class richness ===\n")
cr <- cor.test(bgc_bin$Biosurf_Total, bgc_bin$n_bgc_classes,
               method="spearman", exact=FALSE)
cat(sprintf("Spearman rho=%.3f, p=%.2e\n", cr$estimate, cr$p.value))

# Co-occurrence matrix as wide table
cooc_matrix <- matrix(0, length(key_classes), length(key_classes),
  dimnames=list(key_classes, key_classes))
for(r in 1:nrow(cooc_pairs)) {
  i <- cooc_pairs$class1[r]; j <- cooc_pairs$class2[r]
  cooc_matrix[i,j] <- cooc_pairs$jaccard[r]
  cooc_matrix[j,i] <- cooc_pairs$jaccard[r]
}
diag(cooc_matrix) <- 1

write_csv(as.data.frame(cooc_matrix) %>% tibble::rownames_to_column("class"),
  file.path(BASE, "results/bgc_cooccurrence_matrix.csv"))
write_csv(cooc_pairs,
  file.path(BASE, "results/bgc_cooccurrence_pairs.csv"))
write_csv(eco_v,
  file.path(BASE, "results/bgc_class_versatility_by_ecosystem.csv"))
write_csv(bgc_bin %>% select(IMG.Genome.ID, Ecosystem_Category,
  C_Phylum, Biosurf_Total, n_bgc_classes),
  file.path(BASE, "results/genome_bgc_class_count.csv"))
cat("\nA4 complete. Files saved.\n")
