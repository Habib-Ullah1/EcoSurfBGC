#!/usr/bin/env Rscript
# 02_subset_genomes_base.R
# EcoSurfBGC — Phase 2: Strategic Genome Subsetting (base R)
# Author: Habib Ullah | Version: 1.0 | March 2026
set.seed(42)

BASE   <- '/data/habib/EcoSurfBGC'
META   <- file.path(BASE, 'metadata',
          'Merged_raw_CheckM2_Filtered_JGI-info-table_newTax_ProfileAdded_acc_removed.txt')
OUTDIR <- file.path(BASE, 'genome_lists')

target_genes <- c('fen_pps','fab','itu','rhl','rml','wza','wzb','wzc',
                  'wzy','lic_srf','srf','adh','alg','arf','pso','est','psw','cyp')

ecosystem_categories <- list(
  'Mammalian'=c('Mammals','Mammals: Human','Humans','Human','Humans-Associated','human-associated','Mammalian','Livestock-associated','livestock-associated','livestock','Cats','murine cecal','Camel','Camel-Associated','Human skin','Host, Human skin','Fetus'),
  'Terrestrial'=c('Terrestrial','Soil','Desert','Composting','Agricultural field','Urban waste','Vermicompost','Peat moss','Plant litter','Wood','Rock-dwelling (subaerial biofilms)','Rock-dwelling (endoliths)','Monument','Cement wall','farm-isolated','Metal contaminated soil'),
  'Marine'=c('Marine','marine','Aquatic, Marine','Bay','Sargassum','Marine Mat','Mat','Sediment','Deep subsurface','Subsurface','Porifera','Cnidaria','Ascidians','Tunicates','Echinodermata','Bivalves','Oyster','Seafood','Seafood product','Fish products','fish-associated','Aquaculture'),
  'Plant-Associated'=c('Plants','Plant-Associated','Phyllosphere','Rhizosphere','Endosphere','Rhizome','Seeds','Whole plant body','Host, Plant root, Plants','Plant products','Nodule','Fruiting body','Mycelium','Lichen','Roots','Mesocosm','Feedstock','Spore'),
  'Microbial/Fermentation'=c('Microbial/Fermentation','Microbial','Biofilm','Fermentation','Fermented beverages','Fermented vegetables','Silage fermentation','Fermented food','Fermented seafood','Fermentation starter','Fermentation cellar','Biotransformation','Lab synthesis','Lab enrichment','Defined media','Culture media','Continuous culture','Bioreactor','Photobioreactor (PBR)','Anaerobic digester'),
  'Freshwater'=c('Freshwater','Fresh water','Aquatic, Fresh water','Lake Mat','Floodplain','Drinking water','drinking water'),
  'Food-Products'=c('Food production','Food-Products','Food-Associated','Food sample','Food','Food ','Food waste','Agri-food','Dairy products','Meat products','Bread production','Baby formula','Spices','Nuts','Egg products','Vegetable','Dairy processing facility','Beverages','Fermented food','Fish products'),
  'Animal-Associated (Non-Mammalian)'=c('Host-associated','Animal-associated','Birds','Fish','Reptilia','Amphibia','Arthropoda: Insects','Arthropoda: Crustaceans','Arthropoda: Chelicerates','Arthropoda: Myriapoda','Nematoda','Annelida','Mollusca','Protozoa','Protists','Amoebozoa','Dinoflagellata','Dinoflagellates','Invertebrates','Caenorhabditis elegans','Cockroach','Paenibacillus larvae','Larva','Larvae','Larva: Nauplius','Argopecten purpuratus larvae','Ootheca/Egg mass','Nymph/Instar','Poultry','poultry','poultry sources','Poultry-Associated','Chicken','chicken','Duck','Skeletal system','Muscular system','Fat body','Shell','Bryozoans','Cephalochordata','Nest','Sponge','Rumen','Abdomen','Genome of termite-associated Trabulsiella odontotermitis strain'),
  'Aquatic (Unspecified)'=c('Aquatic','Mat, Aquatic','Aquatic, Biofilm'),
  'Industrial/Waste'=c('Industrial production','Industrial wastewater','Industrial waste','Industrial/Waste','Engineered','Engineered product','Genetically modified','Chemical products','Tailings pond','Bagasse','Household waste','Zoo waste','Animal waste','Solid waste','Landfill','Agricultural waste','Brown waste','Cellulose associated waste'),
  'Wastewater'=c('Wastewater','Sewage','WWTP','Activated Sludge','Activated sludge','UASB (Upflow anaerobic sludge blanket)','MBR (Membrane bioreactor)','Water treatment plant','Drinking water treatment plant','Percolator','Nutrient removal','Influent','Sludge','Drinking Water Filter'),
  'Built-Environment'=c('Built environment','International Space Station','Spacecraft Assembly Cleanrooms','House','Building','City','Vivarium'),
  'Geological'=c('Geologic','Geological','Volcanic','Cave','halite crust','Subsurface'),
  'Saline Environments'=c('Saline','Non-marine Saline and Alkaline','Saline_environments','High-salinity/high-pH','Salt'),
  'Algae'=c('Algae','Green algae','Brown Algae','Red algae','Microalgae','Diatoms'),
  'Air'=c('Air','Outdoor Air','Indoor Air'),
  'Fungi'=c('Fungi'),
  'Clinical'=c('Clinical','clinical','Hospital','Hospital-Associated','Digestive system','Circulatory system','Respiratory system','Urinary system','Integumentary system','Reproductive system','Nervous system','Lymphatic system','Excretory system','Visual system','Auditory/Hearing system','Sensory organs','Gastrointestinal tract','Intestinal tract','Abdominal/Peritoneal cavity','Unspecified system','Multiple systems','Head','Skin','Integument','Tissue','Whole body','Remains','Acinetobacter larvae BRTC-1 Genome sequencing'),
  'Oil-Related Environments'=c('Oil reservoir','Oil refinery','Hydrocarbon'),
  'Extreme Environments'=c('Thermal springs','Hot spring','Hot Lake','Acidic','Anaerobic','Aerobic','Mud volcano')
)

categorize <- function(x) {
  for (cat in names(ecosystem_categories)) {
    if (!is.na(x) && x %in% ecosystem_categories[[cat]]) return(cat)
  }
  return('EXCLUDE')
}

# ── Load data ─────────────────────────────────────────────
cat('Loading metadata...\n')
df <- read.delim(META, sep='\t', stringsAsFactors=FALSE, check.names=FALSE)
for (g in target_genes) {
  if (g %in% colnames(df)) df[[g]] <- ifelse(as.numeric(df[[g]]) > 0, 1, 0)
}
df$Biosurf_Total <- rowSums(df[, target_genes], na.rm=TRUE)
cat('Total genomes:', nrow(df), '\n')

# ── LAYER 1: Quality filter ───────────────────────────────
cat('\n[Layer 1] Quality filtering (>=90% complete, <5% contamination)...\n')
L1 <- df[!is.na(df$CheckM2.Completeness) & !is.na(df$CheckM2.Contamination) &
         df$CheckM2.Completeness >= 90 & df$CheckM2.Contamination < 5, ]
cat('After quality filter:', nrow(L1), '\n')
write.csv(data.frame(IMG.Genome.ID=L1$IMG.Genome.ID),
          file.path(OUTDIR,'quality_filtered.csv'), row.names=FALSE)

# ── LAYER 2: ANI de-replication ───────────────────────────
cat('\n[Layer 2] ANI cluster de-replication...\n')
singletons <- L1[is.na(L1$ANI.Cluster.ID) | L1$ANI.Cluster.Type == 'singleton', ]
clustered  <- L1[!is.na(L1$ANI.Cluster.ID) & L1$ANI.Cluster.Type != 'singleton', ]

# Keep one representative per cluster (highest completeness)
cluster_ids <- unique(clustered$ANI.Cluster.ID)
reps <- do.call(rbind, lapply(cluster_ids, function(cid) {
  sub <- clustered[clustered$ANI.Cluster.ID == cid, ]
  sub[which.max(sub$CheckM2.Completeness), ]
}))

L2 <- rbind(singletons, reps)
cat('After ANI de-replication:', nrow(L2), '\n')
write.csv(data.frame(IMG.Genome.ID=L2$IMG.Genome.ID),
          file.path(OUTDIR,'ani_representatives.csv'), row.names=FALSE)

# ── LAYER 3: Ecosystem categorization and stratified sampling ──
cat('\n[Layer 3] Ecosystem categorization...\n')
L2$Ecosystem_Category <- sapply(L2$Ecosystem.Type, categorize)
L2_clean <- L2[L2$Ecosystem_Category != 'EXCLUDE', ]
cat('After excluding Unclassified/Uncertain:', nrow(L2_clean), '\n')

# Sampling targets: proportional with floor=50, cap=400
eco_counts <- table(L2_clean$Ecosystem_Category)
total_clean <- nrow(L2_clean)
TARGET_TOTAL <- 3000
FLOOR <- 50
CAP   <- 400

cat('\nEcosystem counts and sampling targets:\n')
targets <- data.frame(
  Ecosystem_Category = names(eco_counts),
  available          = as.integer(eco_counts)
)
targets$prop_target <- round(targets$available / total_clean * TARGET_TOTAL)
targets$target      <- pmax(pmin(targets$prop_target, CAP), pmin(FLOOR, targets$available))
print(targets, row.names=FALSE)

# ── LAYER 4: Gene profile diversity sampling ───────────────
cat('\n[Layer 4] Gene profile diversity sampling (50% high, 30% medium, 20% low)...\n')

sample_eco <- function(eco_name, n_target) {
  eco_df <- L2_clean[L2_clean$Ecosystem_Category == eco_name, ]
  n_high   <- round(n_target * 0.50)
  n_medium <- round(n_target * 0.30)
  n_low    <- n_target - n_high - n_medium
  high   <- eco_df[eco_df$Biosurf_Total >= 7, ]
  medium <- eco_df[eco_df$Biosurf_Total >= 3 & eco_df$Biosurf_Total < 7, ]
  low    <- eco_df[eco_df$Biosurf_Total >= 1 & eco_df$Biosurf_Total < 3, ]
  s_high   <- if (nrow(high)   > 0) high[sample(nrow(high),   min(n_high,   nrow(high))),   ] else high
  s_medium <- if (nrow(medium) > 0) medium[sample(nrow(medium), min(n_medium, nrow(medium))), ] else medium
  s_low    <- if (nrow(low)    > 0) low[sample(nrow(low),    min(n_low,    nrow(low))),    ] else low
  rbind(s_high, s_medium, s_low)
}

final_list <- lapply(seq_len(nrow(targets)), function(i) {
  sample_eco(targets$Ecosystem_Category[i], targets$target[i])
})
final_subset <- do.call(rbind, final_list)

cat('\nFinal subset size:', nrow(final_subset), 'genomes\n')
cat('\nFinal subset by ecosystem:\n')
print(sort(table(final_subset$Ecosystem_Category), decreasing=TRUE))
cat('\nGene count distribution in subset:\n')
print(table(final_subset$Biosurf_Total))

# ── Save outputs ──────────────────────────────────────────
write.csv(final_subset,
          file.path(OUTDIR,'final_subset_full_metadata.csv'), row.names=FALSE)
writeLines(as.character(final_subset$IMG.Genome.ID),
           file.path(OUTDIR,'final_subset_ids.txt'))

cat('\nPhase 2 complete.\n')
cat('Final genome IDs saved to: genome_lists/final_subset_ids.txt\n')
cat('Full metadata saved to:    genome_lists/final_subset_full_metadata.csv\n')
