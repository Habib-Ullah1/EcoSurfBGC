#!/usr/bin/env Rscript
# 02_subset_genomes.R
# EcoSurfBGC Project — Phase 2: Strategic Genome Subsetting
# Author: Habib Ullah
# Version: 1.0 | March 2026

library(tidyverse)
set.seed(42)  # Reproducibility

# ── Paths ────────────────────────────────────────────────
BASE    <- '/data/habib/EcoSurfBGC'
META    <- file.path(BASE, 'metadata',
           'Merged_raw_CheckM2_Filtered_JGI-info-table_newTax_ProfileAdded_acc_removed.txt')
OUTDIR  <- file.path(BASE, 'genome_lists')

# ── Target genes ─────────────────────────────────────────
target_genes <- c('fen_pps','fab','itu','rhl','rml','wza','wzb','wzc',
                  'wzy','lic_srf','srf','adh','alg','arf','pso','est','psw','cyp')

# ── Ecosystem mapping ─────────────────────────────────────
# (Same mapping as 01_explore_metadata.R)
ecosystem_categories <- list(
  'Mammalian'=c('Mammals','Mammals: Human','Humans','Human','Humans-Associated',
    'human-associated','Mammalian','Livestock-associated','livestock-associated',
    'livestock','Cats','murine cecal','Camel','Camel-Associated','Human skin',
    'Host, Human skin','Fetus'),
  'Terrestrial'=c('Terrestrial','Soil','Desert','Composting','Agricultural field',
    'Urban waste','Vermicompost','Peat moss','Plant litter','Wood',
    'Rock-dwelling (subaerial biofilms)','Rock-dwelling (endoliths)','Monument',
    'Cement wall','farm-isolated','Metal contaminated soil'),
  'Marine'=c('Marine','marine','Aquatic, Marine','Bay','Sargassum','Marine Mat',
    'Mat','Sediment','Deep subsurface','Subsurface','Porifera','Cnidaria',
    'Ascidians','Tunicates','Echinodermata','Bivalves','Oyster','Seafood',
    'Seafood product','Fish products','fish-associated','Aquaculture'),
  'Plant-Associated'=c('Plants','Plant-Associated','Phyllosphere','Rhizosphere',
    'Endosphere','Rhizome','Seeds','Whole plant body','Host, Plant root, Plants',
    'Plant products','Nodule','Fruiting body','Mycelium','Lichen','Roots',
    'Mesocosm','Feedstock','Spore'),
  'Microbial/Fermentation'=c('Microbial/Fermentation','Microbial','Biofilm',
    'Fermentation','Fermented beverages','Fermented vegetables',
    'Silage fermentation','Fermented food','Fermented seafood',
    'Fermentation starter','Fermentation cellar','Biotransformation',
    'Lab synthesis','Lab enrichment','Defined media','Culture media',
    'Continuous culture','Bioreactor','Photobioreactor (PBR)','Anaerobic digester'),
  'Freshwater'=c('Freshwater','Fresh water','Aquatic, Fresh water','Lake Mat',
    'Floodplain','Drinking water','drinking water'),
  'Food-Products'=c('Food production','Food-Products','Food-Associated',
    'Food sample','Food','Food ','Food waste','Agri-food','Dairy products',
    'Meat products','Bread production','Baby formula','Spices','Nuts',
    'Egg products','Vegetable','Dairy processing facility','Beverages',
    'Fermented food','Fish products'),
  'Animal-Associated (Non-Mammalian)'=c('Host-associated','Animal-associated',
    'Birds','Fish','Reptilia','Amphibia','Arthropoda: Insects',
    'Arthropoda: Crustaceans','Arthropoda: Chelicerates','Arthropoda: Myriapoda',
    'Nematoda','Annelida','Mollusca','Protozoa','Protists','Amoebozoa',
    'Dinoflagellata','Dinoflagellates','Invertebrates','Caenorhabditis elegans',
    'Cockroach','Paenibacillus larvae','Larva','Larvae','Larva: Nauplius',
    'Argopecten purpuratus larvae','Ootheca/Egg mass','Nymph/Instar','Poultry',
    'poultry','poultry sources','Poultry-Associated','Chicken','chicken','Duck',
    'Skeletal system','Muscular system','Fat body','Shell','Bryozoans',
    'Cephalochordata','Nest','Sponge','Rumen','Abdomen',
    'Genome of termite-associated Trabulsiella odontotermitis strain'),
  'Aquatic (Unspecified)'=c('Aquatic','Mat, Aquatic','Aquatic, Biofilm'),
  'Industrial/Waste'=c('Industrial production','Industrial wastewater',
    'Industrial waste','Industrial/Waste','Engineered','Engineered product',
    'Genetically modified','Chemical products','Tailings pond','Bagasse',
    'Household waste','Zoo waste','Animal waste','Solid waste','Landfill',
    'Agricultural waste','Brown waste','Cellulose associated waste'),
  'Wastewater'=c('Wastewater','Sewage','WWTP','Activated Sludge',
    'Activated sludge','UASB (Upflow anaerobic sludge blanket)',
    'MBR (Membrane bioreactor)','Water treatment plant',
    'Drinking water treatment plant','Percolator','Nutrient removal',
    'Influent','Sludge','Drinking Water Filter'),
  'Built-Environment'=c('Built environment','International Space Station',
    'Spacecraft Assembly Cleanrooms','House','Building','City','Vivarium'),
  'Geological'=c('Geologic','Geological','Volcanic','Cave','halite crust','Subsurface'),
  'Saline Environments'=c('Saline','Non-marine Saline and Alkaline',
    'Saline_environments','High-salinity/high-pH','Salt'),
  'Algae'=c('Algae','Green algae','Brown Algae','Red algae','Microalgae','Diatoms'),
  'Air'=c('Air','Outdoor Air','Indoor Air'),
  'Fungi'=c('Fungi'),
  'Clinical'=c('Clinical','clinical','Hospital','Hospital-Associated',
    'Digestive system','Circulatory system','Respiratory system','Urinary system',
    'Integumentary system','Reproductive system','Nervous system',
    'Lymphatic system','Excretory system','Visual system',
    'Auditory/Hearing system','Sensory organs','Gastrointestinal tract',
    'Intestinal tract','Abdominal/Peritoneal cavity','Unspecified system',
    'Multiple systems','Head','Skin','Integument','Tissue','Whole body',
    'Remains','Acinetobacter larvae BRTC-1 Genome sequencing'),
  'Oil-Related Environments'=c('Oil reservoir','Oil refinery','Hydrocarbon'),
  'Extreme Environments'=c('Thermal springs','Hot spring','Hot Lake','Acidic',
    'Anaerobic','Aerobic','Mud volcano')
)

categorize <- function(x) {
  for (cat in names(ecosystem_categories)) {
    if (x %in% ecosystem_categories[[cat]]) return(cat)
  }
  return('EXCLUDE')  # Unknown/Uncertain excluded
}

# ── Load and prepare data ─────────────────────────────────
cat('Loading metadata...\n')
df <- read.delim(META, sep='\t', stringsAsFactors=FALSE, check.names=FALSE) %>%
  mutate(across(all_of(target_genes), ~ ifelse(as.numeric(.) > 0, 1, 0))) %>%
  mutate(Biosurf_Total = rowSums(across(all_of(target_genes)), na.rm=TRUE)) %>%
  mutate(Ecosystem_Category = sapply(Ecosystem.Type, categorize))

cat('Total genomes:', nrow(df), '\n')

# ═══════════════════════════════════════════════════════
# LAYER 1: Quality Filter
# ═══════════════════════════════════════════════════════
cat('\n[Layer 1] Quality filtering...\n')
L1 <- df %>%
  filter(CheckM2.Completeness >= 90,
         CheckM2.Contamination < 5)
cat('After quality filter:', nrow(L1), 'genomes\n')
write_csv(L1 %>% select(IMG.Genome.ID),
  file.path(OUTDIR, 'quality_filtered.txt'))

# ═══════════════════════════════════════════════════════
# LAYER 2: ANI Cluster Non-redundancy
# ═══════════════════════════════════════════════════════
cat('\n[Layer 2] ANI cluster representative selection...\n')

# For singletons: keep all
# For cliques/clique-groups: keep one with highest completeness
L2_singletons <- L1 %>%
  filter(ANI.Cluster.Type == 'singleton' | is.na(ANI.Cluster.ID))

L2_clusters <- L1 %>%
  filter(ANI.Cluster.Type != 'singleton' & !is.na(ANI.Cluster.ID)) %>%
  group_by(ANI.Cluster.ID) %>%
  slice_max(CheckM2.Completeness, n=1, with_ties=FALSE) %>%
  ungroup()

L2 <- bind_rows(L2_singletons, L2_clusters)
cat('After ANI de-replication:', nrow(L2), 'genomes\n')
write_csv(L2 %>% select(IMG.Genome.ID),
  file.path(OUTDIR, 'ani_representatives.txt'))

# ═══════════════════════════════════════════════════════
# LAYER 3: Ecosystem Stratified Sampling
# ═══════════════════════════════════════════════════════
cat('\n[Layer 3] Ecosystem stratified sampling...\n')

# Exclude Unknown/Uncertain
L2_clean <- L2 %>% filter(Ecosystem_Category != 'EXCLUDE')

# Calculate target per ecosystem using hybrid approach
# Proportional with floor=50, cap=400
eco_counts <- L2_clean %>% count(Ecosystem_Category)
total_after_L2 <- nrow(L2_clean)
TARGET_TOTAL  <- 3000  # Adjust based on Phase 1 output
FLOOR         <- 50
CAP           <- 400

eco_targets <- eco_counts %>%
  mutate(
    prop_target = round(n / total_after_L2 * TARGET_TOTAL),
    target = pmax(pmin(prop_target, CAP), pmin(FLOOR, n))
  )

cat('Ecosystem sampling targets:\n')
print(eco_targets)

# ═══════════════════════════════════════════════════════
# LAYER 4: Gene Profile Diversity Within Each Ecosystem
# ═══════════════════════════════════════════════════════
cat('\n[Layer 4] Gene profile diversity sampling...\n')

sample_with_gene_diversity <- function(eco_df, target_n) {
  n_high   <- round(target_n * 0.50)
  n_medium <- round(target_n * 0.30)
  n_low    <- target_n - n_high - n_medium

  high   <- eco_df %>% filter(Biosurf_Total >= 7)
  medium <- eco_df %>% filter(Biosurf_Total >= 3, Biosurf_Total < 7)
  low    <- eco_df %>% filter(Biosurf_Total >= 1, Biosurf_Total < 3)

  s_high   <- high   %>% slice_sample(n=min(n_high,   nrow(high)))
  s_medium <- medium %>% slice_sample(n=min(n_medium, nrow(medium)))
  s_low    <- low    %>% slice_sample(n=min(n_low,    nrow(low)))

  bind_rows(s_high, s_medium, s_low)
}

final_subset <- eco_targets %>%
  purrr::pmap_dfr(function(Ecosystem_Category, n, prop_target, target) {
    eco_df <- L2_clean %>% filter(Ecosystem_Category == !!Ecosystem_Category)
    sample_with_gene_diversity(eco_df, target)
  })

cat('\nFinal subset size:', nrow(final_subset), 'genomes\n')

# Final summary
cat('\nFinal subset by ecosystem:\n')
print(final_subset %>% count(Ecosystem_Category) %>% arrange(desc(n)))

cat('\nGene count distribution in subset:\n')
print(final_subset %>% count(Biosurf_Total) %>% arrange(Biosurf_Total))

# ── Save outputs ──────────────────────────────────────────
write_csv(final_subset, file.path(OUTDIR, 'final_subset_full_metadata.csv'))
write_lines(as.character(final_subset$IMG.Genome.ID),
            file.path(OUTDIR, 'final_subset_ids.txt'))

cat('\nPhase 2 complete.\n')
cat('Genome IDs saved to: genome_lists/final_subset_ids.txt\n')

