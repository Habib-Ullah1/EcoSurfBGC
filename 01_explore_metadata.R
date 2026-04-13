#!/usr/bin/env Rscript
# 01_explore_metadata.R
# EcoSurfBGC Project — Phase 1: Metadata Exploration
# Author: Habib Ullah
# Version: 1.0 | March 2026

library(tidyverse)

# ── Paths ────────────────────────────────────────────────
BASE   <- '/data/habib/EcoSurfBGC'
META   <- file.path(BASE, 'metadata',
          'Merged_raw_CheckM2_Filtered_JGI-info-table_newTax_ProfileAdded_acc_removed.txt')
OUTDIR <- file.path(BASE, 'results')

# ── Target genes ─────────────────────────────────────────
target_genes <- c('fen_pps','fab','itu','rhl','rml','wza','wzb','wzc',
                  'wzy','lic_srf','srf','adh','alg','arf','pso','est','psw','cyp')

# ── Load data ─────────────────────────────────────────────
cat('Loading metadata...\n')
df <- read.delim(META, sep='\t', stringsAsFactors=FALSE, check.names=FALSE) %>%
  mutate(across(all_of(target_genes), ~ ifelse(as.numeric(.) > 0, 1, 0))) %>%
  mutate(Biosurf_Total = rowSums(across(all_of(target_genes)), na.rm=TRUE))

cat('Total genomes loaded:', nrow(df), '\n')

# ── 1. Completeness distribution ─────────────────────────
cat('\n--- CheckM2 Completeness Distribution ---\n')
print(summary(df$CheckM2.Completeness))
cat('Genomes >= 90% complete:', sum(df$CheckM2.Completeness >= 90, na.rm=TRUE), '\n')
cat('Genomes >= 80% complete:', sum(df$CheckM2.Completeness >= 80, na.rm=TRUE), '\n')
cat('Genomes >= 70% complete:', sum(df$CheckM2.Completeness >= 70, na.rm=TRUE), '\n')

# ── 2. Contamination distribution ────────────────────────
cat('\n--- CheckM2 Contamination Distribution ---\n')
print(summary(df$CheckM2.Contamination))
cat('Genomes < 5% contamination:', sum(df$CheckM2.Contamination < 5,  na.rm=TRUE), '\n')
cat('Genomes < 10% contamination:', sum(df$CheckM2.Contamination < 10, na.rm=TRUE), '\n')

# ── 3. Combined quality filter ───────────────────────────
cat('\n--- Combined Quality Filter (>=90% complete, <5% contamination) ---\n')
q_filtered <- df %>% filter(CheckM2.Completeness >= 90, CheckM2.Contamination < 5)
cat('Genomes passing quality filter:', nrow(q_filtered), '\n')

# ── 4. ANI cluster distribution ──────────────────────────
cat('\n--- ANI Cluster Distribution (quality-filtered genomes) ---\n')
print(table(q_filtered$ANI.Cluster.Type, useNA='always'))
cat('Unique ANI cluster IDs:', length(unique(q_filtered$ANI.Cluster.ID)), '\n')
cat('Singletons:', sum(q_filtered$ANI.Cluster.Type == 'singleton', na.rm=TRUE), '\n')
cat('Clique members:', sum(q_filtered$ANI.Cluster.Type %in% c('clique','clique-group'), na.rm=TRUE), '\n')

# ── 5. Ecosystem distribution ────────────────────────────
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
  return('Unclassified/Uncertain')
}

q_filtered <- q_filtered %>%
  mutate(Ecosystem_Category = sapply(Ecosystem.Type, categorize))

eco_summary <- q_filtered %>%
  group_by(Ecosystem_Category) %>%
  summarise(
    n_genomes  = n(),
    mean_genes = round(mean(Biosurf_Total), 2),
    .groups    = 'drop'
  ) %>%
  arrange(desc(n_genomes))

cat('\n--- Ecosystem Distribution (quality-filtered genomes) ---\n')
print(eco_summary, n=25)

# ── 6. Biosurfactant gene count distribution ─────────────
cat('\n--- Biosurfactant Gene Count Distribution (quality-filtered) ---\n')
gene_dist <- q_filtered %>%
  count(Biosurf_Total) %>%
  mutate(pct = round(n / sum(n) * 100, 2))
print(gene_dist, n=30)
cat('Genomes with >= 7 genes:', sum(q_filtered$Biosurf_Total >= 7), '\n')
cat('Genomes with 3-6 genes: ', sum(q_filtered$Biosurf_Total >= 3 & q_filtered$Biosurf_Total < 7), '\n')
cat('Genomes with 1-2 genes: ', sum(q_filtered$Biosurf_Total >= 1 & q_filtered$Biosurf_Total < 3), '\n')
cat('Genomes with 0 genes:   ', sum(q_filtered$Biosurf_Total == 0), '\n')

# ── Save outputs ──────────────────────────────────────────
write_csv(eco_summary,   file.path(OUTDIR, 'ecosystem_summary.csv'))
write_csv(gene_dist,     file.path(OUTDIR, 'gene_count_distribution.csv'))
write_csv(q_filtered %>%
  select(IMG.Genome.ID, CheckM2.Completeness, CheckM2.Contamination,
         ANI.Cluster.ID, ANI.Cluster.Type, Ecosystem_Category, Biosurf_Total),
  file.path(OUTDIR, 'metadata_summary_qfiltered.csv'))

cat('\nPhase 1 complete. Results saved to:', OUTDIR, '\n')
