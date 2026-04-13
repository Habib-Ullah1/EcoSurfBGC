#!/usr/bin/env Rscript
# 01_explore_metadata_base.R
# EcoSurfBGC — Phase 1: Metadata Exploration (base R, no dependencies)
# Author: Habib Ullah | Version: 1.0 | March 2026

BASE   <- '/data/habib/EcoSurfBGC'
META   <- file.path(BASE, 'metadata',
          'Merged_raw_CheckM2_Filtered_JGI-info-table_newTax_ProfileAdded_acc_removed.txt')
OUTDIR <- file.path(BASE, 'results')

target_genes <- c('fen_pps','fab','itu','rhl','rml','wza','wzb','wzc',
                  'wzy','lic_srf','srf','adh','alg','arf','pso','est','psw','cyp')

cat('Loading metadata...\n')
df <- read.delim(META, sep='\t', stringsAsFactors=FALSE, check.names=FALSE)

for (g in target_genes) {
  if (g %in% colnames(df)) {
    df[[g]] <- ifelse(as.numeric(df[[g]]) > 0, 1, 0)
  }
}
df$Biosurf_Total <- rowSums(df[, target_genes], na.rm=TRUE)
cat('Total genomes loaded:', nrow(df), '\n')

cat('\n--- CheckM2 Completeness Distribution ---\n')
print(summary(df$CheckM2.Completeness))
cat('Genomes >= 90% complete:', sum(df$CheckM2.Completeness >= 90, na.rm=TRUE), '\n')
cat('Genomes >= 80% complete:', sum(df$CheckM2.Completeness >= 80, na.rm=TRUE), '\n')
cat('Genomes >= 70% complete:', sum(df$CheckM2.Completeness >= 70, na.rm=TRUE), '\n')

cat('\n--- CheckM2 Contamination Distribution ---\n')
print(summary(df$CheckM2.Contamination))
cat('Genomes < 5% contamination:',  sum(df$CheckM2.Contamination < 5,  na.rm=TRUE), '\n')
cat('Genomes < 10% contamination:', sum(df$CheckM2.Contamination < 10, na.rm=TRUE), '\n')

cat('\n--- Combined Quality Filter (>=90% complete, <5% contamination) ---\n')
q_filtered <- df[!is.na(df$CheckM2.Completeness) & !is.na(df$CheckM2.Contamination) &
                 df$CheckM2.Completeness >= 90 & df$CheckM2.Contamination < 5, ]
cat('Genomes passing quality filter:', nrow(q_filtered), '\n')

cat('\n--- ANI Cluster Distribution (quality-filtered) ---\n')
print(table(q_filtered$ANI.Cluster.Type, useNA='always'))
cat('Unique ANI cluster IDs:', length(unique(q_filtered$ANI.Cluster.ID)), '\n')
cat('Singletons:',     sum(q_filtered$ANI.Cluster.Type == 'singleton',                 na.rm=TRUE), '\n')
cat('Clique members:', sum(q_filtered$ANI.Cluster.Type %in% c('clique','clique-group'), na.rm=TRUE), '\n')
cat('NA cluster type:',sum(is.na(q_filtered$ANI.Cluster.Type)), '\n')

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
    if (!is.na(x) && x %in% ecosystem_categories[[cat]]) return(cat)
  }
  return('Unclassified/Uncertain')
}

cat('Categorizing ecosystems (this may take a minute)...\n')
q_filtered$Ecosystem_Category <- sapply(q_filtered$Ecosystem.Type, categorize)

eco_tab   <- sort(table(q_filtered$Ecosystem_Category), decreasing=TRUE)
eco_means <- tapply(q_filtered$Biosurf_Total, q_filtered$Ecosystem_Category, mean)
eco_summary <- data.frame(
  Ecosystem_Category = names(eco_tab),
  n_genomes          = as.integer(eco_tab),
  mean_genes         = round(eco_means[names(eco_tab)], 2),
  row.names          = NULL
)

cat('\n--- Ecosystem Distribution (quality-filtered genomes) ---\n')
print(eco_summary, row.names=FALSE)

cat('\n--- Biosurfactant Gene Count Distribution ---\n')
gene_tab <- table(q_filtered$Biosurf_Total)
gene_df  <- data.frame(
  Biosurf_Total = as.integer(names(gene_tab)),
  n             = as.integer(gene_tab),
  pct           = round(as.integer(gene_tab) / nrow(q_filtered) * 100, 2),
  row.names     = NULL
)
print(gene_df, row.names=FALSE)
cat('Genomes with >= 7 genes:', sum(q_filtered$Biosurf_Total >= 7), '\n')
cat('Genomes with 3-6 genes: ', sum(q_filtered$Biosurf_Total >= 3 & q_filtered$Biosurf_Total < 7), '\n')
cat('Genomes with 1-2 genes: ', sum(q_filtered$Biosurf_Total >= 1 & q_filtered$Biosurf_Total < 3), '\n')
cat('Genomes with 0 genes:   ', sum(q_filtered$Biosurf_Total == 0), '\n')

write.csv(eco_summary, file.path(OUTDIR, 'ecosystem_summary.csv'),       row.names=FALSE)
write.csv(gene_df,     file.path(OUTDIR, 'gene_count_distribution.csv'), row.names=FALSE)
write.csv(
  q_filtered[, c('IMG.Genome.ID','CheckM2.Completeness','CheckM2.Contamination',
                 'ANI.Cluster.ID','ANI.Cluster.Type','Ecosystem_Category','Biosurf_Total')],
  file.path(OUTDIR, 'metadata_summary_qfiltered.csv'), row.names=FALSE)

cat('\nPhase 1 complete. Results saved to:', OUTDIR, '\n')
