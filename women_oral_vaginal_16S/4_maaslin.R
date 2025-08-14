####LOAD####
library(phyloseq)
library(reshape2)
library(dplyr)
library(Maaslin2)
library(ggplot2)

# Set working directory
# Load phyloseq object

#Uncomment the stratification of interest to run maasling for all taxa levels
load("/Users/danielavargasrobles/Dropbox (UMass Medical School)/Ana_projects/2_PREECLAMPSIA/P16S/R/phy_objects/PHYrar_alpha_ggext_subset02_newseqs.RData")
phy <- PHYrar_alpha_ggext_subset02_newseqs
# phy <- subset_samples(phy,   ethnicity == "Hispanic/Latina")
phy <- subset_samples(
  phy, race == "White/Caucasian" & ethnicity == "Not Hispanic/Latina")


# Filter samples
phy <- subset_samples(phy, body_site == "MomTongue")
phy <- subset_samples(phy, study_id != "401")

# Define taxonomic levels to analyze
tax_levels <- c("Species", "Genus", "Family", "Class", "Order", "Phylum")

# Date for output folder
today <- Sys.Date()

# Process ASV level (no taxonomic agglomeration)
cat("Processing ASV\n")
phy_level <- phy
phy_prop <- transform_sample_counts(phy_level, function(x) x / sum(x) * 100)
perc_samples <- 0.20
phy_fil <- phyloseq::filter_taxa(phy_prop, function(x) sum(x > 0) > (perc_samples * length(x)), TRUE)
abundance_threshold <- 1
phy_fil <- prune_taxa(taxa_sums(phy_fil) > abundance_threshold, phy_fil)
taxa_names(phy_fil) <- make.unique(taxa_names(phy_fil))

# Extract OTU table, metadata, and taxonomy
otu_table <- as.data.frame(t(otu_table(phy_fil)))
sample_metadata <- as(sample_data(phy_fil), "data.frame")

# Define output directory
output_dir <- paste0("results/maaslin/", today, "/ASV_study_group_20prev_sig0_25_ab1per")
dir.create(output_dir, recursive = TRUE)

# Run MaAsLin2
fit_data <- Maaslin2(
  input_data = otu_table, 
  input_metadata = sample_metadata, 
  output = output_dir,
  fixed_effects = c("study_group", "bmi",
                    "antibiotics_during_pregnancy", "ethnicity"),
  random_effects = NULL,
  normalization = "TSS",
  max_significance = 0.25,
  transform = "LOG", 
  analysis_method = "LM", 
  correction = "BH", 
  standardize = TRUE,
  min_abundance = 0,
  min_prevalence = 0
)

for (tax_level in tax_levels) {
  cat("Processing", tax_level, "\n")
  
  # Aggregate data if necessary
  phy_level <- tax_glom(phy, taxrank = tax_level)
  
  # Transform to relative abundance
  phy_prop <- transform_sample_counts(phy_level, function(x) x / sum(x) * 100)
  
  # Filter by prevalence and abundance
  phy_fil <- phyloseq::filter_taxa(phy_prop, function(x) sum(x > 0) > (perc_samples * length(x)), TRUE)
  phy_fil <- prune_taxa(taxa_sums(phy_fil) > abundance_threshold, phy_fil)
  
  # Ensure unique taxa names
  taxa_names(phy_fil) <- make.unique(tax_table(phy_fil)[, tax_level])
  
  # Extract OTU table, metadata, and taxonomy
  otu_table <- as.data.frame(t(otu_table(phy_fil)))
  sample_metadata <- as(sample_data(phy_fil), "data.frame")
  
  # Define output directory
  output_dir <- paste0("results/maaslin/", today, "/", tax_level, "_study_group_20prev_sig0_25_ab1per")
  dir.create(output_dir, recursive = TRUE)
  
  # Run MaAsLin2
  fit_data <- Maaslin2(
    input_data = otu_table, 
    input_metadata = sample_metadata, 
    output = output_dir,
    fixed_effects = c("study_group", "bmi",
                      "antibiotics_during_pregnancy", "ethnicity"),
  random_effects = NULL,
  normalization = "TSS",
  max_significance = 0.25,
  transform = "LOG", 
  analysis_method = "LM", 
  correction = "BH", 
  standardize = TRUE,
  min_abundance = 0,
  min_prevalence = 0
  )

# Load results
results_file <- paste0(output_dir, "/all_results.tsv")
if (file.exists(results_file)) {
  maaslin_results <- read.table(results_file, header = TRUE, sep = "\t")
  
  # Filter significant results
  significant_results <- subset(maaslin_results, pval < 0.05 & metadata == 'study_group' & value == 'PE')
  significant_results <- significant_results[order(significant_results$coef, decreasing = TRUE), ]
  significant_results$feature <- factor(significant_results$feature, levels = significant_results$feature)
  
  # Create coefficient plot
  p <- ggplot(significant_results, aes(x = feature, y = coef)) +
    geom_point() +
    geom_errorbar(aes(ymin = coef - stderr, ymax = coef + stderr), width = 0.2, color = "red") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "blue") +
    labs(title = paste("Significant Coeff (", tax_level, " Study Group: PE, pval < 0.25)", sep = ""),
         x = "Feature", y = "Coefficient") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  # Save plot
  output_plot <- paste0(output_dir, "/coef_plot_PE_study_group_significant_ordered.pdf")
  ggsave(output_plot, plot = p, width = 6, height = 4)
}
}

cat("Analysis completed!\n")


####identifying particular ASV####

phy <- PHYrar_alpha_ggext_subset02_newseqs

# Asegurar que el ASV que buscas está en el objeto phyloseq
asv_id <- "d891a80eccfa02001332b42560424f90"

# Extraer la tabla de taxonomía
tax_table_df <- as.data.frame(tax_table(phy))

# Verificar si el ASV está presente en la tabla
if (asv_id %in% rownames(tax_table_df)) {
  # Extraer la información taxonómica del ASV
  species_info <- tax_table_df[asv_id, ]
  
  # Mostrar el resultado
  print(species_info)
} else {
  print("El ASV especificado no se encuentra en el objeto phyloseq.")
}


