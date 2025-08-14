# Generation of 16S infant stool phyloseq object from QIIME2 output and metadata
# This script processes the 16S rRNA sequencing data from infant stool samples
# and prepares it for further analysis in R using the phyloseq package.


# Load necessary libraries
library(phyloseq)
library(ggplot2)
library(RColorBrewer)
library(vegan)
library(OTUtable)
library(reshape2)
library(ape)
library(qiime2R)
library(tidyverse)
library(readr)
library(tidyr)
library(PERFect)
library(devtools)
library(remotes)
library(QsRutils)
library(readxl)
library(tibble)

####LOAD####
setwd("~/Dropbox (UMass Medical School)/Ana_projects/2_PREECLAMPSIA/P16S/meconium/")

phy <- qza_to_phyloseq(features = "qiime2_files/table_demux.qza")
ASV <- otu_table(phy)

data<-as.data.frame(t(ASV))

# load taxonomies
taxonomy <- read_tsv('qiime2_files/taxonomy/taxonomy.tsv', show_col_types = FALSE)

# clean taxonomy
taxonomy_cleaned <- taxonomy %>%
  dplyr::select(-Confidence) %>%  
  separate(Taxon, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
           sep = ";", fill = "right") %>%
  mutate(across(everything(), ~ gsub("^\\s+", "", .))) %>%  
  mutate(across(everything(), ~ gsub(" ", "_", .))) %>%
  mutate(
    Kingdom = ifelse(is.na(Kingdom) | Kingdom == "" | Kingdom == "d__", "d__Unclassified", Kingdom),
    Phylum = ifelse(is.na(Phylum) | Phylum == "" | Phylum == "p__", "p__Unclassified", Phylum),
    Class = ifelse(is.na(Class) | Class == "" | Class == "c__", "c__Unclassified", Class),
    Order = ifelse(is.na(Order) | Order == "" | Order == "o__", "o__Unclassified", Order),
    Family = ifelse(is.na(Family) | Family == "" | Family == "f__", "f__Unclassified", Family),
    Genus = ifelse(is.na(Genus) | Genus == "" | Genus == "g__", "g__Unclassified", Genus),
    Species = ifelse(
      is.na(Species) | Species == "" | Species == "s__" | Species == "s__Unclassified",
      ifelse(Genus == "g__Unclassified", "g__Unclassified_s__Unclassified", paste0(Genus, "_s__Unclassified")),
      Species  
    )
  ) %>%
  column_to_rownames(var = "Feature ID") 

# Filter out unwanted taxa
taxonomy_filtered <- taxonomy_cleaned %>%
  filter(Kingdom != "d__Archaea") %>%
  filter(Kingdom != "Unassigned") %>%
  filter(Phylum != "p__Unclassified") %>%
  filter(Class != "c__Unclassified") %>%
  filter(Order != "o__Unclassified") %>%
  filter(!if_any(everything(), ~ grepl("mitochondria", ., ignore.case = TRUE)))

#write.csv(taxonomy_filtered, "tables_clean/tax_cleaned_meconium.csv")

# Filter genera with multiple classifications

genera_multiple_classifications <- taxonomy_filtered %>%
  filter(Genus != "g__Unclassified") %>%  # Excluir 'g__Unclassified'
  group_by(Genus) %>%
  summarise(
    unique_kingdoms = n_distinct(Kingdom),
    unique_phylums = n_distinct(Phylum),
    unique_classes = n_distinct(Class),
    unique_orders = n_distinct(Order),
    unique_families = n_distinct(Family)
  ) %>%
  filter(unique_kingdoms > 1 | unique_phylums > 1 | unique_classes > 1 |
           unique_orders > 1 | unique_families > 1)

# Extact the genera with multiple classifications
genera_to_verify <- unique(genera_multiple_classifications$Genus)

taxonomy_cleaned <- tax_table(as.matrix(taxonomy_filtered))

# PRUNE taxonomy after ASV table filtering
common_taxa <- intersect(taxa_names(taxonomy_cleaned), taxa_names(ASV))
TAXA_clean <- prune_taxa(common_taxa, taxonomy_cleaned)
OTU_clean <- prune_taxa(common_taxa, ASV)

# Create the phyloseq object
PHYraw <- phyloseq(TAXA_clean, OTU_clean)
PHYraw
# Save the phyloseq object
save(PHYraw, file = "phy_obj/PHYraw.RData")

#### Filtering ab and prev ####

# We apply the same filters as for maternal fecal WGS

# Transform the phyloseq object to relative abundances
phy1_rel <- transform_sample_counts(PHYrar, function(x) x / sum(x))


# Filter taxa with relative abundance > 0.01
phy1_rel_filtered <- prune_taxa(taxa_sums(phy1_rel) > 0.01, phy1_rel)

# Filter taxa with prevalence > 1
taxa_presence <- apply(X = otu_table(phy1_rel_filtered), MARGIN = 1, FUN = function(x) sum(x > 0))
phy1_rel_filtered_prev <- prune_taxa(taxa_presence > 1, phy1_rel_filtered)

# Since we want to return to counts (not relative abundances), we apply the filters to the original phyloseq object:
# we keep only the taxa that survived both filters
phy1_filtered_final <- prune_taxa(taxa_names(phy1_rel_filtered_prev), phy1)

phy1_filtered_final

#extrat otu table to see the colnames
OTU_table <- as.data.frame(otu_table(phy1_filtered_final))
OTU_table

####METADATA####
# load
met <- read.csv("~/Dropbox (UMass Medical School)/Ana_projects/2_PREECLAMPSIA/met/met_clean0.csv")


# create formatted study_id including "B" prefix
met$study_id_formatted <- paste0("B", met$study_id)


# extract tax_table and otu_table from the filtered phyloseq object
TAXA_clean <- as.matrix(tax_table(phy1_filtered_final))
OTU_clean <- as.matrix(otu_table(phy1_filtered_final, taxa_are_rows = TRUE))

# Sort the metadata based on the formatted study_id
# Ensure that the metadata is ordered according to the columns of OTU_clean
met_ordered <- met[match(colnames(OTU_clean), met$study_id_formatted), ]

# assign the rownames of the metadata to match the column names of OTU_clean
rownames(met_ordered) <- colnames(OTU_clean)

# Create the phyloseq object with the filtered OTU table, taxonomy table, and ordered metadata
PHY_fil <- phyloseq(
  otu_table(OTU_clean, taxa_are_rows = TRUE),
  tax_table(TAXA_clean),
  sample_data(met_ordered)
)

#extract metadata
met_phylo <- data.frame(sample_data(PHY_fil))
met_phylo
save(PHY_fil,file="phy_obj/PHY_fil.RData")
#load("R/phy_objects/PHYraw.RData")

#### Rarefaction ####

# Calculate the sum of sequences per sample
sumatoria <- sample_sums(PHY_fil)
min(sumatoria)

sumatoria_ordenada <- sort(sumatoria)
print(head(sumatoria_ordenada, 30))

# Add the sum of sequences as a new column in the metadata
sample_data(PHY_fil)$num_seqs <- sumatoria

# Show the first few rows of the metadata with the new column
print(sample_data(PHY_fil))

seq_to_check=data.frame(sample_data(PHY_fil))

seq_to_check %>% 
  dplyr::select(num_seqs) %>% 
  arrange(-num_seqs)


# rarefaction
PHY_fil_rar=rarefy_even_depth(PHY_fil, sample.size = 43539,
                         rngseed = 711, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

#One note that can be useful, I also use the estimation of Good's indexes (library(QsRutils)) to evaluate, for instance, if the rarefaction approach is adequate, by looking at the median and mean values
#you can see that both mean and median values dropped a bit after rarefaction. If they remain >90% its ok otherwise the rarefied samples are not reflecting adequately the original taxa diversity
not_rare=summary(goods(otu_table(PHY_fil_rar)))
not_rare
rare=summary(goods(otu_table(PHY_fil_rar)))
rare

#save rarefied object
save(PHY_fil_rar,file="phy_obj/PHY_fil_rar.RData")
met=sample_data(PHY_fil_rar)

#### ALPHA diversity indices#####
#1. extract ASV table 
OTU1 = as(otu_table(PHY_fil_rar), "matrix")
OTUdf = as.data.frame(OTU1)

#Calculating 4 different alpha diversity metrics
#at asv level
shannon = diversity(OTUdf, index = "shannon", MARGIN = 2, base = exp(1))
simpson = diversity(OTUdf, index = "simpson", MARGIN = 2)
#chao = apply(OTUdf, 2, chao1)
chao=estimate_richness(PHY_fil_rar, split = TRUE,measures = "Chao1")
obs = apply(OTUdf, 2, function(x){
  nws = length(which(x>0))
  return(nws)
})

div_asv = data.frame(shannon = shannon, 
                     simpson = simpson,
                     chao1 = chao,
                     observed = obs)

#merging columns of alpha metrics with the rest of the metadata by SampleID
md1 = sample_data(PHY_fil_rar)
md0=as.data.frame(as(md1, "matrix"))
write.csv(md0, "md0.temp", row.names = TRUE)
md = read.csv("md0.temp")

asv=cbind(md,div_asv)
rownames(asv) = sample_names(PHY_fil_rar)

sampledata = sample_data(asv)

#Re doing phyloseq object with new added columns (alpha div metrics) in metadata
PHY_fil_rar_alpha_meco = phyloseq(otu_table(PHY_fil_rar),tax_table(PHY_fil_rar), sampledata)#phy_tree(PHY_fil_rar) )
PHY_fil_rar_alpha_meco

#### Save ####
save(PHY_fil_rar_alpha_meco,file="phy_obj/PHY_fil_rar_alpha_meco.RData")
