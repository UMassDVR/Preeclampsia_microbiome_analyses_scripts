# Generating the phyloseq object
#write all comments of this script in english

#load the necessary libraries
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

####LOAD####
setwd("~/Dropbox (UMass Medical School)/Ana_projects/2_PREECLAMPSIA/")

phy <- qza_to_phyloseq(features = "P16S/qiime2_qza_vag_tong_ggext_newseqs/table_vag_tong_subsample02.qza")

ASV <- otu_table(phy)

write.csv(ASV, "P16S/qiime2_qza_vag_tong_ggext_newseqs/ASV_table_raw.csv")

data<-as.data.frame(t(ASV))
dim(data)

####Clean Taxonomy ####
taxonomy <- read_tsv('P16S/qiime2_qza_vag_tong_ggext_newseqs/taxonomy_full_length_subsample02_gg_ext/taxonomy.tsv', show_col_types = FALSE)

taxonomy_cleaned0 <- taxonomy %>%
  dplyr::select(-Confidence) %>% 
  separate(Taxon, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";", fill = "right") %>%
  mutate(across(everything(), ~ gsub("^\\s+", "", .))) %>%
  mutate(across(everything(), ~ gsub(" ", "_", .))) %>%  
  mutate(
    Kingdom = ifelse(is.na(Kingdom) | Kingdom == '' | Kingdom == 'k__', 'k__Unclassified', Kingdom),
    Phylum = ifelse(is.na(Phylum) | Phylum == '' | Phylum == 'p__', 'p__Unclassified', Phylum),
    Class = ifelse(is.na(Class) | Class == '' | Class == 'c__', 'c__Unclassified', Class),
    Order = ifelse(is.na(Order) | Order == '' | Order == 'o__', 'o__Unclassified', Order),
    Family = ifelse(is.na(Family) | Family == '' | Family == 'f__', 'f__Unclassified', Family),
    Genus = ifelse(is.na(Genus) | Genus == '' | Genus == 'g__', 'g__Unclassified', Genus),
    Species = ifelse(Genus == 'g__Finegoldia', 'g__Finegoldia_magna', Species),  
    Species = ifelse(is.na(Species) | Species == '' | Species == 's__', 's__Unclassified', Species),
    Species = ifelse(Genus != 'g__Unclassified' & Species == 's__Unclassified', paste0(Genus, "_s__Unclassified"), 
                     ifelse(Genus != 'g__Unclassified' & Species != 's__Unclassified', paste0(Genus, " ", sub('s__', '', Species)), 
                            "g__Unclassified_s__Unclassified"))
  ) %>%
  column_to_rownames(var = "Feature ID") %>% 
  filter(Kingdom != "k__Archaea" & Kingdom != "Unassigned" & Phylum != "p__Unclassified", 
         Class != "c__Unclassified", Order != "o__Unclassified", Family!= "f__mitochondria")

write.csv(taxonomy_cleaned0, "taxonomy_cleaned0.csv")

# Generate a phyloseq object with the cleaned taxonomy woth modern classifications
modern_classifications <- tibble::tribble(
  ~Genus, ~Kingdom, ~Phylum, ~Class, ~Order, ~Family,
  "g__Anaerococcus", "k__Bacteria", "p__Firmicutes", "c__Clostridia", "o__Clostridiales", "f__Peptoniphilaceae",
  "g__Atopobium", "k__Bacteria", "p__Actinobacteria", "c__Actinobacteria", "o__Coriobacteriales", "f__Coriobacteriaceae",
  "g__Catonella", "k__Bacteria", "p__Firmicutes", "c__Clostridia", "o__Clostridiales", "f__Lachnospiraceae",
  "g__Corynebacterium", "k__Bacteria", "p__Actinobacteria", "c__Actinobacteria", "o__Corynebacteriales", "f__Corynebacteriaceae",
  "g__Dialister", "k__Bacteria", "p__Firmicutes", "c__Negativicutes", "o__Selenomonadales", "f__Veillonellaceae",
  "g__Finegoldia", "k__Bacteria", "p__Firmicutes", "c__Clostridia", "o__Clostridiales", "f__Peptoniphilaceae",
  "g__Fusobacterium", "k__Bacteria", "p__Fusobacteria", "c__Fusobacteriia", "o__Fusobacteriales", "f__Fusobacteriaceae",
  "g__Gemella", "k__Bacteria", "p__Firmicutes", "c__Bacilli", "o__Gemellales", "f__Gemellaceae",
  "g__Helcococcus", "k__Bacteria", "p__Firmicutes", "c__Clostridia", "o__Clostridiales", "f__Peptoniphilaceae",
  "g__Megasphaera", "k__Bacteria", "p__Firmicutes", "c__Negativicutes", "o__Selenomonadales", "f__Veillonellaceae",
  "g__Mogibacterium", "k__Bacteria", "p__Firmicutes", "c__Clostridia", "o__Clostridiales", "f__Peptostreptococcaceae",
  "g__Mycoplasma", "k__Bacteria", "p__Tenericutes", "c__Mollicutes", "o__Mycoplasmatales", "f__Mycoplasmataceae",
  "g__Parvimonas", "k__Bacteria", "p__Firmicutes", "c__Clostridia", "o__Clostridiales", "f__Peptoniphilaceae",
  "g__Peptoniphilus", "k__Bacteria", "p__Firmicutes", "c__Clostridia", "o__Clostridiales", "f__Peptoniphilaceae",
  "g__Peptostreptococcus", "k__Bacteria", "p__Firmicutes", "c__Clostridia", "o__Clostridiales", "f__Peptostreptococcaceae",
  "g__Propionibacterium", "k__Bacteria", "p__Actinobacteria", "c__Actinobacteria", "o__Propionibacteriales", "f__Propionibacteriaceae",
  "g__Selenomonas", "k__Bacteria", "p__Firmicutes", "c__Negativicutes", "o__Selenomonadales", "f__Veillonellaceae",
  "g__Treponema", "k__Bacteria", "p__Spirochaetes", "c__Spirochaetia", "o__Spirochaetales", "f__Spirochaetaceae",
  "g__Veillonella", "k__Bacteria", "p__Firmicutes", "c__Negativicutes", "o__Selenomonadales", "f__Veillonellaceae",
  "g__Streptococcus", "k__Bacteria", "p__Firmicutes", "c__Bacilli", "o__Lactobacillales", "f__Streptococcaceae",  # Caso 1
  "g__Ureaplasma", "k__Bacteria", "p__Tenericutes", "c__Mollicutes", "o__Mycoplasmatales", "f__Mycoplasmataceae"  # Caso 2
)

taxonomy_cleaned0 <- taxonomy_cleaned0 %>%
  rownames_to_column(var = "ASV_ID")

# merge the modern classifications with the cleaned taxonomy
taxonomy_cleaned1 <- taxonomy_cleaned0 %>%
  left_join(modern_classifications, by = "Genus", suffix = c("_old", "_modern")) %>%
  mutate(
    Kingdom = if_else(Genus != "g__Unclassified" & !is.na(Kingdom_modern), Kingdom_modern, Kingdom_old),
    Phylum = if_else(Genus != "g__Unclassified" & !is.na(Phylum_modern), Phylum_modern, Phylum_old),
    Class = if_else(Genus != "g__Unclassified" & !is.na(Class_modern), Class_modern, Class_old),
    Order = if_else(Genus != "g__Unclassified" & !is.na(Order_modern), Order_modern, Order_old),
    Family = if_else(Genus != "g__Unclassified" & !is.na(Family_modern), Family_modern, Family_old)
  ) %>%
  dplyr::select(-ends_with("_old"), -ends_with("_modern"))

#recover the original rownames (ASV_ID)
taxonomy_cleaned1 <- taxonomy_cleaned1 %>%
  column_to_rownames(var = "ASV_ID")

# We use apply to check each row if the pattern "mitochondria" is present in any column
# taxonomy_cleaned1[apply(taxonomy_cleaned1, 1, function(row) any(grepl("chloroplast", row))), ]

##### Manual Taxonomy ####
# Load the manual taxonomy table
# This table contains the ASV IDs, Genus, and Species for those ASVs that were manually classified
# This table was created based on the ASVs with more than 5,000 sequences that were classified as "Unclassified" at the species level.
# The sequences were extracted from the rep_seq file and manually blasted in NCBI
# The blast results were selected based on 100% identity and the highest total score. If no species classification was available, "Unclassified" was retained.
# This improved the bar plots significantly.

manual_taxonomy <- data.frame(
  ASV_ID = c(
    "5f32e0d308736382d1e4a2cc3a800f03", "3314486729173c20661f139588a86762", 
    "f8d97633268d0cb33a4961ac1d1991e3", "c577d055d016e6c2dfaa9ca155f9f51c", 
    "5a697167dd760a129253902fa6390f95", "7d062509323aa7529dacc301efa91a1e", 
    "eed45c46ab8412632a530e3b2615671e", "f1a1b898667f658b6ba0abe243575021", 
    "61b219d4610c015b5faa10bb8e305700", "821045e2195d26f66a1fe7605b57a6c1", 
    "30f38eea16b97b9f820b34ad88223408", "132c6edec789232615eb69bf39db8b2a", 
    "a019d0954fc18f6c7cff143043b12a87", "176d59c6cca071c80c2b60e9393aa6fc", 
    "2a8e9f5d1713e964b80a0806830fe42b", "969f4202b63057ecc32b1233c3ae2f7c", 
    "4968584ae5a1560f360301ebdffd04f1", "2c5b4b7e04c1cc6fdec353a011a50795", 
    "2bc2f4ce95b919a17db546dd04aa5b09", "e1bf1a1bd4c69c93972eb2210d26ca09", 
    "70c88fdd526ea3b46145647fcddd51dd", "8b5799bb6d57eb3c327a2b21e3a0c3fe", 
    "86e514302b29b2a3bf7b6e62827fdf67", "7da5699cd3ff61097747c640b596cb85", 
    "2569ad3a1d65c0b7e46b029fedc81740",
    "a6f1a209e0780747e50ccf8117c0c854", "04e96f8349e34826431baa6190c35683",
    "06790a632d34be8db2eb09712082abd0", "0cab42d78e6dbabb304a0d5dfe88f356",
    "5158b864829c2ec1c77c13e2e993885a", "163d0544a78051aa4fa36b65d02c17b6",
    "a48572566c42efef9aa4d0c6641561eb", "514b75f20ccc23afa486f434df15f3fe",
    "08968b4577cd35855818f4fe99ac1859", "521d07ab040f6421c93b514d6750fb65",
    "58fd461434d936e5aab4c275d6f6cde7", "2bbcd6a04228ca5ccf105b904828c529",
    "843344d09a7c84470c990f1b6170032b", "7faa1e175e9fb5b005644773c97346ac",
    "930e445523e9c5d5aad858b2e20ba1cf", "1452b3e007e3bab260efe9108618bdf6",
    "9260ffbad675235f1f37dff12ba86038", "4c31e2d1733ca82e1cf384498344ec3b",
    "2239edacecc6a136f597624c94c63086"
  ),
  Genus = c(
    "g__Streptococcus", "g__Streptococcus", "g__Schaalia", "g__Streptococcus", 
    "g__Lactobacillus", "g__Actinomyces", "g__Streptococcus", "g__Rothia", 
    "g__Streptococcus", "g__Neisseria", "g__Veillonella", "g__Schaalia", 
    "g__Neisseria", "g__Granulicatella", "g__Streptococcus", "g__Neisseria", 
    "g__Streptococcus", "g__Streptococcus", "g__Ureaplasma", "g__Streptococcus", 
    "g__Veillonella", "g__Ureaplasma", "g__Streptococcus", "g__Neisseria", 
    "g__Streptococcus",
    "g__Gemella", "g__Neisseria", "g__Streptococcus", "g__Streptococcus", 
    "g__Streptococcus", "g__Lautropia", "g__Streptococcus", "g__Actinomyces", 
    "g__Streptococcus", "g__Streptococcus", "g__Streptococcus", "g__Streptococcus", 
    "g__Finegoldia", "g__Streptococcus", "g__Streptococcus", "g__Streptococcus", 
    "g__Unclassified", "g__Streptococcus", "g__Streptococcus"
  ),
  Species = c(
    "mitis", "Unclassified", "odontolytica", "gwangjuense", 
    "johnsonii", "Unclassified", "Unclassified", "dentocariosa", 
    "Unclassified", "mucosa", "Unclassified", "odontolytica", 
    "subflava", "adiacens", "pseudopneumoniae", "RH3002v2g", 
    "mitis", "sp. oral taxon 431", "parvum", "mitis", 
    "Unclassified", "parvum", "mitis", "Unclassified", 
    "Unclassified",
    "Unclassified", "perflava", "ARUP UnID 635", "Unclassified", 
    "Unclassified", "mirabilis", "australis", "Unclassified", 
    "oralis", "Unclassified", "Unclassified", "Unclassified", 
    "magna", "oralis", "sp. Marseille", "australis", 
    "Unclassified", "Unclassified", "australis"
  )
)

# Ensure that the manual taxonomy table has the same format as the cleaned taxonomy table
taxonomy_cleaned1 <- taxonomy_cleaned1 %>%
  rownames_to_column(var = "ASV_ID")

# Merge the manual taxonomy with the cleaned taxonomy
taxonomy_cleaned2 <- taxonomy_cleaned1 %>%
  left_join(manual_taxonomy, by = "ASV_ID", suffix = c("", ".manual")) %>%
  mutate(
    Genus = ifelse(!is.na(Genus.manual), Genus.manual, Genus),
    Species = ifelse(!is.na(Species.manual), Species.manual, Species)
  ) %>%
  dplyr::select(-Genus.manual, -Species.manual)

# Ensure that the Genus and Species columns are in the correct format
taxonomy_cleaned3 <- taxonomy_cleaned2 %>%
  mutate(
    Species = ifelse(
      Species == "Unclassified",
      paste0(Genus, "_Unclassified"),  # Incluir el género antes de Unclassified
      ifelse(!grepl("^g__", Species), paste0(Genus, " ", Species), Species)
    ),
    Species = gsub(" ", "_", Species)  # Reemplaza espacios por guiones bajos en el nombre de la especie
  )

# This step is to ensure that the Species column is formatted correctly
modern_classifications <- tibble::tribble(
  ~Genus, ~Kingdom, ~Phylum, ~Class, ~Order, ~Family,
  "g__Anaerococcus", "k__Bacteria", "p__Firmicutes", "c__Clostridia", "o__Clostridiales", "f__Peptoniphilaceae",
  "g__Atopobium", "k__Bacteria", "p__Actinobacteria", "c__Actinobacteria", "o__Coriobacteriales", "f__Coriobacteriaceae",
  "g__Catonella", "k__Bacteria", "p__Firmicutes", "c__Clostridia", "o__Clostridiales", "f__Lachnospiraceae",
  "g__Corynebacterium", "k__Bacteria", "p__Actinobacteria", "c__Actinobacteria", "o__Corynebacteriales", "f__Corynebacteriaceae",
  "g__Dialister", "k__Bacteria", "p__Firmicutes", "c__Negativicutes", "o__Selenomonadales", "f__Veillonellaceae",
  "g__Finegoldia", "k__Bacteria", "p__Firmicutes", "c__Clostridia", "o__Clostridiales", "f__Peptoniphilaceae",
  "g__Fusobacterium", "k__Bacteria", "p__Fusobacteria", "c__Fusobacteriia", "o__Fusobacteriales", "f__Fusobacteriaceae",
  "g__Gemella", "k__Bacteria", "p__Firmicutes", "c__Bacilli", "o__Gemellales", "f__Gemellaceae",
  "g__Helcococcus", "k__Bacteria", "p__Firmicutes", "c__Clostridia", "o__Clostridiales", "f__Peptoniphilaceae",
  "g__Megasphaera", "k__Bacteria", "p__Firmicutes", "c__Negativicutes", "o__Selenomonadales", "f__Veillonellaceae",
  "g__Mogibacterium", "k__Bacteria", "p__Firmicutes", "c__Clostridia", "o__Clostridiales", "f__Peptostreptococcaceae",
  "g__Mycoplasma", "k__Bacteria", "p__Tenericutes", "c__Mollicutes", "o__Mycoplasmatales", "f__Mycoplasmataceae",
  "g__Parvimonas", "k__Bacteria", "p__Firmicutes", "c__Clostridia", "o__Clostridiales", "f__Peptoniphilaceae",
  "g__Peptoniphilus", "k__Bacteria", "p__Firmicutes", "c__Clostridia", "o__Clostridiales", "f__Peptoniphilaceae",
  "g__Peptostreptococcus", "k__Bacteria", "p__Firmicutes", "c__Clostridia", "o__Clostridiales", "f__Peptostreptococcaceae",
  "g__Propionibacterium", "k__Bacteria", "p__Actinobacteria", "c__Actinobacteria", "o__Propionibacteriales", "f__Propionibacteriaceae",
  "g__Selenomonas", "k__Bacteria", "p__Firmicutes", "c__Negativicutes", "o__Selenomonadales", "f__Veillonellaceae",
  "g__Treponema", "k__Bacteria", "p__Spirochaetes", "c__Spirochaetia", "o__Spirochaetales", "f__Spirochaetaceae",
  "g__Veillonella", "k__Bacteria", "p__Firmicutes", "c__Negativicutes", "o__Selenomonadales", "f__Veillonellaceae",
  "g__Streptococcus", "k__Bacteria", "p__Firmicutes", "c__Bacilli", "o__Lactobacillales", "f__Streptococcaceae",  # Caso 1
  "g__Ureaplasma", "k__Bacteria", "p__Tenericutes", "c__Mollicutes", "o__Mycoplasmatales", "f__Mycoplasmataceae"  # Caso 2
)

# Ensure that the modern classifications table has the same format as the cleaned taxonomy table
taxonomy_cleaned3 <- taxonomy_cleaned3 %>%
  left_join(modern_classifications, by = "Genus", suffix = c("_old", "_modern")) %>%
  mutate(
    Kingdom = if_else(Genus != "g__Unclassified" & !is.na(Kingdom_modern), Kingdom_modern, Kingdom_old),
    Phylum = if_else(Genus != "g__Unclassified" & !is.na(Phylum_modern), Phylum_modern, Phylum_old),
    Class = if_else(Genus != "g__Unclassified" & !is.na(Class_modern), Class_modern, Class_old),
    Order = if_else(Genus != "g__Unclassified" & !is.na(Order_modern), Order_modern, Order_old),
    Family = if_else(Genus != "g__Unclassified" & !is.na(Family_modern), Family_modern, Family_old)
  ) %>%
  dplyr::select(-ends_with("_old"), -ends_with("_modern"))

# Restore the original rownames (ASV_ID)
taxonomy_cleaned3 <- taxonomy_cleaned3 %>%
  column_to_rownames(var = "ASV_ID")

# Check if there are any ASVs that are not classified at the species level
genera_multiple_classifications <- taxonomy_cleaned3 %>%
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

# Check if there are any genera with multiple classifications
genera_to_verify <- unique(genera_multiple_classifications$Genus)

# If there are genera with multiple classifications, we can filter the taxonomy table to show these genera
filtered_genera_taxonomy <- taxonomy_cleaned3[taxonomy_cleaned3$Genus %in% genera_to_verify, ]

# Display the filtered taxonomy for genera with multiple classifications
unique_classifications <- unique(filtered_genera_taxonomy[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")])

# arrange the taxonomy table to have the correct order of columns
taxonomy_cleaned <- taxonomy_cleaned3 %>%
  dplyr::select(Kingdom, Phylum, Class, Order, Family, Genus, Species)  # Reordenar las columnas

# remove the ASV (489ffd6307c1aef63be92187df3c15c8) that is not classified at the species level
# Directly remove the ASV using row names
taxonomy_cleaned <- taxonomy_cleaned[rownames(taxonomy_cleaned) != "489ffd6307c1aef63be92187df3c15c8", ]
print(tail(taxonomy_cleaned))

write_csv(taxonomy_cleaned, 'P16S/qiime2_qza_vag_tong_ggext_newseqs/taxonomy_cleaned1sep24_ggext_wUnclassified.csv', quote = "none")
# Convert taxonomy_cleaned to tax_table to phyloseq objects
taxonomy_cleaned <- tax_table(as.matrix(taxonomy_cleaned))

# PRUNE taxonomy after ASV table filtering
common_taxa <- intersect(taxa_names(taxonomy_cleaned), taxa_names(ASV))
TAXA_clean <- prune_taxa(common_taxa, taxonomy_cleaned)
OTU_clean <- prune_taxa(common_taxa, ASV)

# Create the phyloseq object
phy1 <- phyloseq(TAXA_clean, OTU_clean)
phy1

# list of samples to remove, "M469MTNG" and "M469MVAG" will be removed since it seems they were flipped in the labeling, ordination shows them as vaginal the oral samples and as oral the vaginal samples. To avoid confusion we will eliminate them
samples_to_remove <- c("NEGCNTRLINXA", "NEGCNTRLINXD", 
                       "Undetermined", "M469MTNG","M469MVAG")

phy1 <- prune_samples(!(sample_names(phy1) %in% samples_to_remove), phy1)

write.csv(tax_table(phy1), "P16S/qiime2_qza_vag_tong_ggext_newseqs/tax_table_cleaned1sep24_ggext_wUnclassified.csv")

#### Filtering ab and prev ####
# Filter taxa based on abundance and prevalence
phy_filtered <- prune_taxa(taxa_sums(phy1) > 10, phy1) 

asv_presence <- apply(otu_table(phy_filtered), 1, function(x) sum(x > 0))

phy1 <- prune_taxa(asv_presence > 2, phy_filtered)

tax_table <- as.data.frame(tax_table(phy1))

##### Verify if there are any ASVs that are not classified at the species level ####
# Check if there are any ASVs that are not classified at the species level

unclassified_asvs <- tax_table(phy1) %>% 
  as.data.frame() %>%
  filter(grepl("Unclassified", Species)) %>%
  rownames()

otu_table <- otu_table(phy1)

unclassified_otu_table <- otu_table[unclassified_asvs, ]

unclassified_taxonomy <- tax_table(phy1) %>%
  as.data.frame() %>%
  filter(rownames(.) %in% unclassified_asvs) %>%
  dplyr::select(Family, Genus)

# Calculate the number of samples that contain each ASV without classification
unclassified_counts <- rowSums(unclassified_otu_table)

# Calculate the total number of reads for each ASV without classification
unclassified_sample_counts <- apply(unclassified_otu_table, 1, function(x) sum(x > 0))

# Convert the counts and sample counts to a data frame
unclassified_summary <- data.frame(
  ASV_ID = unclassified_asvs,
  Family = unclassified_taxonomy$Family,
  Genus = unclassified_taxonomy$Genus,
  Total_Reads = unclassified_counts,
  Sample_Count = unclassified_sample_counts
)

# Sort the summary by Total_Reads in descending order
unclassified_summary %>% 
  arrange(desc(Total_Reads))

# list of ASVs that are not classified at the species level
ASV_IDs <- c(
  "5f32e0d308736382d1e4a2cc3a800f03", "3314486729173c20661f139588a86762", 
  "f8d97633268d0cb33a4961ac1d1991e3", "c577d055d016e6c2dfaa9ca155f9f51c", 
  "5a697167dd760a129253902fa6390f95", "7d062509323aa7529dacc301efa91a1e", 
  "eed45c46ab8412632a530e3b2615671e", "f1a1b898667f658b6ba0abe243575021", 
  "61b219d4610c015b5faa10bb8e305700", "821045e2195d26f66a1fe7605b57a6c1", 
  "30f38eea16b97b9f820b34ad88223408", "132c6edec789232615eb69bf39db8b2a", 
  "a019d0954fc18f6c7cff143043b12a87", "176d59c6cca071c80c2b60e9393aa6fc", 
  "2a8e9f5d1713e964b80a0806830fe42b", "969f4202b63057ecc32b1233c3ae2f7c", 
  "4968584ae5a1560f360301ebdffd04f1", "2c5b4b7e04c1cc6fdec353a011a50795", 
  "2bc2f4ce95b919a17db546dd04aa5b09", "e1bf1a1bd4c69c93972eb2210d26ca09", 
  "70c88fdd526ea3b46145647fcddd51dd", "8b5799bb6d57eb3c327a2b21e3a0c3fe", 
  "86e514302b29b2a3bf7b6e62827fdf67", "7da5699cd3ff61097747c640b596cb85", 
  "2569ad3a1d65c0b7e46b029fedc81740",
  "a6f1a209e0780747e50ccf8117c0c854", "04e96f8349e34826431baa6190c35683",
  "06790a632d34be8db2eb09712082abd0", "0cab42d78e6dbabb304a0d5dfe88f356",
  "5158b864829c2ec1c77c13e2e993885a", "163d0544a78051aa4fa36b65d02c17b6",
  "a48572566c42efef9aa4d0c6641561eb", "514b75f20ccc23afa486f434df15f3fe",
  "08968b4577cd35855818f4fe99ac1859", "521d07ab040f6421c93b514d6750fb65",
  "58fd461434d936e5aab4c275d6f6cde7", "2bbcd6a04228ca5ccf105b904828c529",
  "843344d09a7c84470c990f1b6170032b", "7faa1e175e9fb5b005644773c97346ac",
  "930e445523e9c5d5aad858b2e20ba1cf", "1452b3e007e3bab260efe9108618bdf6",
  "9260ffbad675235f1f37dff12ba86038", "4c31e2d1733ca82e1cf384498344ec3b",
  "2239edacecc6a136f597624c94c63086"
)

# Extract the ASV IDs from the phyloseq object
taxonomy_subset <- tax_table(phy1)[ASV_IDs, ]

####METADATA####

met <- read.csv("met/met_clean_vag_tong.csv", header = TRUE, row.names = 1)
cst <- read.csv("P16S/R/CST/results/vag_pree1.csv", header = TRUE, row.names = 1)
IL6 <- read_excel("docs/UPDATED_IL6_Concentrations_PreE.xlsx", sheet = "SamplesRnd1", col_names = TRUE)

# Add the column 'Conc' from IL6 to met data for Vaginal

met$Conc[met$body_site == "Vaginal"] <- IL6$Conc[match(met$study_id[met$body_site == "Vaginal"], IL6$Sample)]
met$Conc <- as.numeric(met$Conc)

#change name Conc by vag_IL6
colnames(met)[colnames(met) == "Conc"] <- "vag_IL6"

cst <- cst %>%
  mutate(CST5 = case_when(
    CST %in% c("IV-A", "IV-B", "IV-C") ~ "IV",
    TRUE ~ CST
  ))

# Convert rownames of met to a new column for merging
cst$study_id_formatted <- rownames(cst)

cst_subset <- cst[, c("study_id_formatted","CST5", "CST", "subCST")]

merged_data <- merge(met, cst_subset, by = "study_id_formatted", all.x = TRUE)

TAXA_clean=as.matrix(tax_table(phy1))
OTU_clean=as.matrix(otu_table(phy1, taxa_are_rows = TRUE))

# Arrange the metadata to match the order of the ASV table
met_ordered <- merged_data[match(colnames(OTU_clean), merged_data$study_id_formatted), ]

# Assign the row names of met_ordered to the sample names of phy1
rownames(met_ordered) <- colnames(OTU_clean)

# Add the metadata to the phyloseq object
rownames(met_ordered) <- sample_names(phy1)
met_ordered1 <- sample_data(met_ordered)

# Create phyloseq object with the cleaned taxonomy, OTU table, and metadata
PHYraw <- phyloseq(TAXA_clean, OTU_clean, met_ordered1)

save(PHYraw,file="P16S/R/phy_objects/PHYraw.RData")
load("~/Dropbox (UMass Medical School)/Ana_projects/2_PREECLAMPSIA/R/phy_objects/PHYraw.RData")
load("/Users/danielavargasrobles/UMass Medical School Dropbox/Daniela Vargas Robles/Ana_projects/2_PREECLAMPSIA/P16S/R/phy_objects/PHYraw.RData")
#### Rarefaction ####
# Calculate the number of sequences per sample
sumatoria <- sample_sums(PHYraw)
min(sumatoria)

sumatoria_ordenada <- sort(sumatoria)
print(head(sumatoria_ordenada, 30))

# Add the sum of sequences as a new column in the metadata
sample_data(PHYraw)$num_seqs <- sumatoria

seq_to_check=data.frame(sample_data(PHYraw))

seq_to_check %>% 
  dplyr::select(num_seqs) %>% 
  arrange(-num_seqs)

#write.csv(seq_to_check, "metadata/met_13oct22_9_clean_24april24_casos_control_column_nseqs.csv", row.names = TRUE)

#rarefaction
PHYrar=rarefy_even_depth(PHYraw, sample.size = 8634,#14567,#71815,#38854,#16320,#5430,#3489,
                         rngseed = 711, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

#One note that can be useful, I also use the estimation of Good's indexes (library(QsRutils)) to evaluate, for instance, if the rarefaction approach is adequate, by looking at the median and mean values
#you can see that both mean and median values dropped a bit after rarefaction. If they remain >90% its ok otherwise the rarefied samples are not reflecting adequately the original taxa diversity
not_rare=summary(goods(otu_table(PHYraw)))
not_rare
rare=summary(goods(otu_table(PHYrar)))
rare

#save rarefied object
save(PHYrar,file="P16S/R/phy_objects/PHYrar.RData")
met=sample_data(PHYrar)

#### ALPHA diversity indices#####
#1. extract ASV table 
OTU1 = as(otu_table(PHYrar), "matrix")
OTUdf = as.data.frame(OTU1)

#Calculating 4 different alpha diversity metrics
#at asv level
shannon = diversity(OTUdf, index = "shannon", MARGIN = 2, base = exp(1))
simpson = diversity(OTUdf, index = "simpson", MARGIN = 2)
#chao = apply(OTUdf, 2, chao1)
chao=estimate_richness(PHYrar, split = TRUE,measures = "Chao1")
obs = apply(OTUdf, 2, function(x){
  nws = length(which(x>0))
  return(nws)
})

div_asv = data.frame(shannon = shannon, 
                     simpson = simpson,
                     chao1 = chao,
                     observed = obs)



#merging columns of alpha metrics with the rest of the metadata by SampleID
md1 = sample_data(PHYrar)
md0=as.data.frame(as(md1, "matrix"))
write.csv(md0, "md0.temp", row.names = TRUE)
md = read.csv("md0.temp")

asv=cbind(md,div_asv)
rownames(asv) = sample_names(PHYrar)

sampledata = sample_data(asv)

#Re doing phyloseq object with new added columns (alpha div metrics) in metadata
PHYrar_alpha_ggext_subset02_newseqs = phyloseq(otu_table(PHYrar),tax_table(PHYrar), sampledata)#phy_tree(PHYrar) )
PHYrar_alpha_ggext_subset02_newseqs
#### Save ####
save(PHYrar_alpha_ggext_subset02_newseqs,file="P16S/R/phy_objects/PHYrar_alpha_ggext_subset02_newseqs.RData")

#extract metadata as data frame
met_phylo <- data.frame(sample_data(PHYrar_alpha_ggext_subset02_newseqs))
write.csv(met_phylo, "met/met_PHYrar_alpha_ggext_subset02_newseqs.csv")

