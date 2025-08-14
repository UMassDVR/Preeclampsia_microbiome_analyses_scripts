# Genaration of phyloseq object from the women stool sample from WGS 

# LOAD LIBRARIES
library(forcats)
library(dplyr)
library(data.table) #for fread function
library(purrr)
library(anthro)
library(tidyr)

####LOAD####
setwd("~/Dropbox (UMass Medical School)/Ana_projects/2_PREECLAMPSIA/")

#### metadata REDCap
met=read.csv("met/met_clean0.csv")

# remove duplicated study_id
met <- met[!duplicated(met$study_id), ]

# Load alpha  div indices using relative abundance from metaphlan
obs <- read.csv("/Users/danielavargasrobles/Dropbox (UMass Medical School)/Ana_projects/2_PREECLAMPSIA/WGS/seqs/mom/rel_ab/gtdb/diversity_analysis/gtdb_merged_rel_ab_table_richness.tsv", sep="\t", header=TRUE)
sha <- read.csv("/Users/danielavargasrobles/Dropbox (UMass Medical School)/Ana_projects/2_PREECLAMPSIA/WGS/seqs/mom/rel_ab/gtdb/diversity_analysis/gtdb_merged_rel_ab_table_shannon.tsv", sep="\t", header=TRUE)
simp <- read.csv("/Users/danielavargasrobles/Dropbox (UMass Medical School)/Ana_projects/2_PREECLAMPSIA/WGS/seqs/mom/rel_ab/gtdb/diversity_analysis/gtdb_merged_rel_ab_table_simpson.tsv", sep="\t", header=TRUE)


obs$study_id <- rownames(obs)
sha$study_id <- rownames(sha)
simp$study_id <- rownames(simp)

#Combine 
div <- merge(obs, sha, by = "study_id")
div <- merge(div, simp, by = "study_id")


#### VANNI phyloseq object############

#filter 0.01% ab y 2 prev
library("phyloseq")
library("gtools")
library("stringr")

# Import Metaphlan
# Read metaphlan and skip the first row
mtph_dat <- read.csv("/Users/danielavargasrobles/Dropbox (UMass Medical School)/Ana_projects/2_PREECLAMPSIA/WGS/seqs/mom/rel_ab/gtdb/gtdb_merged_rel_ab_table.txt",sep="\t",skip = 1)

mtph_dat <-  mtph_dat[grep("s__",mtph_dat$clade_name),]

# Split into taxonomy table
tax_dt <-  str_split_fixed(mtph_dat$clade_name,pattern = "\\;", n = 7) %>% data.frame()
colnames(tax_dt) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

istk <- which(tax_dt$Species!="s__") 
tax_dt <- tax_dt[istk,]
mtph_dat <- mtph_dat[istk,]

# Remove the *__  from all taxa levels
tax_dt[] <- lapply(tax_dt, function(x) gsub("^.{0,3}", "", x)) 
rownames(tax_dt) <-  tax_dt$Species

# Relative abundance table
rel_dt <- mtph_dat[,2:ncol(mtph_dat)] 

#colSums(rel_dt)
rownames(rel_dt) <-  tax_dt$Species
colnames(rel_dt) <- gsub("^M|Stool_mph.txt_gtdb","",colnames(rel_dt))

meta_dt <- met_div_rel2
#meta_dt <- read.csv("/Users/danielavargasrobles/Dropbox (UMass Medical School)/Ana_projects/2_PREECLAMPSIA/met/met_clean0.csv")
#met$study_id <- paste0("X", met$study_id)
#met_g$study_id <- paste0("X", met_g$study_id)


head(meta_dt)
rownames(meta_dt) <-  meta_dt$study_id

match(rownames(meta_dt),names(rel_dt))

phy_mtph <-  phyloseq(otu_table(rel_dt,taxa_are_rows = T),
                      tax_table(as.matrix(tax_dt)),
                      sample_data(meta_dt))

sample_sums(phy_mtph) 
table(sample_sums(phy_mtph) > 0)


phy_mtph <-  prune_samples(sample_sums(phy_mtph)>0, phy_mtph)
phy_mtph <-  prune_taxa(taxa_sums(phy_mtph)>0, phy_mtph)

# remove from phy_mtph taxa with an abundance of less than 0.01%
phy_mtph1 <- prune_taxa(taxa_sums(phy_mtph) > 0.01, phy_mtph)
print(phy_mtph1)

# Calculate the number of samples where each taxon is present
taxa_presence <- apply(X = phy_mtph1@otu_table, MARGIN = 1, FUN = function(x) sum(x > 0))

# Remove taxa present in only 1 sample
phy_mtph_gtdb <- prune_taxa(taxa_presence > 1, phy_mtph1)
phy_mtph_gtdb

# update the phyloseq object metadata with this one met_div


saveRDS(phy_mtph_gtdb,"/Users/danielavargasrobles/Dropbox (UMass Medical School)/Ana_projects/2_PREECLAMPSIA/WGS/obj/phy_mtph_gtdb.rds")

#extract metadata
meta_phy <- phy_mtph_gtdb@sam_data

