# Cargar librerías necesarias
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
# Configurar el directorio de trabajo
setwd("~/Dropbox (UMass Medical School)/Ana_projects/2_PREECLAMPSIA/P16S/meconium/")

# Cargar datos de ASV
#phy <- qza_to_phyloseq(features = "qiime2_qza_vag_tong_gg2/02vag_tong/table_vag_tong_subsample02.qza")
#new seqs
phy <- qza_to_phyloseq(features = "qiime2_files/table_demux.qza")
phy
dim(phy)
ASV <- otu_table(phy)

write.csv(ASV, "qiime2_files/ASV_table_raw.csv")

data<-as.data.frame(t(ASV))

dim(data)


# Cargar datos de taxonomía gg2
taxonomy <- read_tsv('qiime2_files/taxonomy/taxonomy.tsv', show_col_types = FALSE)

# 2. Limpiar tabla de taxonomía
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
      Species  # Si ya hay Species definida, no la toques
    )
  ) %>%
  column_to_rownames(var = "Feature ID") 

# 3. Filtro paso a paso
taxonomy_filtered <- taxonomy_cleaned %>%
  filter(Kingdom != "d__Archaea") %>%
  filter(Kingdom != "Unassigned") %>%
  filter(Phylum != "p__Unclassified") %>%
  filter(Class != "c__Unclassified") %>%
  filter(Order != "o__Unclassified") %>%
  filter(!if_any(everything(), ~ grepl("mitochondria", ., ignore.case = TRUE)))

# 4. Vista previa
head(taxonomy_filtered)
dim(taxonomy_filtered)
#ver mas filas
write.csv(taxonomy_filtered, "tables_clean/tax_cleaned_meconium.csv")

#verificar si hay géneros con múltiples clasificaciones

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

# Extraer los nombres de los géneros con múltiples clasificaciones
genera_to_verify <- unique(genera_multiple_classifications$Genus)

# # Filtrar la tabla original de taxonomía para mostrar las clasificaciones de los géneros con múltiples clasificaciones
# filtered_genera_taxonomy <- taxonomy_filtered[taxonomy_filtered$Genus %in% genera_to_verify, ]
# 
# # Obtener las clasificaciones únicas para estos géneros
# unique_classifications <- unique(filtered_genera_taxonomy[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")])
# 
# # Ordenar por el nombre del género
# unique_classifications %>%
#   arrange(Genus)
# 
# # Revisar los resultados
# print(tail(taxonomy_filtered))


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
# Lista de muestras a eliminar
# samples_to_remove <- c("NEGCNTRLINXA", "NEGCNTRLINXD", 
#                        "Undetermined", "M469MTNG","M469MVAG")

# Eliminar las muestras del objeto phyloseq
#phy1 <- prune_samples(!(sample_names(phy1) %in% samples_to_remove), phy1)


#### Filtering ab and prev ####
#Se aplican los mismo filtros qye para WGS fecal de la madre

# 1. Tu objeto phyloseq se llama phy1
# 2. Primero pasamos a abundancias relativas
phy1_rel <- transform_sample_counts(PHYrar, function(x) x / sum(x))

# 3. Filtro por abundancia relativa > 0.01
phy1_rel_filtered <- prune_taxa(taxa_sums(phy1_rel) > 0.01, phy1_rel)

# 4. Filtro por prevalencia > 1 muestra
taxa_presence <- apply(X = otu_table(phy1_rel_filtered), MARGIN = 1, FUN = function(x) sum(x > 0))
phy1_rel_filtered_prev <- prune_taxa(taxa_presence > 1, phy1_rel_filtered)

# 5. Como quieres regresar a conteos (no abundancias relativas), aplicamos los filtros a phy1 original:
# Mantén solo los taxones que sobrevivieron ambos filtros
phy1_filtered_final <- prune_taxa(taxa_names(phy1_rel_filtered_prev), phy1)

# 6. Resultado final
phy1_filtered_final
#extrat otu table to see the colnames
OTU_table <- as.data.frame(otu_table(phy1_filtered_final))
OTU_table

####METADATA####
# 1. Cargar la metadata
met <- read.csv("~/Dropbox (UMass Medical School)/Ana_projects/2_PREECLAMPSIA/met/met_clean0.csv")

# 2. Crear el ID formateado (agregar 'B' al principio)
met$study_id_formatted <- paste0("B", met$study_id)

# 3. Extraer tablas de taxonomía y ASV del phyloseq filtrado
TAXA_clean <- as.matrix(tax_table(phy1_filtered_final))
OTU_clean <- as.matrix(otu_table(phy1_filtered_final, taxa_are_rows = TRUE))

# 4. Reordenar la metadata para que coincida con el orden de las columnas de la tabla de ASVs
# (colnames(OTU_clean) son los nombres de las muestras)
met_ordered <- met[match(colnames(OTU_clean), met$study_id_formatted), ]

# 5. Asignar los nombres de las muestras (columnas de OTU_clean) como rownames de la metadata
rownames(met_ordered) <- colnames(OTU_clean)

# 6. Crear nuevo objeto phyloseq
PHY_fil <- phyloseq(
  otu_table(OTU_clean, taxa_are_rows = TRUE),
  tax_table(TAXA_clean),
  sample_data(met_ordered)
)

# 7. Resultado final
PHYraw
#extract metadata
met_phylo <- data.frame(sample_data(PHY_fil))
met_phylo
save(PHY_fil,file="phy_obj/PHY_fil.RData")
#load("R/phy_objects/PHYraw.RData")

#### Rarefaction ####
# Calcular la suma de secuencias por muestra
sumatoria <- sample_sums(PHY_fil)
min(sumatoria)

sumatoria_ordenada <- sort(sumatoria)
print(head(sumatoria_ordenada, 30))

# Agregar la suma de secuencias como una nueva columna en la metadata
sample_data(PHY_fil)$num_seqs <- sumatoria

# Mostrar la metadata con la nueva columna
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

#extract metadata as data frame
met_phylo <- data.frame(sample_data(PHY_fil_rar_alpha_meco))
write.csv(met_phylo, "tables_clean/met_PHY_fil_rar_alpha_meco.csv")
#extract asv table and save
ASV_table <- as.data.frame(otu_table(PHY_fil_rar_alpha_meco))
ASV_table
write.csv(ASV_table, "tables_clean/ASV_table_meconium.csv")
