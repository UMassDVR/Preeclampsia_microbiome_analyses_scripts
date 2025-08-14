####LOAD####
# This script creates demographic tables.
# The tables are created using the metadata from the phyloseq object and saved as CSV files.
# The script includes functions to create tables for all samples, Hispanic samples, and Non-Hispanic White samples.
# The tables include various demographic variables and are stratified by study group or severity.
# The script also includes a section for creating tables based on severity with the actual number of infant stool sequenced samples

# Load libraries
library(phyloseq)
library(dplyr)
library(tableone)

# Set working directory (ensure this path is correct for your system)
setwd("~/Dropbox (UMass Medical School)/Ana_projects/2_PREECLAMPSIA/P16S/")

met <- read.csv("/Users/danielavargasrobles/UMass Medical School Dropbox/Daniela Vargas Robles/Ana_projects/2_PREECLAMPSIA/met/met_clean_vag_tong.csv", header = TRUE, row.names = 1)

### Create table1 ####

# 2. FUNCTION DEFINITION
# A function to process a data frame subset, create a table, and save it.
create_and_save_table <- function(data_subset, variables, file_name) {
  
  # Ensure the directory exists
  dir.create(dirname(file_name), showWarnings = FALSE, recursive = TRUE)
  
  # Reorder factor levels by frequency
  reorder_by_freq <- function(df, col_name) {
    if (!col_name %in% names(df)) return(df[[col_name]]) # Return NULL if col doesn't exist
    counts <- df %>% count(.data[[col_name]]) %>% arrange(desc(n))
    factor(df[[col_name]], levels = counts[[col_name]])
  }
  
  # Apply reordering
  metadata <- data_subset
  metadata$race <- reorder_by_freq(metadata, "race")
  metadata$mod <- reorder_by_freq(metadata, "mod")
  metadata$ethnicity <- reorder_by_freq(metadata, "ethnicity")
  
  # Rename columns for the final table
  metadata_renamed <- metadata %>%
    rename(
      "Age" = age,
      "Race" = race,
      "Ethnicity" = ethnicity,
      "BMI" = bmi,
      "Diabetes" = diabetes,
      "Delivery mode" = mod,
      "Gestational age (weeks)" = gestational_age,
      "Term status (<37 weeks)" = term_status,
      "Mother ICU" = mom_icu,
      "Antibiotics during pregnancy" = antibiotics_during_pregnancy,
      "Antibiotics during c-section" = antibiotic_during_cs,
      "Newborn sex" = sex_of_baby,
      "Newborn weight" = baby_weight,
      "Newborn ICU" = baby_nicu
    )
  
  # Create and save the table
  summary_table <- CreateTableOne(vars = variables, strata = "study_group", data = metadata_renamed)
  table_df <- print(summary_table, quote = FALSE, noSpaces = TRUE, printToggle = FALSE, showAllLevels = FALSE, includeNA = TRUE)
  
  # Add asterisk to significant p-values
  if ("p" %in% colnames(table_df)) {
    p_values_num <- suppressWarnings(as.numeric(table_df[, "p"]))
    significant_rows <- (!is.na(p_values_num) & p_values_num < 0.050) | grepl("<", table_df[, "p"])
    original_p <- table_df[significant_rows, "p"]
    table_df[significant_rows, "p"] <- paste0(original_p, "*")
  }
  
  # Save the modified data frame
  write.csv(table_df, file = file_name)
  cat("Successfully created and saved table to:", file_name, "\n")
}


# 3. EXECUTION
# Define the variables for each table
vars_all <- c("Age", "Race", "Ethnicity", "BMI",  "Diabetes", "Gestational age (weeks)", "Term status (<37 weeks)", "Delivery mode", "Antibiotics during pregnancy", "Mother ICU", "Newborn sex", "Newborn weight", "Newborn ICU")
vars_hispanic <- c("Age", "Race", "BMI",  "Diabetes", "Gestational age (weeks)", "Term status (<37 weeks)","Delivery mode","Antibiotics during c-section", "Mother ICU", "Newborn sex", "Newborn weight", "Newborn ICU")
vars_nhw <- c("Age", "Race", "BMI", "Diabetes", "Gestational age (weeks)","Term status (<37 weeks)", "Delivery mode", "Antibiotics during c-section", "Mother ICU", "Newborn sex", "Newborn weight", "Newborn ICU")

# --- Table 1: All Vaginal Samples ---
data_subset_1 <- filter(met, body_site == "Vaginal")
create_and_save_table(data_subset_1, vars_all, "results/tables/table1_mom_baby.csv")

# --- Table 2: Hispanic Vaginal Samples ---
data_subset_2 <- filter(met, body_site == "Vaginal" & ethnicity == "Hispanic/Latina")
create_and_save_table(data_subset_2, vars_hispanic, "results/tables/table1_mom_baby_hispanics.csv")

# --- Table 3: Non-Hispanic White Vaginal Samples ---
data_subset_3 <- filter(met, body_site == "Vaginal" & race == "White/Caucasian" & ethnicity == "Not Hispanic/Latina")
create_and_save_table(data_subset_3, vars_nhw, "results/tables/table1_mom_baby_non_hispanics_whites.csv")






#### create table by SEVERITY- all samples #####
met <- read.csv("/Users/danielavargasrobles/UMass Medical School Dropbox/Daniela Vargas Robles/Ana_projects/2_PREECLAMPSIA/met/met_clean_vag_tong.csv", header = TRUE, row.names = 1)


# 2. FUNCTION DEFINITION
# A function to process a data frame subset, create a table, and save it.
create_and_save_table <- function(data_subset, variables, file_name) {
  
  # Ensure the directory exists
  dir.create(dirname(file_name), showWarnings = FALSE, recursive = TRUE)
  
  # Reorder factor levels by frequency
  reorder_by_freq <- function(df, col_name) {
    if (!col_name %in% names(df)) return(df[[col_name]])
    counts <- df %>% count(.data[[col_name]]) %>% arrange(desc(n))
    factor(df[[col_name]], levels = counts[[col_name]])
  }
  
  # Apply reordering
  metadata <- data_subset
  metadata$race <- reorder_by_freq(metadata, "race")
  metadata$mod <- reorder_by_freq(metadata, "mod")
  metadata$ethnicity <- reorder_by_freq(metadata, "ethnicity")
  
  # Rename columns for the final table
  metadata_renamed <- metadata %>%
    rename(
      "Age" = age,
      "Race" = race,
      "Ethnicity" = ethnicity,
      "BMI" = bmi,
      "Diabetes" = diabetes,
      "Delivery mode" = mod,
      "Gestational age (weeks)" = gestational_age,
      "Term status (<37 weeks)" = term_status,
      "Mother ICU" = mom_icu,
      "Antibiotics during pregnancy" = antibiotics_during_pregnancy,
      "Antibiotics during c-section" = antibiotic_during_cs,
      "Newborn sex" = sex_of_baby,
      "Newborn weight" = baby_weight,
      "Newborn ICU" = baby_nicu
    )
  
  # Create and save the table, stratified by "severe"
  summary_table <- CreateTableOne(vars = variables, strata = "severe", data = metadata_renamed)
  table_df <- print(summary_table, quote = FALSE, noSpaces = TRUE, printToggle = FALSE, showAllLevels = FALSE, includeNA = TRUE)
  
  # Add asterisk to significant p-values
  if ("p" %in% colnames(table_df)) {
    p_values_num <- suppressWarnings(as.numeric(table_df[, "p"]))
    significant_rows <- (!is.na(p_values_num) & p_values_num < 0.050) | grepl("<", table_df[, "p"])
    original_p <- table_df[significant_rows, "p"]
    table_df[significant_rows, "p"] <- paste0(original_p, "*")
  }
  
  # Save the modified data frame
  write.csv(table_df, file = file_name)
  cat("Successfully created and saved table to:", file_name, "\n")
}


# 3. EXECUTION
# Define the variables for the table
vars_sev <- c("Age", "Race", "Ethnicity", "BMI", "Diabetes", "Gestational age (weeks)", "Term status (<37 weeks)", "Delivery mode", "Antibiotics during pregnancy", "Mother ICU", "Newborn sex", "Newborn weight", "Newborn ICU")

# --- Table: Vaginal Samples, stratified by severity (excluding controls) ---
# The filter assumes your 'met' dataframe has a 'severe_3cat' column
data_subset_4 <- filter(met, body_site == "Vaginal" & severe_3cat != "control")
create_and_save_table(data_subset_4, vars_sev, "results/tables/table1_mom_baby_severe.csv")






#### Create table for the number Sequenced samples of infant stool #####
setwd("~/Dropbox (UMass Medical School)/Ana_projects/2_PREECLAMPSIA/P16S/")

# Load the phyloseq object
load("R/phy_objects/PHYrar_alpha_ggext_subset02_newseqs.RData")
base_phy_object <- PHYrar_alpha_ggext_subset02_newseqs

base_phy_object <- subset_samples(base_phy_object, study_id %in% c(401, 402, 403, 405, 406, 407, 408, 409, 410, 411, 413, 414, 
                                                                   419, 421, 422, 423, 424, 425, 426, 427, 428, 429, 430, 431, 
                                                                   432, 433, 434, 435, 436, 437, 438, 439, 440, 442, 443, 444, 
                                                                   445, 446, 447, 448, 450, 451, 452, 454, 455, 457, 458, 459, 
                                                                   460, 461, 463, 464, 465, 467, 470, 471, 472, 473, 474, 477, 
                                                                   478, 479, 480, 481, 482, 483, 484, 485, 486, 487, 488, 489, 
                                                                   490, 492, 493, 494, 495, 498, 500))                                                                  
                                                                  


# 2. FUNCTION DEFINITION
# A function to process a phyloseq subset, create a table, and save it.
create_and_save_table <- function(phy_subset, variables, file_name) {
  
  # Extract metadata
  metadata <- data.frame(sample_data(phy_subset))
  
  # Reorder factor levels by frequency
  reorder_by_freq <- function(df, col_name) {
    counts <- df %>% count(.data[[col_name]]) %>% arrange(desc(n))
    factor(df[[col_name]], levels = counts[[col_name]])
  }
  metadata$race <- reorder_by_freq(metadata, "race")
  metadata$mod <- reorder_by_freq(metadata, "mod")
  metadata$ethnicity <- reorder_by_freq(metadata, "ethnicity")
  
  # Rename columns for the final table
  metadata_renamed <- metadata %>%
    rename(
      "Age" = age,
      "Race" = race,
      "Ethnicity" = ethnicity,
      "BMI" = bmi,
      "Diabetes" = diabetes,
      "Delivery mode" = mod,
      "Gestational age (weeks)" = gestational_age,
      "Term status (<37 weeks)" = term_status,
      "Mother ICU" = mom_icu,
      "Antibiotics during pregnancy" = antibiotics_during_pregnancy,
      "Antibiotics during c-section" = antibiotic_during_cs,
      "Newborn sex" = sex_of_baby,
      "Newborn weight" = baby_weight,
      "Newborn ICU" = baby_nicu    )
  
  # Create and save the table
  summary_table <- CreateTableOne(vars = variables, strata = "study_group", data = metadata_renamed)
  table_df <- print(summary_table, quote = FALSE, noSpaces = TRUE, printToggle = FALSE, showAllLevels = FALSE, includeNA = TRUE)
  
  # --- Add asterisk to significant p-values ---
  # Check if the 'p' column exists
  if ("p" %in% colnames(table_df)) {
    # Convert p-value column to numeric, suppressing warnings for non-numeric entries like "<0.001"
    p_values_num <- suppressWarnings(as.numeric(table_df[, "p"]))
    
    # Identify rows where p < 0.05 or where the text contains "<"
    significant_rows <- (p_values_num < 0.050 & !is.na(p_values_num)) | grepl("<", table_df[, "p"])
    
    # Add an asterisk to those rows. Ensure we don't add to empty p-value cells.
    original_p <- table_df[significant_rows, "p"]
    table_df[significant_rows, "p"] <- paste0(original_p, "*")
  }
  
  # Save the modified data frame
  write.csv(table_df, file = file_name)
  
  cat("Successfully created and saved table with significance marks to:", file_name, "\n")
}


# 3. EXECUTION
# Define the variables for each table
vars_all <- c("Age", "Race", "Ethnicity", "BMI",  "Diabetes", "Gestational age (weeks)", "Term status (<37 weeks)", "Delivery mode", "Antibiotics during pregnancy", "Mother ICU", "Newborn sex", "Newborn weight", "Newborn ICU")
vars_hispanic <- c("Age", "Race", "BMI",  "Diabetes", "Gestational age (weeks)", "Term status (<37 weeks)","Delivery mode","Antibiotics during c-section", "Mother ICU", "Newborn sex", "Newborn weight", "Newborn ICU")
vars_nhw <- c("Age", "Race", "BMI", "Diabetes", "Gestational age (weeks)","Term status (<37 weeks)", "Delivery mode", "Antibiotics during c-section", "Mother ICU", "Newborn sex", "Newborn weight", "Newborn ICU")



# --- Table 1: All Vaginal Samples (Non-control) ---
phy_subset_1 <- subset_samples(base_phy_object, body_site == "Vaginal")
create_and_save_table(phy_subset_1, vars_all, "results/tables/table1_mom_baby_meconium_sequenced.csv")

# --- Table 2: Hispanic Vaginal Samples ---
phy_subset_2 <- subset_samples(base_phy_object, body_site == "Vaginal" & ethnicity == "Hispanic/Latina")
create_and_save_table(phy_subset_2, vars_hispanic, "results/tables/table1_mom_baby_hispanics_meconium_sequenced.csv")

# --- Table 3: Non-Hispanic White Vaginal Samples ---
phy_subset_3 <- subset_samples(base_phy_object, body_site == "Vaginal" & race == "White/Caucasian" & ethnicity == "Not Hispanic/Latina")
create_and_save_table(phy_subset_3, vars_nhw, "results/tables/table1_mom_baby_non_hispanics_whites_meconium_sequenced.csv")


#### Create table by SEVERITY for the number Sequenced samples of infant stool #####
# Set working directory (ensure this path is correct for your system)
setwd("~/Dropbox (UMass Medical School)/Ana_projects/2_PREECLAMPSIA/P16S/")

# Load the phyloseq object
load("R/phy_objects/PHYrar_alpha_ggext_subset02_newseqs.RData")
base_phy_object <- PHYrar_alpha_ggext_subset02_newseqs

base_phy_object <- subset_samples(base_phy_object, study_id %in% c(401, 402, 403, 405, 406, 407, 408, 409, 410, 411, 413, 414, 
                                                                   419, 421, 422, 423, 424, 425, 426, 427, 428, 429, 430, 431, 
                                                                   432, 433, 434, 435, 436, 437, 438, 439, 440, 442, 443, 444, 
                                                                   445, 446, 447, 448, 450, 451, 452, 454, 455, 457, 458, 459, 
                                                                   460, 461, 463, 464, 465, 467, 470, 471, 472, 473, 474, 477, 
                                                                   478, 479, 480, 481, 482, 483, 484, 485, 486, 487, 488, 489, 
                                                                   490, 492, 493, 494, 495, 498, 500))                                                                  

# 2. FUNCTION DEFINITION
# A function to process a phyloseq subset, create a table, and save it.
create_and_save_table <- function(phy_subset, variables, file_name) {
  
  # Extract metadata
  metadata <- data.frame(sample_data(phy_subset))
  
  # Reorder factor levels by frequency
  reorder_by_freq <- function(df, col_name) {
    counts <- df %>% count(.data[[col_name]]) %>% arrange(desc(n))
    factor(df[[col_name]], levels = counts[[col_name]])
  }
  metadata$race <- reorder_by_freq(metadata, "race")
  metadata$mod <- reorder_by_freq(metadata, "mod")
  metadata$ethnicity <- reorder_by_freq(metadata, "ethnicity")
  # Rename columns for the final table
  metadata_renamed <- metadata %>%
    rename(
      "Age" = age,
      "Race" = race,
      "Ethnicity" = ethnicity,
      "BMI" = bmi,
      "Diabetes" = diabetes,
      "Delivery mode" = mod,
      "Gestational age (weeks)" = gestational_age,
      "Term status (<37 weeks)" = term_status,
      "Mother ICU" = mom_icu,
      "Antibiotics during pregnancy" = antibiotics_during_pregnancy,
      "Antibiotics during c-section" = antibiotic_during_cs,
      "Newborn sex" = sex_of_baby,
      "Newborn weight" = baby_weight,
      "Newborn ICU" = baby_nicu    )
  
  # Create and save the table
  summary_table <- CreateTableOne(vars = variables, strata = "severe", data = metadata_renamed)
  table_df <- print(summary_table, quote = FALSE, noSpaces = TRUE, printToggle = FALSE, showAllLevels = FALSE, includeNA = TRUE)
  
  # --- Add asterisk to significant p-values ---
  # Check if the 'p' column exists
  if ("p" %in% colnames(table_df)) {
    # Convert p-value column to numeric, suppressing warnings for non-numeric entries like "<0.001"
    p_values_num <- suppressWarnings(as.numeric(table_df[, "p"]))
    
    # Identify rows where p < 0.05 or where the text contains "<"
    significant_rows <- (p_values_num < 0.050 & !is.na(p_values_num)) | grepl("<", table_df[, "p"])
    
    # Add an asterisk to those rows. Ensure we don't add to empty p-value cells.
    original_p <- table_df[significant_rows, "p"]
    table_df[significant_rows, "p"] <- paste0(original_p, "*")
  }
  
  # Save the modified data frame
  write.csv(table_df, file = file_name)
  
  cat("Successfully created and saved table with significance marks to:", file_name, "\n")
}


# 3. EXECUTION
# Define the variables for each table
vars_sev <- c("Age", "Race", "Ethnicity","BMI", "Diabetes", "Gestational age (weeks)", "Term status (<37 weeks)","Delivery mode", "Antibiotics during pregnancy", "Mother ICU", "Newborn sex", "Newborn weight", "Newborn ICU")



phy_subset_4 <- subset_samples(base_phy_object, body_site == "Vaginal" & severe_3cat != "control")
create_and_save_table(phy_subset_4, vars_sev, "results/tables/table1_mom_baby_severe_meconium_sequenced.csv")


#### CST table ####
a1 <- read.csv("/Users/danielavargasrobles/UMass Medical School Dropbox/Daniela Vargas Robles/Ana_projects/2_PREECLAMPSIA/met/met_clean_vag_tong.csv", header = TRUE, row.names = 1)
a2 <- a1 %>%
  filter(body_site=="Vaginal") %>%
  dplyr::select(CST5, CST, subCST, study_group)

# Define the variables to include in the table
vars <- colnames(a2)
vars <- vars[!vars %in% c("study_group")]


# Create the table
table <- CreateTableOne(vars = vars, strata = "study_group", data = a2)
print(table)

table1_df <- print(table, quote = FALSE, noSpaces = TRUE, printToggle = FALSE, 
                   showAllLevels = TRUE, includeNA = TRUE)

# Guardar la tabla como un archivo CSV
write.csv(table1_df, file = "results/tables/table1_cst_all.csv")


#### CST table by severity #####

a2 <- a1 %>%
  filter(severe_3cat!="control") %>% 
  dplyr::select(CST5, CST, subCST, study_group,severe)

# Define the variables to include in the table
vars <- colnames(a2)
vars <- vars[!vars %in% c("severe")]


# Create the table
table <- CreateTableOne(vars = vars, strata = "severe", data = a2)
print(table)

table1_df <- print(table, quote = FALSE, noSpaces = TRUE, printToggle = FALSE, 
                   showAllLevels = TRUE, includeNA = TRUE)

# Guardar la tabla como un archivo CSV
write.csv(table1_df, file = "results/tables/table1_cst_severe.csv")

a2 <- a1 %>%
  #filter(severe_3cat!="control") %>% 
  dplyr::select(CST5, CST, subCST, study_group,severe_3cat)


# Define the variables to include in the table
vars <- colnames(a2)
vars <- vars[!vars %in% c("severe_3cat")]


# Create the table
table <- CreateTableOne(vars = vars, strata = "severe_3cat", data = a2)
print(table)

table1_df <- print(table, quote = FALSE, noSpaces = TRUE, printToggle = FALSE, 
                   showAllLevels = TRUE, includeNA = TRUE)

# Guardar la tabla como un archivo CSV
write.csv(table1_df, file = "results/tables/table1_cst_severe_3cat.csv")


#### Stratifying by Ethnicity #####

setwd("~/Dropbox (UMass Medical School)/Ana_projects/2_PREECLAMPSIA/")
a1=read.csv("met/met_clean.csv")


a1=subset(a1,race == "White/Caucasian" & ethnicity == "Not Hispanic/Latina" | ethnicity == "Hispanic/Latina")
#phy2=subset_samples(phy2, study_group == "Control" )

# Calculate the frequency of each level
race_counts <- a1 %>% count(race)
mod_counts <- a1 %>% count(mod)
ethnicity_counts <- a1 %>% count(ethnicity)

# Reorder the levels of 'race' based on these frequencies
a1$race <- factor(a1$race, levels = race_counts %>% arrange(desc(n)) %>% pull(race))
a1$mod <- factor(a1$mod, levels = mod_counts %>% arrange(desc(n)) %>% pull(mod))
a1$ethnicity <- factor(a1$ethnicity, levels = ethnicity_counts %>% arrange(desc(n)) %>% pull(ethnicity))

a2 <- a1 %>%
  rename(
    study_group = "study_group",
    "Age" = age,
    "Race" = race,
    "Ethnicity" = ethnicity,
    "BMI" = bmi,
    "Severe preeclampsia" = severe,
    "Diabetes" = diabetes,
    "Delivery mode" = mod,
    "Antibiotics during pregnancy" = antibiotics_during_pregnancy,
    # "Antibiotics during labor" = antibiotic_during_labor,
    #"Antibiotics during c-section" = antibiotic_during_cs,
    "Mother ICU"= mom_icu,
    # "Baby sex" = sex_of_baby,
    # "Baby ICU" = baby_nicu,
    "Term status" = term_status
  )

# Define the variables
vars <- c("study_group","Age","Race","Ethnicity","BMI","Severe preeclampsia","Diabetes","Delivery mode",
          "Antibiotics during pregnancy","Antibiotics during labor","Antibiotics during c-section",
          "Mother ICU","Baby sex","Baby ICU","Term status")

# Create the table
table <- CreateTableOne(vars = vars, strata = "Ethnicity", data = a2)
#table <- CreateTableOne(vars = vars, strata = "study_group", data = a2)

table
#save table in P16S/results/table
table1_df <- print(table, quote = FALSE, noSpaces = TRUE, printToggle = FALSE, 
                   showAllLevels = TRUE, includeNA = TRUE)


# Guardar la tabla como un archivo CSV
write.csv(table1_df, file = "results/tables/table_all_var.csv")

#### Blood Pressure table ####

#table1: ppsystolic, ppdiastolic,
met <- read.csv("/Users/danielavargasrobles/UMass Medical School Dropbox/Daniela Vargas Robles/Ana_projects/2_PREECLAMPSIA/met/met_clean_vag_tong.csv", header = TRUE, row.names = 1)

# Crear una tabla de resumen para las variables de interés
vars <- c("ppsystolic", "ppdiastolic")
table1 <- CreateTableOne(vars = vars, strata = "severe_3cat", data = met_plot, factorVars = "study_group")
# Imprimir la tabla
print(table1, showAllLevels = TRUE)

#test normality
shapiro.test(met_plot$ppsystolic)
shapiro.test(met_plot$ppdiastolic)

#anove y pairwise
anova_ppsys <- aov(ppsystolic ~ severe_3cat, data = met_plot)
summary(anova_ppsys)
pairwise.t.test(met_plot$ppsystolic, met_plot$severe_3cat, p.adjust.method = "bonferroni")

#diastolic
anova_ppdia <- aov(ppdiastolic ~ severe_3cat, data = met_plot)
summary(anova_ppdia)
pairwise.t.test(met_plot$ppdiastolic, met_plot$severe_3cat, p.adjust.method = "bonferroni")

####Table1 BP####
#table1: ppsystolic, ppdiastolic,
library(tableone)
# Crear una tabla de resumen para las variables de interés
vars <- c("ppsystolic", "ppdiastolic")
table1 <- CreateTableOne(vars = vars, strata = "severe_3cat", data = met_plot, factorVars = "study_group")
# Imprimir la tabla
print(table1, showAllLevels = TRUE)
#number of records what are not NA in ppsystolic and ppdiastolic 
nrow(met_plot[!is.na(met_plot$ppsystolic) & !is.na(met_plot$ppdiastolic), ])
#number of records what are not NA in ppsystolic and ppdiastolic by severity_3cat
table(met_plot$severe_3cat, !is.na(met_plot$ppsystolic) & !is.na(met_plot$ppdiastolic))

#test normality
shapiro.test(met_plot$ppsystolic)
shapiro.test(met_plot$ppdiastolic)

#anove y pairwise
anova_ppsys <- aov(ppsystolic ~ severe_3cat, data = met_plot)
summary(anova_ppsys)
pairwise.t.test(met_plot$ppsystolic, met_plot$severe_3cat, p.adjust.method = "bonferroni")

#diastolic
anova_ppdia <- aov(ppdiastolic ~ severe_3cat, data = met_plot)
summary(anova_ppdia)
pairwise.t.test(met_plot$ppdiastolic, met_plot$severe_3cat, p.adjust.method = "bonferroni")



