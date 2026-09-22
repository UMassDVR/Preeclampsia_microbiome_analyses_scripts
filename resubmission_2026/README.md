# Resubmission 2026 analysis scripts

This directory contains the scripts used for the final resubmission analyses.
It supplements the original analysis notebooks in the parent directories,
which are retained unchanged as the original workflow record.

## Scripts

1. `01_postpartum_BP_alpha_diversity.R` — maternal oral, vaginal, and gut
   alpha diversity at delivery versus 6-week postpartum systolic and diastolic
   blood pressure. The primary models use HC3 robust standard errors;
   gestational-age and inverse-probability-of-follow-up sensitivity analyses,
   as well as racial/ethnic interaction and stratified analyses, are included.
   Benjamini–Hochberg adjusted p values are calculated within each body site,
   blood-pressure outcome, and model across Observed richness, Shannon, and
   Simpson diversity. Stratified simple-slope p values are additionally
   adjusted within body site, blood-pressure outcome, and racial/ethnic group.

2. `02_maternal_gut_alpha_beta_diversity.R` — maternal gut alpha and beta
   diversity analyses, racial/ethnic interaction and stratified analyses, and
   maternal gut–infant resemblance analyses.

3. `03_HUMAnN3_functional_unstratified.R` — maternal gut HUMAnN3 pathway and
   UniRef90 gene-family analyses using unstratified CPM-normalized profiles.

4. `04_HUMAnN3_functional_stratified.R` — targeted taxonomic-contributor
   analyses for PWY-5265 (peptidoglycan biosynthesis II). These analyses are
   exploratory because the parent pathway did not survive pathway-wide FDR
   correction.

5. `05_Figure4_gut_functional_panels.R` — source code for the exploratory
   maternal gut functional panels assembled in Figure 4.

## Running the scripts

The scripts expect the study data directories (`P16S`, `WGS`, and `met`) to
be located in the project root. They can be run from that root, or the root
can be supplied explicitly before execution:

```r
Sys.setenv(PREECLAMPSIA_PROJECT_DIR = "/path/to/2_PREECLAMPSIA")
source("resubmission_2026/01_postpartum_BP_alpha_diversity.R")
```

Each script writes derived files to a dated results directory and does not
overwrite its input data.
