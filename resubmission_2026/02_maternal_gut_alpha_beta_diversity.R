# Maternal stool alpha and beta diversity analyses####
#
# Set PREECLAMPSIA_PROJECT_DIR before running from a different location. When
# run from this repository, the default resolves to the parent project folder.
script_arg <- grep(
  "^--file=",
  commandArgs(trailingOnly = FALSE),
  value = TRUE
)
script_dir <- if (length(script_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1]), mustWork = FALSE))
} else {
  getwd()
}
default_project_dir <- if (dir.exists(file.path(getwd(), "WGS"))) {
  getwd()
} else {
  normalizePath(file.path(script_dir, "..", ".."), mustWork = FALSE)
}
PROJECT_DIR <- Sys.getenv(
  "PREECLAMPSIA_PROJECT_DIR",
  unset = default_project_dir
)

# PE vs Control:
#   Prespecified minimal = study_group + ethnicity + bmi + antibiotics_during_pregnancy
#   Prespecified full = minimal + gestational_age + gest_diabetes + meds_antihypertensives
#   Only the minimally adjusted model repeated among term deliveries (gestational_age >= 37)
#
# PE severity among women with PE:
#   Primary parsimonious = severe + ethnicity + bmi + antibiotics_during_pregnancy
#   GA-adjusted sensitivity = primary + gestational_age
#   Both repeated among term deliveries where applicable
#
# Alpha diversity reproduces the original maternal stool pipeline:
#   sha -> Box-Cox
#   obs -> untransformed
#   simp -> untransformed
#
# Beta diversity:
#   Bray-Curtis
#   Binary Jaccard
#   PERMANOVA with adonis2(..., by = "margin")
#   Beta dispersion with betadisper/permutest

library(phyloseq)
library(vegan)
library(MASS)
library(dplyr)
library(tidyr)
library(broom)
library(readr)
library(tibble)

# Load maternal stool WGS phyloseq object

phy <- readRDS(
  file.path(PROJECT_DIR, "WGS", "obj", "phy_mtph_gtdb.rds")
)

# Output directory

project_dir <- file.path(PROJECT_DIR, "WGS")

out_dir <- file.path(
  project_dir,
  "results",
  "reviewer_adjusted_models",
  as.character(Sys.Date())
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

cat(
  "\nResults will be saved in:\n",
  out_dir,
  "\n"
)

# Extract metadata

met <- data.frame(
  sample_data(phy)
)

# Check required variables

required_vars <- c(
  "study_group",
  "severe",
  "ethnicity",
  "bmi",
  "antibiotics_during_pregnancy",
  "gestational_age",
  "gest_diabetes",
  "meds_antihypertensives",
  "sha",
  "obs",
  "simp"
)

missing_vars <- setdiff(
  required_vars,
  colnames(met)
)

if (length(missing_vars) > 0) {
  stop(
    paste(
      "Missing variables:",
      paste(
        missing_vars,
        collapse = ", "
      )
    )
  )
}

# Inspect sample distributions

cat("\nStudy group:\n")
print(table(met$study_group, useNA = "ifany"))

cat("\nSeverity:\n")
print(table(met$severe, useNA = "ifany"))

cat("\nEthnicity:\n")
print(table(met$ethnicity, useNA = "ifany"))

cat("\nAntibiotics during pregnancy:\n")
print(table(met$antibiotics_during_pregnancy, useNA = "ifany"))

cat("\nGestational diabetes:\n")
print(table(met$gest_diabetes, useNA = "ifany"))

cat("\nAntihypertensive medication:\n")
print(table(met$meds_antihypertensives, useNA = "ifany"))

cat("\nGestational age summary:\n")
print(summary(met$gestational_age))

cat("\nAlpha-diversity summaries:\n")
print(summary(met[, c("sha", "obs", "simp")]))

# Factor coding

met$study_group <- relevel(
  factor(met$study_group),
  ref = "Control"
)

met$ethnicity <- factor(
  met$ethnicity
)

met$antibiotics_during_pregnancy <- factor(
  met$antibiotics_during_pregnancy
)

met$gest_diabetes <- factor(
  met$gest_diabetes
)

met$meds_antihypertensives <- factor(
  met$meds_antihypertensives
)

met$severe <- factor(
  met$severe,
  levels = c("No", "Yes")
)

# Define term status

met$term_restricted <- ifelse(
  !is.na(met$gestational_age) & met$gestational_age >= 37,
  "Term",
  ifelse(
    !is.na(met$gestational_age),
    "Preterm",
    NA
  )
)

met$term_restricted <- factor(
  met$term_restricted
)

cat("\nTerm status:\n")
print(table(met$term_restricted, useNA = "ifany"))

# Box-Cox transformation for Shannon
# This matches the original maternal stool alpha-diversity analysis

shannon_nonmissing <- met %>%
  filter(!is.na(sha))

if (any(shannon_nonmissing$sha <= 0)) {
  stop(
    "The variable sha contains zero or negative values. Box-Cox requires positive values."
  )
}

boxcox_shannon <- MASS::boxcox(
  sha ~ 1,
  data = shannon_nonmissing,
  plotit = FALSE
)

lambda_shannon <- boxcox_shannon$x[
  which.max(boxcox_shannon$y)
]

cat(
  "\nShannon Box-Cox lambda:",
  lambda_shannon,
  "\n"
)

if (abs(lambda_shannon) < 0.001) {
  met$shannon_boxcox <- log(
    met$sha
  )
} else {
  met$shannon_boxcox <- (
    met$sha^lambda_shannon - 1
  ) / lambda_shannon
}

# Keep Observed and Simpson on their original scales

met$observed_original <- met$obs
met$simpson_original <- met$simp

# Add standardized/transformed metadata back to phyloseq

sample_data(phy) <- sample_data(
  met
)

# Safe Shapiro-Wilk helper

safe_shapiro <- function(x) {
  x <- x[is.finite(x)]
  
  if (length(x) < 3 || length(x) > 5000) {
    return(NA_real_)
  }
  
  shapiro.test(x)$p.value
}

# Extract linear-model coefficients and 95% CI

extract_lm_results <- function(
    model,
    metric,
    comparison,
    analysis,
    model_type) {
  
  broom::tidy(
    model,
    conf.int = TRUE,
    conf.level = 0.95
  ) %>%
    mutate(
      body_site = "Maternal_stool",
      metric = metric,
      comparison = comparison,
      analysis = analysis,
      model_type = model_type,
      n = nobs(model),
      adjusted_R2 = summary(model)$adj.r.squared,
      AIC = AIC(model)
    ) %>%
    relocate(
      body_site,
      metric,
      comparison,
      analysis,
      model_type,
      n,
      adjusted_R2,
      AIC
    )
}

# PE vs Control alpha-diversity models

run_alpha_PE <- function(
    met,
    metric,
    outcome,
    term_only = FALSE) {
  
  if (term_only) {
    met_analysis <- met %>%
      filter(
        !is.na(gestational_age),
        gestational_age >= 37
      )
    
    analysis_name <- "Term_only"
  } else {
    met_analysis <- met
    analysis_name <- "All_deliveries"
  }
  
  # Prespecified minimally adjusted model
  # Complete cases are based only on variables in this model,
  # preserving the original minimal-model analysis population.
  
  vars_minimal <- c(
    outcome,
    "study_group",
    "ethnicity",
    "bmi",
    "antibiotics_during_pregnancy"
  )
  
  dat_minimal <- met_analysis %>%
    dplyr::select(all_of(vars_minimal)) %>%
    drop_na()
  
  dat_minimal$study_group <- relevel(
    factor(dat_minimal$study_group),
    ref = "Control"
  )
  
  dat_minimal$ethnicity <- droplevels(
    factor(dat_minimal$ethnicity)
  )
  
  dat_minimal$antibiotics_during_pregnancy <- droplevels(
    factor(dat_minimal$antibiotics_during_pregnancy)
  )
  
  formula_minimal <- as.formula(
    paste(
      outcome,
      "~ study_group + ethnicity + bmi +",
      "antibiotics_during_pregnancy"
    )
  )
  
  mod_minimal <- lm(
    formula_minimal,
    data = dat_minimal
  )
  
  # Prespecified fully adjusted model
  
  vars_full <- c(
    outcome,
    "study_group",
    "ethnicity",
    "bmi",
    "antibiotics_during_pregnancy",
    "gestational_age",
    "gest_diabetes",
    "meds_antihypertensives"
  )
  
  dat_full <- met_analysis %>%
    dplyr::select(all_of(vars_full)) %>%
    drop_na()
  
  dat_full$study_group <- relevel(
    factor(dat_full$study_group),
    ref = "Control"
  )
  
  dat_full$ethnicity <- droplevels(
    factor(dat_full$ethnicity)
  )
  
  dat_full$antibiotics_during_pregnancy <- droplevels(
    factor(dat_full$antibiotics_during_pregnancy)
  )
  
  dat_full$gest_diabetes <- droplevels(
    factor(dat_full$gest_diabetes)
  )
  
  dat_full$meds_antihypertensives <- droplevels(
    factor(dat_full$meds_antihypertensives)
  )
  
  formula_full <- as.formula(
    paste(
      outcome,
      "~ study_group + ethnicity + bmi +",
      "antibiotics_during_pregnancy +",
      "gestational_age + gest_diabetes +",
      "meds_antihypertensives"
    )
  )
  
  mod_full <- lm(
    formula_full,
    data = dat_full
  )
  
  cat(
    "\nPE vs Control alpha",
    "\nMetric:", metric,
    "\nAnalysis:", analysis_name,
    "\nMinimal N:", nrow(dat_minimal),
    "\nFull N:", nrow(dat_full),
    "\n"
  )
  
  cat("\nMinimal model study-group counts:\n")
  print(table(dat_minimal$study_group))
  
  cat("\nFull model study-group counts:\n")
  print(table(dat_full$study_group))
  
  result_minimal <- extract_lm_results(
    mod_minimal,
    metric,
    "PE_vs_Control",
    analysis_name,
    "Prespecified_minimal"
  )
  
  result_full <- extract_lm_results(
    mod_full,
    metric,
    "PE_vs_Control",
    analysis_name,
    "Prespecified_full"
  )
  
  diagnostics <- tibble(
    body_site = "Maternal_stool",
    metric = metric,
    comparison = "PE_vs_Control",
    analysis = analysis_name,
    model_type = c(
      "Prespecified_minimal",
      "Prespecified_full"
    ),
    n = c(
      nrow(dat_minimal),
      nrow(dat_full)
    ),
    shapiro_p = c(
      safe_shapiro(residuals(mod_minimal)),
      safe_shapiro(residuals(mod_full))
    )
  )
  
  formulas <- tibble(
    body_site = "Maternal_stool",
    metric = metric,
    comparison = "PE_vs_Control",
    analysis = analysis_name,
    model_type = c(
      "Prespecified_minimal",
      "Prespecified_full"
    ),
    n = c(
      nrow(dat_minimal),
      nrow(dat_full)
    ),
    formula = c(
      paste(
        deparse(formula_minimal),
        collapse = ""
      ),
      paste(
        deparse(formula_full),
        collapse = ""
      )
    )
  )
  
  list(
    results = bind_rows(
      result_minimal,
      result_full
    ),
    diagnostics = diagnostics,
    formulas = formulas
  )
}

# Run PE vs Control alpha-diversity analyses

alpha_PE_objects <- list(
  run_alpha_PE(
    met = met,
    metric = "Observed",
    outcome = "observed_original",
    term_only = FALSE
  ),
  run_alpha_PE(
    met = met,
    metric = "Shannon",
    outcome = "shannon_boxcox",
    term_only = FALSE
  ),
  run_alpha_PE(
    met = met,
    metric = "Simpson",
    outcome = "simpson_original",
    term_only = FALSE
  ),
  run_alpha_PE(
    met = met,
    metric = "Observed",
    outcome = "observed_original",
    term_only = TRUE
  ),
  run_alpha_PE(
    met = met,
    metric = "Shannon",
    outcome = "shannon_boxcox",
    term_only = TRUE
  ),
  run_alpha_PE(
    met = met,
    metric = "Simpson",
    outcome = "simpson_original",
    term_only = TRUE
  )
)

alpha_PE_results <- bind_rows(
  lapply(
    alpha_PE_objects,
    `[[`,
    "results"
  )
)

alpha_PE_diagnostics <- bind_rows(
  lapply(
    alpha_PE_objects,
    `[[`,
    "diagnostics"
  )
)

alpha_PE_formulas <- bind_rows(
  lapply(
    alpha_PE_objects,
    `[[`,
    "formulas"
  )
)

alpha_PE_effect <- alpha_PE_results %>%
  filter(
    grepl(
      "^study_group",
      term
    )
  ) %>%
  dplyr::select(
    body_site,
    metric,
    comparison,
    analysis,
    model_type,
    n,
    term,
    estimate,
    conf.low,
    conf.high,
    p.value
  ) %>%
  arrange(
    metric,
    analysis,
    model_type
  )

# PE severity alpha-diversity models

run_alpha_severity <- function(
    met,
    metric,
    outcome,
    term_only = FALSE) {
  
  met_analysis <- met %>%
    filter(
      study_group == "PE",
      !is.na(severe)
    )
  
  if (term_only) {
    met_analysis <- met_analysis %>%
      filter(
        !is.na(gestational_age),
        gestational_age >= 37
      )
    
    analysis_name <- "Term_only"
  } else {
    analysis_name <- "All_deliveries"
  }
  
  # Use the same complete-case population for the primary
  # and GA-adjusted sensitivity models.
  
  vars_needed <- c(
    outcome,
    "severe",
    "ethnicity",
    "bmi",
    "antibiotics_during_pregnancy",
    "gestational_age"
  )
  
  dat <- met_analysis %>%
    dplyr::select(all_of(vars_needed)) %>%
    drop_na()
  
  dat$severe <- factor(
    dat$severe,
    levels = c("No", "Yes")
  )
  
  dat$ethnicity <- droplevels(
    factor(dat$ethnicity)
  )
  
  dat$antibiotics_during_pregnancy <- droplevels(
    factor(dat$antibiotics_during_pregnancy)
  )
  
  if (nlevels(droplevels(dat$severe)) < 2) {
    stop(
      paste(
        "Severity has fewer than two levels in",
        metric,
        analysis_name
      )
    )
  }
  
  formula_primary <- as.formula(
    paste(
      outcome,
      "~ severe + ethnicity + bmi +",
      "antibiotics_during_pregnancy"
    )
  )
  
  formula_ga <- as.formula(
    paste(
      outcome,
      "~ severe + ethnicity + bmi +",
      "antibiotics_during_pregnancy +",
      "gestational_age"
    )
  )
  
  mod_primary <- lm(
    formula_primary,
    data = dat
  )
  
  mod_ga <- lm(
    formula_ga,
    data = dat
  )
  
  cat(
    "\nPE severity alpha",
    "\nMetric:", metric,
    "\nAnalysis:", analysis_name,
    "\nN:", nrow(dat),
    "\n"
  )
  
  cat("\nSeverity counts:\n")
  print(table(dat$severe))
  
  result_primary <- extract_lm_results(
    mod_primary,
    metric,
    "PE_SF_vs_PE_NSF",
    analysis_name,
    "Primary_minimal"
  )
  
  result_ga <- extract_lm_results(
    mod_ga,
    metric,
    "PE_SF_vs_PE_NSF",
    analysis_name,
    "GA_adjusted_sensitivity"
  )
  
  diagnostics <- tibble(
    body_site = "Maternal_stool",
    metric = metric,
    comparison = "PE_SF_vs_PE_NSF",
    analysis = analysis_name,
    model_type = c(
      "Primary_minimal",
      "GA_adjusted_sensitivity"
    ),
    n = nrow(dat),
    shapiro_p = c(
      safe_shapiro(residuals(mod_primary)),
      safe_shapiro(residuals(mod_ga))
    )
  )
  
  formulas <- tibble(
    body_site = "Maternal_stool",
    metric = metric,
    comparison = "PE_SF_vs_PE_NSF",
    analysis = analysis_name,
    model_type = c(
      "Primary_minimal",
      "GA_adjusted_sensitivity"
    ),
    n = nrow(dat),
    formula = c(
      paste(
        deparse(formula_primary),
        collapse = ""
      ),
      paste(
        deparse(formula_ga),
        collapse = ""
      )
    )
  )
  
  list(
    results = bind_rows(
      result_primary,
      result_ga
    ),
    diagnostics = diagnostics,
    formulas = formulas
  )
}

# Run PE severity alpha-diversity analyses

alpha_severity_objects <- list(
  run_alpha_severity(
    met = met,
    metric = "Observed",
    outcome = "observed_original",
    term_only = FALSE
  ),
  run_alpha_severity(
    met = met,
    metric = "Shannon",
    outcome = "shannon_boxcox",
    term_only = FALSE
  ),
  run_alpha_severity(
    met = met,
    metric = "Simpson",
    outcome = "simpson_original",
    term_only = FALSE
  ),
  run_alpha_severity(
    met = met,
    metric = "Observed",
    outcome = "observed_original",
    term_only = TRUE
  ),
  run_alpha_severity(
    met = met,
    metric = "Shannon",
    outcome = "shannon_boxcox",
    term_only = TRUE
  ),
  run_alpha_severity(
    met = met,
    metric = "Simpson",
    outcome = "simpson_original",
    term_only = TRUE
  )
)

alpha_severity_results <- bind_rows(
  lapply(
    alpha_severity_objects,
    `[[`,
    "results"
  )
)

alpha_severity_diagnostics <- bind_rows(
  lapply(
    alpha_severity_objects,
    `[[`,
    "diagnostics"
  )
)

alpha_severity_formulas <- bind_rows(
  lapply(
    alpha_severity_objects,
    `[[`,
    "formulas"
  )
)

alpha_severity_effect <- alpha_severity_results %>%
  filter(
    grepl(
      "^severe",
      term
    )
  ) %>%
  dplyr::select(
    body_site,
    metric,
    comparison,
    analysis,
    model_type,
    n,
    term,
    estimate,
    conf.low,
    conf.high,
    p.value
  ) %>%
  arrange(
    metric,
    analysis,
    model_type
  )

# Save alpha-diversity results

write_csv(
  alpha_PE_effect,
  file.path(
    out_dir,
    "stool_alpha_PE_effect_models_95CI.csv"
  )
)

write_csv(
  alpha_PE_results,
  file.path(
    out_dir,
    "stool_alpha_PE_all_coefficients_95CI.csv"
  )
)

write_csv(
  alpha_PE_diagnostics,
  file.path(
    out_dir,
    "stool_alpha_PE_residual_diagnostics.csv"
  )
)

write_csv(
  alpha_PE_formulas,
  file.path(
    out_dir,
    "stool_alpha_PE_model_formulas.csv"
  )
)

write_csv(
  alpha_severity_effect,
  file.path(
    out_dir,
    "stool_alpha_severity_effect_95CI_parsimonious.csv"
  )
)

write_csv(
  alpha_severity_results,
  file.path(
    out_dir,
    "stool_alpha_severity_all_coefficients_95CI_parsimonious.csv"
  )
)

write_csv(
  alpha_severity_diagnostics,
  file.path(
    out_dir,
    "stool_alpha_severity_residual_diagnostics_parsimonious.csv"
  )
)

write_csv(
  alpha_severity_formulas,
  file.path(
    out_dir,
    "stool_alpha_severity_model_formulas_parsimonious.csv"
  )
)

# PE vs Control beta-diversity models

run_beta_PE <- function(
    phy,
    distance_method,
    term_only = FALSE,
    permutations = 9999,
    seed = 711) {
  
  met_all <- data.frame(
    sample_data(phy)
  )
  
  if (term_only) {
    term_ids <- rownames(met_all)[
      !is.na(met_all$gestational_age) &
        met_all$gestational_age >= 37
    ]
    
    phy_analysis <- prune_samples(
      term_ids,
      phy
    )
    
    analysis_name <- "Term_only"
  } else {
    phy_analysis <- phy
    analysis_name <- "All_deliveries"
  }
  
  phy_analysis <- prune_taxa(
    taxa_sums(phy_analysis) > 0,
    phy_analysis
  )
  
  met_analysis <- data.frame(
    sample_data(phy_analysis)
  )
  
  # Use the same complete-case population for minimal and full
  # beta-diversity models, matching the reviewer-adjusted
  # oral/vaginal beta-diversity analysis.
  
  vars_full <- c(
    "study_group",
    "ethnicity",
    "bmi",
    "antibiotics_during_pregnancy",
    "gestational_age",
    "gest_diabetes",
    "meds_antihypertensives"
  )
  
  complete_ids <- rownames(met_analysis)[
    complete.cases(
      met_analysis[
        ,
        vars_full,
        drop = FALSE
      ]
    )
  ]
  
  phy_cc <- prune_samples(
    complete_ids,
    phy_analysis
  )
  
  phy_cc <- prune_taxa(
    taxa_sums(phy_cc) > 0,
    phy_cc
  )
  
  met_cc <- data.frame(
    sample_data(phy_cc)
  )
  
  met_cc$study_group <- relevel(
    factor(met_cc$study_group),
    ref = "Control"
  )
  
  met_cc$ethnicity <- droplevels(
    factor(met_cc$ethnicity)
  )
  
  met_cc$antibiotics_during_pregnancy <- droplevels(
    factor(met_cc$antibiotics_during_pregnancy)
  )
  
  met_cc$gest_diabetes <- droplevels(
    factor(met_cc$gest_diabetes)
  )
  
  met_cc$meds_antihypertensives <- droplevels(
    factor(met_cc$meds_antihypertensives)
  )
  
  if (distance_method == "jaccard") {
    dist_obj <- phyloseq::distance(
      phy_cc,
      method = "jaccard",
      binary = TRUE
    )
  } else if (distance_method == "bray") {
    dist_obj <- phyloseq::distance(
      phy_cc,
      method = "bray"
    )
  } else {
    stop(
      "distance_method must be either 'bray' or 'jaccard'"
    )
  }
  
  formula_minimal <-
    dist_obj ~
    study_group +
    ethnicity +
    bmi +
    antibiotics_during_pregnancy
  
  formula_full <-
    dist_obj ~
    study_group +
    ethnicity +
    bmi +
    antibiotics_during_pregnancy +
    gestational_age +
    gest_diabetes +
    meds_antihypertensives
  
  set.seed(seed)
  
  permanova_minimal <- vegan::adonis2(
    formula_minimal,
    data = met_cc,
    permutations = permutations,
    by = "margin"
  )
  
  set.seed(seed)
  
  permanova_full <- vegan::adonis2(
    formula_full,
    data = met_cc,
    permutations = permutations,
    by = "margin"
  )
  
  minimal_table <- as.data.frame(
    permanova_minimal
  ) %>%
    rownames_to_column(
      "term"
    ) %>%
    mutate(
      body_site = "Maternal_stool",
      comparison = "PE_vs_Control",
      distance = distance_method,
      analysis = analysis_name,
      model_type = "Prespecified_minimal",
      n = nsamples(phy_cc)
    ) %>%
    relocate(
      body_site,
      comparison,
      distance,
      analysis,
      model_type,
      n,
      term
    )
  
  full_table <- as.data.frame(
    permanova_full
  ) %>%
    rownames_to_column(
      "term"
    ) %>%
    mutate(
      body_site = "Maternal_stool",
      comparison = "PE_vs_Control",
      distance = distance_method,
      analysis = analysis_name,
      model_type = "Prespecified_full",
      n = nsamples(phy_cc)
    ) %>%
    relocate(
      body_site,
      comparison,
      distance,
      analysis,
      model_type,
      n,
      term
    )
  
  dispersion_model <- vegan::betadisper(
    dist_obj,
    met_cc$study_group
  )
  
  set.seed(seed)
  
  dispersion_test <- vegan::permutest(
    dispersion_model,
    permutations = permutations
  )
  
  dispersion_table <- tibble(
    body_site = "Maternal_stool",
    comparison = "PE_vs_Control",
    distance = distance_method,
    analysis = analysis_name,
    n = nsamples(phy_cc),
    F = dispersion_test$tab[
      1,
      "F"
    ],
    p_value = dispersion_test$tab[
      1,
      "Pr(>F)"
    ]
  )
  
  dispersion_values <- tibble(
    sample_id = names(
      dispersion_model$distances
    ),
    group = as.character(
      met_cc[
        names(dispersion_model$distances),
        "study_group"
      ]
    ),
    distance_to_centroid = as.numeric(
      dispersion_model$distances
    ),
    body_site = "Maternal_stool",
    comparison = "PE_vs_Control",
    distance = distance_method,
    analysis = analysis_name
  )
  
  formulas <- tibble(
    body_site = "Maternal_stool",
    comparison = "PE_vs_Control",
    distance = distance_method,
    analysis = analysis_name,
    model_type = c(
      "Prespecified_minimal",
      "Prespecified_full"
    ),
    n = nsamples(phy_cc),
    formula = c(
      "study_group + ethnicity + bmi + antibiotics_during_pregnancy",
      "study_group + ethnicity + bmi + antibiotics_during_pregnancy + gestational_age + gest_diabetes + meds_antihypertensives"
    )
  )
  
  cat(
    "\nPE vs Control beta",
    "\nDistance:", distance_method,
    "\nAnalysis:", analysis_name,
    "\nN:", nsamples(phy_cc),
    "\n"
  )
  
  cat("\nStudy-group counts:\n")
  print(table(met_cc$study_group))
  
  list(
    permanova = bind_rows(
      minimal_table,
      full_table
    ),
    dispersion = dispersion_table,
    dispersion_values = dispersion_values,
    formulas = formulas
  )
}

# Run PE vs Control beta-diversity analyses

beta_PE_objects <- list(
  run_beta_PE(
    phy = phy,
    distance_method = "bray",
    term_only = FALSE
  ),
  run_beta_PE(
    phy = phy,
    distance_method = "jaccard",
    term_only = FALSE
  ),
  run_beta_PE(
    phy = phy,
    distance_method = "bray",
    term_only = TRUE
  ),
  run_beta_PE(
    phy = phy,
    distance_method = "jaccard",
    term_only = TRUE
  )
)

beta_PE_results <- bind_rows(
  lapply(
    beta_PE_objects,
    `[[`,
    "permanova"
  )
)

beta_PE_dispersion <- bind_rows(
  lapply(
    beta_PE_objects,
    `[[`,
    "dispersion"
  )
)

beta_PE_dispersion_values <- bind_rows(
  lapply(
    beta_PE_objects,
    `[[`,
    "dispersion_values"
  )
)

beta_PE_formulas <- bind_rows(
  lapply(
    beta_PE_objects,
    `[[`,
    "formulas"
  )
)

beta_PE_effect <- beta_PE_results %>%
  filter(
    term == "study_group"
  ) %>%
  dplyr::select(
    body_site,
    comparison,
    distance,
    analysis,
    model_type,
    n,
    R2,
    F,
    `Pr(>F)`
  ) %>%
  rename(
    p_value = `Pr(>F)`
  ) %>%
  arrange(
    distance,
    analysis,
    model_type
  )

# PE severity beta-diversity models

run_beta_severity <- function(
    phy,
    distance_method,
    term_only = FALSE,
    permutations = 9999,
    seed = 711) {
  
  met_all <- data.frame(
    sample_data(phy)
  )
  
  pe_ids <- rownames(met_all)[
    !is.na(met_all$study_group) &
      met_all$study_group == "PE" &
      !is.na(met_all$severe)
  ]
  
  phy_pe <- prune_samples(
    pe_ids,
    phy
  )
  
  phy_pe <- prune_taxa(
    taxa_sums(phy_pe) > 0,
    phy_pe
  )
  
  if (term_only) {
    met_pe <- data.frame(
      sample_data(phy_pe)
    )
    
    term_ids <- rownames(met_pe)[
      !is.na(met_pe$gestational_age) &
        met_pe$gestational_age >= 37
    ]
    
    phy_pe <- prune_samples(
      term_ids,
      phy_pe
    )
    
    phy_pe <- prune_taxa(
      taxa_sums(phy_pe) > 0,
      phy_pe
    )
    
    analysis_name <- "Term_only"
  } else {
    analysis_name <- "All_deliveries"
  }
  
  met_pe <- data.frame(
    sample_data(phy_pe)
  )
  
  # Same complete-case population for the primary
  # and GA-adjusted severity models.
  
  vars_needed <- c(
    "severe",
    "ethnicity",
    "bmi",
    "antibiotics_during_pregnancy",
    "gestational_age"
  )
  
  complete_ids <- rownames(met_pe)[
    complete.cases(
      met_pe[
        ,
        vars_needed,
        drop = FALSE
      ]
    )
  ]
  
  phy_cc <- prune_samples(
    complete_ids,
    phy_pe
  )
  
  phy_cc <- prune_taxa(
    taxa_sums(phy_cc) > 0,
    phy_cc
  )
  
  met_cc <- data.frame(
    sample_data(phy_cc)
  )
  
  met_cc$severe <- factor(
    met_cc$severe,
    levels = c("No", "Yes")
  )
  
  met_cc$ethnicity <- droplevels(
    factor(met_cc$ethnicity)
  )
  
  met_cc$antibiotics_during_pregnancy <- droplevels(
    factor(met_cc$antibiotics_during_pregnancy)
  )
  
  if (nlevels(droplevels(met_cc$severe)) < 2) {
    stop(
      paste(
        "Severity has fewer than two levels in",
        distance_method,
        analysis_name
      )
    )
  }
  
  if (distance_method == "jaccard") {
    dist_obj <- phyloseq::distance(
      phy_cc,
      method = "jaccard",
      binary = TRUE
    )
  } else if (distance_method == "bray") {
    dist_obj <- phyloseq::distance(
      phy_cc,
      method = "bray"
    )
  } else {
    stop(
      "distance_method must be either 'bray' or 'jaccard'"
    )
  }
  
  formula_primary <-
    dist_obj ~
    severe +
    ethnicity +
    bmi +
    antibiotics_during_pregnancy
  
  formula_ga <-
    dist_obj ~
    severe +
    ethnicity +
    bmi +
    antibiotics_during_pregnancy +
    gestational_age
  
  set.seed(seed)
  
  permanova_primary <- vegan::adonis2(
    formula_primary,
    data = met_cc,
    permutations = permutations,
    by = "margin"
  )
  
  set.seed(seed)
  
  permanova_ga <- vegan::adonis2(
    formula_ga,
    data = met_cc,
    permutations = permutations,
    by = "margin"
  )
  
  primary_table <- as.data.frame(
    permanova_primary
  ) %>%
    rownames_to_column(
      "term"
    ) %>%
    mutate(
      body_site = "Maternal_stool",
      comparison = "PE_SF_vs_PE_NSF",
      distance = distance_method,
      analysis = analysis_name,
      model_type = "Primary_minimal",
      n = nsamples(phy_cc)
    ) %>%
    relocate(
      body_site,
      comparison,
      distance,
      analysis,
      model_type,
      n,
      term
    )
  
  ga_table <- as.data.frame(
    permanova_ga
  ) %>%
    rownames_to_column(
      "term"
    ) %>%
    mutate(
      body_site = "Maternal_stool",
      comparison = "PE_SF_vs_PE_NSF",
      distance = distance_method,
      analysis = analysis_name,
      model_type = "GA_adjusted_sensitivity",
      n = nsamples(phy_cc)
    ) %>%
    relocate(
      body_site,
      comparison,
      distance,
      analysis,
      model_type,
      n,
      term
    )
  
  dispersion_model <- vegan::betadisper(
    dist_obj,
    met_cc$severe
  )
  
  set.seed(seed)
  
  dispersion_test <- vegan::permutest(
    dispersion_model,
    permutations = permutations
  )
  
  dispersion_table <- tibble(
    body_site = "Maternal_stool",
    comparison = "PE_SF_vs_PE_NSF",
    distance = distance_method,
    analysis = analysis_name,
    n = nsamples(phy_cc),
    F = dispersion_test$tab[
      1,
      "F"
    ],
    p_value = dispersion_test$tab[
      1,
      "Pr(>F)"
    ]
  )
  
  dispersion_values <- tibble(
    sample_id = names(
      dispersion_model$distances
    ),
    group = as.character(
      met_cc[
        names(dispersion_model$distances),
        "severe"
      ]
    ),
    distance_to_centroid = as.numeric(
      dispersion_model$distances
    ),
    body_site = "Maternal_stool",
    comparison = "PE_SF_vs_PE_NSF",
    distance = distance_method,
    analysis = analysis_name
  )
  
  formulas <- tibble(
    body_site = "Maternal_stool",
    comparison = "PE_SF_vs_PE_NSF",
    distance = distance_method,
    analysis = analysis_name,
    model_type = c(
      "Primary_minimal",
      "GA_adjusted_sensitivity"
    ),
    n = nsamples(phy_cc),
    formula = c(
      "severe + ethnicity + bmi + antibiotics_during_pregnancy",
      "severe + ethnicity + bmi + antibiotics_during_pregnancy + gestational_age"
    )
  )
  
  cat(
    "\nPE severity beta",
    "\nDistance:", distance_method,
    "\nAnalysis:", analysis_name,
    "\nN:", nsamples(phy_cc),
    "\n"
  )
  
  cat("\nSeverity counts:\n")
  print(table(met_cc$severe))
  
  list(
    permanova = bind_rows(
      primary_table,
      ga_table
    ),
    dispersion = dispersion_table,
    dispersion_values = dispersion_values,
    formulas = formulas
  )
}

# Run PE severity beta-diversity analyses

beta_severity_objects <- list(
  run_beta_severity(
    phy = phy,
    distance_method = "bray",
    term_only = FALSE
  ),
  run_beta_severity(
    phy = phy,
    distance_method = "jaccard",
    term_only = FALSE
  ),
  run_beta_severity(
    phy = phy,
    distance_method = "bray",
    term_only = TRUE
  ),
  run_beta_severity(
    phy = phy,
    distance_method = "jaccard",
    term_only = TRUE
  )
)

beta_severity_results <- bind_rows(
  lapply(
    beta_severity_objects,
    `[[`,
    "permanova"
  )
)

beta_severity_dispersion <- bind_rows(
  lapply(
    beta_severity_objects,
    `[[`,
    "dispersion"
  )
)

beta_severity_dispersion_values <- bind_rows(
  lapply(
    beta_severity_objects,
    `[[`,
    "dispersion_values"
  )
)

beta_severity_formulas <- bind_rows(
  lapply(
    beta_severity_objects,
    `[[`,
    "formulas"
  )
)

beta_severity_effect <- beta_severity_results %>%
  filter(
    term == "severe"
  ) %>%
  dplyr::select(
    body_site,
    comparison,
    distance,
    analysis,
    model_type,
    n,
    R2,
    F,
    `Pr(>F)`
  ) %>%
  rename(
    p_value = `Pr(>F)`
  ) %>%
  arrange(
    distance,
    analysis,
    model_type
  )

# Save PE vs Control beta-diversity results

write_csv(
  beta_PE_effect,
  file.path(
    out_dir,
    "stool_beta_PE_effect_R2.csv"
  )
)

write_csv(
  beta_PE_results,
  file.path(
    out_dir,
    "stool_beta_PE_PERMANOVA_all_models.csv"
  )
)

write_csv(
  beta_PE_dispersion,
  file.path(
    out_dir,
    "stool_beta_PE_dispersion.csv"
  )
)

write_csv(
  beta_PE_dispersion_values,
  file.path(
    out_dir,
    "stool_beta_PE_dispersion_values.csv"
  )
)

write_csv(
  beta_PE_formulas,
  file.path(
    out_dir,
    "stool_beta_PE_model_formulas.csv"
  )
)

# Save PE severity beta-diversity results

write_csv(
  beta_severity_effect,
  file.path(
    out_dir,
    "stool_beta_severity_effect_R2_parsimonious.csv"
  )
)

write_csv(
  beta_severity_results,
  file.path(
    out_dir,
    "stool_beta_severity_PERMANOVA_all_models_parsimonious.csv"
  )
)

write_csv(
  beta_severity_dispersion,
  file.path(
    out_dir,
    "stool_beta_severity_dispersion_parsimonious.csv"
  )
)

write_csv(
  beta_severity_dispersion_values,
  file.path(
    out_dir,
    "stool_beta_severity_dispersion_values_parsimonious.csv"
  )
)

write_csv(
  beta_severity_formulas,
  file.path(
    out_dir,
    "stool_beta_severity_model_formulas_parsimonious.csv"
  )
)

# Print reviewer-focused results

cat("\nMaternal stool PE vs Control alpha results\n")
print(alpha_PE_effect)

cat("\nMaternal stool PE vs Control alpha residual diagnostics\n")
print(alpha_PE_diagnostics)

cat("\nMaternal stool PE vs Control beta results\n")
print(beta_PE_effect)

cat("\nMaternal stool PE vs Control beta dispersion\n")
print(beta_PE_dispersion)

cat("\nMaternal stool PE severity alpha results\n")
print(alpha_severity_effect)

cat("\nMaternal stool PE severity alpha residual diagnostics\n")
print(alpha_severity_diagnostics)

cat("\nMaternal stool PE severity beta results\n")
print(beta_severity_effect)

cat("\nMaternal stool PE severity beta dispersion\n")
print(beta_severity_dispersion)

# List generated files

cat("\nFiles generated:\n")

print(
  list.files(
    out_dir,
    pattern = "^stool_",
    full.names = TRUE
  )
)

####simpson LRM ####
# Maternal stool Simpson diversity: robust linear models
# PE vs Control:
#   Prespecified minimal = study_group + ethnicity + bmi + antibiotics_during_pregnancy
#   Prespecified full = minimal + gestational_age + gest_diabetes + meds_antihypertensives
#   Both repeated among term deliveries (gestational_age >= 37)
#
# PE severity among women with PE:
#   Primary parsimonious = severe + ethnicity + bmi + antibiotics_during_pregnancy
#   GA-adjusted sensitivity = primary + gestational_age
#   Both repeated among term deliveries
#
# Simpson is analyzed using MASS::rlm because residuals from the corresponding
# ordinary linear models were non-normal.

library(phyloseq)
library(MASS)
library(dplyr)
library(tidyr)
library(readr)
library(tibble)

# Load maternal stool phyloseq object

phy <- readRDS(
  file.path(PROJECT_DIR, "WGS", "obj", "phy_mtph_gtdb.rds")
)

met <- data.frame(
  sample_data(phy)
)

# Output directory

project_dir <- file.path(PROJECT_DIR, "WGS")

out_dir <- file.path(
  project_dir,
  "results",
  "reviewer_adjusted_models",
  as.character(Sys.Date())
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

cat(
  "\nResults will be saved in:\n",
  out_dir,
  "\n"
)

# Check required variables

required_vars <- c(
  "simp",
  "study_group",
  "severe",
  "ethnicity",
  "bmi",
  "antibiotics_during_pregnancy",
  "gestational_age",
  "gest_diabetes",
  "meds_antihypertensives"
)

missing_vars <- setdiff(
  required_vars,
  colnames(met)
)

if (length(missing_vars) > 0) {
  stop(
    paste(
      "Missing variables:",
      paste(
        missing_vars,
        collapse = ", "
      )
    )
  )
}

# Inspect variables

cat("\nStudy group:\n")
print(
  table(
    met$study_group,
    useNA = "ifany"
  )
)

cat("\nSeverity:\n")
print(
  table(
    met$severe,
    useNA = "ifany"
  )
)

cat("\nEthnicity:\n")
print(
  table(
    met$ethnicity,
    useNA = "ifany"
  )
)

cat("\nAntibiotics during pregnancy:\n")
print(
  table(
    met$antibiotics_during_pregnancy,
    useNA = "ifany"
  )
)

cat("\nGestational diabetes:\n")
print(
  table(
    met$gest_diabetes,
    useNA = "ifany"
  )
)

cat("\nAntihypertensive medication:\n")
print(
  table(
    met$meds_antihypertensives,
    useNA = "ifany"
  )
)

cat("\nGestational age:\n")
print(
  summary(
    met$gestational_age
  )
)

cat("\nSimpson diversity:\n")
print(
  summary(
    met$simp
  )
)

# Factor coding

met$study_group <- relevel(
  factor(
    met$study_group
  ),
  ref = "Control"
)

met$severe <- factor(
  met$severe,
  levels = c(
    "No",
    "Yes"
  )
)

met$ethnicity <- factor(
  met$ethnicity
)

met$antibiotics_during_pregnancy <- factor(
  met$antibiotics_during_pregnancy
)

met$gest_diabetes <- factor(
  met$gest_diabetes
)

met$meds_antihypertensives <- factor(
  met$meds_antihypertensives
)

# Helper to extract RLM coefficients, approximate p-values and 95% CI

extract_rlm_results <- function(
    model,
    comparison,
    analysis,
    model_type) {
  
  coef_table <- as.data.frame(
    summary(model)$coefficients
  )
  
  coef_table$term <- rownames(
    coef_table
  )
  
  colnames(coef_table)[1:3] <- c(
    "estimate",
    "std.error",
    "statistic"
  )
  
  coef_table %>%
    mutate(
      conf.low = estimate - 1.96 * std.error,
      conf.high = estimate + 1.96 * std.error,
      p.value = 2 * pnorm(
        abs(statistic),
        lower.tail = FALSE
      ),
      body_site = "Maternal_stool",
      metric = "Simpson",
      comparison = comparison,
      analysis = analysis,
      model_type = model_type,
      n = length(
        model$residuals
      )
    ) %>%
    dplyr::select(
      body_site,
      metric,
      comparison,
      analysis,
      model_type,
      n,
      term,
      estimate,
      std.error,
      statistic,
      conf.low,
      conf.high,
      p.value
    )
}

# PE vs Control Simpson RLM

run_simpson_PE_rlm <- function(
    met,
    term_only = FALSE) {
  
  if (term_only) {
    
    met_analysis <- met %>%
      filter(
        !is.na(
          gestational_age
        ),
        gestational_age >= 37
      )
    
    analysis_name <- "Term_only"
    
  } else {
    
    met_analysis <- met
    analysis_name <- "All_deliveries"
  }
  
  # Prespecified minimally adjusted model
  # Complete cases are based only on variables used in the minimal model.
  
  vars_minimal <- c(
    "simp",
    "study_group",
    "ethnicity",
    "bmi",
    "antibiotics_during_pregnancy"
  )
  
  dat_minimal <- met_analysis %>%
    dplyr::select(
      all_of(
        vars_minimal
      )
    ) %>%
    drop_na()
  
  dat_minimal$study_group <- relevel(
    factor(
      dat_minimal$study_group
    ),
    ref = "Control"
  )
  
  dat_minimal$ethnicity <- droplevels(
    factor(
      dat_minimal$ethnicity
    )
  )
  
  dat_minimal$antibiotics_during_pregnancy <- droplevels(
    factor(
      dat_minimal$antibiotics_during_pregnancy
    )
  )
  
  formula_minimal <-
    simp ~
    study_group +
    ethnicity +
    bmi +
    antibiotics_during_pregnancy
  
  mod_minimal <- MASS::rlm(
    formula_minimal,
    data = dat_minimal,
    maxit = 100
  )
  
  # Prespecified fully adjusted model
  
  vars_full <- c(
    "simp",
    "study_group",
    "ethnicity",
    "bmi",
    "antibiotics_during_pregnancy",
    "gestational_age",
    "gest_diabetes",
    "meds_antihypertensives"
  )
  
  dat_full <- met_analysis %>%
    dplyr::select(
      all_of(
        vars_full
      )
    ) %>%
    drop_na()
  
  dat_full$study_group <- relevel(
    factor(
      dat_full$study_group
    ),
    ref = "Control"
  )
  
  dat_full$ethnicity <- droplevels(
    factor(
      dat_full$ethnicity
    )
  )
  
  dat_full$antibiotics_during_pregnancy <- droplevels(
    factor(
      dat_full$antibiotics_during_pregnancy
    )
  )
  
  dat_full$gest_diabetes <- droplevels(
    factor(
      dat_full$gest_diabetes
    )
  )
  
  dat_full$meds_antihypertensives <- droplevels(
    factor(
      dat_full$meds_antihypertensives
    )
  )
  
  formula_full <-
    simp ~
    study_group +
    ethnicity +
    bmi +
    antibiotics_during_pregnancy +
    gestational_age +
    gest_diabetes +
    meds_antihypertensives
  
  mod_full <- MASS::rlm(
    formula_full,
    data = dat_full,
    maxit = 100
  )
  
  cat(
    "\nPE vs Control - Simpson RLM",
    "\nAnalysis:", analysis_name,
    "\nMinimal N:", nrow(
      dat_minimal
    ),
    "\nFull N:", nrow(
      dat_full
    ),
    "\n"
  )
  
  cat("\nMinimal model group counts:\n")
  print(
    table(
      dat_minimal$study_group
    )
  )
  
  cat("\nFull model group counts:\n")
  print(
    table(
      dat_full$study_group
    )
  )
  
  result_minimal <- extract_rlm_results(
    mod_minimal,
    comparison = "PE_vs_Control",
    analysis = analysis_name,
    model_type = "Prespecified_minimal_RLM"
  )
  
  result_full <- extract_rlm_results(
    mod_full,
    comparison = "PE_vs_Control",
    analysis = analysis_name,
    model_type = "Prespecified_full_RLM"
  )
  
  formulas <- tibble(
    body_site = "Maternal_stool",
    metric = "Simpson",
    comparison = "PE_vs_Control",
    analysis = analysis_name,
    model_type = c(
      "Prespecified_minimal_RLM",
      "Prespecified_full_RLM"
    ),
    n = c(
      nrow(
        dat_minimal
      ),
      nrow(
        dat_full
      )
    ),
    formula = c(
      paste(
        deparse(
          formula_minimal
        ),
        collapse = ""
      ),
      paste(
        deparse(
          formula_full
        ),
        collapse = ""
      )
    )
  )
  
  list(
    results = bind_rows(
      result_minimal,
      result_full
    ),
    formulas = formulas
  )
}

# Run PE vs Control RLM analyses

simpson_PE_all <- run_simpson_PE_rlm(
  met = met,
  term_only = FALSE
)

simpson_PE_term <- run_simpson_PE_rlm(
  met = met,
  term_only = TRUE
)

simpson_PE_results <- bind_rows(
  simpson_PE_all$results,
  simpson_PE_term$results
)

simpson_PE_formulas <- bind_rows(
  simpson_PE_all$formulas,
  simpson_PE_term$formulas
)

simpson_PE_effect <- simpson_PE_results %>%
  filter(
    grepl(
      "^study_group",
      term
    )
  ) %>%
  arrange(
    analysis,
    model_type
  )

# PE severity Simpson RLM

run_simpson_severity_rlm <- function(
    met,
    term_only = FALSE) {
  
  met_analysis <- met %>%
    filter(
      study_group == "PE",
      !is.na(
        severe
      )
    )
  
  if (term_only) {
    
    met_analysis <- met_analysis %>%
      filter(
        !is.na(
          gestational_age
        ),
        gestational_age >= 37
      )
    
    analysis_name <- "Term_only"
    
  } else {
    
    analysis_name <- "All_deliveries"
  }
  
  # Use the same complete-case population for the primary
  # and GA-adjusted severity models.
  
  vars_needed <- c(
    "simp",
    "severe",
    "ethnicity",
    "bmi",
    "antibiotics_during_pregnancy",
    "gestational_age"
  )
  
  dat <- met_analysis %>%
    dplyr::select(
      all_of(
        vars_needed
      )
    ) %>%
    drop_na()
  
  dat$severe <- factor(
    dat$severe,
    levels = c(
      "No",
      "Yes"
    )
  )
  
  dat$ethnicity <- droplevels(
    factor(
      dat$ethnicity
    )
  )
  
  dat$antibiotics_during_pregnancy <- droplevels(
    factor(
      dat$antibiotics_during_pregnancy
    )
  )
  
  if (
    nlevels(
      droplevels(
        dat$severe
      )
    ) < 2
  ) {
    stop(
      paste(
        "Severity has fewer than two levels in",
        analysis_name
      )
    )
  }
  
  formula_primary <-
    simp ~
    severe +
    ethnicity +
    bmi +
    antibiotics_during_pregnancy
  
  formula_ga <-
    simp ~
    severe +
    ethnicity +
    bmi +
    antibiotics_during_pregnancy +
    gestational_age
  
  mod_primary <- MASS::rlm(
    formula_primary,
    data = dat,
    maxit = 100
  )
  
  mod_ga <- MASS::rlm(
    formula_ga,
    data = dat,
    maxit = 100
  )
  
  cat(
    "\nPE severity - Simpson RLM",
    "\nAnalysis:", analysis_name,
    "\nN:", nrow(
      dat
    ),
    "\n"
  )
  
  cat("\nSeverity counts:\n")
  print(
    table(
      dat$severe
    )
  )
  
  result_primary <- extract_rlm_results(
    mod_primary,
    comparison = "PE_SF_vs_PE_NSF",
    analysis = analysis_name,
    model_type = "Primary_minimal_RLM"
  )
  
  result_ga <- extract_rlm_results(
    mod_ga,
    comparison = "PE_SF_vs_PE_NSF",
    analysis = analysis_name,
    model_type = "GA_adjusted_sensitivity_RLM"
  )
  
  formulas <- tibble(
    body_site = "Maternal_stool",
    metric = "Simpson",
    comparison = "PE_SF_vs_PE_NSF",
    analysis = analysis_name,
    model_type = c(
      "Primary_minimal_RLM",
      "GA_adjusted_sensitivity_RLM"
    ),
    n = nrow(
      dat
    ),
    formula = c(
      paste(
        deparse(
          formula_primary
        ),
        collapse = ""
      ),
      paste(
        deparse(
          formula_ga
        ),
        collapse = ""
      )
    )
  )
  
  list(
    results = bind_rows(
      result_primary,
      result_ga
    ),
    formulas = formulas
  )
}

# Run PE severity RLM analyses

simpson_severity_all <- run_simpson_severity_rlm(
  met = met,
  term_only = FALSE
)

simpson_severity_term <- run_simpson_severity_rlm(
  met = met,
  term_only = TRUE
)

simpson_severity_results <- bind_rows(
  simpson_severity_all$results,
  simpson_severity_term$results
)

simpson_severity_formulas <- bind_rows(
  simpson_severity_all$formulas,
  simpson_severity_term$formulas
)

simpson_severity_effect <- simpson_severity_results %>%
  filter(
    grepl(
      "^severe",
      term
    )
  ) %>%
  arrange(
    analysis,
    model_type
  )

# Save PE vs Control results

write_csv(
  simpson_PE_effect,
  file.path(
    out_dir,
    "stool_alpha_PE_Simpson_RLM_effect_95CI.csv"
  )
)

write_csv(
  simpson_PE_results,
  file.path(
    out_dir,
    "stool_alpha_PE_Simpson_RLM_all_coefficients_95CI.csv"
  )
)

write_csv(
  simpson_PE_formulas,
  file.path(
    out_dir,
    "stool_alpha_PE_Simpson_RLM_formulas.csv"
  )
)

# Save PE severity results

write_csv(
  simpson_severity_effect,
  file.path(
    out_dir,
    "stool_alpha_severity_Simpson_RLM_effect_95CI.csv"
  )
)

write_csv(
  simpson_severity_results,
  file.path(
    out_dir,
    "stool_alpha_severity_Simpson_RLM_all_coefficients_95CI.csv"
  )
)

write_csv(
  simpson_severity_formulas,
  file.path(
    out_dir,
    "stool_alpha_severity_Simpson_RLM_formulas.csv"
  )
)

# Print main reviewer-focused results

cat(
  "\nSimpson RLM - PE vs Control\n"
)

print(
  simpson_PE_effect
)

cat(
  "\nSimpson RLM - PE severity\n"
)

print(
  simpson_severity_effect
)

cat(
  "\nFiles generated:\n"
)

print(
  list.files(
    out_dir,
    pattern = "Simpson_RLM",
    full.names = TRUE
  )
)


# Ethnicity-stratified maternal stool alpha and beta diversity
# Exploratory analyses matching the oral/vaginal strategy:
#   Hispanic and White women analyzed separately
#   PE vs Control:
#     Primary = study_group + bmi + antibiotics_during_pregnancy
#     GA-adjusted sensitivity = primary + gestational_age
#   PE severity among women with PE:
#     Primary = severe + bmi + antibiotics_during_pregnancy
#     GA-adjusted sensitivity = primary + gestational_age
#   No additional restriction to term deliveries
#   Alpha:
#     Observed = LM, original scale
#     Shannon = LM, Box-Cox transformed as in the overall stool analysis
#     Simpson = RLM because the corresponding ordinary LM residuals were non-normal
#   Beta:
#     Bray-Curtis and binary Jaccard
#     PERMANOVA with adonis2(..., by = "margin")
#     Beta dispersion with betadisper/permutest

library(phyloseq)
library(vegan)
library(MASS)
library(dplyr)
library(tidyr)
library(broom)
library(readr)
library(tibble)

# Re-extract metadata from the maternal stool phyloseq object

met_strat <- data.frame(
  sample_data(phy)
)

required_strat_vars <- c(
  "study_group",
  "severe",
  "ethnicity",
  "race",
  "bmi",
  "antibiotics_during_pregnancy",
  "gestational_age",
  "sha",
  "obs",
  "simp"
)

missing_strat_vars <- setdiff(
  required_strat_vars,
  colnames(met_strat)
)

if (length(missing_strat_vars) > 0) {
  stop(
    paste(
      "Missing variables required for ethnicity-stratified analyses:",
      paste(
        missing_strat_vars,
        collapse = ", "
      )
    )
  )
}

# Define ethnicity strata exactly as used for oral/vaginal analyses

met_strat$ethnicity_analysis <- case_when(
  met_strat$ethnicity %in% c(
    "Hispanic/Latina",
    "Hispanic or Latina",
    "Hispanic"
  ) ~ "Hispanic",
  met_strat$ethnicity %in% c(
    "Not Hispanic/Latina",
    "Not Hispanic or Latina",
    "Not Hispanic"
  ) &
    met_strat$race %in% c(
      "White/Caucasian",
      "White"
    ) ~ "White",
  TRUE ~ "Other"
)

met_strat$ethnicity_analysis <- factor(
  met_strat$ethnicity_analysis,
  levels = c(
    "Hispanic",
    "White",
    "Other"
  )
)

cat(
  "\nEthnicity-analysis strata:\n"
)

print(
  table(
    met_strat$ethnicity_analysis,
    useNA = "ifany"
  )
)

cat(
  "\nStudy-group counts by ethnicity stratum:\n"
)

print(
  table(
    met_strat$ethnicity_analysis,
    met_strat$study_group,
    useNA = "ifany"
  )
)

cat(
  "\nSeverity counts among PE by ethnicity stratum:\n"
)

print(
  table(
    met_strat$ethnicity_analysis[
      met_strat$study_group == "PE"
    ],
    met_strat$severe[
      met_strat$study_group == "PE"
    ],
    useNA = "ifany"
  )
)

# Factor coding

met_strat$study_group <- relevel(
  factor(
    met_strat$study_group
  ),
  ref = "Control"
)

met_strat$severe <- factor(
  met_strat$severe,
  levels = c(
    "No",
    "Yes"
  )
)

met_strat$antibiotics_during_pregnancy <- factor(
  met_strat$antibiotics_during_pregnancy
)

# Recompute the Shannon Box-Cox transformation using the complete maternal stool dataset

shannon_nonmissing_strat <- met_strat %>%
  filter(
    !is.na(sha)
  )

if (any(shannon_nonmissing_strat$sha <= 0)) {
  stop(
    "The variable sha contains zero or negative values. Box-Cox requires positive values."
  )
}

boxcox_shannon_strat <- MASS::boxcox(
  sha ~ 1,
  data = shannon_nonmissing_strat,
  plotit = FALSE
)

lambda_shannon_strat <- boxcox_shannon_strat$x[
  which.max(
    boxcox_shannon_strat$y
  )
]

cat(
  "\nShannon Box-Cox lambda used for ethnicity-stratified analyses:",
  lambda_shannon_strat,
  "\n"
)

if (abs(lambda_shannon_strat) < 0.001) {
  met_strat$shannon_boxcox <- log(
    met_strat$sha
  )
} else {
  met_strat$shannon_boxcox <- (
    met_strat$sha^lambda_shannon_strat - 1
  ) / lambda_shannon_strat
}

met_strat$observed_original <- met_strat$obs
met_strat$simpson_original <- met_strat$simp

# Add the analysis variables back to phyloseq for beta-diversity subsetting

sample_data(phy) <- sample_data(
  met_strat
)

# Helper for user-friendly model-term labels

pretty_stool_term <- function(x) {
  case_when(
    x == "study_groupPE" ~ "Preeclampsia",
    x == "study_group" ~ "Preeclampsia",
    x == "severeYes" ~ "Preeclampsia with severe features",
    x == "severe" ~ "Preeclampsia with severe features",
    grepl(
      "^antibiotics_during_pregnancy",
      x
    ) ~ "Antibiotic use during pregnancy: Yes",
    x == "bmi" ~ "BMI",
    x == "gestational_age" ~ "Gestational age at delivery",
    TRUE ~ x
  )
}

# Helper to extract LM results

extract_strat_lm <- function(
    model_object,
    ethnicity_stratum,
    metric,
    comparison,
    model_type) {
  
  model_n <- nobs(model_object)
  
  broom::tidy(
    model_object,
    conf.int = TRUE,
    conf.level = 0.95
  ) %>%
    mutate(
      body_site = "Maternal_stool",
      ethnicity_stratum = ethnicity_stratum,
      analysis_method = "Alpha diversity",
      analysis_focus = comparison,
      population = "All_deliveries",
      metric = metric,
      model = model_type,
      statistical_method = "LM",
      n = model_n,
      model_term = pretty_stool_term(term)
    ) %>%
    dplyr::select(
      body_site,
      ethnicity_stratum,
      analysis_method,
      analysis_focus,
      population,
      metric,
      model,
      statistical_method,
      n,
      term,
      model_term,
      estimate,
      conf.low,
      conf.high,
      p.value
    )
}

# Helper to extract RLM results

extract_strat_rlm <- function(
    model_object,
    ethnicity_stratum,
    metric,
    comparison,
    model_type) {
  
  model_n <- length(model_object$residuals)
  
  coef_table <- as.data.frame(
    summary(model_object)$coefficients
  )
  
  coef_table$term <- rownames(
    coef_table
  )
  
  colnames(coef_table)[1:3] <- c(
    "estimate",
    "std.error",
    "statistic"
  )
  
  coef_table %>%
    mutate(
      conf.low = estimate - 1.96 * std.error,
      conf.high = estimate + 1.96 * std.error,
      p.value = 2 * pnorm(
        abs(statistic),
        lower.tail = FALSE
      ),
      body_site = "Maternal_stool",
      ethnicity_stratum = ethnicity_stratum,
      analysis_method = "Alpha diversity",
      analysis_focus = comparison,
      population = "All_deliveries",
      metric = metric,
      model = model_type,
      statistical_method = "RLM",
      n = model_n,
      model_term = pretty_stool_term(term)
    ) %>%
    dplyr::select(
      body_site,
      ethnicity_stratum,
      analysis_method,
      analysis_focus,
      population,
      metric,
      model,
      statistical_method,
      n,
      term,
      model_term,
      estimate,
      conf.low,
      conf.high,
      p.value
    )
}

# Alpha-diversity function for PE vs Control within ethnicity

run_alpha_ethnicity_status <- function(
    met,
    ethnicity_stratum,
    metric,
    outcome,
    method = c(
      "LM",
      "RLM"
    )) {
  
  method <- match.arg(
    method
  )
  
  dat0 <- met %>%
    filter(
      ethnicity_analysis == ethnicity_stratum
    )
  
  vars_primary <- c(
    outcome,
    "study_group",
    "bmi",
    "antibiotics_during_pregnancy"
  )
  
  vars_ga <- c(
    vars_primary,
    "gestational_age"
  )
  
  dat_primary <- dat0 %>%
    dplyr::select(
      all_of(
        vars_primary
      )
    ) %>%
    drop_na()
  
  dat_ga <- dat0 %>%
    dplyr::select(
      all_of(
        vars_ga
      )
    ) %>%
    drop_na()
  
  dat_primary$study_group <- relevel(
    factor(
      dat_primary$study_group
    ),
    ref = "Control"
  )
  
  dat_ga$study_group <- relevel(
    factor(
      dat_ga$study_group
    ),
    ref = "Control"
  )
  
  dat_primary$antibiotics_during_pregnancy <- droplevels(
    factor(
      dat_primary$antibiotics_during_pregnancy
    )
  )
  
  dat_ga$antibiotics_during_pregnancy <- droplevels(
    factor(
      dat_ga$antibiotics_during_pregnancy
    )
  )
  
  if (
    nlevels(
      droplevels(
        dat_primary$study_group
      )
    ) < 2
  ) {
    stop(
      paste(
        "study_group has fewer than two levels for",
        ethnicity_stratum,
        metric
      )
    )
  }
  
  formula_primary <- as.formula(
    paste(
      outcome,
      "~ study_group + bmi + antibiotics_during_pregnancy"
    )
  )
  
  formula_ga <- as.formula(
    paste(
      outcome,
      "~ study_group + bmi + antibiotics_during_pregnancy + gestational_age"
    )
  )
  
  if (method == "LM") {
    mod_primary <- lm(
      formula_primary,
      data = dat_primary
    )
    
    mod_ga <- lm(
      formula_ga,
      data = dat_ga
    )
    
    result_primary <- extract_strat_lm(
      mod_primary,
      ethnicity_stratum,
      metric,
      "Preeclampsia status",
      "Primary"
    )
    
    result_ga <- extract_strat_lm(
      mod_ga,
      ethnicity_stratum,
      metric,
      "Preeclampsia status",
      "Gestational age-adjusted"
    )
    
  } else {
    
    mod_primary <- MASS::rlm(
      formula_primary,
      data = dat_primary,
      maxit = 100
    )
    
    mod_ga <- MASS::rlm(
      formula_ga,
      data = dat_ga,
      maxit = 100
    )
    
    result_primary <- extract_strat_rlm(
      mod_primary,
      ethnicity_stratum,
      metric,
      "Preeclampsia status",
      "Primary"
    )
    
    result_ga <- extract_strat_rlm(
      mod_ga,
      ethnicity_stratum,
      metric,
      "Preeclampsia status",
      "Gestational age-adjusted"
    )
  }
  
  formulas <- tibble(
    body_site = "Maternal_stool",
    ethnicity_stratum = ethnicity_stratum,
    analysis_method = "Alpha diversity",
    analysis_focus = "Preeclampsia status",
    population = "All_deliveries",
    metric = metric,
    model = c(
      "Primary",
      "Gestational age-adjusted"
    ),
    statistical_method = method,
    n = c(
      nrow(
        dat_primary
      ),
      nrow(
        dat_ga
      )
    ),
    formula = c(
      paste(
        deparse(
          formula_primary
        ),
        collapse = ""
      ),
      paste(
        deparse(
          formula_ga
        ),
        collapse = ""
      )
    )
  )
  
  cat(
    "\nEthnicity-stratified alpha - PE vs Control",
    "\nEthnicity:", ethnicity_stratum,
    "\nMetric:", metric,
    "\nMethod:", method,
    "\nPrimary N:", nrow(dat_primary),
    "\nGA-adjusted N:", nrow(dat_ga),
    "\n"
  )
  
  cat(
    "\nPrimary study-group counts:\n"
  )
  print(
    table(
      dat_primary$study_group
    )
  )
  
  cat(
    "\nGA-adjusted study-group counts:\n"
  )
  print(
    table(
      dat_ga$study_group
    )
  )
  
  list(
    results = bind_rows(
      result_primary,
      result_ga
    ),
    formulas = formulas
  )
}

# Alpha-diversity function for PE severity within ethnicity

run_alpha_ethnicity_severity <- function(
    met,
    ethnicity_stratum,
    metric,
    outcome,
    method = c(
      "LM",
      "RLM"
    )) {
  
  method <- match.arg(
    method
  )
  
  dat0 <- met %>%
    filter(
      ethnicity_analysis == ethnicity_stratum,
      study_group == "PE",
      !is.na(
        severe
      )
    )
  
  vars_primary <- c(
    outcome,
    "severe",
    "bmi",
    "antibiotics_during_pregnancy"
  )
  
  vars_ga <- c(
    vars_primary,
    "gestational_age"
  )
  
  dat_primary <- dat0 %>%
    dplyr::select(
      all_of(
        vars_primary
      )
    ) %>%
    drop_na()
  
  dat_ga <- dat0 %>%
    dplyr::select(
      all_of(
        vars_ga
      )
    ) %>%
    drop_na()
  
  dat_primary$severe <- factor(
    dat_primary$severe,
    levels = c(
      "No",
      "Yes"
    )
  )
  
  dat_ga$severe <- factor(
    dat_ga$severe,
    levels = c(
      "No",
      "Yes"
    )
  )
  
  dat_primary$antibiotics_during_pregnancy <- droplevels(
    factor(
      dat_primary$antibiotics_during_pregnancy
    )
  )
  
  dat_ga$antibiotics_during_pregnancy <- droplevels(
    factor(
      dat_ga$antibiotics_during_pregnancy
    )
  )
  
  if (
    nlevels(
      droplevels(
        dat_primary$severe
      )
    ) < 2
  ) {
    stop(
      paste(
        "Severity has fewer than two levels for",
        ethnicity_stratum,
        metric
      )
    )
  }
  
  formula_primary <- as.formula(
    paste(
      outcome,
      "~ severe + bmi + antibiotics_during_pregnancy"
    )
  )
  
  formula_ga <- as.formula(
    paste(
      outcome,
      "~ severe + bmi + antibiotics_during_pregnancy + gestational_age"
    )
  )
  
  if (method == "LM") {
    mod_primary <- lm(
      formula_primary,
      data = dat_primary
    )
    
    mod_ga <- lm(
      formula_ga,
      data = dat_ga
    )
    
    result_primary <- extract_strat_lm(
      mod_primary,
      ethnicity_stratum,
      metric,
      "Preeclampsia severity",
      "Primary"
    )
    
    result_ga <- extract_strat_lm(
      mod_ga,
      ethnicity_stratum,
      metric,
      "Preeclampsia severity",
      "Gestational age-adjusted"
    )
    
  } else {
    
    mod_primary <- MASS::rlm(
      formula_primary,
      data = dat_primary,
      maxit = 100
    )
    
    mod_ga <- MASS::rlm(
      formula_ga,
      data = dat_ga,
      maxit = 100
    )
    
    result_primary <- extract_strat_rlm(
      mod_primary,
      ethnicity_stratum,
      metric,
      "Preeclampsia severity",
      "Primary"
    )
    
    result_ga <- extract_strat_rlm(
      mod_ga,
      ethnicity_stratum,
      metric,
      "Preeclampsia severity",
      "Gestational age-adjusted"
    )
  }
  
  formulas <- tibble(
    body_site = "Maternal_stool",
    ethnicity_stratum = ethnicity_stratum,
    analysis_method = "Alpha diversity",
    analysis_focus = "Preeclampsia severity",
    population = "All_deliveries",
    metric = metric,
    model = c(
      "Primary",
      "Gestational age-adjusted"
    ),
    statistical_method = method,
    n = c(
      nrow(
        dat_primary
      ),
      nrow(
        dat_ga
      )
    ),
    formula = c(
      paste(
        deparse(
          formula_primary
        ),
        collapse = ""
      ),
      paste(
        deparse(
          formula_ga
        ),
        collapse = ""
      )
    )
  )
  
  cat(
    "\nEthnicity-stratified alpha - PE severity",
    "\nEthnicity:", ethnicity_stratum,
    "\nMetric:", metric,
    "\nMethod:", method,
    "\nPrimary N:", nrow(dat_primary),
    "\nGA-adjusted N:", nrow(dat_ga),
    "\n"
  )
  
  cat(
    "\nPrimary severity counts:\n"
  )
  print(
    table(
      dat_primary$severe
    )
  )
  
  cat(
    "\nGA-adjusted severity counts:\n"
  )
  print(
    table(
      dat_ga$severe
    )
  )
  
  list(
    results = bind_rows(
      result_primary,
      result_ga
    ),
    formulas = formulas
  )
}

# Run all ethnicity-stratified alpha-diversity analyses

ethnicity_levels_to_run <- c(
  "Hispanic",
  "White"
)

alpha_ethnicity_objects <- list()

for (eth in ethnicity_levels_to_run) {
  
  alpha_ethnicity_objects[[paste0(
    eth,
    "_status_observed"
  )]] <- run_alpha_ethnicity_status(
    met = met_strat,
    ethnicity_stratum = eth,
    metric = "Observed richness",
    outcome = "observed_original",
    method = "LM"
  )
  
  alpha_ethnicity_objects[[paste0(
    eth,
    "_status_shannon"
  )]] <- run_alpha_ethnicity_status(
    met = met_strat,
    ethnicity_stratum = eth,
    metric = "Shannon",
    outcome = "shannon_boxcox",
    method = "LM"
  )
  
  alpha_ethnicity_objects[[paste0(
    eth,
    "_status_simpson"
  )]] <- run_alpha_ethnicity_status(
    met = met_strat,
    ethnicity_stratum = eth,
    metric = "Simpson",
    outcome = "simpson_original",
    method = "RLM"
  )
  
  alpha_ethnicity_objects[[paste0(
    eth,
    "_severity_observed"
  )]] <- run_alpha_ethnicity_severity(
    met = met_strat,
    ethnicity_stratum = eth,
    metric = "Observed richness",
    outcome = "observed_original",
    method = "LM"
  )
  
  alpha_ethnicity_objects[[paste0(
    eth,
    "_severity_shannon"
  )]] <- run_alpha_ethnicity_severity(
    met = met_strat,
    ethnicity_stratum = eth,
    metric = "Shannon",
    outcome = "shannon_boxcox",
    method = "LM"
  )
  
  alpha_ethnicity_objects[[paste0(
    eth,
    "_severity_simpson"
  )]] <- run_alpha_ethnicity_severity(
    met = met_strat,
    ethnicity_stratum = eth,
    metric = "Simpson",
    outcome = "simpson_original",
    method = "RLM"
  )
}

alpha_ethnicity_results <- bind_rows(
  lapply(
    alpha_ethnicity_objects,
    `[[`,
    "results"
  )
)

alpha_ethnicity_formulas <- bind_rows(
  lapply(
    alpha_ethnicity_objects,
    `[[`,
    "formulas"
  )
)

alpha_ethnicity_exposure_effects <- alpha_ethnicity_results %>%
  filter(
    term %in% c(
      "study_groupPE",
      "severeYes"
    )
  ) %>%
  arrange(
    ethnicity_stratum,
    analysis_focus,
    metric,
    model
  )

# Helper to run one ethnicity-stratified PERMANOVA model

run_one_strat_permanova <- function(
    phy,
    ids,
    distance_method,
    formula_type,
    ethnicity_stratum,
    comparison,
    model_type,
    permutations = 9999,
    seed = 711) {
  
  phy_cc <- prune_samples(
    ids,
    phy
  )
  
  phy_cc <- prune_taxa(
    taxa_sums(
      phy_cc
    ) > 0,
    phy_cc
  )
  
  met_cc <- data.frame(
    sample_data(
      phy_cc
    )
  )
  
  if (comparison == "Preeclampsia status") {
    met_cc$study_group <- relevel(
      factor(
        met_cc$study_group
      ),
      ref = "Control"
    )
  } else {
    met_cc$severe <- factor(
      met_cc$severe,
      levels = c(
        "No",
        "Yes"
      )
    )
  }
  
  met_cc$antibiotics_during_pregnancy <- droplevels(
    factor(
      met_cc$antibiotics_during_pregnancy
    )
  )
  
  if (distance_method == "jaccard") {
    dist_obj <- phyloseq::distance(
      phy_cc,
      method = "jaccard",
      binary = TRUE
    )
  } else if (distance_method == "bray") {
    dist_obj <- phyloseq::distance(
      phy_cc,
      method = "bray"
    )
  } else {
    stop(
      "distance_method must be either 'bray' or 'jaccard'"
    )
  }
  
  if (
    comparison == "Preeclampsia status" &&
    formula_type == "primary"
  ) {
    model_formula <-
      dist_obj ~
      study_group +
      bmi +
      antibiotics_during_pregnancy
  }
  
  if (
    comparison == "Preeclampsia status" &&
    formula_type == "ga"
  ) {
    model_formula <-
      dist_obj ~
      study_group +
      bmi +
      antibiotics_during_pregnancy +
      gestational_age
  }
  
  if (
    comparison == "Preeclampsia severity" &&
    formula_type == "primary"
  ) {
    model_formula <-
      dist_obj ~
      severe +
      bmi +
      antibiotics_during_pregnancy
  }
  
  if (
    comparison == "Preeclampsia severity" &&
    formula_type == "ga"
  ) {
    model_formula <-
      dist_obj ~
      severe +
      bmi +
      antibiotics_during_pregnancy +
      gestational_age
  }
  
  set.seed(
    seed
  )
  
  perm <- vegan::adonis2(
    model_formula,
    data = met_cc,
    permutations = permutations,
    by = "margin"
  )
  
  perm_table <- as.data.frame(
    perm
  ) %>%
    rownames_to_column(
      "term"
    ) %>%
    mutate(
      body_site = "Maternal_stool",
      ethnicity_stratum = ethnicity_stratum,
      analysis_method = "PERMANOVA",
      analysis_focus = comparison,
      population = "All_deliveries",
      metric = ifelse(
        distance_method == "bray",
        "Bray-Curtis",
        "Jaccard"
      ),
      model = model_type,
      statistical_method = "PERMANOVA",
      n = nsamples(
        phy_cc
      ),
      model_term = pretty_stool_term(
        term
      )
    ) %>%
    relocate(
      body_site,
      ethnicity_stratum,
      analysis_method,
      analysis_focus,
      population,
      metric,
      model,
      statistical_method,
      n,
      term,
      model_term
    )
  
  list(
    table = perm_table,
    dist_obj = dist_obj,
    metadata = met_cc
  )
}

# Beta-diversity function within ethnicity

run_beta_ethnicity <- function(
    phy,
    ethnicity_stratum,
    comparison,
    distance_method,
    permutations = 9999,
    seed = 711) {
  
  met_all <- data.frame(
    sample_data(
      phy
    )
  )
  
  if (comparison == "Preeclampsia status") {
    
    eligible <- met_all %>%
      mutate(
        sample_id = rownames(
          met_all
        )
      ) %>%
      filter(
        ethnicity_analysis == ethnicity_stratum
      )
    
    exposure_var <- "study_group"
    
  } else {
    
    eligible <- met_all %>%
      mutate(
        sample_id = rownames(
          met_all
        )
      ) %>%
      filter(
        ethnicity_analysis == ethnicity_stratum,
        study_group == "PE",
        !is.na(
          severe
        )
      )
    
    exposure_var <- "severe"
  }
  
  primary_vars <- c(
    exposure_var,
    "bmi",
    "antibiotics_during_pregnancy"
  )
  
  ga_vars <- c(
    primary_vars,
    "gestational_age"
  )
  
  primary_ids <- eligible %>%
    filter(
      complete.cases(
        .[
          ,
          primary_vars,
          drop = FALSE
        ]
      )
    ) %>%
    pull(
      sample_id
    )
  
  ga_ids <- eligible %>%
    filter(
      complete.cases(
        .[
          ,
          ga_vars,
          drop = FALSE
        ]
      )
    ) %>%
    pull(
      sample_id
    )
  
  primary_result <- run_one_strat_permanova(
    phy = phy,
    ids = primary_ids,
    distance_method = distance_method,
    formula_type = "primary",
    ethnicity_stratum = ethnicity_stratum,
    comparison = comparison,
    model_type = "Primary",
    permutations = permutations,
    seed = seed
  )
  
  ga_result <- run_one_strat_permanova(
    phy = phy,
    ids = ga_ids,
    distance_method = distance_method,
    formula_type = "ga",
    ethnicity_stratum = ethnicity_stratum,
    comparison = comparison,
    model_type = "Gestational age-adjusted",
    permutations = permutations,
    seed = seed
  )
  
  # Beta dispersion is evaluated separately on the primary analysis population
  
  if (comparison == "Preeclampsia status") {
    grouping <- primary_result$metadata$study_group
    dispersion_term <- "Preeclampsia"
  } else {
    grouping <- primary_result$metadata$severe
    dispersion_term <- "Preeclampsia with severe features"
  }
  
  dispersion_model <- vegan::betadisper(
    primary_result$dist_obj,
    grouping
  )
  
  set.seed(
    seed
  )
  
  dispersion_test <- vegan::permutest(
    dispersion_model,
    permutations = permutations
  )
  
  dispersion_table <- tibble(
    body_site = "Maternal_stool",
    ethnicity_stratum = ethnicity_stratum,
    analysis_method = "Beta dispersion",
    analysis_focus = comparison,
    population = "All_deliveries",
    metric = ifelse(
      distance_method == "bray",
      "Bray-Curtis",
      "Jaccard"
    ),
    model = "Dispersion test",
    statistical_method = "Betadisper",
    n = length(
      grouping
    ),
    model_term = dispersion_term,
    F = dispersion_test$tab[
      1,
      "F"
    ],
    p_value = dispersion_test$tab[
      1,
      "Pr(>F)"
    ]
  )
  
  dispersion_values <- tibble(
    sample_id = names(
      dispersion_model$distances
    ),
    group = as.character(
      grouping
    ),
    distance_to_centroid = as.numeric(
      dispersion_model$distances
    ),
    body_site = "Maternal_stool",
    ethnicity_stratum = ethnicity_stratum,
    analysis_focus = comparison,
    metric = ifelse(
      distance_method == "bray",
      "Bray-Curtis",
      "Jaccard"
    ),
    population = "All_deliveries"
  )
  
  formulas <- tibble(
    body_site = "Maternal_stool",
    ethnicity_stratum = ethnicity_stratum,
    analysis_method = "PERMANOVA",
    analysis_focus = comparison,
    population = "All_deliveries",
    metric = ifelse(
      distance_method == "bray",
      "Bray-Curtis",
      "Jaccard"
    ),
    model = c(
      "Primary",
      "Gestational age-adjusted"
    ),
    statistical_method = "PERMANOVA",
    n = c(
      length(
        primary_ids
      ),
      length(
        ga_ids
      )
    ),
    formula = if (
      comparison == "Preeclampsia status"
    ) {
      c(
        "study_group + bmi + antibiotics_during_pregnancy",
        "study_group + bmi + antibiotics_during_pregnancy + gestational_age"
      )
    } else {
      c(
        "severe + bmi + antibiotics_during_pregnancy",
        "severe + bmi + antibiotics_during_pregnancy + gestational_age"
      )
    }
  )
  
  cat(
    "\nEthnicity-stratified beta",
    "\nEthnicity:", ethnicity_stratum,
    "\nComparison:", comparison,
    "\nDistance:", distance_method,
    "\nPrimary N:", length(primary_ids),
    "\nGA-adjusted N:", length(ga_ids),
    "\n"
  )
  
  list(
    permanova = bind_rows(
      primary_result$table,
      ga_result$table
    ),
    dispersion = dispersion_table,
    dispersion_values = dispersion_values,
    formulas = formulas
  )
}

# Run all ethnicity-stratified beta-diversity analyses

beta_ethnicity_objects <- list()

for (eth in ethnicity_levels_to_run) {
  
  for (dist_method in c(
    "bray",
    "jaccard"
  )) {
    
    beta_ethnicity_objects[[paste(
      eth,
      "status",
      dist_method,
      sep = "_"
    )]] <- run_beta_ethnicity(
      phy = phy,
      ethnicity_stratum = eth,
      comparison = "Preeclampsia status",
      distance_method = dist_method
    )
    
    beta_ethnicity_objects[[paste(
      eth,
      "severity",
      dist_method,
      sep = "_"
    )]] <- run_beta_ethnicity(
      phy = phy,
      ethnicity_stratum = eth,
      comparison = "Preeclampsia severity",
      distance_method = dist_method
    )
  }
}

beta_ethnicity_results <- bind_rows(
  lapply(
    beta_ethnicity_objects,
    `[[`,
    "permanova"
  )
)

beta_ethnicity_dispersion <- bind_rows(
  lapply(
    beta_ethnicity_objects,
    `[[`,
    "dispersion"
  )
)

beta_ethnicity_dispersion_values <- bind_rows(
  lapply(
    beta_ethnicity_objects,
    `[[`,
    "dispersion_values"
  )
)

beta_ethnicity_formulas <- bind_rows(
  lapply(
    beta_ethnicity_objects,
    `[[`,
    "formulas"
  )
)

beta_ethnicity_exposure_effects <- beta_ethnicity_results %>%
  filter(
    term %in% c(
      "study_group",
      "severe"
    )
  ) %>%
  transmute(
    body_site,
    ethnicity_stratum,
    analysis_method,
    analysis_focus,
    population,
    metric,
    model,
    statistical_method,
    n,
    term,
    model_term,
    R2,
    F,
    p_value = `Pr(>F)`
  ) %>%
  arrange(
    ethnicity_stratum,
    analysis_focus,
    metric,
    model
  )

# Save ethnicity-stratified alpha-diversity outputs

write_csv(
  alpha_ethnicity_results,
  file.path(
    out_dir,
    "stool_alpha_ethnicity_stratified_all_coefficients_95CI.csv"
  )
)

write_csv(
  alpha_ethnicity_exposure_effects,
  file.path(
    out_dir,
    "stool_alpha_ethnicity_stratified_exposure_effects_95CI.csv"
  )
)

write_csv(
  alpha_ethnicity_formulas,
  file.path(
    out_dir,
    "stool_alpha_ethnicity_stratified_model_formulas.csv"
  )
)

# Save ethnicity-stratified beta-diversity outputs

write_csv(
  beta_ethnicity_results,
  file.path(
    out_dir,
    "stool_beta_ethnicity_stratified_PERMANOVA_all_models.csv"
  )
)

write_csv(
  beta_ethnicity_exposure_effects,
  file.path(
    out_dir,
    "stool_beta_ethnicity_stratified_exposure_effects_R2.csv"
  )
)

write_csv(
  beta_ethnicity_dispersion,
  file.path(
    out_dir,
    "stool_beta_ethnicity_stratified_dispersion.csv"
  )
)

write_csv(
  beta_ethnicity_dispersion_values,
  file.path(
    out_dir,
    "stool_beta_ethnicity_stratified_dispersion_values.csv"
  )
)

write_csv(
  beta_ethnicity_formulas,
  file.path(
    out_dir,
    "stool_beta_ethnicity_stratified_model_formulas.csv"
  )
)

# Build a Table S4-style combined table for the maternal gut results

alpha_tableS4_ready <- alpha_ethnicity_results %>%
  transmute(
    `Body site` = "Gut",
    `Ethnicity stratum` = ethnicity_stratum,
    `Analysis method` = analysis_method,
    `Analysis focus` = analysis_focus,
    `Population` = "All deliveries",
    `Metric` = metric,
    `Model` = model,
    `Statistical method` = statistical_method,
    `N` = n,
    `Model term` = model_term,
    `Estimate` = estimate,
    `95% CI lower` = conf.low,
    `95% CI upper` = conf.high,
    `R2` = NA_real_,
    `F` = NA_real_,
    `P value` = p.value
  )

beta_tableS4_ready <- beta_ethnicity_results %>%
  filter(
    !term %in% c(
      "Residual",
      "Total"
    )
  ) %>%
  transmute(
    `Body site` = "Gut",
    `Ethnicity stratum` = ethnicity_stratum,
    `Analysis method` = analysis_method,
    `Analysis focus` = analysis_focus,
    `Population` = "All deliveries",
    `Metric` = metric,
    `Model` = model,
    `Statistical method` = statistical_method,
    `N` = n,
    `Model term` = model_term,
    `Estimate` = NA_real_,
    `95% CI lower` = NA_real_,
    `95% CI upper` = NA_real_,
    `R2` = R2,
    `F` = F,
    `P value` = `Pr(>F)`
  )

dispersion_tableS4_ready <- beta_ethnicity_dispersion %>%
  transmute(
    `Body site` = "Gut",
    `Ethnicity stratum` = ethnicity_stratum,
    `Analysis method` = analysis_method,
    `Analysis focus` = analysis_focus,
    `Population` = "All deliveries",
    `Metric` = metric,
    `Model` = model,
    `Statistical method` = statistical_method,
    `N` = n,
    `Model term` = model_term,
    `Estimate` = NA_real_,
    `95% CI lower` = NA_real_,
    `95% CI upper` = NA_real_,
    `R2` = NA_real_,
    `F` = F,
    `P value` = p_value
  )

stool_ethnicity_tableS4_ready <- bind_rows(
  alpha_tableS4_ready,
  beta_tableS4_ready,
  dispersion_tableS4_ready
) %>%
  arrange(
    `Ethnicity stratum`,
    `Analysis focus`,
    `Analysis method`,
    `Metric`,
    `Model`,
    `Model term`
  )

write_csv(
  stool_ethnicity_tableS4_ready,
  file.path(
    out_dir,
    "stool_ethnicity_stratified_TableS4_ready.csv"
  )
)

# Print the key ethnicity-stratified effects

cat(
  "\nMaternal stool ethnicity-stratified alpha exposure effects\n"
)

print(
  alpha_ethnicity_exposure_effects
)

cat(
  "\nMaternal stool ethnicity-stratified beta exposure effects\n"
)

print(
  beta_ethnicity_exposure_effects
)

cat(
  "\nMaternal stool ethnicity-stratified beta dispersion\n"
)

print(
  beta_ethnicity_dispersion
)

cat(
  "\nEthnicity-stratified files generated:\n"
)

print(
  list.files(
    out_dir,
    pattern = "^stool_.*ethnicity_stratified",
    full.names = TRUE
  )
)

#### Interaction study_group * racial/ethnic group ####
# Maternal gut microbiome:
# Preeclampsia × racial/ethnic group interaction analyses
# Hispanic/Latina vs Non-Hispanic White/Caucasian
# Alpha and beta diversity

library(phyloseq)
library(vegan)
library(MASS)
library(dplyr)
library(tidyr)
library(broom)
library(readr)
library(tibble)

# Load maternal gut phyloseq object

phy_path <- file.path(PROJECT_DIR, "WGS", "obj", "phy_mtph_gtdb.rds")

phy <- readRDS(
  phy_path
)

# Output directory

project_dir <- dirname(
  dirname(
    phy_path
  )
)

out_dir <- file.path(
  project_dir,
  "results",
  "reviewer_adjusted_models",
  as.character(
    Sys.Date()
  ),
  "PE_race_ethnicity_interaction"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

cat(
  "\nGut interaction results will be saved in:\n",
  out_dir,
  "\n"
)

# Create racial/ethnic group

create_race_ethnicity_group <- function(dat) {
  
  dat %>%
    mutate(
      race = as.character(race),
      ethnicity = as.character(ethnicity),
      
      race_ethnicity_group = case_when(
        
        ethnicity == "Hispanic/Latina" ~
          "Hispanic/Latina",
        
        ethnicity == "Not Hispanic/Latina" &
          race == "White/Caucasian" ~
          "Non-Hispanic White/Caucasian",
        
        TRUE ~ NA_character_
      )
    )
}

# Check eligible gut cohort

met_all_check <- data.frame(
  sample_data(phy)
) %>%
  create_race_ethnicity_group()

cat(
  "\nGut samples eligible for racial/ethnic subgroup analysis:\n"
)

print(
  met_all_check %>%
    filter(
      !is.na(race_ethnicity_group)
    ) %>%
    count(
      race_ethnicity_group,
      study_group
    )
)

if ("study_id" %in% names(met_all_check)) {
  
  cat(
    "\nUnique eligible women with gut metadata:",
    met_all_check %>%
      filter(
        !is.na(race_ethnicity_group)
      ) %>%
      summarise(
        n = n_distinct(study_id)
      ) %>%
      pull(n),
    "\n"
  )
}

# Prepare metadata

prepare_gut_metadata <- function(
    phy) {
  
  met <- data.frame(
    sample_data(
      phy
    )
  )
  
  vars_required <- c(
    "study_group",
    "race",
    "ethnicity",
    "bmi",
    "antibiotics_during_pregnancy",
    "obs",
    "sha",
    "simp"
  )
  
  missing_vars <- setdiff(
    vars_required,
    colnames(
      met
    )
  )
  
  if (
    length(
      missing_vars
    ) > 0
  ) {
    
    stop(
      paste(
        "Missing variables:",
        paste(
          missing_vars,
          collapse = ", "
        )
      )
    )
  }
  
  # Create racial/ethnic group and
  # retain only the two prespecified groups
  
  met <- met %>%
    create_race_ethnicity_group() %>%
    filter(
      !is.na(race_ethnicity_group)
    )
  
  met$study_group <- relevel(
    factor(
      met$study_group
    ),
    ref = "Control"
  )
  
  met$race_ethnicity_group <- relevel(
    factor(
      met$race_ethnicity_group
    ),
    ref = "Non-Hispanic White/Caucasian"
  )
  
  met$antibiotics_during_pregnancy <- factor(
    met$antibiotics_during_pregnancy
  )
  
  # Box-Cox transformation for Shannon
  
  sha_nonmissing <- met %>%
    filter(
      !is.na(
        sha
      )
    )
  
  bc <- MASS::boxcox(
    sha ~ 1,
    data = sha_nonmissing,
    plotit = FALSE
  )
  
  lambda_shannon <- bc$x[
    which.max(
      bc$y
    )
  ]
  
  if (
    abs(
      lambda_shannon
    ) < 0.001
  ) {
    
    met$shannon_transformed <- log(
      met$sha
    )
    
  } else {
    
    met$shannon_transformed <-
      (
        met$sha^lambda_shannon - 1
      ) /
      lambda_shannon
  }
  
  attr(
    met,
    "lambda_shannon"
  ) <- lambda_shannon
  
  cat(
    "\nShannon Box-Cox lambda:",
    lambda_shannon,
    "\n"
  )
  
  cat(
    "\nStudy group × racial/ethnic group:\n"
  )
  
  print(
    table(
      met$study_group,
      met$race_ethnicity_group,
      useNA = "ifany"
    )
  )
  
  cat(
    "\nNumber of gut samples before complete-case filtering:",
    nrow(met),
    "\n"
  )
  
  return(
    met
  )
}

met <- prepare_gut_metadata(
  phy
)

# Identify interaction coefficient

find_interaction_term <- function(
    model) {
  
  model_terms <- names(
    coef(
      model
    )
  )
  
  interaction_term <- model_terms[
    grepl(
      "study_groupPE",
      model_terms,
      fixed = TRUE
    ) &
      grepl(
        "race_ethnicity_group",
        model_terms,
        fixed = TRUE
      ) &
      grepl(
        ":",
        model_terms,
        fixed = TRUE
      )
  ]
  
  if (
    length(
      interaction_term
    ) != 1
  ) {
    
    stop(
      paste(
        "Could not uniquely identify the study_group × racial/ethnic group interaction term.",
        "Terms:",
        paste(
          model_terms,
          collapse = ", "
        )
      )
    )
  }
  
  return(
    interaction_term
  )
}

# Extract simple PE effects within each racial/ethnic group

extract_simple_effects <- function(
    model,
    metric,
    method) {
  
  b <- coef(
    model
  )
  
  V <- vcov(
    model
  )
  
  group_term <- "study_groupPE"
  
  interaction_term <- find_interaction_term(
    model
  )
  
  # Reference group:
  # Non-Hispanic White/Caucasian
  
  est_white <- b[
    group_term
  ]
  
  se_white <- sqrt(
    V[
      group_term,
      group_term
    ]
  )
  
  # Hispanic/Latina:
  # PE main effect + interaction
  
  est_hispanic <-
    b[
      group_term
    ] +
    b[
      interaction_term
    ]
  
  var_hispanic <-
    V[
      group_term,
      group_term
    ] +
    V[
      interaction_term,
      interaction_term
    ] +
    2 *
    V[
      group_term,
      interaction_term
    ]
  
  se_hispanic <- sqrt(
    var_hispanic
  )
  
  if (
    method == "LM"
  ) {
    
    df_model <- df.residual(
      model
    )
    
    critical_value <- qt(
      0.975,
      df = df_model
    )
    
    stat_white <-
      est_white /
      se_white
    
    stat_hispanic <-
      est_hispanic /
      se_hispanic
    
    p_white <- 2 * pt(
      abs(
        stat_white
      ),
      df = df_model,
      lower.tail = FALSE
    )
    
    p_hispanic <- 2 * pt(
      abs(
        stat_hispanic
      ),
      df = df_model,
      lower.tail = FALSE
    )
    
    model_n <- nobs(
      model
    )
    
  } else {
    
    critical_value <- 1.96
    
    stat_white <-
      est_white /
      se_white
    
    stat_hispanic <-
      est_hispanic /
      se_hispanic
    
    p_white <- 2 * pnorm(
      abs(
        stat_white
      ),
      lower.tail = FALSE
    )
    
    p_hispanic <- 2 * pnorm(
      abs(
        stat_hispanic
      ),
      lower.tail = FALSE
    )
    
    model_n <- length(
      model$residuals
    )
  }
  
  tibble(
    body_site = "Maternal gut",
    metric = metric,
    method = method,
    
    race_ethnicity_group = c(
      "Non-Hispanic White/Caucasian",
      "Hispanic/Latina"
    ),
    
    contrast = "Preeclampsia - Control",
    
    n = model_n,
    
    estimate = c(
      est_white,
      est_hispanic
    ),
    
    std_error = c(
      se_white,
      se_hispanic
    ),
    
    statistic = c(
      stat_white,
      stat_hispanic
    ),
    
    p_value = c(
      p_white,
      p_hispanic
    ),
    
    conf_low = c(
      est_white -
        critical_value *
        se_white,
      
      est_hispanic -
        critical_value *
        se_hispanic
    ),
    
    conf_high = c(
      est_white +
        critical_value *
        se_white,
      
      est_hispanic +
        critical_value *
        se_hispanic
    )
  )
}

# Run one alpha-diversity interaction model
#
# Observed richness: raw scale
# Shannon: Box-Cox transformed
# Simpson: robust linear regression

run_alpha_interaction <- function(
    met,
    metric,
    outcome,
    robust = FALSE) {
  
  vars_needed <- c(
    outcome,
    "study_group",
    "race_ethnicity_group",
    "bmi",
    "antibiotics_during_pregnancy"
  )
  
  dat <- met %>%
    dplyr::select(
      all_of(
        vars_needed
      )
    ) %>%
    drop_na()
  
  dat$study_group <- relevel(
    droplevels(
      factor(
        dat$study_group
      )
    ),
    ref = "Control"
  )
  
  dat$race_ethnicity_group <- relevel(
    droplevels(
      factor(
        dat$race_ethnicity_group
      )
    ),
    ref = "Non-Hispanic White/Caucasian"
  )
  
  dat$antibiotics_during_pregnancy <- droplevels(
    factor(
      dat$antibiotics_during_pregnancy
    )
  )
  
  cat(
    "\n",
    metric,
    "complete-case distribution:\n"
  )
  
  print(
    table(
      dat$study_group,
      dat$race_ethnicity_group
    )
  )
  
  formula_no_interaction <- as.formula(
    paste(
      outcome,
      "~ study_group + race_ethnicity_group + bmi + antibiotics_during_pregnancy"
    )
  )
  
  formula_interaction <- as.formula(
    paste(
      outcome,
      "~ study_group * race_ethnicity_group + bmi + antibiotics_during_pregnancy"
    )
  )
  
  # LM version for residual diagnostics
  
  mod_lm <- lm(
    formula_interaction,
    data = dat
  )
  
  shapiro_p <- shapiro.test(
    residuals(
      mod_lm
    )
  )$p.value
  
  if (
    robust
  ) {
    
    mod <- MASS::rlm(
      formula_interaction,
      data = dat,
      maxit = 100
    )
    
    method <- "RLM"
    
    s <- summary(
      mod
    )
    
    coef_table <- as.data.frame(
      s$coefficients
    )
    
    coef_table$term <- rownames(
      coef_table
    )
    
    names(
      coef_table
    )[1:3] <- c(
      "estimate",
      "std_error",
      "statistic"
    )
    
    interaction_term <- find_interaction_term(
      mod
    )
    
    interaction_row <- coef_table %>%
      filter(
        term == interaction_term
      )
    
    interaction_estimate <-
      interaction_row$estimate
    
    interaction_se <-
      interaction_row$std_error
    
    interaction_statistic <-
      interaction_row$statistic
    
    interaction_p <- 2 * pnorm(
      abs(
        interaction_statistic
      ),
      lower.tail = FALSE
    )
    
    conf_low <-
      interaction_estimate -
      1.96 *
      interaction_se
    
    conf_high <-
      interaction_estimate +
      1.96 *
      interaction_se
    
    interaction_F <- NA_real_
    
  } else {
    
    mod0 <- lm(
      formula_no_interaction,
      data = dat
    )
    
    mod <- mod_lm
    
    method <- "LM"
    
    comparison <- anova(
      mod0,
      mod
    )
    
    interaction_p <- comparison$`Pr(>F)`[
      2
    ]
    
    interaction_F <- comparison$F[
      2
    ]
    
    interaction_term <- find_interaction_term(
      mod
    )
    
    tidy_mod <- broom::tidy(
      mod,
      conf.int = TRUE,
      conf.level = 0.95
    )
    
    interaction_row <- tidy_mod %>%
      filter(
        term == interaction_term
      )
    
    interaction_estimate <-
      interaction_row$estimate
    
    interaction_se <-
      interaction_row$std.error
    
    interaction_statistic <-
      interaction_row$statistic
    
    conf_low <-
      interaction_row$conf.low
    
    conf_high <-
      interaction_row$conf.high
  }
  
  interaction_result <- tibble(
    body_site = "Maternal gut",
    analysis = "Alpha diversity",
    metric = metric,
    method = method,
    n = nrow(
      dat
    ),
    
    interaction =
      "Preeclampsia × racial/ethnic group",
    
    interaction_term =
      interaction_term,
    
    estimate =
      interaction_estimate,
    
    std_error =
      interaction_se,
    
    conf_low =
      conf_low,
    
    conf_high =
      conf_high,
    
    R2 =
      NA_real_,
    
    F =
      interaction_F,
    
    statistic =
      interaction_statistic,
    
    p_value =
      interaction_p
  )
  
  diagnostics <- tibble(
    body_site = "Maternal gut",
    metric = metric,
    method = method,
    n = nrow(
      dat
    ),
    LM_residual_shapiro_p =
      shapiro_p
  )
  
  simple_effects <- extract_simple_effects(
    model = mod,
    metric = metric,
    method = method
  )
  
  formulas <- tibble(
    body_site = "Maternal gut",
    analysis = "Alpha diversity",
    metric = metric,
    method = method,
    n = nrow(
      dat
    ),
    formula = paste(
      deparse(
        formula_interaction
      ),
      collapse = ""
    )
  )
  
  list(
    interaction = interaction_result,
    simple_effects = simple_effects,
    diagnostics = diagnostics,
    formulas = formulas
  )
}

# Run alpha-diversity interaction analyses

alpha_objects <- list(
  
  run_alpha_interaction(
    met = met,
    metric = "Observed richness",
    outcome = "obs",
    robust = FALSE
  ),
  
  run_alpha_interaction(
    met = met,
    metric = "Shannon",
    outcome = "shannon_transformed",
    robust = FALSE
  ),
  
  run_alpha_interaction(
    met = met,
    metric = "Simpson",
    outcome = "simp",
    robust = TRUE
  )
)

alpha_interaction_results <- bind_rows(
  lapply(
    alpha_objects,
    `[[`,
    "interaction"
  )
) %>%
  mutate(
    q_value = p.adjust(
      p_value,
      method = "BH"
    )
  )

alpha_simple_effects <- bind_rows(
  lapply(
    alpha_objects,
    `[[`,
    "simple_effects"
  )
)

alpha_diagnostics <- bind_rows(
  lapply(
    alpha_objects,
    `[[`,
    "diagnostics"
  )
)

alpha_formulas <- bind_rows(
  lapply(
    alpha_objects,
    `[[`,
    "formulas"
  )
)

# Beta-diversity interaction analysis

run_beta_interaction <- function(
    phy,
    distance_method,
    permutations = 9999,
    seed = 711) {
  
  met <- data.frame(
    sample_data(
      phy
    )
  )
  
  vars_needed <- c(
    "study_group",
    "race",
    "ethnicity",
    "bmi",
    "antibiotics_during_pregnancy"
  )
  
  missing_vars <- setdiff(
    vars_needed,
    colnames(
      met
    )
  )
  
  if (
    length(
      missing_vars
    ) > 0
  ) {
    
    stop(
      paste(
        "Missing variables:",
        paste(
          missing_vars,
          collapse = ", "
        )
      )
    )
  }
  
  # Create racial/ethnic group
  
  met <- met %>%
    create_race_ethnicity_group()
  
  # Retain only Hispanic/Latina and
  # Non-Hispanic White/Caucasian
  
  analysis_ids <- rownames(
    met
  )[
    !is.na(
      met$race_ethnicity_group
    )
  ]
  
  phy_sub <- prune_samples(
    analysis_ids,
    phy
  )
  
  phy_sub <- prune_taxa(
    taxa_sums(
      phy_sub
    ) > 0,
    phy_sub
  )
  
  met_sub <- data.frame(
    sample_data(
      phy_sub
    )
  ) %>%
    create_race_ethnicity_group()
  
  vars_complete <- c(
    "study_group",
    "race_ethnicity_group",
    "bmi",
    "antibiotics_during_pregnancy"
  )
  
  complete_ids <- rownames(
    met_sub
  )[
    complete.cases(
      met_sub[
        ,
        vars_complete,
        drop = FALSE
      ]
    )
  ]
  
  phy_cc <- prune_samples(
    complete_ids,
    phy_sub
  )
  
  phy_cc <- prune_taxa(
    taxa_sums(
      phy_cc
    ) > 0,
    phy_cc
  )
  
  met_cc <- data.frame(
    sample_data(
      phy_cc
    )
  ) %>%
    create_race_ethnicity_group()
  
  met_cc$study_group <- relevel(
    droplevels(
      factor(
        met_cc$study_group
      )
    ),
    ref = "Control"
  )
  
  met_cc$race_ethnicity_group <- relevel(
    droplevels(
      factor(
        met_cc$race_ethnicity_group
      )
    ),
    ref = "Non-Hispanic White/Caucasian"
  )
  
  met_cc$antibiotics_during_pregnancy <- droplevels(
    factor(
      met_cc$antibiotics_during_pregnancy
    )
  )
  
  cat(
    "\nBeta diversity:",
    distance_method,
    "\n"
  )
  
  print(
    table(
      met_cc$study_group,
      met_cc$race_ethnicity_group
    )
  )
  
  cat(
    "N used:",
    nsamples(
      phy_cc
    ),
    "\n"
  )
  
  if (
    distance_method == "jaccard"
  ) {
    
    dist_obj <- phyloseq::distance(
      phy_cc,
      method = "jaccard",
      binary = TRUE
    )
    
    metric_label <- "Jaccard"
    
  } else if (
    distance_method == "bray"
  ) {
    
    dist_obj <- phyloseq::distance(
      phy_cc,
      method = "bray"
    )
    
    metric_label <- "Bray-Curtis"
    
  } else {
    
    stop(
      "distance_method must be either 'bray' or 'jaccard'."
    )
  }
  
  formula_interaction <-
    dist_obj ~
    study_group *
    race_ethnicity_group +
    bmi +
    antibiotics_during_pregnancy
  
  set.seed(
    seed
  )
  
  permanova <- vegan::adonis2(
    formula_interaction,
    data = met_cc,
    permutations = permutations,
    by = "margin"
  )
  
  permanova_table <- as.data.frame(
    permanova
  ) %>%
    tibble::rownames_to_column(
      "term"
    ) %>%
    mutate(
      body_site = "Maternal gut",
      analysis = "Beta diversity",
      metric = metric_label,
      method = "PERMANOVA",
      n = nsamples(
        phy_cc
      )
    ) %>%
    relocate(
      body_site,
      analysis,
      metric,
      method,
      n,
      term
    )
  
  interaction_row <- permanova_table %>%
    filter(
      grepl(
        "study_group:race_ethnicity_group|race_ethnicity_group:study_group",
        term
      )
    )
  
  if (
    nrow(
      interaction_row
    ) != 1
  ) {
    
    stop(
      paste(
        "Could not uniquely identify the PERMANOVA interaction row for",
        metric_label
      )
    )
  }
  
  interaction_result <- tibble(
    body_site = "Maternal gut",
    analysis = "Beta diversity",
    metric = metric_label,
    method = "PERMANOVA",
    n = nsamples(
      phy_cc
    ),
    
    interaction =
      "Preeclampsia × racial/ethnic group",
    
    interaction_term =
      interaction_row$term,
    
    estimate =
      NA_real_,
    
    std_error =
      NA_real_,
    
    conf_low =
      NA_real_,
    
    conf_high =
      NA_real_,
    
    R2 =
      interaction_row$R2,
    
    F =
      interaction_row$F,
    
    statistic =
      NA_real_,
    
    p_value =
      interaction_row$`Pr(>F)`
  )
  
  formulas <- tibble(
    body_site = "Maternal gut",
    analysis = "Beta diversity",
    metric = metric_label,
    method = "PERMANOVA",
    n = nsamples(
      phy_cc
    ),
    
    formula =
      "study_group * race_ethnicity_group + bmi + antibiotics_during_pregnancy"
  )
  
  list(
    interaction = interaction_result,
    full_permanova = permanova_table,
    formulas = formulas
  )
}

# Run beta-diversity analyses

beta_objects <- list(
  
  run_beta_interaction(
    phy = phy,
    distance_method = "bray"
  ),
  
  run_beta_interaction(
    phy = phy,
    distance_method = "jaccard"
  )
)

beta_interaction_results <- bind_rows(
  lapply(
    beta_objects,
    `[[`,
    "interaction"
  )
) %>%
  mutate(
    q_value = p.adjust(
      p_value,
      method = "BH"
    )
  )

beta_full_permanova <- bind_rows(
  lapply(
    beta_objects,
    `[[`,
    "full_permanova"
  )
)

beta_formulas <- bind_rows(
  lapply(
    beta_objects,
    `[[`,
    "formulas"
  )
)

# Combined table for supplementary reporting

interaction_table <- bind_rows(
  alpha_interaction_results,
  beta_interaction_results
) %>%
  arrange(
    analysis,
    metric
  )

# Save results

write_csv(
  alpha_interaction_results,
  file.path(
    out_dir,
    "gut_alpha_PE_race_ethnicity_interaction_tests.csv"
  )
)

write_csv(
  alpha_simple_effects,
  file.path(
    out_dir,
    "gut_alpha_PE_race_ethnicity_simple_effects.csv"
  )
)

write_csv(
  alpha_diagnostics,
  file.path(
    out_dir,
    "gut_alpha_PE_race_ethnicity_interaction_diagnostics.csv"
  )
)

write_csv(
  beta_interaction_results,
  file.path(
    out_dir,
    "gut_beta_PE_race_ethnicity_interaction_tests.csv"
  )
)

write_csv(
  beta_full_permanova,
  file.path(
    out_dir,
    "gut_beta_PE_race_ethnicity_interaction_full_PERMANOVA.csv"
  )
)

write_csv(
  bind_rows(
    alpha_formulas,
    beta_formulas
  ),
  file.path(
    out_dir,
    "gut_PE_race_ethnicity_interaction_model_formulas.csv"
  )
)

write_csv(
  interaction_table,
  file.path(
    out_dir,
    "Table_S13_gut_PE_race_ethnicity_interaction_tests.csv"
  )
)

# Print results

cat(
  "\nGut alpha-diversity interaction tests\n"
)

print(
  alpha_interaction_results
)

cat(
  "\nGut alpha-diversity simple effects\n"
)

print(
  alpha_simple_effects
)

cat(
  "\nGut alpha residual diagnostics\n"
)

print(
  alpha_diagnostics
)

cat(
  "\nGut beta-diversity interaction tests\n"
)

print(
  beta_interaction_results
)

cat(
  "\nFiles generated:\n"
)

print(
  list.files(
    out_dir,
    pattern = "gut|Table_S13_gut",
    full.names = TRUE
  )
)
#mom-infant-pairing diagnostic 16S-WGS####
# Diagnostic check before maternal stool–meconium similarity analysis

library(phyloseq)
library(dplyr)
library(tibble)
library(stringr)

meconium_path <- file.path(PROJECT_DIR, "P16S", "meconium", "phy_obj", "PHY_fil_rar_alpha_meco.RData")

maternal_stool_path <- file.path(PROJECT_DIR, "WGS", "obj", "phy_mtph_gtdb.rds")

output_dir <- file.path(PROJECT_DIR, "P16S", "results", "maternal_stool_meconium_similarity", "diagnostic")

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# Load meconium phyloseq object from RData

meco_env <- new.env()

load(
  meconium_path,
  envir = meco_env
)

meco_objects <- ls(
  meco_env
)

cat("\nObjects found in meconium RData:\n")
print(meco_objects)

meco_phyloseq_objects <- meco_objects[
  vapply(
    meco_objects,
    function(x) {
      inherits(
        get(
          x,
          envir = meco_env
        ),
        "phyloseq"
      )
    },
    logical(1)
  )
]

cat("\nPhyloseq objects found in meconium RData:\n")
print(meco_phyloseq_objects)

if (
  length(
    meco_phyloseq_objects
  ) != 1
) {
  
  stop(
    paste(
      "Expected exactly one phyloseq object in the meconium RData.",
      "Found:",
      paste(
        meco_phyloseq_objects,
        collapse = ", "
      ),
      "\nPlease send me this output before proceeding."
    )
  )
}

phy_meco <- get(
  meco_phyloseq_objects[1],
  envir = meco_env
)

# Load maternal stool WGS phyloseq object

phy_mom_stool <- readRDS(
  maternal_stool_path
)

if (
  !inherits(
    phy_mom_stool,
    "phyloseq"
  )
) {
  
  stop(
    "The maternal stool RDS object is not a phyloseq object."
  )
}

# Basic object summary

summarize_phyloseq <- function(
    phy,
    label) {
  
  cat(
    "\n",
    label,
    "\n",
    sep = ""
  )
  
  cat(
    "Samples:",
    nsamples(
      phy
    ),
    "\n"
  )
  
  cat(
    "Taxa:",
    ntaxa(
      phy
    ),
    "\n"
  )
  
  cat(
    "Taxonomic ranks:\n"
  )
  
  print(
    rank_names(
      phy
    )
  )
  
  cat(
    "Sample-data variables:\n"
  )
  
  print(
    sample_variables(
      phy
    )
  )
}

summarize_phyloseq(
  phy_meco,
  "Meconium 16S object"
)

summarize_phyloseq(
  phy_mom_stool,
  "Maternal stool WGS object"
)

# Inspect taxonomy tables

inspect_taxonomy <- function(
    phy,
    label) {
  
  tt <- as.data.frame(
    tax_table(
      phy
    )
  )
  
  tt <- tibble::rownames_to_column(
    tt,
    "feature_id"
  )
  
  cat(
    "\nTaxonomy columns for ",
    label,
    ":\n",
    sep = ""
  )
  
  print(
    colnames(
      tt
    )
  )
  
  cat(
    "\nFirst taxonomy rows for ",
    label,
    ":\n",
    sep = ""
  )
  
  print(
    head(
      tt,
      20
    )
  )
  
  rank_summary <- tibble(
    rank = setdiff(
      colnames(
        tt
      ),
      "feature_id"
    ),
    n_nonmissing = vapply(
      tt[
        setdiff(
          colnames(
            tt
          ),
          "feature_id"
        )
      ],
      function(x) {
        sum(
          !is.na(
            x
          ) &
            trimws(
              as.character(
                x
              )
            ) != ""
        )
      },
      integer(1)
    ),
    n_unique_nonmissing = vapply(
      tt[
        setdiff(
          colnames(
            tt
          ),
          "feature_id"
        )
      ],
      function(x) {
        length(
          unique(
            as.character(
              x[
                !is.na(
                  x
                ) &
                  trimws(
                    as.character(
                      x
                    )
                  ) != ""
              ]
            )
          )
        )
      },
      integer(1)
    )
  )
  
  write.csv(
    rank_summary,
    file.path(
      output_dir,
      paste0(
        label,
        "_taxonomy_rank_summary.csv"
      )
    ),
    row.names = FALSE
  )
  
  list(
    taxonomy = tt,
    rank_summary = rank_summary
  )
}

meco_tax <- inspect_taxonomy(
  phy_meco,
  "meconium"
)

mom_tax <- inspect_taxonomy(
  phy_mom_stool,
  "maternal_stool"
)

# Detect genus column

find_genus_column <- function(
    taxonomy_df,
    label) {
  
  tax_cols <- setdiff(
    colnames(
      taxonomy_df
    ),
    "feature_id"
  )
  
  exact_match <- tax_cols[
    tolower(
      tax_cols
    ) == "genus"
  ]
  
  if (
    length(
      exact_match
    ) == 1
  ) {
    
    cat(
      "\nGenus column detected for ",
      label,
      ": ",
      exact_match,
      "\n",
      sep = ""
    )
    
    return(
      exact_match
    )
  }
  
  partial_match <- tax_cols[
    grepl(
      "genus",
      tax_cols,
      ignore.case = TRUE
    )
  ]
  
  if (
    length(
      partial_match
    ) > 0
  ) {
    
    cat(
      "\nPossible genus column(s) for ",
      label,
      ": ",
      paste(
        partial_match,
        collapse = ", "
      ),
      "\n",
      sep = ""
    )
    
    return(
      NA_character_
    )
  }
  
  cat(
    "\nNo column explicitly named Genus was detected for ",
    label,
    ".\n",
    sep = ""
  )
  
  return(
    NA_character_
  )
}

meco_genus_col <- find_genus_column(
  meco_tax$taxonomy,
  "meconium"
)

mom_genus_col <- find_genus_column(
  mom_tax$taxonomy,
  "maternal stool"
)

# Genus-name harmonization for diagnostic overlap
#
# The meconium taxonomy contains Greengenes2-style genus labels such as
# g__Bifidobacterium_388775, whereas the maternal WGS GTDB object may contain
# Bifidobacterium. We therefore report overlap using two levels of cleaning:
#
# 1. Conservative: remove only rank prefixes such as g__
# 2. GG2-compatible: additionally remove terminal numeric identifiers
#    (for example, Bifidobacterium_388775 -> Bifidobacterium)
#
# Letter suffixes used by GTDB to distinguish genera (for example Blautia_A,
# Ruminococcus_E, Campylobacter_B) are retained.

clean_genus_conservative <- function(x) {
  
  x <- as.character(x)
  x <- trimws(x)
  
  x <- sub(
    "^g__",
    "",
    x
  )
  
  x <- sub(
    "^g_",
    "",
    x
  )
  
  x <- sub(
    "^Genus__",
    "",
    x,
    ignore.case = TRUE
  )
  
  x[
    is.na(x) |
      x == "" |
      tolower(x) %in% c(
        "na",
        "nan",
        "unknown",
        "unclassified",
        "uncultured"
      )
  ] <- NA_character_
  
  x
}

clean_genus_gg2_compatible <- function(x) {
  
  x <- clean_genus_conservative(x)
  
  # Remove terminal Greengenes2 numeric identifiers only.
  # Keep GTDB letter suffixes such as _A, _B, _E, etc.
  x <- sub(
    "_[0-9]+$",
    "",
    x
  )
  
  x[
    is.na(x) |
      x == "" |
      tolower(x) %in% c(
        "na",
        "nan",
        "unknown",
        "unclassified",
        "uncultured"
      )
  ] <- NA_character_
  
  x
}

if (
  !is.na(meco_genus_col) &&
  !is.na(mom_genus_col)
) {
  
  meco_genera_raw <- unique(
    as.character(
      meco_tax$taxonomy[[meco_genus_col]]
    )
  )
  
  mom_genera_raw <- unique(
    as.character(
      mom_tax$taxonomy[[mom_genus_col]]
    )
  )
  
  meco_cons <- unique(
    na.omit(
      clean_genus_conservative(
        meco_genera_raw
      )
    )
  )
  
  mom_cons <- unique(
    na.omit(
      clean_genus_conservative(
        mom_genera_raw
      )
    )
  )
  
  meco_gg2 <- unique(
    na.omit(
      clean_genus_gg2_compatible(
        meco_genera_raw
      )
    )
  )
  
  mom_gg2 <- unique(
    na.omit(
      clean_genus_gg2_compatible(
        mom_genera_raw
      )
    )
  )
  
  shared_exact <- intersect(
    na.omit(meco_genera_raw),
    na.omit(mom_genera_raw)
  )
  
  shared_cons <- intersect(
    meco_cons,
    mom_cons
  )
  
  shared_gg2 <- intersect(
    meco_gg2,
    mom_gg2
  )
  
  genus_overlap_summary <- tibble(
    meconium_unique_genera_conservative =
      length(meco_cons),
    maternal_stool_unique_genera_conservative =
      length(mom_cons),
    shared_exact_genera =
      length(shared_exact),
    shared_after_prefix_cleaning =
      length(shared_cons),
    meconium_unique_genera_GG2_normalized =
      length(meco_gg2),
    maternal_stool_unique_genera_GG2_normalized =
      length(mom_gg2),
    shared_after_GG2_numeric_suffix_cleaning =
      length(shared_gg2),
    percent_meconium_genera_shared_GG2_normalized =
      100 *
      length(shared_gg2) /
      length(meco_gg2),
    percent_maternal_genera_shared_GG2_normalized =
      100 *
      length(shared_gg2) /
      length(mom_gg2)
  )
  
  cat(
    "\nGenus overlap summary:\n"
  )
  
  print(
    genus_overlap_summary
  )
  
  cat(
    "\nFirst shared genera after conservative cleaning:\n"
  )
  
  print(
    head(
      sort(shared_cons),
      100
    )
  )
  
  cat(
    "\nFirst shared genera after GG2-compatible numeric-suffix cleaning:\n"
  )
  
  print(
    head(
      sort(shared_gg2),
      100
    )
  )
  
  write.csv(
    genus_overlap_summary,
    file.path(
      output_dir,
      "genus_overlap_summary.csv"
    ),
    row.names = FALSE
  )
  
  write.csv(
    tibble(
      shared_genus =
        sort(shared_gg2)
    ),
    file.path(
      output_dir,
      "shared_genera.csv"
    ),
    row.names = FALSE
  )
  
  write.csv(
    tibble(
      meconium_only_genus =
        sort(
          setdiff(
            meco_gg2,
            mom_gg2
          )
        )
    ),
    file.path(
      output_dir,
      "meconium_only_genera.csv"
    ),
    row.names = FALSE
  )
  
  write.csv(
    tibble(
      maternal_stool_only_genus =
        sort(
          setdiff(
            mom_gg2,
            meco_gg2
          )
        )
    ),
    file.path(
      output_dir,
      "maternal_stool_only_genera.csv"
    ),
    row.names = FALSE
  )
  
  genus_mapping_preview <- tibble(
    meconium_raw = meco_genera_raw,
    meconium_prefix_cleaned =
      clean_genus_conservative(
        meco_genera_raw
      ),
    meconium_GG2_normalized =
      clean_genus_gg2_compatible(
        meco_genera_raw
      )
  )
  
  write.csv(
    genus_mapping_preview,
    file.path(
      output_dir,
      "meconium_genus_cleaning_preview.csv"
    ),
    row.names = FALSE
  )
}

# Inspect participant matching

meco_meta <- as.data.frame(
  sample_data(
    phy_meco
  )
) %>%
  tibble::rownames_to_column(
    "sample_name"
  )

mom_meta <- as.data.frame(
  sample_data(
    phy_mom_stool
  )
) %>%
  tibble::rownames_to_column(
    "sample_name"
  )

candidate_id_vars <- c(
  "study_id",
  "Study_ID",
  "studyid",
  "participant_id",
  "subject_id",
  "subject"
)

find_id_variable <- function(
    meta,
    label) {
  
  found <- candidate_id_vars[
    candidate_id_vars %in%
      colnames(
        meta
      )
  ]
  
  cat(
    "\nCandidate participant ID variable(s) for ",
    label,
    ": ",
    paste(
      found,
      collapse = ", "
    ),
    "\n",
    sep = ""
  )
  
  if (
    length(
      found
    ) == 1
  ) {
    
    return(
      found
    )
  }
  
  return(
    NA_character_
  )
}

meco_id_var <- find_id_variable(
  meco_meta,
  "meconium"
)

mom_id_var <- find_id_variable(
  mom_meta,
  "maternal stool"
)

if (
  !is.na(
    meco_id_var
  ) &&
  !is.na(
    mom_id_var
  )
) {
  
  meco_ids <- unique(
    as.character(
      meco_meta[[meco_id_var]]
    )
  )
  
  mom_ids <- unique(
    as.character(
      mom_meta[[mom_id_var]]
    )
  )
  
  meco_ids <- meco_ids[
    !is.na(
      meco_ids
    ) &
      meco_ids != ""
  ]
  
  mom_ids <- mom_ids[
    !is.na(
      mom_ids
    ) &
      mom_ids != ""
  ]
  
  shared_ids <- intersect(
    meco_ids,
    mom_ids
  )
  
  id_overlap_summary <- tibble(
    meconium_participants =
      length(
        meco_ids
      ),
    maternal_stool_participants =
      length(
        mom_ids
      ),
    shared_mother_infant_ids =
      length(
        shared_ids
      )
  )
  
  cat(
    "\nParticipant overlap summary:\n"
  )
  
  print(
    id_overlap_summary
  )
  
  cat(
    "\nMeconium IDs without maternal stool:\n"
  )
  
  print(
    setdiff(
      meco_ids,
      mom_ids
    )
  )
  
  cat(
    "\nMaternal stool IDs without meconium:\n"
  )
  
  print(
    setdiff(
      mom_ids,
      meco_ids
    )
  )
  
  write.csv(
    id_overlap_summary,
    file.path(
      output_dir,
      "participant_overlap_summary.csv"
    ),
    row.names = FALSE
  )
}

# Check variables needed for the adjusted Control vs preeclampsia model

required_model_vars <- c(
  "study_group",
  "days_stool1",
  "mod",
  "antibiotics_during_pregnancy",
  "gestational_age"
)

variable_presence <- bind_rows(
  lapply(
    required_model_vars,
    function(v) {
      
      tibble(
        variable = v,
        present_in_meconium =
          v %in%
          colnames(
            meco_meta
          ),
        present_in_maternal_stool =
          v %in%
          colnames(
            mom_meta
          ),
        meconium_nonmissing =
          if (
            v %in%
            colnames(
              meco_meta
            )
          ) {
            sum(
              !is.na(
                meco_meta[[v]]
              )
            )
          } else {
            NA_integer_
          },
        maternal_stool_nonmissing =
          if (
            v %in%
            colnames(
              mom_meta
            )
          ) {
            sum(
              !is.na(
                mom_meta[[v]]
              )
            )
          } else {
            NA_integer_
          }
      )
    }
  )
)

cat(
  "\nVariables needed for adjusted similarity model:\n"
)

print(
  variable_presence
)

write.csv(
  variable_presence,
  file.path(
    output_dir,
    "required_model_variable_presence.csv"
  ),
  row.names = FALSE
)

cat(
  "\nDiagnostic files saved to:\n",
  output_dir,
  "\n"
)

cat(
  "\nPlease send me the console output plus these files if they were created:\n",
  "genus_overlap_summary.csv\n",
  "shared_genera.csv\n",
  "participant_overlap_summary.csv\n",
  "required_model_variable_presence.csv\n"
)


#PAIRING analysis####
# Maternal stool–meconium microbiome resemblance
#
# Maternal stool: shotgun metagenomics / GTDB
# Infant meconium: 16S rRNA gene sequencing / Greengenes2
#
# Analyses:
# 1. True mother–infant pairs vs randomized pairs
# 2. Direct comparison of true-pair similarity between Control and preeclampsia
# 3. Same adjusted Control vs preeclampsia model restricted to term deliveries
#
# Analyses are performed at the genus level using Bray-Curtis and binary Jaccard
# similarity. Genus labels are harmonized across Greengenes2 and GTDB before
# identifying shared genera.

library(phyloseq)
library(vegan)
library(dplyr)
library(tidyr)
library(readr)
library(MASS)
library(broom)

# Paths

meco_file <- file.path(PROJECT_DIR, "P16S", "meconium", "phy_obj", "PHY_fil_rar_alpha_meco.RData")

maternal_stool_file <- file.path(PROJECT_DIR, "WGS", "obj", "phy_mtph_gtdb.rds")

p16s_project_dir <- file.path(PROJECT_DIR, "P16S")

analysis_date <- as.character(
  Sys.Date()
)

outdir <- file.path(
  p16s_project_dir,
  "results",
  "maternal_stool_meconium_similarity",
  analysis_date
)

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

# Load meconium phyloseq

meco_env <- new.env()

load(
  meco_file,
  envir = meco_env
)

meco_objects <- ls(
  meco_env
)

meco_phyloseq_objects <- meco_objects[
  vapply(
    meco_objects,
    function(x) {
      inherits(
        get(
          x,
          envir = meco_env
        ),
        "phyloseq"
      )
    },
    logical(1)
  )
]

if (
  length(
    meco_phyloseq_objects
  ) != 1
) {
  stop(
    paste(
      "Expected exactly one phyloseq object in the meconium RData. Found:",
      paste(
        meco_phyloseq_objects,
        collapse = ", "
      )
    )
  )
}

phy_meco <- get(
  meco_phyloseq_objects[1],
  envir = meco_env
)

# Load maternal stool phyloseq

phy_mom_stool <- readRDS(
  maternal_stool_file
)

if (
  !inherits(
    phy_mom_stool,
    "phyloseq"
  )
) {
  stop(
    "The maternal stool RDS object is not a phyloseq object."
  )
}

# Remove empty samples and taxa

clean_phyloseq <- function(phy) {
  
  phy <- prune_samples(
    sample_sums(phy) > 0,
    phy
  )
  
  phy <- prune_taxa(
    taxa_sums(phy) > 0,
    phy
  )
  
  phy
}

phy_meco <- clean_phyloseq(
  phy_meco
)

phy_mom_stool <- clean_phyloseq(
  phy_mom_stool
)

cat(
  "\nMeconium samples:",
  nsamples(
    phy_meco
  ),
  "\n"
)

cat(
  "Maternal stool samples:",
  nsamples(
    phy_mom_stool
  ),
  "\n"
)

# Harmonize genus labels
#
# Greengenes2 genus labels can contain a g__ prefix and terminal numeric
# identifiers, for example g__Bifidobacterium_388775.
#
# The maternal WGS GTDB object uses labels such as Bifidobacterium.
#
# We remove:
# - taxonomic rank prefixes such as g__
# - terminal numeric identifiers preceded by an underscore
#
# We retain GTDB letter suffixes such as Blautia_A and Ruminococcus_E.

clean_genus_name <- function(x) {
  
  x <- as.character(
    x
  )
  
  x <- trimws(
    x
  )
  
  x <- sub(
    "^g__",
    "",
    x
  )
  
  x <- sub(
    "^g_",
    "",
    x
  )
  
  x <- sub(
    "^Genus__",
    "",
    x,
    ignore.case = TRUE
  )
  
  x <- sub(
    "_[0-9]+$",
    "",
    x
  )
  
  x[
    is.na(
      x
    ) |
      x == "" |
      tolower(
        x
      ) %in% c(
        "na",
        "nan",
        "unknown",
        "unclassified",
        "uncultured"
      )
  ] <- NA_character_
  
  x
}

# Prepare genus-level abundance matrix

prepare_genus_data <- function(
    phy,
    dataset_name) {
  
  if (
    !"Genus" %in%
    rank_names(
      phy
    )
  ) {
    stop(
      paste(
        "Genus rank is not available in",
        dataset_name
      )
    )
  }
  
  phy_genus <- tax_glom(
    phy,
    taxrank = "Genus",
    NArm = FALSE
  )
  
  tax_df <- as.data.frame(
    tax_table(
      phy_genus
    ),
    stringsAsFactors = FALSE
  )
  
  genus_raw <- as.character(
    tax_df$Genus
  )
  
  genus_clean <- clean_genus_name(
    genus_raw
  )
  
  keep <- !is.na(
    genus_clean
  )
  
  phy_genus <- prune_taxa(
    keep,
    phy_genus
  )
  
  genus_clean <- genus_clean[
    keep
  ]
  
  mat <- as(
    otu_table(
      phy_genus
    ),
    "matrix"
  )
  
  if (
    taxa_are_rows(
      phy_genus
    )
  ) {
    mat <- t(
      mat
    )
  }
  
  colnames(
    mat
  ) <- genus_clean
  
  # Collapse genera that become identical after taxonomy harmonization
  
  if (
    anyDuplicated(
      colnames(
        mat
      )
    )
  ) {
    mat <- t(
      rowsum(
        t(
          mat
        ),
        group = colnames(
          mat
        ),
        reorder = FALSE
      )
    )
  }
  
  # Convert to relative abundance within each original dataset
  
  mat <- mat /
    rowSums(
      mat
    )
  
  meta <- data.frame(
    sample_data(
      phy_genus
    ),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  # Ensure metadata is a plain data.frame. phyloseq::sample_data is an
  # S4 class derived from data.frame and can otherwise propagate through
  # single-column subsetting into dplyr.
  meta <- as.data.frame(
    meta,
    stringsAsFactors = FALSE
  )
  
  meta$sample_name <- rownames(
    meta
  )
  
  if (
    !"study_id" %in%
    colnames(
      meta
    )
  ) {
    stop(
      paste(
        "study_id is missing from",
        dataset_name
      )
    )
  }
  
  meta$study_id <- as.character(
    meta$study_id
  )
  
  if (
    anyDuplicated(
      meta$study_id
    )
  ) {
    duplicate_ids <- unique(
      meta$study_id[
        duplicated(
          meta$study_id
        )
      ]
    )
    
    stop(
      paste(
        "Duplicated study_id values in",
        dataset_name,
        ":",
        paste(
          duplicate_ids,
          collapse = ", "
        )
      )
    )
  }
  
  mat <- mat[
    meta$sample_name,
    ,
    drop = FALSE
  ]
  
  rownames(
    mat
  ) <- meta$study_id
  
  list(
    matrix = mat,
    metadata = meta
  )
}

meco_data <- prepare_genus_data(
  phy_meco,
  "Infant meconium"
)

mom_stool_data <- prepare_genus_data(
  phy_mom_stool,
  "Maternal stool"
)

# Check genus overlap

shared_genera_initial <- intersect(
  colnames(
    meco_data$matrix
  ),
  colnames(
    mom_stool_data$matrix
  )
)

genus_overlap_summary <- tibble::tibble(
  meconium_genera =
    ncol(
      meco_data$matrix
    ),
  maternal_stool_genera =
    ncol(
      mom_stool_data$matrix
    ),
  shared_genera =
    length(
      shared_genera_initial
    ),
  percent_meconium_shared =
    100 *
    length(
      shared_genera_initial
    ) /
    ncol(
      meco_data$matrix
    ),
  percent_maternal_stool_shared =
    100 *
    length(
      shared_genera_initial
    ) /
    ncol(
      mom_stool_data$matrix
    )
)

print(
  genus_overlap_summary
)

readr::write_csv(
  genus_overlap_summary,
  file.path(
    outdir,
    "genus_overlap_summary.csv"
  )
)

readr::write_csv(
  tibble::tibble(
    shared_genus =
      sort(
        shared_genera_initial
      )
  ),
  file.path(
    outdir,
    "shared_genera_used.csv"
  )
)

# Retrieve metadata from infant or maternal datasets
#
# Always return an ordinary vector. This is important because phyloseq
# sample_data objects can retain their S4 class after single-column
# subsetting, which causes dplyr::filter() to fail.

get_metadata_variable <- function(
    ids,
    infant_meta,
    maternal_meta,
    candidates) {
  
  infant_meta <- as.data.frame(
    infant_meta,
    stringsAsFactors = FALSE
  )
  
  maternal_meta <- as.data.frame(
    maternal_meta,
    stringsAsFactors = FALSE
  )
  
  ids <- as.character(
    ids
  )
  
  for (
    var in candidates
  ) {
    
    if (
      var %in%
      colnames(
        infant_meta
      )
    ) {
      
      idx <- match(
        ids,
        as.character(
          infant_meta[["study_id"]]
        )
      )
      
      values <- infant_meta[[var]][
        idx
      ]
      
      # Drop any residual dimensions/classes from sample_data-derived columns.
      values <- as.vector(
        values
      )
      
      if (
        !all(
          is.na(
            values
          )
        )
      ) {
        return(
          values
        )
      }
    }
    
    if (
      var %in%
      colnames(
        maternal_meta
      )
    ) {
      
      idx <- match(
        ids,
        as.character(
          maternal_meta[["study_id"]]
        )
      )
      
      values <- maternal_meta[[var]][
        idx
      ]
      
      values <- as.vector(
        values
      )
      
      if (
        !all(
          is.na(
            values
          )
        )
      ) {
        return(
          values
        )
      }
    }
  }
  
  rep(
    NA,
    length(
      ids
    )
  )
}

# Derangement for randomized mother-infant pairs

make_derangement <- function(x) {
  
  n <- length(
    x
  )
  
  if (
    n < 2
  ) {
    stop(
      "Randomization stratum contains fewer than 2 pairs."
    )
  }
  
  if (
    n == 2
  ) {
    return(
      rev(
        x
      )
    )
  }
  
  for (
    i in 1:1000
  ) {
    
    y <- sample(
      x,
      size = n,
      replace = FALSE
    )
    
    if (
      all(
        y != x
      )
    ) {
      return(
        y
      )
    }
  }
  
  shift <- sample(
    1:(
      n - 1
    ),
    1
  )
  
  x[
    (
      (
        seq_along(
          x
        ) -
          1 +
          shift
      ) %%
        n
    ) +
      1
  ]
}

# Permutation p-value

permutation_p <- function(
    observed,
    null) {
  
  p_upper <- (
    sum(
      null >= observed,
      na.rm = TRUE
    ) +
      1
  ) /
    (
      length(
        null
      ) +
        1
    )
  
  p_lower <- (
    sum(
      null <= observed,
      na.rm = TRUE
    ) +
      1
  ) /
    (
      length(
        null
      ) +
        1
    )
  
  p_two <- min(
    1,
    2 *
      min(
        p_upper,
        p_lower
      )
  )
  
  c(
    p_more_similar =
      p_upper,
    p_two_sided =
      p_two
  )
}

# Extract RLM preeclampsia coefficient

extract_rlm_pe <- function(
    model,
    distance_method,
    population_label) {
  
  ct <- as.data.frame(
    summary(
      model
    )$coefficients
  )
  
  ct$term <- rownames(
    ct
  )
  
  rownames(
    ct
  ) <- NULL
  
  names(
    ct
  )[1:3] <- c(
    "estimate",
    "std_error",
    "statistic"
  )
  
  pe <- ct %>%
    dplyr::filter(
      term ==
        "study_groupPE"
    )
  
  if (
    nrow(
      pe
    ) != 1
  ) {
    stop(
      paste(
        "Could not identify study_groupPE in",
        distance_method,
        population_label
      )
    )
  }
  
  pe %>%
    dplyr::mutate(
      site =
        "Maternal stool",
      distance_method =
        distance_method,
      analysis_population =
        population_label,
      method =
        "RLM",
      n_model =
        length(
          model$residuals
        ),
      conf_low =
        estimate -
        1.96 *
        std_error,
      conf_high =
        estimate +
        1.96 *
        std_error,
      p_value =
        2 *
        pnorm(
          abs(
            statistic
          ),
          lower.tail = FALSE
        )
    ) %>%
    dplyr::select(
      site,
      distance_method,
      analysis_population,
      method,
      n_model,
      term,
      estimate,
      std_error,
      conf_low,
      conf_high,
      statistic,
      p_value
    )
}

# Extract LM preeclampsia coefficient

extract_lm_pe <- function(
    model,
    distance_method,
    population_label) {
  
  broom::tidy(
    model,
    conf.int = TRUE,
    conf.level = 0.95
  ) %>%
    dplyr::filter(
      term ==
        "study_groupPE"
    ) %>%
    dplyr::transmute(
      site =
        "Maternal stool",
      distance_method =
        distance_method,
      analysis_population =
        population_label,
      method =
        "LM",
      n_model =
        stats::nobs(
          model
        ),
      term =
        term,
      estimate =
        estimate,
      std_error =
        std.error,
      conf_low =
        conf.low,
      conf_high =
        conf.high,
      statistic =
        statistic,
      p_value =
        p.value
    )
}

# Main analysis function

run_stool_meconium_similarity <- function(
    maternal_data,
    infant_data,
    distance_method = c(
      "bray",
      "jaccard"
    ),
    B = 5000,
    seed = 20260817) {
  
  distance_method <- match.arg(
    distance_method
  )
  
  set.seed(
    seed
  )
  
  ids <- intersect(
    rownames(
      maternal_data$matrix
    ),
    rownames(
      infant_data$matrix
    )
  )
  
  if (
    length(
      ids
    ) <
    10
  ) {
    stop(
      "Too few matched maternal stool–meconium pairs."
    )
  }
  
  maternal_meta <- maternal_data$metadata
  infant_meta <- infant_data$metadata
  
  pair_meta <- data.frame(
    study_id =
      ids,
    stringsAsFactors = FALSE
  )
  
  pair_meta$study_group <-
    get_metadata_variable(
      ids,
      infant_meta,
      maternal_meta,
      "study_group"
    )
  
  pair_meta$ethnicity <-
    get_metadata_variable(
      ids,
      infant_meta,
      maternal_meta,
      "ethnicity"
    )
  
  pair_meta$gestational_age <-
    get_metadata_variable(
      ids,
      infant_meta,
      maternal_meta,
      c(
        "gestational_age",
        "gestational_age_at_delivery"
      )
    )
  
  pair_meta$days_stool1 <-
    get_metadata_variable(
      ids,
      infant_meta,
      maternal_meta,
      "days_stool1"
    )
  
  pair_meta$mod <-
    get_metadata_variable(
      ids,
      infant_meta,
      maternal_meta,
      c(
        "mod",
        "delivery_mode"
      )
    )
  
  pair_meta$antibiotics_during_pregnancy <-
    get_metadata_variable(
      ids,
      infant_meta,
      maternal_meta,
      c(
        "antibiotics_during_pregnancy",
        "antibiotic_during_pregnancy"
      )
    )
  
  # Defensive coercion: all model metadata must be ordinary atomic vectors.
  pair_meta$study_group <- as.vector(
    pair_meta$study_group
  )
  pair_meta$ethnicity <- as.vector(
    pair_meta$ethnicity
  )
  pair_meta$gestational_age <- as.numeric(
    as.vector(
      pair_meta$gestational_age
    )
  )
  pair_meta$days_stool1 <- as.numeric(
    as.vector(
      pair_meta$days_stool1
    )
  )
  pair_meta$mod <- as.vector(
    pair_meta$mod
  )
  pair_meta$antibiotics_during_pregnancy <- as.vector(
    pair_meta$antibiotics_during_pregnancy
  )
  
  if (
    all(
      is.na(
        pair_meta$study_group
      )
    )
  ) {
    stop(
      "study_group was not found in either phyloseq metadata table."
    )
  }
  
  pair_meta <- pair_meta %>%
    dplyr::filter(
      !is.na(
        study_group
      )
    )
  
  ids <- pair_meta$study_id
  
  # Restrict to genera represented in both sequencing platforms
  
  common_genera <- intersect(
    colnames(
      maternal_data$matrix
    ),
    colnames(
      infant_data$matrix
    )
  )
  
  maternal_matrix <- maternal_data$matrix[
    ids,
    common_genera,
    drop = FALSE
  ]
  
  infant_matrix <- infant_data$matrix[
    ids,
    common_genera,
    drop = FALSE
  ]
  
  # Require each genus to occur at least once in each dataset among matched dyads
  
  keep_genus <-
    colSums(
      maternal_matrix
    ) >
    0 &
    colSums(
      infant_matrix
    ) >
    0
  
  maternal_matrix <- maternal_matrix[
    ,
    keep_genus,
    drop = FALSE
  ]
  
  infant_matrix <- infant_matrix[
    ,
    keep_genus,
    drop = FALSE
  ]
  
  n_shared_genera <- ncol(
    maternal_matrix
  )
  
  # Remove dyads with no abundance among the shared genus set
  
  maternal_sums <- rowSums(
    maternal_matrix
  )
  
  infant_sums <- rowSums(
    infant_matrix
  )
  
  zero_ids <- ids[
    maternal_sums ==
      0 |
      infant_sums ==
      0
  ]
  
  if (
    length(
      zero_ids
    ) >
    0
  ) {
    
    warning(
      paste(
        length(
          zero_ids
        ),
        "dyads had no abundance among shared genera and were removed:",
        paste(
          zero_ids,
          collapse = ", "
        )
      )
    )
    
    ids <- ids[
      !ids %in%
        zero_ids
    ]
    
    pair_meta <- pair_meta[
      match(
        ids,
        pair_meta$study_id
      ),
      ,
      drop = FALSE
    ]
    
    maternal_matrix <- maternal_matrix[
      ids,
      ,
      drop = FALSE
    ]
    
    infant_matrix <- infant_matrix[
      ids,
      ,
      drop = FALSE
    ]
  }
  
  # Renormalize after restriction to the shared genus set
  
  maternal_matrix <-
    maternal_matrix /
    rowSums(
      maternal_matrix
    )
  
  infant_matrix <-
    infant_matrix /
    rowSums(
      infant_matrix
    )
  
  rownames(
    maternal_matrix
  ) <- paste0(
    "M_",
    ids
  )
  
  rownames(
    infant_matrix
  ) <- paste0(
    "I_",
    ids
  )
  
  combined_matrix <- rbind(
    maternal_matrix,
    infant_matrix
  )
  
  if (
    distance_method ==
    "bray"
  ) {
    
    dist_matrix <- as.matrix(
      vegan::vegdist(
        combined_matrix,
        method = "bray"
      )
    )
    
    distance_label <- "Bray-Curtis"
    
  } else {
    
    dist_matrix <- as.matrix(
      vegan::vegdist(
        combined_matrix,
        method = "jaccard",
        binary = TRUE
      )
    )
    
    distance_label <- "Jaccard"
  }
  
  cross_distance <- dist_matrix[
    paste0(
      "M_",
      ids
    ),
    paste0(
      "I_",
      ids
    ),
    drop = FALSE
  ]
  
  rownames(
    cross_distance
  ) <- ids
  
  colnames(
    cross_distance
  ) <- ids
  
  cross_similarity <-
    1 -
    cross_distance
  
  own_similarity <- diag(
    cross_similarity
  )
  
  # Randomization strata
  
  pair_meta$study_group <- factor(
    pair_meta$study_group
  )
  
  if (
    !all(
      is.na(
        pair_meta$ethnicity
      )
    ) &&
    all(
      !is.na(
        pair_meta$ethnicity
      )
    )
  ) {
    
    candidate_strata <- interaction(
      pair_meta$study_group,
      pair_meta$ethnicity,
      drop = TRUE
    )
    
    strata_sizes <- table(
      candidate_strata
    )
    
    if (
      all(
        strata_sizes >=
        2
      )
    ) {
      
      random_strata <- candidate_strata
      strata_definition <- "study_group × ethnicity"
      
    } else {
      
      random_strata <- pair_meta$study_group
      strata_definition <- "study_group"
    }
    
  } else {
    
    random_strata <- pair_meta$study_group
    strata_definition <- "study_group"
  }
  
  if (
    any(
      table(
        random_strata
      ) <
      2
    )
  ) {
    stop(
      "At least one randomization stratum has fewer than 2 participants."
    )
  }
  
  group_levels <- levels(
    droplevels(
      pair_meta$study_group
    )
  )
  
  null_overall <- numeric(
    B
  )
  
  null_group <- matrix(
    NA_real_,
    nrow = B,
    ncol =
      length(
        group_levels
      )
  )
  
  colnames(
    null_group
  ) <- group_levels
  
  for (
    b in seq_len(
      B
    )
  ) {
    
    random_infant <- character(
      length(
        ids
      )
    )
    
    for (
      stratum in unique(
        as.character(
          random_strata
        )
      )
    ) {
      
      idx <- which(
        as.character(
          random_strata
        ) ==
          stratum
      )
      
      random_infant[
        idx
      ] <- make_derangement(
        ids[
          idx
        ]
      )
    }
    
    random_similarity <- cross_similarity[
      cbind(
        ids,
        random_infant
      )
    ]
    
    null_overall[
      b
    ] <- mean(
      random_similarity,
      na.rm = TRUE
    )
    
    for (
      g in group_levels
    ) {
      
      idx_g <- pair_meta$study_group ==
        g
      
      null_group[
        b,
        g
      ] <- mean(
        random_similarity[
          idx_g
        ],
        na.rm = TRUE
      )
    }
  }
  
  observed_overall <- mean(
    own_similarity,
    na.rm = TRUE
  )
  
  observed_group <- sapply(
    group_levels,
    function(g) {
      mean(
        own_similarity[
          pair_meta$study_group ==
            g
        ],
        na.rm = TRUE
      )
    }
  )
  
  p_overall <- permutation_p(
    observed_overall,
    null_overall
  )
  
  real_vs_random <- data.frame(
    site =
      "Maternal stool",
    distance_method =
      distance_method,
    comparison =
      "All",
    n_pairs =
      length(
        ids
      ),
    n_shared_genera =
      n_shared_genera,
    observed_real_similarity =
      observed_overall,
    mean_random_similarity =
      mean(
        null_overall
      ),
    real_minus_random =
      observed_overall -
      mean(
        null_overall
      ),
    p_more_similar =
      p_overall[
        "p_more_similar"
      ],
    p_two_sided =
      p_overall[
        "p_two_sided"
      ],
    randomization_strata =
      strata_definition,
    stringsAsFactors = FALSE
  )
  
  for (
    g in group_levels
  ) {
    
    pg <- permutation_p(
      observed_group[
        g
      ],
      null_group[
        ,
        g
      ]
    )
    
    row_g <- data.frame(
      site =
        "Maternal stool",
      distance_method =
        distance_method,
      comparison =
        g,
      n_pairs =
        sum(
          pair_meta$study_group ==
            g
        ),
      n_shared_genera =
        n_shared_genera,
      observed_real_similarity =
        observed_group[
          g
        ],
      mean_random_similarity =
        mean(
          null_group[
            ,
            g
          ],
          na.rm = TRUE
        ),
      real_minus_random =
        observed_group[
          g
        ] -
        mean(
          null_group[
            ,
            g
          ],
          na.rm = TRUE
        ),
      p_more_similar =
        pg[
          "p_more_similar"
        ],
      p_two_sided =
        pg[
          "p_two_sided"
        ],
      randomization_strata =
        strata_definition,
      stringsAsFactors = FALSE
    )
    
    real_vs_random <- dplyr::bind_rows(
      real_vs_random,
      row_g
    )
  }
  
  # Pair-level data for direct Control vs preeclampsia comparison
  
  pair_results <- pair_meta
  
  pair_results$site <-
    "Maternal stool"
  
  pair_results$distance_method <-
    distance_label
  
  pair_results$own_similarity <-
    own_similarity
  
  # Prepare model variables
  
  pair_results$study_group <- factor(
    pair_results$study_group
  )
  
  if (
    "Control" %in%
    levels(
      pair_results$study_group
    )
  ) {
    pair_results$study_group <- relevel(
      pair_results$study_group,
      ref = "Control"
    )
  }
  
  pair_results$mod <- factor(
    pair_results$mod
  )
  
  pair_results$antibiotics_during_pregnancy <- factor(
    pair_results$antibiotics_during_pregnancy
  )
  
  if (
    "No" %in%
    levels(
      pair_results$antibiotics_during_pregnancy
    )
  ) {
    pair_results$antibiotics_during_pregnancy <- relevel(
      pair_results$antibiotics_during_pregnancy,
      ref = "No"
    )
  }
  
  model_vars <- c(
    "study_id",
    "own_similarity",
    "study_group",
    "days_stool1",
    "mod",
    "antibiotics_during_pregnancy",
    "gestational_age"
  )
  
  model_data_cc <- pair_results %>%
    dplyr::select(
      dplyr::all_of(
        model_vars
      )
    ) %>%
    tidyr::drop_na() %>%
    droplevels()
  
  model_data_term <- model_data_cc %>%
    dplyr::filter(
      gestational_age >=
        37
    ) %>%
    droplevels()
  
  if (
    length(
      unique(
        model_data_cc$study_group
      )
    ) <
    2
  ) {
    stop(
      "Both Control and PE are required for the adjusted all-delivery model."
    )
  }
  
  if (
    length(
      unique(
        model_data_term$study_group
      )
    ) <
    2
  ) {
    stop(
      "Both Control and PE are required for the adjusted term-only model."
    )
  }
  
  adjusted_formula <-
    own_similarity ~
    study_group +
    days_stool1 +
    mod +
    antibiotics_during_pregnancy +
    gestational_age
  
  # RLM primary analysis
  
  rlm_all <- MASS::rlm(
    adjusted_formula,
    data = model_data_cc,
    maxit = 100
  )
  
  rlm_term <- MASS::rlm(
    adjusted_formula,
    data = model_data_term,
    maxit = 100
  )
  
  # LM sensitivity analysis
  
  lm_all <- lm(
    adjusted_formula,
    data = model_data_cc
  )
  
  lm_term <- lm(
    adjusted_formula,
    data = model_data_term
  )
  
  rlm_results <- dplyr::bind_rows(
    extract_rlm_pe(
      rlm_all,
      distance_label,
      "All_deliveries"
    ),
    extract_rlm_pe(
      rlm_term,
      distance_label,
      "Term_only"
    )
  )
  
  lm_results <- dplyr::bind_rows(
    extract_lm_pe(
      lm_all,
      distance_label,
      "All_deliveries"
    ),
    extract_lm_pe(
      lm_term,
      distance_label,
      "Term_only"
    )
  )
  
  descriptives <- dplyr::bind_rows(
    model_data_cc %>%
      dplyr::group_by(
        study_group
      ) %>%
      dplyr::summarise(
        n =
          dplyr::n(),
        mean_similarity =
          mean(
            own_similarity,
            na.rm = TRUE
          ),
        sd_similarity =
          sd(
            own_similarity,
            na.rm = TRUE
          ),
        median_similarity =
          median(
            own_similarity,
            na.rm = TRUE
          ),
        IQR_similarity =
          IQR(
            own_similarity,
            na.rm = TRUE
          ),
        .groups =
          "drop"
      ) %>%
      dplyr::mutate(
        site =
          "Maternal stool",
        distance_method =
          distance_label,
        analysis_population =
          "All_deliveries"
      ),
    model_data_term %>%
      dplyr::group_by(
        study_group
      ) %>%
      dplyr::summarise(
        n =
          dplyr::n(),
        mean_similarity =
          mean(
            own_similarity,
            na.rm = TRUE
          ),
        sd_similarity =
          sd(
            own_similarity,
            na.rm = TRUE
          ),
        median_similarity =
          median(
            own_similarity,
            na.rm = TRUE
          ),
        IQR_similarity =
          IQR(
            own_similarity,
            na.rm = TRUE
          ),
        .groups =
          "drop"
      ) %>%
      dplyr::mutate(
        site =
          "Maternal stool",
        distance_method =
          distance_label,
        analysis_population =
          "Term_only"
      )
  ) %>%
    dplyr::select(
      site,
      distance_method,
      analysis_population,
      study_group,
      n,
      mean_similarity,
      sd_similarity,
      median_similarity,
      IQR_similarity
    )
  
  means_wide <- descriptives %>%
    dplyr::select(
      site,
      distance_method,
      analysis_population,
      study_group,
      mean_similarity
    ) %>%
    tidyr::pivot_wider(
      names_from =
        study_group,
      values_from =
        mean_similarity,
      names_prefix =
        "mean_"
    )
  
  rlm_results <- rlm_results %>%
    dplyr::left_join(
      means_wide,
      by = c(
        "site",
        "distance_method",
        "analysis_population"
      )
    ) %>%
    dplyr::mutate(
      interpretation =
        dplyr::case_when(
          estimate <
            0 ~
            "PE dyads have lower adjusted own-pair similarity than Controls",
          estimate >
            0 ~
            "PE dyads have higher adjusted own-pair similarity than Controls",
          TRUE ~
            "No direction"
        )
    )
  
  lm_results <- lm_results %>%
    dplyr::left_join(
      means_wide,
      by = c(
        "site",
        "distance_method",
        "analysis_population"
      )
    )
  
  diagnostics <- tibble::tibble(
    site =
      "Maternal stool",
    distance_method =
      distance_label,
    analysis_population =
      c(
        "All_deliveries",
        "Term_only"
      ),
    n_model =
      c(
        nrow(
          model_data_cc
        ),
        nrow(
          model_data_term
        )
      ),
    LM_residual_shapiro_p =
      c(
        shapiro.test(
          residuals(
            lm_all
          )
        )$p.value,
        shapiro.test(
          residuals(
            lm_term
          )
        )$p.value
      )
  )
  
  cat(
    "\nMaternal stool–meconium:",
    distance_label,
    "\nMatched pairs:",
    length(
      ids
    ),
    "\nShared genera used:",
    n_shared_genera,
    "\nRandomization strata:",
    strata_definition,
    "\n"
  )
  
  list(
    pair_results =
      pair_results,
    real_vs_random =
      real_vs_random,
    rlm_results =
      rlm_results,
    lm_results =
      lm_results,
    descriptives =
      descriptives,
    diagnostics =
      diagnostics,
    shared_genera =
      colnames(
        maternal_matrix
      )
  )
}

# Run Bray-Curtis and Jaccard analyses

res_bray <- run_stool_meconium_similarity(
  maternal_data =
    mom_stool_data,
  infant_data =
    meco_data,
  distance_method =
    "bray",
  B =
    5000,
  seed =
    20260817
)

res_jaccard <- run_stool_meconium_similarity(
  maternal_data =
    mom_stool_data,
  infant_data =
    meco_data,
  distance_method =
    "jaccard",
  B =
    5000,
  seed =
    20260818
)

# Combine stool–meconium results

real_vs_random_stool <- dplyr::bind_rows(
  res_bray$real_vs_random,
  res_jaccard$real_vs_random
)

# Temporary stool-only FDR.
# Final manuscript FDR should be recalculated across all three maternal sites
# (oral, vaginal, stool) and both metrics.

real_vs_random_stool <- real_vs_random_stool %>%
  dplyr::group_by(
    comparison
  ) %>%
  dplyr::mutate(
    q_two_sided_stool_only =
      p.adjust(
        p_two_sided,
        method = "BH"
      )
  ) %>%
  dplyr::ungroup()

rlm_stool <- dplyr::bind_rows(
  res_bray$rlm_results,
  res_jaccard$rlm_results
) %>%
  dplyr::group_by(
    analysis_population
  ) %>%
  dplyr::mutate(
    q_value_stool_only =
      p.adjust(
        p_value,
        method = "BH"
      )
  ) %>%
  dplyr::ungroup()

lm_stool <- dplyr::bind_rows(
  res_bray$lm_results,
  res_jaccard$lm_results
) %>%
  dplyr::group_by(
    analysis_population
  ) %>%
  dplyr::mutate(
    q_value_stool_only =
      p.adjust(
        p_value,
        method = "BH"
      )
  ) %>%
  dplyr::ungroup()

descriptives_stool <- dplyr::bind_rows(
  res_bray$descriptives,
  res_jaccard$descriptives
)

diagnostics_stool <- dplyr::bind_rows(
  res_bray$diagnostics,
  res_jaccard$diagnostics
)

pair_results_stool <- dplyr::bind_rows(
  res_bray$pair_results,
  res_jaccard$pair_results
)

# Save stool–meconium outputs

readr::write_csv(
  real_vs_random_stool,
  file.path(
    outdir,
    "STOOL_MECONIUM_real_vs_random_summary.csv"
  )
)

readr::write_csv(
  rlm_stool,
  file.path(
    outdir,
    "STOOL_MECONIUM_own_similarity_Control_vs_PE_adjusted_RLM.csv"
  )
)

readr::write_csv(
  lm_stool,
  file.path(
    outdir,
    "STOOL_MECONIUM_own_similarity_Control_vs_PE_adjusted_LM_sensitivity.csv"
  )
)

readr::write_csv(
  descriptives_stool,
  file.path(
    outdir,
    "STOOL_MECONIUM_own_similarity_descriptive_Control_vs_PE.csv"
  )
)

readr::write_csv(
  diagnostics_stool,
  file.path(
    outdir,
    "STOOL_MECONIUM_own_similarity_LM_diagnostics.csv"
  )
)

readr::write_csv(
  pair_results_stool,
  file.path(
    outdir,
    "STOOL_MECONIUM_pair_level_similarity.csv"
  )
)

# Recalculate FDR across oral, vaginal, and stool if previous outputs are available

previous_similarity_dir <- file.path(
  p16s_project_dir,
  "results",
  "maternal_infant_similarity",
  analysis_date
)

previous_real_random_file <- file.path(
  previous_similarity_dir,
  "ALL_METRICS_real_vs_random_summary.csv"
)

previous_rlm_file <- file.path(
  previous_similarity_dir,
  "direct_Control_vs_PE_own_similarity",
  "ALL_METRICS_own_similarity_Control_vs_PE_adjusted_RLM.csv"
)

previous_lm_file <- file.path(
  previous_similarity_dir,
  "direct_Control_vs_PE_own_similarity",
  "ALL_METRICS_own_similarity_Control_vs_PE_adjusted_LM_sensitivity.csv"
)

if (
  file.exists(
    previous_real_random_file
  )
) {
  
  previous_rr <- readr::read_csv(
    previous_real_random_file,
    show_col_types = FALSE
  ) %>%
    dplyr::filter(
      site !=
        "Maternal stool"
    ) %>%
    dplyr::select(
      -dplyr::any_of(
        c(
          "q_two_sided",
          "q_two_sided_stool_only"
        )
      )
    )
  
  combined_rr <- dplyr::bind_rows(
    previous_rr,
    real_vs_random_stool %>%
      dplyr::select(
        -dplyr::any_of(
          "q_two_sided_stool_only"
        )
      )
  ) %>%
    dplyr::group_by(
      comparison
    ) %>%
    dplyr::mutate(
      q_two_sided =
        p.adjust(
          p_two_sided,
          method = "BH"
        )
    ) %>%
    dplyr::ungroup()
  
  readr::write_csv(
    combined_rr,
    file.path(
      outdir,
      "ALL_SITES_real_vs_random_summary_FDR6.csv"
    )
  )
  
  cat(
    "\nCombined oral + vaginal + stool real-vs-random FDR was recalculated across 6 site-by-metric tests.\n"
  )
}

if (
  file.exists(
    previous_rlm_file
  )
) {
  
  previous_rlm <- readr::read_csv(
    previous_rlm_file,
    show_col_types = FALSE
  ) %>%
    dplyr::filter(
      site !=
        "Maternal stool"
    ) %>%
    dplyr::select(
      -dplyr::any_of(
        c(
          "q_value",
          "q_value_stool_only"
        )
      )
    )
  
  combined_rlm <- dplyr::bind_rows(
    previous_rlm,
    rlm_stool %>%
      dplyr::select(
        -dplyr::any_of(
          "q_value_stool_only"
        )
      )
  ) %>%
    dplyr::group_by(
      analysis_population
    ) %>%
    dplyr::mutate(
      q_value =
        p.adjust(
          p_value,
          method = "BH"
        )
    ) %>%
    dplyr::ungroup()
  
  readr::write_csv(
    combined_rlm,
    file.path(
      outdir,
      "ALL_SITES_own_similarity_Control_vs_PE_adjusted_RLM_FDR6.csv"
    )
  )
  
  cat(
    "Combined oral + vaginal + stool RLM FDR was recalculated across 6 site-by-metric tests.\n"
  )
}

if (
  file.exists(
    previous_lm_file
  )
) {
  
  previous_lm <- readr::read_csv(
    previous_lm_file,
    show_col_types = FALSE
  ) %>%
    dplyr::filter(
      site !=
        "Maternal stool"
    ) %>%
    dplyr::select(
      -dplyr::any_of(
        c(
          "q_value",
          "q_value_stool_only"
        )
      )
    )
  
  combined_lm <- dplyr::bind_rows(
    previous_lm,
    lm_stool %>%
      dplyr::select(
        -dplyr::any_of(
          "q_value_stool_only"
        )
      )
  ) %>%
    dplyr::group_by(
      analysis_population
    ) %>%
    dplyr::mutate(
      q_value =
        p.adjust(
          p_value,
          method = "BH"
        )
    ) %>%
    dplyr::ungroup()
  
  readr::write_csv(
    combined_lm,
    file.path(
      outdir,
      "ALL_SITES_own_similarity_Control_vs_PE_adjusted_LM_sensitivity_FDR6.csv"
    )
  )
}

# Console summary

cat(
  "\nMaternal stool–meconium true pairs vs randomized pairs\n"
)

print(
  real_vs_random_stool
)

cat(
  "\nAdjusted direct Control vs preeclampsia comparison of true-pair similarity\n"
)

print(
  rlm_stool
)

cat(
  "\nLM sensitivity analysis\n"
)

print(
  lm_stool
)


# Figures for maternal gut alpha diversity by preeclampsia status and ethnicity####
#
# These plots are intended for manual assembly in Illustrator.  Each PDF has
# two panels (Hispanic and White) comparing Control with PE.  Shannon is shown
# on its original scale, whereas the p value is taken from the prespecified
# minimally adjusted Box-Cox model.  Simpson is shown on its original scale.
# The interaction p value tests whether the PE association differs by ethnicity.


if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package ggplot2 is required for the figure-generation section.")
}

figure_dir <- file.path(out_dir, "figures_gut_alpha_ethnicity")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# Harmonize the two ethnicity strata used in the manuscript.
met_plot <- met %>%
  mutate(
    ethnicity_plot = case_when(
      race_ethnicity_group %in% c("Hispanic", "Hispanic/Latina", "Hispanic-Latina") ~ "Hispanic",
      race_ethnicity_group %in% c("Non-Hispanic White/Caucasian", "Not Hispanic/Latina", "White") ~ "White",
      TRUE ~ NA_character_
    ),
    ethnicity_plot = factor(ethnicity_plot, levels = c("Hispanic", "White")),
    study_group = factor(study_group, levels = c("Control", "PE"))
  ) %>%
  filter(!is.na(ethnicity_plot), study_group %in% c("Control", "PE"))

if (nrow(met_plot) == 0 || nlevels(droplevels(met_plot$ethnicity_plot)) < 2) {
  stop("Could not identify both Hispanic and White maternal stool strata.")
}

metric_specs <- list(
  Shannon = list(display = "Shannon diversity", outcome = "shannon_transformed", raw = "sha"),
  Simpson = list(display = "Simpson diversity", outcome = "simp", raw = "simp")
)

interaction_results <- list()

for (metric_name in names(metric_specs)) {
  spec <- metric_specs[[metric_name]]

  # Use the same covariates as the prespecified minimally adjusted model.
  dat <- met_plot %>%
    dplyr::select(all_of(c(spec$outcome, spec$raw, "study_group", "ethnicity_plot", "bmi", "antibiotics_during_pregnancy"))) %>%
    drop_na(all_of(c(spec$outcome, "study_group", "ethnicity_plot", "bmi", "antibiotics_during_pregnancy")))

  dat$study_group <- relevel(factor(dat$study_group), ref = "Control")
  dat$ethnicity_plot <- relevel(factor(dat$ethnicity_plot), ref = "Hispanic")
  dat$antibiotics_during_pregnancy <- factor(dat$antibiotics_during_pregnancy)

  interaction_formula <- as.formula(
    paste(spec$outcome, "~ study_group * ethnicity_plot + bmi + antibiotics_during_pregnancy")
  )
  interaction_model <- if (metric_name == "Simpson") {
    MASS::rlm(interaction_formula, data = dat)
  } else {
    lm(interaction_formula, data = dat)
  }
  interaction_term <- "study_groupPE:ethnicity_plotWhite"
  if (metric_name == "Simpson") {
    coef_tab <- summary(interaction_model)$coefficients
    if (interaction_term %in% rownames(coef_tab)) {
      interaction_statistic <- coef_tab[interaction_term, "t value"]
      interaction_p <- 2 * pnorm(abs(interaction_statistic), lower.tail = FALSE)
    } else {
      interaction_p <- NA_real_
    }
  } else {
    interaction_p <- broom::tidy(interaction_model) %>%
      filter(term == interaction_term) %>%
      pull(p.value)
    if (length(interaction_p) == 0) interaction_p <- NA_real_
  }

  interaction_results[[metric_name]] <- tibble(
    metric = metric_name,
    n = nrow(dat),
    interaction_term = interaction_term,
    interaction_p = interaction_p
  )

  p_label <- if (is.na(interaction_p)) {
    "Interaction p = NA"
  } else {
    paste0("Interaction p = ", format.pval(interaction_p, digits = 3, eps = 0.001))
  }

  # p values for the PE-versus-Control comparison within each ethnicity panel.
  panel_p <- dat %>%
    group_by(ethnicity_plot) %>%
    group_modify(~ {
      m <- if (metric_name == "Simpson") {
        MASS::rlm(as.formula(paste(spec$outcome, "~ study_group + bmi + antibiotics_during_pregnancy")), data = .x)
      } else {
        lm(as.formula(paste(spec$outcome, "~ study_group + bmi + antibiotics_during_pregnancy")), data = .x)
      }
      if (metric_name == "Simpson") {
        tab <- summary(m)$coefficients
        pv <- if ("study_groupPE" %in% rownames(tab)) 2 * pnorm(abs(tab["study_groupPE", "t value"]), lower.tail = FALSE) else NA_real_
      } else {
        pv <- broom::tidy(m) %>% filter(term == "study_groupPE") %>% pull(p.value)
        if (length(pv) == 0) pv <- NA_real_
      }
      tibble(p_value = pv)
    }) %>%
    ungroup() %>%
    mutate(
      x = 1.5,
      p_label = paste0("PE vs Control: p = ", format.pval(p_value, digits = 3, eps = 0.001))
    )

  p <- ggplot2::ggplot(
    dat,
    ggplot2::aes(x = study_group, y = .data[[spec$raw]], fill = study_group)
  ) +
    ggplot2::geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.75, colour = "black") +
    ggplot2::geom_jitter(width = 0.10, height = 0, size = 1.7, alpha = 0.75, colour = "black") +
    ggplot2::geom_text(
      data = panel_p,
      ggplot2::aes(x = x, y = Inf, label = p_label),
      inherit.aes = FALSE,
      vjust = 1.8,
      size = 3.1
    ) +
    ggplot2::facet_wrap(~ ethnicity_plot, nrow = 1) +
    ggplot2::labs(
      title = spec$display,
      subtitle = p_label,
      x = NULL,
      y = spec$display,
      fill = NULL
    ) +
    ggplot2::scale_fill_manual(values = c(Control = "#5a9bd4", PE = "#ff6347")) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      legend.position = "none",
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 10)
    )

  ggplot2::ggsave(
    filename = file.path(figure_dir, paste0("Figure5_", tolower(metric_name), "_PE_vs_Control_by_ethnicity.pdf")),
    plot = p,
    width = 3.5,
    height = 3.8,
    units = "in",
    device = grDevices::cairo_pdf
  )
}

readr::write_csv(
  bind_rows(interaction_results),
  file.path(figure_dir, "Figure5_gut_alpha_ethnicity_interaction_pvalues.csv")
)

cat(
  "\nGut alpha-diversity PDFs saved to:\n",
  figure_dir,
  "\n"
)

cat(
  "\nResults saved to:\n",
  outdir,
  "\n"
)
#end ####
