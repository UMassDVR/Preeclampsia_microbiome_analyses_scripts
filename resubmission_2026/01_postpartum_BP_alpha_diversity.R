#### Blood pressure analysis####
# Postpartum blood pressure analysis across maternal microbiome sites
# Reviewer-focused analysis for oral, vaginal, and maternal gut alpha diversity

options(stringsAsFactors = FALSE)

required_packages <- c(
  "phyloseq",
  "vegan",
  "dplyr",
  "tidyr",
  "readr",
  "tibble",
  "stringr",
  "sandwich",
  "lmtest",
  "emmeans"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Install these packages before running: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(phyloseq)
  library(vegan)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(stringr)
  library(sandwich)
  library(lmtest)
  library(emmeans)
})

# Paths
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

P16S_DIR <- file.path(
  PROJECT_DIR,
  "P16S"
)

METADATA_FILE <- file.path(
  PROJECT_DIR,
  "met",
  "met_clean0.csv"
)

MATERNAL_16S_FILE <- file.path(
  P16S_DIR,
  "R",
  "phy_objects",
  "PHYrar_alpha_ggext_subset02_newseqs.RData"
)

MATERNAL_GUT_FILE <- file.path(
  PROJECT_DIR,
  "WGS",
  "obj",
  "phy_mtph_gtdb.rds"
)

OUT_DIR <- file.path(
  P16S_DIR,
  "results",
  "reviewer_postpartum_BP",
  as.character(Sys.Date())
)

dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

# Codebook-confirmed definitions

codebook_definitions <- tibble::tibble(
  variable = c(
    "bmi",
    "bomeddischarge",
    "ppsystolic",
    "ppdiastolic",
    "ppmed"
  ),
  definition = c(
    "Pre-pregnancy BMI",
    "Taking blood-pressure medication at hospital discharge",
    "Systolic BP at the 6-week postpartum visit",
    "Diastolic BP at the 6-week postpartum visit",
    "Taking blood-pressure medication at the 6-week postpartum visit"
  )
)

readr::write_csv(
  codebook_definitions,
  file.path(
    OUT_DIR,
    "BP_codebook_variable_definitions.csv"
  )
)

# Read and standardize metadata

met <- readr::read_csv(
  METADATA_FILE,
  show_col_types = FALSE,
  name_repair = "minimal"
)

# Remove unnamed/blank columns before any dplyr transformation.
# The metadata export contains an unnamed first column (typically a row index),
# and dplyr::mutate() cannot operate on a data frame with NA or "" column names.
valid_names <- !is.na(names(met)) & nzchar(trimws(names(met)))
met <- met[, valid_names, drop = FALSE]

# Also remove common exported row-index columns if present.
index_cols <- names(met)[
  names(met) %in% c("X", "...1") |
    grepl("^Unnamed", names(met), ignore.case = TRUE)
]

if (length(index_cols) > 0) {
  met <- met |>
    dplyr::select(-dplyr::all_of(index_cols))
}

required_metadata <- c(
  "study_id",
  "study_group",
  "ethnicity",
  "race",
  "age",
  "bmi",
  "antibiotics_during_pregnancy",
  "gestational_age",
  "gest_diabetes",
  "meds_antihypertensives",
  "bomeddischarge",
  "ppsystolic",
  "ppdiastolic",
  "ppmed"
)

missing_metadata <- setdiff(
  required_metadata,
  names(met)
)

if (length(missing_metadata) > 0) {
  stop(
    "Missing metadata variables: ",
    paste(missing_metadata, collapse = ", ")
  )
}

clean_study_id <- function(x) {
  x <- as.character(x)
  x <- stringr::str_trim(x)
  x <- stringr::str_remove(x, "^M")
  x <- stringr::str_remove(x, "Stool.*$")
  x <- stringr::str_remove(x, "^X")
  x
}

normalize_yes_no <- function(x) {
  x <- as.character(x)
  dplyr::case_when(
    x %in% c("Yes", "YES", "yes", "1", "TRUE", "True") ~ "Yes",
    x %in% c("No", "NO", "no", "0", "FALSE", "False") ~ "No",
    TRUE ~ NA_character_
  )
}

met <- met |>
  dplyr::mutate(
    study_id = clean_study_id(study_id),
    study_group = factor(
      as.character(study_group),
      levels = c("Control", "PE")
    ),
    ethnicity = factor(
      as.character(ethnicity),
      levels = c(
        "Not Hispanic/Latina",
        "Hispanic/Latina"
      )
    ),
    antibiotics_during_pregnancy = factor(
      normalize_yes_no(antibiotics_during_pregnancy),
      levels = c("No", "Yes")
    ),
    gest_diabetes = factor(
      normalize_yes_no(gest_diabetes),
      levels = c("No", "Yes")
    ),
    meds_antihypertensives = factor(
      normalize_yes_no(meds_antihypertensives),
      levels = c("No", "Yes")
    ),
    bomeddischarge = factor(
      normalize_yes_no(bomeddischarge),
      levels = c("No", "Yes")
    ),
    ppmed = factor(
      normalize_yes_no(ppmed),
      levels = c("No", "Yes")
    ),
    bp_followup = !is.na(ppsystolic) &
      !is.na(ppdiastolic),
    ethnicity_analysis = dplyr::case_when(
      ethnicity == "Hispanic/Latina" ~ "Hispanic",
      ethnicity == "Not Hispanic/Latina" &
        race == "White/Caucasian" ~ "White",
      TRUE ~ "Other"
    ),
    ethnicity_analysis = factor(
      ethnicity_analysis,
      levels = c("White", "Hispanic", "Other")
    ),
    term_delivery = factor(
      dplyr::if_else(
        gestational_age >= 37,
        "Term",
        "Preterm"
      ),
      levels = c("Preterm", "Term")
    )
  )

# Follow-up completeness

flow <- tibble::tibble(
  stage = c(
    "Enrolled",
    "Postpartum systolic and diastolic BP available",
    "Postpartum BP unavailable",
    "BP available and 6-week BP medication documented"
  ),
  n = c(
    nrow(met),
    sum(met$bp_followup),
    sum(!met$bp_followup),
    sum(met$bp_followup & !is.na(met$ppmed))
  )
)

readr::write_csv(
  flow,
  file.path(
    OUT_DIR,
    "BP_followup_flow.csv"
  )
)

flow_by_group <- met |>
  dplyr::count(
    study_group,
    bp_followup,
    name = "n"
  ) |>
  dplyr::group_by(study_group) |>
  dplyr::mutate(
    group_total = sum(n),
    percent = 100 * n / group_total
  ) |>
  dplyr::ungroup()

readr::write_csv(
  flow_by_group,
  file.path(
    OUT_DIR,
    "BP_followup_flow_by_study_group.csv"
  )
)

followup_table <- table(
  met$study_group,
  met$bp_followup
)

followup_fisher <- fisher.test(
  followup_table
)

readr::write_csv(
  tibble::tibble(
    comparison = "BP follow-up completeness by preeclampsia status",
    p_value = followup_fisher$p.value
  ),
  file.path(
    OUT_DIR,
    "BP_followup_by_study_group_Fisher.csv"
  )
)

# Compare participants with and without follow-up

numeric_followup_comparison <- function(data, variable) {
  
  x1 <- data[
    data$bp_followup,
    variable,
    drop = TRUE
  ]
  
  x0 <- data[
    !data$bp_followup,
    variable,
    drop = TRUE
  ]
  
  x1 <- x1[!is.na(x1)]
  x0 <- x0[!is.na(x0)]
  
  test <- t.test(
    x1,
    x0
  )
  
  pooled_sd <- sqrt(
    (
      stats::var(x1) +
        stats::var(x0)
    ) / 2
  )
  
  smd <- ifelse(
    pooled_sd > 0,
    (
      mean(x1) -
        mean(x0)
    ) / pooled_sd,
    NA_real_
  )
  
  tibble::tibble(
    variable = variable,
    type = "Continuous",
    followup = sprintf(
      "%.2f (%.2f)",
      mean(x1),
      sd(x1)
    ),
    no_followup = sprintf(
      "%.2f (%.2f)",
      mean(x0),
      sd(x0)
    ),
    p_value = test$p.value,
    standardized_difference = smd
  )
}

categorical_followup_comparison <- function(
    data,
    variable,
    positive_level) {
  
  x <- data[[variable]]
  
  follow_yes <- sum(
    data$bp_followup &
      !is.na(x) &
      as.character(x) == positive_level
  )
  
  follow_n <- sum(
    data$bp_followup &
      !is.na(x)
  )
  
  nofollow_yes <- sum(
    !data$bp_followup &
      !is.na(x) &
      as.character(x) == positive_level
  )
  
  nofollow_n <- sum(
    !data$bp_followup &
      !is.na(x)
  )
  
  if (follow_n == 0 || nofollow_n == 0) {
    return(
      tibble::tibble(
        variable = variable,
        type = "Categorical",
        followup = NA_character_,
        no_followup = NA_character_,
        p_value = NA_real_,
        standardized_difference = NA_real_
      )
    )
  }
  
  mat <- matrix(
    c(
      follow_yes,
      follow_n - follow_yes,
      nofollow_yes,
      nofollow_n - nofollow_yes
    ),
    nrow = 2,
    byrow = TRUE
  )
  
  p_value <- fisher.test(mat)$p.value
  
  p1 <- follow_yes / follow_n
  p0 <- nofollow_yes / nofollow_n
  
  denom <- sqrt(
    (
      p1 * (1 - p1) +
        p0 * (1 - p0)
    ) / 2
  )
  
  smd <- ifelse(
    denom > 0,
    (p1 - p0) / denom,
    NA_real_
  )
  
  tibble::tibble(
    variable = variable,
    type = "Categorical",
    followup = sprintf(
      "%d/%d (%.1f%%)",
      follow_yes,
      follow_n,
      100 * p1
    ),
    no_followup = sprintf(
      "%d/%d (%.1f%%)",
      nofollow_yes,
      nofollow_n,
      100 * p0
    ),
    p_value = p_value,
    standardized_difference = smd
  )
}

followup_comparison <- dplyr::bind_rows(
  numeric_followup_comparison(
    met,
    "age"
  ),
  numeric_followup_comparison(
    met,
    "bmi"
  ),
  numeric_followup_comparison(
    met,
    "gestational_age"
  ),
  categorical_followup_comparison(
    met,
    "study_group",
    "PE"
  ),
  categorical_followup_comparison(
    met,
    "ethnicity",
    "Hispanic/Latina"
  ),
  categorical_followup_comparison(
    met,
    "antibiotics_during_pregnancy",
    "Yes"
  ),
  categorical_followup_comparison(
    met,
    "gest_diabetes",
    "Yes"
  ),
  categorical_followup_comparison(
    met,
    "meds_antihypertensives",
    "Yes"
  ),
  categorical_followup_comparison(
    met,
    "term_delivery",
    "Term"
  )
)

readr::write_csv(
  followup_comparison,
  file.path(
    OUT_DIR,
    "BP_followup_participant_comparison.csv"
  )
)

# BP and medication summaries

bp_summary_by_group <- met |>
  dplyr::filter(bp_followup) |>
  dplyr::group_by(study_group) |>
  dplyr::summarise(
    n = dplyr::n(),
    systolic_mean = mean(ppsystolic),
    systolic_sd = sd(ppsystolic),
    diastolic_mean = mean(ppdiastolic),
    diastolic_sd = sd(ppdiastolic),
    .groups = "drop"
  )

readr::write_csv(
  bp_summary_by_group,
  file.path(
    OUT_DIR,
    "BP_summary_by_study_group.csv"
  )
)

postpartum_medication_summary <- met |>
  dplyr::filter(bp_followup) |>
  dplyr::count(
    study_group,
    ppmed,
    name = "n"
  ) |>
  dplyr::group_by(study_group) |>
  dplyr::mutate(
    denominator_nonmissing = sum(
      n[!is.na(ppmed)]
    ),
    percent_nonmissing = dplyr::if_else(
      !is.na(ppmed),
      100 * n / denominator_nonmissing,
      NA_real_
    )
  ) |>
  dplyr::ungroup()

readr::write_csv(
  postpartum_medication_summary,
  file.path(
    OUT_DIR,
    "postpartum_medication_summary.csv"
  )
)

ppmed_completeness <- met |>
  dplyr::filter(bp_followup) |>
  dplyr::summarise(
    n_with_BP = dplyr::n(),
    n_with_ppmed_documented = sum(!is.na(ppmed)),
    n_missing_ppmed = sum(is.na(ppmed)),
    percent_ppmed_documented = 100 * mean(!is.na(ppmed))
  )

readr::write_csv(
  ppmed_completeness,
  file.path(
    OUT_DIR,
    "postpartum_medication_completeness.csv"
  )
)

medication_at_discharge_summary <- met |>
  dplyr::count(
    study_group,
    bomeddischarge,
    name = "n"
  ) |>
  dplyr::group_by(study_group) |>
  dplyr::mutate(
    denominator_nonmissing = sum(
      n[!is.na(bomeddischarge)]
    ),
    percent_nonmissing = dplyr::if_else(
      !is.na(bomeddischarge),
      100 * n / denominator_nonmissing,
      NA_real_
    )
  ) |>
  dplyr::ungroup()

readr::write_csv(
  medication_at_discharge_summary,
  file.path(
    OUT_DIR,
    "BP_medication_at_discharge_summary.csv"
  )
)

# Load oral and vaginal phyloseq object

maternal_16s_env <- new.env()

load(
  MATERNAL_16S_FILE,
  envir = maternal_16s_env
)

if (
  !"PHYrar_alpha_ggext_subset02_newseqs" %in%
  ls(maternal_16s_env)
) {
  stop(
    "PHYrar_alpha_ggext_subset02_newseqs was not found in the RData file."
  )
}

phy_16s <- maternal_16s_env$PHYrar_alpha_ggext_subset02_newseqs

# Load maternal gut phyloseq object

if (!file.exists(MATERNAL_GUT_FILE)) {
  stop(
    "Maternal gut phyloseq file not found: ",
    MATERNAL_GUT_FILE
  )
}

phy_gut <- readRDS(
  MATERNAL_GUT_FILE
)

# Helper to find pre-computed alpha-diversity columns

find_metric_column <- function(
    data,
    candidates) {
  
  available_lower <- tolower(
    names(data)
  )
  
  candidate_lower <- tolower(
    candidates
  )
  
  hit <- match(
    candidate_lower,
    available_lower
  )
  
  hit <- hit[
    !is.na(hit)
  ]
  
  if (length(hit) == 0) {
    return(NA_character_)
  }
  
  names(data)[
    hit[1]
  ]
}

# Prepare oral/vaginal alpha metrics

prepare_16s_site <- function(
    phy,
    body_site_value,
    site_label) {
  
  # Do not use subset_samples() with a function argument because
  # subset_samples evaluates expressions inside sample_data and cannot
  # reliably see the external object body_site_value.
  sd_all <- data.frame(
    phyloseq::sample_data(phy),
    check.names = FALSE
  )
  
  if (!"body_site" %in% names(sd_all)) {
    stop("body_site is not present in sample_data.")
  }
  
  keep_samples <- rownames(sd_all)[
    as.character(sd_all$body_site) == body_site_value
  ]
  
  if (length(keep_samples) == 0) {
    stop(
      "No samples found for body_site = ",
      body_site_value
    )
  }
  
  phy_site <- phyloseq::prune_samples(
    keep_samples,
    phy
  )
  
  phy_site <- phyloseq::prune_samples(
    phyloseq::sample_sums(phy_site) > 0,
    phy_site
  )
  
  phy_site <- phyloseq::prune_taxa(
    phyloseq::taxa_sums(phy_site) > 0,
    phy_site
  )
  
  sd <- data.frame(
    phyloseq::sample_data(phy_site),
    check.names = FALSE
  ) |>
    tibble::rownames_to_column(
      "sample_id"
    )
  
  observed_col <- find_metric_column(
    sd,
    c(
      "observed",
      "Observed",
      "obs",
      "observed_rar_asv"
    )
  )
  
  shannon_col <- find_metric_column(
    sd,
    c(
      "shannon",
      "Shannon",
      "shannon_rar_asv"
    )
  )
  
  simpson_col <- find_metric_column(
    sd,
    c(
      "simpson",
      "Simpson",
      "simpson_rar_asv"
    )
  )
  
  if (
    any(
      is.na(
        c(
          observed_col,
          shannon_col,
          simpson_col
        )
      )
    )
  ) {
    
    alpha_calc <- phyloseq::estimate_richness(
      phy_site,
      measures = c(
        "Observed",
        "Shannon",
        "Simpson"
      )
    ) |>
      tibble::rownames_to_column(
        "sample_id"
      )
    
    sd <- sd |>
      dplyr::left_join(
        alpha_calc,
        by = "sample_id"
      )
    
    observed_col <- ifelse(
      is.na(observed_col),
      "Observed",
      observed_col
    )
    
    shannon_col <- ifelse(
      is.na(shannon_col),
      "Shannon",
      shannon_col
    )
    
    simpson_col <- ifelse(
      is.na(simpson_col),
      "Simpson",
      simpson_col
    )
  }
  
  if ("study_id" %in% names(sd)) {
    ids <- clean_study_id(
      sd$study_id
    )
  } else {
    ids <- clean_study_id(
      sd$sample_id
    )
  }
  
  tibble::tibble(
    study_id = ids,
    sample_id = sd$sample_id,
    body_site = site_label,
    Observed = as.numeric(
      sd[[observed_col]]
    ),
    Shannon = as.numeric(
      sd[[shannon_col]]
    ),
    Simpson = as.numeric(
      sd[[simpson_col]]
    ),
    alpha_source = paste(
      observed_col,
      shannon_col,
      simpson_col,
      sep = ";"
    )
  )
}

# Prepare maternal gut alpha metrics

prepare_gut_site <- function(
    phy,
    site_label = "Gut") {
  
  phy <- phyloseq::prune_samples(
    phyloseq::sample_sums(phy) > 0,
    phy
  )
  
  phy <- phyloseq::prune_taxa(
    phyloseq::taxa_sums(phy) > 0,
    phy
  )
  
  sd <- data.frame(
    phyloseq::sample_data(phy),
    check.names = FALSE
  ) |>
    tibble::rownames_to_column(
      "sample_id"
    )
  
  observed_col <- find_metric_column(
    sd,
    c(
      "richness",
      "Richness",
      "observed",
      "Observed",
      "obs",
      "species_richness"
    )
  )
  
  shannon_col <- find_metric_column(
    sd,
    c(
      "shannon",
      "Shannon"
    )
  )
  
  simpson_col <- find_metric_column(
    sd,
    c(
      "simpson",
      "Simpson"
    )
  )
  
  alpha_source <- "precomputed_sample_data"
  
  if (
    any(
      is.na(
        c(
          observed_col,
          shannon_col,
          simpson_col
        )
      )
    )
  ) {
    
    otu <- as(
      phyloseq::otu_table(phy),
      "matrix"
    )
    
    if (
      phyloseq::taxa_are_rows(phy)
    ) {
      otu <- t(otu)
    }
    
    alpha_calc <- tibble::tibble(
      sample_id = rownames(otu),
      Observed_calc = rowSums(
        otu > 0
      ),
      Shannon_calc = vegan::diversity(
        otu,
        index = "shannon"
      ),
      Simpson_calc = vegan::diversity(
        otu,
        index = "simpson"
      )
    )
    
    sd <- sd |>
      dplyr::left_join(
        alpha_calc,
        by = "sample_id"
      )
    
    observed_col <- ifelse(
      is.na(observed_col),
      "Observed_calc",
      observed_col
    )
    
    shannon_col <- ifelse(
      is.na(shannon_col),
      "Shannon_calc",
      shannon_col
    )
    
    simpson_col <- ifelse(
      is.na(simpson_col),
      "Simpson_calc",
      simpson_col
    )
    
    alpha_source <- paste0(
      "fallback_from_phyloseq_abundance_table; ",
      "preferred precomputed columns not all found"
    )
  }
  
  if ("study_id" %in% names(sd)) {
    ids <- clean_study_id(
      sd$study_id
    )
  } else {
    ids <- clean_study_id(
      sd$sample_id
    )
  }
  
  tibble::tibble(
    study_id = ids,
    sample_id = sd$sample_id,
    body_site = site_label,
    Observed = as.numeric(
      sd[[observed_col]]
    ),
    Shannon = as.numeric(
      sd[[shannon_col]]
    ),
    Simpson = as.numeric(
      sd[[simpson_col]]
    ),
    alpha_source = paste(
      alpha_source,
      observed_col,
      shannon_col,
      simpson_col,
      sep = ";"
    )
  )
}

oral_alpha <- prepare_16s_site(
  phy_16s,
  "MomTongue",
  "Oral"
)

vaginal_alpha <- prepare_16s_site(
  phy_16s,
  "Vaginal",
  "Vaginal"
)

gut_alpha <- prepare_gut_site(
  phy_gut,
  "Gut"
)

alpha_all <- dplyr::bind_rows(
  oral_alpha,
  vaginal_alpha,
  gut_alpha
)

# Check duplicate participants within site

duplicates <- alpha_all |>
  dplyr::count(
    body_site,
    study_id,
    name = "n"
  ) |>
  dplyr::filter(
    n > 1
  )

if (nrow(duplicates) > 0) {
  print(duplicates)
  stop(
    "Duplicate study_id values were found within at least one body site."
  )
}

readr::write_csv(
  alpha_all |>
    dplyr::distinct(
      body_site,
      alpha_source
    ),
  file.path(
    OUT_DIR,
    "alpha_metric_source_by_body_site.csv"
  )
)

# Merge with master metadata

analysis_all <- alpha_all |>
  dplyr::left_join(
    met,
    by = "study_id"
  ) |>
  dplyr::group_by(
    body_site
  ) |>
  dplyr::mutate(
    Observed_z = as.numeric(
      scale(Observed)
    ),
    Shannon_z = as.numeric(
      scale(Shannon)
    ),
    Simpson_z = as.numeric(
      scale(Simpson)
    )
  ) |>
  dplyr::ungroup()

readr::write_csv(
  analysis_all |>
    dplyr::select(
      study_id,
      sample_id,
      body_site,
      study_group,
      ethnicity,
      ethnicity_analysis,
      race,
      bmi,
      antibiotics_during_pregnancy,
      gestational_age,
      ppmed,
      ppsystolic,
      ppdiastolic,
      Observed,
      Shannon,
      Simpson,
      Observed_z,
      Shannon_z,
      Simpson_z,
      alpha_source
    ),
  file.path(
    OUT_DIR,
    "maternal_all_sites_alpha_BP_analysis_dataset.csv"
  )
)

# Sample counts available for BP analysis by body site

site_counts <- analysis_all |>
  dplyr::group_by(
    body_site
  ) |>
  dplyr::summarise(
    n_microbiome = dplyr::n(),
    n_with_BP = sum(
      !is.na(ppsystolic) &
        !is.na(ppdiastolic)
    ),
    n_with_BP_and_ppmed = sum(
      !is.na(ppsystolic) &
        !is.na(ppdiastolic) &
        !is.na(ppmed)
    ),
    .groups = "drop"
  )

readr::write_csv(
  site_counts,
  file.path(
    OUT_DIR,
    "BP_analysis_sample_counts_by_body_site.csv"
  )
)

# Robust-SE model helper

run_hc3_model <- function(
    data,
    site,
    outcome,
    alpha_metric,
    model_type,
    add_gestational_age = FALSE) {
  
  alpha_z <- paste0(
    alpha_metric,
    "_z"
  )
  
  covariates <- c(
    "study_group",
    "ethnicity",
    "bmi",
    "antibiotics_during_pregnancy",
    "ppmed"
  )
  
  if (add_gestational_age) {
    covariates <- c(
      covariates,
      "gestational_age"
    )
  }
  
  model_variables <- c(
    outcome,
    alpha_z,
    covariates
  )
  
  d <- data |>
    dplyr::filter(
      body_site == site
    ) |>
    dplyr::select(
      dplyr::all_of(
        model_variables
      )
    ) |>
    tidyr::drop_na() |>
    droplevels()
  
  if (nrow(d) < 20) {
    return(
      tibble::tibble(
        body_site = site,
        outcome = outcome,
        alpha_metric = alpha_metric,
        model_type = model_type,
        n = nrow(d),
        estimate_mmHg_per_SD = NA_real_,
        std_error_HC3 = NA_real_,
        conf_low = NA_real_,
        conf_high = NA_real_,
        p_value = NA_real_,
        formula = NA_character_,
        status = "Too few complete cases"
      )
    )
  }
  
  formula_text <- paste(
    outcome,
    "~",
    paste(
      c(
        alpha_z,
        covariates
      ),
      collapse = " + "
    )
  )
  
  fit <- lm(
    stats::as.formula(
      formula_text
    ),
    data = d
  )
  
  robust_vcov <- sandwich::vcovHC(
    fit,
    type = "HC3"
  )
  
  ct <- lmtest::coeftest(
    fit,
    vcov. = robust_vcov
  )
  
  if (
    !alpha_z %in%
    rownames(ct)
  ) {
    return(
      tibble::tibble(
        body_site = site,
        outcome = outcome,
        alpha_metric = alpha_metric,
        model_type = model_type,
        n = nrow(d),
        estimate_mmHg_per_SD = NA_real_,
        std_error_HC3 = NA_real_,
        conf_low = NA_real_,
        conf_high = NA_real_,
        p_value = NA_real_,
        formula = formula_text,
        status = "Alpha term not estimable"
      )
    )
  }
  
  term_row <- ct[
    alpha_z,
    ,
    drop = FALSE
  ]
  
  estimate <- unname(
    term_row[1, "Estimate"]
  )
  
  se <- unname(
    term_row[1, "Std. Error"]
  )
  
  p_value <- unname(
    term_row[1, "Pr(>|t|)"]
  )
  
  critical_value <- qt(
    0.975,
    df = df.residual(fit)
  )
  
  tibble::tibble(
    body_site = site,
    outcome = outcome,
    alpha_metric = alpha_metric,
    model_type = model_type,
    n = nrow(d),
    estimate_mmHg_per_SD = estimate,
    std_error_HC3 = se,
    conf_low = estimate -
      critical_value * se,
    conf_high = estimate +
      critical_value * se,
    p_value = p_value,
    formula = formula_text,
    status = "Completed"
  )
}

# Run the primary tests.
#
# Benjamini-Hochberg correction is applied separately within each body site,
# blood-pressure outcome, and model, across the three prespecified alpha
# diversity metrics (Observed, Shannon, and Simpson). This is the final
# multiple-testing family used for Table S9 and Figures 1 and S2.

sites <- c(
  "Oral",
  "Vaginal",
  "Gut"
)

alpha_metrics <- c(
  "Observed",
  "Shannon",
  "Simpson"
)

bp_outcomes <- c(
  "ppsystolic",
  "ppdiastolic"
)

primary_results <- list()

counter <- 1

for (site in sites) {
  for (metric in alpha_metrics) {
    for (outcome in bp_outcomes) {
      
      primary_results[[counter]] <-
        run_hc3_model(
          analysis_all,
          site = site,
          outcome = outcome,
          alpha_metric = metric,
          model_type = "Primary adjusted",
          add_gestational_age = FALSE
        )
      
      counter <- counter + 1
    }
  }
}

primary_results <- dplyr::bind_rows(
  primary_results
) |>
  dplyr::group_by(
    body_site,
    outcome
  ) |>
  dplyr::mutate(
    adjusted_p = p.adjust(
      p_value,
      method = "BH"
    ),
    nominal_significance = dplyr::case_when(
      p_value < 0.05 ~ "p<0.05",
      p_value >= 0.05 &
        p_value <= 0.08 ~ "p=0.05-0.08",
      TRUE ~ "NS"
    ),
    fdr_significance = dplyr::if_else(
      adjusted_p < 0.05,
      "adjusted p<0.05",
      "NS"
    )
  ) |>
  dplyr::ungroup() |>
  dplyr::arrange(
    body_site,
    outcome,
    alpha_metric
  )

readr::write_csv(
  primary_results,
  file.path(
    OUT_DIR,
    "ALL_SITES_alpha_BP_primary_BH_within_site_BP_outcome.csv"
  )
)

# Gestational-age sensitivity. The same BH family is applied independently
# within each site and BP outcome.

ga_results <- list()

counter <- 1

for (site in sites) {
  for (metric in alpha_metrics) {
    for (outcome in bp_outcomes) {
      
      ga_results[[counter]] <-
        run_hc3_model(
          analysis_all,
          site = site,
          outcome = outcome,
          alpha_metric = metric,
          model_type = "GA sensitivity",
          add_gestational_age = TRUE
        )
      
      counter <- counter + 1
    }
  }
}

ga_results <- dplyr::bind_rows(
  ga_results
) |>
  dplyr::group_by(
    body_site,
    outcome
  ) |>
  dplyr::mutate(
    adjusted_p = p.adjust(
      p_value,
      method = "BH"
    ),
    nominal_significance = dplyr::case_when(
      p_value < 0.05 ~ "p<0.05",
      p_value >= 0.05 &
        p_value <= 0.08 ~ "p=0.05-0.08",
      TRUE ~ "NS"
    ),
    fdr_significance = dplyr::if_else(
      adjusted_p < 0.05,
      "adjusted p<0.05",
      "NS"
    )
  ) |>
  dplyr::ungroup() |>
  dplyr::arrange(
    body_site,
    outcome,
    alpha_metric
  )

readr::write_csv(
  ga_results,
  file.path(
    OUT_DIR,
    "ALL_SITES_alpha_BP_GA_sensitivity_BH_within_site_BP_outcome.csv"
  )
)

# Formal ethnicity-interaction analyses
# Hispanic women are compared with non-Hispanic White women.
# These tests determine whether the diversity-BP slope differs by ethnicity.

run_ethnicity_interaction <- function(
    data,
    site,
    outcome,
    alpha_metric) {
  
  alpha_z <- paste0(
    alpha_metric,
    "_z"
  )
  
  required <- c(
    outcome,
    alpha_z,
    "ethnicity_analysis",
    "study_group",
    "bmi",
    "antibiotics_during_pregnancy",
    "ppmed"
  )
  
  d <- data |>
    dplyr::filter(
      body_site == site,
      ethnicity_analysis %in%
        c(
          "White",
          "Hispanic"
        )
    ) |>
    dplyr::select(
      dplyr::all_of(
        required
      )
    ) |>
    tidyr::drop_na() |>
    droplevels()
  
  if (nrow(d) < 20) {
    return(
      list(
        interaction = tibble::tibble(
          body_site = site,
          outcome = outcome,
          alpha_metric = alpha_metric,
          n = nrow(d),
          estimate = NA_real_,
          std_error_HC3 = NA_real_,
          conf_low = NA_real_,
          conf_high = NA_real_,
          p_value = NA_real_,
          formula = NA_character_,
          status = "Too few complete cases"
        ),
        simple_slopes = tibble::tibble()
      )
    )
  }
  
  formula_text <- paste0(
    outcome,
    " ~ ",
    alpha_z,
    " * ethnicity_analysis + ",
    "study_group + bmi + ",
    "antibiotics_during_pregnancy + ppmed"
  )
  
  fit <- lm(
    stats::as.formula(
      formula_text
    ),
    data = d
  )
  
  robust_vcov <- sandwich::vcovHC(
    fit,
    type = "HC3"
  )
  
  ct <- lmtest::coeftest(
    fit,
    vcov. = robust_vcov
  )
  
  interaction_term <- paste0(
    alpha_z,
    ":ethnicity_analysisHispanic"
  )
  
  if (
    !interaction_term %in%
    rownames(ct)
  ) {
    interaction_term <- paste0(
      "ethnicity_analysisHispanic:",
      alpha_z
    )
  }
  
  if (
    !interaction_term %in%
    rownames(ct)
  ) {
    return(
      list(
        interaction = tibble::tibble(
          body_site = site,
          outcome = outcome,
          alpha_metric = alpha_metric,
          n = nrow(d),
          estimate = NA_real_,
          std_error_HC3 = NA_real_,
          conf_low = NA_real_,
          conf_high = NA_real_,
          p_value = NA_real_,
          formula = formula_text,
          status = "Interaction term not estimable"
        ),
        simple_slopes = tibble::tibble()
      )
    )
  }
  
  interaction_row <- ct[
    interaction_term,
    ,
    drop = FALSE
  ]
  
  estimate <- unname(
    interaction_row[1, "Estimate"]
  )
  
  se <- unname(
    interaction_row[1, "Std. Error"]
  )
  
  critical_value <- qt(
    0.975,
    df = df.residual(fit)
  )
  
  interaction_result <- tibble::tibble(
    body_site = site,
    outcome = outcome,
    alpha_metric = alpha_metric,
    n = nrow(d),
    estimate = estimate,
    std_error_HC3 = se,
    conf_low = estimate -
      critical_value * se,
    conf_high = estimate +
      critical_value * se,
    p_value = unname(
      interaction_row[1, "Pr(>|t|)"]
    ),
    formula = formula_text,
    status = "Completed"
  )
  
  trend <- emmeans::emtrends(
    fit,
    ~ ethnicity_analysis,
    var = alpha_z,
    vcov. = robust_vcov
  )
  
  trend_result <- as.data.frame(
    summary(
      trend,
      infer = c(TRUE, TRUE)
    )
  ) |>
    tibble::as_tibble() |>
    dplyr::mutate(
      body_site = site,
      outcome = outcome,
      alpha_metric = alpha_metric,
      n = nrow(d),
      .before = 1
    )
  
  list(
    interaction = interaction_result,
    simple_slopes = trend_result
  )
}

interaction_objects <- list()

counter <- 1

for (site in sites) {
  for (metric in alpha_metrics) {
    for (outcome in bp_outcomes) {
      
      interaction_objects[[counter]] <-
        run_ethnicity_interaction(
          analysis_all,
          site = site,
          outcome = outcome,
          alpha_metric = metric
        )
      
      counter <- counter + 1
    }
  }
}

ethnicity_interactions <- dplyr::bind_rows(
  lapply(
    interaction_objects,
    `[[`,
    "interaction"
  )
) |>
  dplyr::group_by(
    body_site,
    outcome
  ) |>
  dplyr::mutate(
    adjusted_p = p.adjust(
      p_value,
      method = "BH"
    ),
    nominal_significance = dplyr::case_when(
      p_value < 0.05 ~ "p<0.05",
      p_value >= 0.05 &
        p_value <= 0.08 ~ "p=0.05-0.08",
      TRUE ~ "NS"
    ),
    fdr_significance = dplyr::if_else(
      adjusted_p < 0.05,
      "adjusted p<0.05",
      "NS"
    )
  ) |>
  dplyr::ungroup() |>
  dplyr::arrange(
    body_site,
    outcome,
    alpha_metric
  )

readr::write_csv(
  ethnicity_interactions,
  file.path(
    OUT_DIR,
    "ALL_SITES_alpha_BP_ethnicity_interaction_BH_within_site_BP_outcome.csv"
  )
)

ethnicity_simple_slopes <- dplyr::bind_rows(
  lapply(
    interaction_objects,
    `[[`,
    "simple_slopes"
  )
) |>
  dplyr::group_by(
    body_site,
    outcome,
    ethnicity_analysis
  ) |>
  dplyr::mutate(
    adjusted_p = p.adjust(
      p.value,
      method = "BH"
    )
  ) |>
  dplyr::ungroup()

readr::write_csv(
  ethnicity_simple_slopes,
  file.path(
    OUT_DIR,
    "ALL_SITES_alpha_BP_ethnicity_simple_slopes_BH_within_site_BP_outcome_group.csv"
  )
)

# Inverse-probability-of-follow-up sensitivity
# Weights are estimated in all 100 women, before restricting to follow-up.

selection_model <- glm(
  bp_followup ~
    study_group +
    bmi +
    gestational_age +
    ethnicity,
  family = binomial(),
  data = met
)

met$followup_probability <- predict(
  selection_model,
  type = "response"
)

overall_followup_rate <- mean(
  met$bp_followup
)

met$stabilized_followup_weight <- ifelse(
  met$bp_followup,
  overall_followup_rate /
    met$followup_probability,
  (1 - overall_followup_rate) /
    (1 - met$followup_probability)
)

weight_summary <- met |>
  dplyr::filter(
    bp_followup
  ) |>
  dplyr::summarise(
    n = dplyr::n(),
    mean_weight = mean(
      stabilized_followup_weight
    ),
    sd_weight = sd(
      stabilized_followup_weight
    ),
    min_weight = min(
      stabilized_followup_weight
    ),
    median_weight = median(
      stabilized_followup_weight
    ),
    max_weight = max(
      stabilized_followup_weight
    )
  )

readr::write_csv(
  weight_summary,
  file.path(
    OUT_DIR,
    "followup_IPW_weight_summary.csv"
  )
)

analysis_ipw <- analysis_all |>
  dplyr::left_join(
    met |>
      dplyr::select(
        study_id,
        stabilized_followup_weight
      ),
    by = "study_id"
  )

run_ipw_model <- function(
    data,
    site,
    outcome,
    alpha_metric) {
  
  alpha_z <- paste0(
    alpha_metric,
    "_z"
  )
  
  required <- c(
    outcome,
    alpha_z,
    "study_group",
    "ethnicity",
    "bmi",
    "antibiotics_during_pregnancy",
    "ppmed",
    "stabilized_followup_weight"
  )
  
  d <- data |>
    dplyr::filter(
      body_site == site
    ) |>
    dplyr::select(
      dplyr::all_of(
        required
      )
    ) |>
    tidyr::drop_na() |>
    droplevels()
  
  if (nrow(d) < 20) {
    return(
      tibble::tibble(
        body_site = site,
        outcome = outcome,
        alpha_metric = alpha_metric,
        n = nrow(d),
        estimate_mmHg_per_SD = NA_real_,
        std_error_HC3 = NA_real_,
        conf_low = NA_real_,
        conf_high = NA_real_,
        p_value = NA_real_,
        status = "Too few complete cases"
      )
    )
  }
  
  formula_text <- paste0(
    outcome,
    " ~ ",
    alpha_z,
    " + study_group + ethnicity + bmi + ",
    "antibiotics_during_pregnancy + ppmed"
  )
  
  fit <- lm(
    stats::as.formula(
      formula_text
    ),
    data = d,
    weights = stabilized_followup_weight
  )
  
  robust_vcov <- sandwich::vcovHC(
    fit,
    type = "HC3"
  )
  
  ct <- lmtest::coeftest(
    fit,
    vcov. = robust_vcov
  )
  
  row <- ct[
    alpha_z,
    ,
    drop = FALSE
  ]
  
  estimate <- unname(
    row[1, "Estimate"]
  )
  
  se <- unname(
    row[1, "Std. Error"]
  )
  
  critical_value <- qt(
    0.975,
    df = df.residual(fit)
  )
  
  tibble::tibble(
    body_site = site,
    outcome = outcome,
    alpha_metric = alpha_metric,
    n = nrow(d),
    estimate_mmHg_per_SD = estimate,
    std_error_HC3 = se,
    conf_low = estimate -
      critical_value * se,
    conf_high = estimate +
      critical_value * se,
    p_value = unname(
      row[1, "Pr(>|t|)"]
    ),
    status = "Completed"
  )
}

ipw_results <- list()

counter <- 1

for (site in sites) {
  for (metric in alpha_metrics) {
    for (outcome in bp_outcomes) {
      
      ipw_results[[counter]] <-
        run_ipw_model(
          analysis_ipw,
          site = site,
          outcome = outcome,
          alpha_metric = metric
        )
      
      counter <- counter + 1
    }
  }
}

ipw_results <- dplyr::bind_rows(
  ipw_results
) |>
  dplyr::group_by(
    body_site,
    outcome
  ) |>
  dplyr::mutate(
    adjusted_p = p.adjust(
      p_value,
      method = "BH"
    )
  ) |>
  dplyr::ungroup() |>
  dplyr::arrange(
    body_site,
    outcome,
    alpha_metric
  )

readr::write_csv(
  ipw_results,
  file.path(
    OUT_DIR,
    "ALL_SITES_alpha_BP_IPW_sensitivity_BH_within_site_BP_outcome.csv"
  )
)

# Model manifest

model_manifest <- tibble::tibble(
  analysis = c(
    "Primary",
    "Gestational-age sensitivity",
    "Ethnicity interaction",
    "Follow-up IPW sensitivity"
  ),
  model = c(
    paste(
      "6-week BP ~ alpha diversity z-score + preeclampsia status + ethnicity +",
      "pre-pregnancy BMI + antibiotics during pregnancy + 6-week BP medication"
    ),
    paste(
      "Primary model + gestational age at delivery"
    ),
    paste(
      "6-week BP ~ alpha diversity z-score * ethnicity + preeclampsia status +",
      "pre-pregnancy BMI + antibiotics during pregnancy + 6-week BP medication;",
      "Hispanic and non-Hispanic White women only"
    ),
    paste(
      "Primary model weighted by inverse probability of BP follow-up"
    )
  ),
  multiple_testing_family = c(
    "BH within each body site and BP outcome across 3 alpha metrics",
    "BH within each body site and BP outcome across 3 alpha metrics",
    "BH within each body site and BP outcome across 3 alpha-metric interaction tests",
    "BH within each body site and BP outcome across 3 alpha metrics"
  )
)

readr::write_csv(
  model_manifest,
  file.path(
    OUT_DIR,
    "BP_model_manifest.csv"
  )
)

# Console summary

cat("\nPostpartum BP reviewer analysis completed.\n")
cat("Results directory:\n", OUT_DIR, "\n\n")

cat("Primary results (BH within site and BP outcome):\n")
print(
  primary_results |>
    dplyr::select(
      body_site,
      outcome,
      alpha_metric,
      n,
      estimate_mmHg_per_SD,
      conf_low,
      conf_high,
      p_value,
      adjusted_p
    )
)

cat("\nFormal ethnicity interactions:\n")
print(
  ethnicity_interactions |>
    dplyr::select(
      body_site,
      outcome,
      alpha_metric,
      n,
      estimate,
      conf_low,
      conf_high,
      p_value,
      adjusted_p
    )
)

cat("\nReview these output files first:\n")
cat("  BP_followup_participant_comparison.csv\n")
cat("  BP_followup_flow_by_study_group.csv\n")
cat("  postpartum_medication_summary.csv\n")
cat("  BP_analysis_sample_counts_by_body_site.csv\n")
cat("  alpha_metric_source_by_body_site.csv\n")
cat("  ALL_SITES_alpha_BP_primary_BH_within_site_BP_outcome.csv\n")
cat("  ALL_SITES_alpha_BP_GA_sensitivity_BH_within_site_BP_outcome.csv\n")
cat("  ALL_SITES_alpha_BP_ethnicity_interaction_BH_within_site_BP_outcome.csv\n")
cat("  ALL_SITES_alpha_BP_ethnicity_simple_slopes_BH_within_site_BP_outcome_group.csv\n")
cat("  ALL_SITES_alpha_BP_IPW_sensitivity_BH_within_site_BP_outcome.csv\n")

capture.output(
  sessionInfo(),
  file = file.path(
    OUT_DIR,
    "sessionInfo.txt"
  )
)
