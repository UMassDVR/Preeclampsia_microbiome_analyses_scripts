#species_contributors####
# 03_PWY5265_species_contributors_UPDATED.R
#
# PREECLAMPSIA maternal gut WGS / HUMAnN3
#
# Targeted exploratory follow-up of:
#   PWY-5265: peptidoglycan biosynthesis II (staphylococci)
#
# This script evaluates which HUMAnN taxonomic strata contribute to PWY-5265
# and tests taxon-specific PWY-5265 abundance using the same reviewer-adjusted
# strategy used in the main maternal functional analysis.
#
# PE vs Control:
#   Prespecified minimal =
#     study_group + ethnicity + bmi + antibiotics_during_pregnancy
#
#   Fully adjusted =
#     study_group + ethnicity + bmi + antibiotics_during_pregnancy +
#     gestational_age + gest_diabetes + meds_antihypertensives
#
#   Both models are repeated among term deliveries (gestational_age >= 37).
#
# PE severity:
#   Primary parsimonious =
#     severe + ethnicity + bmi + antibiotics_during_pregnancy
#
#   Gestational-age sensitivity =
#     severe + ethnicity + bmi + antibiotics_during_pregnancy +
#     gestational_age
#
#   Both are also repeated among term deliveries.
#
# Ethnicity-stratified PE vs Control and severity:
#   Primary =
#     exposure + bmi + antibiotics_during_pregnancy
#
#   Gestational-age sensitivity =
#     exposure + bmi + antibiotics_during_pregnancy + gestational_age
#
# Important interpretation:
#   - PWY-5265 is a MetaCyc pathway name; "(staphylococci)" does not identify
#     the organism contributing the signal in these samples.
#   - The outcome here is taxon-stratified PWY-5265 abundance, not total
#     genomic abundance of the taxon.
#   - HUMAnN computes pathway abundance independently at the community and
#     taxonomic-stratum levels. Therefore, stratified pathway abundances do
#     not necessarily sum to the unstratified community pathway abundance.
#   - Taxonomic-stratum results are used to identify contributors and must not
#     be interpreted as an exact additive decomposition of the total pathway.
#   - This is a targeted exploratory follow-up. The parent PWY-5265 association
#     was nominal and did not survive FDR correction in the broad pathway scan.
#   - No stepwise variable selection is used.

# 0. Packages
required_packages <- c(
  "tidyverse",
  "Maaslin2"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_packages) > 0) {
  stop(
    "Install the following packages before running this script: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(Maaslin2)
})

# 1. Paths
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

HUMANN_DIR <- file.path(
  PROJECT_DIR,
  "WGS",
  "HumaNN3"
)

CLEAN_DIR <- file.path(
  HUMANN_DIR,
  "cleaned_100_samples"
)

METADATA_FILE <- file.path(
  PROJECT_DIR,
  "met",
  "met_clean0.csv"
)

UNSTRATIFIED_PATHWAY_FILE <- file.path(
  CLEAN_DIR,
  "PREECLAMPSIA_pathabundance_unstratified_cpm_clean100.tsv.gz"
)

STRATIFIED_PATHWAY_FILE <- file.path(
  CLEAN_DIR,
  "PREECLAMPSIA_pathabundance_stratified_cpm_clean100.tsv.gz"
)

OUTPUT_DIR <- file.path(
  HUMANN_DIR,
  "PWY5265_species_contributors"
)

dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

# 2. Analysis parameters
TARGET_PATHWAY_ID <- "PWY-5265"

# Targeted contributor filter used before MaAsLin2.
# A taxon-specific pathway feature must be detected at >=1 CPM in at least
# 10% of the analytical population and in at least 3 samples.
DETECTION_CPM <- 1
MIN_PREVALENCE <- 0.10
MIN_DETECTED_SAMPLES <- 3

MIN_MODEL_N <- 10
MIN_LEVEL_N <- 3

PRIMARY_FDR <- 0.05
EXPLORATORY_FDR <- 0.25
MARGINAL_P <- 0.08

available_cores <- parallel::detectCores(
  logical = TRUE
)

N_CORES <- max(
  1,
  min(
    4,
    ifelse(
      is.na(available_cores),
      1,
      available_cores - 1
    )
  )
)

TOP_N_CONTRIBUTORS <- 10

# 3. Helper functions
check_file <- function(path) {
  if (!file.exists(path)) {
    stop("Required file not found: ", path)
  }
}

safe_name <- function(x) {
  x |>
    stringr::str_replace_all("[^A-Za-z0-9._-]+", "_") |>
    stringr::str_replace_all("_+", "_") |>
    stringr::str_remove("^_") |>
    stringr::str_remove("_$")
}

# Match the raw HUMAnN ID. This is intentionally based on the pathway ID,
# rather than the descriptive pathway name.
is_target_pathway <- function(x) {
  stringr::str_detect(
    x,
    stringr::regex(
      "^PWY[-.]5265(?=[:.|]|$)",
      ignore_case = TRUE
    )
  )
}

extract_sample_ids <- function(sample_columns) {
  
  ids <- stringr::str_extract(
    sample_columns,
    "M[0-9]+Stool"
  )
  
  if (any(is.na(ids))) {
    bad <- sample_columns[is.na(ids)]
    stop(
      "Could not extract sample IDs from HUMAnN columns: ",
      paste(bad, collapse = ", ")
    )
  }
  
  if (anyDuplicated(ids) > 0) {
    stop("Duplicated HUMAnN sample IDs detected.")
  }
  
  ids
}

taxon_display_label <- function(x) {
  
  species <- stringr::str_extract(
    x,
    "s__[^|]+$"
  )
  
  genus <- stringr::str_extract(
    x,
    "g__[^.|\t]+"
  )
  
  out <- dplyr::case_when(
    !is.na(species) ~
      species |>
      stringr::str_remove("^s__") |>
      stringr::str_replace_all("_", " "),
    !is.na(genus) ~
      genus |>
      stringr::str_remove("^g__") |>
      stringr::str_replace_all("_", " "),
    TRUE ~ x
  )
  
  out
}

read_humann_table <- function(path) {
  
  check_file(path)
  
  x <- readr::read_tsv(
    path,
    show_col_types = FALSE,
    name_repair = "minimal",
    progress = TRUE
  )
  
  names(x)[1] <- "feature"
  
  x <- x |>
    dplyr::mutate(
      feature = as.character(feature)
    )
  
  sample_columns <- names(x)[-1]
  sample_ids <- extract_sample_ids(sample_columns)
  
  names(x)[-1] <- sample_ids
  
  x
}

matrix_from_humann_rows <- function(df, row_ids) {
  
  m <- as.matrix(
    df[
      ,
      -1,
      drop = FALSE
    ]
  )
  
  storage.mode(m) <- "numeric"
  
  if (any(!is.finite(m))) {
    stop("Non-finite abundance values found.")
  }
  
  if (any(m < 0)) {
    stop("Negative abundance values found.")
  }
  
  rownames(m) <- row_ids
  
  m
}

summarize_contributors_by_group <- function(
    long_data,
    group_var
) {
  
  group_sym <- rlang::sym(group_var)
  
  out <- long_data |>
    dplyr::filter(
      !is.na(!!group_sym)
    ) |>
    dplyr::group_by(
      !!group_sym,
      feature_id,
      contributor_raw,
      taxon_label
    ) |>
    dplyr::summarise(
      n_samples = dplyr::n(),
      n_detected = sum(cpm >= DETECTION_CPM),
      prevalence = mean(cpm >= DETECTION_CPM),
      mean_cpm = mean(cpm),
      median_cpm = median(cpm),
      sum_contributor_cpm = sum(cpm),
      .groups = "drop"
    ) |>
    dplyr::rename(
      group = !!group_sym
    ) |>
    dplyr::group_by(group) |>
    dplyr::mutate(
      stratified_share_percent =
        if (sum(sum_contributor_cpm) > 0) {
          100 * sum_contributor_cpm / sum(sum_contributor_cpm)
        } else {
          rep(NA_real_, dplyr::n())
        }
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      grouping_variable = group_var,
      .before = 1
    )
  
  out
}

compare_taxon_specific_group_means <- function(
    long_data,
    group_var,
    comparison_level,
    reference_level
) {
  
  group_sym <- rlang::sym(group_var)
  
  d <- long_data |>
    dplyr::filter(
      !is.na(!!group_sym),
      as.character(!!group_sym) %in%
        c(reference_level, comparison_level)
    )
  
  taxon_means <- d |>
    dplyr::group_by(
      feature_id,
      contributor_raw,
      taxon_label,
      !!group_sym
    ) |>
    dplyr::summarise(
      mean_taxon_specific_cpm = mean(cpm),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      group_character =
        as.character(!!group_sym)
    ) |>
    dplyr::select(
      -!!group_sym
    ) |>
    tidyr::pivot_wider(
      names_from = group_character,
      values_from = mean_taxon_specific_cpm,
      values_fill = 0
    )
  
  if (
    !all(
      c(reference_level, comparison_level) %in%
      names(taxon_means)
    )
  ) {
    stop(
      "Could not compare ",
      comparison_level,
      " vs ",
      reference_level,
      " for ",
      group_var
    )
  }
  
  total_means <- d |>
    dplyr::distinct(
      sample_id,
      !!group_sym,
      total_pathway_cpm
    ) |>
    dplyr::group_by(
      !!group_sym
    ) |>
    dplyr::summarise(
      mean_unstratified_pathway_cpm = mean(total_pathway_cpm),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      group_character =
        as.character(!!group_sym)
    )
  
  ref_total <- total_means |>
    dplyr::filter(
      group_character == reference_level
    ) |>
    dplyr::pull(
      mean_unstratified_pathway_cpm
    )
  
  comp_total <- total_means |>
    dplyr::filter(
      group_character == comparison_level
    ) |>
    dplyr::pull(
      mean_unstratified_pathway_cpm
    )
  
  taxon_means |>
    dplyr::mutate(
      reference_group = reference_level,
      comparison_group = comparison_level,
      mean_taxon_cpm_reference =
        .data[[reference_level]],
      mean_taxon_cpm_comparison =
        .data[[comparison_level]],
      delta_taxon_specific_cpm =
        mean_taxon_cpm_comparison -
        mean_taxon_cpm_reference,
      mean_unstratified_cpm_reference = ref_total,
      mean_unstratified_cpm_comparison = comp_total,
      delta_unstratified_pathway_cpm =
        comp_total - ref_total
    ) |>
    dplyr::select(
      feature_id,
      contributor_raw,
      taxon_label,
      reference_group,
      comparison_group,
      mean_taxon_cpm_reference,
      mean_taxon_cpm_comparison,
      delta_taxon_specific_cpm,
      mean_unstratified_cpm_reference,
      mean_unstratified_cpm_comparison,
      delta_unstratified_pathway_cpm
    ) |>
    dplyr::arrange(
      dplyr::desc(
        delta_taxon_specific_cpm
      )
    )
}

# 4. Read and prepare metadata
check_file(METADATA_FILE)

metadata <- readr::read_csv(
  METADATA_FILE,
  show_col_types = FALSE,
  name_repair = "minimal"
)

empty_names <- which(
  names(metadata) == "" |
    is.na(names(metadata))
)

if (length(empty_names) > 0) {
  names(metadata)[empty_names] <- paste0(
    "unnamed_",
    seq_along(empty_names)
  )
}

required_metadata_variables <- c(
  "study_id",
  "study_group",
  "severe",
  "severe_3cat",
  "bmi",
  "antibiotics_during_pregnancy",
  "ethnicity",
  "race",
  "gestational_age",
  "gest_diabetes",
  "meds_antihypertensives"
)

missing_metadata_variables <- setdiff(
  required_metadata_variables,
  names(metadata)
)

if (length(missing_metadata_variables) > 0) {
  stop(
    "Metadata is missing required variables: ",
    paste(
      missing_metadata_variables,
      collapse = ", "
    )
  )
}

metadata <- metadata |>
  dplyr::mutate(
    study_id = as.character(study_id),
    sample_id = paste0(
      "M",
      study_id,
      "Stool"
    ),
    study_group = factor(
      study_group,
      levels = c(
        "Control",
        "PE"
      )
    ),
    severe = factor(
      severe,
      levels = c(
        "No",
        "Yes"
      )
    ),
    severe_3cat = factor(
      severe_3cat,
      levels = c(
        "control",
        "No",
        "Yes"
      ),
      labels = c(
        "Control",
        "PE-NSF",
        "PE-SF"
      )
    ),
    antibiotics_during_pregnancy = factor(
      antibiotics_during_pregnancy,
      levels = c(
        "No",
        "Yes"
      )
    ),
    ethnicity = factor(
      ethnicity,
      levels = c(
        "Not Hispanic/Latina",
        "Hispanic/Latina"
      )
    ),
    race = factor(race),
    ethnicity_analysis = dplyr::case_when(
      ethnicity == "Hispanic/Latina" ~
        "Hispanic",
      ethnicity == "Not Hispanic/Latina" &
        race == "White/Caucasian" ~
        "White",
      TRUE ~
        "Other"
    ),
    ethnicity_analysis = factor(
      ethnicity_analysis,
      levels = c(
        "White",
        "Hispanic",
        "Other"
      )
    ),
    bmi = as.numeric(bmi),
    gestational_age = as.numeric(gestational_age),
    gest_diabetes = factor(
      gest_diabetes,
      levels = c(
        "No",
        "Yes"
      )
    ),
    meds_antihypertensives = factor(
      meds_antihypertensives,
      levels = c(
        "No",
        "Yes"
      )
    ),
    term_delivery = dplyr::case_when(
      is.na(gestational_age) ~ NA_character_,
      gestational_age >= 37 ~ "Term",
      TRUE ~ "Preterm"
    ),
    term_delivery = factor(
      term_delivery,
      levels = c(
        "Preterm",
        "Term"
      )
    )
  )

if (anyDuplicated(metadata$sample_id) > 0) {
  stop("Duplicated metadata sample IDs detected.")
}

# 5. Read HUMAnN pathway tables and identify PWY-5265
message("Reading unstratified pathway table...")
unstratified <- read_humann_table(
  UNSTRATIFIED_PATHWAY_FILE
)

message("Reading stratified pathway table...")
stratified <- read_humann_table(
  STRATIFIED_PATHWAY_FILE
)

# Require identical samples between the two pathway tables.
unstratified_ids <- names(unstratified)[-1]
stratified_ids <- names(stratified)[-1]

if (!setequal(unstratified_ids, stratified_ids)) {
  stop(
    "Unstratified and stratified pathway tables contain different samples."
  )
}

# Require exact metadata/HUMAnN reconciliation.
humann_ids <- unstratified_ids

humann_only <- setdiff(
  humann_ids,
  metadata$sample_id
)

metadata_only <- setdiff(
  metadata$sample_id,
  humann_ids
)

if (
  length(humann_only) > 0 ||
  length(metadata_only) > 0
) {
  stop(
    "Sample mismatch. HUMAnN-only: ",
    paste(humann_only, collapse = ", "),
    "; metadata-only: ",
    paste(metadata_only, collapse = ", ")
  )
}

# Reorder metadata to HUMAnN sample order.
metadata <- metadata[
  match(
    humann_ids,
    metadata$sample_id
  ),
  ,
  drop = FALSE
]

# Community-level, unstratified PWY-5265 row.
target_total_rows <- unstratified |>
  dplyr::filter(
    is_target_pathway(feature),
    !stringr::str_detect(
      feature,
      stringr::fixed("|")
    )
  )

if (nrow(target_total_rows) != 1) {
  candidate_rows <- unstratified |>
    dplyr::filter(
      stringr::str_detect(
        feature,
        stringr::regex(
          "5265",
          ignore_case = TRUE
        )
      )
    ) |>
    dplyr::pull(feature)
  
  stop(
    "Expected exactly one unstratified PWY-5265 row, found ",
    nrow(target_total_rows),
    ". Candidate rows: ",
    paste(candidate_rows, collapse = " ; ")
  )
}

TARGET_PATHWAY_FULL_NAME <- target_total_rows$feature[[1]]

message(
  "Target pathway identified as: ",
  TARGET_PATHWAY_FULL_NAME
)

total_matrix <- matrix_from_humann_rows(
  target_total_rows,
  row_ids = "PWY5265_total"
)

total_pathway_cpm <- as.numeric(
  total_matrix[1, ]
)

names(total_pathway_cpm) <- colnames(
  total_matrix
)

# Taxonomically stratified contributor rows.
target_stratified_rows <- stratified |>
  dplyr::filter(
    is_target_pathway(feature),
    stringr::str_detect(
      feature,
      stringr::fixed("|")
    )
  )

if (nrow(target_stratified_rows) == 0) {
  candidate_rows <- stratified |>
    dplyr::filter(
      stringr::str_detect(
        feature,
        stringr::regex(
          "5265",
          ignore_case = TRUE
        )
      )
    ) |>
    dplyr::pull(feature)
  
  stop(
    "No taxonomically stratified PWY-5265 rows were found. ",
    "Candidate rows: ",
    paste(candidate_rows, collapse = " ; ")
  )
}

contributor_map <- target_stratified_rows |>
  dplyr::transmute(
    original_feature = feature,
    contributor_raw =
      stringr::str_remove(
        feature,
        "^[^|]+\\|"
      ),
    taxon_label =
      taxon_display_label(
        contributor_raw
      )
  ) |>
  dplyr::mutate(
    feature_id =
      sprintf(
        "taxon_%03d",
        dplyr::row_number()
      ),
    .before = 1
  )

if (anyDuplicated(contributor_map$contributor_raw) > 0) {
  warning(
    "Duplicated contributor labels detected. ",
    "They will be retained as separate HUMAnN rows."
  )
}

contributor_matrix <- matrix_from_humann_rows(
  target_stratified_rows,
  row_ids = contributor_map$feature_id
)

# Reorder samples to metadata order.
contributor_matrix <- contributor_matrix[
  ,
  metadata$sample_id,
  drop = FALSE
]

total_pathway_cpm <- total_pathway_cpm[
  metadata$sample_id
]

readr::write_csv(
  contributor_map,
  file.path(
    OUTPUT_DIR,
    "PWY5265_contributor_feature_map.csv"
  )
)

# 6. QC: compare stratified and unstratified PWY-5265 abundance
sum_stratified_cpm <- colSums(
  contributor_matrix
)

qc_by_sample <- tibble::tibble(
  sample_id = metadata$sample_id,
  total_unstratified_cpm =
    as.numeric(total_pathway_cpm),
  summed_stratified_cpm =
    as.numeric(
      sum_stratified_cpm[
        metadata$sample_id
      ]
    )
) |>
  dplyr::mutate(
    difference_cpm =
      summed_stratified_cpm -
      total_unstratified_cpm,
    absolute_difference_cpm =
      abs(difference_cpm),
    relative_difference_percent =
      ifelse(
        total_unstratified_cpm > 0,
        100 *
          difference_cpm /
          total_unstratified_cpm,
        NA_real_
      )
  )

readr::write_csv(
  qc_by_sample,
  file.path(
    OUTPUT_DIR,
    "PWY5265_stratified_vs_total_QC_by_sample.csv"
  )
)

qc_summary <- tibble::tibble(
  target_pathway =
    TARGET_PATHWAY_FULL_NAME,
  n_samples =
    nrow(qc_by_sample),
  n_stratified_contributors =
    nrow(contributor_matrix),
  max_absolute_difference_cpm =
    max(
      qc_by_sample$absolute_difference_cpm
    ),
  mean_absolute_difference_cpm =
    mean(
      qc_by_sample$absolute_difference_cpm
    ),
  median_absolute_difference_cpm =
    median(
      qc_by_sample$absolute_difference_cpm
    )
)

readr::write_csv(
  qc_summary,
  file.path(
    OUTPUT_DIR,
    "PWY5265_stratified_vs_total_QC_summary.csv"
  )
)

# 7. Build sample-level long-format contributor dataset
contributor_wide <- as.data.frame(
  t(contributor_matrix),
  check.names = FALSE
) |>
  tibble::rownames_to_column(
    "sample_id"
  ) |>
  tibble::as_tibble(
    .name_repair = "minimal"
  )

contributor_long <- contributor_wide |>
  tidyr::pivot_longer(
    cols = -sample_id,
    names_to = "feature_id",
    values_to = "cpm"
  ) |>
  dplyr::left_join(
    contributor_map |>
      dplyr::select(
        feature_id,
        contributor_raw,
        taxon_label
      ),
    by = "feature_id"
  ) |>
  dplyr::left_join(
    metadata,
    by = "sample_id"
  ) |>
  dplyr::mutate(
    total_pathway_cpm =
      .env$total_pathway_cpm[
        match(
          sample_id,
          names(.env$total_pathway_cpm)
        )
      ],
    contribution_percent =
      ifelse(
        total_pathway_cpm > 0,
        100 *
          cpm /
          total_pathway_cpm,
        NA_real_
      )
  )

readr::write_csv(
  contributor_long,
  file.path(
    OUTPUT_DIR,
    "PWY5265_species_abundance_by_sample.csv"
  )
)

total_by_sample <- metadata |>
  dplyr::mutate(
    total_PWY5265_cpm =
      as.numeric(
        .env$total_pathway_cpm[
          sample_id
        ]
      )
  )

readr::write_csv(
  total_by_sample,
  file.path(
    OUTPUT_DIR,
    "PWY5265_total_abundance_by_sample.csv"
  )
)

# 8. Descriptive contributor summaries

overall_contributor_summary <- contributor_long |>
  dplyr::group_by(
    feature_id,
    contributor_raw,
    taxon_label
  ) |>
  dplyr::summarise(
    n_samples = dplyr::n(),
    n_detected =
      sum(
        cpm >= DETECTION_CPM
      ),
    prevalence =
      mean(
        cpm >= DETECTION_CPM
      ),
    mean_cpm =
      mean(cpm),
    median_cpm =
      median(cpm),
    sum_contributor_cpm =
      sum(cpm),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    stratified_share_percent =
      if (sum(sum_contributor_cpm) > 0) {
        100 * sum_contributor_cpm / sum(sum_contributor_cpm)
      } else {
        rep(NA_real_, dplyr::n())
      }
  ) |>
  dplyr::arrange(
    dplyr::desc(
      stratified_share_percent
    )
  ) |>
  dplyr::mutate(
    rank_overall = dplyr::row_number(),
    .before = 1
  )

readr::write_csv(
  overall_contributor_summary,
  file.path(
    OUTPUT_DIR,
    "PWY5265_top_contributors_overall.csv"
  )
)

group_summaries <- dplyr::bind_rows(
  summarize_contributors_by_group(
    contributor_long,
    "study_group"
  ),
  summarize_contributors_by_group(
    contributor_long,
    "severe_3cat"
  ),
  summarize_contributors_by_group(
    contributor_long,
    "ethnicity_analysis"
  )
)

readr::write_csv(
  group_summaries,
  file.path(
    OUTPUT_DIR,
    "PWY5265_contributors_by_group.csv"
  )
)

# 9. Compare taxon-specific contributor abundance between groups

# HUMAnN pathway abundance is non-additive across taxonomic strata.
# Therefore, these taxon-specific differences are shown alongside the
# community-level pathway difference but are not expressed as percentages
# of the total community-level difference.

taxon_group_differences <- dplyr::bind_rows(
  compare_taxon_specific_group_means(
    contributor_long,
    group_var = "study_group",
    comparison_level = "PE",
    reference_level = "Control"
  ) |>
    dplyr::mutate(
      contrast = "PE_vs_Control",
      .before = 1
    ),
  
  compare_taxon_specific_group_means(
    contributor_long,
    group_var = "severe_3cat",
    comparison_level = "PE-NSF",
    reference_level = "Control"
  ) |>
    dplyr::mutate(
      contrast = "PE-NSF_vs_Control",
      .before = 1
    ),
  
  compare_taxon_specific_group_means(
    contributor_long,
    group_var = "severe_3cat",
    comparison_level = "PE-SF",
    reference_level = "Control"
  ) |>
    dplyr::mutate(
      contrast = "PE-SF_vs_Control",
      .before = 1
    )
) |>
  dplyr::group_by(contrast) |>
  dplyr::arrange(
    dplyr::desc(
      delta_taxon_specific_cpm
    ),
    .by_group = TRUE
  ) |>
  dplyr::mutate(
    rank_positive_delta =
      dplyr::row_number()
  ) |>
  dplyr::ungroup()

readr::write_csv(
  taxon_group_differences,
  file.path(
    OUTPUT_DIR,
    "PWY5265_taxon_specific_group_differences.csv"
  )
)

# 10. Figures: total pathway and taxonomic contributors

plot_total_status <- total_by_sample |>
  dplyr::filter(
    !is.na(study_group)
  ) |>
  dplyr::mutate(
    log2_total_PWY5265 =
      log2(
        total_PWY5265_cpm + 1
      )
  ) |>
  ggplot2::ggplot(
    ggplot2::aes(
      x = study_group,
      y = log2_total_PWY5265
    )
  ) +
  ggplot2::geom_boxplot(
    outlier.shape = NA
  ) +
  ggplot2::geom_jitter(
    width = 0.15,
    alpha = 0.55,
    size = 1.8
  ) +
  ggplot2::labs(
    x = NULL,
    y = "log2(PWY-5265 CPM + 1)",
    title =
      "Peptidoglycan biosynthesis II by preeclampsia status"
  ) +
  ggplot2::theme_classic(
    base_size = 12
  )

ggplot2::ggsave(
  file.path(
    OUTPUT_DIR,
    "PWY5265_total_abundance_by_status.pdf"
  ),
  plot_total_status,
  width = 5.5,
  height = 5
)

ggplot2::ggsave(
  file.path(
    OUTPUT_DIR,
    "PWY5265_total_abundance_by_status.png"
  ),
  plot_total_status,
  width = 5.5,
  height = 5,
  dpi = 300
)

plot_total_severity <- total_by_sample |>
  dplyr::filter(
    !is.na(severe_3cat)
  ) |>
  dplyr::mutate(
    log2_total_PWY5265 =
      log2(
        total_PWY5265_cpm + 1
      )
  ) |>
  ggplot2::ggplot(
    ggplot2::aes(
      x = severe_3cat,
      y = log2_total_PWY5265
    )
  ) +
  ggplot2::geom_boxplot(
    outlier.shape = NA
  ) +
  ggplot2::geom_jitter(
    width = 0.15,
    alpha = 0.55,
    size = 1.8
  ) +
  ggplot2::labs(
    x = NULL,
    y = "log2(PWY-5265 CPM + 1)",
    title =
      "Peptidoglycan biosynthesis II by preeclampsia severity"
  ) +
  ggplot2::theme_classic(
    base_size = 12
  )

ggplot2::ggsave(
  file.path(
    OUTPUT_DIR,
    "PWY5265_total_abundance_by_severity.pdf"
  ),
  plot_total_severity,
  width = 6,
  height = 5
)

ggplot2::ggsave(
  file.path(
    OUTPUT_DIR,
    "PWY5265_total_abundance_by_severity.png"
  ),
  plot_total_severity,
  width = 6,
  height = 5,
  dpi = 300
)

top_feature_ids <- overall_contributor_summary |>
  dplyr::slice_head(
    n = TOP_N_CONTRIBUTORS
  ) |>
  dplyr::pull(feature_id)

make_contributor_plot_data <- function(grouping_name) {
  
  group_summaries |>
    dplyr::filter(
      grouping_variable == grouping_name,
      !is.na(group)
    ) |>
    dplyr::mutate(
      contributor_plot =
        ifelse(
          feature_id %in%
            top_feature_ids,
          taxon_label,
          "Other contributors"
        )
    ) |>
    dplyr::group_by(
      group,
      contributor_plot
    ) |>
    dplyr::summarise(
      stratified_share_percent =
        sum(
          stratified_share_percent,
          na.rm = TRUE
        ),
      .groups = "drop"
    )
}

stacked_status <- make_contributor_plot_data(
  "study_group"
) |>
  dplyr::mutate(
    group = factor(
      group,
      levels = c(
        "Control",
        "PE"
      )
    )
  )

plot_contributors_status <- ggplot2::ggplot(
  stacked_status,
  ggplot2::aes(
    x = group,
    y = stratified_share_percent,
    fill = contributor_plot
  )
) +
  ggplot2::geom_col() +
  ggplot2::labs(
    x = NULL,
    y = "% of summed stratified PWY-5265 contributor abundance",
    fill = "Contributor",
    title =
      "Taxonomic contributors to peptidoglycan biosynthesis II"
  ) +
  ggplot2::theme_classic(
    base_size = 12
  ) +
  ggplot2::theme(
    legend.position = "right"
  )

ggplot2::ggsave(
  file.path(
    OUTPUT_DIR,
    "PWY5265_taxonomic_contributors_by_status.pdf"
  ),
  plot_contributors_status,
  width = 9,
  height = 5.5
)

ggplot2::ggsave(
  file.path(
    OUTPUT_DIR,
    "PWY5265_taxonomic_contributors_by_status.png"
  ),
  plot_contributors_status,
  width = 9,
  height = 5.5,
  dpi = 300
)

stacked_severity <- make_contributor_plot_data(
  "severe_3cat"
) |>
  dplyr::mutate(
    group = factor(
      group,
      levels = c(
        "Control",
        "PE-NSF",
        "PE-SF"
      )
    )
  )

plot_contributors_severity <- ggplot2::ggplot(
  stacked_severity,
  ggplot2::aes(
    x = group,
    y = stratified_share_percent,
    fill = contributor_plot
  )
) +
  ggplot2::geom_col() +
  ggplot2::labs(
    x = NULL,
    y = "% of summed stratified PWY-5265 contributor abundance",
    fill = "Contributor",
    title =
      "Taxonomic contributors to peptidoglycan biosynthesis II"
  ) +
  ggplot2::theme_classic(
    base_size = 12
  ) +
  ggplot2::theme(
    legend.position = "right"
  )

ggplot2::ggsave(
  file.path(
    OUTPUT_DIR,
    "PWY5265_taxonomic_contributors_by_severity.pdf"
  ),
  plot_contributors_severity,
  width = 9,
  height = 5.5
)

ggplot2::ggsave(
  file.path(
    OUTPUT_DIR,
    "PWY5265_taxonomic_contributors_by_severity.png"
  ),
  plot_contributors_severity,
  width = 9,
  height = 5.5,
  dpi = 300
)

# 11. Define MaAsLin2 analytical populations

minimal_covariates <- c(
  "ethnicity",
  "bmi",
  "antibiotics_during_pregnancy"
)

full_covariates <- c(
  "ethnicity",
  "bmi",
  "antibiotics_during_pregnancy",
  "gestational_age",
  "gest_diabetes",
  "meds_antihypertensives"
)

severity_primary_covariates <- c(
  "ethnicity",
  "bmi",
  "antibiotics_during_pregnancy"
)

severity_ga_covariates <- c(
  severity_primary_covariates,
  "gestational_age"
)

ethnicity_primary_covariates <- c(
  "bmi",
  "antibiotics_during_pregnancy"
)

ethnicity_ga_covariates <- c(
  ethnicity_primary_covariates,
  "gestational_age"
)

analysis_definitions <- list(
  
  overall_status_minimal = list(
    analysis_family = "PE_vs_Control",
    model_type = "Prespecified minimal",
    population = "All participants",
    description =
      "PE versus Control, all participants, prespecified minimal model",
    subset_function = function(x) {
      rep(TRUE, nrow(x))
    },
    exposure = "study_group",
    covariates = minimal_covariates,
    references = c(
      "study_group,Control",
      "ethnicity,Not Hispanic/Latina",
      "antibiotics_during_pregnancy,No"
    ),
    exploratory = TRUE
  ),
  
  overall_status_full = list(
    analysis_family = "PE_vs_Control",
    model_type = "Fully adjusted",
    population = "All participants",
    description =
      "PE versus Control, all participants, fully adjusted model",
    subset_function = function(x) {
      rep(TRUE, nrow(x))
    },
    exposure = "study_group",
    covariates = full_covariates,
    references = c(
      "study_group,Control",
      "ethnicity,Not Hispanic/Latina",
      "antibiotics_during_pregnancy,No",
      "gest_diabetes,No",
      "meds_antihypertensives,No"
    ),
    exploratory = TRUE
  ),
  
  overall_status_term_minimal = list(
    analysis_family = "PE_vs_Control",
    model_type = "Prespecified minimal",
    population = "Term deliveries only",
    description =
      "PE versus Control, term deliveries only, prespecified minimal model",
    subset_function = function(x) {
      !is.na(x$gestational_age) &
        x$gestational_age >= 37
    },
    exposure = "study_group",
    covariates = minimal_covariates,
    references = c(
      "study_group,Control",
      "ethnicity,Not Hispanic/Latina",
      "antibiotics_during_pregnancy,No"
    ),
    exploratory = TRUE
  ),
  
  overall_status_term_full = list(
    analysis_family = "PE_vs_Control",
    model_type = "Fully adjusted",
    population = "Term deliveries only",
    description =
      "PE versus Control, term deliveries only, fully adjusted model",
    subset_function = function(x) {
      !is.na(x$gestational_age) &
        x$gestational_age >= 37
    },
    exposure = "study_group",
    covariates = full_covariates,
    references = c(
      "study_group,Control",
      "ethnicity,Not Hispanic/Latina",
      "antibiotics_during_pregnancy,No",
      "gest_diabetes,No",
      "meds_antihypertensives,No"
    ),
    exploratory = TRUE
  ),
  
  overall_severity_3cat_primary = list(
    analysis_family = "Severity_3cat",
    model_type = "Parsimonious primary",
    population = "All participants",
    description =
      paste(
        "Control versus PE without severe features versus",
        "PE with severe features, parsimonious primary model"
      ),
    subset_function = function(x) {
      rep(TRUE, nrow(x))
    },
    exposure = "severe_3cat",
    covariates = severity_primary_covariates,
    references = c(
      "severe_3cat,Control",
      "ethnicity,Not Hispanic/Latina",
      "antibiotics_during_pregnancy,No"
    ),
    exploratory = TRUE
  ),
  
  overall_severity_3cat_ga = list(
    analysis_family = "Severity_3cat",
    model_type = "Gestational-age sensitivity",
    population = "All participants",
    description =
      paste(
        "Control versus PE without severe features versus",
        "PE with severe features, gestational-age sensitivity model"
      ),
    subset_function = function(x) {
      rep(TRUE, nrow(x))
    },
    exposure = "severe_3cat",
    covariates = severity_ga_covariates,
    references = c(
      "severe_3cat,Control",
      "ethnicity,Not Hispanic/Latina",
      "antibiotics_during_pregnancy,No"
    ),
    exploratory = TRUE
  ),
  
  overall_severity_3cat_term_primary = list(
    analysis_family = "Severity_3cat",
    model_type = "Parsimonious primary",
    population = "Term deliveries only",
    description =
      paste(
        "Control versus PE without severe features versus",
        "PE with severe features, term deliveries only,",
        "parsimonious primary model"
      ),
    subset_function = function(x) {
      !is.na(x$gestational_age) &
        x$gestational_age >= 37
    },
    exposure = "severe_3cat",
    covariates = severity_primary_covariates,
    references = c(
      "severe_3cat,Control",
      "ethnicity,Not Hispanic/Latina",
      "antibiotics_during_pregnancy,No"
    ),
    exploratory = TRUE
  ),
  
  overall_severity_3cat_term_ga = list(
    analysis_family = "Severity_3cat",
    model_type = "Gestational-age sensitivity",
    population = "Term deliveries only",
    description =
      paste(
        "Control versus PE without severe features versus",
        "PE with severe features, term deliveries only,",
        "gestational-age sensitivity model"
      ),
    subset_function = function(x) {
      !is.na(x$gestational_age) &
        x$gestational_age >= 37
    },
    exposure = "severe_3cat",
    covariates = severity_ga_covariates,
    references = c(
      "severe_3cat,Control",
      "ethnicity,Not Hispanic/Latina",
      "antibiotics_during_pregnancy,No"
    ),
    exploratory = TRUE
  ),
  
  pe_severity_primary = list(
    analysis_family = "PE_severity",
    model_type = "Parsimonious primary",
    population = "PE participants",
    description =
      "Severe versus non-severe PE, parsimonious primary model",
    subset_function = function(x) {
      x$study_group == "PE"
    },
    exposure = "severe",
    covariates = severity_primary_covariates,
    references = c(
      "severe,No",
      "ethnicity,Not Hispanic/Latina",
      "antibiotics_during_pregnancy,No"
    ),
    exploratory = TRUE
  ),
  
  pe_severity_ga = list(
    analysis_family = "PE_severity",
    model_type = "Gestational-age sensitivity",
    population = "PE participants",
    description =
      "Severe versus non-severe PE, gestational-age sensitivity model",
    subset_function = function(x) {
      x$study_group == "PE"
    },
    exposure = "severe",
    covariates = severity_ga_covariates,
    references = c(
      "severe,No",
      "ethnicity,Not Hispanic/Latina",
      "antibiotics_during_pregnancy,No"
    ),
    exploratory = TRUE
  ),
  
  pe_severity_term_primary = list(
    analysis_family = "PE_severity",
    model_type = "Parsimonious primary",
    population = "Term PE participants",
    description =
      paste(
        "Severe versus non-severe PE among term PE participants,",
        "parsimonious primary model"
      ),
    subset_function = function(x) {
      x$study_group == "PE" &
        !is.na(x$gestational_age) &
        x$gestational_age >= 37
    },
    exposure = "severe",
    covariates = severity_primary_covariates,
    references = c(
      "severe,No",
      "ethnicity,Not Hispanic/Latina",
      "antibiotics_during_pregnancy,No"
    ),
    exploratory = TRUE
  ),
  
  pe_severity_term_ga = list(
    analysis_family = "PE_severity",
    model_type = "Gestational-age sensitivity",
    population = "Term PE participants",
    description =
      paste(
        "Severe versus non-severe PE among term PE participants,",
        "gestational-age sensitivity model"
      ),
    subset_function = function(x) {
      x$study_group == "PE" &
        !is.na(x$gestational_age) &
        x$gestational_age >= 37
    },
    exposure = "severe",
    covariates = severity_ga_covariates,
    references = c(
      "severe,No",
      "ethnicity,Not Hispanic/Latina",
      "antibiotics_during_pregnancy,No"
    ),
    exploratory = TRUE
  ),
  
  hispanic_status_primary = list(
    analysis_family = "Ethnicity_stratified_status",
    model_type = "Parsimonious primary",
    population = "Hispanic participants",
    description =
      "PE versus Control among Hispanic participants, primary model",
    subset_function = function(x) {
      x$ethnicity_analysis == "Hispanic"
    },
    exposure = "study_group",
    covariates = ethnicity_primary_covariates,
    references = c(
      "study_group,Control",
      "antibiotics_during_pregnancy,No"
    ),
    exploratory = TRUE
  ),
  
  hispanic_status_ga = list(
    analysis_family = "Ethnicity_stratified_status",
    model_type = "Gestational-age sensitivity",
    population = "Hispanic participants",
    description =
      "PE versus Control among Hispanic participants, GA sensitivity",
    subset_function = function(x) {
      x$ethnicity_analysis == "Hispanic"
    },
    exposure = "study_group",
    covariates = ethnicity_ga_covariates,
    references = c(
      "study_group,Control",
      "antibiotics_during_pregnancy,No"
    ),
    exploratory = TRUE
  ),
  
  white_status_primary = list(
    analysis_family = "Ethnicity_stratified_status",
    model_type = "Parsimonious primary",
    population = "Non-Hispanic White participants",
    description =
      "PE versus Control among non-Hispanic White participants, primary model",
    subset_function = function(x) {
      x$ethnicity_analysis == "White"
    },
    exposure = "study_group",
    covariates = ethnicity_primary_covariates,
    references = c(
      "study_group,Control",
      "antibiotics_during_pregnancy,No"
    ),
    exploratory = TRUE
  ),
  
  white_status_ga = list(
    analysis_family = "Ethnicity_stratified_status",
    model_type = "Gestational-age sensitivity",
    population = "Non-Hispanic White participants",
    description =
      "PE versus Control among non-Hispanic White participants, GA sensitivity",
    subset_function = function(x) {
      x$ethnicity_analysis == "White"
    },
    exposure = "study_group",
    covariates = ethnicity_ga_covariates,
    references = c(
      "study_group,Control",
      "antibiotics_during_pregnancy,No"
    ),
    exploratory = TRUE
  ),
  
  hispanic_pe_severity_primary = list(
    analysis_family = "Ethnicity_stratified_severity",
    model_type = "Parsimonious primary",
    population = "Hispanic PE participants",
    description =
      "Severe versus non-severe PE among Hispanic participants, primary model",
    subset_function = function(x) {
      x$ethnicity_analysis == "Hispanic" &
        x$study_group == "PE"
    },
    exposure = "severe",
    covariates = ethnicity_primary_covariates,
    references = c(
      "severe,No",
      "antibiotics_during_pregnancy,No"
    ),
    exploratory = TRUE
  ),
  
  hispanic_pe_severity_ga = list(
    analysis_family = "Ethnicity_stratified_severity",
    model_type = "Gestational-age sensitivity",
    population = "Hispanic PE participants",
    description =
      "Severe versus non-severe PE among Hispanic participants, GA sensitivity",
    subset_function = function(x) {
      x$ethnicity_analysis == "Hispanic" &
        x$study_group == "PE"
    },
    exposure = "severe",
    covariates = ethnicity_ga_covariates,
    references = c(
      "severe,No",
      "antibiotics_during_pregnancy,No"
    ),
    exploratory = TRUE
  ),
  
  white_pe_severity_primary = list(
    analysis_family = "Ethnicity_stratified_severity",
    model_type = "Parsimonious primary",
    population = "Non-Hispanic White PE participants",
    description =
      paste(
        "Severe versus non-severe PE among non-Hispanic White participants,",
        "primary model"
      ),
    subset_function = function(x) {
      x$ethnicity_analysis == "White" &
        x$study_group == "PE"
    },
    exposure = "severe",
    covariates = ethnicity_primary_covariates,
    references = c(
      "severe,No",
      "antibiotics_during_pregnancy,No"
    ),
    exploratory = TRUE
  ),
  
  white_pe_severity_ga = list(
    analysis_family = "Ethnicity_stratified_severity",
    model_type = "Gestational-age sensitivity",
    population = "Non-Hispanic White PE participants",
    description =
      paste(
        "Severe versus non-severe PE among non-Hispanic White participants,",
        "GA sensitivity"
      ),
    subset_function = function(x) {
      x$ethnicity_analysis == "White" &
        x$study_group == "PE"
    },
    exposure = "severe",
    covariates = ethnicity_ga_covariates,
    references = c(
      "severe,No",
      "antibiotics_during_pregnancy,No"
    ),
    exploratory = TRUE
  )
)

# 12. Run targeted MaAsLin2 models on taxon-specific PWY-5265 abundance

variable_has_variation <- function(x) {
  
  x <- x[!is.na(x)]
  
  if (length(x) == 0) {
    return(FALSE)
  }
  
  if (is.numeric(x)) {
    return(
      length(unique(x)) > 1 &&
        is.finite(stats::var(x)) &&
        stats::var(x) > 0
    )
  }
  
  length(unique(as.character(x))) > 1
}

run_contributor_model <- function(
    definition,
    analysis_name
) {
  
  message("")
  message("Running: ", analysis_name)
  message(definition$description)
  
  analysis_dir <- file.path(
    OUTPUT_DIR,
    "MaAsLin2_contributor_models",
    analysis_name
  )
  
  dir.create(
    analysis_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  subset_index <-
    definition$subset_function(
      metadata
    )
  
  subset_index[
    is.na(subset_index)
  ] <- FALSE
  
  requested_model_variables <- c(
    definition$exposure,
    definition$covariates
  )
  
  model_metadata <- metadata[
    subset_index,
    ,
    drop = FALSE
  ] |>
    dplyr::filter(
      dplyr::if_all(
        dplyr::all_of(
          requested_model_variables
        ),
        ~ !is.na(.x)
      )
    ) |>
    droplevels()
  
  if (
    nrow(model_metadata) <
    MIN_MODEL_N
  ) {
    warning(
      analysis_name,
      ": fewer than ",
      MIN_MODEL_N,
      " complete cases. Skipping."
    )
    
    return(
      list(
        analysis = analysis_name,
        status = "skipped_small_n",
        results = NULL,
        filter_summary = NULL,
        manifest = NULL
      )
    )
  }
  
  exposure_counts <- model_metadata |>
    dplyr::count(
      .data[[
        definition$exposure
      ]],
      name = "n"
    )
  
  names(exposure_counts)[1] <-
    "exposure_level"
  
  readr::write_csv(
    exposure_counts,
    file.path(
      analysis_dir,
      "exposure_counts.csv"
    )
  )
  
  if (
    nrow(exposure_counts) < 2 ||
    any(
      exposure_counts$n <
      MIN_LEVEL_N
    )
  ) {
    warning(
      analysis_name,
      ": insufficient exposure group size. Skipping."
    )
    
    return(
      list(
        analysis = analysis_name,
        status = "skipped_small_group",
        results = NULL,
        filter_summary = NULL,
        manifest = NULL
      )
    )
  }
  
  varying_covariates <- definition$covariates[
    vapply(
      definition$covariates,
      function(v) {
        variable_has_variation(
          model_metadata[[v]]
        )
      },
      logical(1)
    )
  ]
  
  dropped_invariant_covariates <- setdiff(
    definition$covariates,
    varying_covariates
  )
  
  if (
    length(
      dropped_invariant_covariates
    ) > 0
  ) {
    warning(
      analysis_name,
      ": dropping covariate(s) with no variation in this subset: ",
      paste(
        dropped_invariant_covariates,
        collapse = ", "
      )
    )
  }
  
  model_variables <- c(
    definition$exposure,
    varying_covariates
  )
  
  reference_variables <- sub(
    ",.*$",
    "",
    definition$references
  )
  
  effective_references <- definition$references[
    reference_variables %in%
      model_variables
  ]
  
  model_manifest <- tibble::tibble(
    analysis = analysis_name,
    analysis_family =
      definition$analysis_family,
    model_type =
      definition$model_type,
    population =
      definition$population,
    description =
      definition$description,
    exposure =
      definition$exposure,
    requested_covariates =
      paste(
        definition$covariates,
        collapse = " + "
      ),
    covariates_used =
      paste(
        varying_covariates,
        collapse = " + "
      ),
    dropped_invariant_covariates =
      paste(
        dropped_invariant_covariates,
        collapse = "; "
      ),
    model_formula =
      paste(
        "taxon_specific_PWY5265_abundance ~",
        paste(
          model_variables,
          collapse = " + "
        )
      ),
    n_complete_cases =
      nrow(model_metadata)
  )
  
  readr::write_csv(
    model_manifest,
    file.path(
      analysis_dir,
      "model_manifest.csv"
    )
  )
  
  sample_ids <- model_metadata$sample_id
  
  model_matrix <- contributor_matrix[
    ,
    sample_ids,
    drop = FALSE
  ]
  
  filter_stats <- tibble::tibble(
    feature_id =
      rownames(model_matrix),
    n_samples =
      ncol(model_matrix),
    n_detected =
      rowSums(
        model_matrix >=
          DETECTION_CPM
      ),
    prevalence =
      rowMeans(
        model_matrix >=
          DETECTION_CPM
      ),
    mean_cpm =
      rowMeans(
        model_matrix
      ),
    variance =
      apply(
        model_matrix,
        1,
        stats::var
      )
  ) |>
    dplyr::left_join(
      contributor_map |>
        dplyr::select(
          feature_id,
          contributor_raw,
          taxon_label
        ),
      by = "feature_id"
    ) |>
    dplyr::mutate(
      passes_filter =
        n_detected >=
        MIN_DETECTED_SAMPLES &
        prevalence >=
        MIN_PREVALENCE &
        is.finite(variance) &
        variance > 0
    )
  
  readr::write_csv(
    filter_stats,
    file.path(
      analysis_dir,
      "contributor_prefilter_statistics.csv"
    )
  )
  
  kept_features <- filter_stats |>
    dplyr::filter(
      passes_filter
    ) |>
    dplyr::pull(feature_id)
  
  filter_summary <- tibble::tibble(
    analysis =
      analysis_name,
    analysis_family =
      definition$analysis_family,
    model_type =
      definition$model_type,
    population =
      definition$population,
    n_complete_cases =
      nrow(model_metadata),
    n_contributors_total =
      nrow(model_matrix),
    n_contributors_tested =
      length(kept_features)
  )
  
  readr::write_csv(
    filter_summary,
    file.path(
      analysis_dir,
      "contributor_prefilter_summary.csv"
    )
  )
  
  if (
    length(kept_features) == 0
  ) {
    warning(
      analysis_name,
      ": no contributor features passed filtering."
    )
    
    return(
      list(
        analysis = analysis_name,
        status = "skipped_no_features",
        results = NULL,
        filter_summary = filter_summary,
        manifest = model_manifest
      )
    )
  }
  
  filtered_matrix <- model_matrix[
    kept_features,
    ,
    drop = FALSE
  ]
  
  input_data <- as.data.frame(
    t(filtered_matrix),
    check.names = FALSE
  )
  
  input_metadata <- model_metadata |>
    dplyr::select(
      sample_id,
      dplyr::all_of(
        model_variables
      )
    ) |>
    as.data.frame(
      check.names = FALSE
    )
  
  rownames(input_data) <-
    model_metadata$sample_id
  
  rownames(input_metadata) <-
    model_metadata$sample_id
  
  input_metadata$sample_id <-
    NULL
  
  maaslin_output <- file.path(
    analysis_dir,
    "Maaslin2_output"
  )
  
  if (
    dir.exists(
      maaslin_output
    )
  ) {
    unlink(
      maaslin_output,
      recursive = TRUE,
      force = TRUE
    )
  }
  
  Maaslin2::Maaslin2(
    input_data =
      input_data,
    input_metadata =
      input_metadata,
    output =
      maaslin_output,
    fixed_effects =
      model_variables,
    random_effects =
      NULL,
    reference =
      effective_references,
    normalization =
      "NONE",
    transform =
      "LOG",
    analysis_method =
      "LM",
    standardize =
      FALSE,
    min_abundance =
      0,
    min_prevalence =
      0,
    min_variance =
      0,
    max_significance =
      EXPLORATORY_FDR,
    correction =
      "BH",
    cores =
      N_CORES,
    plot_heatmap =
      FALSE,
    plot_scatter =
      FALSE
  )
  
  results_file <- file.path(
    maaslin_output,
    "all_results.tsv"
  )
  
  if (!file.exists(results_file)) {
    stop(
      "MaAsLin2 result file not found for ",
      analysis_name
    )
  }
  
  all_results <- readr::read_tsv(
    results_file,
    show_col_types = FALSE,
    name_repair = "minimal"
  )
  
  required_result_columns <- c(
    "feature",
    "metadata",
    "value",
    "coef",
    "stderr",
    "pval",
    "qval"
  )
  
  missing_result_columns <- setdiff(
    required_result_columns,
    names(all_results)
  )
  
  if (
    length(
      missing_result_columns
    ) > 0
  ) {
    stop(
      "Unexpected MaAsLin2 output columns for ",
      analysis_name,
      ": missing ",
      paste(
        missing_result_columns,
        collapse = ", "
      )
    )
  }
  
  exposure_results <- all_results |>
    dplyr::filter(
      metadata ==
        definition$exposure
    ) |>
    dplyr::group_by(
      value
    ) |>
    dplyr::mutate(
      q_exposure =
        stats::p.adjust(
          pval,
          method = "BH"
        )
    ) |>
    dplyr::ungroup() |>
    dplyr::rename(
      feature_id = feature
    ) |>
    dplyr::left_join(
      contributor_map |>
        dplyr::select(
          feature_id,
          contributor_raw,
          taxon_label
        ),
      by = "feature_id"
    ) |>
    dplyr::mutate(
      analysis =
        analysis_name,
      analysis_family =
        definition$analysis_family,
      model_type =
        definition$model_type,
      population =
        definition$population,
      analysis_description =
        definition$description,
      exploratory_analysis =
        definition$exploratory,
      n_complete_cases =
        nrow(model_metadata),
      conf_low =
        coef - 1.96 * stderr,
      conf_high =
        coef + 1.96 * stderr,
      direction =
        dplyr::case_when(
          coef > 0 ~
            "Higher in displayed level",
          coef < 0 ~
            "Lower in displayed level",
          TRUE ~
            "No direction"
        ),
      evidence_category =
        dplyr::case_when(
          q_exposure <
            PRIMARY_FDR ~
            "FDR q < 0.05",
          q_exposure <
            EXPLORATORY_FDR ~
            "Exploratory q < 0.25",
          pval <
            0.05 ~
            "Nominal p < 0.05",
          pval <=
            MARGINAL_P ~
            "Marginal p = 0.05-0.08",
          TRUE ~
            "Not prioritized"
        )
    ) |>
    dplyr::arrange(
      value,
      q_exposure,
      pval
    )
  
  readr::write_csv(
    exposure_results,
    file.path(
      analysis_dir,
      "PWY5265_contributor_exposure_results.csv"
    )
  )
  
  prioritized <- exposure_results |>
    dplyr::filter(
      q_exposure <
        EXPLORATORY_FDR |
        pval <=
        MARGINAL_P
    )
  
  readr::write_csv(
    prioritized,
    file.path(
      analysis_dir,
      "PWY5265_contributor_prioritized_results.csv"
    )
  )
  
  list(
    analysis =
      analysis_name,
    status =
      "completed",
    results =
      exposure_results,
    filter_summary =
      filter_summary,
    manifest =
      model_manifest
  )
}

model_outputs <- purrr::imap(
  analysis_definitions,
  run_contributor_model
)

# 13. Combine targeted model results

all_contributor_results <- dplyr::bind_rows(
  purrr::map(
    model_outputs,
    "results"
  )
)

if (nrow(all_contributor_results) == 0) {
  stop(
    "No contributor MaAsLin2 models produced results. ",
    "Inspect RUN_STATUS.csv and the model-specific folders."
  )
}

all_filter_summaries <- dplyr::bind_rows(
  purrr::map(
    model_outputs,
    "filter_summary"
  )
)

all_model_manifests <- dplyr::bind_rows(
  purrr::map(
    model_outputs,
    "manifest"
  )
)

run_status <- tibble::tibble(
  analysis =
    names(model_outputs),
  status =
    purrr::map_chr(
      model_outputs,
      "status"
    )
)

readr::write_csv(
  run_status,
  file.path(
    OUTPUT_DIR,
    "RUN_STATUS.csv"
  )
)

readr::write_csv(
  all_filter_summaries,
  file.path(
    OUTPUT_DIR,
    "ALL_CONTRIBUTOR_PREFILTER_SUMMARIES.csv"
  )
)

readr::write_csv(
  all_model_manifests,
  file.path(
    OUTPUT_DIR,
    "ALL_MODEL_MANIFESTS.csv"
  )
)

readr::write_csv(
  all_contributor_results,
  file.path(
    OUTPUT_DIR,
    "ALL_PWY5265_CONTRIBUTOR_RESULTS.csv"
  )
)

all_prioritized <- all_contributor_results |>
  dplyr::filter(
    q_exposure <
      EXPLORATORY_FDR |
      pval <=
      MARGINAL_P
  ) |>
  dplyr::arrange(
    analysis,
    value,
    q_exposure,
    pval
  )

readr::write_csv(
  all_prioritized,
  file.path(
    OUTPUT_DIR,
    "ALL_PWY5265_CONTRIBUTOR_PRIORITIZED_RESULTS.csv"
  )
)

all_fdr_results <- all_contributor_results |>
  dplyr::filter(
    q_exposure <
      PRIMARY_FDR
  ) |>
  dplyr::arrange(
    q_exposure,
    pval
  )

readr::write_csv(
  all_fdr_results,
  file.path(
    OUTPUT_DIR,
    "ALL_PWY5265_CONTRIBUTOR_FDR05_RESULTS.csv"
  )
)

all_exploratory_fdr <- all_contributor_results |>
  dplyr::filter(
    q_exposure <
      EXPLORATORY_FDR
  ) |>
  dplyr::arrange(
    q_exposure,
    pval
  )

readr::write_csv(
  all_exploratory_fdr,
  file.path(
    OUTPUT_DIR,
    "ALL_PWY5265_CONTRIBUTOR_FDR25_RESULTS.csv"
  )
)

status_model_results <- all_contributor_results |>
  dplyr::filter(
    analysis_family ==
      "PE_vs_Control",
    value ==
      "PE"
  ) |>
  dplyr::arrange(
    taxon_label,
    factor(
      analysis,
      levels = c(
        "overall_status_minimal",
        "overall_status_full",
        "overall_status_term_minimal",
        "overall_status_term_full"
      )
    )
  )

readr::write_csv(
  status_model_results,
  file.path(
    OUTPUT_DIR,
    "PWY5265_PE_VS_CONTROL_ALL_MODELS.csv"
  )
)

# 14. Focused PE vs Control results

focused_status_results <- all_contributor_results |>
  dplyr::filter(
    analysis %in%
      c(
        "overall_status_minimal",
        "overall_status_full",
        "overall_status_term_minimal",
        "overall_status_term_full"
      ),
    value ==
      "PE"
  ) |>
  dplyr::arrange(
    analysis,
    q_exposure,
    pval
  )

readr::write_csv(
  focused_status_results,
  file.path(
    OUTPUT_DIR,
    "FOCUSED_PWY5265_PE_VS_CONTROL_RESULTS.csv"
  )
)

secondary_results <- all_contributor_results |>
  dplyr::filter(
    !analysis %in%
      c(
        "overall_status_minimal",
        "overall_status_full",
        "overall_status_term_minimal",
        "overall_status_term_full"
      )
  ) |>
  dplyr::arrange(
    analysis,
    value,
    q_exposure,
    pval
  )

readr::write_csv(
  secondary_results,
  file.path(
    OUTPUT_DIR,
    "PWY5265_SECONDARY_SEVERITY_ETHNICITY_RESULTS.csv"
  )
)

# 15. Plot PE vs Control coefficients for contributors prioritized in any model

prioritized_status_taxa <- focused_status_results |>
  dplyr::filter(
    q_exposure <
      EXPLORATORY_FDR |
      pval <=
      MARGINAL_P
  ) |>
  dplyr::distinct(feature_id) |>
  dplyr::pull(feature_id)

if (
  length(
    prioritized_status_taxa
  ) > 0
) {
  
  plot_results <- focused_status_results |>
    dplyr::filter(
      feature_id %in%
        prioritized_status_taxa
    ) |>
    dplyr::mutate(
      model_label =
        dplyr::recode(
          analysis,
          overall_status_minimal =
            "Minimal",
          overall_status_full =
            "Fully adjusted",
          overall_status_term_minimal =
            "Term minimal",
          overall_status_term_full =
            "Term fully adjusted"
        ),
      model_label = factor(
        model_label,
        levels = c(
          "Minimal",
          "Fully adjusted",
          "Term minimal",
          "Term fully adjusted"
        )
      ),
      taxon_label =
        forcats::fct_reorder(
          taxon_label,
          coef
        )
    )
  
  coefficient_plot <- ggplot2::ggplot(
    plot_results,
    ggplot2::aes(
      x = coef,
      y = taxon_label
    )
  ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linetype = 2
    ) +
    ggplot2::geom_errorbarh(
      ggplot2::aes(
        xmin = conf_low,
        xmax = conf_high
      ),
      height = 0
    ) +
    ggplot2::geom_point(
      size = 2.3
    ) +
    ggplot2::facet_wrap(
      ~ model_label,
      scales = "free_y"
    ) +
    ggplot2::labs(
      x =
        "Adjusted MaAsLin2 coefficient after LOG transformation",
      y = NULL,
      title =
        "Taxon-specific contributors to PWY-5265 in preeclampsia"
    ) +
    ggplot2::theme_classic(
      base_size = 11
    )
  
  ggplot2::ggsave(
    file.path(
      OUTPUT_DIR,
      "PWY5265_PE_VS_CONTROL_contributor_coefficients.pdf"
    ),
    coefficient_plot,
    width = 10,
    height = 7
  )
  
  ggplot2::ggsave(
    file.path(
      OUTPUT_DIR,
      "PWY5265_PE_VS_CONTROL_contributor_coefficients.png"
    ),
    coefficient_plot,
    width = 10,
    height = 7,
    dpi = 300
  )
}

# 16. Analysis settings and notes

analysis_settings <- tibble::tibble(
  setting = c(
    "Target pathway",
    "Stratified pathway input",
    "Normalization",
    "MaAsLin2 transform",
    "MaAsLin2 model",
    "Detection threshold CPM",
    "Minimum prevalence",
    "Minimum detected samples",
    "Primary FDR",
    "Exploratory FDR",
    "Marginal p threshold",
    "Maximum cores"
  ),
  value = c(
    TARGET_PATHWAY_FULL_NAME,
    STRATIFIED_PATHWAY_FILE,
    "NONE (HUMAnN CPM already normalized)",
    "LOG",
    "LM",
    DETECTION_CPM,
    MIN_PREVALENCE,
    MIN_DETECTED_SAMPLES,
    PRIMARY_FDR,
    EXPLORATORY_FDR,
    MARGINAL_P,
    N_CORES
  )
)

readr::write_csv(
  analysis_settings,
  file.path(
    OUTPUT_DIR,
    "analysis_settings.csv"
  )
)

notes <- c(
  paste0(
    "Target pathway: ",
    TARGET_PATHWAY_FULL_NAME
  ),
  "",
  paste0(
    "Taxon-specific PWY-5265 abundances were extracted from the ",
    "HUMAnN3 stratified CPM pathway table."
  ),
  "",
  paste0(
    "Overall PE vs Control minimal model: study_group + ethnicity + bmi + ",
    "antibiotics_during_pregnancy."
  ),
  "",
  paste0(
    "Overall PE vs Control full model: study_group + ethnicity + bmi + ",
    "antibiotics_during_pregnancy + gestational_age + gest_diabetes + ",
    "meds_antihypertensives."
  ),
  "",
  paste0(
    "Both PE vs Control models were repeated among term deliveries ",
    "(gestational_age >= 37 weeks)."
  ),
  "",
  paste0(
    "PE severity primary models are parsimonious and adjust for ethnicity, ",
    "BMI, and antibiotics during pregnancy; GA-adjusted sensitivity models ",
    "add gestational_age. Both are also repeated among term deliveries."
  ),
  "",
  paste0(
    "Ethnicity-stratified models adjust for BMI and antibiotics during ",
    "pregnancy, with a separate gestational-age sensitivity model."
  ),
  "",
  "No stepwise variable selection is used.",
  "",
  paste0(
    "If a covariate has no variation within a restricted subgroup, it is ",
    "dropped because it cannot be estimated; this is recorded in ",
    "ALL_MODEL_MANIFESTS.csv."
  ),
  "",
  paste0(
    "For each exposure contrast, q_exposure is BH-adjusted across the ",
    "taxon-specific PWY-5265 contributors tested in that targeted model."
  ),
  "",
  paste0(
    "HUMAnN computes pathway abundance independently at the community and ",
    "taxonomic-stratum levels. Consequently, taxonomic-stratum pathway ",
    "abundances do not necessarily sum to the unstratified community-level ",
    "pathway abundance."
  ),
  "",
  paste0(
    "PWY5265_stratified_vs_total_QC_by_sample.csv documents this ",
    "non-additivity. PWY5265_taxon_specific_group_differences.csv therefore ",
    "reports taxon-specific group differences but does not present them as ",
    "percentages of the total community-level pathway difference."
  ),
  "",
  paste0(
    "Contributor-composition plots express each taxon's share of the summed ",
    "stratified contributor abundance, not its percentage of the unstratified ",
    "community pathway abundance."
  ),
  "",
  paste0(
    "This targeted analysis is exploratory because the parent PWY-5265 ",
    "association was nominal and did not survive FDR correction in the broad ",
    "functional screen."
  ),
  "",
  paste0(
    "These results do not demonstrate increased peptidoglycan production, ",
    "circulating peptidoglycan, or causal activity."
  )
)

writeLines(
  notes,
  file.path(
    OUTPUT_DIR,
    "ANALYSIS_NOTES.txt"
  )
)

capture.output(
  sessionInfo(),
  file =
    file.path(
      OUTPUT_DIR,
      "sessionInfo.txt"
    )
)

# 17. Console summary

message("\nAnalysis completed.")
message("Output directory: ", OUTPUT_DIR)
message(
  "PWY-5265 contributors detected in stratified table: ",
  nrow(contributor_matrix)
)
message("\nStart by inspecting:")
message(
  "  1. RUN_STATUS.csv"
)
message(
  "  2. ALL_MODEL_MANIFESTS.csv"
)
message(
  "  3. PWY5265_PE_VS_CONTROL_ALL_MODELS.csv"
)
message(
  "  4. ALL_PWY5265_CONTRIBUTOR_PRIORITIZED_RESULTS.csv"
)
message(
  "  5. PWY5265_top_contributors_overall.csv"
)
message(
  "  6. PWY5265_taxon_specific_group_differences.csv"
)
message(
  "  7. PWY5265_PE_VS_CONTROL_contributor_coefficients.pdf"
)


# Additional sparse-contributor analyses
#
# Species-level taxon-specific PWY-5265 contributions are extremely sparse in
# this cohort. We therefore do not lower the prespecified 10% prevalence filter
# to force continuous MaAsLin2 models. Instead, this section adds two targeted
# exploratory analyses:
#   1) exact presence/absence comparisons using Fisher's exact test;
#   2) aggregation of species-stratified contributor abundance at the genus
#      level, plus a combined classified-contributor feature.
#
# The genus-level quantities below are sums of HUMAnN species-stratified
# pathway abundances assigned to members of a genus. They are NOT equivalent
# to a de novo HUMAnN genus-level pathway reconstruction.

extract_genus_label <- function(x) {
  dplyr::case_when(
    x == "unclassified" ~ "unclassified",
    stringr::str_detect(x, "g__") ~
      stringr::str_extract(x, "g__[^.|]+") |>
      stringr::str_remove("^g__") |>
      stringr::str_replace_all("_", " "),
    TRUE ~ "other"
  )
}

species_presence_fisher <- function(long_data, term_only = FALSE) {
  d <- long_data
  
  if (term_only) {
    d <- d |>
      dplyr::filter(
        !is.na(gestational_age),
        gestational_age >= 37
      )
    population_label <- "Term deliveries only"
  } else {
    population_label <- "All participants"
  }
  
  d <- d |>
    dplyr::filter(!is.na(study_group)) |>
    dplyr::mutate(detected = cpm >= DETECTION_CPM)
  
  out <- d |>
    dplyr::group_by(feature_id, contributor_raw, taxon_label) |>
    dplyr::group_modify(function(.x, .y) {
      tab <- table(
        factor(.x$study_group, levels = c("Control", "PE")),
        factor(.x$detected, levels = c(FALSE, TRUE))
      )
      
      n_control <- sum(.x$study_group == "Control")
      n_pe <- sum(.x$study_group == "PE")
      det_control <- sum(.x$study_group == "Control" & .x$detected)
      det_pe <- sum(.x$study_group == "PE" & .x$detected)
      
      ft <- stats::fisher.test(tab)
      
      tibble::tibble(
        population = population_label,
        n_control = n_control,
        n_pe = n_pe,
        detected_control = det_control,
        detected_pe = det_pe,
        prevalence_control = det_control / n_control,
        prevalence_pe = det_pe / n_pe,
        odds_ratio = unname(ft$estimate),
        pval = ft$p.value
      )
    }) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      q_fisher = stats::p.adjust(pval, method = "BH")
    ) |>
    dplyr::arrange(pval)
  
  out
}

species_fisher_all <- species_presence_fisher(
  contributor_long,
  term_only = FALSE
)

species_fisher_term <- species_presence_fisher(
  contributor_long,
  term_only = TRUE
)

species_fisher_results <- dplyr::bind_rows(
  species_fisher_all,
  species_fisher_term
)

readr::write_csv(
  species_fisher_results,
  file.path(
    OUTPUT_DIR,
    "PWY5265_species_presence_Fisher_PE_vs_Control.csv"
  )
)

# Aggregate species-stratified contributor abundance by genus.

genus_long <- contributor_long |>
  dplyr::mutate(
    genus_label = extract_genus_label(contributor_raw)
  ) |>
  dplyr::group_by(sample_id, genus_label) |>
  dplyr::summarise(
    cpm = sum(cpm),
    .groups = "drop"
  ) |>
  tidyr::complete(
    sample_id = metadata$sample_id,
    genus_label,
    fill = list(cpm = 0)
  ) |>
  dplyr::left_join(
    metadata,
    by = "sample_id"
  )

# Add a combined feature representing all taxonomically classified contributors.
classified_total <- contributor_long |>
  dplyr::filter(contributor_raw != "unclassified") |>
  dplyr::group_by(sample_id) |>
  dplyr::summarise(
    cpm = sum(cpm),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    genus_label = "All classified contributors"
  ) |>
  dplyr::left_join(
    metadata,
    by = "sample_id"
  )

genus_long_augmented <- dplyr::bind_rows(
  genus_long,
  classified_total
)

readr::write_csv(
  genus_long_augmented |>
    dplyr::select(
      sample_id,
      genus_label,
      cpm,
      study_group,
      severe,
      severe_3cat,
      ethnicity,
      ethnicity_analysis,
      bmi,
      antibiotics_during_pregnancy,
      gestational_age,
      gest_diabetes,
      meds_antihypertensives
    ),
  file.path(
    OUTPUT_DIR,
    "PWY5265_summed_species_stratified_contributors_by_genus.csv"
  )
)

genus_summary <- genus_long_augmented |>
  dplyr::group_by(genus_label) |>
  dplyr::summarise(
    n_samples = dplyr::n(),
    n_detected = sum(cpm >= DETECTION_CPM),
    prevalence = mean(cpm >= DETECTION_CPM),
    mean_cpm = mean(cpm),
    median_cpm = median(cpm),
    sum_cpm = sum(cpm),
    .groups = "drop"
  ) |>
  dplyr::arrange(dplyr::desc(prevalence), dplyr::desc(sum_cpm))

readr::write_csv(
  genus_summary,
  file.path(
    OUTPUT_DIR,
    "PWY5265_genus_aggregated_contributor_summary.csv"
  )
)

# Exact detection analysis for genus-aggregated contributor features.

genus_presence_fisher <- function(long_data, term_only = FALSE) {
  d <- long_data
  
  if (term_only) {
    d <- d |>
      dplyr::filter(
        !is.na(gestational_age),
        gestational_age >= 37
      )
    population_label <- "Term deliveries only"
  } else {
    population_label <- "All participants"
  }
  
  d <- d |>
    dplyr::filter(!is.na(study_group)) |>
    dplyr::mutate(detected = cpm >= DETECTION_CPM)
  
  d |>
    dplyr::group_by(genus_label) |>
    dplyr::group_modify(function(.x, .y) {
      tab <- table(
        factor(.x$study_group, levels = c("Control", "PE")),
        factor(.x$detected, levels = c(FALSE, TRUE))
      )
      
      n_control <- sum(.x$study_group == "Control")
      n_pe <- sum(.x$study_group == "PE")
      det_control <- sum(.x$study_group == "Control" & .x$detected)
      det_pe <- sum(.x$study_group == "PE" & .x$detected)
      ft <- stats::fisher.test(tab)
      
      tibble::tibble(
        population = population_label,
        n_control = n_control,
        n_pe = n_pe,
        detected_control = det_control,
        detected_pe = det_pe,
        prevalence_control = det_control / n_control,
        prevalence_pe = det_pe / n_pe,
        odds_ratio = unname(ft$estimate),
        pval = ft$p.value
      )
    }) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      q_fisher = stats::p.adjust(pval, method = "BH")
    ) |>
    dplyr::arrange(pval)
}

genus_fisher_results <- dplyr::bind_rows(
  genus_presence_fisher(genus_long_augmented, term_only = FALSE),
  genus_presence_fisher(genus_long_augmented, term_only = TRUE)
)

readr::write_csv(
  genus_fisher_results,
  file.path(
    OUTPUT_DIR,
    "PWY5265_genus_presence_Fisher_PE_vs_Control.csv"
  )
)

# Adjusted MaAsLin2 analysis of genus-aggregated contributor features.
# Only the four primary PE-vs-Control models are evaluated here.

primary_status_analysis_names <- c(
  "overall_status_minimal",
  "overall_status_full",
  "overall_status_term_minimal",
  "overall_status_term_full"
)

genus_feature_levels <- unique(genus_long_augmented$genus_label)
genus_feature_ids <- paste0(
  "genus_",
  seq_along(genus_feature_levels)
)
genus_feature_map <- tibble::tibble(
  feature_id = genus_feature_ids,
  genus_label = genus_feature_levels
)

genus_matrix_df <- genus_long_augmented |>
  dplyr::select(sample_id, genus_label, cpm) |>
  dplyr::left_join(genus_feature_map, by = "genus_label") |>
  dplyr::select(sample_id, feature_id, cpm) |>
  tidyr::pivot_wider(
    names_from = sample_id,
    values_from = cpm,
    values_fill = 0
  )

genus_matrix <- as.matrix(
  genus_matrix_df[, -1, drop = FALSE]
)
storage.mode(genus_matrix) <- "numeric"
rownames(genus_matrix) <- genus_matrix_df$feature_id

genus_model_outputs <- list()

for (analysis_name in primary_status_analysis_names) {
  definition <- analysis_definitions[[analysis_name]]
  
  subset_index <- definition$subset_function(metadata)
  subset_index[is.na(subset_index)] <- FALSE
  
  requested_covariates <- definition$covariates
  subset_metadata <- metadata[subset_index, , drop = FALSE]
  
  model_variables <- c(definition$exposure, requested_covariates)
  model_metadata <- subset_metadata |>
    dplyr::select(sample_id, dplyr::all_of(model_variables)) |>
    tidyr::drop_na(dplyr::all_of(model_variables)) |>
    dplyr::mutate(dplyr::across(where(is.factor), droplevels))
  
  analysis_dir <- file.path(
    OUTPUT_DIR,
    "Genus_aggregated_models",
    analysis_name
  )
  dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)
  
  model_matrix <- genus_matrix[, model_metadata$sample_id, drop = FALSE]
  
  prevalence_stats <- tibble::tibble(
    feature_id = rownames(model_matrix),
    n_detected = rowSums(model_matrix >= DETECTION_CPM),
    prevalence = rowMeans(model_matrix >= DETECTION_CPM),
    variance = apply(model_matrix, 1, stats::var)
  ) |>
    dplyr::left_join(genus_feature_map, by = "feature_id") |>
    dplyr::mutate(
      passes_filter =
        n_detected >= MIN_DETECTED_SAMPLES &
        prevalence >= MIN_PREVALENCE &
        variance > 0
    )
  
  readr::write_csv(
    prevalence_stats,
    file.path(analysis_dir, "genus_prefilter_statistics.csv")
  )
  
  kept_features <- prevalence_stats |>
    dplyr::filter(passes_filter) |>
    dplyr::pull(feature_id)
  
  if (length(kept_features) == 0) {
    genus_model_outputs[[analysis_name]] <- tibble::tibble(
      analysis = analysis_name,
      status = "skipped_no_features"
    )
    next
  }
  
  input_data <- as.data.frame(
    t(model_matrix[kept_features, , drop = FALSE]),
    check.names = FALSE
  )
  input_metadata <- as.data.frame(
    model_metadata |>
      dplyr::select(-sample_id),
    check.names = FALSE
  )
  rownames(input_data) <- model_metadata$sample_id
  rownames(input_metadata) <- model_metadata$sample_id
  
  maaslin_output <- file.path(analysis_dir, "Maaslin2_output")
  if (dir.exists(maaslin_output)) {
    unlink(maaslin_output, recursive = TRUE, force = TRUE)
  }
  
  Maaslin2::Maaslin2(
    input_data = input_data,
    input_metadata = input_metadata,
    output = maaslin_output,
    fixed_effects = model_variables,
    random_effects = NULL,
    reference = definition$references,
    normalization = "NONE",
    transform = "LOG",
    analysis_method = "LM",
    standardize = FALSE,
    min_abundance = 0,
    min_prevalence = 0,
    min_variance = 0,
    max_significance = EXPLORATORY_FDR,
    correction = "BH",
    cores = N_CORES
  )
  
  results_file <- file.path(maaslin_output, "all_results.tsv")
  res <- readr::read_tsv(results_file, show_col_types = FALSE) |>
    dplyr::filter(metadata == definition$exposure) |>
    dplyr::rename(feature_id = feature) |>
    dplyr::left_join(genus_feature_map, by = "feature_id") |>
    dplyr::group_by(value) |>
    dplyr::mutate(
      q_exposure = stats::p.adjust(pval, method = "BH")
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      analysis = analysis_name,
      model_type = definition$model_type,
      population = definition$population,
      n_complete_cases = nrow(model_metadata),
      conf_low = coef - 1.96 * stderr,
      conf_high = coef + 1.96 * stderr,
      interpretation =
        "Summed HUMAnN species-stratified PWY-5265 contributor abundance within genus"
    )
  
  readr::write_csv(
    res,
    file.path(analysis_dir, "genus_contributor_exposure_results.csv")
  )
  
  genus_model_outputs[[analysis_name]] <- res
}

genus_adjusted_results <- dplyr::bind_rows(
  genus_model_outputs[purrr::map_lgl(genus_model_outputs, ~ is.data.frame(.x) && "coef" %in% names(.x))]
)

if (nrow(genus_adjusted_results) > 0) {
  readr::write_csv(
    genus_adjusted_results,
    file.path(
      OUTPUT_DIR,
      "PWY5265_GENUS_AGGREGATED_PE_VS_CONTROL_ALL_MODELS.csv"
    )
  )
}

genus_run_status <- purrr::imap_dfr(
  genus_model_outputs,
  function(x, nm) {
    if (is.data.frame(x) && "status" %in% names(x)) {
      tibble::tibble(analysis = nm, status = x$status[[1]])
    } else {
      tibble::tibble(analysis = nm, status = "completed")
    }
  }
)

readr::write_csv(
  genus_run_status,
  file.path(
    OUTPUT_DIR,
    "PWY5265_GENUS_AGGREGATED_RUN_STATUS.csv"
  )
)

message("")
message("Additional sparse-contributor analyses completed.")
message("Inspect:")
message("  PWY5265_species_presence_Fisher_PE_vs_Control.csv")
message("  PWY5265_genus_aggregated_contributor_summary.csv")
message("  PWY5265_genus_presence_Fisher_PE_vs_Control.csv")
message("  PWY5265_GENUS_AGGREGATED_RUN_STATUS.csv")
message("  PWY5265_GENUS_AGGREGATED_PE_VS_CONTROL_ALL_MODELS.csv, if created")
