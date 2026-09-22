
# PREECLAMPSIA WGS FUNCTIONAL DIFFERENTIAL ABUNDANCE ANALYSI####
# PREECLAMPSIA WGS FUNCTIONAL DIFFERENTIAL ABUNDANCE ANALYSIS
#
# Inputs:
#   1) HUMAnN3 unstratified pathway abundance, CPM-normalized
#   2) HUMAnN3 unstratified gene-family abundance, CPM-normalized
#   3) Clinical metadata: met_clean0.csv
#
# Primary exposure:
#   study_group (PE vs Control; Control is the reference)
#
# Overall PE vs Control models:
#   - Prespecified minimal: study_group + ethnicity + bmi + antibiotics_during_pregnancy
#   - Fully adjusted: minimal + gestational_age + gest_diabetes + meds_antihypertensives
#   - Both models repeated after restriction to term deliveries (>=37 weeks)
#
# Severity analyses:
#   - Primary parsimonious: severe + ethnicity + bmi + antibiotics_during_pregnancy
#   - Gestational-age sensitivity: primary + gestational_age
#   - Both repeated among term deliveries
#
# Ethnicity-stratified analyses:
#   - Primary parsimonious: study_group/severe + bmi + antibiotics_during_pregnancy
#   - Gestational-age sensitivity: primary + gestational_age
#
# Important statistical decisions:
#   - No rarefaction is performed.
#   - CPM tables are already library-size normalized.
#   - MaAsLin2 is run with normalization = "NONE" and transform = "LOG".
#   - Confounders are retained; no stepwise covariate selection is performed.
#   - Rare features are filtered independently within each analytical population.
#   - FDR is recalculated for the exposure of interest within each contrast.
#   - Unstratified profiles are analyzed here. Species-stratified contributions
#     should be examined later only for prioritized functions.


options(stringsAsFactors = FALSE)
set.seed(20260806)


# 0. Packages


required_packages <- c(
  "tidyverse",
  "Maaslin2",
  "ggrepel"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
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
  library(ggrepel)
})


# 1. User-editable paths
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

PATHWAY_FILE <- file.path(
  CLEAN_DIR,
  "PREECLAMPSIA_pathabundance_unstratified_cpm_clean100.tsv.gz"
)

GENEFAMILY_FILE <- file.path(
  CLEAN_DIR,
  "PREECLAMPSIA_genefamilies_unstratified_cpm_clean100.tsv.gz"
)

OUTPUT_DIR <- file.path(
  HUMANN_DIR,
  "functional_differential_abundance_adjusted_models"
)

dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# 2. Analysis parameters


# A feature is considered detected when its abundance is at least 1 CPM.
DETECTION_CPM <- 1

# A feature must be detected in at least 20% of the analytical population.
MIN_PREVALENCE <- 0.20

# For small subgroups, require detection in at least three participants.
MIN_DETECTED_SAMPLES <- 3

# Remove features with extremely small average abundance after the prevalence
# criterion. A mean of 1 CPM corresponds to 0.0001% of a CPM-normalized table.
MIN_MEAN_CPM <- 1

# Require at least this many complete cases to run a model.
MIN_MODEL_N <- 10

# Require at least this many complete cases per exposure level.
MIN_LEVEL_N <- 3

# MaAsLin2 parallel processing.
available_cores <- parallel::detectCores(logical = TRUE)
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

# Visualization thresholds.
PRIMARY_FDR <- 0.05
EXPLORATORY_FDR <- 0.25
NOMINAL_P <- 0.05
MARGINAL_P <- 0.08
MAX_LABELS_VOLCANO <- 10

settings <- tibble(
  parameter = c(
    "DETECTION_CPM",
    "MIN_PREVALENCE",
    "MIN_DETECTED_SAMPLES",
    "MIN_MEAN_CPM",
    "MIN_MODEL_N",
    "MIN_LEVEL_N",
    "N_CORES",
    "PRIMARY_FDR",
    "EXPLORATORY_FDR",
    "NOMINAL_P",
    "MARGINAL_P"
  ),
  value = c(
    DETECTION_CPM,
    MIN_PREVALENCE,
    MIN_DETECTED_SAMPLES,
    MIN_MEAN_CPM,
    MIN_MODEL_N,
    MIN_LEVEL_N,
    N_CORES,
    PRIMARY_FDR,
    EXPLORATORY_FDR,
    NOMINAL_P,
    MARGINAL_P
  )
)

write_csv(
  settings,
  file.path(
    OUTPUT_DIR,
    "analysis_settings.csv"
  )
)


# 3. Utility functions


safe_name <- function(x) {
  x |>
    stringr::str_replace_all("[^A-Za-z0-9._-]+", "_") |>
    stringr::str_replace_all("_+", "_") |>
    stringr::str_remove("^_") |>
    stringr::str_remove("_$")
}

make_short_label <- function(x, width = 65) {
  x |>
    stringr::str_replace("^UniRef90_", "") |>
    stringr::str_replace("^UniRef50_", "") |>
    stringr::str_trunc(width = width)
}

is_technical_feature <- function(x) {
  stringr::str_detect(
    stringr::str_to_upper(x),
    "^(UNMAPPED|UNINTEGRATED|UNGROUPED)(\\b|\\||:|$)"
  )
}

check_file <- function(path) {
  if (!file.exists(path)) {
    stop("Required file not found: ", path)
  }
}


# 4. Read and prepare clinical metadata


check_file(METADATA_FILE)

metadata <- readr::read_csv(
  METADATA_FILE,
  show_col_types = FALSE,
  name_repair = "minimal"
)

# Repair an unnamed first column without changing the biological variables.
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
    "The metadata is missing required variables: ",
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
    gestational_age = as.numeric(gestational_age),
    race = factor(race),
    ethnicity_analysis = case_when(
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
    bmi = as.numeric(bmi)
  )

if (anyDuplicated(metadata$sample_id) > 0) {
  duplicated_ids <- metadata$sample_id[
    duplicated(metadata$sample_id)
  ]
  
  stop(
    "Duplicated metadata sample IDs: ",
    paste(
      unique(duplicated_ids),
      collapse = ", "
    )
  )
}

metadata_summary <- metadata |>
  dplyr::summarise(
    n_participants = n(),
    n_control = sum(
      study_group == "Control",
      na.rm = TRUE
    ),
    n_pe = sum(
      study_group == "PE",
      na.rm = TRUE
    ),
    n_pe_nsf = sum(
      severe == "No",
      na.rm = TRUE
    ),
    n_pe_sf = sum(
      severe == "Yes",
      na.rm = TRUE
    ),
    n_hispanic = sum(
      ethnicity_analysis == "Hispanic",
      na.rm = TRUE
    ),
    n_white = sum(
      ethnicity_analysis == "White",
      na.rm = TRUE
    ),
    missing_bmi = sum(is.na(bmi)),
    missing_antibiotics = sum(
      is.na(
        antibiotics_during_pregnancy
      )
    ),
    missing_ethnicity = sum(
      is.na(ethnicity)
    )
  )

write_csv(
  metadata_summary,
  file.path(
    OUTPUT_DIR,
    "metadata_summary.csv"
  )
)

metadata_counts <- metadata |>
  dplyr::count(
    ethnicity_analysis,
    study_group,
    severe,
    severe_3cat,
    name = "n",
    .drop = FALSE
  )

write_csv(
  metadata_counts,
  file.path(
    OUTPUT_DIR,
    "metadata_group_counts.csv"
  )
)


# 5. Read HUMAnN3 unstratified CPM tables


read_humann_cpm <- function(
    input_file,
    feature_type
) {
  
  check_file(input_file)
  
  message(
    "Reading ",
    feature_type,
    ": ",
    input_file
  )
  
  humann <- readr::read_tsv(
    input_file,
    show_col_types = FALSE,
    name_repair = "minimal",
    progress = TRUE
  )
  
  names(humann)[1] <- "feature"
  
  humann <- humann |>
    dplyr::mutate(
      feature = as.character(feature)
    )
  
  # The requested analysis uses unstratified profiles only.
  if (any(
    stringr::str_detect(
      humann$feature,
      fixed("|")
    )
  )) {
    stop(
      feature_type,
      " input contains stratified rows with '|'. ",
      "Use the unstratified HUMAnN3 table."
    )
  }
  
  sample_columns <- names(humann)[-1]
  
  extracted_ids <- stringr::str_extract(
    sample_columns,
    "M[0-9]+Stool"
  )
  
  if (any(is.na(extracted_ids))) {
    bad_columns <- sample_columns[
      is.na(extracted_ids)
    ]
    
    stop(
      "Could not identify sample IDs in ",
      feature_type,
      " columns: ",
      paste(
        bad_columns,
        collapse = ", "
      )
    )
  }
  
  if (anyDuplicated(extracted_ids) > 0) {
    stop(
      "Duplicated HUMAnN sample IDs in ",
      feature_type
    )
  }
  
  removed_technical <- humann |>
    dplyr::filter(
      is_technical_feature(feature)
    ) |>
    dplyr::select(feature)
  
  write_csv(
    removed_technical,
    file.path(
      OUTPUT_DIR,
      paste0(
        "removed_technical_features_",
        feature_type,
        ".csv"
      )
    )
  )
  
  humann <- humann |>
    dplyr::filter(
      !is_technical_feature(feature)
    )
  
  feature_names <- make.unique(
    humann$feature
  )
  
  feature_by_sample <- as.matrix(
    humann[
      ,
      -1,
      drop = FALSE
    ]
  )
  
  storage.mode(
    feature_by_sample
  ) <- "numeric"
  
  rownames(
    feature_by_sample
  ) <- feature_names
  
  colnames(
    feature_by_sample
  ) <- extracted_ids
  
  if (any(
    !is.finite(
      feature_by_sample
    )
  )) {
    stop(
      "Non-finite abundance values found in ",
      feature_type
    )
  }
  
  if (any(
    feature_by_sample < 0
  )) {
    stop(
      "Negative abundance values found in ",
      feature_type
    )
  }
  
  feature_by_sample
}

pathway_cpm <- read_humann_cpm(
  PATHWAY_FILE,
  "pathways"
)

genefamily_cpm <- read_humann_cpm(
  GENEFAMILY_FILE,
  "gene_families"
)

feature_tables <- list(
  pathways = pathway_cpm,
  gene_families = genefamily_cpm
)


# 6. Confirm exact sample matching


matching_results <- map_dfr(
  names(feature_tables),
  function(feature_type) {
    
    feature_table <- feature_tables[[
      feature_type
    ]]
    
    humann_ids <- colnames(
      feature_table
    )
    
    metadata_ids <- metadata$sample_id
    
    tibble(
      feature_type = feature_type,
      metric = c(
        "HUMAnN samples",
        "Metadata participants",
        "Duplicated HUMAnN IDs",
        "HUMAnN IDs without metadata",
        "Metadata IDs without HUMAnN"
      ),
      value = c(
        length(humann_ids),
        length(metadata_ids),
        anyDuplicated(humann_ids),
        length(
          setdiff(
            humann_ids,
            metadata_ids
          )
        ),
        length(
          setdiff(
            metadata_ids,
            humann_ids
          )
        )
      )
    )
  }
)

write_csv(
  matching_results,
  file.path(
    OUTPUT_DIR,
    "sample_matching_summary.csv"
  )
)

for (feature_type in names(feature_tables)) {
  
  humann_ids <- colnames(
    feature_tables[[feature_type]]
  )
  
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
      "Sample mismatch for ",
      feature_type,
      ". HUMAnN-only: ",
      paste(
        humann_only,
        collapse = ", "
      ),
      "; metadata-only: ",
      paste(
        metadata_only,
        collapse = ", "
      )
    )
  }
}


# 7. Define analytical populations and models


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
    description = "PE versus Control, all participants, prespecified minimal model",
    subset_function = function(x) rep(TRUE, nrow(x)),
    exposure = "study_group",
    covariates = minimal_covariates,
    references = c(
      "study_group,Control",
      "ethnicity,Not Hispanic/Latina",
      "antibiotics_during_pregnancy,No"
    ),
    exploratory = FALSE
  ),
  
  overall_status_full = list(
    description = "PE versus Control, all participants, fully adjusted model",
    subset_function = function(x) rep(TRUE, nrow(x)),
    exposure = "study_group",
    covariates = full_covariates,
    references = c(
      "study_group,Control",
      "ethnicity,Not Hispanic/Latina",
      "antibiotics_during_pregnancy,No",
      "gest_diabetes,No",
      "meds_antihypertensives,No"
    ),
    exploratory = FALSE
  ),
  
  overall_status_term_minimal = list(
    description = "PE versus Control, term deliveries only, prespecified minimal model",
    subset_function = function(x) {
      !is.na(x$gestational_age) & x$gestational_age >= 37
    },
    exposure = "study_group",
    covariates = minimal_covariates,
    references = c(
      "study_group,Control",
      "ethnicity,Not Hispanic/Latina",
      "antibiotics_during_pregnancy,No"
    ),
    exploratory = FALSE
  ),
  
  overall_status_term_full = list(
    description = "PE versus Control, term deliveries only, fully adjusted model",
    subset_function = function(x) {
      !is.na(x$gestational_age) & x$gestational_age >= 37
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
    exploratory = FALSE
  ),
  
  overall_severity_3cat_primary = list(
    description = paste(
      "Control versus PE without severe features",
      "versus PE with severe features, parsimonious primary model"
    ),
    subset_function = function(x) rep(TRUE, nrow(x)),
    exposure = "severe_3cat",
    covariates = severity_primary_covariates,
    references = c(
      "severe_3cat,Control",
      "ethnicity,Not Hispanic/Latina",
      "antibiotics_during_pregnancy,No"
    ),
    exploratory = FALSE
  ),
  
  overall_severity_3cat_ga = list(
    description = paste(
      "Control versus PE without severe features",
      "versus PE with severe features, gestational-age sensitivity model"
    ),
    subset_function = function(x) rep(TRUE, nrow(x)),
    exposure = "severe_3cat",
    covariates = severity_ga_covariates,
    references = c(
      "severe_3cat,Control",
      "ethnicity,Not Hispanic/Latina",
      "antibiotics_during_pregnancy,No"
    ),
    exploratory = FALSE
  ),
  
  pe_severity_primary = list(
    description = "Severe versus non-severe PE among PE participants, parsimonious primary model",
    subset_function = function(x) x$study_group == "PE",
    exposure = "severe",
    covariates = severity_primary_covariates,
    references = c(
      "severe,No",
      "ethnicity,Not Hispanic/Latina",
      "antibiotics_during_pregnancy,No"
    ),
    exploratory = FALSE
  ),
  
  pe_severity_ga = list(
    description = "Severe versus non-severe PE among PE participants, gestational-age sensitivity model",
    subset_function = function(x) x$study_group == "PE",
    exposure = "severe",
    covariates = severity_ga_covariates,
    references = c(
      "severe,No",
      "ethnicity,Not Hispanic/Latina",
      "antibiotics_during_pregnancy,No"
    ),
    exploratory = FALSE
  ),
  
  pe_severity_term_primary = list(
    description = "Severe versus non-severe PE among term PE participants, parsimonious primary model",
    subset_function = function(x) {
      x$study_group == "PE" & !is.na(x$gestational_age) & x$gestational_age >= 37
    },
    exposure = "severe",
    covariates = severity_primary_covariates,
    references = c(
      "severe,No",
      "ethnicity,Not Hispanic/Latina",
      "antibiotics_during_pregnancy,No"
    ),
    exploratory = FALSE
  ),
  
  pe_severity_term_ga = list(
    description = "Severe versus non-severe PE among term PE participants, gestational-age sensitivity model",
    subset_function = function(x) {
      x$study_group == "PE" & !is.na(x$gestational_age) & x$gestational_age >= 37
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
    description = "PE versus Control among Hispanic participants, parsimonious primary model",
    subset_function = function(x) x$ethnicity_analysis == "Hispanic",
    exposure = "study_group",
    covariates = ethnicity_primary_covariates,
    references = c(
      "study_group,Control",
      "antibiotics_during_pregnancy,No"
    ),
    exploratory = TRUE
  ),
  
  hispanic_status_ga = list(
    description = "PE versus Control among Hispanic participants, gestational-age sensitivity model",
    subset_function = function(x) x$ethnicity_analysis == "Hispanic",
    exposure = "study_group",
    covariates = ethnicity_ga_covariates,
    references = c(
      "study_group,Control",
      "antibiotics_during_pregnancy,No"
    ),
    exploratory = TRUE
  ),
  
  white_status_primary = list(
    description = "PE versus Control among non-Hispanic White participants, parsimonious primary model",
    subset_function = function(x) x$ethnicity_analysis == "White",
    exposure = "study_group",
    covariates = ethnicity_primary_covariates,
    references = c(
      "study_group,Control",
      "antibiotics_during_pregnancy,No"
    ),
    exploratory = TRUE
  ),
  
  white_status_ga = list(
    description = "PE versus Control among non-Hispanic White participants, gestational-age sensitivity model",
    subset_function = function(x) x$ethnicity_analysis == "White",
    exposure = "study_group",
    covariates = ethnicity_ga_covariates,
    references = c(
      "study_group,Control",
      "antibiotics_during_pregnancy,No"
    ),
    exploratory = TRUE
  ),
  
  hispanic_pe_severity_primary = list(
    description = "Severe versus non-severe PE among Hispanic PE participants, parsimonious primary model",
    subset_function = function(x) {
      x$ethnicity_analysis == "Hispanic" & x$study_group == "PE"
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
    description = "Severe versus non-severe PE among Hispanic PE participants, gestational-age sensitivity model",
    subset_function = function(x) {
      x$ethnicity_analysis == "Hispanic" & x$study_group == "PE"
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
    description = "Severe versus non-severe PE among non-Hispanic White PE participants, parsimonious primary model",
    subset_function = function(x) {
      x$ethnicity_analysis == "White" & x$study_group == "PE"
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
    description = "Severe versus non-severe PE among non-Hispanic White PE participants, gestational-age sensitivity model",
    subset_function = function(x) {
      x$ethnicity_analysis == "White" & x$study_group == "PE"
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


# 8. Prepare metadata for one analysis


prepare_model_metadata <- function(
    full_metadata,
    analysis_name,
    definition
) {
  
  subset_index <- definition$subset_function(
    full_metadata
  )
  
  subset_index[
    is.na(subset_index)
  ] <- FALSE
  
  model_variables <- c(
    definition$exposure,
    definition$covariates
  )
  
  model_metadata <- full_metadata[
    subset_index,
    ,
    drop = FALSE
  ] |>
    dplyr::select(
      sample_id,
      all_of(
        model_variables
      )
    ) |>
    drop_na(
      all_of(
        model_variables
      )
    ) |>
    dplyr::mutate(
      across(
        where(is.factor),
        droplevels
      )
    )
  
  exposure_counts <- model_metadata |>
    dplyr::count(
      .data[[
        definition$exposure
      ]],
      name = "n"
    )
  
  names(exposure_counts)[1] <-
    "exposure_level"
  
  exposure_counts <- exposure_counts |>
    dplyr::mutate(
      analysis = analysis_name,
      exposure =
        definition$exposure,
      description =
        definition$description,
      exploratory =
        definition$exploratory,
      .before = 1
    )
  
  list(
    metadata = model_metadata,
    exposure_counts =
      exposure_counts
  )
}


# 9. Prefilter one functional table


prefilter_features <- function(
    feature_by_sample,
    sample_ids,
    feature_type,
    analysis_name
) {
  
  subset_table <- feature_by_sample[
    ,
    sample_ids,
    drop = FALSE
  ]
  
  n_samples <- ncol(
    subset_table
  )
  
  required_detected_n <- max(
    MIN_DETECTED_SAMPLES,
    ceiling(
      MIN_PREVALENCE *
        n_samples
    )
  )
  
  detected_n <- rowSums(
    subset_table >=
      DETECTION_CPM
  )
  
  prevalence <- detected_n /
    n_samples
  
  mean_cpm <- rowMeans(
    subset_table
  )
  
  variance_cpm <- apply(
    subset_table,
    1,
    stats::var
  )
  
  keep <- (
    detected_n >=
      required_detected_n &
      mean_cpm >=
      MIN_MEAN_CPM &
      is.finite(
        variance_cpm
      ) &
      variance_cpm > 0
  )
  
  filter_table <- tibble(
    feature = rownames(
      subset_table
    ),
    detected_n = detected_n,
    prevalence = prevalence,
    mean_cpm = mean_cpm,
    variance_cpm =
      variance_cpm,
    retained = keep
  ) |>
    dplyr::arrange(
      desc(retained),
      desc(prevalence),
      desc(mean_cpm)
    )
  
  analysis_dir <- file.path(
    OUTPUT_DIR,
    feature_type,
    analysis_name
  )
  
  dir.create(
    analysis_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  write_csv(
    filter_table,
    file.path(
      analysis_dir,
      "prefilter_feature_statistics.csv"
    )
  )
  
  filter_summary <- tibble(
    feature_type =
      feature_type,
    analysis =
      analysis_name,
    n_samples =
      n_samples,
    detection_cpm =
      DETECTION_CPM,
    minimum_prevalence =
      MIN_PREVALENCE,
    minimum_detected_n =
      required_detected_n,
    minimum_mean_cpm =
      MIN_MEAN_CPM,
    features_before_filter =
      nrow(subset_table),
    features_after_filter =
      sum(keep)
  )
  
  write_csv(
    filter_summary,
    file.path(
      analysis_dir,
      "prefilter_summary.csv"
    )
  )
  
  list(
    filtered_table =
      subset_table[
        keep,
        ,
        drop = FALSE
      ],
    filter_summary =
      filter_summary
  )
}


# 10. Volcano plot


create_volcano_plot <- function(
    exposure_results,
    feature_type,
    analysis_name,
    contrast_value,
    output_directory
) {
  
  contrast_results <- exposure_results |>
    dplyr::filter(
      value ==
        contrast_value
    ) |>
    dplyr::mutate(
      p_for_plot = pmax(
        pval,
        .Machine$double.xmin
      ),
      neg_log10_p =
        -log10(p_for_plot),
      evidence_category =
        case_when(
          q_exposure <
            PRIMARY_FDR ~
            "FDR q < 0.05",
          q_exposure <
            EXPLORATORY_FDR ~
            "Exploratory q < 0.25",
          pval <
            NOMINAL_P ~
            "Nominal p < 0.05",
          pval >=
            NOMINAL_P &
            pval <=
            MARGINAL_P ~
            "Marginal p = 0.05–0.08",
          TRUE ~
            "Not prioritized"
        ),
      short_label =
        make_short_label(feature)
    )
  
  labels <- contrast_results |>
    dplyr::filter(
      q_exposure <
        EXPLORATORY_FDR |
        pval <=
        MARGINAL_P
    ) |>
    dplyr::arrange(
      q_exposure,
      pval
    ) |>
    slice_head(
      n =
        MAX_LABELS_VOLCANO
    )
  
  if (nrow(labels) == 0) {
    labels <- contrast_results |>
      dplyr::arrange(pval) |>
      slice_head(
        n = 5
      )
  }
  
  plot_title <- paste0(
    feature_type,
    ": ",
    analysis_name,
    " — ",
    contrast_value
  )
  
  p <- ggplot(
    contrast_results,
    aes(
      x = coef,
      y = neg_log10_p,
      color =
        evidence_category
    )
  ) +
    geom_point(
      alpha = 0.70,
      size = 1.8
    ) +
    geom_vline(
      xintercept = 0,
      linetype = 2,
      linewidth = 0.4
    ) +
    geom_hline(
      yintercept =
        -log10(NOMINAL_P),
      linetype = 2,
      linewidth = 0.4
    ) +
    ggrepel::geom_text_repel(
      data = labels,
      aes(label = short_label),
      size = 3,
      max.overlaps = Inf,
      min.segment.length = 0,
      box.padding = 0.35,
      point.padding = 0.20,
      show.legend = FALSE
    ) +
    scale_color_manual(
      values = c(
        "FDR q < 0.05" =
          "#B2182B",
        "Exploratory q < 0.25" =
          "#EF8A62",
        "Nominal p < 0.05" =
          "#2166AC",
        "Marginal p = 0.05–0.08" =
          "#67A9CF",
        "Not prioritized" =
          "grey70"
      ),
      breaks = c(
        "FDR q < 0.05",
        "Exploratory q < 0.25",
        "Nominal p < 0.05",
        "Marginal p = 0.05–0.08",
        "Not prioritized"
      )
    ) +
    labs(
      title = plot_title,
      subtitle = paste(
        "Adjusted MaAsLin2 model;",
        "FDR recalculated for this exposure contrast"
      ),
      x = paste(
        "Adjusted MaAsLin2 coefficient",
        "(positive = higher in displayed level)"
      ),
      y = expression(
        -log[10](
          italic(p)
        )
      ),
      color = "Evidence"
    ) +
    theme_classic(
      base_size = 12
    ) +
    theme(
      legend.position = "right",
      plot.title =
        element_text(
          face = "bold"
        )
    )
  
  contrast_file_name <- safe_name(
    contrast_value
  )
  
  ggsave(
    filename = file.path(
      output_directory,
      paste0(
        "volcano_",
        contrast_file_name,
        ".pdf"
      )
    ),
    plot = p,
    width = 9,
    height = 7
  )
  
  ggsave(
    filename = file.path(
      output_directory,
      paste0(
        "volcano_",
        contrast_file_name,
        ".png"
      )
    ),
    plot = p,
    width = 9,
    height = 7,
    dpi = 300
  )
  
  invisible(p)
}


# 11. Run one MaAsLin2 analysis


run_functional_model <- function(
    feature_by_sample,
    feature_type,
    analysis_name,
    definition,
    full_metadata
) {
  
  message("")
  message(
    "============================================================"
  )
  message(
    "Feature type: ",
    feature_type
  )
  message(
    "Analysis: ",
    analysis_name
  )
  message(
    definition$description
  )
  message(
    "============================================================"
  )
  
  prepared <- prepare_model_metadata(
    full_metadata,
    analysis_name,
    definition
  )
  
  model_metadata <- prepared$metadata
  exposure_counts <-
    prepared$exposure_counts
  
  analysis_dir <- file.path(
    OUTPUT_DIR,
    feature_type,
    analysis_name
  )
  
  dir.create(
    analysis_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  write_csv(
    exposure_counts,
    file.path(
      analysis_dir,
      "model_exposure_counts.csv"
    )
  )
  
  model_variables <- c(
    definition$exposure,
    definition$covariates
  )
  
  model_manifest <- tibble(
    feature_type =
      feature_type,
    analysis =
      analysis_name,
    description =
      definition$description,
    exploratory =
      definition$exploratory,
    exposure =
      definition$exposure,
    covariates =
      paste(
        definition$covariates,
        collapse = " + "
      ),
    model_formula =
      paste(
        paste0(
          feature_type,
          "_abundance"
        ),
        "~",
        paste(
          model_variables,
          collapse = " + "
        )
      ),
    n_complete_cases =
      nrow(model_metadata)
  )
  
  write_csv(
    model_manifest,
    file.path(
      analysis_dir,
      "model_manifest.csv"
    )
  )
  
  if (
    nrow(model_metadata) <
    MIN_MODEL_N
  ) {
    warning(
      "Skipping ",
      feature_type,
      " / ",
      analysis_name,
      ": only ",
      nrow(model_metadata),
      " complete cases."
    )
    
    return(
      list(
        status = "skipped_small_n",
        exposure_results = NULL,
        manifest = model_manifest
      )
    )
  }
  
  exposure_vector <- model_metadata[[
    definition$exposure
  ]]
  
  if (
    length(
      unique(
        exposure_vector
      )
    ) < 2
  ) {
    warning(
      "Skipping ",
      feature_type,
      " / ",
      analysis_name,
      ": exposure has fewer than two levels."
    )
    
    return(
      list(
        status = "skipped_one_level",
        exposure_results = NULL,
        manifest = model_manifest
      )
    )
  }
  
  level_counts <- table(
    exposure_vector
  )
  
  if (
    any(
      level_counts <
      MIN_LEVEL_N
    )
  ) {
    warning(
      "Skipping ",
      feature_type,
      " / ",
      analysis_name,
      ": at least one exposure level has fewer than ",
      MIN_LEVEL_N,
      " samples."
    )
    
    return(
      list(
        status = "skipped_small_level",
        exposure_results = NULL,
        manifest = model_manifest
      )
    )
  }
  
  filtered <- prefilter_features(
    feature_by_sample =
      feature_by_sample,
    sample_ids =
      model_metadata$sample_id,
    feature_type =
      feature_type,
    analysis_name =
      analysis_name
  )
  
  filtered_table <-
    filtered$filtered_table
  
  if (
    nrow(filtered_table) == 0
  ) {
    warning(
      "Skipping ",
      feature_type,
      " / ",
      analysis_name,
      ": no features passed prefiltering."
    )
    
    return(
      list(
        status = "skipped_no_features",
        exposure_results = NULL,
        manifest = model_manifest
      )
    )
  }
  
  # MaAsLin2 expects samples as rows and features as columns.
  input_data <- as.data.frame(
    t(filtered_table),
    check.names = FALSE
  )
  
  input_metadata <- model_metadata |>
    dplyr::select(
      sample_id,
      all_of(
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
  
  message(
    "Running MaAsLin2 with ",
    nrow(input_data),
    " samples and ",
    ncol(input_data),
    " retained features."
  )
  
  Maaslin2::Maaslin2(
    input_data = input_data,
    input_metadata =
      input_metadata,
    output = maaslin_output,
    fixed_effects =
      model_variables,
    random_effects = NULL,
    reference =
      definition$references,
    normalization = "NONE",
    transform = "LOG",
    analysis_method = "LM",
    standardize = FALSE,
    min_abundance = 0,
    min_prevalence = 0,
    min_variance = 0,
    max_significance =
      EXPLORATORY_FDR,
    correction = "BH",
    cores = N_CORES
  )
  
  all_results_file <- file.path(
    maaslin_output,
    "all_results.tsv"
  )
  
  if (
    !file.exists(
      all_results_file
    )
  ) {
    stop(
      "MaAsLin2 did not generate all_results.tsv for ",
      feature_type,
      " / ",
      analysis_name
    )
  }
  
  all_results <- readr::read_tsv(
    all_results_file,
    show_col_types = FALSE
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
      "Unexpected MaAsLin2 result structure. Missing columns: ",
      paste(
        missing_result_columns,
        collapse = ", "
      )
    )
  }
  
  all_results <- all_results |>
    dplyr::mutate(
      feature_type =
        feature_type,
      analysis =
        analysis_name,
      analysis_description =
        definition$description,
      exploratory_analysis =
        definition$exploratory,
      n_complete_cases =
        nrow(model_metadata),
      .before = 1
    ) |>
    dplyr::arrange(
      qval,
      pval
    )
  
  write_csv(
    all_results,
    file.path(
      analysis_dir,
      "all_model_terms_sorted.csv"
    )
  )
  
  exposure_results <- all_results |>
    dplyr::filter(
      metadata ==
        definition$exposure
    ) |>
    dplyr::group_by(value) |>
    dplyr::mutate(
      # Primary q-value for interpretation:
      # BH correction across all tested functions for this exposure contrast.
      q_exposure =
        p.adjust(
          pval,
          method = "BH"
        )
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      direction =
        case_when(
          coef > 0 ~
            "Higher in displayed level",
          coef < 0 ~
            "Lower in displayed level",
          TRUE ~
            "No direction"
        ),
      evidence_category =
        case_when(
          q_exposure <
            PRIMARY_FDR ~
            "FDR q < 0.05",
          q_exposure <
            EXPLORATORY_FDR ~
            "Exploratory q < 0.25",
          pval <
            NOMINAL_P ~
            "Nominal p < 0.05",
          pval >=
            NOMINAL_P &
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
  
  if (
    nrow(exposure_results) == 0
  ) {
    warning(
      "No exposure coefficients were returned for ",
      feature_type,
      " / ",
      analysis_name
    )
    
    return(
      list(
        status = "completed_no_exposure_rows",
        exposure_results = NULL,
        manifest = model_manifest
      )
    )
  }
  
  write_csv(
    exposure_results,
    file.path(
      analysis_dir,
      "exposure_results_all.csv"
    )
  )
  
  prioritized_results <- exposure_results |>
    dplyr::filter(
      q_exposure <
        EXPLORATORY_FDR |
        pval <=
        MARGINAL_P
    )
  
  write_csv(
    prioritized_results,
    file.path(
      analysis_dir,
      "exposure_results_prioritized.csv"
    )
  )
  
  significance_summary <- exposure_results |>
    dplyr::count(
      value,
      evidence_category,
      name = "n_features"
    ) |>
    complete(
      value,
      evidence_category = c(
        "FDR q < 0.05",
        "Exploratory q < 0.25",
        "Nominal p < 0.05",
        "Marginal p = 0.05-0.08",
        "Not prioritized"
      ),
      fill = list(
        n_features = 0
      )
    ) |>
    dplyr::mutate(
      feature_type =
        feature_type,
      analysis =
        analysis_name,
      .before = 1
    )
  
  write_csv(
    significance_summary,
    file.path(
      analysis_dir,
      "exposure_significance_summary.csv"
    )
  )
  
  contrast_values <- unique(
    exposure_results$value
  )
  
  walk(
    contrast_values,
    ~ create_volcano_plot(
      exposure_results =
        exposure_results,
      feature_type =
        feature_type,
      analysis_name =
        analysis_name,
      contrast_value = .x,
      output_directory =
        analysis_dir
    )
  )
  
  list(
    status = "completed",
    exposure_results =
      exposure_results,
    manifest =
      model_manifest,
    filter_summary =
      filtered$filter_summary,
    significance_summary =
      significance_summary
  )
}


# 12. Run all feature types and all analyses


all_run_results <- list()
all_exposure_results <- list()
all_manifests <- list()
all_filter_summaries <- list()
all_significance_summaries <- list()

for (
  feature_type in
  names(feature_tables)
) {
  
  for (
    analysis_name in
    names(analysis_definitions)
  ) {
    
    definition <- analysis_definitions[[
      analysis_name
    ]]
    
    result <- run_functional_model(
      feature_by_sample =
        feature_tables[[
          feature_type
        ]],
      feature_type =
        feature_type,
      analysis_name =
        analysis_name,
      definition =
        definition,
      full_metadata =
        metadata
    )
    
    result_key <- paste(
      feature_type,
      analysis_name,
      sep = "__"
    )
    
    all_run_results[[
      result_key
    ]] <- tibble(
      feature_type =
        feature_type,
      analysis =
        analysis_name,
      status =
        result$status
    )
    
    if (
      !is.null(
        result$exposure_results
      )
    ) {
      all_exposure_results[[
        result_key
      ]] <-
        result$exposure_results
    }
    
    if (
      !is.null(
        result$manifest
      )
    ) {
      all_manifests[[
        result_key
      ]] <-
        result$manifest
    }
    
    if (
      !is.null(
        result$filter_summary
      )
    ) {
      all_filter_summaries[[
        result_key
      ]] <-
        result$filter_summary
    }
    
    if (
      !is.null(
        result$significance_summary
      )
    ) {
      all_significance_summaries[[
        result_key
      ]] <-
        result$significance_summary
    }
  }
}


# 13. Combined result tables


run_status <- bind_rows(
  all_run_results
)

write_csv(
  run_status,
  file.path(
    OUTPUT_DIR,
    "RUN_STATUS.csv"
  )
)

if (
  length(
    all_manifests
  ) > 0
) {
  write_csv(
    bind_rows(
      all_manifests
    ),
    file.path(
      OUTPUT_DIR,
      "ALL_MODEL_MANIFESTS.csv"
    )
  )
}

if (
  length(
    all_filter_summaries
  ) > 0
) {
  write_csv(
    bind_rows(
      all_filter_summaries
    ),
    file.path(
      OUTPUT_DIR,
      "ALL_PREFILTER_SUMMARIES.csv"
    )
  )
}

if (
  length(
    all_significance_summaries
  ) > 0
) {
  write_csv(
    bind_rows(
      all_significance_summaries
    ),
    file.path(
      OUTPUT_DIR,
      "ALL_SIGNIFICANCE_SUMMARIES.csv"
    )
  )
}

if (
  length(
    all_exposure_results
  ) > 0
) {
  
  combined_exposure_results <- bind_rows(
    all_exposure_results
  ) |>
    dplyr::arrange(
      feature_type,
      analysis,
      value,
      q_exposure,
      pval
    )
  
  write_csv(
    combined_exposure_results,
    file.path(
      OUTPUT_DIR,
      "ALL_EXPOSURE_RESULTS.csv"
    )
  )
  
  combined_prioritized <- combined_exposure_results |>
    dplyr::filter(
      q_exposure <
        EXPLORATORY_FDR |
        pval <=
        MARGINAL_P
    )
  
  write_csv(
    combined_prioritized,
    file.path(
      OUTPUT_DIR,
      "ALL_PRIORITIZED_RESULTS.csv"
    )
  )
  
  primary_overall_results <- combined_exposure_results |>
    dplyr::filter(
      analysis %in% c(
        "overall_status_minimal",
        "overall_status_full",
        "overall_status_term_minimal",
        "overall_status_term_full"
      ),
      value ==
        "PE"
    )
  
  write_csv(
    primary_overall_results,
    file.path(
      OUTPUT_DIR,
      "PRIMARY_PE_VS_CONTROL_ALL_MODELS.csv"
    )
  )
  
  write_csv(
    primary_overall_results |> dplyr::filter(analysis == "overall_status_minimal"),
    file.path(OUTPUT_DIR, "PRIMARY_PE_VS_CONTROL_MINIMAL_RESULTS.csv")
  )
  
  write_csv(
    primary_overall_results |> dplyr::filter(analysis == "overall_status_full"),
    file.path(OUTPUT_DIR, "PRIMARY_PE_VS_CONTROL_FULL_RESULTS.csv")
  )
  
  write_csv(
    primary_overall_results |> dplyr::filter(analysis == "overall_status_term_minimal"),
    file.path(OUTPUT_DIR, "PRIMARY_PE_VS_CONTROL_TERM_MINIMAL_RESULTS.csv")
  )
  
  write_csv(
    primary_overall_results |> dplyr::filter(analysis == "overall_status_term_full"),
    file.path(OUTPUT_DIR, "PRIMARY_PE_VS_CONTROL_TERM_FULL_RESULTS.csv")
  )
}


# 14. Reproducibility records


capture.output(
  sessionInfo(),
  file = file.path(
    OUTPUT_DIR,
    "sessionInfo.txt"
  )
)

write_lines(
  c(
    "Overall PE vs Control prespecified minimal model:",
    paste(
      "functional abundance ~ study_group + ethnicity + bmi +",
      "antibiotics_during_pregnancy"
    ),
    "",
    "Overall PE vs Control fully adjusted model:",
    paste(
      "functional abundance ~ study_group + ethnicity + bmi +",
      "antibiotics_during_pregnancy + gestational_age +",
      "gest_diabetes + meds_antihypertensives"
    ),
    "",
    "Both overall models were repeated after restriction to term deliveries (>=37 weeks).",
    "",
    "PE severity primary model:",
    paste(
      "functional abundance ~ severe + ethnicity + bmi +",
      "antibiotics_during_pregnancy"
    ),
    "",
    "PE severity gestational-age sensitivity model:",
    paste(
      "functional abundance ~ severe + ethnicity + bmi +",
      "antibiotics_during_pregnancy + gestational_age"
    ),
    "",
    "Severity models were also repeated among term PE participants.",
    "",
    "Ethnicity-stratified primary model:",
    paste(
      "functional abundance ~ study_group + bmi +",
      "antibiotics_during_pregnancy"
    ),
    "",
    "Ethnicity-stratified gestational-age sensitivity model:",
    paste(
      "functional abundance ~ study_group + bmi +",
      "antibiotics_during_pregnancy + gestational_age"
    ),
    "",
    "No stepwise selection was used. Prespecified covariates were retained.",
    "",
    "The input CPM files were not rarefied and were not normalized again.",
    "",
    "Technical HUMAnN rows were excluded before biological analyses.",
    "",
    paste0(
      "Prefilter: detection >= ",
      DETECTION_CPM,
      " CPM; prevalence >= ",
      MIN_PREVALENCE * 100,
      "%; mean abundance >= ",
      MIN_MEAN_CPM,
      " CPM; variance > 0."
    ),
    "",
    paste0(
      "FDR was recalculated by BH across tested functions ",
      "for each exposure contrast within each analytical model."
    )
  ),
  file.path(
    OUTPUT_DIR,
    "ANALYSIS_NOTES.txt"
  )
)

message("")
message(
  "All functional differential-abundance analyses finished."
)
message(
  "Results directory: ",
  OUTPUT_DIR
)
message(
  "Review RUN_STATUS.csv first."
)
message(
  "Primary PE vs Control comparison: PRIMARY_PE_VS_CONTROL_ALL_MODELS.csv"
)
message(
  "Primary results: PRIMARY_PE_VS_CONTROL_ADJUSTED_RESULTS.csv"
)
