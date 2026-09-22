#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(scales)
  library(ggrepel)
})

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
root <- Sys.getenv(
  "PREECLAMPSIA_PROJECT_DIR",
  unset = default_project_dir
)
pathway_file <- file.path(root, "WGS/HumaNN3/cleaned_100_samples/PREECLAMPSIA_pathabundance_unstratified_cpm_clean100.tsv.gz")
contrib_file <- file.path(root, "WGS/HumaNN3/PWY5265_species_contributors/PWY5265_summed_species_stratified_contributors_by_genus.csv")
out_dir <- file.path(root, "WGS/results/reviewer_adjusted_models/2026-09-04/functional_panel")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pal <- c(Control = "#619CFF", PE = "#F8766D")

# ---- PWY-5265 total pathway abundance ----
pa <- read.delim(gzfile(pathway_file), check.names = FALSE, stringsAsFactors = FALSE)
pathway_col <- names(pa)[1]
pwy_row <- pa[grepl("PWY-5265", pa[[pathway_col]], fixed = TRUE), , drop = FALSE]
if (nrow(pwy_row) != 1) stop("Could not identify exactly one PWY-5265 row")

pwy_long <- pwy_row |>
  pivot_longer(-all_of(pathway_col), names_to = "sample_col", values_to = "pwy5265_cpm") |>
  mutate(sample_id = sub("_Abundance$", "", sample_col),
         pwy5265_cpm = as.numeric(pwy5265_cpm))

meta <- read.csv(contrib_file, stringsAsFactors = FALSE) |>
  select(sample_id, study_group, ethnicity) |>
  distinct()

pwy_long <- pwy_long |>
  left_join(meta, by = "sample_id") |>
  filter(!is.na(study_group)) |>
  mutate(study_group = factor(study_group, levels = c("Control", "PE")),
         detectable = pwy5265_cpm > 1)

pwy_presence <- pwy_long |>
  group_by(study_group) |>
  summarise(detected = sum(detectable), total = n(), .groups = "drop") |>
  mutate(proportion = detected / total, label = paste0(detected, "/", total))

pwy_fisher <- fisher.test(table(pwy_long$detectable, pwy_long$study_group))

pwy_positive <- pwy_long |> filter(detectable)

# ---- Streptococcus contribution and presence/absence ----
strep <- read.csv(contrib_file, stringsAsFactors = FALSE) |>
  filter(genus_label == "Streptococcus") |>
  select(sample_id, study_group, cpm) |>
  group_by(sample_id, study_group) |>
  summarise(strep_cpm = sum(cpm, na.rm = TRUE), .groups = "drop") |>
  mutate(study_group = factor(study_group, levels = c("Control", "PE")),
         # Match the prespecified contributor-detection rule used in Table S19.
         detectable = strep_cpm > 1)

presence <- strep |>
  group_by(study_group) |>
  summarise(detected = sum(detectable), total = n(), .groups = "drop") |>
  mutate(label = paste0(detected, "/", total),
         proportion = detected / total)

theme_pub <- theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 9),
        plot.subtitle = element_text(hjust = 0.5, size = 7.5, color = "black"),
        axis.title = element_text(size = 10),
        axis.text = element_text(color = "black"),
        legend.position = "none",
        plot.margin = margin(5.5, 8, 5.5, 5.5))

p1 <- ggplot(pwy_presence, aes(study_group, proportion, fill = study_group)) +
  geom_col(width = 0.58, color = "black") +
  geom_text(aes(label = label), vjust = -0.35, size = 3.5) +
  scale_fill_manual(values = pal) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1), expand = expansion(mult = c(0, 0.08))) +
  labs(title = "PWY-5265 detection", subtitle = paste0("Detection threshold: >1 CPM; Fisher p = ", format.pval(pwy_fisher$p.value, digits = 3)), x = NULL, y = "Samples detected") +
  theme_pub

p2 <- ggplot(pwy_positive, aes(study_group, log10(pwy5265_cpm + 1), fill = study_group)) +
  geom_boxplot(width = 0.58, outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.09, size = 1.4, alpha = 0.7) +
  scale_fill_manual(values = pal) +
  labs(title = "PWY-5265 abundance among detected samples", subtitle = "MaAsLin2 minimal: p = 0.0029; adjusted p = 0.957", x = NULL, y = "log10(CPM + 1)") +
  theme_pub

p_strep_abund <- ggplot(strep |> filter(detectable), aes(study_group, log10(strep_cpm + 1), fill = study_group)) +
  geom_boxplot(width = 0.58, outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.09, size = 1.4, alpha = 0.7) +
  scale_fill_manual(values = pal) +
  labs(title = "Streptococcus abundance among detected samples", subtitle = "Targeted minimal: p = 0.0058; adjusted p = 0.0058", x = NULL, y = "log10(CPM + 1)") +
  theme_pub

p3 <- ggplot(presence, aes(study_group, proportion, fill = study_group)) +
  geom_col(width = 0.58, color = "black") +
  geom_text(aes(label = label), vjust = -0.35, size = 3.5) +
  scale_fill_manual(values = pal) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 0.30), expand = expansion(mult = c(0, 0.08))) +
  labs(title = "Detectable Streptococcus", subtitle = "Fisher p = 0.0035; adjusted p = 0.0070", x = NULL, y = "Participants") +
  theme_pub

final_plot <- ((p1 | p2) / (p_strep_abund | p3)) + plot_annotation(tag_levels = "A")
ggsave(file.path(out_dir, "Figure4C_functional_PWY5265_Streptococcus.pdf"), final_plot, width = 8.8, height = 7.0, device = cairo_pdf, bg = "white")
write.csv(pwy_presence, file.path(out_dir, "Figure4C_PWY5265_presence_summary.csv"), row.names = FALSE)
write.csv(presence, file.path(out_dir, "Figure4C_Streptococcus_presence_summary.csv"), row.names = FALSE)

# ---- Exploratory forest plot of all nominally associated pathways ----
result_files <- c(Minimal = file.path(root, "WGS/HumaNN3/functional_differential_abundance_adjusted_models/pathways/overall_status_minimal/exposure_results_all.csv"))
nominal_ids <- c("METHGLYUT-PWY", "PWY-5265", "PWY-5367", "PWY-822")
pathway_labels <- c(
  "METHGLYUT-PWY" = "METHGLYUT-PWY\nMethylglyoxal degradation",
  "PWY-5265" = "PWY-5265\nPeptidoglycan biosynthesis II",
  "PWY-5367" = "PWY-5367\nPetroselinate biosynthesis",
  "PWY-822" = "PWY-822\nFructan biosynthesis"
)

forest <- lapply(names(result_files), function(model_name) {
  z <- read.csv(result_files[[model_name]], stringsAsFactors = FALSE)
  z <- z[z$feature_type == "pathways" & z$metadata == "study_group" & z$value == "PE", , drop = FALSE]
  z$pathway_id <- NA_character_
  z$pathway_id[grepl("METHGLYUT", z$feature, fixed = TRUE)] <- "METHGLYUT-PWY"
  z$pathway_id[grepl("PWY.5265", z$feature, fixed = TRUE)] <- "PWY-5265"
  z$pathway_id[grepl("PWY.5367", z$feature, fixed = TRUE)] <- "PWY-5367"
  z$pathway_id[grepl("PWY.822", z$feature, fixed = TRUE)] <- "PWY-822"
  z <- z[z$pathway_id %in% nominal_ids & z$pval < 0.05, , drop = FALSE]
  z$model <- model_name
  z
}) |> bind_rows()

forest <- forest |>
  mutate(pathway_id = factor(pathway_id, levels = rev(nominal_ids)),
         pathway = factor(unname(pathway_labels[as.character(pathway_id)]), levels = rev(unname(pathway_labels))),
         model = factor(model, levels = names(result_files)),
         lower = coef - 1.96 * stderr,
         upper = coef + 1.96 * stderr,
         p_label = paste0("p=", formatC(pval, format = "f", digits = 3), "\nq=", formatC(q_exposure, format = "f", digits = 3)))

pd <- position_dodge(width = 0.62)
p_forest <- ggplot(forest, aes(x = coef, y = pathway, color = model)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey55") +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0, position = pd, linewidth = 0.5) +
  geom_point(position = pd, size = 2.2) +
  scale_color_manual(values = c(Minimal = "#1F78B4")) +
  labs(title = "Nominally associated maternal gut pathways",
       subtitle = "Prespecified minimal model; no pathway survived FDR correction",
       x = "Effect estimate (PE vs Control, 95% CI)", y = NULL, color = "Model") +
  theme_pub + theme(legend.position = "none", plot.subtitle = element_text(size = 7))

ggsave(file.path(out_dir, "Figure4C_exploratory_pathway_forest.pdf"), p_forest, width = 7.4, height = 4.2, device = cairo_pdf, bg = "white")
write.csv(forest, file.path(out_dir, "Figure4C_exploratory_pathway_forest_data.csv"), row.names = FALSE)

# ---- Minimal-model volcano plot with nominal pathways labelled ----
vdat <- read.csv(result_files[[1]], stringsAsFactors = FALSE) |>
  filter(feature_type == "pathways", metadata == "study_group", value == "PE") |>
  mutate(pathway_id = case_when(
    grepl("METHGLYUT", feature, fixed = TRUE) ~ "METHGLYUT-PWY",
    grepl("PWY.5265", feature, fixed = TRUE) ~ "PWY-5265",
    grepl("PWY.5367", feature, fixed = TRUE) ~ "PWY-5367",
    grepl("PWY.822", feature, fixed = TRUE) ~ "PWY-822",
    TRUE ~ NA_character_),
    nominal = !is.na(pathway_id) & pval < 0.05,
    neglog10p = -log10(pval))

p_volcano <- ggplot(vdat, aes(coef, neglog10p)) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, color = "grey55") +
  geom_point(aes(color = nominal), alpha = 0.65, size = 1.7) +
  geom_text_repel(data = subset(vdat, nominal), aes(label = unname(pathway_labels[pathway_id])), size = 2.7,
                  max.overlaps = Inf, box.padding = 0.45, min.segment.length = 0,
                  seed = 1, color = "#1F78B4") +
  scale_color_manual(values = c(`FALSE` = "grey70", `TRUE` = "#1F78B4"), guide = "none") +
  labs(title = "Maternal gut functional profiling: minimal model",
       subtitle = "Nominal pathways labelled; no pathway survived FDR correction",
       x = "MaAsLin2 coefficient (PE vs Control)", y = expression(-log[10](p))) +
  theme_pub

ggsave(file.path(out_dir, "Figure4C_minimal_pathway_volcano.pdf"), p_volcano, width = 5.8, height = 4.8, device = cairo_pdf, bg = "white")
write.csv(vdat, file.path(out_dir, "Figure4C_minimal_pathway_volcano_data.csv"), row.names = FALSE)
message("Wrote functional Panel C outputs to: ", normalizePath(out_dir))
