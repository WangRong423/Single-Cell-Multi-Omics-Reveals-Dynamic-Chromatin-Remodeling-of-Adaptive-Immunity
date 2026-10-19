# Locate the companion file for source() and Rscript.
.local_source_file <- tryCatch(sys.frame(1)$ofile, error=function(e) NULL)
.local_args <- grep("^--file=", commandArgs(FALSE), value=TRUE)
.local_script_dir <- if (!is.null(.local_source_file)) {
  dirname(normalizePath(.local_source_file))
} else if (length(.local_args)) {
  dirname(normalizePath(sub("^--file=", "", .local_args[1])))
} else getwd()
source(file.path(.local_script_dir, "Tcell_updated_labels_helpers.R"))

## ========================================================================
## TCellSI Prime/Booster baseline-referenced z-scores
## with donor-paired statistical testing
##
## Statistical inference is performed on RAW TCellSI scores:
##   Score ~ 0 + timepoint + donor
##
## Prime contrasts:   day1/day3/day6/day14/day30 versus day0
## Booster contrasts: day31/day33/day36/day44 versus day30
## ========================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(limma)
})

## ========================================================================
## 0. Configuration
## ========================================================================

INPUT_DIR <- file.path(UPDATED_BASE, "TCellSI_TB_Subtype_reanalysis_updated_labels")
INPUT_CSV <- file.path(
  INPUT_DIR,
  "03_TCellSI",
  "TCellSI_scores_long_TB_Subtype.csv"
)
OUTDIR <- file.path(INPUT_DIR, "06_Phase_Zscore_Statistics")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

ZSCORE_METHOD <- "paired_change"
MIN_PAIRED_DONORS <- 3
DISPLAY_FDR_THRESHOLD <- 0.10
DISPLAY_FDR_COLUMN <- "FDR_phase"

PRIME_TIMES <- c("day0", "day1", "day3", "day6", "day14", "day30")
BOOSTER_TIMES <- c("day30", "day31", "day33", "day36", "day44")

PRIME_CONTRASTS <- list(
  Prime_D1_vs_D0 = c("day1", "day0"),
  Prime_D3_vs_D0 = c("day3", "day0"),
  Prime_D6_vs_D0 = c("day6", "day0"),
  Prime_D14_vs_D0 = c("day14", "day0"),
  Prime_D30_vs_D0 = c("day30", "day0")
)

BOOSTER_CONTRASTS <- list(
  Booster_D31_vs_D30 = c("day31", "day30"),
  Booster_D33_vs_D30 = c("day33", "day30"),
  Booster_D36_vs_D30 = c("day36", "day30"),
  Booster_D44_vs_D30 = c("day44", "day30")
)

STATE_ORDER <- c(
  "Quiescence", "Regulating", "Proliferation", "Helper",
  "Cytotoxicity", "Progenitor exhaustion",
  "Terminal exhaustion", "Senescence"
)

T_SUBTYPE_ORDER <- c(
  "CD4_NaiveT", "CD4_Naive_CCR4", "CD8_Memory_GZMK", "Treg", "CD4_Memory_ZEB2", "CD4_Memory",
  "CD8_NaiveT", "CD8_Memory", "gdT", "CD4_CD8_DP"
)

SUBTYPE_COLORS <- c(
  "CD4_NaiveT" = "#009E73",
  "CD4_Naive_CCR4" = "#CC79A7",
  "CD8_Memory_GZMK" = "#D55E00",
  "Treg" = "#F781BF",
  "CD4_Memory_ZEB2" = "#E41A1C",
  "CD4_Memory" = "#4DAF4A",
  "CD8_NaiveT" = "#D8A500",
  "CD8_Memory" = "#0072B2",
  "gdT" = "#A6761D",
  "CD4_CD8_DP" = "#6A3D9A"
)

save_pdf <- function(plot_object, filename, width = 11.5, height = 7.2) {
  ggplot2::ggsave(
    filename = filename,
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    device = grDevices::cairo_pdf
  )
}

safe_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(NA_real_)
  value <- stats::sd(x)
  if (!is.finite(value) || value <= 0) NA_real_ else value
}

significance_label <- function(fdr) {
  dplyr::case_when(
    is.na(fdr) ~ "",
    fdr < 0.001 ~ "***",
    fdr < 0.01 ~ "**",
    fdr < 0.05 ~ "*",
    fdr < 0.10 ~ "\u2020",
    TRUE ~ ""
  )
}

## ========================================================================
## 1. Read TCellSI long-format data
## ========================================================================

# Read the updated scoring output explicitly; do not reuse a stale session object.
if (!file.exists(INPUT_CSV)) stop("Run the updated TCellSI reanalysis script first: ", INPUT_CSV)
score_data <- read.csv(INPUT_CSV, stringsAsFactors=FALSE, check.names=FALSE)

required_columns <- c("donor", "timepoint", "subtype", "State", "Score")
missing_columns <- setdiff(required_columns, colnames(score_data))
if (length(missing_columns) > 0) {
  stop("The input is missing: ", paste(missing_columns, collapse = ", "))
}

score_data$subtype <- rename_tb_labels(score_data$subtype,
  file.path(OUTDIR, "Subtype_label_mapping.csv"))

score_data <- score_data %>%
  transmute(
    donor = as.character(donor),
    timepoint = as.character(timepoint),
    subtype = as.character(subtype),
    State = as.character(State),
    Score = as.numeric(Score)
  ) %>%
  filter(
    !is.na(donor), !is.na(timepoint), !is.na(subtype), !is.na(State),
    is.finite(Score), subtype %in% T_SUBTYPE_ORDER, State %in% STATE_ORDER
  ) %>%
  group_by(donor, timepoint, subtype, State) %>%
  summarise(Score = mean(Score, na.rm = TRUE), .groups = "drop")

write.csv(
  score_data,
  file.path(OUTDIR, "TCellSI_raw_scores_used_for_statistics.csv"),
  row.names = FALSE
)

## ========================================================================
## 2. Calculate phase-specific z-scores
## ========================================================================

phase_data <- bind_rows(
  score_data %>%
    filter(timepoint %in% PRIME_TIMES) %>%
    mutate(
      Phase = "Prime",
      baseline_time = "day0",
      phase_time_index = match(timepoint, PRIME_TIMES)
    ),
  score_data %>%
    filter(timepoint %in% BOOSTER_TIMES) %>%
    mutate(
      Phase = "Booster",
      baseline_time = "day30",
      phase_time_index = match(timepoint, BOOSTER_TIMES)
    )
)

baseline_by_donor <- phase_data %>%
  filter(timepoint == baseline_time) %>%
  transmute(
    donor, subtype, State, Phase, baseline_time,
    baseline_score = Score
  ) %>%
  distinct()

phase_data <- phase_data %>%
  left_join(
    baseline_by_donor,
    by = c("donor", "subtype", "State", "Phase", "baseline_time")
  ) %>%
  mutate(
    has_paired_baseline = is.finite(baseline_score),
    Delta_to_baseline = Score - baseline_score
  )

paired_scale <- phase_data %>%
  filter(
    has_paired_baseline,
    timepoint != baseline_time,
    is.finite(Delta_to_baseline)
  ) %>%
  group_by(subtype, State, Phase) %>%
  summarise(
    paired_change_sd = safe_sd(Delta_to_baseline),
    .groups = "drop"
  )

baseline_stats <- phase_data %>%
  filter(timepoint == baseline_time, is.finite(Score)) %>%
  group_by(subtype, State, Phase) %>%
  summarise(
    baseline_mean = mean(Score, na.rm = TRUE),
    baseline_sd = safe_sd(Score),
    .groups = "drop"
  )

phase_zscore <- phase_data %>%
  left_join(paired_scale, by = c("subtype", "State", "Phase")) %>%
  left_join(baseline_stats, by = c("subtype", "State", "Phase")) %>%
  mutate(
    Z_paired_change = ifelse(
      has_paired_baseline & is.finite(paired_change_sd) & paired_change_sd > 0,
      Delta_to_baseline / paired_change_sd,
      NA_real_
    ),
    Z_paired_change = ifelse(
      timepoint == baseline_time & has_paired_baseline,
      0,
      Z_paired_change
    ),
    Z_baseline_distribution = ifelse(
      is.finite(baseline_sd) & baseline_sd > 0,
      (Score - baseline_mean) / baseline_sd,
      NA_real_
    )
  )

if (ZSCORE_METHOD == "paired_change") {
  phase_zscore$Zscore <- phase_zscore$Z_paired_change
  y_axis_label <- paste0(
    "Paired baseline-referenced z-score\n",
    "(Prime vs day0; Booster vs day30)"
  )
} else if (ZSCORE_METHOD == "baseline_distribution") {
  phase_zscore$Zscore <- phase_zscore$Z_baseline_distribution
  y_axis_label <- paste0(
    "Baseline-distribution z-score\n",
    "(Prime vs day0; Booster vs day30)"
  )
} else {
  stop("ZSCORE_METHOD must be paired_change or baseline_distribution.")
}

write.csv(
  phase_zscore,
  file.path(OUTDIR, "TCellSI_phase_zscores_donor_level.csv"),
  row.names = FALSE
)

zscore_summary <- phase_zscore %>%
  filter(is.finite(Zscore)) %>%
  group_by(Phase, subtype, State, timepoint, phase_time_index) %>%
  summarise(
    mean_z = mean(Zscore, na.rm = TRUE),
    sd_z = stats::sd(Zscore, na.rm = TRUE),
    n_donor = n_distinct(donor),
    se_z = ifelse(n_donor > 1, sd_z / sqrt(n_donor), 0),
    .groups = "drop"
  )

## ========================================================================
## 3. Paired limma test for one subtype and contrast
## ========================================================================

run_paired_limma_test <- function(
    data,
    subtype_value,
    phase_name,
    contrast_name,
    target_time,
    baseline_time,
    min_paired_donors = 3
) {

  data_sub <- data %>%
    filter(
      subtype == subtype_value,
      timepoint %in% c(baseline_time, target_time)
    )

  common_donors <- intersect(
    unique(data_sub$donor[data_sub$timepoint == baseline_time]),
    unique(data_sub$donor[data_sub$timepoint == target_time])
  )

  if (length(common_donors) < min_paired_donors) {
    warning(
      "Skipping ", subtype_value, " / ", contrast_name,
      ": only ", length(common_donors), " paired donors."
    )
    return(NULL)
  }

  data_sub <- data_sub %>%
    filter(donor %in% common_donors) %>%
    mutate(sample_id = paste(donor, timepoint, sep = "__"))

  sample_meta <- data_sub %>%
    distinct(sample_id, donor, timepoint) %>%
    arrange(
      donor,
      factor(timepoint, levels = c(baseline_time, target_time))
    )

  score_matrix_df <- data_sub %>%
    select(State, sample_id, Score) %>%
    pivot_wider(names_from = sample_id, values_from = Score)

  missing_states <- setdiff(STATE_ORDER, score_matrix_df$State)
  if (length(missing_states) > 0) {
    warning(
      "Skipping ", subtype_value, " / ", contrast_name,
      ": missing states ", paste(missing_states, collapse = ", ")
    )
    return(NULL)
  }

  score_matrix_df <- score_matrix_df[
    match(STATE_ORDER, score_matrix_df$State),
    ,
    drop = FALSE
  ]

  score_matrix <- as.matrix(score_matrix_df[, -1, drop = FALSE])
  rownames(score_matrix) <- score_matrix_df$State
  score_matrix <- score_matrix[, sample_meta$sample_id, drop = FALSE]

  sample_meta$timepoint <- factor(
    sample_meta$timepoint,
    levels = c(baseline_time, target_time)
  )
  sample_meta$donor <- factor(sample_meta$donor)

  design <- model.matrix(~ 0 + timepoint + donor, data = sample_meta)
  colnames(design) <- sub("^timepoint", "", colnames(design))

  if (qr(design)$rank < ncol(design)) {
    stop("Rank-deficient design for ", subtype_value, " / ", contrast_name)
  }

  fit <- limma::lmFit(score_matrix, design)
  contrast_matrix <- limma::makeContrasts(
    contrasts = paste0(target_time, "-", baseline_time),
    levels = design
  )
  fit2 <- limma::contrasts.fit(fit, contrast_matrix)
  fit2 <- limma::eBayes(fit2, trend = TRUE, robust = TRUE)

  test_result <- limma::topTable(
    fit2,
    number = Inf,
    sort.by = "none"
  ) %>%
    rownames_to_column("State") %>%
    rename(
      mean_difference_limma = logFC,
      raw_P = P.Value
    ) %>%
    select(State, mean_difference_limma, AveExpr, t, raw_P, B)

  paired_effect <- data_sub %>%
    select(donor, State, timepoint, Score) %>%
    pivot_wider(names_from = timepoint, values_from = Score)

  paired_effect$paired_difference <-
    paired_effect[[target_time]] - paired_effect[[baseline_time]]

  paired_summary <- paired_effect %>%
    group_by(State) %>%
    summarise(
      mean_baseline = mean(.data[[baseline_time]], na.rm = TRUE),
      mean_target = mean(.data[[target_time]], na.rm = TRUE),
      mean_paired_difference = mean(paired_difference, na.rm = TRUE),
      median_paired_difference = median(paired_difference, na.rm = TRUE),
      sd_paired_difference = stats::sd(paired_difference, na.rm = TRUE),
      Cohen_dz = ifelse(
        is.finite(sd_paired_difference) & sd_paired_difference > 0,
        mean_paired_difference / sd_paired_difference,
        NA_real_
      ),
      n_paired_donors = sum(is.finite(paired_difference)),
      n_increased = sum(paired_difference > 0, na.rm = TRUE),
      n_decreased = sum(paired_difference < 0, na.rm = TRUE),
      n_unchanged = sum(paired_difference == 0, na.rm = TRUE),
      direction_consistency = pmax(n_increased, n_decreased) /
        pmax(n_paired_donors, 1),
      .groups = "drop"
    )

  test_result %>%
    left_join(paired_summary, by = "State") %>%
    mutate(
      Phase = phase_name,
      subtype = subtype_value,
      contrast = contrast_name,
      target_time = target_time,
      baseline_time = baseline_time
    )
}

## ========================================================================
## 4. Run all tests
## ========================================================================

all_test_results <- list()
result_counter <- 1

for (subtype_value in T_SUBTYPE_ORDER) {
  for (contrast_name in names(PRIME_CONTRASTS)) {
    values <- PRIME_CONTRASTS[[contrast_name]]
    all_test_results[[result_counter]] <- run_paired_limma_test(
      data = score_data,
      subtype_value = subtype_value,
      phase_name = "Prime",
      contrast_name = contrast_name,
      target_time = values[1],
      baseline_time = values[2],
      min_paired_donors = MIN_PAIRED_DONORS
    )
    result_counter <- result_counter + 1
  }

  for (contrast_name in names(BOOSTER_CONTRASTS)) {
    values <- BOOSTER_CONTRASTS[[contrast_name]]
    all_test_results[[result_counter]] <- run_paired_limma_test(
      data = score_data,
      subtype_value = subtype_value,
      phase_name = "Booster",
      contrast_name = contrast_name,
      target_time = values[1],
      baseline_time = values[2],
      min_paired_donors = MIN_PAIRED_DONORS
    )
    result_counter <- result_counter + 1
  }
}

statistics_all <- bind_rows(all_test_results)
if (nrow(statistics_all) == 0) {
  stop("No statistical comparisons were completed.")
}

statistics_all <- statistics_all %>%
  mutate(FDR_global = p.adjust(raw_P, method = "BH")) %>%
  group_by(Phase) %>%
  mutate(FDR_phase = p.adjust(raw_P, method = "BH")) %>%
  ungroup() %>%
  group_by(Phase, subtype) %>%
  mutate(FDR_subtype_phase = p.adjust(raw_P, method = "BH")) %>%
  ungroup() %>%
  mutate(
    significance_global = significance_label(FDR_global),
    significance_phase = significance_label(FDR_phase),
    significance_subtype_phase = significance_label(FDR_subtype_phase),
    effect_direction = case_when(
      mean_paired_difference > 0 ~ "Increased",
      mean_paired_difference < 0 ~ "Decreased",
      TRUE ~ "Unchanged"
    )
  )

write.csv(
  statistics_all,
  file.path(OUTDIR, "TCellSI_paired_limma_all_statistical_tests.csv"),
  row.names = FALSE
)

statistics_significant <- statistics_all %>%
  filter(.data[[DISPLAY_FDR_COLUMN]] < DISPLAY_FDR_THRESHOLD) %>%
  arrange(Phase, subtype, target_time, .data[[DISPLAY_FDR_COLUMN]])

write.csv(
  statistics_significant,
  file.path(
    OUTDIR,
    paste0(
      "TCellSI_significant_tests_",
      DISPLAY_FDR_COLUMN,
      "_lt_",
      DISPLAY_FDR_THRESHOLD,
      ".csv"
    )
  ),
  row.names = FALSE
)

## ========================================================================
## 5. Join statistics to z-score summaries
## ========================================================================

statistics_for_plot <- statistics_all %>%
  transmute(
    Phase,
    subtype,
    State,
    timepoint = target_time,
    raw_P,
    FDR_global,
    FDR_phase,
    FDR_subtype_phase,
    mean_paired_difference,
    Cohen_dz,
    direction_consistency,
    significance = significance_label(.data[[DISPLAY_FDR_COLUMN]]),
    display_fdr = .data[[DISPLAY_FDR_COLUMN]]
  )

zscore_plot_data <- zscore_summary %>%
  left_join(
    statistics_for_plot,
    by = c("Phase", "subtype", "State", "timepoint")
  )

write.csv(
  zscore_plot_data,
  file.path(OUTDIR, "TCellSI_zscore_summary_with_statistics.csv"),
  row.names = FALSE
)

## ========================================================================
## 6. Z-score line plots with significance rings
## ========================================================================

make_phase_zscore_plot <- function(phase_name, phase_times, phase_labels) {
  plot_data <- zscore_plot_data %>% filter(Phase == phase_name)
  significant_points <- plot_data %>%
    filter(
      is.finite(display_fdr),
      display_fdr < DISPLAY_FDR_THRESHOLD
    )

  ggplot(
    plot_data,
    aes(
      x = phase_time_index,
      y = mean_z,
      group = subtype,
      color = subtype
    )
  ) +
    geom_hline(
      yintercept = 0,
      linewidth = 0.35,
      linetype = 2,
      color = "grey45"
    ) +
    geom_line(linewidth = 0.78, na.rm = TRUE) +
    geom_point(size = 1.5, na.rm = TRUE) +
    geom_errorbar(
      aes(ymin = mean_z - se_z, ymax = mean_z + se_z),
      width = 0.10,
      linewidth = 0.24,
      alpha = 0.60,
      na.rm = TRUE
    ) +
    geom_point(
      data = significant_points,
      aes(x = phase_time_index, y = mean_z),
      inherit.aes = FALSE,
      shape = 21,
      fill = NA,
      color = "black",
      stroke = 0.55,
      size = 2.8
    ) +
    facet_wrap(~ State, ncol = 4, scales = "fixed") +
    scale_x_continuous(
      breaks = seq_along(phase_times),
      labels = phase_labels,
      limits = c(0.75, length(phase_times) + 0.25),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_color_manual(
      values = SUBTYPE_COLORS,
      breaks = T_SUBTYPE_ORDER,
      drop = FALSE
    ) +
    labs(
      x = NULL,
      y = y_axis_label,
      color = "T-cell subtype",
      title = paste0(phase_name, " response"),
      subtitle = paste0(
        "Black rings: ", DISPLAY_FDR_COLUMN, " < ",
        DISPLAY_FDR_THRESHOLD,
        "; paired limma test on raw TCellSI scores"
      )
    ) +
    theme_classic(base_size = 9) +
    theme(
      plot.title = element_text(face = "bold", size = 10),
      plot.subtitle = element_text(size = 7.5),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
      axis.text.y = element_text(size = 7),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 9),
      legend.position = "right",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8),
      panel.spacing = grid::unit(3, "mm")
    ) +
    guides(
      color = guide_legend(
        override.aes = list(linewidth = 1, size = 2)
      )
    )
}

p_prime_zscore_stats <- make_phase_zscore_plot(
  "Prime",
  PRIME_TIMES,
  c("D0", "D1", "D3", "D6", "D14", "D30")
)

p_booster_zscore_stats <- make_phase_zscore_plot(
  "Booster",
  BOOSTER_TIMES,
  c("D30", "D31", "D33", "D36", "D44")
)

save_pdf(
  p_prime_zscore_stats,
  file.path(
    OUTDIR,
    paste0(
      "Fig_TCellSI_Prime_vs_day0_",
      ZSCORE_METHOD,
      "_with_statistics.pdf"
    )
  )
)

save_pdf(
  p_booster_zscore_stats,
  file.path(
    OUTDIR,
    paste0(
      "Fig_TCellSI_Booster_vs_day30_",
      ZSCORE_METHOD,
      "_with_statistics.pdf"
    )
  )
)

## ========================================================================
## 7. Statistical dot plots
## ========================================================================

make_statistics_dotplot <- function(phase_name, target_order) {
  plot_data <- statistics_all %>%
    filter(Phase == phase_name) %>%
    mutate(
      target_time = factor(target_time, levels = target_order),
      subtype = factor(subtype, levels = rev(T_SUBTYPE_ORDER)),
      display_fdr = .data[[DISPLAY_FDR_COLUMN]],
      neg_log10_fdr = -log10(pmax(display_fdr, 1e-300)),
      significance = significance_label(display_fdr)
    )

  ggplot(plot_data, aes(x = target_time, y = subtype)) +
    geom_point(
      aes(
        size = neg_log10_fdr,
        fill = mean_paired_difference
      ),
      shape = 21,
      stroke = 0.30
    ) +
    geom_text(
      aes(label = significance),
      size = 2.2,
      vjust = -1.0,
      color = "black"
    ) +
    facet_wrap(~ State, ncol = 4) +
    scale_fill_gradient2(
      low = "#3B6FB6",
      mid = "white",
      high = "#B64040",
      midpoint = 0,
      name = "Paired mean\nscore change"
    ) +
    scale_size_continuous(
      range = c(0.8, 5.5),
      name = expression(-log[10](FDR))
    ) +
    labs(
      x = NULL,
      y = NULL,
      title = paste0(phase_name, " paired statistical tests"),
      subtitle = paste0(
        "limma: score ~ timepoint + donor; FDR = ",
        DISPLAY_FDR_COLUMN,
        "; \u2020 <0.10, * <0.05, ** <0.01, *** <0.001"
      )
    ) +
    theme_classic(base_size = 8) +
    theme(
      plot.title = element_text(face = "bold", size = 10),
      plot.subtitle = element_text(size = 7.5),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
      axis.text.y = element_text(size = 7),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 8),
      panel.spacing = grid::unit(3, "mm"),
      legend.position = "right"
    )
}

p_prime_statistics <- make_statistics_dotplot(
  "Prime",
  c("day1", "day3", "day6", "day14", "day30")
)

p_booster_statistics <- make_statistics_dotplot(
  "Booster",
  c("day31", "day33", "day36", "day44")
)

save_pdf(
  p_prime_statistics,
  file.path(OUTDIR, "Fig_TCellSI_Prime_statistical_dotplot.pdf"),
  width = 11.5,
  height = 9
)

save_pdf(
  p_booster_statistics,
  file.path(OUTDIR, "Fig_TCellSI_Booster_statistical_dotplot.pdf"),
  width = 11.5,
  height = 9
)

## ========================================================================
## 8. Save all objects and method summary
## ========================================================================

saveRDS(
  list(
    zscore_method = ZSCORE_METHOD,
    display_fdr_column = DISPLAY_FDR_COLUMN,
    display_fdr_threshold = DISPLAY_FDR_THRESHOLD,
    donor_level_zscores = phase_zscore,
    zscore_summary = zscore_summary,
    statistics = statistics_all,
    significant_statistics = statistics_significant,
    plots = list(
      Prime_zscore = p_prime_zscore_stats,
      Booster_zscore = p_booster_zscore_stats,
      Prime_statistics = p_prime_statistics,
      Booster_statistics = p_booster_statistics
    )
  ),
  file.path(OUTDIR, "TCellSI_phase_zscore_statistics_objects.rds")
)

writeLines(
  c(
    paste0("ZSCORE_METHOD = ", ZSCORE_METHOD),
    paste0("MIN_PAIRED_DONORS = ", MIN_PAIRED_DONORS),
    paste0("DISPLAY_FDR_COLUMN = ", DISPLAY_FDR_COLUMN),
    paste0("DISPLAY_FDR_THRESHOLD = ", DISPLAY_FDR_THRESHOLD),
    "Statistical test: limma moderated paired model on raw TCellSI scores",
    "Model: Score ~ 0 + timepoint + donor",
    "Prime baseline: day0",
    "Booster baseline: day30",
    "Effect size: Cohen's dz from paired donor differences",
    "Multiple testing: global, phase-wise, and subtype-within-phase BH FDR"
  ),
  con = file.path(OUTDIR, "Statistical_method_summary.txt")
)

message("TCellSI z-score and statistical analysis completed.")
message("Output directory: ", normalizePath(OUTDIR, mustWork = FALSE))

writeLines(capture.output(sessionInfo()), file.path(OUTDIR, "sessionInfo_updated.txt"))
