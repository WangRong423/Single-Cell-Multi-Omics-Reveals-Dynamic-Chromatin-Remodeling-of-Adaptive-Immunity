#!/usr/bin/env Rscript

# ==============================================================================
# 01 — B-cell variable motif deviations and top-20 longitudinal dynamics
#
# Input
# -----
# proj_TB_deviations_complete.rds
#
# Analysis
# --------
# 1. Subset the completed ArchR project to B-cell subtypes.
# 2. Rank the full motif-deviation names with ArchR::getVarDeviations().
# 3. Preserve the complete ArchR motif names, such as SPIB_336 or FOSL1_142.
# 4. Select the top 20 ranked deviations.
# 5. Aggregate the existing MotifMatrix by B-cell subtype x sample.
# 6. Plot longitudinal mean chromVAR z-scores with donor SEM.
#
# This script does NOT:
# - rerun addDeviationsMatrix()
# - collapse multiple motifs into one TF symbol
# - select TFs from a manually curated list
# ==============================================================================


# ==============================================================================
# 0. Locate and source the helper
# ==============================================================================

script_argument <- grep(
  "^--file=",
  commandArgs(
    trailingOnly = FALSE
  ),
  value = TRUE
)

if (length(script_argument) > 0) {
  SCRIPT_DIR <- dirname(
    normalizePath(
      sub(
        "^--file=",
        "",
        script_argument[1]
      )
    )
  )
} else {
  source_file <- tryCatch(
    sys.frame(1)$ofile,
    error = function(e) NULL
  )

  SCRIPT_DIR <- if (
    !is.null(source_file) &&
      length(source_file) == 1
  ) {
    dirname(
      normalizePath(
        source_file
      )
    )
  } else {
    getwd()
  }
}

HELPER_CANDIDATES <- file.path(
  SCRIPT_DIR, "Bcell_deviations_helpers_updated.R"
)

helper_file <- HELPER_CANDIDATES[
  file.exists(
    HELPER_CANDIDATES
  )
]

if (length(helper_file) == 0) {
  stop(
    "No B-cell helper file was found.\nChecked:\n",
    paste(
      HELPER_CANDIDATES,
      collapse = "\n"
    )
  )
}

source(
  helper_file[1],
  local = .GlobalEnv
)

required_helper_formals <- c(
  "project_b",
  "pseudobulk_metadata",
  "use_matrix",
  "cache_file",
  "divide_n",
  "scale_to",
  "preferred_assays"
)

actual_helper_formals <- names(
  formals(
    get_bcell_pseudobulk
  )
)

if (
  !all(
    required_helper_formals %in%
      actual_helper_formals
  )
) {
  stop(
    "Incorrect get_bcell_pseudobulk() function loaded.\nExpected: ",
    paste(
      required_helper_formals,
      collapse = ", "
    ),
    "\nLoaded: ",
    paste(
      actual_helper_formals,
      collapse = ", "
    )
  )
}


# ==============================================================================
# 1. Parameters
# ==============================================================================

TOP_N_DEVIATIONS <- 20L

MOTIF_MATRIX_NAME <- "MotifMatrix"

MOTIF_PB_CACHE <- file.path(
  BCELL_OBJECT_DIR,
  "Bcell_MotifMatrix_fullname_pseudobulk.rds"
)

REUSE_MOTIF_PB_CACHE <- FALSE

VAR_DEVIATION_TABLE_RDS <- file.path(
  BCELL_OBJECT_DIR,
  "Bcell_getVarDeviations_complete_table.rds"
)

TOP20_OBJECT_RDS <- file.path(
  BCELL_OBJECT_DIR,
  "Bcell_getVarDeviations_top20_dynamics_objects.rds"
)

FIGURE_WIDTH_MM <- 240
FIGURE_HEIGHT_MM <- 210


# ==============================================================================
# 2. Utilities specific to this script
# ==============================================================================

normalize_deviation_name <- function(x) {
  x <- as.character(
    x
  )

  x <- sub(
    "^(z|deviations)[:|_]",
    "",
    x,
    ignore.case = TRUE
  )

  x
}


standardize_var_deviations <- function(
    var_deviation_object
) {
  var_df <- as.data.frame(
    var_deviation_object
  ) |>
    tibble::rownames_to_column(
      "ArchRRowName"
    ) |>
    tibble::as_tibble()

  required_columns <- c(
    "name",
    "combinedVars"
  )

  missing_columns <- setdiff(
    required_columns,
    colnames(
      var_df
    )
  )

  if (length(missing_columns) > 0) {
    stop(
      "getVarDeviations() result lacks: ",
      paste(
        missing_columns,
        collapse = ", "
      ),
      ". Available columns: ",
      paste(
        colnames(
          var_df
        ),
        collapse = ", "
      )
    )
  }

  if ("seqnames" %in% colnames(var_df)) {
    z_rows <- as.character(
      var_df$seqnames
    ) ==
      "z"

    if (any(z_rows)) {
      var_df <- var_df[
        z_rows,
        ,
        drop = FALSE
      ]
    }
  }

  if (!"combinedMeans" %in% colnames(var_df)) {
    var_df$combinedMeans <- NA_real_
  }

  if (!"rank" %in% colnames(var_df)) {
    var_df <- var_df |>
      dplyr::arrange(
        dplyr::desc(
          combinedVars
        )
      ) |>
      dplyr::mutate(
        rank = dplyr::row_number()
      )
  }

  if (!"seqnames" %in% colnames(var_df)) {
    var_df$seqnames <- rep(NA_character_, nrow(var_df))
  }

  var_df |>
    dplyr::transmute(
      ArchRRowName = ArchRRowName,
      Seqnames = as.character(seqnames),
      FullDeviationName = as.character(
        name
      ),
      NormalizedDeviationName =
        normalize_deviation_name(
          name
        ),
      CombinedVariance = as.numeric(
        combinedVars
      ),
      CombinedMean = as.numeric(
        combinedMeans
      ),
      Rank = as.integer(
        rank
      )
    ) |>
    dplyr::filter(
      !is.na(
        FullDeviationName
      ),
      FullDeviationName != "",
      is.finite(
        CombinedVariance
      ),
      is.finite(
        Rank
      )
    ) |>
    dplyr::arrange(
      Rank,
      dplyr::desc(
        CombinedVariance
      )
    ) |>
    dplyr::distinct(
      FullDeviationName,
      .keep_all = TRUE
    )
}


align_top_deviations_to_pseudobulk <- function(
    top_deviation_table,
    pseudobulk_feature_names
) {
  pseudobulk_feature_table <- tibble::tibble(
    PseudobulkRowIndex = seq_along(
      pseudobulk_feature_names
    ),
    PseudobulkFeatureName = as.character(
      pseudobulk_feature_names
    ),
    NormalizedDeviationName =
      normalize_deviation_name(
        pseudobulk_feature_names
      )
  )

  exact_match <- match(
    top_deviation_table$FullDeviationName,
    pseudobulk_feature_table$PseudobulkFeatureName
  )

  normalized_match <- match(
    top_deviation_table$NormalizedDeviationName,
    pseudobulk_feature_table$NormalizedDeviationName
  )

  selected_match <- exact_match

  selected_match[
    is.na(
      selected_match
    )
  ] <- normalized_match[
    is.na(
      selected_match
    )
  ]

  output <- top_deviation_table |>
    dplyr::mutate(
      PseudobulkRowIndex = selected_match,
      MatchType = dplyr::case_when(
        !is.na(
          exact_match
        ) ~
          "Exact",
        is.na(
          exact_match
        ) &
          !is.na(
            normalized_match
          ) ~
          "Normalized",
        TRUE ~
          "Unmatched"
      ),
      PseudobulkFeatureName =
        dplyr::if_else(
          !is.na(
            PseudobulkRowIndex
          ),
          pseudobulk_feature_table$PseudobulkFeatureName[
            PseudobulkRowIndex
          ],
          NA_character_
        )
    )

  duplicated_rows <- output$PseudobulkRowIndex[
    !is.na(
      output$PseudobulkRowIndex
    )
  ]

  if (anyDuplicated(duplicated_rows)) {
    warning(
      "More than one top deviation matched the same pseudobulk row. ",
      "Inspect the alignment table."
    )
  }

  output
}


# ==============================================================================
# 3. Load completed project and subset B cells
# ==============================================================================

project <- load_deviations_project()

bcell_object <- prepare_b_project(
  project
)

project_b <- bcell_object$Project
pseudobulk_metadata <- bcell_object$PBMetadata

message(
  "B cells retained: ",
  format(
    length(
      ArchR::getCellNames(
        project_b
      )
    ),
    big.mark = ","
  )
)

print(
  table(
    as.character(
      project_b$TB_Subtype
    ),
    useNA = "ifany"
  )
)


# ==============================================================================
# 4. Rank complete deviation names with ArchR::getVarDeviations()
# ==============================================================================

var_deviation_data <- ArchR::getVarDeviations(
  ArchRProj = project_b,
  name = MOTIF_MATRIX_NAME,
  plot = FALSE,
  n = TOP_N_DEVIATIONS
)

var_deviation_table <- standardize_var_deviations(
  var_deviation_data
)

if (
  nrow(
    var_deviation_table
  ) <
    TOP_N_DEVIATIONS
) {
  stop(
    "Only ",
    nrow(
      var_deviation_table
    ),
    " variable deviations were returned."
  )
}

saveRDS(
  var_deviation_table,
  VAR_DEVIATION_TABLE_RDS,
  compress = FALSE
)

write.csv(
  var_deviation_table,
  file.path(
    BCELL_TABLE_DIR,
    "Bcell_getVarDeviations_complete_ranking.csv"
  ),
  row.names = FALSE
)

top20_deviation_table <- var_deviation_table |>
  dplyr::slice_head(
    n = TOP_N_DEVIATIONS
  )

write.csv(
  top20_deviation_table,
  file.path(
    BCELL_TABLE_DIR,
    "Bcell_getVarDeviations_top20_full_names.csv"
  ),
  row.names = FALSE
)

message(
  "Top ",
  TOP_N_DEVIATIONS,
  " complete deviation names:\n",
  paste(
    sprintf(
      "%02d. %s",
      top20_deviation_table$Rank,
      top20_deviation_table$FullDeviationName
    ),
    collapse = "\n"
  )
)

p_var_deviations <- ArchR::getVarDeviations(
  ArchRProj = project_b,
  name = MOTIF_MATRIX_NAME,
  plot = TRUE,
  n = TOP_N_DEVIATIONS
)

ggplot2::ggsave(
  filename = file.path(
    BCELL_PLOT_DIR,
    "Fig4A1_Bcell_getVarDeviations_top20.pdf"
  ),
  plot = p_var_deviations,
  width = 120,
  height = 100,
  units = "mm",
  limitsize = FALSE
)

ggplot2::ggsave(
  filename = file.path(
    BCELL_PLOT_DIR,
    "Fig4A1_Bcell_getVarDeviations_top20.png"
  ),
  plot = p_var_deviations,
  width = 120,
  height = 100,
  units = "mm",
  dpi = 600,
  limitsize = FALSE
)


# ==============================================================================
# 5. Aggregate the existing MotifMatrix by subtype x sample
# ==============================================================================

if (
  REUSE_MOTIF_PB_CACHE &&
    file.exists(
      MOTIF_PB_CACHE
    )
) {
  motif_pb <- readRDS(
    MOTIF_PB_CACHE
  )

  message(
    "Loaded cached full-name MotifMatrix pseudobulk:\n",
    MOTIF_PB_CACHE
  )
} else {
  motif_pb <- get_bcell_pseudobulk(
    project_b = project_b,
    pseudobulk_metadata =
      pseudobulk_metadata,
    use_matrix = MOTIF_MATRIX_NAME,
    cache_file = MOTIF_PB_CACHE,
    divide_n = TRUE,
    scale_to = NULL,
    preferred_assays = c(
      "z",
      "deviations"
    )
  )
}

motif_pseudobulk_matrix <- as.matrix(
  motif_pb$Matrix
)

motif_pseudobulk_metadata <- motif_pb$Metadata

motif_pseudobulk_feature_names <-
  extract_feature_names(
    summarized_experiment =
      motif_pb$SE,
    matrix_object =
      motif_pb$Matrix
  )

rownames(
  motif_pseudobulk_matrix
) <- motif_pseudobulk_feature_names

motif_assay_name <- if (
  "z" %in%
    motif_pb$AssayNames
) {
  "z"
} else if (
  "deviations" %in%
    motif_pb$AssayNames
) {
  "deviations"
} else {
  motif_pb$AssayNames[1]
}

message(
  "MotifMatrix assay used for lines: ",
  motif_assay_name
)


# ==============================================================================
# 6. Match the top 20 complete names to pseudobulk matrix rows
# ==============================================================================

top20_alignment <- align_top_deviations_to_pseudobulk(
  top_deviation_table =
    top20_deviation_table,
  pseudobulk_feature_names =
    motif_pseudobulk_feature_names
)

write.csv(
  top20_alignment,
  file.path(
    BCELL_TABLE_DIR,
    "Bcell_getVarDeviations_top20_pseudobulk_alignment.csv"
  ),
  row.names = FALSE
)

unmatched_names <- top20_alignment |>
  dplyr::filter(
    is.na(
      PseudobulkRowIndex
    )
  ) |>
  dplyr::pull(
    FullDeviationName
  )

if (length(unmatched_names) > 0) {
  warning(
    "The following getVarDeviations names were not found in the ",
    "pseudobulk matrix:\n",
    paste(
      unmatched_names,
      collapse = "\n"
    )
  )
}

matched_top20 <- top20_alignment |>
  dplyr::filter(
    !is.na(
      PseudobulkRowIndex
    )
  )

if (
  nrow(
    matched_top20
  ) ==
    0
) {
  stop(
    "None of the top deviations matched the MotifMatrix pseudobulk rows."
  )
}

top20_matrix <- motif_pseudobulk_matrix[
  matched_top20$PseudobulkRowIndex,
  ,
  drop = FALSE
]

rownames(
  top20_matrix
) <- matched_top20$FullDeviationName


# ==============================================================================
# 7. Top-20 longitudinal line plot
# ==============================================================================

top20_sample_long <- as.data.frame(
  t(
    top20_matrix
  )
) |>
  tibble::rownames_to_column(
    "B_PB_Group"
  ) |>
  tidyr::pivot_longer(
    cols = -B_PB_Group,
    names_to = "FullDeviationName",
    values_to = "MotifDeviation"
  ) |>
  dplyr::left_join(
    motif_pseudobulk_metadata,
    by = "B_PB_Group"
  ) |>
  dplyr::left_join(
    matched_top20 |>
      dplyr::select(
        FullDeviationName,
        Rank,
        CombinedVariance,
        CombinedMean
      ),
    by = "FullDeviationName"
  ) |>
  dplyr::mutate(
    TB_Subtype = factor(
      as.character(
        TB_Subtype
      ),
      levels = B_SUBTYPES
    ),
    Timepoint = factor(
      as.character(
        Timepoint
      ),
      levels = TIME_LEVELS
    ),
    FullDeviationName = factor(
      FullDeviationName,
      levels = matched_top20$FullDeviationName
    )
  )

write.csv(
  top20_sample_long,
  gzfile(
    file.path(
      BCELL_TABLE_DIR,
      "Bcell_getVarDeviations_top20_sample_scores.csv.gz"
    )
  ),
  row.names = FALSE
)

top20_time_summary <- mean_sem(
  data = top20_sample_long,
  value = "MotifDeviation",
  groups = c(
    "FullDeviationName",
    "Rank",
    "TB_Subtype",
    "Timepoint"
  )
) |>
  dplyr::mutate(
    FullDeviationName = factor(
      as.character(
        FullDeviationName
      ),
      levels = matched_top20$FullDeviationName
    ),
    Timepoint = factor(
      as.character(
        Timepoint
      ),
      levels = TIME_LEVELS,
      labels = TIME_LABELS
    )
  )

write.csv(
  top20_time_summary,
  file.path(
    BCELL_TABLE_DIR,
    "Bcell_getVarDeviations_top20_time_mean_SEM.csv"
  ),
  row.names = FALSE
)

p_top20_lines <- ggplot2::ggplot(
  top20_time_summary,
  ggplot2::aes(
    x = Timepoint,
    y = Mean,
    color = TB_Subtype,
    fill = TB_Subtype,
    group = TB_Subtype
  )
) +
  ggplot2::geom_hline(
    yintercept = 0,
    linewidth = 0.25,
    linetype = 2,
    color = "grey72"
  ) +
  ggplot2::geom_ribbon(
    ggplot2::aes(
      ymin = Mean -
        SEM,
      ymax = Mean +
        SEM
    ),
    alpha = 0.10,
    color = NA
  ) +
  ggplot2::geom_line(
    linewidth = 0.72
  ) +
  ggplot2::geom_point(
    size = 1.35
  ) +
  ggplot2::facet_wrap(
    ~ FullDeviationName,
    ncol = 5,
    scales = "free_y"
  ) +
  ggplot2::scale_color_manual(
    values = B_SUBTYPE_COLORS,
    drop = FALSE
  ) +
  ggplot2::scale_fill_manual(
    values = B_SUBTYPE_COLORS,
    drop = FALSE
  ) +
  ggplot2::labs(
    title = paste0(
      "Top ",
      TOP_N_DEVIATIONS,
      " variable motif deviations in B-cell subtypes"
    ),
    subtitle = paste0(
      "Full ArchR deviation names retained; MotifMatrix assay = ",
      motif_assay_name,
      "; lines show donor mean ± SEM"
    ),
    x = NULL,
    y = if (motif_assay_name == "z") "Mean chromVAR deviation Z-score" else "Mean chromVAR deviation",
    color = NULL,
    fill = NULL
  ) +
  ggplot2::theme_classic(
    base_size = 8.5
  ) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 45,
      hjust = 1,
      size = 6.3
    ),
    axis.text.y = ggplot2::element_text(
      size = 6.5
    ),
    strip.background = ggplot2::element_rect(
      fill = "grey95",
      color = "grey45",
      linewidth = 0.30
    ),
    strip.text = ggplot2::element_text(
      face = "bold",
      size = 7.0
    ),
    panel.border = ggplot2::element_rect(
      fill = NA,
      color = "grey65",
      linewidth = 0.25
    ),
    panel.spacing = grid::unit(
      0.45,
      "lines"
    ),
    legend.position = "bottom",
    plot.title = ggplot2::element_text(
      face = "bold",
      hjust = 0.5,
      size = 11
    ),
    plot.subtitle = ggplot2::element_text(
      hjust = 0.5,
      size = 7.8
    )
  )

ggplot2::ggsave(
  filename = file.path(
    BCELL_PLOT_DIR,
    "Fig4A2_Bcell_getVarDeviations_top20_longitudinal_lines.pdf"
  ),
  plot = p_top20_lines,
  width = FIGURE_WIDTH_MM,
  height = FIGURE_HEIGHT_MM,
  units = "mm",
  limitsize = FALSE
)

ggplot2::ggsave(
  filename = file.path(
    BCELL_PLOT_DIR,
    "Fig4A2_Bcell_getVarDeviations_top20_longitudinal_lines.png"
  ),
  plot = p_top20_lines,
  width = FIGURE_WIDTH_MM,
  height = FIGURE_HEIGHT_MM,
  units = "mm",
  dpi = 600,
  limitsize = FALSE
)


# ==============================================================================
# 8. Combined figure
# ==============================================================================

p_combined <- p_var_deviations +
  p_top20_lines +
  patchwork::plot_layout(
    widths = c(
      0.75,
      2.25
    ),
    guides = "collect"
  ) +
  patchwork::plot_annotation(
    tag_levels = "A"
  ) &
  ggplot2::theme(
    legend.position = "bottom"
  )

ggplot2::ggsave(
  filename = file.path(
    BCELL_PLOT_DIR,
    "Fig4A_Bcell_getVarDeviations_top20_combined.pdf"
  ),
  plot = p_combined,
  width = 300,
  height = FIGURE_HEIGHT_MM,
  units = "mm",
  limitsize = FALSE
)

ggplot2::ggsave(
  filename = file.path(
    BCELL_PLOT_DIR,
    "Fig4A_Bcell_getVarDeviations_top20_combined.png"
  ),
  plot = p_combined,
  width = 300,
  height = FIGURE_HEIGHT_MM,
  units = "mm",
  dpi = 600,
  limitsize = FALSE
)


# ==============================================================================
# 9. Save reusable objects
# ==============================================================================

saveRDS(
  list(
    Parameters = list(
      TOP_N_DEVIATIONS =
        TOP_N_DEVIATIONS,
      MOTIF_MATRIX_NAME =
        MOTIF_MATRIX_NAME,
      MotifAssayName =
        motif_assay_name
    ),
    VarDeviationData =
      var_deviation_data,
    VarDeviationTable =
      var_deviation_table,
    Top20DeviationTable =
      top20_deviation_table,
    Top20Alignment =
      top20_alignment,
    MotifPseudobulkMatrix =
      motif_pseudobulk_matrix,
    MotifPseudobulkMetadata =
      motif_pseudobulk_metadata,
    Top20Matrix =
      top20_matrix,
    Top20SampleLong =
      top20_sample_long,
    Top20TimeSummary =
      top20_time_summary,
    pVarDeviations =
      p_var_deviations,
    pTop20Lines =
      p_top20_lines,
    pCombined =
      p_combined
  ),
  TOP20_OBJECT_RDS,
  compress = FALSE
)

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    BCELL_OUTDIR,
    "sessionInfo_getVarDeviations_top20.txt"
  )
)

cat(
  "\n============================================================\n",
  "B-cell getVarDeviations top-20 analysis completed.\n",
  "\nArchR variable-deviation plot:\n",
  file.path(
    BCELL_PLOT_DIR,
    "Fig4A1_Bcell_getVarDeviations_top20.pdf"
  ),
  "\n\nTop-20 longitudinal line plot:\n",
  file.path(
    BCELL_PLOT_DIR,
    "Fig4A2_Bcell_getVarDeviations_top20_longitudinal_lines.pdf"
  ),
  "\n\nCombined figure:\n",
  file.path(
    BCELL_PLOT_DIR,
    "Fig4A_Bcell_getVarDeviations_top20_combined.pdf"
  ),
  "\n\nTop-20 full names:\n",
  file.path(
    BCELL_TABLE_DIR,
    "Bcell_getVarDeviations_top20_full_names.csv"
  ),
  "\n============================================================\n",
  sep = ""
)

writeLines(capture.output(sessionInfo()), file.path(BCELL_OUTDIR, "sessionInfo_updated_labels.txt"))
