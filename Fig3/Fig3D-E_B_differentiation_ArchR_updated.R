#!/usr/bin/env Rscript

# ==============================================================================
# Supervised B-cell state continuum — ArchR only
#
# Trajectory:
#   Naive_B -> Intermediate_B -> PlasmaB
#
# Input:
#   proj_TB_deviations_complete.rds
#
# Outputs:
#   1. ArchR trajectory embedding and pseudotime diagnostic plot
#   2. MotifMatrix trajectory heatmap
#   3. GeneExpressionMatrix trajectory heatmap
#   4. Cell-level pseudotime table
#   5. ArchR project containing the fitted trajectory
#
# This script:
#   - uses ArchR::addTrajectory()
#   - does not use Slingshot or Monocle
#   - completely excludes Memory_B from trajectory fitting
#   - does not rerun addDeviationsMatrix()
# ==============================================================================


# ==============================================================================
# 0. Locate shared helper
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
    error = function(e) {
      NULL
    }
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
    "B-cell helper file was not found.\nChecked:\n",
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

suppressPackageStartupMessages({
  library(ArchR)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(SummarizedExperiment)
  library(ComplexHeatmap)
})


# ==============================================================================
# 1. Parameters
# ==============================================================================

PLASMA_TRAJECTORY_GROUPS <- c(
  "Naive_B",
  "Intermediate_B",
  "PlasmaB"
)

TRAJECTORY_NAME <- "B_Plasma_StateContinuum"

GROUP_EVERY <- 2L
SMOOTH_WINDOW <- 5L

MOTIF_VAR_CUTOFF <- 0.90
MOTIF_MAX_FEATURES <- 5000L
MOTIF_LABEL_TOP <- 20L

RNA_VAR_CUTOFF <- 0.85
RNA_MAX_FEATURES <- 10000L
RNA_LABEL_TOP <- 20L

PLASMA_MARKER_GENES <- c(
  "CD79A",
  "MS4A1",
  "CD74",
  "CD37",
  "CD27",
  "IRF4",
  "PRDM1",
  "XBP1",
  "MZB1",
  "JCHAIN",
  "SDC1",
  "DERL3",
  "SEC11C",
  "SSR4",
  "IGKC"
)

PLASMA_TRAJECTORY_OUTDIR <- file.path(
  BCELL_OUTDIR,
  "PlasmaB_ArchR_trajectory"
)

PLASMA_TRAJECTORY_PLOT_DIR <- file.path(
  PLASMA_TRAJECTORY_OUTDIR,
  "Plots"
)

PLASMA_TRAJECTORY_TABLE_DIR <- file.path(
  PLASMA_TRAJECTORY_OUTDIR,
  "Tables"
)

PLASMA_TRAJECTORY_OBJECT_DIR <- file.path(
  PLASMA_TRAJECTORY_OUTDIR,
  "Objects"
)

for (
  output_directory in c(
    PLASMA_TRAJECTORY_PLOT_DIR,
    PLASMA_TRAJECTORY_TABLE_DIR,
    PLASMA_TRAJECTORY_OBJECT_DIR
  )
) {
  dir.create(
    output_directory,
    recursive = TRUE,
    showWarnings = FALSE
  )
}


# ==============================================================================
# 2. Load completed project and retain the predefined plasma-directed state axis
# ==============================================================================

project <- load_deviations_project()

bcell_object <- prepare_b_project(
  project
)

project_b <- bcell_object$Project

bcell_metadata <- bcell_object$CellMetadata

plasma_axis_metadata <- bcell_metadata |>
  dplyr::filter(
    TB_Subtype %in%
      PLASMA_TRAJECTORY_GROUPS
  ) |>
  dplyr::mutate(
    TB_Subtype = factor(
      TB_Subtype,
      levels = PLASMA_TRAJECTORY_GROUPS
    )
  )

if (
  nrow(
    plasma_axis_metadata
  ) ==
    0
) {
  stop(
    "No cells were retained for the plasma differentiation trajectory."
  )
}

plasma_group_counts <- plasma_axis_metadata |>
  dplyr::count(
    TB_Subtype,
    name = "nCells"
  ) |>
  tidyr::complete(
    TB_Subtype = factor(
      PLASMA_TRAJECTORY_GROUPS,
      levels = PLASMA_TRAJECTORY_GROUPS
    ),
    fill = list(
      nCells = 0L
    )
  )

print(
  plasma_group_counts
)

if (
  any(
    plasma_group_counts$nCells <
      50
  )
) {
  warning(
    "At least one trajectory group contains fewer than 50 cells:\n",
    paste(
      plasma_group_counts$TB_Subtype,
      plasma_group_counts$nCells,
      sep = "=",
      collapse = "; "
    )
  )
}

if (
  any(
    plasma_group_counts$nCells ==
      0
  )
) {
  stop(
    "One or more required trajectory groups are absent."
  )
}

project_plasma <- ArchR::subsetCells(
  ArchRProj = project_b,
  cellNames = plasma_axis_metadata$Cell
)

project_plasma <- ArchR::addCellColData(
  ArchRProj = project_plasma,
  data = as.character(
    plasma_axis_metadata$TB_Subtype
  ),
  cells = plasma_axis_metadata$Cell,
  name = "PlasmaTrajectoryGroup",
  force = TRUE
)

project_plasma$PlasmaTrajectoryGroup <- factor(
  project_plasma$PlasmaTrajectoryGroup,
  levels = PLASMA_TRAJECTORY_GROUPS
)

write.csv(
  plasma_group_counts,
  file.path(
    PLASMA_TRAJECTORY_TABLE_DIR,
    "PlasmaB_trajectory_group_cell_counts.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 3. Fit one supervised ArchR trajectory
# ==============================================================================

embedding_name <- detect_embedding(
  project_plasma
)

message(
  "Using embedding: ",
  embedding_name
)

project_plasma <- ArchR::addTrajectory(
  ArchRProj = project_plasma,
  name = TRAJECTORY_NAME,
  trajectory = PLASMA_TRAJECTORY_GROUPS,
  groupBy = "PlasmaTrajectoryGroup",
  embedding = embedding_name,
  preFilterQuantile = 0.95,
  postFilterQuantile = 0.95,
  useAll = FALSE,
  dof = 250,
  spar = 1,
  force = TRUE,
  seed = 123
)

cell_coldata_names <- colnames(
  ArchR::getCellColData(
    project_plasma
  )
)

if (
  !TRAJECTORY_NAME %in%
    cell_coldata_names
) {
  stop(
    "ArchR::addTrajectory() completed without adding ",
    TRAJECTORY_NAME,
    " to cellColData."
  )
}


# ==============================================================================
# 4. Plot trajectory
# ==============================================================================

trajectory_plot_list <- ArchR::plotTrajectory(
  ArchRProj = project_plasma,
  embedding = embedding_name,
  trajectory = TRAJECTORY_NAME,
  colorBy = "cellColData",
  name = TRAJECTORY_NAME,
  addArrow = TRUE,
  size = 0.18,
  baseSize = 9,
  keepAxis = FALSE
)

if (
  !is.list(
    trajectory_plot_list
  ) ||
    length(
      trajectory_plot_list
    ) <
      2
) {
  stop(
    "ArchR::plotTrajectory() did not return its expected two-plot list."
  )
}

p_trajectory_embedding <-
  trajectory_plot_list[[1]] +
  ggplot2::ggtitle(
    "Naive B → Intermediate B → Plasma B"
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold",
      hjust = 0.5
    )
  )

p_trajectory_diagnostic <-
  trajectory_plot_list[[2]] +
  ggplot2::ggtitle(
    "Pseudotime along the predefined B-cell continuum"
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold",
      hjust = 0.5
    )
  )

pdf(
  file = file.path(
    PLASMA_TRAJECTORY_PLOT_DIR,
    "Fig4C1_PlasmaB_ArchR_trajectory_embedding.pdf"
  ),
  width = 10,
  height = 5.5,
  useDingbats = FALSE
)

print(
  p_trajectory_embedding
)

print(
  p_trajectory_diagnostic
)

grDevices::dev.off()


# ==============================================================================
# 5. Export cell-level pseudotime
# ==============================================================================

cell_metadata_with_trajectory <- ArchR::getCellColData(
  project_plasma
) |>
  as.data.frame() |>
  tibble::rownames_to_column(
    "Cell"
  ) |>
  tibble::as_tibble()

required_export_columns <- intersect(
  c(
    "Cell",
    "Sample",
    "Donor",
    "Timepoint",
    "TB_Subtype",
    "PlasmaTrajectoryGroup",
    TRAJECTORY_NAME
  ),
  colnames(
    cell_metadata_with_trajectory
  )
)

pseudotime_table <-
  cell_metadata_with_trajectory |>
  dplyr::select(
    dplyr::all_of(
      required_export_columns
    )
  ) |>
  dplyr::rename(
    Pseudotime =
      dplyr::all_of(
        TRAJECTORY_NAME
      )
  )

write.csv(
  pseudotime_table,
  gzfile(
    file.path(
      PLASMA_TRAJECTORY_TABLE_DIR,
      "PlasmaB_ArchR_cell_pseudotime.csv.gz"
    )
  ),
  row.names = FALSE
)

pseudotime_qc <- pseudotime_table |>
  dplyr::group_by(
    PlasmaTrajectoryGroup
  ) |>
  dplyr::summarise(
    nCells = dplyr::n(),
    nFinitePseudotime = sum(
      is.finite(
        Pseudotime
      )
    ),
    FractionFinite = mean(
      is.finite(
        Pseudotime
      )
    ),
    MedianPseudotime = stats::median(
      Pseudotime,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

write.csv(
  pseudotime_qc,
  file.path(
    PLASMA_TRAJECTORY_TABLE_DIR,
    "PlasmaB_ArchR_pseudotime_QC.csv"
  ),
  row.names = FALSE
)

print(
  pseudotime_qc
)


# ==============================================================================
# 6. MotifMatrix variation along the predefined state continuum
# ==============================================================================

if (
  !"MotifMatrix" %in%
    ArchR::getAvailableMatrices(
      project_plasma
    )
) {
  stop(
    "MotifMatrix is absent from the plasma-axis project."
  )
}

trajectory_motif <- ArchR::getTrajectory(
  ArchRProj = project_plasma,
  name = TRAJECTORY_NAME,
  useMatrix = "MotifMatrix",
  groupEvery = GROUP_EVERY,
  log2Norm = FALSE,
  scaleTo = NULL,
  smoothWindow = SMOOTH_WINDOW,
  trajectoryLabel = "PlasmaTrajectoryGroup",
  threads = 8
)

p_motif_heatmap <- ArchR::plotTrajectoryHeatmap(
  seTrajectory = trajectory_motif,
  varCutOff = MOTIF_VAR_CUTOFF,
  maxFeatures = MOTIF_MAX_FEATURES,
  scaleRows = TRUE,
  limits = c(
    -1.5,
    1.5
  ),
  pal = ArchR::paletteContinuous(
    set = "solarExtra"
  ),
  colorColumns = TRUE,
  columnPal = c(
    "Naive_B" = "#4E79A7",
    "Intermediate_B" = "#E15759",
    "PlasmaB" = "#F28E2B"
  ),
  labelTop = MOTIF_LABEL_TOP,
  labelRows = FALSE,
  useSeqnames = "z",
  force = FALSE
)

pdf(
  file = file.path(
    PLASMA_TRAJECTORY_PLOT_DIR,
    "Fig4C2_PlasmaB_MotifMatrix_trajectory_heatmap.pdf"
  ),
  width = 7,
  height = 9,
  useDingbats = FALSE
)

ComplexHeatmap::draw(
  p_motif_heatmap
)

grDevices::dev.off()


# ==============================================================================
# 7. GeneExpressionMatrix variation along the predefined state continuum
# ==============================================================================

if (
  !"GeneExpressionMatrix" %in%
    ArchR::getAvailableMatrices(
      project_plasma
    )
) {
  stop(
    "GeneExpressionMatrix is absent from the plasma-axis project."
  )
}

trajectory_rna <- ArchR::getTrajectory(
  ArchRProj = project_plasma,
  name = TRAJECTORY_NAME,
  useMatrix = "GeneExpressionMatrix",
  groupEvery = GROUP_EVERY,
  log2Norm = TRUE,
  scaleTo = 10000,
  smoothWindow = SMOOTH_WINDOW,
  trajectoryLabel = "PlasmaTrajectoryGroup",
  threads = 8
)

available_trajectory_genes <- rownames(
  trajectory_rna
)

retained_plasma_markers <- intersect(
  PLASMA_MARKER_GENES,
  available_trajectory_genes
)

message(
  "Plasma-associated markers retained in RNA trajectory: ",
  paste(
    retained_plasma_markers,
    collapse = ", "
  )
)

write.csv(
  tibble::tibble(
    RequestedGene =
      PLASMA_MARKER_GENES,
    Retained =
      PLASMA_MARKER_GENES %in%
        retained_plasma_markers
  ),
  file.path(
    PLASMA_TRAJECTORY_TABLE_DIR,
    "PlasmaB_requested_RNA_markers.csv"
  ),
  row.names = FALSE
)

p_rna_heatmap <- ArchR::plotTrajectoryHeatmap(
  seTrajectory = trajectory_rna,
  varCutOff = RNA_VAR_CUTOFF,
  maxFeatures = RNA_MAX_FEATURES,
  scaleRows = TRUE,
  limits = c(
    -1.5,
    1.5
  ),
  pal = ArchR::paletteContinuous(
    set = "horizonExtra"
  ),
  colorColumns = TRUE,
  columnPal = c(
    "Naive_B" = "#4E79A7",
    "Intermediate_B" = "#E15759",
    "PlasmaB" = "#F28E2B"
  ),
  labelMarkers = retained_plasma_markers,
  labelTop = RNA_LABEL_TOP,
  labelRows = FALSE
)

pdf(
  file = file.path(
    PLASMA_TRAJECTORY_PLOT_DIR,
    "Fig4C3_PlasmaB_GeneExpression_trajectory_heatmap.pdf"
  ),
  width = 7,
  height = 9,
  useDingbats = FALSE
)

ComplexHeatmap::draw(
  p_rna_heatmap
)

grDevices::dev.off()


# ==============================================================================
# 8. Save objects
# ==============================================================================

saveRDS(
  project_plasma,
  file.path(
    PLASMA_TRAJECTORY_OBJECT_DIR,
    "proj_B_Plasma_ArchR_trajectory.rds"
  ),
  compress = FALSE
)

saveRDS(
  list(
    Parameters = list(
      TrajectoryGroups =
        PLASMA_TRAJECTORY_GROUPS,
      TrajectoryName =
        TRAJECTORY_NAME,
      EmbeddingName =
        embedding_name,
      GroupEvery =
        GROUP_EVERY,
      SmoothWindow =
        SMOOTH_WINDOW
    ),
    GroupCounts =
      plasma_group_counts,
    PseudotimeTable =
      pseudotime_table,
    PseudotimeQC =
      pseudotime_qc,
    TrajectoryMotif =
      trajectory_motif,
    TrajectoryRNA =
      trajectory_rna,
    RetainedPlasmaMarkers =
      retained_plasma_markers,
    pTrajectoryEmbedding =
      p_trajectory_embedding,
    pTrajectoryDiagnostic =
      p_trajectory_diagnostic,
    pMotifHeatmap =
      p_motif_heatmap,
    pRNAHeatmap =
      p_rna_heatmap
  ),
  file.path(
    PLASMA_TRAJECTORY_OBJECT_DIR,
    "Fig4C_PlasmaB_ArchR_trajectory_objects.rds"
  ),
  compress = FALSE
)

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    PLASMA_TRAJECTORY_OUTDIR,
    "sessionInfo.txt"
  )
)

cat(
  "\n============================================================\n",
  "ArchR supervised B-cell state continuum completed.\n",
  "\nTrajectory:\n",
  paste(
    PLASMA_TRAJECTORY_GROUPS,
    collapse = " -> "
  ),
  "\n\nEmbedding plot:\n",
  file.path(
    PLASMA_TRAJECTORY_PLOT_DIR,
    "Fig4C1_PlasmaB_ArchR_trajectory_embedding.pdf"
  ),
  "\n\nMotif trajectory heatmap:\n",
  file.path(
    PLASMA_TRAJECTORY_PLOT_DIR,
    "Fig4C2_PlasmaB_MotifMatrix_trajectory_heatmap.pdf"
  ),
  "\n\nRNA trajectory heatmap:\n",
  file.path(
    PLASMA_TRAJECTORY_PLOT_DIR,
    "Fig4C3_PlasmaB_GeneExpression_trajectory_heatmap.pdf"
  ),
  "\n============================================================\n",
  sep = ""
)

writeLines(
  c("Supervised path: Naive_B -> Intermediate_B -> PlasmaB.",
    "Memory_B is excluded from trajectory fitting.",
    "Directionality is predefined, not inferred from lineage tracing or RNA velocity.",
    "Pseudotime describes molecular variation along selected cell states; it is not vaccination time.",
    "The analysis does not establish in vivo differentiation or a validated lineage relationship."),
  file.path(PLASMA_TRAJECTORY_OUTDIR, "trajectory_interpretation.txt")
)
