#!/usr/bin/env Rscript

# ==============================================================================
# Figure 3A — donor-aware differential peak counts and Prime/Booster D1 overlap
#
# Biological questions
# --------------------
# 1. How many differential peaks are induced at:
#      Prime D1     : day1  vs day0
#      Booster D1   : day31 vs day30
#      Prime D14    : day14 vs day0
#      Prime D30    : day30 vs day0
#      Booster D14  : day44 vs day30
#
# 2. At D1, which peaks are:
#      Prime-specific
#      Booster-specific
#      Shared and opened in both stages
#      Shared and closed in both stages
#      Shared but changed in opposite directions
#
# Statistical design
# ------------------
# PeakMatrix pseudobulk counts are aggregated by:
#      TB_Subtype × biological sample
#
# Within every cell subtype:
#      Peak counts ~ Donor + Timepoint
#
# Prime contrasts are fitted together:
#      day1/day14/day30 versus day0
#
# Booster contrasts are fitted together:
#      day31/day44 versus day30
#
# The script uses ArchR::getGroupSE(
#     divideN = FALSE,
#     scaleTo = NULL
# )
# so that edgeR receives summed PeakMatrix counts rather than normalized values.
#
# Output
# ------
# Publication-level count figure containing all T/B subtypes.
# Significant peak tables and GRanges are saved for the next peak-to-gene step.
# ==============================================================================


# ==============================================================================
# 0. Limit hidden parallelism
# ==============================================================================

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

options(mc.cores = 1)


# ==============================================================================
# 1. Packages
# ==============================================================================

required_packages <- c(
  "ArchR",
  "Matrix",
  "edgeR",
  "dplyr",
  "tidyr",
  "tibble",
  "ggplot2",
  "patchwork",
  "scales",
  "SummarizedExperiment",
  "GenomicRanges",
  "S4Vectors"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  stop(
    "Missing packages: ",
    paste(missing_packages, collapse = ", ")
  )
}

if ("package:plyr" %in% search()) {
  detach("package:plyr", unload = FALSE, character.only = TRUE)
}

suppressPackageStartupMessages({
  library(ArchR)
  library(Matrix)
  library(edgeR)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(SummarizedExperiment)
  library(GenomicRanges)
  library(S4Vectors)
})

set.seed(123)

ARCHR_THREADS <- 16L
addArchRThreads(threads = ARCHR_THREADS)
addArchRGenome("hg38")
addArchRLocking(locking = TRUE)


# ==============================================================================
# 2. Input and output paths
# ==============================================================================

BASE_TB_DIR <- paste0(
  "/public/users/wangrong/home/wangrong/results/T_B/",
  "TB_final_downstream"
)

PROJECT_RDS_CANDIDATES <- c(
  file.path(
    BASE_TB_DIR,
    "Checkpoints",
    "99_proj_TB_down_complete.rds"
  ),
  file.path(
    BASE_TB_DIR,
    "Checkpoints",
    "04_proj_TB_down_peak2gene.rds"
  ),
  file.path(
    BASE_TB_DIR,
    "Checkpoints",
    "03_proj_TB_down_peak_matrix.rds"
  ),
  file.path(
    BASE_TB_DIR,
    "proj_TB_only_annotated_with_new_UMAP.rds"
  )
)

OUTDIR <- file.path(
  BASE_TB_DIR,
  "Figure3_DAR",
  "Fig3A_DAR_counts"
)

TABLE_DIR <- file.path(
  OUTDIR,
  "Tables"
)

PEAK_DIR <- file.path(
  OUTDIR,
  "Significant_peaks"
)

PLOT_DIR <- file.path(
  OUTDIR,
  "Plots"
)

dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PEAK_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)


# ==============================================================================
# 3. Analysis parameters
# ==============================================================================

CELLTYPE_COL <- "TB_Subtype"
SAMPLE_COL <- "Sample"
DONOR_COL <- "Donor"
TIME_COL <- "Timepoint"

# PlasmaB is rare, so 10 cells is used as the initial pseudobulk threshold.
# Increase to 20–25 for a stricter sensitivity analysis.
MIN_CELLS_PER_PSEUDOBULK <- 10L

MIN_PAIRED_DONORS <- 3L

# Independent peak filtering before edgeR fitting.
FILTER_MIN_COUNT <- 5
FILTER_MIN_TOTAL_COUNT <- 10

# Main DAR definition.
DAR_FDR <- 0.10
DAR_ABS_LOG2FC <- 0.50

# Save only significant peaks by default. Saving every tested peak for all
# 70 celltype/contrast combinations can produce a very large file.
SAVE_ALL_TESTED_RESULTS <- FALSE

# The pseudobulk matrix can be large and is not required for peak-to-gene.
SAVE_PSEUDOBULK_MATRIX <- FALSE

FIGURE_WIDTH_MM <- 180
FIGURE_HEIGHT_MM <- 245


# ==============================================================================
# 4. Cell subtype and contrast orders
# ==============================================================================

tb_levels <- c(
  "CD8_NaiveT",
  "CD4_NaiveT",
  "CD4_Tcm",
  "Treg",
  "Th2",
  "CD4_CD8_DP",
  "Th1",
  "gdT",
  "CD4_Teff",
  "CD8_Tem",
  "Naive_B",
  "Activated_B",
  "Resting_B",
  "PlasmaB"
)

contrast_table <- tribble(
  ~Comparison,    ~Stage,    ~Target,  ~Baseline, ~DisplayLabel,  ~ContrastOrder,
  "Prime_D1",     "Prime",   "day1",   "day0",    "Prime D1",      1L,
  "Booster_D1",   "Booster", "day31",  "day30",   "Booster D1",    2L,
  "Prime_D14",    "Prime",   "day14",  "day0",    "Prime D14",     3L,
  "Prime_D30",    "Prime",   "day30",  "day0",    "Prime D30",     4L,
  "Booster_D14",  "Booster", "day44",  "day30",   "Booster D14",   5L
)

comparison_levels <- contrast_table |>
  dplyr::arrange(ContrastOrder) |>
  dplyr::pull(DisplayLabel)

stage_timepoints <- list(
  Prime = c(
    "day0",
    "day1",
    "day14",
    "day30"
  ),
  Booster = c(
    "day30",
    "day31",
    "day44"
  )
)


# ==============================================================================
# 5. Helper functions
# ==============================================================================

safe_filename <- function(x) {
  gsub(
    "[^A-Za-z0-9_.-]",
    "_",
    x
  )
}


assign_lineage <- function(x) {
  dplyr::case_when(
    x %in% c(
      "Naive_B",
      "Activated_B",
      "Resting_B",
      "PlasmaB"
    ) ~ "B cells",
    TRUE ~ "T cells"
  )
}


derive_timepoint_from_sample <- function(sample_vector) {
  sample_index <- suppressWarnings(
    as.integer(
      sub(
        "^.*-",
        "",
        sample_vector
      )
    )
  )

  lookup <- c(
    "1" = "day0",
    "2" = "day1",
    "3" = "day3",
    "4" = "day6",
    "5" = "day14",
    "6" = "day30",
    "7" = "day31",
    "8" = "day33",
    "9" = "day36",
    "10" = "day44"
  )

  unname(
    lookup[
      as.character(sample_index)
    ]
  )
}


load_archr_project <- function() {
  object_candidates <- c(
    "proj_TB_down",
    "proj_TB_only"
  )

  for (object_name in object_candidates) {
    if (
      exists(
        object_name,
        envir = .GlobalEnv
      )
    ) {
      message(
        "Using ArchR project from active session: ",
        object_name
      )

      return(
        get(
          object_name,
          envir = .GlobalEnv
        )
      )
    }
  }

  project_file <- PROJECT_RDS_CANDIDATES[
    file.exists(
      PROJECT_RDS_CANDIDATES
    )
  ]

  if (length(project_file) == 0) {
    stop(
      "No ArchR project was found.\nChecked:\n",
      paste(
        PROJECT_RDS_CANDIDATES,
        collapse = "\n"
      )
    )
  }

  message(
    "Loading ArchR project:\n",
    project_file[1]
  )

  readRDS(
    project_file[1]
  )
}


make_peak_id <- function(gr) {
  if (
    is.null(gr) ||
      !inherits(gr, "GRanges") ||
      length(gr) == 0
  ) {
    stop(
      "make_peak_id() requires a non-empty GRanges object."
    )
  }

  paste0(
    as.character(
      GenomicRanges::seqnames(gr)
    ),
    ":",
    BiocGenerics::start(gr),
    "-",
    BiocGenerics::end(gr)
  )
}


extract_peak_ranges <- function(
    group_se,
    ArchRProj,
    count_matrix
) {
  n_features <- nrow(
    count_matrix
  )

  feature_names <- rownames(
    count_matrix
  )

  # ------------------------------------------------------------------------
  # Option 1: rowRanges(group_se)
  # ------------------------------------------------------------------------
  ranges_from_se <- tryCatch(
    SummarizedExperiment::rowRanges(
      group_se
    ),
    error = function(e) {
      NULL
    }
  )

  if (
    !is.null(
      ranges_from_se
    ) &&
      inherits(
        ranges_from_se,
        "GRanges"
      ) &&
      length(
        ranges_from_se
      ) ==
        n_features
  ) {
    message(
      "Peak coordinates obtained from rowRanges(group_se)."
    )

    return(
      ranges_from_se
    )
  }

  # ------------------------------------------------------------------------
  # Option 2: build GRanges from rowData(group_se)
  # ------------------------------------------------------------------------
  row_data <- tryCatch(
    as.data.frame(
      SummarizedExperiment::rowData(
        group_se
      )
    ),
    error = function(e) {
      NULL
    }
  )

  if (
    !is.null(
      row_data
    ) &&
      nrow(
        row_data
      ) ==
        n_features
  ) {
    chromosome_candidates <- c(
      "seqnames",
      "seqname",
      "chr",
      "chromosome",
      "Chr",
      "Chromosome"
    )

    start_candidates <- c(
      "start",
      "Start"
    )

    end_candidates <- c(
      "end",
      "End"
    )

    chromosome_column <- chromosome_candidates[
      chromosome_candidates %in%
        colnames(
          row_data
        )
    ]

    start_column <- start_candidates[
      start_candidates %in%
        colnames(
          row_data
        )
    ]

    end_column <- end_candidates[
      end_candidates %in%
        colnames(
          row_data
        )
    ]

    if (
      length(
        chromosome_column
      ) >
        0 &&
        length(
          start_column
        ) >
          0 &&
        length(
          end_column
        ) >
          0
    ) {
      ranges_from_rowdata <- GenomicRanges::GRanges(
        seqnames =
          as.character(
            row_data[[chromosome_column[1]]]
          ),
        ranges = IRanges::IRanges(
          start =
            as.integer(
              row_data[[start_column[1]]]
            ),
          end =
            as.integer(
              row_data[[end_column[1]]]
            )
        )
      )

      if (
        length(
          ranges_from_rowdata
        ) ==
          n_features
      ) {
        message(
          "Peak coordinates reconstructed from rowData(group_se)."
        )

        return(
          ranges_from_rowdata
        )
      }
    }
  }

  # ------------------------------------------------------------------------
  # Option 3: use the peak set stored in the ArchRProject.
  #
  # For PeakMatrix, the row order should correspond to getPeakSet().
  # We still validate and, when possible, reorder by feature names.
  # ------------------------------------------------------------------------
  project_peak_set <- tryCatch(
    ArchR::getPeakSet(
      ArchRProj
    ),
    error = function(e) {
      NULL
    }
  )

  if (
    is.null(
      project_peak_set
    ) ||
      !inherits(
        project_peak_set,
        "GRanges"
      ) ||
      length(
        project_peak_set
      ) ==
        0
  ) {
    stop(
      "rowRanges(group_se) is NULL and getPeakSet(ArchRProj) ",
      "did not return a valid GRanges peak set."
    )
  }

  if (
    length(
      project_peak_set
    ) !=
      n_features
  ) {
    stop(
      "PeakMatrix feature count does not match the project peak set: ",
      n_features,
      " matrix rows versus ",
      length(
        project_peak_set
      ),
      " peaks in getPeakSet(). ",
      "The PeakMatrix may have been generated from a different peak set."
    )
  }

  project_peak_ids <- make_peak_id(
    project_peak_set
  )

  # If matrix feature names contain genomic coordinates, align explicitly.
  if (
    !is.null(
      feature_names
    ) &&
      length(
        feature_names
      ) ==
        n_features &&
      all(
        feature_names %in%
          project_peak_ids
      )
  ) {
    project_peak_set <- project_peak_set[
      match(
        feature_names,
        project_peak_ids
      )
    ]

    message(
      "Peak coordinates obtained from getPeakSet() and aligned ",
      "to PeakMatrix row names."
    )

    return(
      project_peak_set
    )
  }

  # Some ArchR versions return simple feature indices rather than genomic
  # coordinates as row names. In that case, PeakMatrix and getPeakSet() are
  # both stored in peak-set order, so use the validated order.
  message(
    "rowRanges(group_se) was NULL. Peak coordinates obtained from ",
    "getPeakSet(); feature counts matched and project peak-set order was used."
  )

  project_peak_set
}


non_integer_fraction <- function(x) {
  values <- if (
    inherits(
      x,
      "sparseMatrix"
    )
  ) {
    x@x
  } else {
    as.numeric(x)
  }

  values <- values[
    is.finite(values)
  ]

  if (length(values) == 0) {
    return(0)
  }

  mean(
    abs(
      values -
        round(values)
    ) >
      1e-8
  )
}


make_design_and_contrasts <- function(
    sample_info,
    stage_name
) {
  stage_contrasts <- contrast_table |>
    dplyr::filter(
      Stage == stage_name
    )

  available_timepoints <- stage_timepoints[[stage_name]]

  sample_info <- sample_info |>
    dplyr::filter(
      Timepoint %in%
        available_timepoints
    ) |>
    dplyr::group_by(
      Donor
    ) |>
    dplyr::filter(
      dplyr::n_distinct(
        Timepoint
      ) >=
        2
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      Donor = factor(
        Donor
      ),
      Timepoint = factor(
        Timepoint,
        levels =
          available_timepoints[
            available_timepoints %in%
              unique(
                Timepoint
              )
          ]
      )
    ) |>
    dplyr::arrange(
      Donor,
      Timepoint,
      Sample
    )

  if (
    nrow(
      sample_info
    ) <
      6
  ) {
    return(NULL)
  }

  design <- model.matrix(
    ~ 0 + Timepoint + Donor,
    data = sample_info
  )

  if (
    qr(
      design
    )$rank <
      ncol(
        design
      )
  ) {
    return(
      list(
        status =
          "Rank-deficient design",
        sample_info =
          sample_info
      )
    )
  }

  timepoint_columns <- grep(
    "^Timepoint",
    colnames(
      design
    ),
    value = TRUE
  )

  timepoint_lookup <- setNames(
    timepoint_columns,
    sub(
      "^Timepoint",
      "",
      timepoint_columns
    )
  )

  contrast_vectors <- list()
  contrast_metadata <- list()
  contrast_index <- 1L

  for (
    i in seq_len(
      nrow(
        stage_contrasts
      )
    )
  ) {
    info <- stage_contrasts[
      i,
      ,
      drop = FALSE
    ]

    if (
      !all(
        c(
          info$Baseline,
          info$Target
        ) %in%
          names(
            timepoint_lookup
          )
      )
    ) {
      next
    }

    paired_donors <- sample_info |>
      dplyr::filter(
        as.character(
          Timepoint
        ) %in%
          c(
            info$Baseline,
            info$Target
          )
      ) |>
      dplyr::distinct(
        Donor,
        Timepoint
      ) |>
      dplyr::count(
        Donor,
        name = "n_timepoints"
      ) |>
      dplyr::filter(
        n_timepoints ==
          2
      ) |>
      dplyr::pull(
        Donor
      ) |>
      as.character()

    if (
      length(
        paired_donors
      ) <
        MIN_PAIRED_DONORS
    ) {
      next
    }

    contrast_vector <- rep(
      0,
      ncol(
        design
      )
    )

    names(
      contrast_vector
    ) <- colnames(
      design
    )

    contrast_vector[timepoint_lookup[[info$Target]]] <- 1

    contrast_vector[timepoint_lookup[[info$Baseline]]] <- -1

    contrast_vectors[[contrast_index]] <- contrast_vector

    contrast_metadata[[contrast_index]] <- info |>
      dplyr::mutate(
        nPairedDonors =
          length(
            paired_donors
          ),
        PairedDonors = paste(
          paired_donors,
          collapse = ","
        )
      )

    contrast_index <- contrast_index +
      1L
  }

  if (
    length(
      contrast_vectors
    ) ==
      0
  ) {
    return(NULL)
  }

  contrast_matrix <- do.call(
    cbind,
    contrast_vectors
  )

  metadata <- dplyr::bind_rows(
    contrast_metadata
  )

  colnames(
    contrast_matrix
  ) <- metadata$Comparison

  list(
    status = "OK",
    sample_info = sample_info,
    design = design,
    contrast_matrix = contrast_matrix,
    metadata = metadata
  )
}


run_subtype_stage_edger <- function(
    count_matrix,
    pb_metadata,
    subtype,
    stage_name,
    peak_annotation
) {
  sample_info_initial <- pb_metadata |>
    dplyr::filter(
      TB_Subtype == subtype,
      Timepoint %in%
        stage_timepoints[[stage_name]]
    )

  design_objects <- make_design_and_contrasts(
    sample_info =
      sample_info_initial,
    stage_name =
      stage_name
  )

  if (
    is.null(
      design_objects
    )
  ) {
    return(
      list(
        significant = tibble::tibble(),
        all_tested = tibble::tibble(),
        summary = tibble::tibble(
          TB_Subtype = subtype,
          Stage = stage_name,
          Comparison = NA_character_,
          Status = "Insufficient paired samples"
        )
      )
    )
  }

  if (
    design_objects$status !=
      "OK"
  ) {
    return(
      list(
        significant = tibble::tibble(),
        all_tested = tibble::tibble(),
        summary = tibble::tibble(
          TB_Subtype = subtype,
          Stage = stage_name,
          Comparison = NA_character_,
          Status = design_objects$status
        )
      )
    )
  }

  sample_info <- design_objects$sample_info
  design <- design_objects$design

  counts_sub <- count_matrix[
    ,
    sample_info$PBGroup,
    drop = FALSE
  ]

  stopifnot(
    identical(
      colnames(
        counts_sub
      ),
      sample_info$PBGroup
    )
  )

  positive_library <- Matrix::colSums(
    counts_sub
  ) >
    0

  counts_sub <- counts_sub[
    ,
    positive_library,
    drop = FALSE
  ]

  sample_info <- sample_info[
    positive_library,
    ,
    drop = FALSE
  ]

  design_objects <- make_design_and_contrasts(
    sample_info =
      sample_info,
    stage_name =
      stage_name
  )

  if (
    is.null(
      design_objects
    ) ||
      design_objects$status !=
        "OK"
  ) {
    return(
      list(
        significant = tibble::tibble(),
        all_tested = tibble::tibble(),
        summary = tibble::tibble(
          TB_Subtype = subtype,
          Stage = stage_name,
          Comparison = NA_character_,
          Status = "Insufficient samples after zero-library removal"
        )
      )
    )
  }

  sample_info <- design_objects$sample_info
  design <- design_objects$design

  counts_sub <- counts_sub[
    ,
    sample_info$PBGroup,
    drop = FALSE
  ]

  y <- edgeR::DGEList(
    counts =
      counts_sub
  )

  keep_peak <- edgeR::filterByExpr(
    y,
    design = design,
    min.count =
      FILTER_MIN_COUNT,
    min.total.count =
      FILTER_MIN_TOTAL_COUNT
  )

  if (
    sum(
      keep_peak
    ) <
      100
  ) {
    return(
      list(
        significant = tibble::tibble(),
        all_tested = tibble::tibble(),
        summary = tibble::tibble(
          TB_Subtype = subtype,
          Stage = stage_name,
          Comparison = NA_character_,
          Status = "Fewer than 100 peaks passed filtering"
        )
      )
    )
  }

  y <- y[
    keep_peak,
    ,
    keep.lib.sizes = FALSE
  ]

  # TMMwsp is stable for sparse count matrices with a high zero proportion.
  y <- edgeR::calcNormFactors(
    y,
    method = "TMMwsp"
  )

  y <- edgeR::estimateDisp(
    y,
    design,
    robust = TRUE
  )

  fit <- edgeR::glmQLFit(
    y,
    design,
    robust = TRUE
  )

  significant_list <- list()
  all_tested_list <- list()
  summary_list <- list()

  for (
    contrast_index in seq_len(
      ncol(
        design_objects$contrast_matrix
      )
    )
  ) {
    contrast_info <- design_objects$metadata[
      contrast_index,
      ,
      drop = FALSE
    ]

    contrast_vector <- design_objects$contrast_matrix[
      ,
      contrast_index
    ]

    qlf <- edgeR::glmQLFTest(
      fit,
      contrast =
        contrast_vector
    )

    result <- edgeR::topTags(
      qlf,
      n = Inf,
      sort.by = "none"
    )$table |>
      rownames_to_column(
        "PeakID"
      ) |>
      tibble::as_tibble() |>
      dplyr::left_join(
        peak_annotation,
        by = "PeakID"
      ) |>
      dplyr::mutate(
        TB_Subtype = subtype,
        Stage = stage_name,
        Comparison =
          contrast_info$Comparison,
        DisplayLabel =
          contrast_info$DisplayLabel,
        Baseline =
          contrast_info$Baseline,
        Target =
          contrast_info$Target,
        nPairedDonors =
          contrast_info$nPairedDonors,
        PairedDonors =
          contrast_info$PairedDonors,
        Significant = (
          FDR <=
            DAR_FDR &
            abs(
              logFC
            ) >=
              DAR_ABS_LOG2FC
        ),
        Direction = dplyr::case_when(
          Significant &
            logFC >
              0 ~
            "Opening",
          Significant &
            logFC <
              0 ~
            "Closing",
          TRUE ~
            "Not significant"
        )
      )

    significant_result <- result |>
      dplyr::filter(
        Significant
      )

    significant_list[[contrast_index]] <- significant_result

    if (
      SAVE_ALL_TESTED_RESULTS
    ) {
      all_tested_list[[contrast_index]] <- result
    }

    summary_list[[contrast_index]] <- tibble::tibble(
      TB_Subtype = subtype,
      Stage = stage_name,
      Comparison =
        contrast_info$Comparison,
      DisplayLabel =
        contrast_info$DisplayLabel,
      Baseline =
        contrast_info$Baseline,
      Target =
        contrast_info$Target,
      nPairedDonors =
        contrast_info$nPairedDonors,
      PairedDonors =
        contrast_info$PairedDonors,
      nSamplesInModel =
        nrow(
          sample_info
        ),
      nTestedPeaks =
        nrow(
          result
        ),
      nOpening =
        sum(
          result$Direction ==
            "Opening"
        ),
      nClosing =
        sum(
          result$Direction ==
            "Closing"
        ),
      nDAR =
        sum(
          result$Significant
        ),
      minFDR =
        min(
          result$FDR,
          na.rm = TRUE
        ),
      maxAbsLog2FC =
        max(
          abs(
            result$logFC
          ),
          na.rm = TRUE
        ),
      Status = "OK"
    )

    output_file <- file.path(
      PEAK_DIR,
      paste0(
        safe_filename(
          subtype
        ),
        "__",
        contrast_info$Comparison,
        "__significant_peaks.csv.gz"
      )
    )

    write.csv(
      significant_result,
      gzfile(
        output_file
      ),
      row.names = FALSE
    )
  }

  list(
    significant = dplyr::bind_rows(
      significant_list
    ),
    all_tested = dplyr::bind_rows(
      all_tested_list
    ),
    summary = dplyr::bind_rows(
      summary_list
    )
  )
}


# ==============================================================================
# 6. Load project and validate metadata
# ==============================================================================

proj_TB_fig3 <- load_archr_project()

available_matrices <- getAvailableMatrices(
  proj_TB_fig3
)

if (
  !"PeakMatrix" %in%
    available_matrices
) {
  stop(
    "PeakMatrix is absent. Load proj_TB_down after addPeakMatrix()."
  )
}

cell_colnames <- colnames(
  as.data.frame(
    getCellColData(
      proj_TB_fig3
    )
  )
)

if (
  !CELLTYPE_COL %in%
    cell_colnames
) {
  if (
    "TB_Subtype_All" %in%
      cell_colnames
  ) {
    subtype_values <- as.character(
      proj_TB_fig3$TB_Subtype_All
    )

    subtype_values[
      subtype_values ==
        "Activaed_B"
    ] <- "Activated_B"

    proj_TB_fig3 <- addCellColData(
      ArchRProj = proj_TB_fig3,
      data = subtype_values,
      cells = getCellNames(
        proj_TB_fig3
      ),
      name = CELLTYPE_COL,
      force = TRUE
    )

  } else {
    stop(
      "Neither TB_Subtype nor TB_Subtype_All exists in cellColData."
    )
  }
}

if (
  !SAMPLE_COL %in%
    cell_colnames
) {
  stop(
    "Sample is absent from cellColData."
  )
}

metadata_columns <- unique(
  c(
    CELLTYPE_COL,
    SAMPLE_COL,
    DONOR_COL,
    TIME_COL
  )
)

metadata_columns <- metadata_columns[
  metadata_columns %in%
    colnames(
      as.data.frame(
        getCellColData(
          proj_TB_fig3
        )
      )
    )
]

cell_metadata <- as.data.frame(
  getCellColData(
    proj_TB_fig3,
    select = metadata_columns
  )
) |>
  rownames_to_column(
    "Cell"
  )

if (
  !DONOR_COL %in%
    colnames(
      cell_metadata
    )
) {
  cell_metadata[[DONOR_COL]] <- sub(
    "-.*$",
    "",
    cell_metadata[[SAMPLE_COL]]
  )
}

if (
  !TIME_COL %in%
    colnames(
      cell_metadata
    )
) {
  cell_metadata[[TIME_COL]] <- derive_timepoint_from_sample(
    cell_metadata[[SAMPLE_COL]]
  )
}

cell_metadata <- cell_metadata |>
  dplyr::transmute(
    Cell = as.character(
      Cell
    ),
    TB_Subtype = as.character(
      .data[[CELLTYPE_COL]]
    ),
    Sample = as.character(
      .data[[SAMPLE_COL]]
    ),
    Donor = as.character(
      .data[[DONOR_COL]]
    ),
    Timepoint = as.character(
      .data[[TIME_COL]]
    )
  ) |>
  dplyr::mutate(
    TB_Subtype = dplyr::if_else(
      TB_Subtype ==
        "Activaed_B",
      "Activated_B",
      TB_Subtype
    )
  ) |>
  dplyr::filter(
    TB_Subtype %in%
      tb_levels,
    !is.na(
      Sample
    ),
    Sample !=
      "",
    !is.na(
      Donor
    ),
    Donor !=
      "",
    !is.na(
      Timepoint
    ),
    Timepoint %in%
      unique(
        unlist(
          stage_timepoints
        )
      )
  )

if (
  nrow(
    cell_metadata
  ) ==
    0
) {
  stop(
    "No valid T/B cells remained after metadata filtering."
  )
}

# Remove unused intermediate timepoints before pseudobulk aggregation.
# This also guarantees that every retained cell receives a DAR_PB_Group value.
proj_TB_fig3 <- subsetCells(
  ArchRProj = proj_TB_fig3,
  cellNames = cell_metadata$Cell
)

cell_metadata <- cell_metadata[
  match(
    getCellNames(proj_TB_fig3),
    cell_metadata$Cell
  ),
  ,
  drop = FALSE
]

stopifnot(
  identical(
    cell_metadata$Cell,
    getCellNames(proj_TB_fig3)
  )
)


# ==============================================================================
# 7. Build subtype × sample pseudobulk labels
# ==============================================================================

pb_metadata <- cell_metadata |>
  dplyr::count(
    TB_Subtype,
    Sample,
    Donor,
    Timepoint,
    name = "nCells"
  ) |>
  dplyr::arrange(
    factor(
      TB_Subtype,
      levels = tb_levels
    ),
    Donor,
    Timepoint,
    Sample
  ) |>
  dplyr::mutate(
    PBGroup = sprintf(
      "PB%04d",
      dplyr::row_number()
    ),
    KeepPseudobulk = (
      nCells >=
        MIN_CELLS_PER_PSEUDOBULK
    )
  )

cell_metadata <- cell_metadata |>
  dplyr::left_join(
    pb_metadata |>
      dplyr::select(
        TB_Subtype,
        Sample,
        Donor,
        Timepoint,
        PBGroup
      ),
    by = c(
      "TB_Subtype",
      "Sample",
      "Donor",
      "Timepoint"
    )
  )

stopifnot(
  !anyNA(
    cell_metadata$PBGroup
  )
)

proj_TB_fig3 <- addCellColData(
  ArchRProj = proj_TB_fig3,
  data =
    cell_metadata$PBGroup,
  cells =
    cell_metadata$Cell,
  name =
    "DAR_PB_Group",
  force = TRUE
)

write.csv(
  pb_metadata,
  file.path(
    TABLE_DIR,
    "DAR_pseudobulk_metadata_before_filtering.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 8. Aggregate PeakMatrix counts with ArchR::getGroupSE
# ==============================================================================

message(
  "Aggregating PeakMatrix by TB subtype × sample..."
)

peak_pb_se <- getGroupSE(
  ArchRProj = proj_TB_fig3,
  useMatrix = "PeakMatrix",
  groupBy = "DAR_PB_Group",
  divideN = FALSE,
  scaleTo = NULL,
  threads = ARCHR_THREADS,
  verbose = TRUE
)

peak_assay_name <- assayNames(
  peak_pb_se
)[1]

peak_pb_counts <- assay(
  peak_pb_se,
  peak_assay_name
)

peak_pb_counts <- methods::as(
  Matrix::Matrix(
    peak_pb_counts,
    sparse = TRUE
  ),
  "dgCMatrix"
)

peak_ranges <- extract_peak_ranges(
  group_se = peak_pb_se,
  ArchRProj = proj_TB_fig3,
  count_matrix = peak_pb_counts
)

stopifnot(
  inherits(
    peak_ranges,
    "GRanges"
  ),
  length(
    peak_ranges
  ) ==
    nrow(
      peak_pb_counts
    )
)

peak_ids <- make_peak_id(
  peak_ranges
)

rownames(
  peak_pb_counts
) <- peak_ids

peak_annotation <- tibble::tibble(
  PeakID = peak_ids,
  seqnames = as.character(
    GenomicRanges::seqnames(
      peak_ranges
    )
  ),
  start = BiocGenerics::start(
    peak_ranges
  ),
  end = BiocGenerics::end(
    peak_ranges
  ),
  width = BiocGenerics::width(
    peak_ranges
  )
)

if (
  anyDuplicated(
    peak_annotation$PeakID
  )
) {
  stop(
    "Duplicated PeakID values were found."
  )
}

pb_match <- match(
  colnames(
    peak_pb_counts
  ),
  pb_metadata$PBGroup
)

if (
  anyNA(
    pb_match
  )
) {
  stop(
    "Some getGroupSE columns could not be matched to pseudobulk metadata:\n",
    paste(
      colnames(
        peak_pb_counts
      )[
        is.na(
          pb_match
        )
      ],
      collapse = ", "
    )
  )
}

pb_metadata <- pb_metadata[
  pb_match,
  ,
  drop = FALSE
]

stopifnot(
  identical(
    pb_metadata$PBGroup,
    colnames(
      peak_pb_counts
    )
  )
)

keep_pb <- (
  pb_metadata$KeepPseudobulk &
    Matrix::colSums(
      peak_pb_counts
    ) >
      0
)

peak_pb_counts <- peak_pb_counts[
  ,
  keep_pb,
  drop = FALSE
]

pb_metadata <- pb_metadata[
  keep_pb,
  ,
  drop = FALSE
]

stopifnot(
  identical(
    pb_metadata$PBGroup,
    colnames(
      peak_pb_counts
    )
  )
)

fraction_non_integer <- non_integer_fraction(
  peak_pb_counts
)

validation_table <- tibble::tibble(
  nPeaks = nrow(
    peak_pb_counts
  ),
  nPseudobulks = ncol(
    peak_pb_counts
  ),
  minLibrarySize = min(
    Matrix::colSums(
      peak_pb_counts
    )
  ),
  medianLibrarySize = median(
    Matrix::colSums(
      peak_pb_counts
    )
  ),
  maxLibrarySize = max(
    Matrix::colSums(
      peak_pb_counts
    )
  ),
  nonIntegerFraction =
    fraction_non_integer,
  minimumCellsPerPseudobulk =
    MIN_CELLS_PER_PSEUDOBULK
)

print(
  validation_table
)

write.csv(
  validation_table,
  file.path(
    TABLE_DIR,
    "PeakMatrix_pseudobulk_validation.csv"
  ),
  row.names = FALSE
)

write.csv(
  pb_metadata,
  file.path(
    TABLE_DIR,
    "DAR_pseudobulk_metadata_used.csv"
  ),
  row.names = FALSE
)

if (
  fraction_non_integer >
    1e-6
) {
  stop(
    "The aggregated PeakMatrix contains non-integer values. ",
    "Confirm that getGroupSE used divideN = FALSE and scaleTo = NULL."
  )
}

if (
  SAVE_PSEUDOBULK_MATRIX
) {
  saveRDS(
    list(
      counts = peak_pb_counts,
      metadata = pb_metadata,
      peak_ranges = peak_ranges
    ),
    file.path(
      OUTDIR,
      "PeakMatrix_sample_pseudobulk_counts.rds"
    ),
    compress = FALSE
  )
}


# ==============================================================================
# 9. Differential peak analysis
# ==============================================================================

subtypes_present <- tb_levels[
  tb_levels %in%
    unique(
      pb_metadata$TB_Subtype
    )
]

significant_result_list <- list()
all_tested_result_list <- list()
summary_result_list <- list()

result_index <- 1L
summary_index <- 1L

for (
  subtype in subtypes_present
) {
  for (
    stage_name in c(
      "Prime",
      "Booster"
    )
  ) {
    message(
      "edgeR DAR: ",
      subtype,
      " — ",
      stage_name
    )

    stage_result <- run_subtype_stage_edger(
      count_matrix =
        peak_pb_counts,
      pb_metadata =
        pb_metadata,
      subtype =
        subtype,
      stage_name =
        stage_name,
      peak_annotation =
        peak_annotation
    )

    if (
      nrow(
        stage_result$significant
      ) >
        0
    ) {
      significant_result_list[[result_index]] <- stage_result$significant
    }

    if (
      SAVE_ALL_TESTED_RESULTS &&
        nrow(
          stage_result$all_tested
        ) >
          0
    ) {
      all_tested_result_list[[result_index]] <- stage_result$all_tested
    }

    summary_result_list[[summary_index]] <- stage_result$summary

    result_index <- result_index +
      1L

    summary_index <- summary_index +
      1L

    invisible(
      gc()
    )
  }
}

significant_peaks <- dplyr::bind_rows(
  significant_result_list
)

all_tested_peaks <- dplyr::bind_rows(
  all_tested_result_list
)

dar_summary_raw <- dplyr::bind_rows(
  summary_result_list
)


# ==============================================================================
# 10. Complete missing subtype × contrast combinations
# ==============================================================================

expected_grid <- tidyr::crossing(
  TB_Subtype = subtypes_present,
  contrast_table |>
    dplyr::select(
      Comparison,
      Stage,
      Target,
      Baseline,
      DisplayLabel,
      ContrastOrder
    )
)

dar_summary <- expected_grid |>
  dplyr::left_join(
    dar_summary_raw |>
      dplyr::filter(
        !is.na(
          Comparison
        )
      ),
    by = c(
      "TB_Subtype",
      "Comparison",
      "Stage",
      "Target",
      "Baseline",
      "DisplayLabel"
    )
  ) |>
  dplyr::mutate(
    Status = tidyr::replace_na(
      Status,
      "Not tested"
    ),
    nOpening = dplyr::if_else(
      Status ==
        "OK",
      as.integer(
        tidyr::replace_na(
          nOpening,
          0L
        )
      ),
      NA_integer_
    ),
    nClosing = dplyr::if_else(
      Status ==
        "OK",
      as.integer(
        tidyr::replace_na(
          nClosing,
          0L
        )
      ),
      NA_integer_
    ),
    nDAR = dplyr::if_else(
      Status ==
        "OK",
      as.integer(
        tidyr::replace_na(
          nDAR,
          0L
        )
      ),
      NA_integer_
    ),
    TB_Subtype = factor(
      TB_Subtype,
      levels =
        subtypes_present
    ),
    Lineage = assign_lineage(
      as.character(
        TB_Subtype
      )
    ),
    DisplayLabel = factor(
      DisplayLabel,
      levels =
        comparison_levels
    )
  ) |>
  dplyr::arrange(
    TB_Subtype,
    ContrastOrder
  )

write.csv(
  dar_summary,
  file.path(
    TABLE_DIR,
    "Fig3A_differential_peak_count_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  significant_peaks,
  gzfile(
    file.path(
      TABLE_DIR,
      "Fig3A_all_significant_DARs.csv.gz"
    )
  ),
  row.names = FALSE
)

if (
  SAVE_ALL_TESTED_RESULTS
) {
  write.csv(
    all_tested_peaks,
    gzfile(
      file.path(
        TABLE_DIR,
        "Fig3A_all_tested_peak_results.csv.gz"
      )
    ),
    row.names = FALSE
  )
}


# ==============================================================================
# 11. Prime D1 versus Booster D1 shared and unique peaks
# ==============================================================================

prime_d1 <- significant_peaks |>
  dplyr::filter(
    Comparison ==
      "Prime_D1"
  ) |>
  dplyr::select(
    TB_Subtype,
    PeakID,
    Prime_logFC = logFC,
    Prime_FDR = FDR
  )

booster_d1 <- significant_peaks |>
  dplyr::filter(
    Comparison ==
      "Booster_D1"
  ) |>
  dplyr::select(
    TB_Subtype,
    PeakID,
    Booster_logFC = logFC,
    Booster_FDR = FDR
  )

day1_overlap_details <- dplyr::full_join(
  prime_d1,
  booster_d1,
  by = c(
    "TB_Subtype",
    "PeakID"
  )
) |>
  dplyr::left_join(
    peak_annotation,
    by = "PeakID"
  ) |>
  dplyr::mutate(
    InPrime = !is.na(
      Prime_logFC
    ),
    InBooster = !is.na(
      Booster_logFC
    ),
    OverlapClass = dplyr::case_when(
      InPrime &
        !InBooster ~
        "Prime-specific",
      !InPrime &
        InBooster ~
        "Booster-specific",
      InPrime &
        InBooster &
        Prime_logFC >
          0 &
        Booster_logFC >
          0 ~
        "Shared opening",
      InPrime &
        InBooster &
        Prime_logFC <
          0 &
        Booster_logFC <
          0 ~
        "Shared closing",
      InPrime &
        InBooster &
        sign(
          Prime_logFC
        ) !=
          sign(
            Booster_logFC
          ) ~
        "Shared discordant",
      TRUE ~
        "Unclassified"
    )
  )

overlap_levels <- c(
  "Prime-specific",
  "Shared opening",
  "Shared closing",
  "Shared discordant",
  "Booster-specific"
)

day1_overlap_summary <- day1_overlap_details |>
  dplyr::count(
    TB_Subtype,
    OverlapClass,
    name = "nPeaks"
  ) |>
  tidyr::complete(
    TB_Subtype =
      subtypes_present,
    OverlapClass =
      overlap_levels,
    fill = list(
      nPeaks = 0L
    )
  ) |>
  dplyr::mutate(
    TB_Subtype = factor(
      TB_Subtype,
      levels =
        subtypes_present
    ),
    Lineage = assign_lineage(
      as.character(
        TB_Subtype
      )
    ),
    OverlapClass = factor(
      OverlapClass,
      levels =
        overlap_levels
    )
  )

write.csv(
  day1_overlap_details,
  gzfile(
    file.path(
      TABLE_DIR,
      "Fig3A_PrimeD1_BoosterD1_overlap_details.csv.gz"
    )
  ),
  row.names = FALSE
)

write.csv(
  day1_overlap_summary,
  file.path(
    TABLE_DIR,
    "Fig3A_PrimeD1_BoosterD1_overlap_counts.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 12. Save significant peaks as GRanges for peak-to-gene analysis
# ==============================================================================

significant_peak_granges <- peak_ranges[
  match(
    significant_peaks$PeakID,
    peak_annotation$PeakID
  )
]

mcols(
  significant_peak_granges
) <- S4Vectors::DataFrame(
  significant_peaks[
    ,
    setdiff(
      colnames(
        significant_peaks
      ),
      c(
        "seqnames",
        "start",
        "end",
        "width"
      )
    ),
    drop = FALSE
  ]
)

significant_peak_granges_split <- split(
  significant_peak_granges,
  paste(
    significant_peaks$TB_Subtype,
    significant_peaks$Comparison,
    sep = "__"
  )
)

saveRDS(
  significant_peak_granges_split,
  file.path(
    OUTDIR,
    "Fig3A_significant_DAR_GRanges_by_subtype_comparison.rds"
  ),
  compress = FALSE
)


# ==============================================================================
# 13. Publication-level plots
# ==============================================================================

count_plot_data <- dar_summary |>
  dplyr::filter(
    Status ==
      "OK"
  ) |>
  dplyr::select(
    TB_Subtype,
    Lineage,
    DisplayLabel,
    Opening =
      nOpening,
    Closing =
      nClosing
  ) |>
  tidyr::pivot_longer(
    cols = c(
      Opening,
      Closing
    ),
    names_to =
      "Direction",
    values_to =
      "nPeaks"
  ) |>
  dplyr::mutate(
    Direction = factor(
      Direction,
      levels = c(
        "Opening",
        "Closing"
      )
    ),
    TB_Subtype = factor(
      as.character(
        TB_Subtype
      ),
      levels = rev(
        subtypes_present
      )
    ),
    DisplayLabel = factor(
      as.character(
        DisplayLabel
      ),
      levels =
        comparison_levels
    ),
    Lineage = factor(
      Lineage,
      levels = c(
        "T cells",
        "B cells"
      )
    )
  )

p_dar_counts <- ggplot(
  count_plot_data,
  aes(
    x = DisplayLabel,
    y = TB_Subtype
  )
) +
  geom_tile(
    fill = "grey98",
    color = "white",
    linewidth = 0.25
  ) +
  geom_point(
    aes(
      size = nPeaks,
      fill = Direction
    ),
    shape = 21,
    color = "grey20",
    stroke = 0.30
  ) +
  geom_vline(
    xintercept = 2.5,
    color = "grey30",
    linewidth = 0.55,
    linetype = "dashed"
  ) +
  facet_grid(
    rows = vars(
      Lineage
    ),
    cols = vars(
      Direction
    ),
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  scale_fill_manual(
    values = c(
      "Opening" = "#C7473D",
      "Closing" = "#3F6FA3"
    ),
    guide = "none"
  ) +
  scale_size_area(
    max_size = 9,
    breaks = scales::breaks_pretty(
      n = 4
    ),
    name = "DAR count"
  ) +
  scale_x_discrete(
    drop = FALSE
  ) +
  scale_y_discrete(
    drop = FALSE
  ) +
  labs(
    tag = "A",
    title = "Temporal differential chromatin accessibility across T and B cell subtypes",
    subtitle = paste0(
      "Donor-adjusted pseudobulk edgeR; FDR \u2264 ",
      DAR_FDR,
      " and |log2FC| \u2265 ",
      DAR_ABS_LOG2FC
    ),
    x = NULL,
    y = NULL
  ) +
  theme_classic(
    base_size = 9
  ) +
  theme(
    axis.text.x = element_text(
      angle = 40,
      hjust = 1,
      color = "black",
      size = 7.5
    ),
    axis.text.y = element_text(
      color = "black",
      size = 7.5
    ),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    strip.placement = "outside",
    strip.background = element_rect(
      fill = "grey95",
      color = "grey45",
      linewidth = 0.40
    ),
    strip.text.x = element_text(
      face = "bold",
      color = "black",
      size = 9
    ),
    strip.text.y.left = element_text(
      angle = 0,
      face = "bold",
      color = "black",
      size = 8
    ),
    panel.spacing = grid::unit(
      0.18,
      "lines"
    ),
    panel.border = element_rect(
      color = "grey45",
      fill = NA,
      linewidth = 0.35
    ),
    legend.position = "right",
    legend.title = element_text(
      face = "bold"
    ),
    plot.title = element_text(
      face = "bold",
      size = 11
    ),
    plot.subtitle = element_text(
      size = 8.5
    ),
    plot.tag = element_text(
      face = "bold",
      size = 14
    )
  )


overlap_plot_data <- day1_overlap_summary |>
  dplyr::mutate(
    TB_Subtype = factor(
      as.character(
        TB_Subtype
      ),
      levels = rev(
        subtypes_present
      )
    ),
    Lineage = factor(
      Lineage,
      levels = c(
        "T cells",
        "B cells"
      )
    )
  )

p_day1_overlap <- ggplot(
  overlap_plot_data,
  aes(
    x = nPeaks,
    y = TB_Subtype,
    fill = OverlapClass
  )
) +
  geom_col(
    width = 0.72,
    color = "white",
    linewidth = 0.20
  ) +
  facet_grid(
    rows = vars(
      Lineage
    ),
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  scale_fill_manual(
    values = c(
      "Prime-specific" = "#E69F00",
      "Shared opening" = "#C7473D",
      "Shared closing" = "#2A9D8F",
      "Shared discordant" = "#8E6BBE",
      "Booster-specific" = "#3F6FA3"
    ),
    drop = FALSE,
    name = "D1 peak class"
  ) +
  scale_x_continuous(
    labels = scales::label_number(
      big.mark = ","
    ),
    expand = expansion(
      mult = c(
        0,
        0.04
      )
    )
  ) +
  labs(
    tag = "B",
    title = "Prime D1 and Booster D1 differential peaks: shared and stage-specific components",
    subtitle = "Shared peaks are separated by concordant opening, concordant closing, or discordant direction",
    x = "Number of differential peaks",
    y = NULL
  ) +
  theme_classic(
    base_size = 9
  ) +
  theme(
    axis.text.x = element_text(
      color = "black",
      size = 7.5
    ),
    axis.text.y = element_text(
      color = "black",
      size = 7.5
    ),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    strip.placement = "outside",
    strip.background = element_rect(
      fill = "grey95",
      color = "grey45",
      linewidth = 0.40
    ),
    strip.text.y.left = element_text(
      angle = 0,
      face = "bold",
      color = "black",
      size = 8
    ),
    panel.spacing.y = grid::unit(
      0.18,
      "lines"
    ),
    panel.border = element_rect(
      color = "grey45",
      fill = NA,
      linewidth = 0.35
    ),
    legend.position = "bottom",
    legend.title = element_text(
      face = "bold"
    ),
    plot.title = element_text(
      face = "bold",
      size = 11
    ),
    plot.subtitle = element_text(
      size = 8.5
    ),
    plot.tag = element_text(
      face = "bold",
      size = 14
    )
  )


p_figure3a <- p_dar_counts /
  p_day1_overlap +
  plot_layout(
    heights = c(
      1.25,
      1
    ),
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom"
  )

print(
  p_figure3a
)

ggsave(
  filename = file.path(
    PLOT_DIR,
    "Fig3A_DAR_counts_and_PrimeD1_BoosterD1_overlap.pdf"
  ),
  plot = p_figure3a,
  width = FIGURE_WIDTH_MM,
  height = FIGURE_HEIGHT_MM,
  units = "mm",
  limitsize = FALSE
)

ggsave(
  filename = file.path(
    PLOT_DIR,
    "Fig3A_DAR_counts_and_PrimeD1_BoosterD1_overlap.png"
  ),
  plot = p_figure3a,
  width = FIGURE_WIDTH_MM,
  height = FIGURE_HEIGHT_MM,
  units = "mm",
  dpi = 600,
  limitsize = FALSE
)

ggsave(
  filename = file.path(
    PLOT_DIR,
    "Fig3A_panelA_DAR_counts.pdf"
  ),
  plot = p_dar_counts,
  width = 180,
  height = 135,
  units = "mm",
  limitsize = FALSE
)

ggsave(
  filename = file.path(
    PLOT_DIR,
    "Fig3A_panelB_D1_overlap.pdf"
  ),
  plot = p_day1_overlap,
  width = 180,
  height = 125,
  units = "mm",
  limitsize = FALSE
)


# ==============================================================================
# 14. Save reusable objects for Figure 3B/3C peak-to-gene analysis
# ==============================================================================

saveRDS(
  list(
    parameters = list(
      MIN_CELLS_PER_PSEUDOBULK =
        MIN_CELLS_PER_PSEUDOBULK,
      MIN_PAIRED_DONORS =
        MIN_PAIRED_DONORS,
      FILTER_MIN_COUNT =
        FILTER_MIN_COUNT,
      FILTER_MIN_TOTAL_COUNT =
        FILTER_MIN_TOTAL_COUNT,
      DAR_FDR =
        DAR_FDR,
      DAR_ABS_LOG2FC =
        DAR_ABS_LOG2FC
    ),
    contrast_table =
      contrast_table,
    pseudobulk_metadata =
      pb_metadata,
    peak_annotation =
      peak_annotation,
    peak_ranges =
      peak_ranges,
    dar_summary =
      dar_summary,
    significant_peaks =
      significant_peaks,
    day1_overlap_details =
      day1_overlap_details,
    day1_overlap_summary =
      day1_overlap_summary,
    p_dar_counts =
      p_dar_counts,
    p_day1_overlap =
      p_day1_overlap,
    p_figure3a =
      p_figure3a
  ),
  file.path(
    OUTDIR,
    "Fig3A_DAR_count_and_overlap_objects.rds"
  ),
  compress = FALSE
)

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    OUTDIR,
    "sessionInfo.txt"
  )
)

cat(
  "\n============================================================\n",
  "Figure 3A DAR count analysis completed.\n",
  "\nMain combined figure:\n",
  file.path(
    PLOT_DIR,
    "Fig3A_DAR_counts_and_PrimeD1_BoosterD1_overlap.pdf"
  ),
  "\n\nDAR count table:\n",
  file.path(
    TABLE_DIR,
    "Fig3A_differential_peak_count_summary.csv"
  ),
  "\n\nSignificant DAR table:\n",
  file.path(
    TABLE_DIR,
    "Fig3A_all_significant_DARs.csv.gz"
  ),
  "\n\nGRanges for peak-to-gene:\n",
  file.path(
    OUTDIR,
    "Fig3A_significant_DAR_GRanges_by_subtype_comparison.rds"
  ),
  "\n============================================================\n",
  sep = ""
)
