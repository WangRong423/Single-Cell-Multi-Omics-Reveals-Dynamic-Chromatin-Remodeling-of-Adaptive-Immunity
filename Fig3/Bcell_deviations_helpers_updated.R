# ==============================================================================
# Shared helpers for B-cell analysis using proj_TB_deviations_complete.rds
# ==============================================================================

suppressPackageStartupMessages({
  library(ArchR)
  library(Matrix)
  library(limma)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(SummarizedExperiment)
  library(S4Vectors)
  library(GenomicRanges)
  library(BiocGenerics)
})

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

options(mc.cores = 1)

if ("package:plyr" %in% search()) {
  detach(
    "package:plyr",
    unload = FALSE,
    character.only = TRUE
  )
}

set.seed(123)

ArchR::addArchRThreads(threads = 8)
ArchR::addArchRGenome("hg38")
ArchR::addArchRLocking(locking = TRUE)

# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------

BASE_TB_DIR <- paste0(
  "/public/users/wangrong/home/wangrong/results/T_B/",
  "TB_final_downstream"
)

PROJECT_RDS_CANDIDATES <- c(
  file.path(
    getwd(),
    "proj_TB_deviations_complete.rds"
  ),
  file.path(
    BASE_TB_DIR,
    "proj_TB_deviations_complete.rds"
  ),
  file.path(
    BASE_TB_DIR,
    "Checkpoints",
    "proj_TB_deviations_complete.rds"
  )
)

BCELL_OUTDIR <- file.path(
  BASE_TB_DIR,
  "Figure4_Bcell_dynamics_deviations_updated_labels"
)

BCELL_TABLE_DIR <- file.path(
  BCELL_OUTDIR,
  "Tables"
)

BCELL_PLOT_DIR <- file.path(
  BCELL_OUTDIR,
  "Plots"
)

BCELL_OBJECT_DIR <- file.path(
  BCELL_OUTDIR,
  "Objects"
)

BCELL_TRAJECTORY_DIR <- file.path(
  BCELL_OUTDIR,
  "Trajectory"
)

for (
  output_directory in c(
    BCELL_TABLE_DIR,
    BCELL_PLOT_DIR,
    BCELL_OBJECT_DIR,
    BCELL_TRAJECTORY_DIR
  )
) {
  dir.create(
    output_directory,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

# ------------------------------------------------------------------------------
# Biological order
# ------------------------------------------------------------------------------

B_SUBTYPES <- c(
  "Naive_B",
  "Intermediate_B",
  "Memory_B",
  "PlasmaB"
)

TIME_LEVELS <- c(
  "day0",
  "day1",
  "day3",
  "day6",
  "day14",
  "day30",
  "day31",
  "day33",
  "day36",
  "day44"
)

TIME_LABELS <- c(
  "D0",
  "D1",
  "D3",
  "D6",
  "D14",
  "D30",
  "B+D1",
  "B+D3",
  "B+D6",
  "B+D14"
)

B_SUBTYPE_COLORS <- c(
  "Naive_B" = "#4E79A7",
  "Intermediate_B" = "#E15759",
  "Memory_B" = "#59A14F",
  "PlasmaB" = "#F28E2B"
)

MIN_CELLS_PER_PSEUDOBULK <- 10L

# ------------------------------------------------------------------------------
# Project and metadata
# ------------------------------------------------------------------------------

load_deviations_project <- function() {
  if (
    exists(
      "proj_TB_deviations_complete",
      envir = .GlobalEnv
    )
  ) {
    message(
      "Using active object: proj_TB_deviations_complete"
    )

    project <- get(
      "proj_TB_deviations_complete",
      envir = .GlobalEnv
    )

  } else {
    project_paths <- PROJECT_RDS_CANDIDATES[
      file.exists(
        PROJECT_RDS_CANDIDATES
      )
    ]

    if (length(project_paths) == 0) {
      stop(
        "proj_TB_deviations_complete.rds was not found.\nChecked:\n",
        paste(
          PROJECT_RDS_CANDIDATES,
          collapse = "\n"
        ),
        "\nEdit PROJECT_RDS_CANDIDATES in Bcell_deviations_helpers.R ",
        "if the file is stored elsewhere."
      )
    }

    message(
      "Loading deviations-complete project:\n",
      normalizePath(
        project_paths[1]
      )
    )

    project <- readRDS(
      project_paths[1]
    )
  }

  available_matrices <- ArchR::getAvailableMatrices(
    project
  )

  message(
    "Available matrices: ",
    paste(
      available_matrices,
      collapse = ", "
    )
  )

  if (
    !"MotifMatrix" %in%
      available_matrices
  ) {
    stop(
      "The loaded project does not contain MotifMatrix."
    )
  }

  project
}


derive_timepoint <- function(sample_vector) {
  sample_suffix <- sub(
    "^.*-",
    "",
    as.character(
      sample_vector
    )
  )

  suffix_map <- c(
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
    suffix_map[
      sample_suffix
    ]
  )
}


derive_donor <- function(sample_vector) {
  sub(
    "-[0-9]+$",
    "",
    as.character(
      sample_vector
    )
  )
}


standardize_b_metadata <- function(project) {
  metadata <- ArchR::getCellColData(
    project
  ) |>
    as.data.frame() |>
    tibble::rownames_to_column(
      "Cell"
    ) |>
    tibble::as_tibble()

  if (
    !"TB_Subtype" %in%
      colnames(
        metadata
      )
  ) {
    if (
      "TB_Subtype_All" %in%
        colnames(
          metadata
        )
    ) {
      metadata$TB_Subtype <-
        metadata$TB_Subtype_All
    } else {
      stop(
        "No TB_Subtype or TB_Subtype_All column was found."
      )
    }
  }

  metadata$TB_Subtype <- as.character(
    metadata$TB_Subtype
  )

  metadata$TB_Subtype_original <- metadata$TB_Subtype
  label_map <- c(
    Activaed_B = "Intermediate_B",
    Activated_B = "Intermediate_B",
    Resting_B = "Memory_B"
  )
  hit <- match(metadata$TB_Subtype, names(label_map))
  metadata$TB_Subtype[!is.na(hit)] <- unname(label_map[hit[!is.na(hit)]])
  mapping <- unique(metadata[, c("TB_Subtype_original", "TB_Subtype")])
  bmapping <- mapping[mapping$TB_Subtype %in% B_SUBTYPES, , drop = FALSE]
  if (anyDuplicated(bmapping$TB_Subtype)) {
    print(bmapping)
    stop("Multiple existing labels map to one B subtype; verify membership before merging.")
  }
  write.csv(mapping, file.path(BCELL_TABLE_DIR, "Subtype_label_mapping.csv"), row.names=FALSE)

  if (
    !"Sample" %in%
      colnames(
        metadata
      )
  ) {
    stop(
      "Cell metadata lacks Sample."
    )
  }

  if (
    !"Donor" %in%
      colnames(
        metadata
      )
  ) {
    metadata$Donor <- derive_donor(
      metadata$Sample
    )
  }

  if (
    !"Timepoint" %in%
      colnames(
        metadata
      )
  ) {
    metadata$Timepoint <- derive_timepoint(
      metadata$Sample
    )
  }

  metadata |>
    dplyr::mutate(
      Sample = as.character(
        Sample
      ),
      Donor = as.character(
        Donor
      ),
      Timepoint = as.character(
        Timepoint
      )
    ) |>
    dplyr::filter(
      TB_Subtype %in%
        B_SUBTYPES,
      Timepoint %in%
        TIME_LEVELS
    )
}


prepare_b_project <- function(project) {
  metadata <- standardize_b_metadata(
    project
  )

  if (nrow(metadata) == 0) {
    stop(
      "No B cells were retained."
    )
  }

  project_b <- ArchR::subsetCells(
    ArchRProj = project,
    cellNames = metadata$Cell
  )

  project_b <- ArchR::addCellColData(
    ArchRProj = project_b, data = metadata$TB_Subtype,
    cells = metadata$Cell, name = "TB_Subtype", force = TRUE
  )
  # One sample per donor/time/subtype is required for donor-level SEM.
  sample_check <- metadata |>
    dplyr::distinct(TB_Subtype, Donor, Timepoint, Sample) |>
    dplyr::count(TB_Subtype, Donor, Timepoint, name="nSamples")
  if (any(sample_check$nSamples > 1L)) {
    stop("Multiple samples per donor/time/subtype: aggregate biological replicates explicitly before donor SEM.")
  }

  pseudobulk_group <- paste(
    metadata$TB_Subtype,
    metadata$Sample,
    sep = "__"
  )

  names(
    pseudobulk_group
  ) <- metadata$Cell

  project_b <- ArchR::addCellColData(
    ArchRProj = project_b,
    data = pseudobulk_group,
    cells = metadata$Cell,
    name = "B_PB_Group",
    force = TRUE
  )

  pseudobulk_metadata <- metadata |>
    dplyr::mutate(
      B_PB_Group = paste(
        TB_Subtype,
        Sample,
        sep = "__"
      )
    ) |>
    dplyr::count(
      B_PB_Group,
      TB_Subtype,
      Sample,
      Donor,
      Timepoint,
      name = "nCells"
    ) |>
    dplyr::mutate(
      Keep = nCells >=
        MIN_CELLS_PER_PSEUDOBULK,
      TB_Subtype = factor(
        TB_Subtype,
        levels = B_SUBTYPES
      ),
      Timepoint = factor(
        Timepoint,
        levels = TIME_LEVELS
      )
    ) |>
    dplyr::arrange(
      TB_Subtype,
      Timepoint,
      Donor
    )

  list(
    Project = project_b,
    CellMetadata = metadata,
    PBMetadata = pseudobulk_metadata
  )
}

# ------------------------------------------------------------------------------
# SummarizedExperiment utilities
# ------------------------------------------------------------------------------

extract_assay <- function(
    summarized_experiment,
    preferred_assays = character(0)
) {
  assay_names <- SummarizedExperiment::assayNames(
    summarized_experiment
  )

  if (length(assay_names) == 0) {
    stop(
      "SummarizedExperiment has no assay."
    )
  }

  preferred_hit <- preferred_assays[
    preferred_assays %in%
      assay_names
  ]

  assay_name <- if (
    length(preferred_hit) >
      0
  ) {
    preferred_hit[1]
  } else {
    assay_names[1]
  }

  message(
    "Using assay: ",
    assay_name
  )

  SummarizedExperiment::assay(
    summarized_experiment,
    assay_name
  )
}


extract_feature_names <- function(
    summarized_experiment,
    matrix_object = NULL
) {
  if (is.null(matrix_object)) {
    matrix_object <- extract_assay(
      summarized_experiment
    )
  }

  matrix_names <- rownames(
    matrix_object
  )

  if (
    !is.null(matrix_names) &&
      length(matrix_names) ==
        nrow(
          matrix_object
        ) &&
      all(
        !is.na(matrix_names) &
          matrix_names !=
            ""
      )
  ) {
    return(
      as.character(
        matrix_names
      )
    )
  }

  row_data <- SummarizedExperiment::rowData(
    summarized_experiment
  ) |>
    as.data.frame()

  candidate_columns <- c(
    "name",
    "Name",
    "motifName",
    "MotifName",
    "geneName",
    "GeneName",
    "symbol",
    "Symbol"
  )

  selected_column <- candidate_columns[
    candidate_columns %in%
      colnames(
        row_data
      )
  ]

  if (length(selected_column) == 0) {
    stop(
      "No usable feature-name column was found in rowData(). ",
      "Available columns: ",
      paste(
        colnames(
          row_data
        ),
        collapse = ", "
      )
    )
  }

  as.character(
    row_data[[selected_column[1]]]
  )
}


get_bcell_pseudobulk <- function(
    project_b,
    pseudobulk_metadata,
    use_matrix,
    cache_file,
    divide_n,
    scale_to,
    preferred_assays = character(0)
) {
  # Recompute from the current project; never load stale annotation caches here.

  pseudobulk_se <- ArchR::getGroupSE(
    ArchRProj = project_b,
    useMatrix = use_matrix,
    groupBy = "B_PB_Group",
    divideN = divide_n,
    scaleTo = scale_to,
    threads = 8,
    verbose = TRUE
  )

  pseudobulk_matrix <- extract_assay(
    pseudobulk_se,
    preferred_assays = preferred_assays
  )

  pseudobulk_match <- match(
    colnames(
      pseudobulk_matrix
    ),
    pseudobulk_metadata$B_PB_Group
  )

  if (anyNA(pseudobulk_match)) {
    stop(
      use_matrix,
      " pseudobulk columns do not match B-cell metadata."
    )
  }

  matched_metadata <- pseudobulk_metadata[
    pseudobulk_match,
    ,
    drop = FALSE
  ]

  finite_column <- vapply(
    seq_len(
      ncol(
        pseudobulk_matrix
      )
    ),
    function(column_index) {
      any(
        is.finite(
          pseudobulk_matrix[
            ,
            column_index
          ]
        )
      )
    },
    logical(1)
  )

  keep_columns <- (
    matched_metadata$Keep &
      finite_column
  )

  pseudobulk_matrix <- pseudobulk_matrix[
    ,
    keep_columns,
    drop = FALSE
  ]

  matched_metadata <- matched_metadata[
    keep_columns,
    ,
    drop = FALSE
  ]

  stopifnot(
    identical(
      matched_metadata$B_PB_Group,
      colnames(
        pseudobulk_matrix
      )
    )
  )

  output <- list(
    Matrix = pseudobulk_matrix,
    Metadata = matched_metadata,
    SE = pseudobulk_se,
    AssayNames =
      SummarizedExperiment::assayNames(
        pseudobulk_se
      )
  )

  saveRDS(
    output,
    cache_file,
    compress = FALSE
  )

  output
}

# ------------------------------------------------------------------------------
# Numeric utilities
# ------------------------------------------------------------------------------

row_zscore <- function(matrix_object) {
  matrix_object <- as.matrix(
    matrix_object
  )

  row_mean <- rowMeans(
    matrix_object,
    na.rm = TRUE
  )

  row_sd <- apply(
    matrix_object,
    1,
    stats::sd,
    na.rm = TRUE
  )

  row_sd[
    !is.finite(row_sd) |
      row_sd ==
        0
  ] <- 1

  sweep(
    sweep(
      matrix_object,
      1,
      row_mean,
      "-"
    ),
    1,
    row_sd,
    "/"
  )
}


mean_sem <- function(
    data,
    value,
    groups
) {
  data |>
    dplyr::group_by(
      dplyr::across(
        dplyr::all_of(
          groups
        )
      )
    ) |>
    dplyr::summarise(
      Mean = mean(
        .data[[value]],
        na.rm = TRUE
      ),
      SD = stats::sd(
        .data[[value]],
        na.rm = TRUE
      ),
      N = sum(
        is.finite(
          .data[[value]]
        )
      ),
      SEM = dplyr::if_else(
        N >
          1,
        SD /
          sqrt(
            N
          ),
        NA_real_
      ),
      .groups = "drop"
    )
}


clean_tf_name <- function(feature_names) {
  feature_names <- as.character(
    feature_names
  )

  feature_names <- sub(
    "^(z|deviations|Motif):",
    "",
    feature_names,
    ignore.case = TRUE
  )

  feature_names <- sub(
    "_[0-9]+$",
    "",
    feature_names
  )

  feature_names <- sub(
    "\\.$",
    "",
    feature_names
  )

  feature_names
}


collapse_motifs_to_tf <- function(
    motif_matrix,
    motif_names
) {
  if (
    nrow(
      motif_matrix
    ) !=
      length(
        motif_names
      )
  ) {
    stop(
      "Motif matrix rows and motif names differ in length."
    )
  }

  tf_names <- clean_tf_name(
    motif_names
  )

  valid_rows <- (
    !is.na(tf_names) &
      tf_names !=
        ""
  )

  motif_matrix <- motif_matrix[
    valid_rows,
    ,
    drop = FALSE
  ]

  tf_names <- tf_names[
    valid_rows
  ]

  tf_index <- split(
    seq_len(
      nrow(
        motif_matrix
      )
    ),
    tf_names
  )

  tf_matrix <- do.call(
    rbind,
    lapply(
      tf_index,
      function(row_index) {
        colMeans(
          motif_matrix[
            row_index,
            ,
            drop = FALSE
          ],
          na.rm = TRUE
        )
      }
    )
  )

  rownames(
    tf_matrix
  ) <- names(
    tf_index
  )

  colnames(
    tf_matrix
  ) <- colnames(
    motif_matrix
  )

  tf_matrix
}


run_paired_limma <- function(
    feature_matrix,
    metadata,
    feature_name_column,
    contrast_table
) {
  result_list <- list()
  result_index <- 1L

  for (subtype in B_SUBTYPES) {
    subtype_metadata <- metadata |>
      dplyr::filter(
        as.character(
          TB_Subtype
        ) ==
          subtype
      )

    for (
      contrast_index in seq_len(
        nrow(
          contrast_table
        )
      )
    ) {
      contrast_info <- contrast_table[
        contrast_index,
        ,
        drop = FALSE
      ]

      paired_donors <- subtype_metadata |>
        dplyr::filter(
          as.character(
            Timepoint
          ) %in%
            c(
              contrast_info$Baseline,
              contrast_info$Target
            )
        ) |>
        dplyr::distinct(
          Donor,
          Timepoint
        ) |>
        dplyr::count(
          Donor,
          name = "nTimepoints"
        ) |>
        dplyr::filter(
          nTimepoints ==
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
          3
      ) {
        next
      }

      analysis_metadata <- subtype_metadata |>
        dplyr::filter(
          Donor %in%
            paired_donors,
          as.character(
            Timepoint
          ) %in%
            c(
              contrast_info$Baseline,
              contrast_info$Target
            )
        ) |>
        dplyr::mutate(
          Condition = factor(
            as.character(
              Timepoint
            ),
            levels = c(
              contrast_info$Baseline,
              contrast_info$Target
            )
          ),
          Donor = factor(
            Donor
          )
        ) |>
        dplyr::arrange(
          Donor,
          Condition
        )

      analysis_matrix <- feature_matrix[
        ,
        analysis_metadata$B_PB_Group,
        drop = FALSE
      ]

      design <- stats::model.matrix(
        ~ Donor + Condition,
        data = analysis_metadata
      )

      target_coefficient <- grep(
        "^Condition",
        colnames(
          design
        ),
        value = TRUE
      )

      if (
        length(
          target_coefficient
        ) !=
          1
      ) {
        next
      }

      fit <- limma::lmFit(
        analysis_matrix,
        design
      )

      fit <- limma::eBayes(
        fit,
        trend = TRUE,
        robust = FALSE
      )

      result_table <- limma::topTable(
        fit,
        coef = target_coefficient,
        number = Inf,
        sort.by = "none",
        adjust.method = "BH"
      ) |>
        tibble::rownames_to_column(
          feature_name_column
        ) |>
        tibble::as_tibble() |>
        dplyr::mutate(
          TB_Subtype = subtype,
          Comparison =
            contrast_info$Comparison,
          DisplayLabel =
            contrast_info$DisplayLabel,
          Target =
            contrast_info$Target,
          Baseline =
            contrast_info$Baseline,
          nPairedDonors =
            length(
              paired_donors
            ),
          PairedDonors = paste(
            paired_donors,
            collapse = ","
          )
        )

      result_list[[result_index]] <-
        result_table

      result_index <- result_index + 1L
    }
  }

  dplyr::bind_rows(
    result_list
  )
}


detect_embedding <- function(project) {
  embedding_names <- tryCatch(
    names(
      project@embeddings
    ),
    error = function(e) {
      character(0)
    }
  )

  candidates <- c(
    "UMAP_Combined",
    "Combined_UMAP",
    "UMAPCombined",
    "UMAP"
  )

  selected <- candidates[
    candidates %in%
      embedding_names
  ]

  if (length(selected) == 0) {
    stop(
      "No supported UMAP embedding was found. Available embeddings: ",
      paste(
        embedding_names,
        collapse = ", "
      )
    )
  }

  selected[1]
}
