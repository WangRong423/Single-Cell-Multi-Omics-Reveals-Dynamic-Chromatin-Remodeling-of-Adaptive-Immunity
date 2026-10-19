# Locate the companion file for source() and Rscript.
.local_source_file <- tryCatch(sys.frame(1)$ofile, error=function(e) NULL)
.local_args <- grep("^--file=", commandArgs(FALSE), value=TRUE)
.local_script_dir <- if (!is.null(.local_source_file)) {
  dirname(normalizePath(.local_source_file))
} else if (length(.local_args)) {
  dirname(normalizePath(sub("^--file=", "", .local_args[1])))
} else getwd()
source(file.path(.local_script_dir, "Tcell_updated_labels_helpers.R"))

###############################################################################
# T/B cell longitudinal chromatin remodeling pipeline
#
# Main analyses
#   1. T-cell TCellSI functional-state scoring using donor x time x subtype
#      RNA pseudobulk samples.
#   2. Paired day44-vs-day30 TF motif activity and TF RNA expression analysis.
#   3. Dynamic trajectories of day44-upregulated TFs across all ten timepoints.
#   4. B-cell and T-cell day44-vs-day30 significant peak analysis and display.
#   5. B-cell TF footprinting for selected day44-upregulated motifs.
#
# Input
#   An ArchRProject containing:
#     GeneExpressionMatrix, PeakMatrix, MotifMatrix, Motif annotation positions.
#
# Important
#   Only edit the CONFIGURATION section before running.
###############################################################################

options(stringsAsFactors = FALSE)
set.seed(20260731)

suppressPackageStartupMessages({
  library(ArchR)
  library(Matrix)
  library(SummarizedExperiment)
  library(GenomicRanges)
  library(S4Vectors)
  library(TCellSI)
  library(edgeR)
  library(limma)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
})

###############################################################################
# 0. CONFIGURATION: EDIT ONLY THIS SECTION
###############################################################################

# Existing ArchRProject object in the R environment.
PROJ <- select_tb_project()

# Metadata column names in PROJ@cellColData.
COL_DONOR     <- "Donor"
COL_TIME      <- "Timepoint"
COL_T_SUBTYPE <- "T_subtype"
COL_B_SUBTYPE <- "B_subtype"

# Ten longitudinal timepoints.
TIME_LEVELS <- c(
  "day0", "day1", "day3", "day6", "day14",
  "day30", "day31", "day33", "day36", "day44"
)

# Booster baseline and late booster endpoint.
BASELINE_TIME <- "day30"
TARGET_TIME   <- "day44"

# Minimum cells in each donor x time x subtype pseudobulk.
MIN_CELLS_RNA   <- 30
MIN_CELLS_MOTIF <- 30
MIN_CELLS_PEAK  <- 30

# Significance thresholds.
TF_FDR         <- 0.10
TF_MIN_DELTA_Z <- 0.20
RNA_TF_FDR     <- 0.10
PEAK_FDR       <- 0.10
PEAK_MIN_LOGFC <- 0.50

# Numbers of displayed features.
TOP_TF_PER_SUBTYPE             <- 8
TOP_PEAKS_PER_DIRECTION        <- 20
TOP_PEAK_LINES_PER_SUBTYPE     <- 6
TOP_FOOTPRINT_MOTIFS_PER_BTYPE <- 3

# Footprint settings.
RUN_B_FOOTPRINTS <- TRUE
FOOTPRINT_TIMES <- c(
  "day0", "day1", "day14",
  "day30", "day31", "day36", "day44"
)
MIN_FOOTPRINT_CELLS <- 25

# Output directory.
OUTDIR <- file.path(UPDATED_BASE, "T_B_TCellSI_TF_dynamic_results_updated_labels")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

###############################################################################
# 1. GENERAL UTILITIES
###############################################################################

`%||%` <- function(x, y) {
  if (!is.null(x) && length(x) > 0L) x else y
}

safe_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  gsub("_+", "_", x)
}

check_columns <- function(df, columns, object_name = "metadata") {
  missing_columns <- setdiff(columns, colnames(df))
  if (length(missing_columns) > 0L) {
    stop(
      object_name, " is missing columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }
}

first_existing_column <- function(df, candidates) {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) == 0L) return(NULL)
  hit[[1L]]
}

sem <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) <= 1L) return(NA_real_)
  stats::sd(x) / sqrt(length(x))
}

scale_rows_safely <- function(x) {
  x <- as.numeric(x)
  if (sum(is.finite(x)) < 2L || stats::sd(x, na.rm = TRUE) == 0) {
    return(rep(0, length(x)))
  }
  as.numeric(scale(x))
}

get_project_metadata <- function(proj) {
  meta <- as.data.frame(getCellColData(proj))
  meta$cell <- rownames(meta)
  if (!"TB_Subtype" %in% names(meta)) stop("TB_Subtype is required.")
  meta$TB_Subtype <- rename_tb_labels(meta$TB_Subtype,
    file.path(OUTDIR, "Subtype_label_mapping.csv"))
  meta[[COL_T_SUBTYPE]] <- ifelse(meta$TB_Subtype %in% UPDATED_T_TYPES, meta$TB_Subtype, NA_character_)
  meta[[COL_B_SUBTYPE]] <- ifelse(meta$TB_Subtype %in% UPDATED_B_TYPES, meta$TB_Subtype, NA_character_)

  check_columns(
    meta,
    c(COL_DONOR, COL_TIME, COL_T_SUBTYPE, COL_B_SUBTYPE),
    "ArchR cellColData"
  )

  meta[[COL_TIME]] <- factor(
    as.character(meta[[COL_TIME]]),
    levels = TIME_LEVELS
  )

  meta
}

# Build a donor x time x subtype pseudobulk matrix.
# mode = "sum" for counts and "mean" for chromVAR z-scores.
make_pseudobulk <- function(
    mat,
    meta,
    subtype_col,
    min_cells = 30,
    mode = c("sum", "mean")) {

  mode <- match.arg(mode)

  common_cells <- intersect(colnames(mat), rownames(meta))
  if (length(common_cells) == 0L) {
    stop("No shared cells between matrix and ArchR metadata.")
  }

  mat <- mat[, common_cells, drop = FALSE]
  meta_sub <- meta[common_cells, , drop = FALSE]

  valid <-
    !is.na(meta_sub[[COL_DONOR]]) &
    !is.na(meta_sub[[COL_TIME]]) &
    !is.na(meta_sub[[subtype_col]]) &
    as.character(meta_sub[[COL_DONOR]]) != "" &
    as.character(meta_sub[[subtype_col]]) != ""

  mat <- mat[, valid, drop = FALSE]
  meta_sub <- meta_sub[valid, , drop = FALSE]

  meta_sub$donor <- as.character(meta_sub[[COL_DONOR]])
  meta_sub$time <- as.character(meta_sub[[COL_TIME]])
  meta_sub$subtype <- as.character(meta_sub[[subtype_col]])

  meta_sub$pb_id <- paste(
    meta_sub$donor,
    meta_sub$time,
    meta_sub$subtype,
    sep = "__"
  )

  group_size <- table(meta_sub$pb_id)
  valid_groups <- names(group_size)[group_size >= min_cells]

  keep <- meta_sub$pb_id %in% valid_groups
  mat <- mat[, keep, drop = FALSE]
  meta_sub <- meta_sub[keep, , drop = FALSE]

  if (ncol(mat) == 0L) {
    stop(
      "No pseudobulk groups retained for ", subtype_col,
      ". Reduce min_cells or check metadata."
    )
  }

  group_levels <- unique(meta_sub$pb_id)

  cell_to_group <- Matrix::sparseMatrix(
    i = seq_len(nrow(meta_sub)),
    j = match(meta_sub$pb_id, group_levels),
    x = 1,
    dims = c(nrow(meta_sub), length(group_levels))
  )
  colnames(cell_to_group) <- group_levels

  pb_mat <- mat %*% cell_to_group

  if (mode == "mean") {
    group_n <- Matrix::colSums(cell_to_group)
    pb_mat <- pb_mat %*% Matrix::Diagonal(x = 1 / group_n)
  }

  colnames(pb_mat) <- group_levels

  pb_meta <- meta_sub %>%
    distinct(pb_id, donor, time, subtype)

  pb_meta <- pb_meta[
    match(colnames(pb_mat), pb_meta$pb_id),
    ,
    drop = FALSE
  ]

  pb_meta$time <- factor(pb_meta$time, levels = TIME_LEVELS)
  pb_meta$n_cells <- as.integer(group_size[pb_meta$pb_id])

  list(
    matrix = pb_mat,
    metadata = pb_meta,
    cell_metadata = meta_sub,
    cell_to_group = cell_to_group
  )
}

# Keep donors with both denominator and numerator samples.
get_paired_metadata <- function(
    pb_meta,
    subtype_value,
    denominator = BASELINE_TIME,
    numerator = TARGET_TIME) {

  pair_meta <- pb_meta %>%
    filter(
      subtype == subtype_value,
      as.character(time) %in% c(denominator, numerator)
    )

  paired_donors <- pair_meta %>%
    distinct(donor, time) %>%
    count(donor, name = "n_time") %>%
    filter(n_time == 2L) %>%
    pull(donor)

  pair_meta <- pair_meta %>%
    filter(donor %in% paired_donors)

  pair_meta$time <- factor(
    as.character(pair_meta$time),
    levels = c(denominator, numerator)
  )
  pair_meta$donor <- factor(pair_meta$donor)

  pair_meta
}

# Paired limma analysis for normalized continuous matrices:
# log2(TPM+1), chromVAR z-scores, or TCellSI scores.
run_paired_limma <- function(
    mat,
    pb_meta,
    denominator = BASELINE_TIME,
    numerator = TARGET_TIME) {

  subtype_values <- sort(unique(pb_meta$subtype))

  purrr::map_dfr(subtype_values, function(subtype_value) {

    pair_meta <- get_paired_metadata(
      pb_meta = pb_meta,
      subtype_value = subtype_value,
      denominator = denominator,
      numerator = numerator
    )

    if (n_distinct(pair_meta$donor) < 2L) {
      message("Skipping ", subtype_value, ": fewer than two paired donors.")
      return(NULL)
    }

    sample_ids <- pair_meta$pb_id
    sample_ids <- sample_ids[sample_ids %in% colnames(mat)]
    pair_meta <- pair_meta[match(sample_ids, pair_meta$pb_id), , drop = FALSE]

    y <- as.matrix(mat[, sample_ids, drop = FALSE])

    design <- model.matrix(~ donor + time, data = pair_meta)
    coef_name <- paste0("time", numerator)

    if (!coef_name %in% colnames(design)) {
      stop("Cannot find coefficient: ", coef_name)
    }

    fit <- limma::lmFit(y, design)
    fit <- limma::eBayes(
      fit,
      trend = TRUE,
      robust = nrow(y) >= 50L
    )

    limma::topTable(
      fit,
      coef = coef_name,
      number = Inf,
      sort.by = "P"
    ) %>%
      rownames_to_column("Feature") %>%
      mutate(
        subtype = subtype_value,
        contrast = paste0(numerator, "_vs_", denominator),
        n_donor = n_distinct(pair_meta$donor)
      )
  })
}

# Paired edgeR quasi-likelihood analysis for pseudobulk count matrices.
run_paired_edger <- function(
    count_mat,
    pb_meta,
    denominator = BASELINE_TIME,
    numerator = TARGET_TIME) {

  subtype_values <- sort(unique(pb_meta$subtype))

  purrr::map_dfr(subtype_values, function(subtype_value) {

    pair_meta <- get_paired_metadata(
      pb_meta = pb_meta,
      subtype_value = subtype_value,
      denominator = denominator,
      numerator = numerator
    )

    if (n_distinct(pair_meta$donor) < 2L) {
      message("Skipping peak DE for ", subtype_value,
              ": fewer than two paired donors.")
      return(NULL)
    }

    sample_ids <- pair_meta$pb_id
    sample_ids <- sample_ids[sample_ids %in% colnames(count_mat)]
    pair_meta <- pair_meta[match(sample_ids, pair_meta$pb_id), , drop = FALSE]

    y_counts <- as.matrix(count_mat[, sample_ids, drop = FALSE])
    if (any(!is.finite(y_counts) | y_counts < 0 | abs(y_counts-round(y_counts)) > 1e-8)) {
      stop("PeakMatrix must contain nonnegative raw counts; refusing to truncate values.")
    }

    design <- model.matrix(~ donor + time, data = pair_meta)
    coef_name <- paste0("time", numerator)

    y <- edgeR::DGEList(counts = y_counts)
    keep <- edgeR::filterByExpr(y, design = design, min.count = 10)

    if (sum(keep) < 10L) {
      message("Skipping peak DE for ", subtype_value,
              ": too few peaks passed filterByExpr.")
      return(NULL)
    }

    y <- y[keep, , keep.lib.sizes = FALSE]
    y <- edgeR::calcNormFactors(y, method = "TMM")
    y <- edgeR::estimateDisp(y, design, robust = TRUE)

    fit <- edgeR::glmQLFit(y, design, robust = TRUE)
    qlf <- edgeR::glmQLFTest(
      fit,
      coef = match(coef_name, colnames(design))
    )

    edgeR::topTags(qlf, n = Inf, sort.by = "PValue")$table %>%
      rownames_to_column("Feature") %>%
      mutate(
        subtype = subtype_value,
        contrast = paste0(numerator, "_vs_", denominator),
        n_donor = n_distinct(pair_meta$donor)
      )
  })
}

# Convert selected pseudobulk counts to log2(CPM+1) for visualization.
selected_counts_to_logcpm <- function(count_mat, selected_features) {
  selected_features <- intersect(selected_features, rownames(count_mat))
  selected_counts <- as.matrix(
    count_mat[selected_features, , drop = FALSE]
  )
  lib_size <- Matrix::colSums(count_mat)
  log2(sweep(selected_counts, 2, lib_size, "/") * 1e6 + 1)
}

# Generic dynamic heatmap and donor-line plots.
plot_dynamic_features <- function(
    mat,
    pb_meta,
    selected_table,
    label_col,
    output_prefix,
    baseline_time = BASELINE_TIME,
    max_line_features = 6) {

  if (nrow(selected_table) == 0L) {
    message("No selected features for ", output_prefix)
    return(invisible(NULL))
  }

  selected_table <- selected_table %>%
    distinct(subtype, Feature, .keep_all = TRUE)

  features <- intersect(unique(selected_table$Feature), rownames(mat))
  if (length(features) == 0L) {
    message("No selected features found in matrix for ", output_prefix)
    return(invisible(NULL))
  }

  selected_labels <- selected_table %>%
    transmute(
      subtype,
      Feature,
      Label = .data[[label_col]]
    )

  long_df <- as.data.frame(
    t(as.matrix(mat[features, , drop = FALSE]))
  ) %>%
    rownames_to_column("pb_id") %>%
    left_join(pb_meta, by = "pb_id") %>%
    pivot_longer(
      cols = all_of(features),
      names_to = "Feature",
      values_to = "value"
    ) %>%
    inner_join(
      selected_labels,
      by = c("subtype", "Feature")
    ) %>%
    group_by(donor, subtype, Feature, Label) %>%
    mutate(
      baseline_value = dplyr::first(
        value[as.character(time) == baseline_time],
        default = NA_real_
      ),
      delta_baseline = value - baseline_value
    ) %>%
    ungroup() %>%
    filter(is.finite(delta_baseline))

  write.csv(
    long_df,
    file.path(OUTDIR, paste0(output_prefix, "_donor_trajectories.csv")),
    row.names = FALSE
  )

  summary_df <- long_df %>%
    group_by(subtype, Feature, Label, time) %>%
    summarise(
      mean_delta = mean(delta_baseline, na.rm = TRUE),
      se_delta = sem(delta_baseline),
      .groups = "drop"
    ) %>%
    group_by(subtype, Feature) %>%
    mutate(dynamic_z = scale_rows_safely(mean_delta)) %>%
    ungroup()

  write.csv(
    summary_df,
    file.path(OUTDIR, paste0(output_prefix, "_summary.csv")),
    row.names = FALSE
  )

  p_heatmap <- ggplot(
    summary_df,
    aes(x = time, y = Label, fill = dynamic_z)
  ) +
    geom_tile() +
    facet_wrap(~ subtype, scales = "free_y", ncol = 2) +
    scale_fill_gradient2(
      low = "#3B6FB6",
      mid = "white",
      high = "#B64040",
      midpoint = 0,
      name = "Row-scaled\ndynamics"
    ) +
    labs(x = NULL, y = NULL) +
    theme_classic(base_size = 8) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 5.5),
      strip.background = element_blank(),
      panel.spacing = grid::unit(2.5, "mm")
    )

  ggsave(
    file.path(OUTDIR, paste0(output_prefix, "_heatmap.pdf")),
    p_heatmap,
    width = 9,
    height = max(5, 2.2 * length(unique(summary_df$subtype))),
    limitsize = FALSE
  )

  subtype_values <- unique(selected_table$subtype)

  purrr::walk(subtype_values, function(subtype_value) {

    selected_sub <- selected_table %>%
      filter(subtype == subtype_value) %>%
      slice_head(n = max_line_features)

    line_df <- long_df %>%
      filter(
        subtype == subtype_value,
        Feature %in% selected_sub$Feature
      ) %>%
      mutate(time_index = match(as.character(time), TIME_LEVELS))

    if (nrow(line_df) == 0L) return(invisible(NULL))

    mean_df <- line_df %>%
      group_by(Feature, Label, time, time_index) %>%
      summarise(
        mean_delta = mean(delta_baseline, na.rm = TRUE),
        se_delta = sem(delta_baseline),
        .groups = "drop"
      )

    p_line <- ggplot() +
      geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey70") +
      geom_vline(
        xintercept = 6.5,
        linetype = 2,
        linewidth = 0.35,
        colour = "grey45"
      ) +
      geom_line(
        data = line_df,
        aes(
          x = time_index,
          y = delta_baseline,
          group = donor
        ),
        linewidth = 0.35,
        alpha = 0.35,
        colour = "grey40"
      ) +
      geom_point(
        data = line_df,
        aes(x = time_index, y = delta_baseline),
        size = 0.8,
        alpha = 0.45,
        colour = "grey40"
      ) +
      geom_ribbon(
        data = mean_df,
        aes(
          x = time_index,
          ymin = mean_delta - se_delta,
          ymax = mean_delta + se_delta
        ),
        alpha = 0.16,
        fill = "#B64040"
      ) +
      geom_line(
        data = mean_df,
        aes(x = time_index, y = mean_delta),
        linewidth = 0.8,
        colour = "#B64040"
      ) +
      geom_point(
        data = mean_df,
        aes(x = time_index, y = mean_delta),
        size = 1.3,
        colour = "#B64040"
      ) +
      facet_wrap(~ Label, scales = "free_y", ncol = 3) +
      scale_x_continuous(
        breaks = seq_along(TIME_LEVELS),
        labels = TIME_LEVELS
      ) +
      labs(
        title = subtype_value,
        x = NULL,
        y = paste0("Change from ", baseline_time)
      ) +
      theme_classic(base_size = 8) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank(),
        plot.title = element_text(face = "bold", size = 9)
      )

    ggsave(
      file.path(
        OUTDIR,
        paste0(
          output_prefix,
          "_lines_",
          safe_filename(subtype_value),
          ".pdf"
        )
      ),
      p_line,
      width = 8.5,
      height = 5.5
    )
  })

  invisible(list(long = long_df, summary = summary_df, heatmap = p_heatmap))
}

###############################################################################
# 2. LOAD AND CHECK ARCHR METADATA/MATRICES
###############################################################################

meta_all <- get_project_metadata(PROJ)

available_matrices <- getAvailableMatrices(PROJ)
required_matrices <- c(
  "GeneExpressionMatrix",
  "PeakMatrix",
  "MotifMatrix"
)

missing_matrices <- setdiff(required_matrices, available_matrices)
if (length(missing_matrices) > 0L) {
  stop(
    "ArchRProject is missing matrices: ",
    paste(missing_matrices, collapse = ", ")
  )
}

write.csv(
  meta_all,
  file.path(OUTDIR, "ArchR_cell_metadata_used.csv"),
  row.names = FALSE
)

###############################################################################
# 3. RNA PSEUDOBULK AND TCellSI
###############################################################################

message("[1/8] Extracting GeneExpressionMatrix...")

se_rna <- getMatrixFromProject(
  ArchRProj = PROJ,
  useMatrix = "GeneExpressionMatrix"
)

rna_counts <- assay(se_rna)
rna_rowdata <- as.data.frame(rowData(se_rna))

symbol_col <- first_existing_column(
  rna_rowdata,
  c("name", "symbol", "gene_name", "Gene", "gene")
)

if (is.null(symbol_col)) {
  gene_symbol <- rownames(se_rna)
} else {
  gene_symbol <- as.character(rna_rowdata[[symbol_col]])
}

if (is.null(gene_symbol)) {
  stop("Cannot identify gene symbols in GeneExpressionMatrix rowData.")
}

gene_length <- tryCatch(
  as.numeric(width(rowRanges(se_rna))),
  error = function(e) rep(NA_real_, nrow(se_rna))
)

if (length(gene_length) != nrow(se_rna) || all(!is.finite(gene_length))) {
  start_col <- first_existing_column(rna_rowdata, c("start", "Start"))
  end_col <- first_existing_column(rna_rowdata, c("end", "End"))
  if (!is.null(start_col) && !is.null(end_col)) {
    gene_length <-
      as.numeric(rna_rowdata[[end_col]]) -
      as.numeric(rna_rowdata[[start_col]]) + 1
  }
}

valid_gene <-
  !is.na(gene_symbol) &
  gene_symbol != "" &
  is.finite(gene_length) &
  gene_length > 0

rna_counts <- rna_counts[valid_gene, , drop = FALSE]
gene_symbol <- toupper(gene_symbol[valid_gene])
gene_length <- gene_length[valid_gene]

if (mean(grepl("^ENSG", gene_symbol)) > 0.5) {
  stop(
    "Most RNA feature names are Ensembl IDs. Convert them to HGNC gene symbols ",
    "before running TCellSI."
  )
}

# Collapse duplicate gene symbols by summing counts.
symbol_levels <- unique(gene_symbol)
gene_to_symbol <- Matrix::sparseMatrix(
  i = match(gene_symbol, symbol_levels),
  j = seq_along(gene_symbol),
  x = 1,
  dims = c(length(symbol_levels), length(gene_symbol))
)

rna_symbol_counts <- gene_to_symbol %*% rna_counts
rownames(rna_symbol_counts) <- symbol_levels
colnames(rna_symbol_counts) <- colnames(rna_counts)

gene_length_symbol <- tapply(
  gene_length,
  factor(gene_symbol, levels = symbol_levels),
  median,
  na.rm = TRUE
)
gene_length_symbol <- as.numeric(
  gene_length_symbol[rownames(rna_symbol_counts)]
)
names(gene_length_symbol) <- rownames(rna_symbol_counts)

counts_to_log2tpm <- function(pb_counts, gene_lengths_bp) {
  gene_lengths_bp <- gene_lengths_bp[rownames(pb_counts)]
  gene_lengths_kb <- gene_lengths_bp / 1000

  valid <- is.finite(gene_lengths_kb) & gene_lengths_kb > 0
  pb_counts <- pb_counts[valid, , drop = FALSE]
  gene_lengths_kb <- gene_lengths_kb[valid]

  rpk <- Matrix::Diagonal(x = 1 / gene_lengths_kb) %*% pb_counts
  rpk_sum <- Matrix::colSums(rpk)

  if (any(!is.finite(rpk_sum) | rpk_sum <= 0)) {
    stop("One or more pseudobulk samples have zero/invalid total RPK.")
  }

  tpm <- rpk %*% Matrix::Diagonal(x = 1e6 / rpk_sum)
  result <- log2(as.matrix(tpm) + 1)
  dimnames(result) <- dimnames(pb_counts)
  result
}

message("[2/8] Creating T-cell and B-cell RNA pseudobulk matrices...")

T_rna_pb <- make_pseudobulk(
  mat = rna_symbol_counts,
  meta = meta_all,
  subtype_col = COL_T_SUBTYPE,
  min_cells = MIN_CELLS_RNA,
  mode = "sum"
)

B_rna_pb <- make_pseudobulk(
  mat = rna_symbol_counts,
  meta = meta_all,
  subtype_col = COL_B_SUBTYPE,
  min_cells = MIN_CELLS_RNA,
  mode = "sum"
)

T_log2tpm <- counts_to_log2tpm(
  T_rna_pb$matrix,
  gene_length_symbol
)

B_log2tpm <- counts_to_log2tpm(
  B_rna_pb$matrix,
  gene_length_symbol
)

saveRDS(
  list(
    T_log2tpm = T_log2tpm,
    T_pb_counts = T_rna_pb$matrix,
    T_pb_meta = T_rna_pb$metadata,
    B_log2tpm = B_log2tpm,
    B_pb_counts = B_rna_pb$matrix,
    B_pb_meta = B_rna_pb$metadata
  ),
  file.path(OUTDIR, "RNA_pseudobulk_TPM_objects.rds")
)

message("[3/8] Running TCellSI...")

TCellSI_scores <- TCellSI::TCSS_Calculate(
  object = as.data.frame(T_log2tpm),
  ref = TRUE
)

TCellSI_scores <- as.matrix(TCellSI_scores)

if (!all(colnames(TCellSI_scores) %in% T_rna_pb$metadata$pb_id)) {
  stop("TCellSI output sample names do not match T-cell pseudobulk IDs.")
}

write.csv(
  TCellSI_scores,
  file.path(OUTDIR, "TCellSI_scores_matrix.csv")
)

TCellSI_long <- as.data.frame(t(TCellSI_scores)) %>%
  rownames_to_column("pb_id") %>%
  left_join(T_rna_pb$metadata, by = "pb_id") %>%
  pivot_longer(
    cols = all_of(rownames(TCellSI_scores)),
    names_to = "State",
    values_to = "TCSS"
  ) %>%
  group_by(donor, subtype, State) %>%
  mutate(
    day30_score = dplyr::first(
      TCSS[as.character(time) == BASELINE_TIME],
      default = NA_real_
    ),
    Delta_day30 = TCSS - day30_score
  ) %>%
  ungroup()

write.csv(
  TCellSI_long,
  file.path(OUTDIR, "TCellSI_scores_long.csv"),
  row.names = FALSE
)

TCellSI_DE <- run_paired_limma(
  mat = TCellSI_scores,
  pb_meta = T_rna_pb$metadata
) %>%
  rename(State = Feature)

write.csv(
  TCellSI_DE,
  file.path(OUTDIR, "TCellSI_day44_vs_day30_statistics.csv"),
  row.names = FALSE
)

TCellSI_summary <- TCellSI_long %>%
  filter(is.finite(Delta_day30)) %>%
  group_by(subtype, State, time) %>%
  summarise(
    mean_TCSS = mean(TCSS, na.rm = TRUE),
    se_TCSS = sem(TCSS),
    mean_delta = mean(Delta_day30, na.rm = TRUE),
    se_delta = sem(Delta_day30),
    .groups = "drop"
  )

p_TCellSI_heatmap <- ggplot(
  TCellSI_summary,
  aes(x = time, y = State, fill = mean_delta)
) +
  geom_tile() +
  facet_wrap(~ subtype, ncol = 2) +
  scale_fill_gradient2(
    low = "#3B6FB6",
    mid = "white",
    high = "#B64040",
    midpoint = 0,
    name = paste0("TCSS change\nfrom ", BASELINE_TIME)
  ) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 8) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_blank(),
    panel.spacing = grid::unit(2.5, "mm")
  )

ggsave(
  file.path(OUTDIR, "TCellSI_dynamic_heatmap.pdf"),
  p_TCellSI_heatmap,
  width = 9,
  height = max(5, 2.1 * length(unique(TCellSI_summary$subtype))),
  limitsize = FALSE
)

purrr::walk(unique(TCellSI_long$subtype), function(subtype_value) {

  plot_df <- TCellSI_long %>%
    filter(subtype == subtype_value, is.finite(Delta_day30)) %>%
    mutate(time_index = match(as.character(time), TIME_LEVELS))

  mean_df <- plot_df %>%
    group_by(State, time, time_index) %>%
    summarise(
      mean_delta = mean(Delta_day30, na.rm = TRUE),
      se_delta = sem(Delta_day30),
      .groups = "drop"
    )

  p <- ggplot() +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey70") +
    geom_vline(
      xintercept = 6.5,
      linetype = 2,
      linewidth = 0.35,
      colour = "grey45"
    ) +
    geom_line(
      data = plot_df,
      aes(x = time_index, y = Delta_day30, group = donor),
      linewidth = 0.35,
      alpha = 0.35,
      colour = "grey40"
    ) +
    geom_point(
      data = plot_df,
      aes(x = time_index, y = Delta_day30),
      size = 0.8,
      alpha = 0.45,
      colour = "grey40"
    ) +
    geom_ribbon(
      data = mean_df,
      aes(
        x = time_index,
        ymin = mean_delta - se_delta,
        ymax = mean_delta + se_delta
      ),
      fill = "#B64040",
      alpha = 0.16
    ) +
    geom_line(
      data = mean_df,
      aes(x = time_index, y = mean_delta),
      linewidth = 0.8,
      colour = "#B64040"
    ) +
    geom_point(
      data = mean_df,
      aes(x = time_index, y = mean_delta),
      size = 1.2,
      colour = "#B64040"
    ) +
    facet_wrap(~ State, scales = "free_y", ncol = 4) +
    scale_x_continuous(
      breaks = seq_along(TIME_LEVELS),
      labels = TIME_LEVELS
    ) +
    labs(
      title = subtype_value,
      x = NULL,
      y = paste0("TCSS change from ", BASELINE_TIME)
    ) +
    theme_classic(base_size = 8) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.background = element_blank(),
      plot.title = element_text(face = "bold", size = 9)
    )

  ggsave(
    file.path(
      OUTDIR,
      paste0("TCellSI_lines_", safe_filename(subtype_value), ".pdf")
    ),
    p,
    width = 10,
    height = 6
  )
})

###############################################################################
# 4. TF RNA DIFFERENTIAL EXPRESSION: day44 vs day30
###############################################################################

message("[4/8] Testing TF RNA expression at day44 vs day30...")

T_RNA_DE <- run_paired_limma(
  mat = T_log2tpm,
  pb_meta = T_rna_pb$metadata
) %>%
  rename(Gene = Feature)

B_RNA_DE <- run_paired_limma(
  mat = B_log2tpm,
  pb_meta = B_rna_pb$metadata
) %>%
  rename(Gene = Feature)

write.csv(
  T_RNA_DE,
  file.path(OUTDIR, "Tcell_RNA_day44_vs_day30_all_genes.csv"),
  row.names = FALSE
)

write.csv(
  B_RNA_DE,
  file.path(OUTDIR, "Bcell_RNA_day44_vs_day30_all_genes.csv"),
  row.names = FALSE
)

###############################################################################
# 5. MOTIF ACTIVITY DIFFERENTIAL ANALYSIS AND TF DYNAMICS
###############################################################################

message("[5/8] Extracting chromVAR motif z-scores...")

se_motif_z <- getMatrixFromProject(
  ArchRProj = PROJ,
  useMatrix = "MotifMatrix",
  useSeqnames = "z"
)

motif_z <- assay(se_motif_z)
motif_rowdata <- as.data.frame(rowData(se_motif_z))

motif_name_col <- first_existing_column(
  motif_rowdata,
  c("name", "Name", "motif", "Motif")
)

if (is.null(motif_name_col)) {
  motif_name <- rownames(se_motif_z)
} else {
  motif_name <- as.character(motif_rowdata[[motif_name_col]])
}

if (is.null(motif_name)) {
  stop("Cannot identify motif names from MotifMatrix rowData.")
}

motif_feature <- make.unique(motif_name)
rownames(motif_z) <- motif_feature

motif_map <- data.frame(
  Feature = motif_feature,
  Motif = motif_name,
  TF = toupper(sub("_[0-9]+$", "", motif_name)),
  stringsAsFactors = FALSE
)

T_motif_pb <- make_pseudobulk(
  mat = motif_z,
  meta = meta_all,
  subtype_col = COL_T_SUBTYPE,
  min_cells = MIN_CELLS_MOTIF,
  mode = "mean"
)

B_motif_pb <- make_pseudobulk(
  mat = motif_z,
  meta = meta_all,
  subtype_col = COL_B_SUBTYPE,
  min_cells = MIN_CELLS_MOTIF,
  mode = "mean"
)

T_motif_DE <- run_paired_limma(
  mat = T_motif_pb$matrix,
  pb_meta = T_motif_pb$metadata
) %>%
  left_join(motif_map, by = "Feature")

B_motif_DE <- run_paired_limma(
  mat = B_motif_pb$matrix,
  pb_meta = B_motif_pb$metadata
) %>%
  left_join(motif_map, by = "Feature")

write.csv(
  T_motif_DE,
  file.path(OUTDIR, "Tcell_motif_day44_vs_day30_all.csv"),
  row.names = FALSE
)

write.csv(
  B_motif_DE,
  file.path(OUTDIR, "Bcell_motif_day44_vs_day30_all.csv"),
  row.names = FALSE
)

# Merge motif activity with TF RNA expression.
prepare_tf_candidates <- function(motif_de, rna_de) {

  rna_support <- rna_de %>%
    transmute(
      subtype,
      TF = toupper(Gene),
      RNA_logFC = logFC,
      RNA_AveExpr = AveExpr,
      RNA_PValue = P.Value,
      RNA_FDR = adj.P.Val
    )

  motif_de %>%
    filter(
      adj.P.Val <= TF_FDR,
      logFC >= TF_MIN_DELTA_Z
    ) %>%
    mutate(
      TF = toupper(TF),
      single_TF_motif = !grepl("::|,|/", TF)
    ) %>%
    left_join(
      rna_support,
      by = c("subtype", "TF")
    ) %>%
    mutate(
      RNA_up = !is.na(RNA_logFC) & RNA_logFC > 0,
      RNA_significant_up =
        !is.na(RNA_FDR) &
        RNA_FDR <= RNA_TF_FDR &
        RNA_logFC > 0,
      Evidence = case_when(
        RNA_significant_up ~ "Motif+RNA significant",
        RNA_up ~ "Motif+RNA concordant",
        TRUE ~ "Motif only"
      )
    ) %>%
    arrange(
      subtype,
      desc(RNA_significant_up),
      adj.P.Val,
      desc(logFC)
    ) %>%
    group_by(subtype, TF) %>%
    slice_head(n = 1) %>%
    ungroup()
}

T_TF_candidates <- prepare_tf_candidates(T_motif_DE, T_RNA_DE)
B_TF_candidates <- prepare_tf_candidates(B_motif_DE, B_RNA_DE)

write.csv(
  T_TF_candidates,
  file.path(OUTDIR, "Tcell_day44_up_TF_candidates.csv"),
  row.names = FALSE
)

write.csv(
  B_TF_candidates,
  file.path(OUTDIR, "Bcell_day44_up_TF_candidates.csv"),
  row.names = FALSE
)

T_selected_TF <- T_TF_candidates %>%
  group_by(subtype) %>%
  slice_head(n = TOP_TF_PER_SUBTYPE) %>%
  ungroup() %>%
  mutate(Label = TF)

B_selected_TF <- B_TF_candidates %>%
  group_by(subtype) %>%
  slice_head(n = TOP_TF_PER_SUBTYPE) %>%
  ungroup() %>%
  mutate(Label = TF)

T_TF_dynamic <- plot_dynamic_features(
  mat = T_motif_pb$matrix,
  pb_meta = T_motif_pb$metadata,
  selected_table = T_selected_TF,
  label_col = "Label",
  output_prefix = "Tcell_day44_up_TF_dynamic",
  baseline_time = BASELINE_TIME,
  max_line_features = min(6, TOP_TF_PER_SUBTYPE)
)

B_TF_dynamic <- plot_dynamic_features(
  mat = B_motif_pb$matrix,
  pb_meta = B_motif_pb$metadata,
  selected_table = B_selected_TF,
  label_col = "Label",
  output_prefix = "Bcell_day44_up_TF_dynamic",
  baseline_time = BASELINE_TIME,
  max_line_features = min(6, TOP_TF_PER_SUBTYPE)
)

# TF differential bubble plot.
plot_tf_bubble <- function(candidate_table, output_name) {
  if (nrow(candidate_table) == 0L) return(invisible(NULL))

  plot_df <- candidate_table %>%
    group_by(subtype) %>%
    slice_head(n = TOP_TF_PER_SUBTYPE) %>%
    ungroup() %>%
    mutate(
      minus_log10_fdr = -log10(pmax(adj.P.Val, 1e-300)),
      Evidence = factor(
        Evidence,
        levels = c(
          "Motif only",
          "Motif+RNA concordant",
          "Motif+RNA significant"
        )
      )
    )

  p <- ggplot(
    plot_df,
    aes(
      x = subtype,
      y = TF,
      size = minus_log10_fdr,
      colour = logFC,
      shape = Evidence
    )
  ) +
    geom_point(alpha = 0.9) +
    scale_colour_gradient2(
      low = "#3B6FB6",
      mid = "white",
      high = "#B64040",
      midpoint = 0,
      name = "Motif Δz"
    ) +
    scale_size_continuous(name = "-log10 FDR", range = c(1.5, 5)) +
    labs(x = NULL, y = NULL) +
    theme_classic(base_size = 8) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right"
    )

  ggsave(
    file.path(OUTDIR, output_name),
    p,
    width = 8,
    height = max(4.5, 0.25 * length(unique(plot_df$TF)) + 2)
  )

  invisible(p)
}

plot_tf_bubble(
  T_TF_candidates,
  "Tcell_day44_up_TF_bubble.pdf"
)

plot_tf_bubble(
  B_TF_candidates,
  "Bcell_day44_up_TF_bubble.pdf"
)

###############################################################################
# 6. PEAK DIFFERENTIAL ANALYSIS: B CELLS AND T CELLS
###############################################################################

message("[6/8] Extracting PeakMatrix and testing day44 vs day30 peaks...")

se_peak <- getMatrixFromProject(
  ArchRProj = PROJ,
  useMatrix = "PeakMatrix"
)

peak_counts <- assay(se_peak)
peak_rowdata <- as.data.frame(rowData(se_peak))

peak_chr <- if ("seqnames" %in% colnames(peak_rowdata)) {
  as.character(peak_rowdata$seqnames)
} else {
  as.character(seqnames(rowRanges(se_peak)))
}

peak_start <- if ("start" %in% colnames(peak_rowdata)) {
  as.integer(peak_rowdata$start)
} else {
  start(rowRanges(se_peak))
}

peak_end <- if ("end" %in% colnames(peak_rowdata)) {
  as.integer(peak_rowdata$end)
} else {
  end(rowRanges(se_peak))
}

peak_label <- paste0(peak_chr, ":", peak_start, "-", peak_end)
peak_feature <- make.unique(peak_label)
rownames(peak_counts) <- peak_feature

peak_map <- data.frame(
  Feature = peak_feature,
  Peak = peak_label,
  chr = peak_chr,
  start = peak_start,
  end = peak_end,
  stringsAsFactors = FALSE
)

T_peak_pb <- make_pseudobulk(
  mat = peak_counts,
  meta = meta_all,
  subtype_col = COL_T_SUBTYPE,
  min_cells = MIN_CELLS_PEAK,
  mode = "sum"
)

B_peak_pb <- make_pseudobulk(
  mat = peak_counts,
  meta = meta_all,
  subtype_col = COL_B_SUBTYPE,
  min_cells = MIN_CELLS_PEAK,
  mode = "sum"
)

T_peak_DE <- run_paired_edger(
  count_mat = T_peak_pb$matrix,
  pb_meta = T_peak_pb$metadata
) %>%
  left_join(peak_map, by = "Feature")

B_peak_DE <- run_paired_edger(
  count_mat = B_peak_pb$matrix,
  pb_meta = B_peak_pb$metadata
) %>%
  left_join(peak_map, by = "Feature")

write.csv(
  T_peak_DE,
  file.path(OUTDIR, "Tcell_peaks_day44_vs_day30_all.csv"),
  row.names = FALSE
)

write.csv(
  B_peak_DE,
  file.path(OUTDIR, "Bcell_peaks_day44_vs_day30_all.csv"),
  row.names = FALSE
)

select_dynamic_peaks <- function(peak_de) {
  peak_de %>%
    filter(
      FDR <= PEAK_FDR,
      abs(logFC) >= PEAK_MIN_LOGFC
    ) %>%
    mutate(Direction = if_else(logFC > 0, "Up", "Down")) %>%
    group_by(subtype, Direction) %>%
    arrange(FDR, desc(abs(logFC)), .by_group = TRUE) %>%
    slice_head(n = TOP_PEAKS_PER_DIRECTION) %>%
    ungroup() %>%
    mutate(Label = Peak)
}

T_selected_peaks <- select_dynamic_peaks(T_peak_DE)
B_selected_peaks <- select_dynamic_peaks(B_peak_DE)

write.csv(
  T_selected_peaks,
  file.path(OUTDIR, "Tcell_selected_day44_significant_peaks.csv"),
  row.names = FALSE
)

write.csv(
  B_selected_peaks,
  file.path(OUTDIR, "Bcell_selected_day44_significant_peaks.csv"),
  row.names = FALSE
)

# Compute selected-feature logCPM only, avoiding a dense all-peak matrix.
if (nrow(T_selected_peaks) > 0L) {
  T_peak_logcpm <- selected_counts_to_logcpm(
    T_peak_pb$matrix,
    unique(T_selected_peaks$Feature)
  )

  T_peak_dynamic <- plot_dynamic_features(
    mat = T_peak_logcpm,
    pb_meta = T_peak_pb$metadata,
    selected_table = T_selected_peaks,
    label_col = "Label",
    output_prefix = "Tcell_day44_significant_peak_dynamic",
    baseline_time = BASELINE_TIME,
    max_line_features = TOP_PEAK_LINES_PER_SUBTYPE
  )
}

if (nrow(B_selected_peaks) > 0L) {
  B_peak_logcpm <- selected_counts_to_logcpm(
    B_peak_pb$matrix,
    unique(B_selected_peaks$Feature)
  )

  B_peak_dynamic <- plot_dynamic_features(
    mat = B_peak_logcpm,
    pb_meta = B_peak_pb$metadata,
    selected_table = B_selected_peaks,
    label_col = "Label",
    output_prefix = "Bcell_day44_significant_peak_dynamic",
    baseline_time = BASELINE_TIME,
    max_line_features = TOP_PEAK_LINES_PER_SUBTYPE
  )
}

###############################################################################
# 7. B-CELL TF FOOTPRINTING
###############################################################################

message("[7/8] Preparing B-cell TF footprints...")

run_bcell_footprints <- function(
    proj,
    meta,
    candidate_table,
    footprint_times = FOOTPRINT_TIMES) {

  if (nrow(candidate_table) == 0L) {
    message("No B-cell TF candidates available for footprinting.")
    return(invisible(NULL))
  }

  motif_positions <- getPositions(
    ArchRProj = proj,
    name = "Motif"
  )

  available_motifs <- names(motif_positions)
  if (is.null(available_motifs)) {
    stop("Motif positions have no names. Check ArchR motif annotations.")
  }

  candidate_use <- candidate_table %>%
    filter(Motif %in% available_motifs) %>%
    group_by(subtype) %>%
    slice_head(n = TOP_FOOTPRINT_MOTIFS_PER_BTYPE) %>%
    ungroup()

  if (nrow(candidate_use) == 0L) {
    message(
      "No candidate motif names matched getPositions(PROJ, name='Motif')."
    )
    return(invisible(NULL))
  }

  write.csv(
    candidate_use,
    file.path(OUTDIR, "Bcell_footprint_motifs_selected.csv"),
    row.names = FALSE
  )

  purrr::walk(unique(candidate_use$subtype), function(subtype_value) {

    motifs_now <- candidate_use %>%
      filter(subtype == subtype_value) %>%
      pull(Motif) %>%
      unique()

    cells_now <- meta %>%
      filter(
        as.character(.data[[COL_B_SUBTYPE]]) == subtype_value,
        as.character(.data[[COL_TIME]]) %in% footprint_times
      ) %>%
      pull(cell)

    if (length(cells_now) < MIN_FOOTPRINT_CELLS) {
      message("Skipping footprint for ", subtype_value, ": too few cells.")
      return(invisible(NULL))
    }

    keep_flag <- getCellNames(proj) %in% cells_now
    proj_sub <- proj[keep_flag, ]

    sub_meta <- as.data.frame(getCellColData(proj_sub))
    proj_sub$FootprintTime <- factor(
      as.character(sub_meta[[COL_TIME]]),
      levels = footprint_times
    )

    se_foot <- getFootprints(
      ArchRProj = proj_sub,
      positions = motif_positions[motifs_now],
      groupBy = "FootprintTime",
      minCells = MIN_FOOTPRINT_CELLS
    )

    plotFootprints(
      seFoot = se_foot,
      names = motifs_now,
      flank = 250,
      flankNorm = 50,
      normMethod = "Subtract",
      smoothWindow = 5,
      baseSize = 8,
      ArchRProj = proj_sub,
      plotName = paste0(
        "Bcell_Footprints_",
        safe_filename(subtype_value)
      ),
      width = 5.5,
      height = 4.5,
      addDOC = FALSE,
      force = TRUE
    )
  })

  invisible(candidate_use)
}

if (RUN_B_FOOTPRINTS) {
  B_footprint_candidates <- run_bcell_footprints(
    proj = PROJ,
    meta = meta_all,
    candidate_table = B_TF_candidates
  )
}

###############################################################################
# 8. SAVE FINAL OBJECTS AND SESSION INFORMATION
###############################################################################

message("[8/8] Saving final R objects...")

saveRDS(
  list(
    TCellSI_scores = TCellSI_scores,
    TCellSI_long = TCellSI_long,
    TCellSI_DE = TCellSI_DE,
    T_RNA_DE = T_RNA_DE,
    B_RNA_DE = B_RNA_DE,
    T_motif_DE = T_motif_DE,
    B_motif_DE = B_motif_DE,
    T_TF_candidates = T_TF_candidates,
    B_TF_candidates = B_TF_candidates,
    T_peak_DE = T_peak_DE,
    B_peak_DE = B_peak_DE,
    T_selected_peaks = T_selected_peaks,
    B_selected_peaks = B_selected_peaks
  ),
  file.path(OUTDIR, "T_B_dynamic_analysis_final_objects.rds")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(OUTDIR, "sessionInfo.txt")
)

message("Pipeline completed. Results are in: ", normalizePath(OUTDIR))

writeLines(capture.output(sessionInfo()), file.path(OUTDIR, "sessionInfo_updated.txt"))
