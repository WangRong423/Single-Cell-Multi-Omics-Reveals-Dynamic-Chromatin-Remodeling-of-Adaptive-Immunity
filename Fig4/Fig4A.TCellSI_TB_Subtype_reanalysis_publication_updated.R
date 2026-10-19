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
## Re-analysis of TCellSI using PROJ$TB_Subtype
##
## ArchR object:
##   proj_TB_fig3
##
## Cell metadata:
##   Donor
##   Timepoint
##   TB_Subtype
##
## T-cell subtypes retained:
##   CD4_NaiveT, CD4_Naive_CCR4, CD8_Memory_GZMK, Treg, CD4_Memory_ZEB2, CD4_Memory,
##   CD8_NaiveT, CD8_Memory, gdT, CD4_CD8_DP
##
## Analysis unit:
##   Donor × Timepoint × TB_Subtype RNA pseudobulk
##
## Main outputs:
##   1. TCellSI scores for the unified 8 functional states
##   2. Mean ± SE longitudinal line plot
##   3. Donor trajectories plus subtype mean trajectories
##   4. CD4-oriented and CD8-oriented line plots
## ========================================================================

suppressPackageStartupMessages({
  library(ArchR)
  library(SummarizedExperiment)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(edgeR)
  library(TCellSI)
})

## ========================================================================
## 0. Configuration
## ========================================================================

PROJ <- select_tb_project()

COL_DONOR <- "Donor"
COL_TIME <- "Timepoint"
COL_SUBTYPE <- "TB_Subtype"

TIME_LEVELS <- c(
  "day0", "day1", "day3", "day6", "day14",
  "day30", "day31", "day33", "day36", "day44"
)

TIME_LABELS <- c(
  "D0", "D1", "D3", "D6", "D14",
  "D30", "D31", "D33", "D36", "D44"
)

## Explicitly retain only T-cell subtypes from TB_Subtype
T_SUBTYPES <- c(
  "CD4_NaiveT",
  "CD4_Naive_CCR4",
  "CD8_Memory_GZMK",
  "Treg",
  "CD4_Memory_ZEB2",
  "CD4_Memory",
  "CD8_NaiveT",
  "CD8_Memory",
  "gdT",
  "CD4_CD8_DP"
)

## Order used in the figure legend
T_SUBTYPE_ORDER <- c(
  "CD4_NaiveT",
  "CD4_Naive_CCR4",
  "CD8_Memory_GZMK",
  "Treg",
  "CD4_Memory_ZEB2",
  "CD4_Memory",
  "CD8_NaiveT",
  "CD8_Memory",
  "gdT",
  "CD4_CD8_DP"
)

## Unified TCellSI state order
STATE_ORDER <- c(
  "Quiescence",
  "Regulating",
  "Proliferation",
  "Helper",
  "Cytotoxicity",
  "Progenitor exhaustion",
  "Terminal exhaustion",
  "Senescence"
)

## Minimum cells in each Donor × Timepoint × TB_Subtype pseudobulk
MIN_CELLS_PER_PB <- 30

## Use TPM as the primary TCellSI input when gene length is available.
## A log2CPM sensitivity result is also calculated.
PRIMARY_METHOD <- "log2TPM"

OUTDIR <- file.path(UPDATED_BASE, "TCellSI_TB_Subtype_reanalysis_updated_labels")

for (subdir in c(
  "01_QC",
  "02_Pseudobulk",
  "03_TCellSI",
  "04_Figures"
)) {
  dir.create(
    file.path(OUTDIR, subdir),
    recursive = TRUE,
    showWarnings = FALSE
  )
}

## Publication palette: one fixed color per T-cell subtype
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

save_pdf <- function(
    plot_object,
    filename,
    width,
    height
) {
  dir.create(
    dirname(filename),
    recursive = TRUE,
    showWarnings = FALSE
  )

  ggplot2::ggsave(
    filename = filename,
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    device = grDevices::cairo_pdf
  )
}

## ========================================================================
## 1. Check metadata and retain TB_Subtype-defined T cells
## ========================================================================

meta_all <- as.data.frame(
  ArchR::getCellColData(PROJ)
)

meta_all$cell <- rownames(meta_all)

required_columns <- c(
  COL_DONOR,
  COL_TIME,
  COL_SUBTYPE
)

missing_columns <- setdiff(
  required_columns,
  colnames(meta_all)
)

if (length(missing_columns) > 0) {
  stop(
    "The ArchR metadata is missing: ",
    paste(missing_columns, collapse = ", ")
  )
}

meta_all[[COL_SUBTYPE]] <- rename_tb_labels(meta_all[[COL_SUBTYPE]],
  file.path(OUTDIR, "01_QC", "Subtype_label_mapping.csv"))

## Save the complete TB_Subtype count table
tb_subtype_counts <- as.data.frame(
  table(
    as.character(meta_all[[COL_SUBTYPE]]),
    useNA = "ifany"
  )
)

colnames(tb_subtype_counts) <- c(
  "TB_Subtype",
  "n_cells"
)

write.csv(
  tb_subtype_counts,
  file.path(
    OUTDIR,
    "01_QC",
    "TB_Subtype_all_cell_counts.csv"
  ),
  row.names = FALSE
)

meta_T <- meta_all %>%
  transmute(
    cell = cell,
    donor = as.character(.data[[COL_DONOR]]),
    timepoint = as.character(.data[[COL_TIME]]),
    subtype = as.character(.data[[COL_SUBTYPE]])
  ) %>%
  filter(
    !is.na(cell),
    !is.na(donor),
    !is.na(timepoint),
    !is.na(subtype),
    subtype %in% T_SUBTYPES,
    timepoint %in% TIME_LEVELS
  ) %>%
  mutate(
    timepoint = factor(
      timepoint,
      levels = TIME_LEVELS
    ),
    subtype = factor(
      subtype,
      levels = T_SUBTYPE_ORDER
    )
  )

if (nrow(meta_T) == 0) {
  stop(
    "No T cells were retained using PROJ$TB_Subtype."
  )
}

message(
  "T cells retained from TB_Subtype: ",
  format(nrow(meta_T), big.mark = ",")
)

tcell_count_qc <- meta_T %>%
  count(
    donor,
    timepoint,
    subtype,
    name = "n_cells",
    .drop = FALSE
  ) %>%
  filter(
    n_cells > 0
  )

write.csv(
  tcell_count_qc,
  file.path(
    OUTDIR,
    "01_QC",
    "Tcell_counts_by_donor_time_TB_Subtype.csv"
  ),
  row.names = FALSE
)

## ========================================================================
## 2. Extract GeneExpressionMatrix
## ========================================================================

if (
  !"GeneExpressionMatrix" %in%
    ArchR::getAvailableMatrices(PROJ)
) {
  stop(
    "GeneExpressionMatrix is not available in the ArchR project."
  )
}

se_rna <- ArchR::getMatrixFromProject(
  ArchRProj = PROJ,
  useMatrix = "GeneExpressionMatrix"
)

rna_cell_matrix <- SummarizedExperiment::assay(
  se_rna
)

rna_rowdata <- as.data.frame(
  SummarizedExperiment::rowData(
    se_rna
  )
)

## Align T-cell metadata with the RNA matrix
common_cells <- intersect(
  meta_T$cell,
  colnames(rna_cell_matrix)
)

if (length(common_cells) == 0) {
  stop(
    "No common cells between TB_Subtype T cells and GeneExpressionMatrix."
  )
}

rna_cell_matrix <- rna_cell_matrix[
  ,
  common_cells,
  drop = FALSE
]

meta_T <- meta_T[
  match(common_cells, meta_T$cell),
  ,
  drop = FALSE
]

stopifnot(
  identical(
    colnames(rna_cell_matrix),
    meta_T$cell
  )
)

## ========================================================================
## 3. Select the correct gene-symbol field
## ========================================================================

known_human_genes <- c(
  "BTG2", "KLF2", "FOXP1", "RUNX1", "MYC",
  "CD4", "CD8A", "CD8B", "TCF7", "IL7R",
  "CCR7", "LTB", "MALAT1", "NKG7", "GZMB",
  "PRF1", "IFNG", "TBX21", "FOXP3", "CTLA4",
  "IL2RA", "BATF", "JUNB", "FOS", "STAT5A"
)

candidate_columns <- intersect(
  c(
    "symbol",
    "gene_symbol",
    "gene_name",
    "gene",
    "name",
    "feature",
    "gene_id",
    "id"
  ),
  colnames(rna_rowdata)
)

gene_id_candidates <- list(
  rownames = rownames(se_rna)
)

for (candidate in candidate_columns) {
  gene_id_candidates[[paste0("rowData$", candidate)]] <- rna_rowdata[[candidate]]
}

gene_identifier_qc <- do.call(
  rbind,
  lapply(
    names(gene_id_candidates),
    function(source_name) {

      values <- toupper(
        trimws(
          as.character(
            gene_id_candidates[[source_name]]
          )
        )
      )

      data.frame(
        source = source_name,
        n_non_missing = sum(
          !is.na(values) &
            values != ""
        ),
        n_unique = length(
          unique(
            values[
              !is.na(values) &
                values != ""
            ]
          )
        ),
        known_gene_overlap = sum(
          known_human_genes %in% values
        ),
        CD4_present = "CD4" %in% values,
        CD8A_present = "CD8A" %in% values,
        stringsAsFactors = FALSE
      )
    }
  )
)

gene_identifier_qc <- gene_identifier_qc[
  order(
    -gene_identifier_qc$known_gene_overlap,
    -gene_identifier_qc$n_unique
  ),
  ,
  drop = FALSE
]

write.csv(
  gene_identifier_qc,
  file.path(
    OUTDIR,
    "01_QC",
    "GeneExpressionMatrix_gene_identifier_candidates.csv"
  ),
  row.names = FALSE
)

best_gene_source <- gene_identifier_qc$source[1]

if (
  gene_identifier_qc$known_gene_overlap[1] < 3
) {
  stop(
    "No reliable HGNC gene-symbol column was identified. ",
    "Inspect GeneExpressionMatrix_gene_identifier_candidates.csv."
  )
}

if (best_gene_source == "rownames") {

  gene_symbol_raw <- rownames(se_rna)

} else {

  best_column <- sub(
    "^rowData\\$",
    "",
    best_gene_source
  )

  gene_symbol_raw <- rna_rowdata[[best_column]]
}

gene_symbol_raw <- toupper(
  trimws(
    as.character(gene_symbol_raw)
  )
)

gene_symbol_raw[
  is.na(gene_symbol_raw) |
    gene_symbol_raw == "" |
    gene_symbol_raw == "NA"
] <- NA_character_

writeLines(
  c(
    paste0(
      "Selected gene-symbol source: ",
      best_gene_source
    ),
    paste0(
      "Known-gene overlap: ",
      gene_identifier_qc$known_gene_overlap[1]
    )
  ),
  con = file.path(
    OUTDIR,
    "01_QC",
    "Selected_gene_identifier.txt"
  )
)

message(
  "Selected gene-symbol source: ",
  best_gene_source
)

## ========================================================================
## 4. Collapse duplicate gene symbols
## ========================================================================

collapse_sparse_rows <- function(
    mat,
    feature_names
) {

  valid <- !is.na(feature_names) &
    feature_names != ""

  mat <- mat[
    valid,
    ,
    drop = FALSE
  ]

  feature_names <- feature_names[
    valid
  ]

  unique_features <- unique(
    feature_names
  )

  row_map <- Matrix::sparseMatrix(
    i = match(
      feature_names,
      unique_features
    ),
    j = seq_along(
      feature_names
    ),
    x = 1,
    dims = c(
      length(unique_features),
      length(feature_names)
    )
  )

  collapsed <- row_map %*% mat

  rownames(collapsed) <- unique_features
  colnames(collapsed) <- colnames(mat)

  collapsed
}

rna_gene_counts <- collapse_sparse_rows(
  mat = rna_cell_matrix,
  feature_names = gene_symbol_raw
)

## Check whether GeneExpressionMatrix is count-like
nonzero_values <- rna_gene_counts@x

integer_fraction <- mean(
  abs(
    nonzero_values -
      round(nonzero_values)
  ) < 1e-8
)

writeLines(
  paste0(
    "Fraction of nonzero matrix values that are integer-like: ",
    signif(integer_fraction, 5)
  ),
  con = file.path(
    OUTDIR,
    "01_QC",
    "GeneExpressionMatrix_count_like_check.txt"
  )
)

if (integer_fraction < 0.95) {
  warning(
    "GeneExpressionMatrix contains many non-integer values. ",
    "Confirm that it represents raw UMI counts before using pseudobulk."
  )
}

## ========================================================================
## 5. Build Donor × Timepoint × TB_Subtype RNA pseudobulk
## ========================================================================

pb_id <- paste(
  meta_T$donor,
  as.character(meta_T$timepoint),
  as.character(meta_T$subtype),
  sep = "__"
)

pb_ncells <- table(pb_id)

valid_pb_ids <- names(pb_ncells)[
  pb_ncells >= MIN_CELLS_PER_PB
]

keep_cells <- pb_id %in% valid_pb_ids

meta_T_pb <- meta_T[
  keep_cells,
  ,
  drop = FALSE
]

rna_gene_counts_pb_cells <- rna_gene_counts[
  ,
  keep_cells,
  drop = FALSE
]

pb_id <- pb_id[
  keep_cells
]

pb_levels <- unique(
  pb_id
)

cell_to_pb <- Matrix::sparseMatrix(
  i = seq_along(pb_id),
  j = match(pb_id, pb_levels),
  x = 1,
  dims = c(
    length(pb_id),
    length(pb_levels)
  )
)

colnames(cell_to_pb) <- pb_levels

T_counts <- rna_gene_counts_pb_cells %*%
  cell_to_pb

rownames(T_counts) <- rownames(
  rna_gene_counts_pb_cells
)

colnames(T_counts) <- pb_levels

T_pb_meta <- meta_T_pb %>%
  mutate(
    pb_id = pb_id
  ) %>%
  distinct(
    pb_id,
    donor,
    timepoint,
    subtype
  )

T_pb_meta <- T_pb_meta[
  match(
    colnames(T_counts),
    T_pb_meta$pb_id
  ),
  ,
  drop = FALSE
]

T_pb_meta$n_cells <- as.integer(
  pb_ncells[T_pb_meta$pb_id]
)

T_pb_meta$timepoint <- factor(
  as.character(T_pb_meta$timepoint),
  levels = TIME_LEVELS
)

T_pb_meta$subtype <- factor(
  as.character(T_pb_meta$subtype),
  levels = T_SUBTYPE_ORDER
)

write.csv(
  T_pb_meta,
  file.path(
    OUTDIR,
    "02_Pseudobulk",
    "Tcell_RNA_pseudobulk_metadata.csv"
  ),
  row.names = FALSE
)

saveRDS(
  list(
    counts = T_counts,
    metadata = T_pb_meta
  ),
  file.path(
    OUTDIR,
    "02_Pseudobulk",
    "Tcell_RNA_pseudobulk_counts.rds"
  )
)

## ========================================================================
## 6. Calculate log2(CPM + 1)
## ========================================================================

T_log2CPM <- edgeR::cpm(
  T_counts,
  log = TRUE,
  prior.count = 1
)

dimnames(T_log2CPM) <- dimnames(
  T_counts
)

## ========================================================================
## 7. Extract gene lengths and calculate log2(TPM + 1)
## ========================================================================

extract_gene_lengths <- function(
    se_object,
    rowdata
) {

  rr <- tryCatch(
    SummarizedExperiment::rowRanges(
      se_object
    ),
    error = function(e) NULL
  )

  if (
    inherits(rr, "GRanges") &&
      length(rr) == nrow(se_object)
  ) {

    return(
      list(
        length_bp = as.numeric(
          GenomicRanges::width(rr)
        ),
        source = "rowRanges"
      )
    )
  }

  if (
    all(
      c("start", "end") %in%
        colnames(rowdata)
    )
  ) {

    return(
      list(
        length_bp = abs(
          as.numeric(rowdata$end) -
            as.numeric(rowdata$start)
        ) + 1,
        source = "rowData start/end"
      )
    )
  }

  if ("width" %in% colnames(rowdata)) {

    return(
      list(
        length_bp = as.numeric(
          rowdata$width
        ),
        source = "rowData width"
      )
    )
  }

  list(
    length_bp = rep(
      NA_real_,
      nrow(se_object)
    ),
    source = "unavailable"
  )
}

gene_length_info <- extract_gene_lengths(
  se_object = se_rna,
  rowdata = rna_rowdata
)

gene_length_raw <- gene_length_info$length_bp

writeLines(
  paste0(
    "Gene length source: ",
    gene_length_info$source
  ),
  con = file.path(
    OUTDIR,
    "01_QC",
    "Gene_length_source.txt"
  )
)

valid_length <- !is.na(gene_symbol_raw) &
  gene_symbol_raw != "" &
  is.finite(gene_length_raw) &
  gene_length_raw > 0

length_table <- data.frame(
  gene = gene_symbol_raw[
    valid_length
  ],
  length_bp = gene_length_raw[
    valid_length
  ]
) %>%
  group_by(gene) %>%
  summarise(
    length_bp = median(
      length_bp,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

gene_lengths <- length_table$length_bp[
  match(
    rownames(T_counts),
    length_table$gene
  )
]

valid_gene_lengths <- is.finite(
  gene_lengths
) &
  gene_lengths > 0

T_log2TPM <- NULL

if (sum(valid_gene_lengths) > 0) {

  fallback_length <- median(
    gene_lengths[
      valid_gene_lengths
    ],
    na.rm = TRUE
  )

  gene_lengths[
    !valid_gene_lengths
  ] <- fallback_length

  gene_scaler <- Matrix::Diagonal(
    x = 1000 / gene_lengths
  )

  dimnames(gene_scaler) <- list(
    rownames(T_counts),
    rownames(T_counts)
  )

  rpk <- gene_scaler %*%
    T_counts

  dimnames(rpk) <- dimnames(
    T_counts
  )

  rpk_library_size <- Matrix::colSums(
    rpk
  )

  sample_scaler <- Matrix::Diagonal(
    x = 1e6 /
      pmax(
        rpk_library_size,
        1e-12
      )
  )

  dimnames(sample_scaler) <- list(
    colnames(T_counts),
    colnames(T_counts)
  )

  tpm <- rpk %*%
    sample_scaler

  dimnames(tpm) <- dimnames(
    T_counts
  )

  T_log2TPM <- log2(
    as.matrix(tpm) + 1
  )

  dimnames(T_log2TPM) <- dimnames(
    T_counts
  )
}

saveRDS(
  list(
    log2TPM = T_log2TPM,
    log2CPM = T_log2CPM,
    metadata = T_pb_meta
  ),
  file.path(
    OUTDIR,
    "02_Pseudobulk",
    "Tcell_TCellSI_input_matrices.rds"
  )
)

## ========================================================================
## 8. TCellSI preflight and scoring
## ========================================================================

normalize_state_names <- function(x) {

  x_original <- as.character(x)

  x_clean <- tolower(
    trimws(
      gsub(
        "[._-]+",
        " ",
        x_original
      )
    )
  )

  dplyr::case_when(
    grepl("^quies", x_clean) ~
      "Quiescence",

    grepl("^regulat", x_clean) ~
      "Regulating",

    grepl("^prolifer", x_clean) ~
      "Proliferation",

    grepl("^helper", x_clean) ~
      "Helper",

    grepl("^cytotox", x_clean) ~
      "Cytotoxicity",

    grepl("progenitor.*exhaust", x_clean) ~
      "Progenitor exhaustion",

    grepl("terminal.*exhaust", x_clean) ~
      "Terminal exhaustion",

    grepl("^senesc", x_clean) ~
      "Senescence",

    TRUE ~ x_original
  )
}

tcellsi_input_check <- function(
    expression_matrix,
    input_name
) {

  if (is.null(expression_matrix)) {
    return(NULL)
  }

  if (is.null(rownames(expression_matrix))) {
    stop(
      input_name,
      " has no gene row names."
    )
  }

  required_check <- c(
    "BTG2", "KLF2", "FOXP1",
    "RUNX1", "MYC", "CD4"
  )

  check_table <- data.frame(
    gene = required_check,
    present = required_check %in%
      rownames(expression_matrix)
  )

  write.csv(
    check_table,
    file.path(
      OUTDIR,
      "01_QC",
      paste0(
        "TCellSI_key_gene_check_",
        input_name,
        ".csv"
      )
    ),
    row.names = FALSE
  )

  if (
    sum(
      check_table$present[
        check_table$gene != "CD4"
      ]
    ) < 3
  ) {
    stop(
      input_name,
      " contains fewer than three Quiescence markers."
    )
  }

  if (!"CD4" %in% rownames(expression_matrix)) {
    stop(
      input_name,
      " does not contain CD4; ref = TRUE cannot be used."
    )
  }

  expression_matrix
}

TCellSI_scores <- list()

if (!is.null(T_log2TPM)) {

  T_log2TPM <- tcellsi_input_check(
    T_log2TPM,
    "log2TPM"
  )

  TCellSI_scores$log2TPM <-
    TCellSI::TCSS_Calculate(
      as.data.frame(
        T_log2TPM,
        check.names = FALSE
      ),
      ref = TRUE
    )
}

T_log2CPM <- tcellsi_input_check(
  T_log2CPM,
  "log2CPM"
)

TCellSI_scores$log2CPM <-
  TCellSI::TCSS_Calculate(
    as.data.frame(
      T_log2CPM,
      check.names = FALSE
    ),
    ref = TRUE
  )

## Standardize output orientation as state × pseudobulk sample
standardize_tcellsi_output <- function(
    score_object,
    sample_ids
) {

  score_matrix <- as.matrix(
    score_object
  )

  row_hits <- sum(
    normalize_state_names(
      rownames(score_matrix)
    ) %in% STATE_ORDER
  )

  col_hits <- sum(
    normalize_state_names(
      colnames(score_matrix)
    ) %in% STATE_ORDER
  )

  if (row_hits >= col_hits) {

    out <- score_matrix

  } else {

    out <- t(score_matrix)
  }

  rownames(out) <- normalize_state_names(
    rownames(out)
  )

  missing_states <- setdiff(
    STATE_ORDER,
    rownames(out)
  )

  if (length(missing_states) > 0) {
    stop(
      "TCellSI output is missing states: ",
      paste(missing_states, collapse = ", ")
    )
  }

  out <- out[
    STATE_ORDER,
    ,
    drop = FALSE
  ]

  missing_samples <- setdiff(
    sample_ids,
    colnames(out)
  )

  if (length(missing_samples) > 0) {
    stop(
      "TCellSI output is missing pseudobulk samples: ",
      paste(
        head(missing_samples, 10),
        collapse = ", "
      )
    )
  }

  out[
    ,
    sample_ids,
    drop = FALSE
  ]
}

for (method_name in names(TCellSI_scores)) {
  TCellSI_scores[[method_name]] <-
    standardize_tcellsi_output(
      score_object =
        TCellSI_scores[[method_name]],
      sample_ids = T_pb_meta$pb_id
    )
}

saveRDS(
  TCellSI_scores,
  file.path(
    OUTDIR,
    "03_TCellSI",
    "TCellSI_scores_TB_Subtype.rds"
  )
)

## Select primary result
if (
  PRIMARY_METHOD %in%
    names(TCellSI_scores)
) {

  primary_method_used <- PRIMARY_METHOD

} else {

  primary_method_used <- "log2CPM"

  warning(
    "log2TPM was unavailable; log2CPM is used as the primary result."
  )
}

TCellSI_primary <- TCellSI_scores[[primary_method_used]]

writeLines(
  paste0(
    "Primary TCellSI input method: ",
    primary_method_used
  ),
  con = file.path(
    OUTDIR,
    "03_TCellSI",
    "Primary_TCellSI_method.txt"
  )
)

## ========================================================================
## 9. Convert scores to long format
## ========================================================================

tcellsi_long <- as.data.frame(
  t(TCellSI_primary),
  check.names = FALSE
) %>%
  rownames_to_column(
    "pb_id"
  ) %>%
  left_join(
    T_pb_meta %>%
      transmute(
        pb_id = as.character(pb_id),
        donor = as.character(donor),
        timepoint = as.character(timepoint),
        subtype = as.character(subtype),
        n_cells = n_cells
      ),
    by = "pb_id"
  ) %>%
  pivot_longer(
    cols = all_of(
      STATE_ORDER
    ),
    names_to = "State",
    values_to = "Score"
  ) %>%
  filter(
    !is.na(donor),
    !is.na(timepoint),
    !is.na(subtype)
  ) %>%
  mutate(
    timepoint = factor(
      timepoint,
      levels = TIME_LEVELS
    ),
    time_index = match(
      as.character(timepoint),
      TIME_LEVELS
    ),
    subtype = factor(
      subtype,
      levels = T_SUBTYPE_ORDER
    ),
    State = factor(
      State,
      levels = STATE_ORDER
    )
  )

write.csv(
  tcellsi_long,
  file.path(
    OUTDIR,
    "03_TCellSI",
    "TCellSI_scores_long_TB_Subtype.csv"
  ),
  row.names = FALSE
)

tcellsi_summary <- tcellsi_long %>%
  group_by(
    subtype,
    timepoint,
    time_index,
    State
  ) %>%
  summarise(
    mean_score = mean(
      Score,
      na.rm = TRUE
    ),
    sd_score = stats::sd(
      Score,
      na.rm = TRUE
    ),
    n_donor = sum(
      is.finite(Score)
    ),
    se_score = ifelse(
      n_donor > 1,
      sd_score /
        sqrt(n_donor),
      0
    ),
    .groups = "drop"
  )

write.csv(
  tcellsi_summary,
  file.path(
    OUTDIR,
    "03_TCellSI",
    "TCellSI_scores_summary_TB_Subtype.csv"
  ),
  row.names = FALSE
)

## ========================================================================
## 10. Main figure: 8 states, all 10 T-cell subtypes
## ========================================================================

p_TCellSI_all <- ggplot(
  tcellsi_summary,
  aes(
    x = time_index,
    y = mean_score,
    group = subtype,
    color = subtype
  )
) +
  annotate(
    "rect",
    xmin = 5.5,
    xmax = 10.5,
    ymin = -Inf,
    ymax = Inf,
    fill = "grey96",
    color = NA
  ) +
  geom_vline(
    xintercept = 5.5,
    linewidth = 0.35,
    linetype = 2,
    color = "grey45"
  ) +
  geom_line(
    linewidth = 0.75,
    na.rm = TRUE
  ) +
  geom_point(
    size = 1.45,
    na.rm = TRUE
  ) +
  geom_errorbar(
    aes(
      ymin = mean_score - se_score,
      ymax = mean_score + se_score
    ),
    width = 0.12,
    linewidth = 0.25,
    alpha = 0.65,
    na.rm = TRUE
  ) +
  facet_wrap(
    ~ State,
    ncol = 4,
    scales = "fixed"
  ) +
  scale_x_continuous(
    breaks = seq_along(
      TIME_LEVELS
    ),
    labels = TIME_LABELS,
    limits = c(0.75, 10.25),
    expand = expansion(
      mult = c(0, 0)
    )
  ) +
  scale_color_manual(
    values = SUBTYPE_COLORS,
    breaks = T_SUBTYPE_ORDER,
    drop = FALSE
  ) +
  labs(
    x = NULL,
    y = "TCellSI score",
    color = "T-cell subtype"
  ) +
  theme_classic(
    base_size = 9
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 7
    ),
    axis.text.y = element_text(
      size = 7
    ),
    axis.title.y = element_text(
      size = 9
    ),
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 9
    ),
    legend.position = "right",
    legend.title = element_text(
      size = 9
    ),
    legend.text = element_text(
      size = 8
    ),
    legend.key.height = grid::unit(
      3.6,
      "mm"
    ),
    panel.spacing = grid::unit(
      3,
      "mm"
    )
  ) +
  guides(
    color = guide_legend(
      override.aes = list(
        linewidth = 1,
        size = 2
      ),
      ncol = 1
    )
  )

save_pdf(
  p_TCellSI_all,
  file.path(
    OUTDIR,
    "04_Figures",
    "Fig_TCellSI_TB_Subtype_8states_all_Tcells.pdf"
  ),
  width = 11.5,
  height = 7.2
)

## ========================================================================
## 11. Donor-level trajectories plus subtype mean
## ========================================================================

p_TCellSI_donor <- ggplot() +
  annotate(
    "rect",
    xmin = 5.5,
    xmax = 10.5,
    ymin = -Inf,
    ymax = Inf,
    fill = "grey96",
    color = NA
  ) +
  geom_vline(
    xintercept = 5.5,
    linewidth = 0.35,
    linetype = 2,
    color = "grey45"
  ) +
  geom_line(
    data = tcellsi_long,
    aes(
      x = time_index,
      y = Score,
      group = interaction(
        donor,
        subtype,
        drop = TRUE
      ),
      color = subtype
    ),
    linewidth = 0.28,
    alpha = 0.20,
    na.rm = TRUE
  ) +
  geom_line(
    data = tcellsi_summary,
    aes(
      x = time_index,
      y = mean_score,
      group = subtype,
      color = subtype
    ),
    linewidth = 0.85,
    na.rm = TRUE
  ) +
  geom_point(
    data = tcellsi_summary,
    aes(
      x = time_index,
      y = mean_score,
      color = subtype
    ),
    size = 1.4,
    na.rm = TRUE
  ) +
  facet_wrap(
    ~ State,
    ncol = 4,
    scales = "fixed"
  ) +
  scale_x_continuous(
    breaks = seq_along(
      TIME_LEVELS
    ),
    labels = TIME_LABELS,
    limits = c(0.75, 10.25),
    expand = expansion(
      mult = c(0, 0)
    )
  ) +
  scale_color_manual(
    values = SUBTYPE_COLORS,
    breaks = T_SUBTYPE_ORDER,
    drop = FALSE
  ) +
  labs(
    x = NULL,
    y = "TCellSI score",
    color = "T-cell subtype"
  ) +
  theme_classic(
    base_size = 9
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 7
    ),
    axis.text.y = element_text(
      size = 7
    ),
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 9
    ),
    legend.position = "right",
    legend.title = element_text(
      size = 9
    ),
    legend.text = element_text(
      size = 8
    ),
    panel.spacing = grid::unit(
      3,
      "mm"
    )
  ) +
  guides(
    color = guide_legend(
      override.aes = list(
        linewidth = 1,
        alpha = 1,
        size = 2
      )
    )
  )

save_pdf(
  p_TCellSI_donor,
  file.path(
    OUTDIR,
    "04_Figures",
    "Fig_TCellSI_TB_Subtype_8states_donor_trajectories.pdf"
  ),
  width = 11.5,
  height = 7.2
)

## ========================================================================
## 12. CD4-oriented figure
## ========================================================================

CD4_SUBTYPES <- c(
  "CD4_NaiveT",
  "CD4_Naive_CCR4",
  "Treg",
  "CD4_Memory_ZEB2",
  "CD4_Memory"
)

cd4_summary <- tcellsi_summary %>%
  filter(
    subtype %in% CD4_SUBTYPES
  ) %>%
  droplevels()

p_TCellSI_CD4 <- ggplot(
  cd4_summary,
  aes(
    x = time_index,
    y = mean_score,
    group = subtype,
    color = subtype
  )
) +
  annotate(
    "rect",
    xmin = 5.5,
    xmax = 10.5,
    ymin = -Inf,
    ymax = Inf,
    fill = "grey96",
    color = NA
  ) +
  geom_vline(
    xintercept = 5.5,
    linewidth = 0.35,
    linetype = 2,
    color = "grey45"
  ) +
  geom_line(
    linewidth = 0.85
  ) +
  geom_point(
    size = 1.55
  ) +
  geom_errorbar(
    aes(
      ymin = mean_score - se_score,
      ymax = mean_score + se_score
    ),
    width = 0.12,
    linewidth = 0.25,
    alpha = 0.65
  ) +
  facet_wrap(
    ~ State,
    ncol = 4,
    scales = "fixed"
  ) +
  scale_x_continuous(
    breaks = seq_along(
      TIME_LEVELS
    ),
    labels = TIME_LABELS,
    limits = c(0.75, 10.25),
    expand = expansion(
      mult = c(0, 0)
    )
  ) +
  scale_color_manual(
    values = SUBTYPE_COLORS[
      CD4_SUBTYPES
    ],
    breaks = CD4_SUBTYPES,
    drop = FALSE
  ) +
  labs(
    x = NULL,
    y = "TCellSI score",
    color = "CD4 T-cell subtype"
  ) +
  theme_classic(
    base_size = 9
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 7
    ),
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold"
    ),
    legend.position = "right",
    panel.spacing = grid::unit(
      3,
      "mm"
    )
  )

save_pdf(
  p_TCellSI_CD4,
  file.path(
    OUTDIR,
    "04_Figures",
    "Fig_TCellSI_TB_Subtype_CD4_8states.pdf"
  ),
  width = 10.5,
  height = 7
)

## ========================================================================
## 13. CD8 / gdT / DP-oriented figure
## ========================================================================

CD8_OTHER_SUBTYPES <- c(
  "CD8_NaiveT",
  "CD8_Memory_GZMK",
  "CD8_Memory",
  "gdT",
  "CD4_CD8_DP"
)

cd8_summary <- tcellsi_summary %>%
  filter(
    subtype %in% CD8_OTHER_SUBTYPES
  ) %>%
  droplevels()

p_TCellSI_CD8 <- ggplot(
  cd8_summary,
  aes(
    x = time_index,
    y = mean_score,
    group = subtype,
    color = subtype
  )
) +
  annotate(
    "rect",
    xmin = 5.5,
    xmax = 10.5,
    ymin = -Inf,
    ymax = Inf,
    fill = "grey96",
    color = NA
  ) +
  geom_vline(
    xintercept = 5.5,
    linewidth = 0.35,
    linetype = 2,
    color = "grey45"
  ) +
  geom_line(
    linewidth = 0.9
  ) +
  geom_point(
    size = 1.65
  ) +
  geom_errorbar(
    aes(
      ymin = mean_score - se_score,
      ymax = mean_score + se_score
    ),
    width = 0.12,
    linewidth = 0.25,
    alpha = 0.65
  ) +
  facet_wrap(
    ~ State,
    ncol = 4,
    scales = "fixed"
  ) +
  scale_x_continuous(
    breaks = seq_along(
      TIME_LEVELS
    ),
    labels = TIME_LABELS,
    limits = c(0.75, 10.25),
    expand = expansion(
      mult = c(0, 0)
    )
  ) +
  scale_color_manual(
    values = SUBTYPE_COLORS[
      CD8_OTHER_SUBTYPES
    ],
    breaks = CD8_OTHER_SUBTYPES,
    drop = FALSE
  ) +
  labs(
    x = NULL,
    y = "TCellSI score",
    color = "T-cell subtype"
  ) +
  theme_classic(
    base_size = 9
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 7
    ),
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold"
    ),
    legend.position = "right",
    panel.spacing = grid::unit(
      3,
      "mm"
    )
  )

save_pdf(
  p_TCellSI_CD8,
  file.path(
    OUTDIR,
    "04_Figures",
    "Fig_TCellSI_TB_Subtype_CD8_gdT_DP_8states.pdf"
  ),
  width = 10.5,
  height = 7
)

## ========================================================================
## 14. Save plotting objects
## ========================================================================

saveRDS(
  list(
    primary_method = primary_method_used,
    TCellSI_scores = TCellSI_scores,
    long_data = tcellsi_long,
    summary_data = tcellsi_summary,
    subtype_colors = SUBTYPE_COLORS,
    plots = list(
      all_T_subtypes = p_TCellSI_all,
      donor_trajectories = p_TCellSI_donor,
      CD4 = p_TCellSI_CD4,
      CD8_gdT_DP = p_TCellSI_CD8
    )
  ),
  file.path(
    OUTDIR,
    "03_TCellSI",
    "TCellSI_TB_Subtype_analysis_objects.rds"
  )
)

message(
  "TCellSI re-analysis using PROJ$TB_Subtype completed."
)

message(
  "Output directory: ",
  normalizePath(
    OUTDIR,
    mustWork = FALSE
  )
)

writeLines(capture.output(sessionInfo()), file.path(OUTDIR, "sessionInfo_updated.txt"))
