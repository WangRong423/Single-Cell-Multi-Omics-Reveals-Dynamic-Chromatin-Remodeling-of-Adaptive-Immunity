#!/usr/bin/env Rscript

# ==============================================================================
# 02 — RNA-based B-cell functional scoring
#
# Uses GeneExpressionMatrix from proj_TB_deviations_complete.rds.
# Does not use GeneScoreMatrix.
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
      length(source_file) ==
        1
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

source(
  file.path(
    SCRIPT_DIR,
    "Bcell_deviations_helpers_updated.R"
  )
)

suppressPackageStartupMessages({
  library(msigdbr)
})

RNA_PB_CACHE <- file.path(
  BCELL_OBJECT_DIR,
  "Bcell_GeneExpression_pseudobulk.rds"
)

MIN_GENES_PER_MODULE <- 5L
FUNCTION_FDR_DISPLAY <- 0.10

contrast_table <- tibble::tribble(
  ~Comparison,    ~DisplayLabel, ~Target,  ~Baseline,
  "Prime_D1",     "Prime D1",     "day1",   "day0",
  "Prime_D3",     "Prime D3",     "day3",   "day0",
  "Prime_D6",     "Prime D6",     "day6",   "day0",
  "Prime_D14",    "Prime D14",    "day14",  "day0",
  "Prime_D30",    "Prime D30",    "day30",  "day0",
  "Booster_D1",   "Booster D1",   "day31",  "day30",
  "Booster_D3",   "Booster D3",   "day33",  "day30",
  "Booster_D6",   "Booster D6",   "day36",  "day30",
  "Booster_D14",  "Booster D14",  "day44",  "day30"
)

contrast_levels <- contrast_table$DisplayLabel

get_msigdbr_collection <- function(
    collection,
    subcollection = NULL
) {
  argument_names <- names(
    formals(
      msigdbr::msigdbr
    )
  )

  arguments <- list(
    species = "Homo sapiens"
  )

  if (
    "db_species" %in%
      argument_names
  ) {
    arguments$db_species <- "HS"
  }

  if (
    "collection" %in%
      argument_names
  ) {
    arguments$collection <- collection
  } else {
    arguments$category <- collection
  }

  if (!is.null(subcollection)) {
    if (
      "subcollection" %in%
        argument_names
    ) {
      arguments$subcollection <-
        subcollection
    } else {
      arguments$subcategory <-
        subcollection
    }
  }

  result <- tryCatch(
    do.call(
      msigdbr::msigdbr,
      arguments
    ),
    error = function(e) NULL
  )

  if (
    collection ==
      "C5" &&
      !is.null(
        subcollection
      ) &&
      (
        is.null(result) ||
          nrow(result) ==
            0
      )
  ) {
    alternative <- if (
      subcollection ==
        "GO:BP"
    ) {
      "BP"
    } else {
      "GO:BP"
    }

    if (
      "subcollection" %in%
        argument_names
    ) {
      arguments$subcollection <-
        alternative
    } else {
      arguments$subcategory <-
        alternative
    }

    result <- do.call(
      msigdbr::msigdbr,
      arguments
    )
  }

  if (
    is.null(result) ||
      nrow(result) ==
        0
  ) {
    stop(
      "No msigdbr records were returned for ",
      collection,
      "/",
      subcollection
    )
  }

  tibble::as_tibble(
    result
  )
}


select_gene_set <- function(
    gene_set_table,
    exact_names = character(0),
    regex_patterns = character(0)
) {
  selected <- gene_set_table |>
    dplyr::filter(
      gs_name %in%
        exact_names
    )

  if (
    nrow(selected) ==
      0 &&
      length(regex_patterns) >
        0
  ) {
    selected <- gene_set_table |>
      dplyr::filter(
        grepl(
          paste(
            regex_patterns,
            collapse = "|"
          ),
          gs_name,
          perl = TRUE
        )
      )
  }

  unique(
    as.character(
      selected$gene_symbol
    )
  )
}


build_bcell_modules <- function() {
  hallmark <- get_msigdbr_collection(
    "H"
  )

  go_bp <- get_msigdbr_collection(
    "C5",
    "GO:BP"
  )

  immune_signatures <- get_msigdbr_collection(
    "C7"
  )

  list(
    "BCR signaling" = select_gene_set(
      go_bp,
      exact_names =
        "GOBP_B_CELL_RECEPTOR_SIGNALING_PATHWAY",
      regex_patterns =
        "B_CELL_RECEPTOR_SIGNALING_PATHWAY"
    ),
    "B-cell activation" = select_gene_set(
      go_bp,
      exact_names =
        "GOBP_B_CELL_ACTIVATION",
      regex_patterns =
        "^GOBP_B_CELL_ACTIVATION$"
    ),
    "Antigen presentation MHC-II" =
      select_gene_set(
        go_bp,
        regex_patterns =
          "ANTIGEN_PROCESSING_AND_PRESENTATION.*MHC_CLASS_II"
      ),
    "Germinal-center program" =
      select_gene_set(
        go_bp,
        exact_names =
          "GOBP_GERMINAL_CENTER_FORMATION",
        regex_patterns =
          "GERMINAL_CENTER_FORMATION"
      ),
    "Memory B-cell program" =
      select_gene_set(
        immune_signatures,
        exact_names =
          "GSE13411_NAIVE_VS_MEMORY_BCELL_DN",
        regex_patterns =
          "GSE13411_NAIVE_VS_MEMORY_BCELL_DN"
      ),
    "Plasma-cell differentiation" =
      select_gene_set(
        go_bp,
        exact_names =
          "GOBP_PLASMA_CELL_DIFFERENTIATION",
        regex_patterns =
          "PLASMA_CELL_DIFFERENTIATION"
      ),
    "Plasma vs memory program" =
      select_gene_set(
        immune_signatures,
        exact_names =
          "GSE13411_PLASMA_CELL_VS_MEMORY_BCELL_UP",
        regex_patterns =
          "GSE13411_PLASMA_CELL_VS_MEMORY_BCELL_UP"
      ),
    "Immunoglobulin production" =
      select_gene_set(
        go_bp,
        exact_names =
          "GOBP_IMMUNOGLOBULIN_PRODUCTION",
        regex_patterns =
          "IMMUNOGLOBULIN_PRODUCTION"
      ),
    "Secretory / UPR" =
      select_gene_set(
        hallmark,
        exact_names =
          "HALLMARK_UNFOLDED_PROTEIN_RESPONSE"
      ),
    "Proliferation" = unique(
      c(
        select_gene_set(
          hallmark,
          exact_names =
            "HALLMARK_E2F_TARGETS"
        ),
        select_gene_set(
          hallmark,
          exact_names =
            "HALLMARK_G2M_CHECKPOINT"
        )
      )
    ),
    "Type-I interferon" =
      select_gene_set(
        hallmark,
        exact_names =
          "HALLMARK_INTERFERON_ALPHA_RESPONSE"
      ),
    "NF-kB activation" =
      select_gene_set(
        hallmark,
        exact_names =
          "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
      )
  )
}


score_modules <- function(
    expression_matrix,
    module_list
) {
  expression_matrix <- as.matrix(
    expression_matrix
  )

  expression_z <- row_zscore(
    expression_matrix
  )

  score_list <- list()
  retained_list <- list()
  retained_index <- 1L

  for (module_name in names(module_list)) {
    retained_genes <- intersect(
      unique(
        module_list[[module_name]]
      ),
      rownames(
        expression_z
      )
    )

    if (
      length(
        retained_genes
      ) <
        MIN_GENES_PER_MODULE
    ) {
      message(
        "Skipping module ",
        module_name,
        ": ",
        length(
          retained_genes
        ),
        " genes retained."
      )

      next
    }

    score_list[[module_name]] <- colMeans(
      expression_z[
        retained_genes,
        ,
        drop = FALSE
      ],
      na.rm = TRUE
    )

    retained_list[[retained_index]] <-
      tibble::tibble(
        Module = module_name,
        Gene = retained_genes
      )

    retained_index <- retained_index + 1L
  }

  if (length(score_list) == 0) {
    stop(
      "No B-cell functional module retained enough genes."
    )
  }

  list(
    ScoreMatrix = do.call(
      rbind,
      score_list
    ),
    RetainedGenes = dplyr::bind_rows(
      retained_list
    )
  )
}

project <- load_deviations_project()

bcell_object <- prepare_b_project(
  project
)

project_b <- bcell_object$Project
pseudobulk_metadata <-
  bcell_object$PBMetadata

if (
  !"GeneExpressionMatrix" %in%
    ArchR::getAvailableMatrices(
      project_b
    )
) {
  stop(
    "GeneExpressionMatrix is absent from the project."
  )
}

rna_pb <- get_bcell_pseudobulk(
  project_b = project_b,
  pseudobulk_metadata =
    pseudobulk_metadata,
  use_matrix = "GeneExpressionMatrix",
  cache_file = RNA_PB_CACHE,
  divide_n = TRUE,
  scale_to = 10000,
  preferred_assays = character(0)
)

expression_matrix <- log2(
  as.matrix(
    rna_pb$Matrix
  ) +
    1
)

expression_metadata <- rna_pb$Metadata

gene_names <- extract_feature_names(
  summarized_experiment =
    rna_pb$SE,
  matrix_object =
    rna_pb$Matrix
)

rownames(
  expression_matrix
) <- gene_names

if (
  anyDuplicated(
    rownames(
      expression_matrix
    )
  )
) {
  duplicate_index <- split(
    seq_len(
      nrow(
        expression_matrix
      )
    ),
    rownames(
      expression_matrix
    )
  )

  expression_matrix <- do.call(
    rbind,
    lapply(
      duplicate_index,
      function(row_index) {
        colMeans(
          expression_matrix[
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
    expression_matrix
  ) <- names(
    duplicate_index
  )
}

module_list <- build_bcell_modules()

score_object <- score_modules(
  expression_matrix =
    expression_matrix,
  module_list =
    module_list
)

score_matrix <- score_object$ScoreMatrix
retained_genes <- score_object$RetainedGenes

write.csv(
  retained_genes,
  file.path(
    BCELL_TABLE_DIR,
    "Bcell_function_modules_retained_genes.csv"
  ),
  row.names = FALSE
)

score_long <- as.data.frame(
  t(
    score_matrix
  )
) |>
  tibble::rownames_to_column(
    "B_PB_Group"
  ) |>
  tidyr::pivot_longer(
    cols = -B_PB_Group,
    names_to = "Module",
    values_to = "FunctionScore"
  ) |>
  dplyr::left_join(
    expression_metadata,
    by = "B_PB_Group"
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
    )
  )

score_time_mean <- mean_sem(
  data = score_long,
  value = "FunctionScore",
  groups = c(
    "Module",
    "TB_Subtype",
    "Timepoint"
  )
)

module_order <- c(
  "BCR signaling",
  "B-cell activation",
  "Antigen presentation MHC-II",
  "Germinal-center program",
  "Memory B-cell program",
  "Plasma-cell differentiation",
  "Plasma vs memory program",
  "Immunoglobulin production",
  "Secretory / UPR",
  "Proliferation",
  "Type-I interferon",
  "NF-kB activation"
)

module_order <- module_order[
  module_order %in%
    unique(
      score_time_mean$Module
    )
]

function_limma <- run_paired_limma(
  feature_matrix = score_matrix,
  metadata = expression_metadata,
  feature_name_column = "Module",
  contrast_table = contrast_table
) |>
  dplyr::transmute(
    Module = Module,
    ScoreDifference = logFC,
    AverageScore = AveExpr,
    TStatistic = t,
    PValue = P.Value,
    FDR = adj.P.Val,
    TB_Subtype = TB_Subtype,
    Comparison = Comparison,
    DisplayLabel = DisplayLabel,
    Target = Target,
    Baseline = Baseline,
    nPairedDonors = nPairedDonors,
    PairedDonors = PairedDonors
  ) |>
  dplyr::mutate(
    DisplayLabel = factor(
      DisplayLabel,
      levels = contrast_levels
    ),
    TB_Subtype = factor(
      TB_Subtype,
      levels = B_SUBTYPES
    ),
    Module = factor(
      Module,
      levels = rev(
        module_order
      )
    ),
    MinusLog10FDR = -log10(
      pmax(
        FDR,
        1e-300
      )
    ),
    FDRClass = dplyr::if_else(
      FDR <=
        FUNCTION_FDR_DISPLAY,
      paste0(
        "FDR ≤ ",
        FUNCTION_FDR_DISPLAY
      ),
      paste0(
        "FDR > ",
        FUNCTION_FDR_DISPLAY
      )
    )
  )

write.csv(
  score_long,
  gzfile(
    file.path(
      BCELL_TABLE_DIR,
      "Bcell_function_score_sample_values.csv.gz"
    )
  ),
  row.names = FALSE
)

write.csv(
  function_limma,
  gzfile(
    file.path(
      BCELL_TABLE_DIR,
      "Bcell_function_score_paired_limma.csv.gz"
    )
  ),
  row.names = FALSE
)

p_lines <- score_time_mean |>
  dplyr::mutate(
    Module = factor(
      Module,
      levels = module_order
    ),
    Timepoint = factor(
      Timepoint,
      levels = TIME_LEVELS,
      labels = TIME_LABELS
    )
  ) |>
  ggplot2::ggplot(
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
    color = "grey70"
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
    linewidth = 0.70
  ) +
  ggplot2::geom_point(
    size = 1.30
  ) +
  ggplot2::facet_wrap(
    ~ Module,
    ncol = 4,
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
    title = "RNA-based functional dynamics of B-cell subtypes",
    subtitle = "GeneExpressionMatrix; score = mean gene-wise z-score",
    x = NULL,
    y = "Functional module score",
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
      size = 6.5
    ),
    strip.background = ggplot2::element_rect(
      fill = "grey95",
      color = "grey45",
      linewidth = 0.30
    ),
    strip.text = ggplot2::element_text(
      face = "bold",
      size = 7
    ),
    panel.border = ggplot2::element_rect(
      fill = NA,
      color = "grey65",
      linewidth = 0.25
    ),
    plot.title = ggplot2::element_text(
      face = "bold",
      hjust = 0.5
    ),
    plot.subtitle = ggplot2::element_text(
      hjust = 0.5,
      size = 7.5
    ),
    legend.position = "bottom"
  )

finite_score_difference <- abs(
  function_limma$ScoreDifference[
    is.finite(
      function_limma$ScoreDifference
    )
  ]
)

score_color_limit <- if (
  length(
    finite_score_difference
  ) >
    0
) {
  max(
    0.20,
    as.numeric(
      stats::quantile(
        finite_score_difference,
        0.98,
        na.rm = TRUE
      )
    )
  )
} else {
  1
}

p_limma <- ggplot2::ggplot(
  function_limma,
  ggplot2::aes(
    x = DisplayLabel,
    y = Module
  )
) +
  ggplot2::geom_tile(
    fill = "grey98",
    color = "white",
    linewidth = 0.20
  ) +
  ggplot2::geom_point(
    ggplot2::aes(
      fill = ScoreDifference,
      size = MinusLog10FDR,
      alpha = FDRClass
    ),
    shape = 21,
    color = "grey20",
    stroke = 0.25
  ) +
  ggplot2::facet_wrap(
    ~ TB_Subtype,
    ncol = 2
  ) +
  ggplot2::scale_fill_gradient2(
    low = "#3B6FB6",
    mid = "white",
    high = "#C83E4D",
    midpoint = 0,
    limits = c(
      -score_color_limit,
      score_color_limit
    ),
    oob = scales::squish,
    name = "Score\ndifference"
  ) +
  ggplot2::scale_size_continuous(
    range = c(
      1.2,
      6
    ),
    name = expression(
      -log[10] *
        "(FDR)"
    )
  ) +
  ggplot2::scale_alpha_manual(
    values = stats::setNames(
      c(
        1,
        0.28
      ),
      c(
        paste0(
          "FDR ≤ ",
          FUNCTION_FDR_DISPLAY
        ),
        paste0(
          "FDR > ",
          FUNCTION_FDR_DISPLAY
        )
      )
    ),
    name = NULL
  ) +
  ggplot2::labs(
    title = "Donor-paired changes in B-cell functional programs",
    x = NULL,
    y = NULL
  ) +
  ggplot2::theme_classic(
    base_size = 8.5
  ) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 45,
      hjust = 1,
      size = 6.4
    ),
    axis.text.y = ggplot2::element_text(
      size = 6.6
    ),
    axis.ticks = ggplot2::element_blank(),
    axis.line = ggplot2::element_blank(),
    strip.background = ggplot2::element_rect(
      fill = "grey95",
      color = "grey45",
      linewidth = 0.30
    ),
    strip.text = ggplot2::element_text(
      face = "bold"
    ),
    panel.border = ggplot2::element_rect(
      fill = NA,
      color = "grey65",
      linewidth = 0.25
    ),
    plot.title = ggplot2::element_text(
      face = "bold",
      hjust = 0.5
    )
  )

p_function <- p_lines /
  p_limma +
  patchwork::plot_layout(
    heights = c(
      1.20,
      1
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
    "Fig4B_Bcell_RNA_function_scores.pdf"
  ),
  plot = p_function,
  width = 240,
  height = 280,
  units = "mm",
  limitsize = FALSE
)

ggplot2::ggsave(
  filename = file.path(
    BCELL_PLOT_DIR,
    "Fig4B_Bcell_RNA_function_scores.png"
  ),
  plot = p_function,
  width = 240,
  height = 280,
  units = "mm",
  dpi = 600,
  limitsize = FALSE
)

saveRDS(
  list(
    ExpressionMatrix =
      expression_matrix,
    Metadata =
      expression_metadata,
    Modules =
      module_list,
    RetainedGenes =
      retained_genes,
    ScoreMatrix =
      score_matrix,
    ScoreLong =
      score_long,
    ScoreTimeMean =
      score_time_mean,
    FunctionLimma =
      function_limma,
    pLines =
      p_lines,
    pLimma =
      p_limma,
    pFunction =
      p_function
  ),
  file.path(
    BCELL_OBJECT_DIR,
    "Fig4B_Bcell_RNA_function_objects.rds"
  ),
  compress = FALSE
)

cat(
  "\nB-cell RNA function-score analysis completed:\n",
  file.path(
    BCELL_PLOT_DIR,
    "Fig4B_Bcell_RNA_function_scores.pdf"
  ),
  "\n"
)

writeLines(capture.output(sessionInfo()), file.path(BCELL_OBJECT_DIR, "Bcell_function_sessionInfo.txt"))
saveRDS(module_list, file.path(BCELL_OBJECT_DIR, "Bcell_function_module_definitions.rds"))
