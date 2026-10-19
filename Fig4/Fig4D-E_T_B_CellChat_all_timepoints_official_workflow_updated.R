# Locate the companion file for source() and Rscript.
.local_source_file <- tryCatch(sys.frame(1)$ofile, error=function(e) NULL)
.local_args <- grep("^--file=", commandArgs(FALSE), value=TRUE)
.local_script_dir <- if (!is.null(.local_source_file)) {
  dirname(normalizePath(.local_source_file))
} else if (length(.local_args)) {
  dirname(normalizePath(sub("^--file=", "", .local_args[1])))
} else getwd()
source(file.path(.local_script_dir, "Tcell_updated_labels_helpers.R"))


## =============================================================================
## Traditional CellChat workflow across all timepoints
## Data source: ArchR project proj_TB_fig3
##
## Strategy:
##   1. Pool all donors.
##   2. Analyze all 10 timepoints separately with standard CellChat workflow.
##   3. Group cells by TB_Subtype.
##   4. Export official CellChat plots for each timepoint.
##   5. Merge all CellChat objects and export cross-timepoint comparison plots.
## =============================================================================

suppressPackageStartupMessages({
  library(ArchR)
  library(SummarizedExperiment)
  library(Matrix)
  library(CellChat)
  library(future)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

proj <- select_tb_project()

TIME_LEVELS <- c(
  "day0", "day1", "day3", "day6", "day14",
  "day30", "day31", "day33", "day36", "day44"
)

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

B_SUBTYPES <- c(
  "Naive_B",
  "Intermediate_B",
  "Memory_B",
  "PlasmaB"
)

ALL_SUBTYPES <- c(T_SUBTYPES, B_SUBTYPES)

MAX_CELLS_PER_SUBTYPE <- 300
MIN_CELLS_PER_SUBTYPE <- 20

NBOOT <- 20
RANDOM_SEED <- 20260802
RECOMPUTE <- TRUE

OUTDIR <- file.path(UPDATED_BASE, "T_B_CellChat_all_timepoints_official_updated_labels")
OBJECT_DIR <- file.path(OUTDIR, "01_objects")
TABLE_DIR <- file.path(OUTDIR, "02_tables")
FIG_DIR <- file.path(OUTDIR, "03_figures")
QC_DIR <- file.path(OUTDIR, "04_qc")

for (x in c(OUTDIR, OBJECT_DIR, TABLE_DIR, FIG_DIR, QC_DIR)) {
  dir.create(x, recursive = TRUE, showWarnings = FALSE)
}

save_pdf_plot <- function(filename, width, height, expr) {
  grDevices::pdf(
    file = filename,
    width = width,
    height = height,
    useDingbats = FALSE
  )
  on.exit(grDevices::dev.off(), add = TRUE)
  eval.parent(substitute(expr))
}

collapse_sparse_rows <- function(mat, feature_names) {
  keep <- !is.na(feature_names) & feature_names != ""
  mat <- mat[keep, , drop = FALSE]
  feature_names <- feature_names[keep]

  unique_features <- unique(feature_names)

  row_map <- Matrix::sparseMatrix(
    i = match(feature_names, unique_features),
    j = seq_along(feature_names),
    x = 1,
    dims = c(length(unique_features), length(feature_names))
  )

  collapsed <- row_map %*% mat
  rownames(collapsed) <- unique_features
  colnames(collapsed) <- colnames(mat)
  collapsed
}

normalize_log1p_1e4 <- function(count_mat) {
  sf <- 1e4 / pmax(Matrix::colSums(count_mat), 1)
  norm_mat <- count_mat %*% Matrix::Diagonal(x = sf)
  dimnames(norm_mat) <- dimnames(count_mat)
  norm_mat@x <- log1p(norm_mat@x)
  norm_mat
}

classify_group <- function(x) {
  dplyr::case_when(
    x %in% T_SUBTYPES ~ "T",
    x %in% B_SUBTYPES ~ "B",
    TRUE ~ "Other"
  )
}

get_top_pathways <- function(cellchat_obj, top_n = 6) {
  pathway_df <- tryCatch(
    {
      CellChat::subsetCommunication(
        cellchat_obj,
        slot.name = "netP"
      )
    },
    error = function(e) {
      NULL
    }
  )

  if (is.null(pathway_df) || nrow(pathway_df) == 0) {
    return(character(0))
  }

  pathway_df %>%
    dplyr::group_by(pathway_name) %>%
    dplyr::summarise(total_prob = sum(prob, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(total_prob)) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::pull(pathway_name)
}

meta <- as.data.frame(ArchR::getCellColData(proj))
meta$cell <- rownames(meta)

required_cols <- c("Timepoint", "TB_Subtype")
missing_cols <- setdiff(required_cols, colnames(meta))
if (length(missing_cols) > 0) {
  stop(
    "Missing metadata columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

meta$TB_Subtype <- rename_tb_labels(meta$TB_Subtype,
  file.path(QC_DIR, "Subtype_label_mapping.csv"))

meta_use <- meta %>%
  dplyr::transmute(
    cell = cell,
    timepoint = as.character(Timepoint),
    subtype = as.character(TB_Subtype)
  ) %>%
  dplyr::filter(
    timepoint %in% TIME_LEVELS,
    subtype %in% ALL_SUBTYPES
  )

if (nrow(meta_use) == 0) {
  stop("No cells retained after filtering metadata.")
}

write.csv(
  meta_use %>%
    dplyr::count(timepoint, subtype, name = "n_cells"),
  file.path(QC_DIR, "cell_counts_before_downsampling.csv"),
  row.names = FALSE
)

if (!"GeneExpressionMatrix" %in% ArchR::getAvailableMatrices(proj)) {
  stop("GeneExpressionMatrix is unavailable in the ArchR project.")
}

se_rna <- ArchR::getMatrixFromProject(
  ArchRProj = proj,
  useMatrix = "GeneExpressionMatrix"
)

rna_counts <- SummarizedExperiment::assay(se_rna)
rna_rowdata <- as.data.frame(SummarizedExperiment::rowData(se_rna))

candidate_gene_cols <- intersect(
  c("symbol", "gene_symbol", "gene_name", "gene", "name", "feature"),
  colnames(rna_rowdata)
)

if (length(candidate_gene_cols) > 0) {
  gene_names <- as.character(rna_rowdata[[candidate_gene_cols[1]]])
} else {
  gene_names <- rownames(se_rna)
}

gene_names <- toupper(trimws(gene_names))
gene_names[gene_names == ""] <- NA_character_

common_cells <- intersect(meta_use$cell, colnames(rna_counts))
if (length(common_cells) == 0) {
  stop("No common cells between metadata and GeneExpressionMatrix.")
}

meta_use <- meta_use[match(common_cells, meta_use$cell), , drop = FALSE]
rna_counts <- rna_counts[, common_cells, drop = FALSE]

stopifnot(identical(colnames(rna_counts), meta_use$cell))

rna_gene <- collapse_sparse_rows(
  mat = rna_counts,
  feature_names = gene_names
)

writeLines(
  paste(
    "Prepared RNA matrix:",
    nrow(rna_gene), "genes x", ncol(rna_gene), "cells"
  ),
  con = file.path(QC_DIR, "RNA_matrix_summary.txt")
)

data("CellChatDB.human", package = "CellChat")
CellChatDB.use <- CellChat::subsetDB(CellChatDB.human)

run_cellchat_one_timepoint <- function(time_value) {

  object_file <- file.path(OBJECT_DIR, paste0("CellChat_", time_value, ".rds"))
  lr_file <- file.path(TABLE_DIR, paste0("LR_", time_value, ".csv"))
  pathway_file <- file.path(TABLE_DIR, paste0("Pathway_", time_value, ".csv"))

  if (file.exists(object_file) && !RECOMPUTE) {
    message("Loading cached CellChat object: ", time_value)
    cellchat <- readRDS(object_file)

    lr_tbl <- CellChat::subsetCommunication(cellchat)
    lr_tbl$timepoint <- time_value

    pathway_tbl <- tryCatch(
      {
        x <- CellChat::subsetCommunication(cellchat, slot.name = "netP")
        x$timepoint <- time_value
        x
      },
      error = function(e) {
        data.frame()
      }
    )

    return(
      list(
        object = cellchat,
        lr = lr_tbl,
        pathway = pathway_tbl
      )
    )
  }

  set.seed(RANDOM_SEED + match(time_value, TIME_LEVELS))

  meta_time <- meta_use %>%
    dplyr::filter(timepoint == time_value) %>%
    dplyr::mutate(rand = stats::runif(dplyr::n())) %>%
    dplyr::group_by(subtype) %>%
    dplyr::arrange(rand, .by_group = TRUE) %>%
    dplyr::slice_head(n = MAX_CELLS_PER_SUBTYPE) %>%
    dplyr::ungroup() %>%
    dplyr::select(-rand)

  subtype_counts <- meta_time %>%
    dplyr::count(subtype, name = "n_cells")

  keep_subtypes <- subtype_counts %>%
    dplyr::filter(n_cells >= MIN_CELLS_PER_SUBTYPE) %>%
    dplyr::pull(subtype)

  meta_time <- meta_time %>%
    dplyr::filter(subtype %in% keep_subtypes)

  if (nrow(meta_time) == 0) {
    stop("No cells retained for timepoint: ", time_value)
  }

  write.csv(
    meta_time %>% dplyr::count(subtype, name = "n_cells"),
    file.path(QC_DIR, paste0("cell_counts_used_", time_value, ".csv")),
    row.names = FALSE
  )

  cells_use <- intersect(meta_time$cell, colnames(rna_gene))
  meta_time <- meta_time[match(cells_use, meta_time$cell), , drop = FALSE]
  expr_time <- rna_gene[, cells_use, drop = FALSE]

  expr_time <- normalize_log1p_1e4(expr_time)

  meta_cellchat <- data.frame(
    labels = meta_time$subtype,
    row.names = meta_time$cell,
    stringsAsFactors = FALSE
  )

  message(
    "Running CellChat for ", time_value, ": ",
    ncol(expr_time), " cells, ",
    length(unique(meta_cellchat$labels)), " subtypes."
  )

  cellchat <- CellChat::createCellChat(
    object = expr_time,
    meta = meta_cellchat,
    group.by = "labels"
  )

  cellchat@DB <- CellChatDB.use
  cellchat <- CellChat::subsetData(cellchat)
  cellchat <- CellChat::identifyOverExpressedGenes(cellchat)
  cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)
  cellchat <- CellChat::computeCommunProb(
    cellchat,
    type = "truncatedMean",
    trim = 0.1,
    raw.use = TRUE,
    population.size = FALSE,
    nboot = NBOOT,
    seed.use = RANDOM_SEED
  )
  cellchat <- CellChat::filterCommunication(
    cellchat,
    min.cells = MIN_CELLS_PER_SUBTYPE
  )
  cellchat <- CellChat::computeCommunProbPathway(cellchat)
  cellchat <- CellChat::aggregateNet(cellchat)
  cellchat <- CellChat::netAnalysis_computeCentrality(
    cellchat,
    slot.name = "netP"
  )

  saveRDS(cellchat, object_file)

  lr_tbl <- CellChat::subsetCommunication(cellchat)
  lr_tbl$timepoint <- time_value
  write.csv(lr_tbl, lr_file, row.names = FALSE)

  pathway_tbl <- tryCatch(
    {
      x <- CellChat::subsetCommunication(cellchat, slot.name = "netP")
      x$timepoint <- time_value
      write.csv(x, pathway_file, row.names = FALSE)
      x
    },
    error = function(e) {
      data.frame()
    }
  )

  list(
    object = cellchat,
    lr = lr_tbl,
    pathway = pathway_tbl
  )
}

future::plan("multisession", workers = 4)
options(future.globals.maxSize = 20 * 1024^3)

cellchat_results <- vector("list", length(TIME_LEVELS))
names(cellchat_results) <- TIME_LEVELS

for (time_value in TIME_LEVELS) {
  cellchat_results[[time_value]] <- run_cellchat_one_timepoint(time_value)
}

future::plan("sequential")

plot_one_timepoint <- function(cellchat_obj, time_value) {

  plot_dir <- file.path(FIG_DIR, time_value)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

  group_size <- as.numeric(table(cellchat_obj@idents))

  save_pdf_plot(
    file.path(plot_dir, paste0(time_value, "_circle_count.pdf")),
    width = 8,
    height = 8,
    CellChat::netVisual_circle(
      net = cellchat_obj@net$count,
      vertex.weight = group_size,
      weight.scale = TRUE,
      label.edge = FALSE,
      title.name = paste0(time_value, " | Number of interactions")
    )
  )

  save_pdf_plot(
    file.path(plot_dir, paste0(time_value, "_circle_weight.pdf")),
    width = 8,
    height = 8,
    CellChat::netVisual_circle(
      net = cellchat_obj@net$weight,
      vertex.weight = group_size,
      weight.scale = TRUE,
      label.edge = FALSE,
      title.name = paste0(time_value, " | Interaction weights")
    )
  )

  save_pdf_plot(
    file.path(plot_dir, paste0(time_value, "_heatmap_count.pdf")),
    width = 9,
    height = 8,
    print(
      CellChat::netVisual_heatmap(
        object = cellchat_obj,
        measure = "count",
        color.heatmap = "Reds"
      )
    )
  )

  save_pdf_plot(
    file.path(plot_dir, paste0(time_value, "_heatmap_weight.pdf")),
    width = 9,
    height = 8,
    print(
      CellChat::netVisual_heatmap(
        object = cellchat_obj,
        measure = "weight",
        color.heatmap = "Reds"
      )
    )
  )

  save_pdf_plot(
    file.path(plot_dir, paste0(time_value, "_signalingRole_scatter.pdf")),
    width = 7,
    height = 6,
    print(
      CellChat::netAnalysis_signalingRole_scatter(
        object = cellchat_obj
      )
    )
  )

  save_pdf_plot(
    file.path(plot_dir, paste0(time_value, "_signalingRole_heatmap_outgoing.pdf")),
    width = 8,
    height = 9,
    print(
      CellChat::netAnalysis_signalingRole_heatmap(
        object = cellchat_obj,
        pattern = "outgoing"
      )
    )
  )

  save_pdf_plot(
    file.path(plot_dir, paste0(time_value, "_signalingRole_heatmap_incoming.pdf")),
    width = 8,
    height = 9,
    print(
      CellChat::netAnalysis_signalingRole_heatmap(
        object = cellchat_obj,
        pattern = "incoming"
      )
    )
  )

  id_levels <- levels(cellchat_obj@idents)

  b_idx <- which(id_levels %in% B_SUBTYPES)
  t_idx <- which(id_levels %in% T_SUBTYPES)

  if (length(b_idx) > 0 && length(t_idx) > 0) {

    save_pdf_plot(
      file.path(plot_dir, paste0(time_value, "_bubble_B_to_T.pdf")),
      width = 10,
      height = 8,
      print(
        CellChat::netVisual_bubble(
          object = cellchat_obj,
          sources.use = b_idx,
          targets.use = t_idx,
          remove.isolate = FALSE
        )
      )
    )

    save_pdf_plot(
      file.path(plot_dir, paste0(time_value, "_bubble_T_to_B.pdf")),
      width = 10,
      height = 8,
      print(
        CellChat::netVisual_bubble(
          object = cellchat_obj,
          sources.use = t_idx,
          targets.use = b_idx,
          remove.isolate = FALSE
        )
      )
    )
  }

  top_pathways <- get_top_pathways(cellchat_obj, top_n = 6)

  if (length(top_pathways) > 0) {
    for (pathway_name in top_pathways) {

      safe_name <- gsub("[^A-Za-z0-9._-]+", "_", pathway_name)

      save_pdf_plot(
        file.path(plot_dir, paste0(time_value, "_contribution_", safe_name, ".pdf")),
        width = 7,
        height = 5,
        print(
          CellChat::netAnalysis_contribution(
            object = cellchat_obj,
            signaling = pathway_name
          )
        )
      )

      save_pdf_plot(
        file.path(plot_dir, paste0(time_value, "_aggregate_circle_", safe_name, ".pdf")),
        width = 8,
        height = 8,
        CellChat::netVisual_aggregate(
          object = cellchat_obj,
          signaling = pathway_name,
          layout = "circle"
        )
      )

      save_pdf_plot(
        file.path(plot_dir, paste0(time_value, "_aggregate_hierarchy_", safe_name, ".pdf")),
        width = 10,
        height = 8,
        CellChat::netVisual_aggregate(
          object = cellchat_obj,
          signaling = pathway_name,
          layout = "hierarchy"
        )
      )
    }
  }
}

TIME_LEVELS <- c(
  "day0", "day1", "day3", "day6", "day14",
  "day30", "day31", "day33", "day36", "day44"
)

T_SUBTYPES <- c(
  "CD4_NaiveT",
  "CD4_Naive_CCR4",
  "CD4_Memory",
  "CD4_Memory_ZEB2",
  "Treg",
  "CD8_NaiveT",
  "CD8_Memory",
  "CD8_Memory_GZMK",
  "gdT",
  "CD4_CD8_DP"
)

B_SUBTYPES <- c(
  "Naive_B",
  "Intermediate_B",
  "Memory_B",
  "PlasmaB"
)

# 保留 TCellSI 中相应亚群的颜色
SUBTYPE_COLORS <- c(
  CD4_NaiveT       = "#009E73",
  CD4_Naive_CCR4   = "#CC79A7",
  CD4_Memory      = "#4DAF4A",
  CD4_Memory_ZEB2 = "#E41A1C",
  Treg            = "#F781BF",
  CD8_NaiveT      = "#D8A500",
  CD8_Memory     = "#0072B2",
  CD8_Memory_GZMK = "#D55E00",
  gdT             = "#A6761D",
  CD4_CD8_DP      = "#6A3D9A",
  Naive_B         = "#4E79A7",
  Intermediate_B  = "#E15759",
  Memory_B        = "#59A14F",
  PlasmaB         = "#F28E2B"
)

OUTDIR <- paste0(
  "/public/users/wangrong/home/wangrong/results/T_B/",
  "TB_final_downstream/",
  "T_B_CellChat_all_timepoints_official_updated_labels"
)

# 单独保存重绘结果
REPLOT_DIR <- file.path(OUTDIR, "05_replot")
FIG_DIR <- file.path(REPLOT_DIR, "Figures")
TABLE_DIR <- file.path(REPLOT_DIR, "Tables")
OBJECT_DIR <- file.path(REPLOT_DIR, "Objects")
QC_DIR <- file.path(REPLOT_DIR, "QC")

for (directory in c(
  REPLOT_DIR, FIG_DIR, TABLE_DIR, OBJECT_DIR, QC_DIR
)) {
  dir.create(
    directory,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

# ==============================================================================
# 2. 检查已有对象
# ==============================================================================

if (!exists("cellchat_results")) {
  stop("当前 R 环境中没有 cellchat_results，请先载入已有分析结果。")
}

object_check <- data.frame(
  Timepoint = TIME_LEVELS,
  HasObject = vapply(
    TIME_LEVELS,
    function(tp) {
      result <- cellchat_results[[tp]]

      !is.null(result) &&
        !is.null(result$object) &&
        inherits(result$object, "CellChat")
    },
    logical(1)
  )
)

print(object_check, row.names = FALSE)

write.csv(
  object_check,
  file.path(QC_DIR, "CellChat_object_availability.csv"),
  row.names = FALSE
)

if (!all(object_check$HasObject)) {
  stop(
    "缺少有效 CellChat 对象：",
    paste(
      object_check$Timepoint[!object_check$HasObject],
      collapse = ", "
    )
  )
}

cellchat_list <- setNames(
  lapply(
    TIME_LEVELS,
    function(tp) cellchat_results[[tp]]$object
  ),
  TIME_LEVELS
)

# 检查对象中的真实标签，避免仅修改图例
observed_labels <- unique(unlist(
  lapply(cellchat_list, function(x) levels(x@idents)),
  use.names = FALSE
))

unknown_labels <- setdiff(
  observed_labels,
  names(SUBTYPE_COLORS)
)

if (length(unknown_labels) > 0L) {
  stop(
    "CellChat 对象存在未配置的标签：",
    paste(unknown_labels, collapse = ", "),
    "。请确认这些对象已使用更新后的 TB_Subtype 计算。"
  )
}

# ==============================================================================
# 3. 错误记录与安全绘图
# ==============================================================================

error_file <- file.path(QC_DIR, "CellChat_plot_errors.csv")

write.csv(
  data.frame(
    File = character(),
    Error = character()
  ),
  error_file,
  row.names = FALSE
)

record_error <- function(filename, error_message) {
  write.table(
    data.frame(
      File = filename,
      Error = error_message,
      stringsAsFactors = FALSE
    ),
    file = error_file,
    sep = ",",
    row.names = FALSE,
    col.names = FALSE,
    append = TRUE,
    quote = TRUE
  )
}

save_pdf_plot <- function(filename, width, height, expr) {

  plot_expression <- substitute(expr)
  calling_environment <- parent.frame()

  dir.create(
    dirname(filename),
    recursive = TRUE,
    showWarnings = FALSE
  )

  device_id <- NULL
  success <- FALSE

  tryCatch(
    {
      grDevices::pdf(
        file = filename,
        width = width,
        height = height,
        useDingbats = FALSE
      )

      device_id <- grDevices::dev.cur()

      eval(
        plot_expression,
        envir = calling_environment
      )

      success <- TRUE
    },
    error = function(e) {
      record_error(filename, conditionMessage(e))

      message(
        "跳过失败图：", basename(filename),
        "\n原因：", conditionMessage(e)
      )
    },
    finally = {
      devices <- grDevices::dev.list()

      if (
        !is.null(device_id) &&
        !is.null(devices) &&
        device_id %in% unname(devices)
      ) {
        grDevices::dev.off(which = device_id)
      }
    }
  )

  if (!success && file.exists(filename)) {
    unlink(filename)
  }

  invisible(success)
}

safe_filename <- function(x) {
  gsub("[^A-Za-z0-9._-]+", "_", x)
}

classify_group <- function(x) {
  dplyr::case_when(
    x %in% T_SUBTYPES ~ "T",
    x %in% B_SUBTYPES ~ "B",
    TRUE ~ "Other"
  )
}

# ==============================================================================
# 4. 提取通信结果
#    提取失败记为 unavailable，不当作零通信
# ==============================================================================

extract_communication <- function(object, timepoint, pathway = FALSE) {

  tryCatch(
    {
      result <- if (pathway) {
        CellChat::subsetCommunication(
          object,
          slot.name = "netP"
        )
      } else {
        CellChat::subsetCommunication(object)
      }

      result <- as.data.frame(result)

      result$timepoint <- rep(
        timepoint,
        nrow(result)
      )

      list(
        data = result,
        status = "OK"
      )
    },
    error = function(e) {
      record_error(
        paste0(
          timepoint,
          if (pathway) "_pathway_extraction" else "_LR_extraction"
        ),
        conditionMessage(e)
      )

      list(
        data = data.frame(),
        status = "Unavailable"
      )
    }
  )
}

lr_results <- setNames(
  lapply(
    TIME_LEVELS,
    function(tp) {
      extract_communication(
        cellchat_list[[tp]],
        timepoint = tp
      )
    }
  ),
  TIME_LEVELS
)

pathway_results <- setNames(
  lapply(
    TIME_LEVELS,
    function(tp) {
      extract_communication(
        cellchat_list[[tp]],
        timepoint = tp,
        pathway = TRUE
      )
    }
  ),
  TIME_LEVELS
)

extraction_qc <- data.frame(
  Timepoint = TIME_LEVELS,
  LRStatus = vapply(
    lr_results,
    function(x) x$status,
    character(1)
  ),
  PathwayStatus = vapply(
    pathway_results,
    function(x) x$status,
    character(1)
  )
)

write.csv(
  extraction_qc,
  file.path(QC_DIR, "Communication_extraction_status.csv"),
  row.names = FALSE
)

# ==============================================================================
# 5. 各时间点绘图函数
# ==============================================================================

plot_one_timepoint <- function(cellchat_obj, time_value) {

  plot_dir <- file.path(FIG_DIR, time_value)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

  id_levels <- levels(cellchat_obj@idents)

  group_counts <- table(cellchat_obj@idents)
  group_size <- as.numeric(group_counts[id_levels])

  colors_now <- unname(SUBTYPE_COLORS[id_levels])

  # --------------------------------------------------------------------------
  # 5.1 全局 circle 图
  # --------------------------------------------------------------------------

  for (measure in c("count", "weight")) {

    network <- cellchat_obj@net[[measure]]

    # 按网络节点名称严格匹配细胞数和颜色
    network_groups <- rownames(network)

    if (is.null(network_groups)) {
      network_groups <- id_levels
    }

    save_pdf_plot(
      file.path(
        plot_dir,
        paste0(time_value, "_circle_", measure, ".pdf")
      ),
      width = 8,
      height = 8,
      CellChat::netVisual_circle(
        net = network,
        vertex.weight = as.numeric(group_counts[network_groups]),
        color.use = unname(SUBTYPE_COLORS[network_groups]),
        weight.scale = TRUE,
        label.edge = FALSE,
        title.name = paste0(
          time_value,
          if (measure == "count") {
            " | Predicted interaction count"
          } else {
            " | Predicted interaction weight"
          }
        )
      )
    )

    # ------------------------------------------------------------------------
    # 5.2 通信热图
    # ------------------------------------------------------------------------

    save_pdf_plot(
      file.path(
        plot_dir,
        paste0(time_value, "_heatmap_", measure, ".pdf")
      ),
      width = 9,
      height = 8,
      print(
        CellChat::netVisual_heatmap(
          object = cellchat_obj,
          measure = measure,
          color.heatmap = "Reds"
        )
      )
    )
  }

  # --------------------------------------------------------------------------
  # 5.3 信号角色图
  # --------------------------------------------------------------------------

  save_pdf_plot(
    file.path(
      plot_dir,
      paste0(time_value, "_signalingRole_scatter.pdf")
    ),
    width = 7,
    height = 6,
    print(
      CellChat::netAnalysis_signalingRole_scatter(
        object = cellchat_obj,
        color.use = colors_now
      )
    )
  )

  for (pattern_value in c("outgoing", "incoming")) {

    save_pdf_plot(
      file.path(
        plot_dir,
        paste0(
          time_value,
          "_signalingRole_heatmap_",
          pattern_value,
          ".pdf"
        )
      ),
      width = 8,
      height = 9,
      print(
        CellChat::netAnalysis_signalingRole_heatmap(
          object = cellchat_obj,
          pattern = pattern_value
        )
      )
    )
  }

  # --------------------------------------------------------------------------
  # 5.4 T→B 和 B→T 气泡图
  # --------------------------------------------------------------------------

  b_index <- which(id_levels %in% B_SUBTYPES)
  t_index <- which(id_levels %in% T_SUBTYPES)

  if (length(b_index) > 0L && length(t_index) > 0L) {

    save_pdf_plot(
      file.path(
        plot_dir,
        paste0(time_value, "_bubble_B_to_T.pdf")
      ),
      width = 12,
      height = 8,
      print(
        CellChat::netVisual_bubble(
          object = cellchat_obj,
          sources.use = b_index,
          targets.use = t_index,
          remove.isolate = FALSE
        )
      )
    )

    save_pdf_plot(
      file.path(
        plot_dir,
        paste0(time_value, "_bubble_T_to_B.pdf")
      ),
      width = 12,
      height = 8,
      print(
        CellChat::netVisual_bubble(
          object = cellchat_obj,
          sources.use = t_index,
          targets.use = b_index,
          remove.isolate = FALSE
        )
      )
    )
  }

  # --------------------------------------------------------------------------
  # 5.5 前六个通路：贡献图和 circle 图
  # --------------------------------------------------------------------------

  pathway_table <- pathway_results[[time_value]]$data

  required_columns <- c("pathway_name", "prob")

  if (
    nrow(pathway_table) > 0L &&
    all(required_columns %in% colnames(pathway_table))
  ) {

    top_pathways <- pathway_table %>%
      group_by(pathway_name) %>%
      summarise(
        total_prob = sum(prob, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(desc(total_prob)) %>%
      slice_head(n = 6) %>%
      pull(pathway_name)

    for (pathway_name in top_pathways) {

      pathway_filename <- safe_filename(pathway_name)

      save_pdf_plot(
        file.path(
          plot_dir,
          paste0(
            time_value,
            "_contribution_",
            pathway_filename,
            ".pdf"
          )
        ),
        width = 7,
        height = 5,
        print(
          CellChat::netAnalysis_contribution(
            object = cellchat_obj,
            signaling = pathway_name
          )
        )
      )

      save_pdf_plot(
        file.path(
          plot_dir,
          paste0(
            time_value,
            "_aggregate_circle_",
            pathway_filename,
            ".pdf"
          )
        ),
        width = 8,
        height = 8,
        CellChat::netVisual_aggregate(
          object = cellchat_obj,
          signaling = pathway_name,
          layout = "circle",
          color.use = colors_now
        )
      )
    }
  }

  invisible(NULL)
}

# ==============================================================================
# 6. 绘制全部时间点
# ==============================================================================

for (time_value in TIME_LEVELS) {

  message("\n正在绘制：", time_value)

  tryCatch(
    {
      plot_one_timepoint(
        cellchat_obj = cellchat_list[[time_value]],
        time_value = time_value
      )
    },
    error = function(e) {
      record_error(
        paste0(time_value, "_plot_one_timepoint"),
        conditionMessage(e)
      )

      message(
        time_value,
        " 绘图中断，继续下一时间点：",
        conditionMessage(e)
      )
    }
  )
}

# ==============================================================================
# 7. 导出全部时间点通信表
# ==============================================================================

lr_all <- bind_rows(
  lapply(lr_results, function(x) x$data)
)

pathway_all <- bind_rows(
  lapply(pathway_results, function(x) x$data)
)

if (nrow(lr_all) > 0L) {
  lr_all <- lr_all %>%
    mutate(
      source = as.character(source),
      target = as.character(target),
      source_class = classify_group(source),
      target_class = classify_group(target),
      direction = paste0(source_class, " to ", target_class)
    )
}

write.csv(
  lr_all,
  file.path(TABLE_DIR, "All_timepoints_LR_communications.csv"),
  row.names = FALSE
)

write.csv(
  pathway_all,
  file.path(TABLE_DIR, "All_timepoints_pathway_communications.csv"),
  row.names = FALSE
)

# ==============================================================================
# 8. 通信方向动态：B→T、T→B、T→T、B→B
# ==============================================================================

DIRECTION_LEVELS <- c(
  "B to T",
  "T to B",
  "T to T",
  "B to B"
)

direction_observed <- if (nrow(lr_all) > 0L) {

  lr_all %>%
    filter(direction %in% DIRECTION_LEVELS) %>%
    group_by(timepoint, direction) %>%
    summarise(
      strength = sum(prob, na.rm = TRUE),
      n_interactions = n(),
      .groups = "drop"
    )

} else {

  data.frame(
    timepoint = character(),
    direction = character(),
    strength = numeric(),
    n_interactions = integer()
  )
}

direction_summary <- tidyr::crossing(
  timepoint = TIME_LEVELS,
  direction = DIRECTION_LEVELS
) %>%
  left_join(
    direction_observed,
    by = c("timepoint", "direction")
  ) %>%
  left_join(
    extraction_qc %>%
      transmute(
        timepoint = Timepoint,
        ExtractionStatus = LRStatus
      ),
    by = "timepoint"
  )

# 确认相应 T/B 群体实际存在，缺失群体不作为零通信
direction_summary$GroupsPresent <- mapply(
  function(tp, direction) {

    labels <- levels(cellchat_list[[tp]]@idents)

    has_t <- any(labels %in% T_SUBTYPES)
    has_b <- any(labels %in% B_SUBTYPES)

    switch(
      direction,
      "B to T" = has_b && has_t,
      "T to B" = has_t && has_b,
      "T to T" = has_t,
      "B to B" = has_b
    )
  },
  direction_summary$timepoint,
  direction_summary$direction
)

direction_summary <- direction_summary %>%
  mutate(
    Available = ExtractionStatus == "OK" & GroupsPresent,
    strength = ifelse(
      Available,
      coalesce(strength, 0),
      NA_real_
    ),
    n_interactions = ifelse(
      Available,
      coalesce(n_interactions, 0L),
      NA_integer_
    ),
    timepoint = factor(timepoint, levels = TIME_LEVELS),
    direction = factor(direction, levels = DIRECTION_LEVELS)
  )

write.csv(
  direction_summary,
  file.path(TABLE_DIR, "Communication_direction_summary.csv"),
  row.names = FALSE
)

comparison_dir <- file.path(
  FIG_DIR,
  "all_timepoints_comparison"
)

dir.create(
  comparison_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

p_strength <- ggplot(
  direction_summary,
  aes(
    x = timepoint,
    y = strength,
    color = direction,
    group = direction
  )
) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.2) +
  labs(
    x = NULL,
    y = "Summed inferred communication probability",
    color = NULL,
    title = "Inferred communication dynamics"
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )

save_pdf_plot(
  file.path(
    comparison_dir,
    "Communication_direction_dynamics.pdf"
  ),
  width = 8.5,
  height = 4.8,
  print(p_strength)
)

# ==============================================================================
# 9. 合并已有对象
# ==============================================================================

cellchat_merged <- tryCatch(
  {
    CellChat::mergeCellChat(
      object.list = cellchat_list,
      add.names = TIME_LEVELS,
      cell.prefix = TRUE
    )
  },
  error = function(e) {
    record_error("mergeCellChat", conditionMessage(e))
    message("对象合并失败：", conditionMessage(e))
    NULL
  }
)

# ==============================================================================
# 10. 跨时间点比较
# ==============================================================================

if (!is.null(cellchat_merged)) {

  saveRDS(
    cellchat_merged,
    file.path(
      OBJECT_DIR,
      "CellChat_merged_all_timepoints.rds"
    )
  )

  for (measure in c("count", "weight")) {

    save_pdf_plot(
      file.path(
        comparison_dir,
        paste0("compareInteractions_", measure, ".pdf")
      ),
      width = 8,
      height = 6,
      print(
        CellChat::compareInteractions(
          object = cellchat_merged,
          measure = measure,
          show.legend = FALSE,
          group = seq_along(TIME_LEVELS)
        )
      )
    )
  }

  for (stacked_value in c(TRUE, FALSE)) {

    save_pdf_plot(
      file.path(
        comparison_dir,
        if (stacked_value) {
          "rankNet_stacked.pdf"
        } else {
          "rankNet_unstacked.pdf"
        }
      ),
      width = 12,
      height = 8,
      print(
        CellChat::rankNet(
          object = cellchat_merged,
          mode = "comparison",
          stacked = stacked_value,
          do.stat = FALSE
        )
      )
    )
  }
}

# ==============================================================================
# 11. 保存汇总对象与运行记录
# ==============================================================================

saveRDS(
  list(
    Timepoints = TIME_LEVELS,
    SubtypeColors = SUBTYPE_COLORS,
    ObjectAvailability = object_check,
    ExtractionQC = extraction_qc,
    LR = lr_all,
    Pathways = pathway_all,
    DirectionSummary = direction_summary,
    DirectionPlot = p_strength
  ),
  file.path(
    OBJECT_DIR,
    "CellChat_replot_summary.rds"
  )
)

writeLines(
  capture.output(sessionInfo()),
  file.path(QC_DIR, "sessionInfo.txt")
)

plot_errors <- read.csv(
  error_file,
  stringsAsFactors = FALSE
)

if (nrow(plot_errors) == 0L) {

  message("\n绘图完成，没有记录到错误。")

} else {

  message(
    "\n处理完成，共记录 ",
    nrow(plot_errors),
    " 项错误："
  )

  print(plot_errors, row.names = FALSE)
}

message(
  "\n输出目录：",
  normalizePath(REPLOT_DIR, mustWork = FALSE)
)