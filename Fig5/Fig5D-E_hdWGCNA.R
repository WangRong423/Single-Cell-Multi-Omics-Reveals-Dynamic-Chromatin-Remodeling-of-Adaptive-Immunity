#!/usr/bin/env Rscript

library(ArchR)
library(SummarizedExperiment)
library(Seurat)
library(hdWGCNA)
library(WGCNA)
library(Matrix)
library(dplyr)
library(tidyr)
library(tibble)
library(readr)
library(stringr)
library(purrr)
library(ggplot2)
library(patchwork)
library(pheatmap)
library(igraph)
library(gtable)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(GO.db)

# # 从 ArchR 项目准备各细胞类型的 Seurat 对象
# set.seed(12345)
# PROJ_RDS <- "/public/users/wangrong/home/wangrong/results/T_B/TB_final_downstream/proj_TB_only_annotated_with_new_UMAP.rds"
# OUT_DIR <- "/public/users/huhuajie/projects/T_B_CoVarNet/results"
# TARGET_CELLTYPES <- c(
#     "CD4_CD8_DP",
#     "CD4_Memory",
#     "CD4_Memory_ZEB2",
#     "CD4_Naive_CCR4",
#     "CD4_NaiveT",
#     "CD8_Memory",
#     "CD8_Memory_GZMK",
#     "CD8_NaiveT",
#     "gdT",
#     "Intermediate_B",
#     "Memory_B",
#     "Naive_B",
#     "PlasmaB",
#     "Treg"
# )
# MIN_TOTAL_CELLS <- 500
# MIN_GENES_PER_CELL <- 100
# MIN_CELLS_PER_GENE <- 3
# N_VARIABLE_FEATURES <- 3000
# N_PCS <- 30
# proj <- readRDS(PROJ_RDS)
# old_labels <- as.character(getCellColData(proj, select = "TB_Subtype_All")[, 1])
# rename_map <- c(
#     Activated_B = "Intermediate_B",
#     Resting_B = "Memory_B",
#     CD4_Teff = "CD8_Memory_GZMK",
#     Th1 = "CD4_Memory_ZEB2",
#     Th2 = "CD4_Memory",
#     CD4_Tcm = "CD4_Naive_CCR4",
#     CD8_Tem = "CD8_Memory"
# )
# new_labels <- old_labels
# idx <- old_labels %in% names(rename_map)
# new_labels[idx] <- unname(rename_map[old_labels[idx]])
# proj <- addCellColData(
#     ArchRProj = proj,
#     data = new_labels,
#     cells = getCellNames(proj),
#     name = "TB_Subtype_All",
#     force = TRUE
# )
# meta_raw <- getCellColData(proj) %>%
#     as.data.frame() %>%
#     rownames_to_column("cell_name") %>%
#     as_tibble()
# meta <- meta_raw %>%
#     transmute(
#         cell_name = as.character(cell_name),
#         Sample = as.character(Sample),
#         Donor = as.character(Donor),
#         Timepoint_original = as.character(Timepoint),
#         TB_Subtype_All = as.character(TB_Subtype_All)
#     )
# valid_days <- c(0, 1, 3, 6, 14, 30, 31, 33, 36, 44)
# timepoint_order <- paste0("day", valid_days)
# day_text <- gsub("[^0-9]", "", meta$Timepoint_original)
# meta$day_number <- suppressWarnings(as.integer(day_text))
# meta$Timepoint_canonical <- paste0("day", meta$day_number)
# meta <- meta %>%
#     mutate(
#         Timepoint = factor(Timepoint_canonical, levels = timepoint_order),
#         dose = factor(
#             if_else(day_number < 30, "Dose1", "Dose2"),
#             levels = c("Dose1", "Dose2")
#         ),
#         relative_day_number = if_else(
#             day_number < 30,
#             day_number,
#             day_number -
#                 30L
#         ),
#         relative_day = factor(relative_day_number, levels = c(0, 1, 3, 6, 14)),
#         library_id = Sample,
#         bio_sample_id = paste(Donor, as.character(Timepoint), sep = "__")
#     ) %>%
#     dplyr::select(-Timepoint_canonical)
# celltype_counts <- meta %>%
#     count(TB_Subtype_All, name = "n_cells") %>%
#     arrange(desc(n_cells))
# write_tsv(celltype_counts, file.path(OUT_DIR, "celltype_counts.tsv"))
# sample_design <- meta %>%
#     distinct(
#         Sample,
#         Donor,
#         Timepoint,
#         day_number,
#         dose,
#         relative_day,
#         relative_day_number,
#         bio_sample_id
#     ) %>%
#     arrange(Donor, day_number, Sample)
# write_tsv(sample_design, file.path(OUT_DIR, "sample_design.tsv"))
# sample_cell_counts <- meta %>%
#     count(
#         TB_Subtype_All,
#         Donor,
#         Timepoint,
#         day_number,
#         Sample,
#         bio_sample_id,
#         name = "n_cells"
#     ) %>%
#     arrange(TB_Subtype_All, Donor, day_number, Sample)
# write_tsv(
#     sample_cell_counts,
#     file.path(OUT_DIR, "celltype_donor_timepoint_counts.tsv")
# )
# celltypes_to_run <- intersect(TARGET_CELLTYPES, unique(meta$TB_Subtype_All))
# as_dgCMatrix <- function(x) {
#     if (!inherits(x, "dgCMatrix")) {
#         x <- as(x, "dgCMatrix")
#     }
#     x
# }
# deduplicate_gene_symbols <- function(expr, gene_names) {
#     valid <- !is.na(gene_names) & trimws(gene_names) != ""
#     expr <- expr[valid, , drop = FALSE]
#     gene_names <- trimws(gene_names[valid])
#     if (!anyDuplicated(gene_names)) {
#         rownames(expr) <- gene_names
#         return(list(expr = expr, n_duplicates_removed = 0L))
#     }
#     totals <- Matrix::rowSums(expr)
#     order_by_expression <- order(totals, decreasing = TRUE, na.last = TRUE)
#     keep_in_rank_order <- order_by_expression[
#         !duplicated(gene_names[order_by_expression])
#     ]
#     keep_original_order <- sort(keep_in_rank_order)
#     n_removed <- length(gene_names) - length(keep_original_order)
#     expr <- expr[keep_original_order, , drop = FALSE]
#     rownames(expr) <- gene_names[keep_original_order]
#     list(expr = expr, n_duplicates_removed = as.integer(n_removed))
# }
# sparse_log1p <- function(x) {
#     out <- x
#     if (length(out@x) > 0) {
#         out@x <- log1p(out@x)
#     }
#     Matrix::drop0(out)
# }
# sparse_row_variance <- function(x) {
#     n <- ncol(x)
#     means <- Matrix::rowMeans(x)
#     x_squared <- x
#     x_squared@x <- x_squared@x^2
#     mean_squares <- Matrix::rowMeans(x_squared)
#     rm(x_squared)
#     variances <- (mean_squares - means^2) * n / (n - 1)
#     variances[variances < 0 & variances > -1e-10] <- 0
#     names(variances) <- rownames(x)
#     variances
# }
# export_summary <- list()
# for (celltype in celltypes_to_run) {
#     cells_use <- meta %>% filter(TB_Subtype_All == celltype) %>% pull(cell_name)
#     if (length(cells_use) < MIN_TOTAL_CELLS) {
#         next
#     }
#     proj_sub <- subsetCells(ArchRProj = proj, cellNames = cells_use)
#     rna_se <- getMatrixFromProject(
#         ArchRProj = proj_sub,
#         useMatrix = "GeneExpressionMatrix",
#         binarize = FALSE,
#         verbose = TRUE
#     )
#     rna_norm <- as_dgCMatrix(assay(rna_se))
#     feature_data <- as.data.frame(rowData(rna_se))
#     dedup_result <- deduplicate_gene_symbols(
#         expr = rna_norm,
#         gene_names = as.character(feature_data$name)
#     )
#     rna_norm <- dedup_result$expr
#     col_sums_before_filter <- Matrix::colSums(rna_norm)
#     positive_col_sums <- col_sums_before_filter[col_sums_before_filter > 0]
#     median_archr_colsum <- median(positive_col_sums)
#     keep_cells <- Matrix::colSums(rna_norm > 0) >= MIN_GENES_PER_CELL
#     rna_norm <- rna_norm[, keep_cells, drop = FALSE]
#     if (ncol(rna_norm) < MIN_TOTAL_CELLS) {
#         rm(proj_sub, rna_se, rna_norm)
#         gc()
#         next
#     }
#     keep_genes <- Matrix::rowSums(rna_norm > 0) >= MIN_CELLS_PER_GENE
#     rna_norm <- rna_norm[keep_genes, , drop = FALSE]
#     rna_norm <- Matrix::drop0(rna_norm)
#     rna_log <- sparse_log1p(rna_norm)
#     match_index <- match(colnames(rna_norm), meta$cell_name)
#     meta_sub <- as.data.frame(meta[match_index, , drop = FALSE])
#     rownames(meta_sub) <- meta_sub$cell_name
#     meta_sub$cell_name <- NULL
#     safe_name <- gsub("[^A-Za-z0-9]+", "_", celltype)
#     # hdWGCNA 在 Seurat v4 中使用 counts 槽构建宏细胞，因此同时保留非对数矩阵和 log1p 矩阵
#     seu <- CreateSeuratObject(
#         counts = rna_norm,
#         assay = "RNA",
#         meta.data = meta_sub,
#         project = safe_name,
#         min.cells = 0,
#         min.features = 0
#     )
#     seu <- SetAssayData(
#         object = seu,
#         assay = "RNA",
#         slot = "data",
#         new.data = rna_log
#     )
#     DefaultAssay(seu) <- "RNA"
#     Key(seu[["RNA"]]) <- "rna_"
#     Idents(seu) <- "TB_Subtype_All"
#     gene_variance <- sparse_row_variance(rna_log)
#     valid_variance <- is.finite(gene_variance) & gene_variance > 0
#     n_variable_use <- min(N_VARIABLE_FEATURES, sum(valid_variance))
#     variable_genes <- names(sort(
#         gene_variance[valid_variance],
#         decreasing = TRUE
#     ))[seq_len(n_variable_use)]
#     VariableFeatures(seu) <- variable_genes
#     seu <- ScaleData(
#         seu,
#         assay = "RNA",
#         features = variable_genes,
#         verbose = FALSE
#     )
#     n_pcs_use <- min(N_PCS, length(variable_genes) - 1L, ncol(seu) - 1L)
#     seu <- RunPCA(
#         seu,
#         assay = "RNA",
#         features = variable_genes,
#         npcs = n_pcs_use,
#         verbose = FALSE,
#         seed.use = 12345
#     )
#     seu@misc$input_matrix <- "ArchR GeneExpressionMatrix"
#     seu@misc$counts_slot_semantics <- paste0(
#         "ArchR non-log normalized expression (NormCounts), not raw UMI counts; ",
#         "median pre-filter column sum=",
#         signif(median_archr_colsum, 6)
#     )
#     seu@misc$data_slot_semantics <- "natural-log log1p of RNA counts slot"
#     seu@misc$celltype <- celltype
#     seu@misc$preparation_script <- "02_prepare_hdWGCNA_rds_fixed.R"
#     seu@misc$preparation_date <- as.character(Sys.Date())
#     celltype_sample_counts <- seu@meta.data %>%
#         rownames_to_column("cell_name") %>%
#         count(
#             Donor,
#             Timepoint,
#             day_number,
#             dose,
#             relative_day,
#             library_id,
#             bio_sample_id,
#             name = "n_cells"
#         ) %>%
#         arrange(Donor, day_number, library_id)
#     count_file <- file.path(
#         OUT_DIR,
#         paste0(safe_name, "_sample_cell_counts.tsv")
#     )
#     write_tsv(celltype_sample_counts, count_file)
#     output_file <- file.path(OUT_DIR, paste0(safe_name, "_RNA_for_hdWGCNA.rds"))
#     saveRDS(seu, output_file, compress = FALSE)
#     export_summary[[celltype]] <- data.frame(
#         celltype = celltype,
#         n_cells = ncol(seu),
#         n_genes = nrow(seu),
#         n_variable_genes = length(VariableFeatures(seu)),
#         n_pcs = n_pcs_use,
#         median_archr_colsum_before_filter = median_archr_colsum,
#         n_duplicate_symbols_removed = dedup_result$n_duplicates_removed,
#         output_file = output_file,
#         stringsAsFactors = FALSE
#     )
#     rm(proj_sub, rna_se, rna_norm, rna_log, meta_sub, seu, gene_variance)
#     gc()
# }
# if (length(export_summary) > 0) {
#     export_summary_df <- bind_rows(export_summary)
#     summary_file <- file.path(OUT_DIR, "hdWGCNA_Seurat_export_summary.tsv")
#     write_tsv(export_summary_df, summary_file)
#     print(export_summary_df, row.names = FALSE)
# } else {
#     NULL
# }


# 构建各细胞类型的 hdWGCNA 网络
set.seed(12345)
options(stringsAsFactors = FALSE)
INPUT_DIR <- "/public/users/huhuajie/projects/T_B_CoVarNet/results"
RESULT_ROOT <- "/public/users/huhuajie/projects/codeCheck/T_B_scMultiomics/Fig5/hdWGCNA_results"
dir.create(RESULT_ROOT, recursive = TRUE, showWarnings = FALSE)
N_THREADS <- 8
GENE_SELECT_FRACTION <- 0.05
MIN_NETWORK_GENES <- 200
DEFAULT_K <- 25
MIN_K <- 10
MIN_CELLS_PER_BIOSAMPLE <- 25
MAX_SHARED <- 10
# 所有生物样本合计的目标宏细胞数，运行时按合格样本数分摊
TARGET_METACELLS <- 1000
MIN_BIOLOGICAL_SAMPLES <- 3
MIN_DATEXPR_SAMPLES <- 20
SOFT_POWERS <- c(1:10, seq(12, 30, by = 2))
SOFT_POWER_BLOCK_SIZE <- 1000
MIN_SOFT_POWER <- 3
PRIMARY_SFT_CUTOFF <- 0.8
FALLBACK_SFT_CUTOFF <- 0.7
NETWORK_TYPE <- "signed"
TOM_TYPE <- "signed"
CORRELATION_FUNCTION <- "bicor"
MIN_MODULE_SIZE <- 50
DEEP_SPLIT <- 4
MERGE_CUT_HEIGHT <- 0.2
PAM_STAGE <- FALSE
MAX_BLOCK_SIZE <- 30000
OVERWRITE_TOM <- FALSE
# 仅在存在独立技术批次时设置；不要使用可能与生物学时间点混杂的 Sample 或 library_id
HME_BATCH_COLUMN <- NULL
N_HUBS <- 20
WGCNA::enableWGCNAThreads(nThreads = N_THREADS)
rds_files <- list.files(
    INPUT_DIR,
    pattern = "_RNA_for_hdWGCNA\\.rds$",
    full.names = TRUE
)
input_stubs <- sub("_RNA_for_hdWGCNA\\.rds$", "", basename(rds_files))
completed_rds <- file.path(RESULT_ROOT, input_stubs, paste0(input_stubs, "_hdWGCNA_final.rds"))
rds_files <- rds_files[!file.exists(completed_rds)]
safe_file_name <- function(x) {
    gsub("[^A-Za-z0-9]+", "_", x)
}
coerce_metadata_vector <- function(x, column_name) {
    if (inherits(x, "Rle")) {
        x <- as.vector(x)
    }
    if (is.data.frame(x)) {
        x <- x[[1L]]
    }
    if (is.list(x)) {
        element_lengths <- lengths(x)
        x <- vapply(
            x,
            FUN = function(z) {
                if (is.null(z) || length(z) == 0L || all(is.na(z))) {
                    return(NA_character_)
                }
                as.character(z[[1L]])
            },
            FUN.VALUE = character(1L),
            USE.NAMES = FALSE
        )
    } else {
        x <- as.vector(x)
    }
    as.character(x)
}
standardize_required_metadata <- function(seu, required_metadata) {
    metadata <- as.data.frame(seu@meta.data, stringsAsFactors = FALSE)
    for (column_name in required_metadata) {
        metadata[[column_name]] <- coerce_metadata_vector(
            metadata[[column_name]],
            column_name = column_name
        )
    }
    seu@meta.data <- metadata
    seu
}
# 优先选择满足无标度拟合阈值且斜率为负的最低软阈值，未达阈值时使用最优拟合回退规则
select_soft_power_auto <- function(
    power_table,
    min_power = 3,
    primary_cutoff = 0.8,
    fallback_cutoff = 0.7
) {
    pt <- power_table %>%
        dplyr::mutate(
            Power = as.numeric(Power),
            SFT.R.sq = as.numeric(SFT.R.sq),
            slope = as.numeric(slope)
        ) %>%
        dplyr::filter(
            is.finite(.data$Power),
            is.finite(.data$SFT.R.sq),
            is.finite(.data$slope),
            .data$Power >= min_power
        ) %>%
        dplyr::arrange(.data$Power)
    negative_slope <- pt %>% dplyr::filter(.data$slope < 0)
    if (nrow(negative_slope) == 0) {
        best_r2 <- max(pt$SFT.R.sq, na.rm = TRUE)
        chosen <- pt %>%
            dplyr::filter(.data$SFT.R.sq >= best_r2 - 0.05) %>%
            dplyr::slice_min(.data$Power, n = 1, with_ties = FALSE)
        return(list(
            power = chosen$Power[[1]],
            method = "fallback_no_negative_slope_within_0.05_of_max_R2",
            cutoff_reached = FALSE,
            selected_row = chosen
        ))
    }
    primary <- negative_slope %>%
        dplyr::filter(.data$SFT.R.sq >= primary_cutoff) %>%
        dplyr::slice_min(.data$Power, n = 1, with_ties = FALSE)
    if (nrow(primary) == 1) {
        return(list(
            power = primary$Power[[1]],
            method = paste0(
                "lowest_power_R2_ge_",
                primary_cutoff,
                "_negative_slope"
            ),
            cutoff_reached = TRUE,
            selected_row = primary
        ))
    }
    fallback_1 <- negative_slope %>%
        dplyr::filter(.data$SFT.R.sq >= fallback_cutoff) %>%
        dplyr::slice_min(.data$Power, n = 1, with_ties = FALSE)
    if (nrow(fallback_1) == 1) {
        return(list(
            power = fallback_1$Power[[1]],
            method = paste0(
                "fallback_lowest_power_R2_ge_",
                fallback_cutoff,
                "_negative_slope"
            ),
            cutoff_reached = FALSE,
            selected_row = fallback_1
        ))
    }
    max_r2 <- max(negative_slope$SFT.R.sq, na.rm = TRUE)
    fallback_2 <- negative_slope %>%
        dplyr::filter(.data$SFT.R.sq >= max_r2 - 0.05) %>%
        dplyr::slice_min(.data$Power, n = 1, with_ties = FALSE)
    list(
        power = fallback_2$Power[[1]],
        method = "fallback_negative_slope_within_0.05_of_max_R2",
        cutoff_reached = FALSE,
        selected_row = fallback_2
    )
}
standardize_seurat_input <- function(seu) {
    required_metadata <- c(
        "TB_Subtype_All",
        "bio_sample_id",
        "Donor",
        "Timepoint",
        "dose",
        "relative_day",
        "library_id"
    )
    seu <- standardize_required_metadata(
        seu,
        required_metadata = required_metadata
    )
    seu
}
save_plot_list_pdf <- function(plot_list, output_file) {
    pdf(output_file, width = 14, height = 10, onefile = TRUE)
    on.exit(dev.off(), add = TRUE)
    if (inherits(plot_list, "ggplot")) {
        print(plot_list)
    } else if (is.list(plot_list)) {
        print(patchwork::wrap_plots(plot_list, ncol = 2))
    } else {
        print(plot_list)
    }
}
run_summary <- list()
run_errors <- list()
for (input_file in rds_files) {
    input_stub <- sub("_RNA_for_hdWGCNA\\.rds$", "", basename(input_file))
    current_stage <- "initialize"
    one_result <- tryCatch(
        {
            current_stage <- "readRDS"
            seu <- readRDS(input_file)
            current_stage <- "standardize_input"
            seu <- standardize_seurat_input(seu)
            current_stage <- "set_default_assay"
            Seurat::DefaultAssay(seu) <- "RNA"
            current_stage <- "extract_celltype"
            celltype_values <- unique(seu@meta.data[["TB_Subtype_All"]])
            celltype_values <- celltype_values[
                !is.na(celltype_values) & celltype_values != ""
            ]
            celltype <- celltype_values[[1]]
            safe_celltype <- safe_file_name(celltype)
            wgcna_name <- paste0(safe_celltype, "_hdWGCNA")
            celltype_dir <- file.path(RESULT_ROOT, safe_celltype)
            tom_dir <- file.path(celltype_dir, "TOM")
            dir.create(celltype_dir, recursive = TRUE, showWarnings = FALSE)
            dir.create(tom_dir, recursive = TRUE, showWarnings = FALSE)
            current_stage <- "count_biological_samples"
            biosample_counts <- seu@meta.data %>%
                tibble::rownames_to_column("cell_name") %>%
                dplyr::count(
                    .data$bio_sample_id,
                    .data$Donor,
                    .data$Timepoint,
                    name = "n_cells"
                ) %>%
                dplyr::arrange(dplyr::desc(.data$n_cells))
            readr::write_tsv(
                biosample_counts,
                file.path(celltype_dir, "biological_sample_cell_counts_all.tsv")
            )
            eligible_biosamples <- biosample_counts %>%
                dplyr::filter(.data$n_cells >= MIN_CELLS_PER_BIOSAMPLE)
            excluded_biosamples <- biosample_counts %>%
                dplyr::filter(.data$n_cells < MIN_CELLS_PER_BIOSAMPLE)
            readr::write_tsv(
                eligible_biosamples,
                file.path(celltype_dir, "biological_samples_used.tsv")
            )
            readr::write_tsv(
                excluded_biosamples,
                file.path(celltype_dir, "biological_samples_excluded.tsv")
            )
            current_stage <- "subset_eligible_biological_samples"
            biosample_vector <- seu@meta.data[["bio_sample_id"]]
            cells_keep <- rownames(seu@meta.data)[
                biosample_vector %in% eligible_biosamples[["bio_sample_id"]]
            ]
            seu <- subset(seu, cells = cells_keep)
            Seurat::Idents(seu) <- "TB_Subtype_All"
            min_group_cells <- min(eligible_biosamples$n_cells)
            k_use <- min(DEFAULT_K, floor(min_group_cells / 2))
            k_use <- max(MIN_K, k_use)
            if (k_use >= min_group_cells) {
                k_use <- min_group_cells - 1L
            }
            max_shared_use <- min(MAX_SHARED, max(1L, floor(k_use / 2)))
            target_metacells_use <- max(20L, ceiling(TARGET_METACELLS / nrow(eligible_biosamples)))
            pca_dims_available <- ncol(Embeddings(seu, reduction = "pca"))
            dims_use <- seq_len(min(30L, pca_dims_available))
            current_stage <- "SetupForWGCNA"
            seu <- SetupForWGCNA(
                seu,
                gene_select = "fraction",
                fraction = GENE_SELECT_FRACTION,
                wgcna_name = wgcna_name
            )
            network_genes <- GetWGCNAGenes(seu, wgcna_name = wgcna_name)
            readr::write_lines(
                network_genes,
                file.path(celltype_dir, "network_genes.txt")
            )
            current_stage <- "MetacellsByGroups"
            seu <- MetacellsByGroups(
                seurat_obj = seu,
                group.by = c("TB_Subtype_All", "bio_sample_id"),
                ident.group = "TB_Subtype_All",
                k = k_use,
                reduction = "pca",
                dims = dims_use,
                assay = "RNA",
                slot = "counts",
                mode = "average",
                min_cells = MIN_CELLS_PER_BIOSAMPLE,
                max_shared = max_shared_use,
                target_metacells = target_metacells_use,
                verbose = TRUE,
                wgcna_name = wgcna_name
            )
            current_stage <- "NormalizeMetacells"
            seu <- NormalizeMetacells(
                seu,
                wgcna_name = wgcna_name,
                normalization.method = "LogNormalize",
                scale.factor = 10000,
                verbose = FALSE
            )
            metacell_obj <- GetMetacellObject(seu, wgcna_name = wgcna_name)
            metacell_metadata <- metacell_obj@meta.data %>%
                tibble::rownames_to_column("metacell")
            readr::write_tsv(
                metacell_metadata,
                file.path(celltype_dir, "metacell_metadata.tsv")
            )
            current_stage <- "SetDatExpr"
            seu <- SetDatExpr(
                seu,
                group_name = celltype,
                use_metacells = TRUE,
                group.by = "TB_Subtype_All",
                assay = "RNA",
                slot = "data",
                wgcna_name = wgcna_name
            )
            dat_expr <- GetDatExpr(seu, wgcna_name = wgcna_name)
            # hdWGCNA 0.4.11 已在内部传递 verbose 参数，此处不重复设置
            current_stage <- "TestSoftPowers"
            seu <- TestSoftPowers(
                seu,
                powers = SOFT_POWERS,
                networkType = NETWORK_TYPE,
                corFnc = CORRELATION_FUNCTION,
                blockSize = SOFT_POWER_BLOCK_SIZE,
                wgcna_name = wgcna_name
            )
            power_table <- GetPowerTable(seu, wgcna_name = wgcna_name)
            readr::write_tsv(
                as.data.frame(power_table),
                file.path(celltype_dir, "soft_power_test.tsv")
            )
            soft_selection <- select_soft_power_auto(
                power_table = power_table,
                min_power = MIN_SOFT_POWER,
                primary_cutoff = PRIMARY_SFT_CUTOFF,
                fallback_cutoff = FALLBACK_SFT_CUTOFF
            )
            selected_power <- as.numeric(soft_selection$power)
            selected_power_row <- as.data.frame(soft_selection$selected_row)
            selected_power_row$selection_method <- soft_selection$method
            selected_power_row$primary_cutoff_reached <- soft_selection$cutoff_reached
            readr::write_tsv(
                selected_power_row,
                file.path(celltype_dir, "selected_soft_power.tsv")
            )
            soft_power_plots <- PlotSoftPowers(
                seu,
                selected_power = selected_power,
                wgcna_name = wgcna_name
            )
            save_plot_list_pdf(
                soft_power_plots,
                file.path(celltype_dir, "soft_power_diagnostics.pdf")
            )
            checkpoint_file <- file.path(
                celltype_dir,
                paste0(safe_celltype, "_before_ConstructNetwork.rds")
            )
            saveRDS(seu, checkpoint_file, compress = FALSE)
            current_stage <- "ConstructNetwork"
            seu <- ConstructNetwork(
                seu,
                soft_power = selected_power,
                min_power = MIN_SOFT_POWER,
                tom_outdir = tom_dir,
                tom_name = safe_celltype,
                overwrite_tom = OVERWRITE_TOM,
                wgcna_name = wgcna_name,
                maxBlockSize = MAX_BLOCK_SIZE,
                randomSeed = 12345,
                corType = CORRELATION_FUNCTION,
                networkType = NETWORK_TYPE,
                TOMType = TOM_TYPE,
                deepSplit = DEEP_SPLIT,
                pamStage = PAM_STAGE,
                minModuleSize = MIN_MODULE_SIZE,
                mergeCutHeight = MERGE_CUT_HEIGHT
            )
            pdf(
                file.path(celltype_dir, "module_dendrogram.pdf"),
                width = 12,
                height = 8
            )
            PlotDendrogram(
                seu,
                main = paste0(celltype, " hdWGCNA dendrogram"),
                wgcna_name = wgcna_name
            )
            dev.off()
            use_harmony <- !is.null(HME_BATCH_COLUMN) &&
                HME_BATCH_COLUMN %in% colnames(seu@meta.data) &&
                length(unique(as.character(seu@meta.data[[
                    HME_BATCH_COLUMN
                ]]))) >=
                    2L
            if (use_harmony) {
                current_stage <- "ModuleEigengenes_harmonized"
                seu <- ModuleEigengenes(
                    seu,
                    group.by.vars = HME_BATCH_COLUMN,
                    assay = "RNA",
                    wgcna_name = wgcna_name
                )
            } else {
                current_stage <- "ModuleEigengenes_unharmonized"
                seu <- ModuleEigengenes(
                    seu,
                    group.by.vars = NULL,
                    assay = "RNA",
                    wgcna_name = wgcna_name
                )
            }
            current_stage <- "ModuleConnectivity"
            seu <- ModuleConnectivity(
                seu,
                group.by = "TB_Subtype_All",
                group_name = celltype,
                corFnc = CORRELATION_FUNCTION,
                harmonized = use_harmony,
                assay = "RNA",
                slot = "data",
                sparse = TRUE,
                reassign_modules = TRUE,
                TOM_use = wgcna_name,
                wgcna_name = wgcna_name
            )
            modules_all <- GetModules(seu, wgcna_name = wgcna_name)
            modules_non_grey <- modules_all %>%
                dplyr::filter(.data$module != "grey")
            readr::write_tsv(
                as.data.frame(modules_all),
                file.path(celltype_dir, "modules_all.tsv")
            )
            readr::write_tsv(
                as.data.frame(modules_non_grey),
                file.path(celltype_dir, "modules_non_grey.tsv")
            )
            hub_genes <- GetHubGenes(
                seu,
                n_hubs = N_HUBS,
                mods = as.character(unique(modules_non_grey$module)),
                wgcna_name = wgcna_name
            )
            readr::write_tsv(
                as.data.frame(hub_genes),
                file.path(
                    celltype_dir,
                    paste0("top_", N_HUBS, "_hub_genes.tsv")
                )
            )
            MEs <- GetMEs(seu, harmonized = FALSE, wgcna_name = wgcna_name)
            readr::write_tsv(
                as.data.frame(MEs) %>% tibble::rownames_to_column("cell_name"),
                file.path(celltype_dir, "module_eigengenes_MEs.tsv")
            )
            if (use_harmony) {
                hMEs <- GetMEs(seu, harmonized = TRUE, wgcna_name = wgcna_name)
                readr::write_tsv(
                    as.data.frame(hMEs) %>%
                        tibble::rownames_to_column("cell_name"),
                    file.path(
                        celltype_dir,
                        "harmonized_module_eigengenes_hMEs.tsv"
                    )
                )
            }
            final_rds <- file.path(
                celltype_dir,
                paste0(safe_celltype, "_hdWGCNA_final.rds")
            )
            saveRDS(seu, final_rds, compress = FALSE)
            module_count <- length(unique(modules_non_grey$module))
            summary_row <- data.frame(
                celltype = celltype,
                input_file = input_file,
                cells_used = ncol(seu),
                biological_samples_used = nrow(eligible_biosamples),
                metacells = ncol(metacell_obj),
                network_genes_selected = length(network_genes),
                datExpr_samples = nrow(dat_expr),
                datExpr_genes = ncol(dat_expr),
                metacell_k = k_use,
                target_metacells_per_sample = target_metacells_use,
                selected_soft_power = selected_power,
                soft_power_method = soft_selection$method,
                primary_SFT_cutoff_reached = soft_selection$cutoff_reached,
                non_grey_modules = module_count,
                hME_harmony_applied = use_harmony,
                final_rds = final_rds,
                stringsAsFactors = FALSE
            )
            summary_row
        },
        error = function(e) {
            error_text <- paste0(
                "[stage: ",
                current_stage,
                "] ",
                conditionMessage(e)
            )
            message(error_text)
            error_df <- data.frame(
                input_file = input_file,
                input_stub = input_stub,
                stage = current_stage,
                error = error_text,
                stringsAsFactors = FALSE
            )
            readr::write_tsv(
                error_df,
                file.path(RESULT_ROOT, paste0(input_stub, "_ERROR.tsv"))
            )
            run_errors[[input_stub]] <<- error_df
        }
    )
    if (!is.null(one_result)) {
        run_summary[[input_stub]] <- one_result
    }
    gc()
}
if (length(run_summary) > 0) {
    run_summary_df <- dplyr::bind_rows(run_summary)
    readr::write_tsv(
        run_summary_df,
        file.path(RESULT_ROOT, "hdWGCNA_run_summary.tsv")
    )
    print(run_summary_df, row.names = FALSE)
}
if (length(run_errors) > 0) {
    run_errors_df <- dplyr::bind_rows(run_errors)
    readr::write_tsv(
        run_errors_df,
        file.path(RESULT_ROOT, "hdWGCNA_run_errors.tsv")
    )
}


# 下游枢纽基因、热图及 GO 富集分析
Sys.setenv(
    OMP_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1",
    NUMEXPR_NUM_THREADS = "1"
)
options(stringsAsFactors = FALSE)
set.seed(12345)
WGCNA::disableWGCNAThreads()
WGCNA::enableWGCNAThreads(nThreads = 8)
RESULT_ROOT <- "/public/users/huhuajie/projects/codeCheck/T_B_scMultiomics/Fig5/hdWGCNA_results"
FINAL_RDS_PATTERN <- "_hdWGCNA_final\\.rds$"
TARGET_CELLTYPES <- NULL
OVERWRITE_OUTPUT <- TRUE
CELLTYPE_COL <- "TB_Subtype_All"
BIOSAMPLE_COL <- "bio_sample_id"
DONOR_COL <- "Donor"
TIMEPOINT_COL <- "Timepoint"
DOSE_COL <- "dose"
RELATIVE_DAY_COL <- "relative_day"
LIBRARY_COL <- "library_id"
USE_HARMONIZED_ME <- FALSE
MODULE_LABEL_SEPARATOR <- "_"
HEATMAP_ROW_Z_LEGEND_TITLE <- "Row Z-score"
HEATMAP_RAW_ME_LEGEND_TITLE <- "Mean module eigengene"
HUB_N_GENES_PER_MODULE <- 25
HUB_NETWORK_TOP_EDGES <- 80
HUB_NETWORK_MIN_EDGES_FALLBACK <- 20
HUB_NETWORK_MIN_SIGNED_ADJACENCY <- 0.05
HUB_NETWORK_LABEL_TOP_N <- 25
HUB_NETWORK_LAYOUT_SEED <- 12345
HUB_HEATMAP_GENES_PER_MODULE <- 10
# 分别保存完整 GO 富集结果和聚焦 B/T 细胞生物过程的筛选结果
GO_ONTOLOGIES <- c("BP")
GO_PADJ_CUTOFF <- 0.05
GO_MIN_GS_SIZE <- 1
GO_MAX_GS_SIZE <- 1e+05
GO_MAX_SIGNIFICANT_TERMS_PER_MODULE <- Inf
GO_OVERVIEW_TOP_N_PER_MODULE_ONTOLOGY <- 10
BT_GO_PADJ_CUTOFF <- 0.05
BT_GO_MIN_COUNT <- 3
BT_GO_MIN_FOLD_ENRICHMENT <- 1.3
BT_GO_MAX_BG_COUNT <- 1000
BT_GO_TOP_N_PER_CATEGORY <- 8
BT_GO_OVERVIEW_TOP_N_PER_MODULE <- 20
# 从下列免疫相关根节点递归纳入子条目，priority 较小者优先归类
BT_GO_ROOTS <- tibble::tribble(
    ~root_go_id         , ~immune_category              , ~priority                      , "GO:0050853"             , "B-cell receptor signaling" ,                                               1 , "GO:0042113"                       ,
    "B-cell activation" ,                             2 , "GO:0030183"                   , "B-cell differentiation" ,                           3 , "GO:0016064"                                    , "Immunoglobulin-mediated immunity" ,
                      4 , "GO:0002377"                  , "Immunoglobulin production"    ,                        5 , "GO:0050852"                , "T-cell receptor signaling"                     ,                                  6 , "GO:0042110"                  , "T-cell activation"                   ,
                      7 , "GO:0030217"                  , "T-cell differentiation"       ,                        8 , "GO:0050851"                , "Antigen-receptor signaling"                    ,                                  9 , "GO:0019882"                  , "Antigen processing and presentation" ,
                     10 , "GO:0046649"                  , "Lymphocyte activation"        ,                       11 , "GO:0030098"                , "Lymphocyte differentiation"                    ,                                 12 , "GO:0002449"                  , "Lymphocyte-mediated immunity"        ,
                     13 , "GO:0002250"                  , "Adaptive immune response"     ,                       14 , "GO:0002460"                , "Somatic-recombination-based adaptive immunity" ,                                 15 ,
    "GO:0019221"        , "Cytokine-mediated signaling" ,                             16 , "GO:0001816"             , "Cytokine production"       ,                                              17 , "GO:0060337"                       , "Type-I interferon signaling" ,
                     18 , "GO:0060333"                  , "Interferon-gamma signaling"   ,                       19 , "GO:0070098"                , "Chemokine-mediated signaling"                  ,                                 20 , "GO:0050900"                  , "Leukocyte migration"                 ,
                     21 , "GO:0007159"                  , "Leukocyte cell-cell adhesion" ,                       22 , "GO:0001906"                , "Cell killing and cytotoxicity"                 ,                                 23
)
BT_GO_EXCLUDE_IDS <- c(
    "GO:0008150",
    "GO:0065007",
    "GO:0050789",
    "GO:0050794",
    "GO:0009987",
    "GO:0050896",
    "GO:0048583",
    "GO:0048584",
    "GO:0048585",
    "GO:0002376",
    "GO:0006955",
    "GO:0002682",
    "GO:0050776",
    "GO:0002684",
    "GO:0002683",
    "GO:0006952"
)
BT_GO_EXCLUDE_DESCRIPTION_REGEX <- paste0(
    "^(",
    paste(
        c(
            "biological process",
            "biological regulation",
            "cellular process",
            "regulation of biological process",
            "regulation of cellular process",
            "response to stimulus",
            "immune system process",
            "immune response",
            "regulation of immune system process",
            "regulation of immune response",
            "positive regulation of immune system process",
            "negative regulation of immune system process",
            "defense response"
        ),
        collapse = "|"
    ),
    ")$"
)
message2 <- function(...) NULL
safe_name <- function(x) {
    x <- as.character(x)
    x <- gsub("[[:space:]]+", "_", x)
    x <- gsub("[/:;|,()\\[\\]\\{\\}]+", "_", x, perl = TRUE)
    x <- gsub("[^A-Za-z0-9_.+-]+", "_", x)
    x <- gsub("_+", "_", x)
    x <- gsub("^_|_$", "", x)
    x
}
first_existing_column <- function(
    df,
    exact_candidates,
    regex_candidates = character()
) {
    cn <- colnames(df)
    for (candidate in exact_candidates) {
        hit <- cn[tolower(cn) == tolower(candidate)]
        if (length(hit) > 0) {
            return(hit[[1]])
        }
    }
    for (candidate in regex_candidates) {
        hit <- cn[grepl(candidate, cn, ignore.case = TRUE)]
        if (length(hit) > 0) {
            return(hit[[1]])
        }
    }
    NA_character_
}
valid_r_color <- function(x, fallback = "skyblue2") {
    x <- as.character(x)[1]
    ok <- tryCatch(
        {
            grDevices::col2rgb(x)
            TRUE
        },
        error = function(e) FALSE
    )
    if (isTRUE(ok)) {
        x
    } else {
        fallback
    }
}
row_zscore <- function(mat) {
    mat <- as.matrix(mat)
    out <- t(scale(t(mat)))
    out[!is.finite(out)] <- 0
    out
}
parse_gene_ratio <- function(x) {
    x <- as.character(x)
    out <- rep(NA_real_, length(x))
    frac_idx <- grepl("^[0-9.]+/[0-9.]+$", x)
    if (any(frac_idx)) {
        pieces <- strsplit(x[frac_idx], "/", fixed = TRUE)
        out[frac_idx] <- vapply(
            pieces,
            function(z) as.numeric(z[[1]]) / as.numeric(z[[2]]),
            numeric(1)
        )
    }
    numeric_idx <- !frac_idx & grepl("^[0-9.]+$", x)
    out[numeric_idx] <- suppressWarnings(as.numeric(x[numeric_idx]))
    out
}
extract_day_number <- function(x) {
    suppressWarnings(as.integer(gsub("[^0-9]", "", as.character(x))))
}
clean_go_term_name <- function(x) {
    x <- as.character(x)
    x <- stringr::str_replace_all(x, "\\s*[\\[(]GO:[0-9]{7}[\\])]\\s*$", "")
    x <- stringr::str_replace_all(x, "\\s*GO:[0-9]{7}\\s*$", "")
    stringr::str_squish(x)
}
module_number_value <- function(x) {
    suppressWarnings(as.integer(gsub("[^0-9]", "", as.character(x))))
}
make_output_dirs <- function(celltype_dir) {
    downstream_dir <- file.path(celltype_dir, "downstream_hub_heatmap_GO")
    dirs <- list(
        root = downstream_dir,
        network = file.path(downstream_dir, "hub_gene_network"),
        heatmap = file.path(downstream_dir, "heatmap"),
        go = file.path(downstream_dir, "GO_enrichment_all_terms"),
        tables = file.path(downstream_dir, "tables")
    )
    invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
    dirs
}
get_celltype <- function(seu, final_rds) {
    if (CELLTYPE_COL %in% colnames(seu@meta.data)) {
        values <- unique(as.character(seu@meta.data[[CELLTYPE_COL]]))
        values <- values[!is.na(values) & values != ""]
        if (length(values) == 1) {
            return(values[[1]])
        }
    }
    safe_name(basename(dirname(final_rds)))
}
get_wgcna_name_safe <- function(seu, celltype) {
    preferred <- paste0(safe_name(celltype), "_hdWGCNA")
    is_valid_wgcna_name <- function(x) {
        if (length(x) != 1L || is.na(x) || x == "") {
            return(FALSE)
        }
        tryCatch(
            {
                mods <- hdWGCNA::GetModules(seu, wgcna_name = x)
                !is.null(mods) && nrow(as.data.frame(mods)) > 0
            },
            error = function(e) FALSE
        )
    }
    if (is_valid_wgcna_name(preferred)) {
        return(preferred)
    }
    active_name <- tryCatch(hdWGCNA::GetActiveWGCNA(seu), error = function(e) {
        NULL
    })
    active_name <- as.character(active_name)[1]
    if (
        !is.na(active_name) &&
            nzchar(active_name) &&
            is_valid_wgcna_name(active_name)
    ) {
        return(active_name)
    }
    active_misc <- tryCatch(seu@misc$active_wgcna, error = function(e) NULL)
    active_misc <- as.character(active_misc)[1]
    if (
        !is.na(active_misc) &&
            nzchar(active_misc) &&
            is_valid_wgcna_name(active_misc)
    ) {
        return(active_misc)
    }
    misc_names <- tryCatch(names(seu@misc), error = function(e) character())
    misc_names <- setdiff(misc_names, c("active_wgcna", "hdWGCNA"))
    misc_names <- misc_names[!is.na(misc_names) & misc_names != ""]
    valid_top_level <- misc_names[vapply(
        misc_names,
        is_valid_wgcna_name,
        logical(1)
    )]
    preferred_hits <- valid_top_level[
        valid_top_level == preferred |
            grepl(
                paste0("^", safe_name(celltype), "(_hdWGCNA)?$"),
                valid_top_level
            )
    ]
    if (length(preferred_hits) == 1L) {
        return(preferred_hits[[1]])
    }
    suffix_hits <- valid_top_level[grepl("_hdWGCNA$", valid_top_level)]
    if (length(suffix_hits) == 1L) {
        return(suffix_hits[[1]])
    }
    if (length(valid_top_level) == 1L) {
        return(valid_top_level[[1]])
    }
    nested_names <- tryCatch(names(seu@misc$hdWGCNA), error = function(e) {
        character()
    })
    nested_names <- nested_names[!is.na(nested_names) & nested_names != ""]
    valid_nested <- nested_names[vapply(
        nested_names,
        is_valid_wgcna_name,
        logical(1)
    )]
    if (preferred %in% valid_nested) {
        return(preferred)
    }
    if (length(valid_nested) == 1L) {
        return(valid_nested[[1]])
    }
    all_valid <- unique(c(valid_top_level, valid_nested))
}
get_wgcna_misc_safe <- function(seu, wgcna_name) {
    direct <- tryCatch(seu@misc[[wgcna_name]], error = function(e) NULL)
    if (!is.null(direct)) {
        return(direct)
    }
    tryCatch(seu@misc$hdWGCNA[[wgcna_name]], error = function(e) NULL)
}
get_selected_soft_power <- function(celltype_dir, seu, wgcna_name) {
    selected_file <- file.path(celltype_dir, "selected_soft_power.tsv")
    if (file.exists(selected_file)) {
        df <- readr::read_tsv(selected_file, show_col_types = FALSE)
        power_col <- first_existing_column(
            df,
            exact_candidates = c(
                "Power",
                "power",
                "soft_power",
                "selected_power"
            ),
            regex_candidates = c("power")
        )
        if (!is.na(power_col)) {
            power_now <- suppressWarnings(as.numeric(df[[power_col]][[1]]))
            if (is.finite(power_now)) {
                return(power_now)
            }
        }
    }
    misc_now <- get_wgcna_misc_safe(seu, wgcna_name)
    if (!is.null(misc_now)) {
        candidate_names <- c("soft_power", "softPower", "power")
        for (candidate in candidate_names) {
            value <- suppressWarnings(as.numeric(misc_now[[candidate]][[1]]))
            if (length(value) == 1 && is.finite(value)) {
                return(value)
            }
        }
    }
}
get_dat_expr_safe <- function(seu, wgcna_name) {
    dat_expr <- tryCatch(
        hdWGCNA::GetDatExpr(seu, wgcna_name = wgcna_name),
        error = function(e) NULL
    )
    if (!is.null(dat_expr)) {
        dat_expr <- as.data.frame(dat_expr)
        if (nrow(dat_expr) > 0 && ncol(dat_expr) > 0) {
            return(dat_expr)
        }
    }
    misc_now <- get_wgcna_misc_safe(seu, wgcna_name)
    if (is.null(misc_now)) {
        return(NULL)
    }
    for (candidate in c("datExpr", "dat_expr", "data_expr", "wgcna_datExpr")) {
        if (!is.null(misc_now[[candidate]])) {
            out <- as.data.frame(misc_now[[candidate]])
            if (nrow(out) > 0 && ncol(out) > 0) {
                return(out)
            }
        }
    }
}
standardize_modules <- function(modules_df) {
    df <- tibble::as_tibble(modules_df)
    gene_col <- first_existing_column(
        df,
        exact_candidates = c("gene_name", "gene", "symbol"),
        regex_candidates = c("gene")
    )
    module_col <- first_existing_column(
        df,
        exact_candidates = c("module"),
        regex_candidates = c("module")
    )
    color_col <- first_existing_column(
        df,
        exact_candidates = c("color", "module_color"),
        regex_candidates = c("color")
    )
    out0 <- tibble::tibble(
        gene_name = as.character(df[[gene_col]]),
        module = as.character(df[[module_col]]),
        color = if (!is.na(color_col)) {
            as.character(df[[color_col]])
        } else {
            as.character(df[[module_col]])
        }
    ) %>%
        dplyr::filter(!is.na(.data$gene_name), .data$gene_name != "") %>%
        dplyr::filter(!is.na(.data$module), .data$module != "") %>%
        dplyr::distinct(.data$module, .data$gene_name, .keep_all = TRUE)
    module_map <- out0 %>%
        dplyr::filter(tolower(.data$module) != "grey") %>%
        dplyr::group_by(.data$module) %>%
        dplyr::summarise(
            color = dplyr::first(.data$color),
            module_size = dplyr::n_distinct(.data$gene_name),
            .groups = "drop"
        ) %>%
        dplyr::arrange(dplyr::desc(.data$module_size), .data$module) %>%
        dplyr::mutate(
            module_number = dplyr::row_number(),
            module_id = paste0("M", .data$module_number),
            module_label = paste0(
                .data$module,
                MODULE_LABEL_SEPARATOR,
                .data$module_id
            )
        )
    grey_map <- out0 %>%
        dplyr::filter(tolower(.data$module) == "grey") %>%
        dplyr::group_by(.data$module) %>%
        dplyr::summarise(
            color = dplyr::first(.data$color),
            module_size = dplyr::n_distinct(.data$gene_name),
            .groups = "drop"
        ) %>%
        dplyr::mutate(
            module_number = NA_integer_,
            module_id = NA_character_,
            module_label = "grey"
        )
    module_map <- dplyr::bind_rows(module_map, grey_map)
    out <- out0 %>%
        dplyr::select(-.data$color) %>%
        dplyr::left_join(module_map, by = "module") %>%
        dplyr::arrange(
            is.na(.data$module_number),
            .data$module_number,
            .data$gene_name
        )
    out
}
get_module_mapping <- function(modules_std) {
    modules_std %>%
        dplyr::select(
            module,
            module_label,
            module_id,
            module_number,
            module_size,
            color
        ) %>%
        dplyr::distinct() %>%
        dplyr::arrange(
            is.na(.data$module_number),
            .data$module_number,
            .data$module
        )
}
standardize_hubs <- function(hub_df, modules_std) {
    df <- tibble::as_tibble(hub_df)
    gene_col <- first_existing_column(
        df,
        exact_candidates = c("gene_name", "gene", "symbol"),
        regex_candidates = c("gene")
    )
    module_col <- first_existing_column(
        df,
        exact_candidates = c("module"),
        regex_candidates = c("module")
    )
    kme_col <- first_existing_column(
        df,
        exact_candidates = c("kME", "kme", "kMEmax", "kme_max"),
        regex_candidates = c("kME")
    )
    module_map <- get_module_mapping(modules_std)
    out <- tibble::tibble(
        gene_name = as.character(df[[gene_col]]),
        module = as.character(df[[module_col]]),
        kME = if (!is.na(kme_col)) {
            suppressWarnings(as.numeric(df[[kme_col]]))
        } else {
            NA_real_
        }
    ) %>%
        dplyr::left_join(module_map, by = "module") %>%
        dplyr::filter(!is.na(.data$gene_name), .data$gene_name != "") %>%
        dplyr::filter(
            !is.na(.data$module),
            .data$module != "",
            tolower(.data$module) != "grey"
        ) %>%
        dplyr::group_by(
            .data$module,
            .data$module_label,
            .data$module_number
        ) %>%
        dplyr::arrange(dplyr::desc(abs(.data$kME)), .by_group = TRUE) %>%
        dplyr::ungroup() %>%
        dplyr::distinct(.data$module, .data$gene_name, .keep_all = TRUE) %>%
        dplyr::arrange(.data$module_number, dplyr::desc(abs(.data$kME)))
    out
}
add_pheatmap_legend_title <- function(pheatmap_result, legend_title) {
    if (
        is.null(legend_title) ||
            length(legend_title) == 0 ||
            is.na(legend_title[[1]]) ||
            legend_title[[1]] == ""
    ) {
        return(pheatmap_result)
    }
    legend_index <- which(pheatmap_result$gtable$layout$name == "legend")
    if (length(legend_index) != 1) {
        return(pheatmap_result)
    }
    legend_grob <- pheatmap_result$gtable$grobs[[legend_index]]
    title_grob <- grid::textGrob(
        label = legend_title[[1]],
        x = grid::unit(0, "npc"),
        just = "left",
        gp = grid::gpar(fontsize = 10, fontface = "bold")
    )
    legend_box <- gtable::gtable(
        widths = grid::unit.pmax(
            grid::grobWidth(title_grob),
            grid::grobWidth(legend_grob)
        ),
        heights = grid::unit.c(
            grid::grobHeight(title_grob) +
                grid::unit(2, "mm"),
            grid::grobHeight(legend_grob)
        )
    )
    legend_box <- gtable::gtable_add_grob(
        legend_box,
        title_grob,
        t = 1,
        l = 1,
        clip = "off"
    )
    legend_box <- gtable::gtable_add_grob(
        legend_box,
        legend_grob,
        t = 2,
        l = 1,
        clip = "off"
    )
    pheatmap_result$gtable$grobs[[legend_index]] <- legend_box
    pheatmap_result
}
draw_pheatmap_with_title <- function(mat, legend_title = NULL, ...) {
    ph <- pheatmap::pheatmap(mat, ..., silent = TRUE)
    ph <- add_pheatmap_legend_title(ph, legend_title)
    grid::grid.newpage()
    grid::grid.draw(ph$gtable)
    invisible(ph)
}
save_pheatmap_both <- function(
    mat,
    pdf_file,
    png_file,
    width,
    height,
    legend_title = NULL,
    ...
) {
    if (nrow(mat) == 0 || ncol(mat) == 0) {
        return(invisible(NULL))
    }
    grDevices::pdf(pdf_file, width = width, height = height)
    draw_pheatmap_with_title(mat = mat, legend_title = legend_title, ...)
    grDevices::dev.off()
    grDevices::png(
        png_file,
        width = width,
        height = height,
        units = "in",
        res = 300
    )
    draw_pheatmap_with_title(mat = mat, legend_title = legend_title, ...)
    grDevices::dev.off()
}
plot_hub_networks <- function(
    seu,
    modules_std,
    hubs_std,
    dat_expr,
    selected_power,
    outdir,
    prefix
) {
    dat_expr <- as.data.frame(dat_expr)
    gene_available <- colnames(dat_expr)
    module_order <- hubs_std %>%
        dplyr::select(module, module_label, module_id, module_number) %>%
        dplyr::distinct() %>%
        dplyr::arrange(.data$module_number)
    all_nodes <- list()
    all_edges <- list()
    combined_pdf <- file.path(
        outdir,
        paste0(prefix, "__all_module_hub_networks.pdf")
    )
    grDevices::pdf(combined_pdf, width = 9, height = 9, onefile = TRUE)
    on.exit(
        {
            if (grDevices::dev.cur() != 1) {
                grDevices::dev.off()
            }
        },
        add = TRUE
    )
    for (i in seq_len(nrow(module_order))) {
        module_now <- module_order$module[[i]]
        module_label_now <- module_order$module_label[[i]]
        hubs_module <- hubs_std %>%
            dplyr::filter(.data$module == module_now) %>%
            dplyr::arrange(dplyr::desc(abs(.data$kME))) %>%
            dplyr::slice_head(n = HUB_N_GENES_PER_MODULE)
        genes_use <- intersect(hubs_module$gene_name, gene_available)
        if (length(genes_use) < 3) {
            message2(
                "Skip hub network for ",
                prefix,
                " ",
                module_label_now,
                ": fewer than 3 genes in datExpr."
            )
            next
        }
        expr_module <- as.matrix(dat_expr[, genes_use, drop = FALSE])
        cor_matrix <- suppressWarnings(WGCNA::bicor(
            expr_module,
            use = "pairwise.complete.obs",
            maxPOutliers = 0.1
        ))
        cor_matrix[!is.finite(cor_matrix)] <- 0
        diag(cor_matrix) <- 0
        adjacency_matrix <- ((1 + cor_matrix) / 2)^selected_power
        adjacency_matrix[!is.finite(adjacency_matrix)] <- 0
        diag(adjacency_matrix) <- 0
        upper_idx <- which(upper.tri(adjacency_matrix), arr.ind = TRUE)
        edge_candidates <- tibble::tibble(
            from = rownames(adjacency_matrix)[upper_idx[, 1]],
            to = colnames(adjacency_matrix)[upper_idx[,
                2
            ]],
            correlation = cor_matrix[upper_idx],
            signed_adjacency = adjacency_matrix[upper_idx]
        ) %>%
            dplyr::arrange(dplyr::desc(.data$signed_adjacency))
        edge_df <- edge_candidates %>%
            dplyr::filter(
                .data$signed_adjacency >= HUB_NETWORK_MIN_SIGNED_ADJACENCY
            ) %>%
            dplyr::slice_head(n = HUB_NETWORK_TOP_EDGES)
        minimum_edges <- min(
            HUB_NETWORK_MIN_EDGES_FALLBACK,
            nrow(edge_candidates)
        )
        if (nrow(edge_df) < minimum_edges) {
            edge_df <- edge_candidates %>%
                dplyr::slice_head(
                    n = min(
                        HUB_NETWORK_TOP_EDGES,
                        max(minimum_edges, nrow(edge_df))
                    )
                )
        }
        edge_df <- edge_df %>%
            dplyr::mutate(
                module_original = module_now,
                module = module_label_now,
                sign = ifelse(.data$correlation >= 0, "positive", "negative")
            )
        node_df <- hubs_module %>%
            dplyr::filter(.data$gene_name %in% genes_use) %>%
            dplyr::arrange(dplyr::desc(abs(.data$kME))) %>%
            dplyr::mutate(
                module_original = .data$module,
                module = .data$module_label,
                label_rank = dplyr::row_number(),
                label = ifelse(
                    .data$label_rank <= HUB_NETWORK_LABEL_TOP_N,
                    .data$gene_name,
                    ""
                ),
                abs_kME = abs(.data$kME)
            )
        graph_now <- igraph::graph_from_data_frame(
            d = edge_df %>%
                dplyr::select(from, to, correlation, signed_adjacency, sign),
            vertices = node_df %>% dplyr::rename(name = gene_name),
            directed = FALSE
        )
        set.seed(HUB_NETWORK_LAYOUT_SEED)
        layout_now <- igraph::layout_with_fr(
            graph_now,
            weights = igraph::E(graph_now)$signed_adjacency
        )
        vertex_kme <- igraph::V(graph_now)$abs_kME
        vertex_kme[!is.finite(vertex_kme)] <- 0
        if (length(unique(vertex_kme)) <= 1) {
            vertex_size <- rep(11, length(vertex_kme))
        } else {
            vertex_size <- 7 +
                9 *
                    (vertex_kme - min(vertex_kme)) /
                    (max(vertex_kme) - min(vertex_kme))
        }
        edge_adj <- igraph::E(graph_now)$signed_adjacency
        if (length(unique(edge_adj)) <= 1) {
            edge_width <- rep(1.5, length(edge_adj))
        } else {
            edge_width <- 0.4 +
                4 * (edge_adj - min(edge_adj)) / (max(edge_adj) - min(edge_adj))
        }
        module_color <- valid_r_color(unique(node_df$color)[[1]])
        edge_color <- ifelse(
            igraph::E(graph_now)$correlation >= 0,
            "grey35",
            "grey75"
        )
        edge_lty <- ifelse(igraph::E(graph_now)$correlation >= 0, 1, 2)
        plot_network <- function() {
            plot(
                graph_now,
                layout = layout_now,
                vertex.size = vertex_size,
                vertex.color = module_color,
                vertex.frame.color = "grey25",
                vertex.label = igraph::V(graph_now)$label,
                vertex.label.cex = 0.75,
                vertex.label.color = "black",
                edge.width = edge_width,
                edge.color = edge_color,
                edge.lty = edge_lty,
                main = paste0(
                    prefix,
                    " | ",
                    module_label_now,
                    "\nsoft power = ",
                    selected_power,
                    "; nodes = ",
                    igraph::vcount(graph_now),
                    "; edges = ",
                    igraph::ecount(graph_now)
                )
            )
        }
        plot_network()
        module_stub <- safe_name(module_label_now)
        grDevices::pdf(
            file.path(
                outdir,
                paste0(prefix, "__", module_stub, "__hub_network.pdf")
            ),
            width = 9,
            height = 9
        )
        plot_network()
        grDevices::dev.off()
        grDevices::png(
            file.path(
                outdir,
                paste0(prefix, "__", module_stub, "__hub_network.png")
            ),
            width = 2700,
            height = 2700,
            res = 300
        )
        plot_network()
        grDevices::dev.off()
        readr::write_tsv(
            edge_df,
            file.path(
                outdir,
                paste0(prefix, "__", module_stub, "__hub_network_edges.tsv")
            )
        )
        readr::write_tsv(
            node_df,
            file.path(
                outdir,
                paste0(prefix, "__", module_stub, "__hub_network_nodes.tsv")
            )
        )
        all_edges[[module_label_now]] <- edge_df
        all_nodes[[module_label_now]] <- node_df
    }
    grDevices::dev.off()
    nodes_all <- dplyr::bind_rows(all_nodes)
    edges_all <- dplyr::bind_rows(all_edges)
    if (nrow(nodes_all) > 0) {
        readr::write_tsv(
            nodes_all,
            file.path(
                outdir,
                paste0(prefix, "__all_module_hub_network_nodes.tsv")
            )
        )
    }
    if (nrow(edges_all) > 0) {
        readr::write_tsv(
            edges_all,
            file.path(
                outdir,
                paste0(prefix, "__all_module_hub_network_edges.tsv")
            )
        )
    }
    list(nodes = nodes_all, edges = edges_all, combined_pdf = combined_pdf)
}
get_me_matrix_safe <- function(seu, wgcna_name, harmonized = FALSE) {
    me <- tryCatch(
        hdWGCNA::GetMEs(seu, harmonized = harmonized, wgcna_name = wgcna_name),
        error = function(e) NULL
    )
    if (is.null(me) && harmonized) {
        message2("Harmonized MEs unavailable; falling back to raw MEs.")
        me <- tryCatch(
            hdWGCNA::GetMEs(seu, harmonized = FALSE, wgcna_name = wgcna_name),
            error = function(e) NULL
        )
    }
    if (is.null(me)) {
        return(NULL)
    }
    as.data.frame(me)
}
strip_me_prefix <- function(x) {
    sub("^(h?ME)[_.-]?", "", as.character(x), ignore.case = TRUE)
}
build_me_module_map <- function(me_cols, modules_std) {
    module_map <- get_module_mapping(modules_std) %>%
        dplyr::filter(tolower(.data$module) != "grey") %>%
        dplyr::mutate(module_key = tolower(.data$module))
    out <- tibble::tibble(
        me_col = as.character(me_cols),
        module_guess = strip_me_prefix(me_cols),
        module_key = tolower(strip_me_prefix(me_cols))
    ) %>%
        dplyr::left_join(module_map, by = "module_key") %>%
        dplyr::mutate(
            module = dplyr::if_else(
                is.na(.data$module) | .data$module == "",
                .data$module_guess,
                .data$module
            ),
            module_label = dplyr::if_else(
                is.na(.data$module_label) | .data$module_label == "",
                .data$me_col,
                .data$module_label
            )
        ) %>%
        dplyr::filter(
            !grepl("grey", .data$me_col, ignore.case = TRUE),
            !grepl("grey", .data$module, ignore.case = TRUE)
        )
    duplicated_labels <- out$module_label[duplicated(out$module_label)]
    out %>%
        dplyr::arrange(
            is.na(.data$module_number),
            .data$module_number,
            .data$me_col
        )
}
rename_me_columns_for_output <- function(df, module_col_map) {
    out <- as.data.frame(df)
    match_index <- match(module_col_map$me_col, colnames(out))
    valid <- !is.na(match_index)
    colnames(out)[match_index[valid]] <- module_col_map$module_label[valid]
    out
}
prepare_sample_me_summary <- function(seu, me_df, modules_std) {
    common_cells <- intersect(rownames(me_df), rownames(seu@meta.data))
    me_df <- me_df[common_cells, , drop = FALSE]
    meta_df <- seu@meta.data[common_cells, , drop = FALSE]
    required_meta <- c(BIOSAMPLE_COL, DONOR_COL, TIMEPOINT_COL)
    missing_meta <- setdiff(required_meta, colnames(meta_df))
    optional_meta <- intersect(
        c(DOSE_COL, RELATIVE_DAY_COL, LIBRARY_COL),
        colnames(meta_df)
    )
    grouping_optional_meta <- intersect(
        c(DOSE_COL, RELATIVE_DAY_COL),
        colnames(meta_df)
    )
    module_col_map <- build_me_module_map(colnames(me_df), modules_std)
    module_cols <- module_col_map$me_col
    me_meta <- dplyr::bind_cols(
        meta_df %>%
            tibble::rownames_to_column("cell_name") %>%
            dplyr::select(
                cell_name,
                dplyr::all_of(required_meta),
                dplyr::all_of(optional_meta)
            ),
        me_df[, module_cols, drop = FALSE] %>% tibble::as_tibble()
    )
    grouping_cols <- unique(c(required_meta, grouping_optional_meta))
    sample_summary <- me_meta %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(grouping_cols))) %>%
        dplyr::summarise(
            dplyr::across(dplyr::all_of(module_cols), ~ mean(.x, na.rm = TRUE)),
            n_cells = dplyr::n(),
            .groups = "drop"
        )
    if (LIBRARY_COL %in% colnames(me_meta)) {
        library_summary <- me_meta %>%
            dplyr::group_by(.data[[BIOSAMPLE_COL]]) %>%
            dplyr::summarise(
                library_ids = paste(
                    sort(unique(as.character(.data[[LIBRARY_COL]]))),
                    collapse = ","
                ),
                .groups = "drop"
            )
        sample_summary <- sample_summary %>%
            dplyr::left_join(library_summary, by = BIOSAMPLE_COL)
    }
    sample_summary <- sample_summary %>%
        dplyr::mutate(
            .timepoint_day = extract_day_number(.data[[TIMEPOINT_COL]])
        ) %>%
        dplyr::arrange(
            .data$.timepoint_day,
            .data[[DONOR_COL]],
            .data[[BIOSAMPLE_COL]]
        )
    list(
        me_meta = me_meta,
        me_meta_labeled = rename_me_columns_for_output(me_meta, module_col_map),
        sample_summary = sample_summary,
        sample_summary_labeled = rename_me_columns_for_output(
            sample_summary,
            module_col_map
        ),
        module_cols = module_cols,
        module_col_map = module_col_map,
        grouping_cols = grouping_cols
    )
}
plot_module_activity_heatmaps <- function(
    sample_me,
    module_cols,
    module_col_map,
    outdir,
    prefix
) {
    sample_me <- as.data.frame(sample_me)
    sample_matrix <- t(as.matrix(sample_me[, module_cols, drop = FALSE]))
    rownames(sample_matrix) <- module_col_map$module_label[match(
        rownames(sample_matrix),
        module_col_map$me_col
    )]
    colnames(sample_matrix) <- as.character(sample_me[[BIOSAMPLE_COL]])
    sample_matrix_z <- row_zscore(sample_matrix)
    annotation_cols <- intersect(
        c(DONOR_COL, TIMEPOINT_COL, DOSE_COL, RELATIVE_DAY_COL),
        colnames(sample_me)
    )
    annotation_col <- sample_me %>%
        dplyr::select(dplyr::all_of(c(BIOSAMPLE_COL, annotation_cols))) %>%
        dplyr::distinct(.data[[BIOSAMPLE_COL]], .keep_all = TRUE) %>%
        tibble::column_to_rownames(BIOSAMPLE_COL)
    annotation_col <- annotation_col[colnames(sample_matrix_z), , drop = FALSE]
    sample_output <- rename_me_columns_for_output(sample_me, module_col_map)
    readr::write_tsv(
        sample_output %>% dplyr::select(-dplyr::any_of(".timepoint_day")),
        file.path(outdir, paste0(prefix, "__module_activity_by_biosample.tsv"))
    )
    sample_width <- max(9, 3.5 + 0.2 * ncol(sample_matrix_z))
    sample_height <- max(6, 3.5 + 0.28 * nrow(sample_matrix_z))
    save_pheatmap_both(
        mat = sample_matrix_z,
        pdf_file = file.path(
            outdir,
            paste0(prefix, "__module_activity_biosample_rowZ_heatmap.pdf")
        ),
        png_file = file.path(
            outdir,
            paste0(prefix, "__module_activity_biosample_rowZ_heatmap.png")
        ),
        width = sample_width,
        height = sample_height,
        legend_title = HEATMAP_ROW_Z_LEGEND_TITLE,
        annotation_col = annotation_col,
        cluster_rows = TRUE,
        cluster_cols = FALSE,
        border_color = NA,
        show_colnames = TRUE,
        fontsize_col = 7,
        angle_col = 45,
        main = paste0(prefix, " module activity by biological sample"),
        na_col = "grey90"
    )
    timepoint_summary <- sample_me %>%
        dplyr::group_by(.data[[TIMEPOINT_COL]]) %>%
        dplyr::summarise(
            dplyr::across(dplyr::all_of(module_cols), ~ mean(.x, na.rm = TRUE)),
            n_biological_samples = dplyr::n(),
            .groups = "drop"
        ) %>%
        dplyr::mutate(
            .timepoint_day = extract_day_number(.data[[TIMEPOINT_COL]])
        ) %>%
        dplyr::arrange(.data$.timepoint_day)
    timepoint_output <- rename_me_columns_for_output(
        timepoint_summary,
        module_col_map
    )
    readr::write_tsv(
        timepoint_output %>% dplyr::select(-dplyr::any_of(".timepoint_day")),
        file.path(outdir, paste0(prefix, "__module_activity_by_timepoint.tsv"))
    )
    time_matrix <- t(as.matrix(timepoint_summary[, module_cols, drop = FALSE]))
    rownames(time_matrix) <- module_col_map$module_label[match(
        rownames(time_matrix),
        module_col_map$me_col
    )]
    colnames(time_matrix) <- as.character(timepoint_summary[[TIMEPOINT_COL]])
    time_matrix_z <- row_zscore(time_matrix)
    time_width <- max(8, 3.5 + 0.65 * ncol(time_matrix_z))
    time_height <- max(6, 3.5 + 0.3 * nrow(time_matrix_z))
    save_pheatmap_both(
        mat = time_matrix_z,
        pdf_file = file.path(
            outdir,
            paste0(prefix, "__module_activity_timepoint_rowZ_heatmap.pdf")
        ),
        png_file = file.path(
            outdir,
            paste0(prefix, "__module_activity_timepoint_rowZ_heatmap.png")
        ),
        width = time_width,
        height = time_height,
        legend_title = HEATMAP_ROW_Z_LEGEND_TITLE,
        cluster_rows = TRUE,
        cluster_cols = FALSE,
        border_color = NA,
        show_colnames = TRUE,
        fontsize_col = 9,
        angle_col = 45,
        main = paste0(prefix, " mean module activity by timepoint"),
        na_col = "grey90"
    )
    save_pheatmap_both(
        mat = time_matrix,
        pdf_file = file.path(
            outdir,
            paste0(prefix, "__module_activity_timepoint_raw_heatmap.pdf")
        ),
        png_file = file.path(
            outdir,
            paste0(prefix, "__module_activity_timepoint_raw_heatmap.png")
        ),
        width = time_width,
        height = time_height,
        legend_title = HEATMAP_RAW_ME_LEGEND_TITLE,
        cluster_rows = TRUE,
        cluster_cols = FALSE,
        border_color = NA,
        show_colnames = TRUE,
        fontsize_col = 9,
        angle_col = 45,
        main = paste0(prefix, " mean module eigengene by timepoint"),
        na_col = "grey90"
    )
    list(
        sample_matrix = sample_matrix,
        sample_matrix_z = sample_matrix_z,
        timepoint_summary = timepoint_summary,
        time_matrix = time_matrix,
        time_matrix_z = time_matrix_z
    )
}
compute_group_average_expression <- function(
    seu,
    genes,
    group_col,
    sample_order = NULL,
    assay = "RNA",
    slot = "data"
) {
    genes <- unique(intersect(as.character(genes), rownames(seu[[assay]])))
    group_vector <- as.character(seu@meta.data[[group_col]])
    valid_cells <- !is.na(group_vector) & group_vector != ""
    cells_use <- colnames(seu)[valid_cells]
    group_vector <- group_vector[valid_cells]
    if (is.null(sample_order)) {
        sample_order <- unique(group_vector)
    } else {
        sample_order <- sample_order[sample_order %in% group_vector]
    }
    group_factor <- factor(group_vector, levels = sample_order)
    design <- Matrix::sparse.model.matrix(~ 0 + group_factor)
    colnames(design) <- levels(group_factor)
    expression_matrix <- Seurat::GetAssayData(seu, assay = assay, slot = slot)[
        genes,
        cells_use,
        drop = FALSE
    ]
    summed <- expression_matrix %*% design
    group_n <- Matrix::colSums(design)
    averaged <- sweep(as.matrix(summed), 2, group_n, "/")
    averaged[!is.finite(averaged)] <- NA_real_
    averaged
}
plot_hub_gene_expression_heatmap <- function(
    seu,
    hubs_std,
    sample_me,
    outdir,
    prefix
) {
    hubs_heat <- hubs_std %>%
        dplyr::group_by(
            .data$module,
            .data$module_label,
            .data$module_number
        ) %>%
        dplyr::arrange(dplyr::desc(abs(.data$kME)), .by_group = TRUE) %>%
        dplyr::slice_head(n = HUB_HEATMAP_GENES_PER_MODULE) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(.data$module_number, dplyr::desc(abs(.data$kME))) %>%
        dplyr::distinct(.data$gene_name, .keep_all = TRUE)
    sample_order <- as.character(sample_me[[BIOSAMPLE_COL]])
    average_expression <- compute_group_average_expression(
        seu = seu,
        genes = hubs_heat$gene_name,
        group_col = BIOSAMPLE_COL,
        sample_order = sample_order,
        assay = "RNA",
        slot = "data"
    )
    ordered_genes <- hubs_heat$gene_name[
        hubs_heat$gene_name %in% rownames(average_expression)
    ]
    average_expression <- average_expression[ordered_genes, , drop = FALSE]
    expression_z <- row_zscore(average_expression)
    hubs_heat <- hubs_heat %>%
        dplyr::filter(.data$gene_name %in% rownames(expression_z)) %>%
        dplyr::mutate(
            gene_name = factor(.data$gene_name, levels = rownames(expression_z))
        ) %>%
        dplyr::arrange(.data$gene_name)
    row_annotation <- hubs_heat %>%
        dplyr::transmute(
            gene_name = as.character(.data$gene_name),
            Module = .data$module_label,
            kME = .data$kME
        ) %>%
        tibble::column_to_rownames("gene_name")
    annotation_cols <- intersect(
        c(DONOR_COL, TIMEPOINT_COL, DOSE_COL, RELATIVE_DAY_COL),
        colnames(sample_me)
    )
    column_annotation <- sample_me %>%
        dplyr::select(dplyr::all_of(c(BIOSAMPLE_COL, annotation_cols))) %>%
        dplyr::distinct(.data[[BIOSAMPLE_COL]], .keep_all = TRUE) %>%
        tibble::column_to_rownames(BIOSAMPLE_COL)
    column_annotation <- column_annotation[
        colnames(expression_z),
        ,
        drop = FALSE
    ]
    if (nrow(hubs_heat) > 1) {
        gap_rows <- which(
            hubs_heat$module_label[-nrow(hubs_heat)] !=
                hubs_heat$module_label[-1]
        )
    } else {
        gap_rows <- integer()
    }
    expression_table <- as.data.frame(average_expression) %>%
        tibble::rownames_to_column("gene_name") %>%
        dplyr::left_join(
            hubs_heat %>%
                dplyr::mutate(gene_name = as.character(.data$gene_name)) %>%
                dplyr::select(
                    gene_name,
                    module_original = module,
                    module = module_label,
                    module_id,
                    module_number,
                    kME
                ),
            by = "gene_name"
        ) %>%
        dplyr::relocate(
            module_original,
            module,
            module_id,
            module_number,
            kME,
            .after = gene_name
        )
    readr::write_tsv(
        expression_table,
        file.path(
            outdir,
            paste0(prefix, "__top_hub_gene_biosample_mean_expression.tsv")
        )
    )
    heat_width <- max(10, 4 + 0.2 * ncol(expression_z))
    heat_height <- max(8, 4 + 0.2 * nrow(expression_z))
    save_pheatmap_both(
        mat = expression_z,
        pdf_file = file.path(
            outdir,
            paste0(
                prefix,
                "__top_hub_gene_expression_biosample_rowZ_heatmap.pdf"
            )
        ),
        png_file = file.path(
            outdir,
            paste0(
                prefix,
                "__top_hub_gene_expression_biosample_rowZ_heatmap.png"
            )
        ),
        width = heat_width,
        height = heat_height,
        legend_title = HEATMAP_ROW_Z_LEGEND_TITLE,
        annotation_col = column_annotation,
        annotation_row = row_annotation,
        cluster_rows = FALSE,
        cluster_cols = FALSE,
        gaps_row = gap_rows,
        border_color = NA,
        show_rownames = TRUE,
        show_colnames = TRUE,
        fontsize_row = 7,
        fontsize_col = 7,
        angle_col = 45,
        main = paste0(prefix, " top hub-gene expression by biological sample"),
        na_col = "grey90"
    )
    list(
        hubs = hubs_heat,
        average_expression = average_expression,
        expression_z = expression_z
    )
}
map_symbols_to_entrez <- function(symbols) {
    symbols <- unique(as.character(symbols))
    symbols <- symbols[!is.na(symbols) & symbols != ""]
    if (length(symbols) == 0) {
        return(tibble::tibble(SYMBOL = character(), ENTREZID = character()))
    }
    mapped <- suppressMessages(AnnotationDbi::select(
        org.Hs.eg.db,
        keys = symbols,
        keytype = "SYMBOL",
        columns = c("SYMBOL", "ENTREZID")
    )) %>%
        tibble::as_tibble() %>%
        dplyr::filter(!is.na(.data$SYMBOL), .data$SYMBOL != "") %>%
        dplyr::filter(!is.na(.data$ENTREZID), .data$ENTREZID != "") %>%
        dplyr::distinct(.data$SYMBOL, .data$ENTREZID)
    mapped
}
get_go_bp_descendants <- function(root_go_id) {
    offspring <- tryCatch(
        {
            value <- mget(root_go_id, GO.db::GOBPOFFSPRING, ifnotfound = NA)[[
                1
            ]]
            as.character(value)
        },
        error = function(e) character()
    )
    unique(c(root_go_id, offspring[!is.na(offspring) & offspring != ""]))
}
build_bt_go_catalog <- function(go_term_names) {
    catalog_long <- purrr::pmap_dfr(
        BT_GO_ROOTS,
        function(root_go_id, immune_category, priority) {
            tibble::tibble(
                ID = get_go_bp_descendants(root_go_id),
                root_go_id = root_go_id,
                immune_category = immune_category,
                priority = as.integer(priority)
            )
        }
    )
    catalog <- catalog_long %>%
        dplyr::inner_join(
            go_term_names %>% dplyr::select(ID, Description),
            by = "ID"
        ) %>%
        dplyr::mutate(Description = clean_go_term_name(.data$Description)) %>%
        dplyr::arrange(.data$priority, .data$ID) %>%
        dplyr::group_by(.data$ID) %>%
        dplyr::slice_head(n = 1) %>%
        dplyr::ungroup()
    catalog
}
select_bt_relevant_go <- function(go_df, bt_catalog) {
    if (nrow(go_df) == 0 || nrow(bt_catalog) == 0) {
        return(tibble::tibble())
    }
    go_df %>%
        dplyr::filter(
            .data$ONTOLOGY == "BP",
            .data$ID %in% bt_catalog$ID,
            !.data$ID %in% BT_GO_EXCLUDE_IDS,
            !is.na(.data$p.adjust),
            .data$p.adjust <= BT_GO_PADJ_CUTOFF,
            .data$Count >= BT_GO_MIN_COUNT,
            .data$FoldEnrichment >= BT_GO_MIN_FOLD_ENRICHMENT,
            .data$BgCount <= BT_GO_MAX_BG_COUNT
        ) %>%
        dplyr::left_join(
            bt_catalog %>%
                dplyr::select(ID, immune_category, root_go_id, priority),
            by = "ID"
        ) %>%
        dplyr::mutate(Description = clean_go_term_name(.data$Description)) %>%
        dplyr::filter(
            !grepl(
                BT_GO_EXCLUDE_DESCRIPTION_REGEX,
                .data$Description,
                ignore.case = TRUE
            )
        ) %>%
        dplyr::mutate(
            GeneRatio_numeric = parse_gene_ratio(.data$GeneRatio),
            neg_log10_padj = -log10(.data$p.adjust + 1e-300),
            specificity_score = .data$neg_log10_padj *
                log2(pmax(.data$FoldEnrichment, 1))
        ) %>%
        dplyr::arrange(
            .data$priority,
            dplyr::desc(.data$specificity_score),
            .data$p.adjust,
            dplyr::desc(.data$Count)
        )
}
plot_bt_go_one_module <- function(bt_df, outdir, prefix, module_label_now) {
    if (nrow(bt_df) == 0) {
        return(invisible(NULL))
    }
    plot_df <- bt_df %>%
        dplyr::group_by(.data$immune_category) %>%
        dplyr::arrange(
            dplyr::desc(.data$specificity_score),
            .data$p.adjust,
            .by_group = TRUE
        ) %>%
        dplyr::slice_head(n = BT_GO_TOP_N_PER_CATEGORY) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(
            term_label = stringr::str_wrap(
                clean_go_term_name(.data$Description),
                width = 55
            ),
            term_plot_id = paste0(
                .data$immune_category,
                "___",
                .data$term_label
            )
        ) %>%
        dplyr::arrange(
            .data$immune_category,
            .data$GeneRatio_numeric,
            .data$p.adjust
        ) %>%
        dplyr::mutate(
            term_plot_id = factor(
                .data$term_plot_id,
                levels = unique(.data$term_plot_id)
            )
        )
    p <- ggplot2::ggplot(
        plot_df,
        ggplot2::aes(x = .data$GeneRatio_numeric, y = .data$term_plot_id)
    ) +
        ggplot2::geom_point(
            ggplot2::aes(size = .data$Count, color = .data$neg_log10_padj),
            alpha = 0.9
        ) +
        ggplot2::facet_grid(
            immune_category ~ .,
            scales = "free_y",
            space = "free_y"
        ) +
        ggplot2::scale_y_discrete(labels = function(x) {
            sub("^.*___", "", x)
        }) +
        ggplot2::scale_color_gradient(low = "grey75", high = "red3") +
        ggplot2::theme_bw(base_size = 10) +
        ggplot2::labs(
            x = "Gene ratio",
            y = NULL,
            size = "Gene count",
            color = "-log10 adjusted P",
            title = paste0(
                prefix,
                " | ",
                module_label_now,
                " | B/T-cell-relevant GO BP"
            )
        ) +
        ggplot2::theme(
            axis.text.y = ggplot2::element_text(size = 7),
            strip.text.y = ggplot2::element_text(angle = 0, size = 8),
            plot.title = ggplot2::element_text(hjust = 0.5),
            panel.grid.minor = ggplot2::element_blank()
        )
    plot_height <- max(6, 3.5 + 0.24 * nrow(plot_df))
    module_stub <- safe_name(module_label_now)
    ggplot2::ggsave(
        file.path(
            outdir,
            paste0(prefix, "__", module_stub, "__BT_relevant_GO_BP_dotplot.pdf")
        ),
        p,
        width = 10,
        height = plot_height,
        limitsize = FALSE
    )
    ggplot2::ggsave(
        file.path(
            outdir,
            paste0(prefix, "__", module_stub, "__BT_relevant_GO_BP_dotplot.png")
        ),
        p,
        width = 10,
        height = plot_height,
        dpi = 300,
        limitsize = FALSE
    )
    invisible(plot_df)
}
plot_bt_go_overview <- function(bt_all, outdir, prefix) {
    if (nrow(bt_all) == 0) {
        return(invisible(NULL))
    }
    overview_df <- bt_all %>%
        dplyr::group_by(.data$module) %>%
        dplyr::arrange(
            dplyr::desc(.data$specificity_score),
            .data$p.adjust,
            .by_group = TRUE
        ) %>%
        dplyr::slice_head(n = BT_GO_OVERVIEW_TOP_N_PER_MODULE) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(
            term_label = stringr::str_wrap(
                clean_go_term_name(.data$Description),
                width = 45
            ),
            term_plot_id = paste0(.data$module, "___", .data$term_label)
        ) %>%
        dplyr::arrange(
            .data$module_number,
            .data$GeneRatio_numeric,
            .data$p.adjust
        ) %>%
        dplyr::mutate(
            term_plot_id = factor(
                .data$term_plot_id,
                levels = unique(.data$term_plot_id)
            )
        )
    p <- ggplot2::ggplot(
        overview_df,
        ggplot2::aes(x = .data$GeneRatio_numeric, y = .data$term_plot_id)
    ) +
        ggplot2::geom_point(
            ggplot2::aes(size = .data$Count, color = .data$immune_category),
            alpha = 0.88
        ) +
        ggplot2::facet_wrap(~module, scales = "free_y") +
        ggplot2::scale_y_discrete(labels = function(x) {
            sub("^.*___", "", x)
        }) +
        ggplot2::theme_bw(base_size = 9) +
        ggplot2::labs(
            x = "Gene ratio",
            y = NULL,
            size = "Gene count",
            color = "Immune process",
            title = paste0(prefix, " B/T-cell-relevant GO BP overview")
        ) +
        ggplot2::theme(
            axis.text.y = ggplot2::element_text(size = 6.5),
            strip.text = ggplot2::element_text(size = 8),
            plot.title = ggplot2::element_text(hjust = 0.5),
            legend.position = "right"
        )
    plot_height <- max(8, 3 + 1.6 * length(unique(overview_df$module)))
    ggplot2::ggsave(
        file.path(outdir, paste0(prefix, "__BT_relevant_GO_BP_overview.pdf")),
        p,
        width = 16,
        height = plot_height,
        limitsize = FALSE
    )
    ggplot2::ggsave(
        file.path(outdir, paste0(prefix, "__BT_relevant_GO_BP_overview.png")),
        p,
        width = 16,
        height = plot_height,
        dpi = 300,
        limitsize = FALSE
    )
    invisible(overview_df)
}
run_go_enrichment <- function(modules_std, dat_expr, outdir, prefix) {
    modules_non_grey <- modules_std %>%
        dplyr::filter(tolower(.data$module) != "grey")
    universe_symbols <- colnames(dat_expr)
    universe_mapping <- map_symbols_to_entrez(universe_symbols)
    universe_entrez <- unique(universe_mapping$ENTREZID)
    readr::write_tsv(
        universe_mapping,
        file.path(outdir, paste0(prefix, "__GO_universe_SYMBOL_to_ENTREZ.tsv"))
    )
    module_mapping <- modules_non_grey %>%
        dplyr::transmute(
            module_original = .data$module,
            module = .data$module_label,
            module_id = .data$module_id,
            module_number = .data$module_number,
            SYMBOL = .data$gene_name
        ) %>%
        dplyr::left_join(universe_mapping, by = "SYMBOL") %>%
        dplyr::filter(!is.na(.data$ENTREZID), .data$ENTREZID != "") %>%
        dplyr::distinct(
            .data$module_original,
            .data$module,
            .data$SYMBOL,
            .data$ENTREZID,
            .keep_all = TRUE
        )
    readr::write_tsv(
        module_mapping,
        file.path(outdir, paste0(prefix, "__module_gene_SYMBOL_to_ENTREZ.tsv"))
    )
    go_universe_raw <- suppressMessages(AnnotationDbi::select(
        org.Hs.eg.db,
        keys = universe_entrez,
        keytype = "ENTREZID",
        columns = c("GOALL", "ONTOLOGYALL")
    )) %>%
        tibble::as_tibble() %>%
        dplyr::transmute(
            ENTREZID = as.character(.data$ENTREZID),
            ID = as.character(.data$GOALL),
            ONTOLOGY = as.character(.data$ONTOLOGYALL)
        ) %>%
        dplyr::filter(
            !is.na(.data$ENTREZID),
            .data$ENTREZID != "",
            !is.na(.data$ID),
            .data$ID != "",
            .data$ONTOLOGY %in% GO_ONTOLOGIES
        ) %>%
        dplyr::distinct(.data$ONTOLOGY, .data$ID, .data$ENTREZID)
    go_term_names <- suppressMessages(AnnotationDbi::select(
        GO.db::GO.db,
        keys = unique(go_universe_raw$ID),
        keytype = "GOID",
        columns = c("TERM", "ONTOLOGY")
    )) %>%
        tibble::as_tibble() %>%
        dplyr::transmute(
            ID = as.character(.data$GOID),
            Description = clean_go_term_name(as.character(.data$TERM)),
            GO_DB_ONTOLOGY = as.character(.data$ONTOLOGY)
        ) %>%
        dplyr::filter(!is.na(.data$ID), .data$ID != "") %>%
        dplyr::distinct(.data$ID, .keep_all = TRUE)
    bt_catalog <- build_bt_go_catalog(go_term_names)
    readr::write_tsv(
        bt_catalog,
        file.path(
            outdir,
            paste0(prefix, "__BT_relevant_GO_BP_catalog_used.tsv")
        )
    )
    all_results <- list()
    all_bt_results <- list()
    module_order <- get_module_mapping(modules_std) %>%
        dplyr::filter(tolower(.data$module) != "grey") %>%
        dplyr::arrange(.data$module_number)
    for (i in seq_len(nrow(module_order))) {
        module_now <- module_order$module[[i]]
        module_label_now <- module_order$module_label[[i]]
        module_id_now <- module_order$module_id[[i]]
        module_number_now <- module_order$module_number[[i]]
        module_map_now <- module_mapping %>%
            dplyr::filter(.data$module_original == module_now)
        module_entrez <- unique(module_map_now$ENTREZID)
        module_entrez <- intersect(module_entrez, universe_entrez)
        if (length(module_entrez) < 3) {
            message2(
                "Skip GO for ",
                prefix,
                " ",
                module_label_now,
                ": fewer than 3 mapped genes."
            )
            next
        }
        module_results <- list()
        for (ontology_now in GO_ONTOLOGIES) {
            term2gene <- go_universe_raw %>%
                dplyr::filter(.data$ONTOLOGY == ontology_now)
            term_sizes <- term2gene %>%
                dplyr::count(.data$ID, name = "BgCount") %>%
                dplyr::filter(
                    .data$BgCount >= GO_MIN_GS_SIZE,
                    .data$BgCount <= GO_MAX_GS_SIZE
                )
            term2gene <- term2gene %>% dplyr::semi_join(term_sizes, by = "ID")
            annotated_universe <- unique(term2gene$ENTREZID)
            module_annotated <- intersect(module_entrez, annotated_universe)
            N <- length(annotated_universe)
            n <- length(module_annotated)
            if (N == 0 || n == 0 || n > N) {
                message2(
                    "No annotated module genes for ",
                    prefix,
                    " ",
                    module_label_now,
                    " ",
                    ontology_now,
                    "."
                )
                next
            }
            overlap_long <- term2gene %>%
                dplyr::filter(.data$ENTREZID %in% module_annotated) %>%
                dplyr::inner_join(
                    module_map_now %>%
                        dplyr::select(.data$SYMBOL, .data$ENTREZID),
                    by = "ENTREZID"
                ) %>%
                dplyr::distinct(.data$ID, .data$ENTREZID, .data$SYMBOL)
            if (nrow(overlap_long) == 0) {
                next
            }
            overlap_summary <- overlap_long %>%
                dplyr::group_by(.data$ID) %>%
                dplyr::summarise(
                    Count = dplyr::n_distinct(.data$ENTREZID),
                    geneID = paste(sort(unique(.data$SYMBOL)), collapse = "/"),
                    ENTREZID_overlap = paste(
                        sort(unique(.data$ENTREZID)),
                        collapse = "/"
                    ),
                    .groups = "drop"
                )
            result_now <- overlap_summary %>%
                dplyr::left_join(term_sizes, by = "ID") %>%
                dplyr::mutate(
                    GeneRatio = paste0(.data$Count, "/", n),
                    BgRatio = paste0(.data$BgCount, "/", N),
                    RichFactor = .data$Count / .data$BgCount,
                    FoldEnrichment = (.data$Count / n) / (.data$BgCount / N),
                    pvalue = stats::phyper(
                        q = .data$Count - 1,
                        m = .data$BgCount,
                        n = N - .data$BgCount,
                        k = n,
                        lower.tail = FALSE
                    ),
                    p.adjust = stats::p.adjust(.data$pvalue, method = "BH")
                ) %>%
                dplyr::left_join(go_term_names, by = "ID") %>%
                dplyr::mutate(
                    Description = dplyr::if_else(
                        is.na(.data$Description) | .data$Description == "",
                        .data$ID,
                        clean_go_term_name(.data$Description)
                    ),
                    module_original = module_now,
                    module = module_label_now,
                    module_id = module_id_now,
                    module_number = module_number_now,
                    ONTOLOGY = ontology_now,
                    n_module_entrez = length(module_entrez),
                    n_module_annotated_ontology = n,
                    n_universe_entrez = length(universe_entrez),
                    n_universe_annotated_ontology = N,
                    qvalue = NA_real_,
                    enrichment_method = "one_sided_hypergeometric_GOALL_BH"
                ) %>%
                dplyr::arrange(
                    .data$p.adjust,
                    .data$pvalue,
                    dplyr::desc(.data$Count)
                ) %>%
                dplyr::select(
                    .data$module_original,
                    .data$module,
                    .data$module_id,
                    .data$module_number,
                    .data$ONTOLOGY,
                    .data$n_module_entrez,
                    .data$n_module_annotated_ontology,
                    .data$n_universe_entrez,
                    .data$n_universe_annotated_ontology,
                    .data$ID,
                    .data$Description,
                    .data$GeneRatio,
                    .data$BgRatio,
                    .data$RichFactor,
                    .data$FoldEnrichment,
                    .data$pvalue,
                    .data$p.adjust,
                    .data$qvalue,
                    .data$geneID,
                    .data$ENTREZID_overlap,
                    .data$Count,
                    .data$BgCount,
                    .data$enrichment_method
                )
            module_results[[ontology_now]] <- result_now
        }
        module_all <- dplyr::bind_rows(module_results)
        if (nrow(module_all) == 0) {
            message2(
                "No GO results returned for ",
                prefix,
                " ",
                module_label_now
            )
            next
        }
        module_sig <- module_all %>%
            dplyr::filter(
                !is.na(.data$p.adjust),
                .data$p.adjust <= GO_PADJ_CUTOFF
            ) %>%
            dplyr::arrange(
                .data$ONTOLOGY,
                .data$p.adjust,
                dplyr::desc(.data$Count)
            )
        module_bt <- select_bt_relevant_go(module_all, bt_catalog)
        module_stub <- safe_name(module_label_now)
        readr::write_tsv(
            module_all,
            file.path(
                outdir,
                paste0(
                    prefix,
                    "__",
                    module_stub,
                    "__GO_all_overlapping_tested_terms.tsv"
                )
            )
        )
        readr::write_tsv(
            module_sig,
            file.path(
                outdir,
                paste0(
                    prefix,
                    "__",
                    module_stub,
                    "__GO_significant_padj_",
                    GO_PADJ_CUTOFF,
                    ".tsv"
                )
            )
        )
        readr::write_tsv(
            module_bt,
            file.path(
                outdir,
                paste0(prefix, "__", module_stub, "__BT_relevant_GO_BP.tsv")
            )
        )
        if (nrow(module_sig) > 0) {
            plot_df <- module_sig %>%
                dplyr::mutate(
                    GeneRatio_numeric = parse_gene_ratio(.data$GeneRatio),
                    neg_log10_padj = -log10(
                        .data$p.adjust +
                            1e-300
                    ),
                    term_label = stringr::str_wrap(
                        clean_go_term_name(.data$Description),
                        width = 60
                    ),
                    term_plot_id = paste0(
                        .data$ONTOLOGY,
                        "___",
                        .data$term_label
                    )
                ) %>%
                dplyr::filter(!is.na(.data$GeneRatio_numeric))
            if (is.finite(GO_MAX_SIGNIFICANT_TERMS_PER_MODULE)) {
                plot_df <- plot_df %>%
                    dplyr::group_by(.data$ONTOLOGY) %>%
                    dplyr::arrange(
                        .data$p.adjust,
                        dplyr::desc(.data$GeneRatio_numeric),
                        .by_group = TRUE
                    ) %>%
                    dplyr::slice_head(
                        n = GO_MAX_SIGNIFICANT_TERMS_PER_MODULE
                    ) %>%
                    dplyr::ungroup()
            }
            plot_df <- plot_df %>%
                dplyr::group_by(.data$ONTOLOGY) %>%
                dplyr::arrange(
                    .data$GeneRatio_numeric,
                    .data$p.adjust,
                    .by_group = TRUE
                ) %>%
                dplyr::mutate(
                    term_plot_id = factor(
                        .data$term_plot_id,
                        levels = unique(.data$term_plot_id)
                    )
                ) %>%
                dplyr::ungroup()
            go_plot <- ggplot2::ggplot(
                plot_df,
                ggplot2::aes(
                    x = .data$GeneRatio_numeric,
                    y = .data$term_plot_id
                )
            ) +
                ggplot2::geom_point(
                    ggplot2::aes(
                        size = .data$Count,
                        color = .data$neg_log10_padj
                    ),
                    alpha = 0.9
                ) +
                ggplot2::facet_grid(
                    ONTOLOGY ~ .,
                    scales = "free_y",
                    space = "free_y"
                ) +
                ggplot2::scale_y_discrete(labels = function(x) {
                    sub("^.*___", "", x)
                }) +
                ggplot2::scale_color_gradient(low = "grey75", high = "red3") +
                ggplot2::theme_bw(base_size = 10) +
                ggplot2::labs(
                    x = "Gene ratio",
                    y = NULL,
                    size = "Gene count",
                    color = "-log10 adjusted P",
                    title = paste0(
                        prefix,
                        " | ",
                        module_label_now,
                        " | all significant GO terms"
                    )
                ) +
                ggplot2::theme(
                    axis.text.y = ggplot2::element_text(size = 7),
                    axis.text.x = ggplot2::element_text(size = 8),
                    strip.text.y = ggplot2::element_text(angle = 0),
                    plot.title = ggplot2::element_text(hjust = 0.5),
                    panel.grid.minor = ggplot2::element_blank()
                )
            plot_height <- max(6, 3.5 + 0.22 * nrow(plot_df))
            ggplot2::ggsave(
                file.path(
                    outdir,
                    paste0(
                        prefix,
                        "__",
                        module_stub,
                        "__GO_all_significant_terms.pdf"
                    )
                ),
                go_plot,
                width = 10,
                height = plot_height,
                limitsize = FALSE
            )
            ggplot2::ggsave(
                file.path(
                    outdir,
                    paste0(
                        prefix,
                        "__",
                        module_stub,
                        "__GO_all_significant_terms.png"
                    )
                ),
                go_plot,
                width = 10,
                height = plot_height,
                dpi = 300,
                limitsize = FALSE
            )
        }
        if (nrow(module_bt) > 0) {
            plot_bt_go_one_module(
                bt_df = module_bt,
                outdir = outdir,
                prefix = prefix,
                module_label_now = module_label_now
            )
        } else {
            message2(
                "No B/T-cell-focused GO BP passed filters for ",
                prefix,
                " ",
                module_label_now,
                "."
            )
        }
        all_results[[module_label_now]] <- module_all
        all_bt_results[[module_label_now]] <- module_bt
    }
    all_tested <- dplyr::bind_rows(all_results)
    if (nrow(all_tested) == 0) {
        message2("No GO enrichment result was generated for ", prefix)
        return(list(
            all_tested = tibble::tibble(),
            significant = tibble::tibble(),
            bt_relevant = tibble::tibble()
        ))
    }
    all_significant <- all_tested %>%
        dplyr::filter(
            !is.na(.data$p.adjust),
            .data$p.adjust <= GO_PADJ_CUTOFF
        ) %>%
        dplyr::arrange(.data$module_number, .data$ONTOLOGY, .data$p.adjust)
    bt_all <- dplyr::bind_rows(all_bt_results)
    if (nrow(bt_all) > 0) {
        bt_all <- bt_all %>%
            dplyr::arrange(
                .data$module_number,
                .data$priority,
                dplyr::desc(.data$specificity_score)
            )
    }
    readr::write_tsv(
        all_tested,
        file.path(
            outdir,
            paste0(prefix, "__all_modules_GO_all_overlapping_tested_terms.tsv")
        )
    )
    readr::write_tsv(
        all_significant,
        file.path(
            outdir,
            paste0(
                prefix,
                "__all_modules_GO_significant_padj_",
                GO_PADJ_CUTOFF,
                ".tsv"
            )
        )
    )
    readr::write_tsv(
        bt_all,
        file.path(outdir, paste0(prefix, "__all_modules_BT_relevant_GO_BP.tsv"))
    )
    if (nrow(all_significant) > 0) {
        overview_df <- all_significant %>%
            dplyr::mutate(
                GeneRatio_numeric = parse_gene_ratio(.data$GeneRatio),
                neg_log10_padj = -log10(
                    .data$p.adjust +
                        1e-300
                ),
                term_label = stringr::str_wrap(
                    stringr::str_trunc(
                        clean_go_term_name(.data$Description),
                        65
                    ),
                    width = 50
                ),
                term_plot_id = paste0(
                    .data$module,
                    "___",
                    .data$ONTOLOGY,
                    "___",
                    .data$term_label
                )
            ) %>%
            dplyr::group_by(.data$module, .data$ONTOLOGY) %>%
            dplyr::arrange(
                .data$p.adjust,
                dplyr::desc(.data$GeneRatio_numeric),
                .by_group = TRUE
            ) %>%
            dplyr::slice_head(n = GO_OVERVIEW_TOP_N_PER_MODULE_ONTOLOGY) %>%
            dplyr::ungroup() %>%
            dplyr::arrange(
                .data$module_number,
                .data$ONTOLOGY,
                .data$GeneRatio_numeric
            ) %>%
            dplyr::mutate(
                term_plot_id = factor(
                    .data$term_plot_id,
                    levels = unique(.data$term_plot_id)
                )
            )
        overview_plot <- ggplot2::ggplot(
            overview_df,
            ggplot2::aes(x = .data$neg_log10_padj, y = .data$term_plot_id)
        ) +
            ggplot2::geom_point(
                ggplot2::aes(size = .data$Count, color = .data$ONTOLOGY),
                alpha = 0.85
            ) +
            ggplot2::facet_wrap(~ module + ONTOLOGY, scales = "free_y") +
            ggplot2::scale_y_discrete(labels = function(x) {
                sub("^.*___", "", x)
            }) +
            ggplot2::theme_bw(base_size = 9) +
            ggplot2::labs(
                x = "-log10 adjusted P value",
                y = NULL,
                size = "Gene count",
                color = "Ontology",
                title = paste0(prefix, " complete GO enrichment overview")
            ) +
            ggplot2::theme(
                axis.text.y = ggplot2::element_text(size = 6),
                strip.text = ggplot2::element_text(size = 7),
                plot.title = ggplot2::element_text(hjust = 0.5)
            )
        overview_height <- max(8, 3 + 1.8 * length(unique(overview_df$module)))
        ggplot2::ggsave(
            file.path(outdir, paste0(prefix, "__GO_enrichment_overview.pdf")),
            overview_plot,
            width = 16,
            height = overview_height,
            limitsize = FALSE
        )
        ggplot2::ggsave(
            file.path(outdir, paste0(prefix, "__GO_enrichment_overview.png")),
            overview_plot,
            width = 16,
            height = overview_height,
            dpi = 300,
            limitsize = FALSE
        )
    }
    if (nrow(bt_all) > 0) {
        plot_bt_go_overview(bt_all = bt_all, outdir = outdir, prefix = prefix)
    }
    list(
        all_tested = all_tested,
        significant = all_significant,
        bt_relevant = bt_all
    )
}
run_one_celltype <- function(final_rds) {
    current_stage <- "readRDS"
    tryCatch(
        {
            message2(
                "\n############################################################"
            )
            message2("Input: ", final_rds)
            message2(
                "############################################################"
            )
            seu <- readRDS(final_rds)
            current_stage <- "identify_celltype_and_wgcna"
            celltype <- get_celltype(seu, final_rds)
            prefix <- safe_name(celltype)
            wgcna_name <- get_wgcna_name_safe(seu, celltype)
            celltype_dir <- dirname(final_rds)
            dirs <- make_output_dirs(celltype_dir)
            completion_file <- file.path(
                dirs$root,
                paste0(prefix, "__DOWNSTREAM_COMPLETED.txt")
            )
            if (file.exists(completion_file) && !OVERWRITE_OUTPUT) {
                message2("Skip completed downstream result: ", celltype)
                return(tibble::tibble(
                    celltype = celltype,
                    final_rds = final_rds,
                    status = "skipped_existing",
                    output_dir = dirs$root
                ))
            }
            message2("Cell type: ", celltype)
            message2("hdWGCNA experiment: ", wgcna_name)
            current_stage <- "retrieve_modules_and_hubs"
            modules_raw <- hdWGCNA::GetModules(seu, wgcna_name = wgcna_name)
            modules_std <- standardize_modules(modules_raw)
            modules_non_grey <- modules_std %>%
                dplyr::filter(tolower(.data$module) != "grey")
            module_mapping <- get_module_mapping(modules_std)
            hub_raw <- hdWGCNA::GetHubGenes(
                seu,
                n_hubs = HUB_N_GENES_PER_MODULE,
                mods = module_mapping$module[
                    tolower(module_mapping$module) != "grey"
                ],
                wgcna_name = wgcna_name
            )
            hubs_std <- standardize_hubs(hub_raw, modules_std)
            readr::write_tsv(
                module_mapping,
                file.path(
                    dirs$tables,
                    paste0(prefix, "__module_color_to_M_number_mapping.tsv")
                )
            )
            readr::write_tsv(
                modules_std,
                file.path(
                    dirs$tables,
                    paste0(prefix, "__modules_standardized.tsv")
                )
            )
            readr::write_tsv(
                hubs_std,
                file.path(
                    dirs$tables,
                    paste0(prefix, "__top_hub_genes_standardized.tsv")
                )
            )
            current_stage <- "retrieve_datExpr_and_power"
            dat_expr <- get_dat_expr_safe(seu, wgcna_name)
            selected_power <- get_selected_soft_power(
                celltype_dir,
                seu,
                wgcna_name
            )
            message2(
                "datExpr: ",
                nrow(dat_expr),
                " profiles x ",
                ncol(dat_expr),
                " genes; selected soft power: ",
                selected_power
            )
            current_stage <- "hub_gene_network"
            network_result <- plot_hub_networks(
                seu = seu,
                modules_std = modules_std,
                hubs_std = hubs_std,
                dat_expr = dat_expr,
                selected_power = selected_power,
                outdir = dirs$network,
                prefix = prefix
            )
            current_stage <- "module_eigengene_heatmaps"
            me_df <- get_me_matrix_safe(
                seu,
                wgcna_name = wgcna_name,
                harmonized = USE_HARMONIZED_ME
            )
            sample_me_result <- prepare_sample_me_summary(
                seu = seu,
                me_df = me_df,
                modules_std = modules_std
            )
            readr::write_tsv(
                sample_me_result$me_meta_labeled,
                file.path(
                    dirs$tables,
                    paste0(
                        prefix,
                        "__cell_level_module_eigengenes_with_metadata.tsv"
                    )
                )
            )
            readr::write_tsv(
                sample_me_result$module_col_map,
                file.path(
                    dirs$tables,
                    paste0(prefix, "__ME_column_to_module_label_mapping.tsv")
                )
            )
            module_heatmap_result <- plot_module_activity_heatmaps(
                sample_me = sample_me_result$sample_summary,
                module_cols = sample_me_result$module_cols,
                module_col_map = sample_me_result$module_col_map,
                outdir = dirs$heatmap,
                prefix = prefix
            )
            current_stage <- "hub_gene_expression_heatmap"
            hub_heatmap_result <- plot_hub_gene_expression_heatmap(
                seu = seu,
                hubs_std = hubs_std,
                sample_me = sample_me_result$sample_summary,
                outdir = dirs$heatmap,
                prefix = prefix
            )
            current_stage <- "GO_enrichment"
            go_result <- run_go_enrichment(
                modules_std = modules_std,
                dat_expr = dat_expr,
                outdir = dirs$go,
                prefix = prefix
            )
            completion_text <- c(
                paste0("celltype\t", celltype),
                paste0("wgcna_name\t", wgcna_name),
                paste0("selected_soft_power\t", selected_power),
                paste0(
                    "non_grey_modules\t",
                    length(unique(modules_non_grey$module))
                ),
                paste0("hub_network_nodes\t", nrow(network_result$nodes)),
                paste0("hub_network_edges\t", nrow(network_result$edges)),
                paste0("GO_all_tested_rows\t", nrow(go_result$all_tested)),
                paste0("GO_significant_rows\t", nrow(go_result$significant)),
                paste0("BT_relevant_GO_BP_rows\t", nrow(go_result$bt_relevant)),
                paste0(
                    "completed_at\t",
                    format(Sys.time(), "%Y-%m-%d %H:%M:%S")
                )
            )
            readr::write_lines(completion_text, completion_file)
            message2("Completed downstream analysis: ", celltype)
            tibble::tibble(
                celltype = celltype,
                final_rds = final_rds,
                status = "completed",
                selected_soft_power = selected_power,
                non_grey_modules = length(unique(modules_non_grey$module)),
                hub_network_nodes = nrow(network_result$nodes),
                hub_network_edges = nrow(network_result$edges),
                GO_all_tested_rows = nrow(go_result$all_tested),
                GO_significant_rows = nrow(go_result$significant),
                BT_relevant_GO_BP_rows = nrow(go_result$bt_relevant),
                output_dir = dirs$root
            )
        },
        error = function(e) {
            error_message <- paste0(
                "[stage: ",
                current_stage,
                "] ",
                conditionMessage(e)
            )
            message2("ERROR: ", error_message)
            error_file <- file.path(
                dirname(final_rds),
                paste0(safe_name(basename(final_rds)), "__DOWNSTREAM_ERROR.tsv")
            )
            error_df <- tibble::tibble(
                final_rds = final_rds,
                stage = current_stage,
                error = error_message
            )
            readr::write_tsv(error_df, error_file)
            tibble::tibble(
                celltype = safe_name(basename(dirname(final_rds))),
                final_rds = final_rds,
                status = "error",
                selected_soft_power = NA_real_,
                non_grey_modules = NA_integer_,
                hub_network_nodes = NA_integer_,
                hub_network_edges = NA_integer_,
                GO_all_tested_rows = NA_integer_,
                GO_significant_rows = NA_integer_,
                BT_relevant_GO_BP_rows = NA_integer_,
                output_dir = dirname(final_rds),
                error = error_message
            )
        }
    )
}
final_rds_files <- list.files(
    RESULT_ROOT,
    pattern = FINAL_RDS_PATTERN,
    recursive = TRUE,
    full.names = TRUE
)
if (!is.null(TARGET_CELLTYPES)) {
    target_safe <- safe_name(TARGET_CELLTYPES)
    final_rds_files <- final_rds_files[
        safe_name(basename(dirname(final_rds_files))) %in% target_safe
    ]
}
final_rds_files <- sort(unique(final_rds_files))
message2("Found completed hdWGCNA final RDS files: ", length(final_rds_files))
summary_list <- lapply(final_rds_files, run_one_celltype)
summary_df <- dplyr::bind_rows(summary_list)
readr::write_tsv(
    summary_df,
    file.path(
        RESULT_ROOT,
        "hdWGCNA_downstream_hub_heatmap_GO_BT_selected_summary.tsv"
    )
)
message2("\nAll available completed cell types have been processed.")
message2(
    "Completed: ",
    sum(summary_df$status == "completed", na.rm = TRUE),
    "; skipped: ",
    sum(summary_df$status == "skipped_existing", na.rm = TRUE),
    "; errors: ",
    sum(summary_df$status == "error", na.rm = TRUE)
)
