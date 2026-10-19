#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(ArchR)
    library(CoVarNet)
    library(NMF)
    library(dplyr)
    library(tidyr)
    library(tibble)
    library(readr)
    library(tidytext)
    library(ggplot2)
    library(patchwork)
    library(ggpubr)
    library(extrafont)
    library(viridis)
    library(ComplexHeatmap)
    library(circlize)
    library(ggsci)
    library(grid)
})


# 准备 CoVarNet 输入数据
proj <- readRDS(
    "/public/users/wangrong/home/wangrong/results/T_B/TB_final_downstream/proj_TB_only_annotated_with_new_UMAP.rds"
)
old_labels <- as.character(getCellColData(proj, select = "TB_Subtype_All")[, 1])
rename_map <- c(
    Activated_B = "Intermediate_B",
    Resting_B = "Memory_B",
    CD4_Teff = "CD8_Memory_GZMK",
    Th1 = "CD4_Memory_ZEB2",
    Th2 = "CD4_Memory",
    CD4_Tcm = "CD4_Naive_CCR4",
    CD8_Tem = "CD8_Memory"
)
new_labels <- old_labels
idx <- old_labels %in% names(rename_map)
new_labels[idx] <- unname(rename_map[old_labels[idx]])
proj <- addCellColData(
    ArchRProj = proj,
    data = new_labels,
    cells = getCellNames(proj),
    name = "TB_Subtype_All",
    force = TRUE
)

cell_meta <- getCellColData(proj)
cell_meta_df <- as.data.frame(cell_meta) %>% rownames_to_column("cell_id")
meta <- cell_meta_df %>%
    dplyr::select(
        cell_id,
        sampleID = Sample,
        Donor,
        subCluster = TB_Subtype_All,
        Timepoint
    ) %>%
    mutate(
        majorCluster = case_when(
            subCluster %in%
                c("Memory_B", "Naive_B", "PlasmaB", "Intermediate_B") ~ "B",
            subCluster %in%
                c(
                    "CD4_NaiveT",
                    "CD4_Naive_CCR4",
                    "CD4_Memory_ZEB2",
                    "CD4_Memory",
                    "Treg"
                ) ~ "CD4_T",
            subCluster %in%
                c("CD8_NaiveT", "CD8_Memory", "CD8_Memory_GZMK") ~ "CD8_T",
            subCluster %in%
                c("gdT", "CD4_CD8_DP") ~ "other_T"
        ),
        Dose = case_when(
            Timepoint %in% c("day0", "day1", "day3", "day6", "day14") ~ "Dose1",
            Timepoint %in%
                c("day30", "day31", "day33", "day36", "day44") ~ "Dose2"
        )
    )
meta <- as.data.frame(meta)
rownames(meta) <- meta$cell_id
meta$cell_id <- NULL
write.table(
    meta,
    file = "/public/users/huhuajie/projects/T_B_CoVarNet/results/meta_CoVarNet.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = TRUE,
    col.names = NA
)

# CoVarNet 分析及图表输出

SUPP_DIR <- "/public/users/huhuajie/projects/T_B_CoVarNet/results/supplementary_tables"
dir.create(SUPP_DIR, recursive = TRUE, showWarnings = FALSE)
meta <- read.table(
    "/public/users/huhuajie/projects/T_B_CoVarNet/results/meta_CoVarNet.tsv",
    header = TRUE,
    sep = "\t",
    row.names = 1,
    quote = "",
    comment.char = "",
    check.names = FALSE
)
rt_sp <- names(table(meta$sampleID))[table(meta$sampleID) > 100]
meta <- meta[meta$sampleID %in% rt_sp, ]
mat_fq_raw <- freq_calculate(meta)
mat_fq_norm <- freq_normalize(mat_fq_raw, normalize = "zscore")
res <- nmf(
    mat_fq_norm,
    rank = 2:20,
    method = "nsNMF",
    seed = rep(123456, 6),
    .options = "vp"
)
# 图 S7A：NMF 聚类数评估
plot(res)
K <- 5
NMF <- nmf(mat_fq_norm, K, method = "nsNMF", seed = rep(123456, 6), nrun = 100)
colnames(basis(NMF)) <- sprintf("CM%02d", seq_len(ncol(basis(NMF))))
rownames(coef(NMF)) <- sprintf("CM%02d", seq_len(nrow(coef(NMF))))
w <- basis(NMF)
# 图 S7B、图 5A 和图 5C：模块权重、时序分布及主要细胞组成
gr.weight_all(NMF)
gr.distribution(NMF, meta = meta, group = "Timepoint")
gr.distribution(NMF, meta = meta, group = "Dose")
h <- coef(NMF)
max_cm <- apply(h, 2, function(x) rownames(h)[which.max(x)])
max_cm <- gsub("CM", "CMT", max_cm)
max_cm
p <- gr.weight_top(NMF, num = 10)
p +
    facet_wrap(~cm, scales = "free", ncol = 5, dir = "v") +
    coord_flip(clip = "off") +
    theme(
        aspect.ratio = NULL,
        axis.text.y = element_text(size = 7),
        strip.text = element_text(size = 8),
        panel.spacing = unit(0.8, "lines")
    )
cor_pair <- pair_correlation(mat_fq_raw, method = "pearson")
top_n <- 10
network <- cm_network(NMF, cor_pair, top_n = top_n, corr = 0.2, pval_fdr = 0.05)
# 输出 CoVarNet 模块组成、网络边和样本模块活性补充表
celltype_annotation <- meta %>%
    as.data.frame() %>%
    distinct(subCluster, majorCluster)
final_network_nodes <- network$filter %>%
    as_tibble() %>%
    transmute(
        CM = as.character(cm),
        cell_subset = as.character(subCluster),
        retained_in_final_network = TRUE
    ) %>%
    distinct()
supp_table_1 <- as.data.frame(w) %>%
    rownames_to_column("cell_subset") %>%
    pivot_longer(
        cols = -cell_subset,
        names_to = "CM",
        values_to = "NMF_weight"
    ) %>%
    left_join(celltype_annotation, by = c(cell_subset = "subCluster")) %>%
    group_by(CM) %>%
    arrange(desc(NMF_weight), .by_group = TRUE) %>%
    mutate(
        rank_within_CM = row_number(),
        NMF_weight_fraction = NMF_weight / sum(NMF_weight),
        top10_NMF_member = rank_within_CM <= top_n
    ) %>%
    ungroup() %>%
    left_join(final_network_nodes, by = c("CM", "cell_subset")) %>%
    mutate(
        retained_in_final_network = replace_na(retained_in_final_network, FALSE)
    ) %>%
    transmute(
        CM,
        cell_subset,
        major_compartment = majorCluster,
        NMF_weight,
        NMF_weight_fraction,
        rank_within_CM,
        top10_NMF_member,
        retained_in_final_network
    ) %>%
    arrange(CM, rank_within_CM)
write_tsv(
    supp_table_1,
    file.path(SUPP_DIR, "Table_CoVarNet_module_composition_and_NMF_weights.tsv")
)
n_subsets <- length(unique(c(
    as.character(cor_pair$subCluster1),
    as.character(cor_pair$subCluster2)
)))
specificity_cutoff <- 1 - ((top_n - 1) * 2 - 1) / ((n_subsets - 1) * 2 - 1)
supp_table_2 <- network$each$edge %>%
    as_tibble() %>%
    transmute(
        CM = as.character(cm),
        cell_subset_1 = as.character(subCluster1),
        major_compartment_1 = as.character(majorCluster1),
        cell_subset_2 = as.character(subCluster2),
        major_compartment_2 = as.character(majorCluster2),
        Pearson_r = as.numeric(correlation),
        P_value = as.numeric(pval),
        FDR = as.numeric(pval_fdr),
        correlation_specificity = as.numeric(spe),
        correlation_cutoff = 0.2,
        FDR_cutoff = 0.05,
        specificity_cutoff = specificity_cutoff
    ) %>%
    arrange(CM, desc(Pearson_r), FDR)
write_tsv(
    supp_table_2,
    file.path(
        SUPP_DIR,
        "Table_CoVarNet_specifically_correlated_network_edges.tsv"
    )
)
each <- network$each
df_cm <- as.data.frame(each$node)
gr_igraph_each_mfrow <- function(each, nrow = 2, ncol = 5, ...) {
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)
    base_par <- graphics::par
    e <- new.env(parent = environment(gr.igraph_each))
    e$par <- function(...) {
        args <- list(...)
        if (!is.null(args$mfrow)) {
            args$mfrow <- c(nrow, ncol)
        }
        do.call(base_par, args)
    }
    f <- gr.igraph_each
    environment(f) <- e
    f(each, ...)
}
gr_igraph_each_mfrow(each, nrow = 1, ncol = 5)
global <- network$global
par(plt = c(0, 1, 0, 1), fig = c(0, 1, 0, 1))
gr.igraph_global(global, Layout = layout_with_fr)
cs <- colSums(h)
h <- sweep(h, 2, cs, "/")
meta_sp <- meta[, c("sampleID", "Timepoint", "Dose", "Donor")]
meta_sp[[1]] <- as.character(meta_sp[[1]])
meta_sp <- meta_sp[!duplicated(meta_sp[[1]]), , drop = FALSE]
rownames(meta_sp) <- meta_sp[[1]]
meta_sp <- meta_sp[, -1, drop = FALSE]
freq <- mat_fq_raw
h <- h[, rownames(meta_sp)]
freq <- freq[, rownames(meta_sp)]
meta_sp <- cbind(meta_sp, t(h))
meta_sp <- cbind(meta_sp, t(freq))
sample_metadata <- meta %>%
    as.data.frame() %>%
    transmute(
        sampleID = as.character(sampleID),
        Donor = as.character(Donor),
        Timepoint = as.character(Timepoint),
        Dose = as.character(Dose)
    ) %>%
    distinct()
cm_columns <- rownames(h)
sample_cm_wide <- as.data.frame(t(h)) %>% rownames_to_column("sampleID")
sample_cm_matrix <- as.matrix(sample_cm_wide[, cm_columns, drop = FALSE])
sample_cm_wide$dominant_CM <- cm_columns[max.col(
    sample_cm_matrix,
    ties.method = "first"
)]
sample_cm_wide$dominant_CM_activity <- apply(
    sample_cm_matrix,
    1,
    max,
    na.rm = TRUE
)
supp_table_3 <- sample_cm_wide %>%
    left_join(sample_metadata, by = "sampleID") %>%
    mutate(
        timepoint_day = suppressWarnings(as.integer(gsub(
            "[^0-9]",
            "",
            Timepoint
        )))
    ) %>%
    pivot_longer(
        cols = all_of(cm_columns),
        names_to = "CM",
        values_to = "CM_activity"
    ) %>%
    dplyr::select(
        sampleID,
        Donor,
        Timepoint,
        timepoint_day,
        Dose,
        CM,
        CM_activity,
        dominant_CM,
        dominant_CM_activity
    ) %>%
    arrange(timepoint_day, Donor, CM)
write_tsv(
    supp_table_3,
    file.path(SUPP_DIR, "Table_CoVarNet_sample_level_module_activities.tsv")
)
writexl::write_xlsx(
    list(
        CM_module_weights = supp_table_1,
        CM_network_edges = supp_table_2,
        Sample_CM_activities = supp_table_3
    ),
    path = file.path(SUPP_DIR, "CoVarNet_supplementary_tables.xlsx")
)
cell_type_groups <- list(
    CM01 = c("Intermediate_B", "CD4_CD8_DP", "CD8_Memory", "PlasmaB"),
    CM03 = c("CD4_NaiveT", "CD8_NaiveT", "Naive_B", "gdT")
)
timepoint_order <- c(
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
for (group_name in names(cell_type_groups)) {
    cell_types <- cell_type_groups[[group_name]]
    meta_sp_long <- meta_sp %>%
        pivot_longer(
            cols = all_of(cell_types),
            names_to = "celltype",
            values_to = "activity"
        ) %>%
        mutate(celltype = factor(celltype, levels = cell_types))
    meta_sp_long$Timepoint <- factor(
        meta_sp_long$Timepoint,
        levels = timepoint_order
    )
    # 图 S7C：指定模块细胞类型的时间动态
    p2 <- ggplot(meta_sp_long, aes(Timepoint, activity, col = Timepoint)) +
        geom_boxplot(outlier.shape = NA, linewidth = 0.3, alpha = 0) +
        geom_jitter(size = 0.5, width = 0.3) +
        stat_summary(
            aes(group = 1),
            fun = median,
            geom = "line",
            color = "black",
            linewidth = 0.5
        ) +
        stat_summary(
            aes(group = 1),
            fun = median,
            geom = "point",
            color = "black",
            size = 1
        ) +
        scale_color_viridis_d(option = "cividis") +
        stat_compare_means(
            size = 2,
            method = "kruskal.test",
            label.x.npc = 0.05,
            label.y.npc = 0.95
        ) +
        facet_wrap(
            ~celltype,
            ncol = 4,
            scales = "free_y",
            axes = "all_x",
            axis.labels = "all_x"
        ) +
        labs(x = "", y = "Cell type activity") +
        theme_classic() +
        theme(
            text = element_text(size = 6, color = "black"),
            axis.text.x = element_text(angle = 45, hjust = 1),
            strip.background = element_blank(),
            legend.position = "none",
            aspect.ratio = 1.2 / 2
        )
    p2
}
h <- h[, rownames(meta_sp), drop = FALSE]
cm_names <- rownames(h)
timepoint_order <- c(
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
donor_order <- sort(unique(as.character(meta_sp$Donor)))
meta_h <- meta_sp %>%
    dplyr::mutate(
        Donor = factor(Donor, levels = donor_order),
        Timepoint = factor(Timepoint, levels = timepoint_order)
    )
timepoint_order <- c(
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
meta_h <- meta_sp %>%
    mutate(Timepoint = factor(Timepoint, levels = timepoint_order))
h_timepoint <- meta_h %>%
    group_by(Timepoint) %>%
    summarise(
        across(all_of(cm_names), ~ mean(.x, na.rm = TRUE)),
        .groups = "drop"
    )
h_tp <- h_timepoint %>% column_to_rownames("Timepoint") %>% as.matrix() %>% t()
h_tp <- h_tp[, timepoint_order, drop = FALSE]
roe <- sweep(h_tp, MARGIN = 1, STATS = rowSums(h_tp), FUN = "/") * sum(h_tp)
ann <- data.frame(
    Dose = ifelse(
        colnames(roe) %in% c("day0", "day1", "day3", "day6", "day14"),
        "Dose1",
        "Dose2"
    ),
    row.names = colnames(roe)
)
ann$Dose <- factor(ann$Dose, levels = c("Dose1", "Dose2"))
dose_col <- c(Dose1 = "#E64B35", Dose2 = "#3C5488")
ann_top <- HeatmapAnnotation(
    Dose = ann$Dose,
    col = list(Dose = dose_col),
    simple_anno_size = unit(3, "mm"),
    border = TRUE,
    annotation_name_gp = gpar(fontsize = 6),
    annotation_legend_param = list(
        labels_gp = gpar(fontsize = 6),
        title_gp = gpar(fontsize = 6),
        grid_width = unit(3, "mm"),
        grid_height = unit(2, "mm")
    )
)
col_fun <- colorRamp2(
    c(0, 1, max(3, max(roe, na.rm = TRUE))),
    c("#F7F7F7", "#F4A582", "#B2182B")
)
# 绘制各时间点的模块富集热图
p3 <- Heatmap(
    roe,
    name = "Propeller",
    col = col_fun,
    border = TRUE,
    cluster_rows = TRUE,
    clustering_method_rows = "ward.D2",
    cluster_columns = FALSE,
    width = unit(10, "cm"),
    height = unit(max(3.5, nrow(roe) * 0.35), "cm"),
    column_title = "CM abundance prevalence across timepoints",
    column_title_gp = gpar(fontsize = 8),
    column_names_rot = 90,
    column_names_gp = gpar(fontsize = 6),
    row_names_gp = gpar(fontsize = 6),
    top_annotation = ann_top,
    row_dend_width = unit(5, "mm"),
    cell_fun = function(j, i, x, y, width, height, fill) {
        if (roe[i, j] > 3) {
            grid.text("+++", x, y, gp = gpar(fontsize = 5))
        } else if (roe[i, j] > 2) {
            grid.text("++", x, y, gp = gpar(fontsize = 5))
        } else if (roe[i, j] > 1) {
            grid.text("+", x, y, gp = gpar(fontsize = 5))
        }
    },
    heatmap_legend_param = list(
        at = c(0, 1, 2, 3),
        labels = c("0", "1", "2", "3"),
        title = "Ro/e",
        grid_width = unit(3, "mm"),
        labels_gp = gpar(fontsize = 6),
        title_gp = gpar(fontsize = 6)
    )
)
draw(p3, heatmap_legend_side = "right", annotation_legend_side = "right")
