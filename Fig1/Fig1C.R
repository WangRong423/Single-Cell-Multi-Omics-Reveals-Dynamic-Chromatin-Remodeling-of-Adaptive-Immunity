library(ArchR)
library(ggplot2)
library(Matrix)
library(SummarizedExperiment)

# ============================================================
# 1. 使用当前 ArchR 对象
# ============================================================
proj <- TB_proj

OUT_DIR <- paste0(
  "/public/users/wangrong/home/wangrong/results/T_B/",
  "TB_final_downstream/Subtype_marker_dotplot"
)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 2. Marker 分组，与参考图一致
# ============================================================
all_marker <- list(
    "B cell" = c("MS4A1", "CD19", "CD79A", "CD79B"),
    "Naive B" = c("FCER2", "TCL1A", "IL4R"),
    "Memory B" = c("AIM2", "TNFRSF13B", "GPR183"),
    "Plasma B" = c("MZB1", "IGHG1", "SDC1", "IGHA1", "CD27", "CD38"),
    "T cell" = c("CD3D", "CD3E", "CD3G", "CD4", "CD8A", "CD8B"),
    "Naive" = c("TCF7", "CCR7", "SELL", "LEF1"),"CCR4"
=c("CCR4"),"Memory T" = c("LMNA", "ANXA1", "S100A4"),
    "GZMK" = c("ZEB2", "GZMK"),
	"Gamma Delta T" = c("TRDC", "TRGC1", "TRDV2", "TRGV9"),
    "CD4 Treg" = c("IKZF2", "FOXP3", "IL2RA", "CTLA4")
    
)

celltype_levels <- c(
  "Naive_B",
  "Intermediate_B",
  "Memory_B",
  "PlasmaB",
  "CD4_NaiveT",
  "CD4_Naive_CCR4",
  "CD4_Memory",
   "CD4_Memory_ZEB2",
  "CD8_NaiveT",
  "CD8_Memory",
  "CD8_Memory_GZMK",
  "gdT", "Treg",
  "CD4_CD8_DP"
)

# ============================================================
# 3. 获取注释，仅统一绘图名称
# ============================================================
meta <- as.data.frame(getCellColData(proj))



cell_types_plot <- as.character(meta$TB_Subtype)
names(cell_types_plot) <- rownames(meta)


cells_use <- names(cell_types_plot)[
  !is.na(cell_types_plot) &
    cell_types_plot %in% celltype_levels
]

if (!length(cells_use)) stop("没有可用于绘图的细胞。")

print(table(factor(
  cell_types_plot[cells_use],
  levels = celltype_levels
)))

# ============================================================
# 4. 提取已归一化的 RNA 表达矩阵
# ============================================================
if (!"GeneExpressionMatrix" %in% getAvailableMatrices(proj)) {
  stop("对象中不存在 GeneExpressionMatrix。")
}

proj_plot <- subsetCells(
  ArchRProj = proj,
  cellNames = cells_use
)

se <- getMatrixFromProject(
  ArchRProj = proj_plot,
  useMatrix = "GeneExpressionMatrix"
)

# ============================================================
# 5. 提取 marker
# ============================================================
marker_info <- data.frame(
  gene = unlist(all_marker, use.names = FALSE),
  marker_group = rep(names(all_marker), lengths(all_marker)),
  stringsAsFactors = FALSE
)

# 允许同一个基因在不同展示组出现
marker_info$marker_id <- paste0(
  "Marker_", seq_len(nrow(marker_info))
)

marker_genes <- unique(marker_info$gene)

rd <- as.data.frame(rowData(se))

candidate_cols <- intersect(
  c("name", "gene_name", "symbol", "geneName"),
  colnames(rd)
)

gene_candidates <- lapply(
  candidate_cols,
  function(k) as.character(rd[[k]])
)
names(gene_candidates) <- candidate_cols

if (!is.null(rownames(se))) {
  gene_candidates[["rownames"]] <- rownames(se)
}

if (!length(gene_candidates)) {
  stop("没有可识别的基因名称字段。")
}

match_counts <- vapply(
  gene_candidates,
  function(x) sum(marker_genes %in% x),
  integer(1)
)

gene_names <- gene_candidates[[which.max(match_counts)]]

idx <- which(
  !is.na(gene_names) & gene_names %in% marker_genes
)

if (!length(idx)) stop("没有匹配到 marker。")

# 仅 marker 小矩阵转为普通矩阵
mat <- as.matrix(assay(se)[idx, , drop = FALSE])

# 合并同名基因行
mat <- rowsum(
  mat,
  group = gene_names[idx],
  reorder = FALSE
)

genes_use <- marker_genes[marker_genes %in% rownames(mat)]
mat <- mat[genes_use, , drop = FALSE]

if (any(!is.finite(mat)) || any(mat < 0)) {
  stop("表达矩阵存在非有限值或负值。")
}

missing_genes <- setdiff(marker_genes, genes_use)

if (length(missing_genes)) {
  message("缺失基因：", paste(missing_genes, collapse = ", "))
}

write.csv(
  data.frame(Gene = missing_genes),
  file.path(OUT_DIR, "Missing_marker_genes.csv"),
  row.names = FALSE
)

marker_info <- marker_info[
  marker_info$gene %in% genes_use,
  ,
  drop = FALSE
]

# ============================================================
# 6. 对齐细胞注释
# ============================================================
if (!all(colnames(mat) %in% names(cell_types_plot))) {
  stop("表达矩阵细胞名与注释不匹配。")
}

cell_group <- unname(cell_types_plot[colnames(mat)])

valid_cells <- !is.na(cell_group) &
  cell_group %in% celltype_levels

mat <- mat[, valid_cells, drop = FALSE]
cell_group <- cell_group[valid_cells]

if (!ncol(mat)) stop("没有匹配的细胞。")

# ============================================================
# 7. 计算平均表达和表达比例
# ============================================================
dot_data <- do.call(
  rbind,
  lapply(celltype_levels, function(ct) {

    cells <- which(cell_group == ct)
    if (!length(cells)) return(NULL)

    m <- mat[, cells, drop = FALSE]

    data.frame(
      celltype = ct,
      gene = rownames(m),
      nCells = length(cells),
      avg_exp = rowMeans(m),
      pct_exp = 100 * rowMeans(m > 0),
      row.names = NULL
    )
  })
)

# 已归一化矩阵直接汇总；
# log1p 和 Z-score 仅用于颜色展示
dot_data$log_avg <- log1p(dot_data$avg_exp)

dot_data$avg_scaled <- ave(
  dot_data$log_avg,
  dot_data$gene,
  FUN = function(x) {
    s <- sd(x)
    if (!is.finite(s) || s == 0) {
      rep(0, length(x))
    } else {
      (x - mean(x)) / s
    }
  }
)

dot_data$avg_scaled_plot <- pmax(
  -2.5, pmin(2.5, dot_data$avg_scaled)
)

# ============================================================
# 8. 合并 marker 分组
# ============================================================
plot_df <- merge(
  dot_data,
  marker_info,
  by = "gene",
  sort = FALSE
)

plot_df$marker_id <- factor(
  plot_df$marker_id,
  levels = marker_info$marker_id
)

plot_df$marker_group <- factor(
  plot_df$marker_group,
  levels = names(all_marker)
)

present_celltypes <- celltype_levels[
  celltype_levels %in% dot_data$celltype
]

plot_df$celltype <- factor(
  plot_df$celltype,
  levels = rev(present_celltypes)
)

# ============================================================
# 9. 绘图
# ============================================================
p_all_marker_dotplot <- ggplot(
  plot_df,
  aes(x = marker_id, y = celltype)
) +
  geom_point(
    aes(size = pct_exp, colour = avg_scaled_plot)
  ) +
  facet_grid(
    . ~ marker_group,
    scales = "free_x",
    space = "free_x"
  ) +
  scale_x_discrete(
    labels = function(x) {
      marker_info$gene[
        match(as.character(x), marker_info$marker_id)
      ]
    }
  ) +
  scale_colour_gradient2(
    low = "#3B6FB6",
    mid = "white",
    high = "#C83E3E",
    midpoint = 0,
    limits = c(-2.5, 2.5),
    name = "Scaled expression"
  ) +
  scale_size(
    range = c(0, 6),
    limits = c(0, 100),
    breaks = c(0, 25, 50, 75, 100),
    name = "Expressing cells (%)"
  ) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    panel.spacing.x = grid::unit(0.3, "lines"),
    panel.border = element_rect(
      colour = "#555555",
      linewidth = 0.5
    ),
    strip.background = element_rect(
      fill = "grey95",
      colour = "#555555"
    ),
    strip.text.x = element_text(
      face = "bold",
      size = 10
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 9
    ),
    axis.text.y = element_text(size = 11),
    legend.position = "right"
  )

print(p_all_marker_dotplot)

# ============================================================
# 10. 保存
# ============================================================
ggsave(
  file.path(OUT_DIR, "All_marker_DotPlot.pdf"),
  plot = p_all_marker_dotplot,
  width = 20,
  height = 9,
  units = "in"
)

ggsave(
  file.path(OUT_DIR, "All_marker_DotPlot.png"),
  plot = p_all_marker_dotplot,
  width = 20,
  height = 9,
  units = "in",
  dpi = 300,
  bg = "white"
)

write.csv(
  plot_df,
  file.path(OUT_DIR, "All_marker_DotPlot_values.csv"),
  row.names = FALSE
)

write.csv(
  marker_info,
  file.path(OUT_DIR, "Marker_groups.csv"),
  row.names = FALSE
)

cat("输出目录：\n", OUT_DIR, "\n")

