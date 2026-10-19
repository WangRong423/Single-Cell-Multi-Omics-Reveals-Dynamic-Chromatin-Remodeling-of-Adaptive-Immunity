
STAGE <- 1L
REBUILD_ROOT <- "/database/wangrong/Results/0712_ATAC+RNA/atac_output_repaired"
TB_OUTDIR <- file.path(REBUILD_ROOT, "TB")
dir.create(TB_OUTDIR, recursive = TRUE, showWarnings = FALSE)
library(dplyr)
library(Matrix)
library(SummarizedExperiment)
set.seed(123)
annotation_from_csv <- function(project, cluster_col, file, label_col) {
  clusters <- as.character(getCellColData(project, select = cluster_col)[[1]])
  if (!file.exists(file)) {
    template <- data.frame(Cluster = sort(unique(clusters)), Label = "")
    names(template)[2] <- label_col
    write.csv(template, file, row.names = FALSE)
    stop("请根据 marker 填写注释表后继续：", file)
  }
  tab <- read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("Cluster", label_col) %in% names(tab))) stop("注释表列名错误")
  if (anyDuplicated(tab$Cluster)) stop("注释表 Cluster 重复")
  labels <- trimws(as.character(tab[[label_col]][match(clusters, tab$Cluster)]))
  if (anyNA(labels) || any(labels == "")) stop("注释表有遗漏/空标签：", file)
  addCellColData(project, data = labels, cells = getCellNames(project),
                 name = label_col, force = TRUE)
}
if (!STAGE %in% 1:3) stop("STAGE 必须为 1、2 或 3")
library(ArchR)
library(pheatmap)
library(Seurat) 
library(ggplot2)
library(BSgenome.Hsapiens.UCSC.hg38)
library(chromVARmotifs)

addArchRThreads(threads = max(1L, min(16L, parallel::detectCores())))
addArchRGenome("hg38")

if (STAGE == 1L) {
batch_list=list("M1-1","M1-2","M1-3","M1-4","M1-5","M1-6","M1-7","M1-8","M1-9","M1-10",
                "M2-1","M2-2","M2-3","M2-4","M2-5","M2-6","M2-7","M2-8","M2-9","M2-10",
                "M3-1","M3-2","M3-3","M3-4","M3-5","M3-6","M3-7","M3-8","M3-9","M3-10",
                "M5-1","M5-2","M5-3","M5-4","M5-5","M5-6","M5-7","M5-8","M5-9","M5-10")

name_list=c("M1-1","M1-2","M1-3","M1-4","M1-5","M1-6","M1-7","M1-8","M1-9","M1-10",
            "M2-1","M2-2","M2-3","M2-4","M2-5","M2-6","M2-7","M2-8","M2-9","M2-10",
            "M3-1","M3-2","M3-3","M3-4","M3-5","M3-6","M3-7","M3-8","M3-9","M3-10",
            "M5-1","M5-2","M5-3","M5-4","M5-5","M5-6","M5-7","M5-8","M5-9","M5-10")

atac_input<-c()
for(i in 1:length(batch_list)){
  atac_input<-c(atac_input,paste0('/database/wangrong/Results/0712_ATAC+RNA/',batch_list[[i]],'/outs/atac_fragments.tsv.gz'))
}
stopifnot(all(file.exists(atac_input)))
setwd(REBUILD_ROOT)
ArrowFiles <- createArrowFiles(inputFiles = atac_input, sampleNames = name_list,
  minTSS = 4, minFrags = 1000, maxFrags = 100000,
  addTileMat = TRUE, addGeneScoreMat = TRUE)

rna_input<-c()
for(i in 1:length(batch_list)){
  rna_input<-c(rna_input,paste0('/database/wangrong/Results/0712_ATAC+RNA/',batch_list[[i]],'/outs/filtered_feature_bc_matrix.h5'))
}
stopifnot(all(file.exists(rna_input)))
seRNA <- import10xFeatureMatrix(input = rna_input,names = name_list,featureType = "Gene Expression")

projpbmc <- ArchRProject(ArrowFiles = ArrowFiles,outputDirectory = file.path(REBUILD_ROOT, "Project"), copyArrows = TRUE)

# 在 addGeneExpressionMatrix 前保存/计算原始 RNA counts 指标。
raw_counts <- SummarizedExperiment::assay(seRNA)
count_values <- if (inherits(raw_counts, "sparseMatrix")) raw_counts@x else as.vector(raw_counts)
if (any(!is.finite(count_values)) || any(count_values < 0) ||
    any(abs(count_values - round(count_values)) > 1e-8)) stop("RNA 非原始计数")
rd <- as.data.frame(SummarizedExperiment::rowData(seRNA))
symbol_col <- intersect(c("name", "symbol", "gene_name"), names(rd))
if (!length(symbol_col)) stop("seRNA 缺少基因符号，无法可靠计算线粒体比例")
symbols <- as.character(rd[[symbol_col[1]]])
mito <- grepl("^MT-", symbols, ignore.case = TRUE)
if (!any(mito)) stop("未匹配到 MT- 基因；请核查基因符号/参考基因组")
rna_qc <- data.frame(Cell = colnames(raw_counts),
  RNA_nUMI_raw = Matrix::colSums(raw_counts),
  RNA_nGenes_raw = Matrix::colSums(raw_counts > 0),
  RNA_mitoRatio_raw = Matrix::colSums(raw_counts[mito,,drop=FALSE]) /
    pmax(Matrix::colSums(raw_counts), 1))
stopifnot(!anyDuplicated(rna_qc$Cell))
saveRDS(seRNA, file.path(REBUILD_ROOT, "seRNA_raw_counts.rds"), compress = FALSE)

new_ArrowFiles<-addGeneExpressionMatrix(input=projpbmc, seRNA=seRNA)
proj <- new_ArrowFiles
idx <- match(getCellNames(proj), rna_qc$Cell)
if (anyNA(idx)) stop("RNA/ATAC sample#barcode 未完整匹配；请检查样本名称与条形码后缀")
for (nm in setdiff(names(rna_qc), "Cell")) {
  proj <- addCellColData(proj, data = rna_qc[[nm]][idx],
    cells = getCellNames(proj), name = nm, force = TRUE)
}
sample_meta <- data.frame(Sample = name_list,
  Donor = sub("-.*$", "", name_list),
  Timepoint = paste0("day", rep(c(0,1,3,6,14,30,31,33,36,44), 4)))
# 上述后缀 1:10 -> 时间点对应关系须与原始样本记录核对。
write.csv(sample_meta, file.path(REBUILD_ROOT, "sample_metadata.csv"), row.names = FALSE)
for (nm in c("Donor", "Timepoint")) {
  proj <- addCellColData(proj, data = sample_meta[[nm]][match(as.character(proj$Sample), sample_meta$Sample)],
    cells = getCellNames(proj), name = nm, force = TRUE)
}

minFeatures <- 200
maxFeatures <- 6000

minUMI <- 1000
maxUMI <- 50000

# 比例值，不是百分数
maxMT <- 0.20

minFrags <- 1000
maxFrags <- 100000
minTSS <- 10

qc_meta <- as.data.frame(
  getCellColData(
    proj,
    select = c(
      "Sample",
      "RNA_nUMI_raw",
      "RNA_nGenes_raw",
      "RNA_mitoRatio_raw",
      "nFrags",
      "TSSEnrichment"
    )
  )
)

qc_meta <- qc_meta[
  match(
    getCellNames(proj),
    rownames(qc_meta)
  ),
  ,
  drop = FALSE
]

keep_joint <- with(
  qc_meta,

  is.finite(RNA_nUMI_raw) &
    RNA_nUMI_raw >= minUMI &
    RNA_nUMI_raw <= maxUMI &

  is.finite(RNA_nGenes_raw) &
    RNA_nGenes_raw >= minFeatures &
    RNA_nGenes_raw <= maxFeatures &

  is.finite(RNA_mitoRatio_raw) &
    RNA_mitoRatio_raw <= maxMT &

  is.finite(nFrags) &
    nFrags >= minFrags &
    nFrags <= maxFrags &

  is.finite(TSSEnrichment) &
    TSSEnrichment >= minTSS
)

keep_joint[
  is.na(keep_joint)
] <- FALSE

qc_retention <- tibble(
  Cell = getCellNames(proj),
  Sample = as.character(proj$Sample),
  Pass_QC = keep_joint
) %>%
  group_by(Sample) %>%
  summarise(
    n_before = n(),
    n_after = sum(Pass_QC),
    retention = mean(Pass_QC),
    .groups = "drop"
  )

print(
  qc_retention,
  n = Inf
)

stopifnot(
  nrow(qc_retention) == 40,
  all(qc_retention$n_after > 0)
)

write.csv(qc_retention, file.path(REBUILD_ROOT, "joint_QC_retention.csv"), row.names = FALSE)
proj <- subsetCells(proj, cellNames = getCellNames(proj)[keep_joint])
proj <- addDoubletScores(
  input = proj,
  useMatrix = "TileMatrix",
  k = 10,
  nTrials = 5,
  dimsToUse = 1:30,
  LSIMethod = 1,
  knnMethod = "UMAP",
  force = TRUE,
  threads = getArchRThreads()
)

proj <- filterDoublets(
  ArchRProj = proj
)

proj <- addIterativeLSI(
  ArchRProj = proj,
  useMatrix = "TileMatrix",
  depthCol = "nFrags",
  name = "LSI_ATAC",
  iterations = 2,
  clusterParams = list(
    resolution = 0.2,
    sampleCells = 10000,
    maxClusters = 30,
    n.start = 10
  ),
  varFeatures = 25000,
  dimsToUse = 1:30,
  firstSelection = "top",
  binarize = TRUE,
  saveIterations = FALSE,
  force = TRUE
)

proj <- addIterativeLSI(
  ArchRProj = proj,
  useMatrix = "GeneExpressionMatrix",
  depthCol = "Gex_nUMI",
  name = "LSI_RNA",
  iterations = 2,
  clusterParams = list(
    resolution = 0.2,
    sampleCells = 10000,
    maxClusters = 30,
    n.start = 10
  ),
  varFeatures = 2500,
  dimsToUse = 1:30,
  firstSelection = "variable",
  binarize = FALSE,
  saveIterations = FALSE,
  force = TRUE
)

proj <- addCombinedDims(
  ArchRProj = proj,
  reducedDims = c(
    "LSI_ATAC",
    "LSI_RNA"
  ),
  name = "LSI_Combined",
  dimsToUse = 1:30,
  dimWeights = c(1, 1),
  corCutOff = 0.75
)


umap_sample_cells <- min(
  200000,
  length(getCellNames(proj))
)

proj <- addUMAP(
  ArchRProj = proj,
  reducedDims = "LSI_ATAC",
  name = "UMAP_ATAC",
  nNeighbors = 40,
  minDist = 0.5,
  metric = "cosine",
  sampleCells = umap_sample_cells,
  force = TRUE,
  threads = 1
)

proj <- addUMAP(
  ArchRProj = proj,
  reducedDims = "LSI_RNA",
  name = "UMAP_RNA",
  nNeighbors = 40,
  minDist = 0.5,
  metric = "cosine",
  sampleCells = umap_sample_cells,
  force = TRUE,
  threads = 1
)

proj <- addUMAP(
  ArchRProj = proj,
  reducedDims = "LSI_Combined",
  name = "UMAP_Combined",
  nNeighbors = 40,
  minDist = 0.5,
  metric = "cosine",
  sampleCells = umap_sample_cells,
  force = TRUE,
  threads = 1
)

proj <- addClusters(
  input = proj,
  reducedDims = "LSI_Combined",
  name = "Clusters",
  resolution = 0.4,
  method = "Seurat",
  sampleCells = umap_sample_cells,
  seed = 123,
  maxClusters = 50,
  force = TRUE
)

p1 <- plotEmbedding(
  ArchRProj = proj,
  name = "Clusters",
  embedding = "UMAP_ATAC",
  colorBy = "cellColData",
  size = 0.1,
  labelAsFactors = FALSE,
  labelMeans = TRUE
)

p2 <- plotEmbedding(
  ArchRProj = proj,
  name = "Clusters",
  embedding = "UMAP_RNA",
  colorBy = "cellColData",
  size = 0.1,
  labelAsFactors = FALSE,
  labelMeans = TRUE
)

p3 <- plotEmbedding(
  ArchRProj = proj,
  name = "Clusters",
  embedding = "UMAP_Combined",
  colorBy = "cellColData",
  size = 0.1,
  labelAsFactors = FALSE,
  labelMeans = TRUE
)

plotPDF(
  plotList = list(
    p1,
    p2,
    p3
  ),
  name = "UMAP_ATAC_RNA_Combined_Clusters",
  ArchRProj = proj,
  addDOC = FALSE,
  width = 6,
  height = 5
)

proj_anno <- addHarmony(
  ArchRProj = proj,
  reducedDims = "LSI_ATAC",
  dimsToUse = 1:30,
  name = "Harmony_ATAC",
  groupBy = "Donor",
  corCutOff = 0.75,
  force = TRUE
)

proj_anno <- addHarmony(
  ArchRProj = proj_anno,
  reducedDims = "LSI_RNA",
  dimsToUse = 1:30,
  name = "Harmony_RNA",
  groupBy = "Donor",
  corCutOff = 0.75,
  force = TRUE
)

names(proj_anno@reducedDims)

dim(
  getReducedDims(
    proj_anno,
    reducedDims = "Harmony_ATAC"
  )
)

dim(
  getReducedDims(
    proj_anno,
    reducedDims = "Harmony_RNA"
  )
)

proj_anno <- addCombinedDims(
  ArchRProj = proj_anno,
  reducedDims = c(
    "Harmony_ATAC",
    "Harmony_RNA"
  ),
  name = "Harmony_Combined",
  dimsToUse = 1:30,
  dimWeights = c(1, 1),
  corCutOff = 0.75
)

proj_anno <- addUMAP(
  ArchRProj = proj_anno,
  reducedDims = "Harmony_ATAC",
  name = "UMAP_ATAC_Harmony",
  nNeighbors = 40,
  minDist = 0.5,
  metric = "cosine",
  force = TRUE
)

proj_anno <- addUMAP(
  ArchRProj = proj_anno,
  reducedDims = "Harmony_RNA",
  name = "UMAP_RNA_Harmony",
  nNeighbors = 40,
  minDist = 0.5,
  metric = "cosine",
  force = TRUE
)

proj_anno <- addUMAP(
  ArchRProj = proj_anno,
  reducedDims = "Harmony_Combined",
  name = "UMAP_Combined_Harmony",
  nNeighbors = 40,
  minDist = 0.5,
  metric = "cosine",
  force = TRUE
)

proj_anno <- addClusters(
  input = proj_anno,
  reducedDims = "Harmony_Combined",
  name = "Clusters_Harmony",
  resolution = 0.4,
  method = "Seurat",
  force = TRUE
)

p_donor <- plotEmbedding(
  ArchRProj = proj_anno,
  colorBy = "cellColData",
  name = "Donor",
  embedding = "UMAP_Combined_Harmony",
  size = 0.15,
  labelMeans = FALSE
)

p_time <- plotEmbedding(
  ArchRProj = proj_anno,
  colorBy = "cellColData",
  name = "Timepoint",
  embedding = "UMAP_Combined_Harmony",
  size = 0.15,
  labelMeans = FALSE
)

p_cluster <- plotEmbedding(
  ArchRProj = proj_anno,
  colorBy = "cellColData",
  name = "Clusters_Harmony",
  embedding = "UMAP_Combined_Harmony",
  size = 0.15,
  labelMeans = TRUE
)

plotPDF(
  plotList = list(
    p_donor,
    p_time,
    p_cluster
  ),
  name = "Harmony_Combined_batch_QC",
  ArchRProj = proj_anno,
  addDOC = FALSE,
  width = 6,
  height = 5
)

markers_rna <- getMarkerFeatures(
  ArchRProj = proj_anno,
  useMatrix = "GeneExpressionMatrix",
  groupBy = "Clusters_Harmony",
  bias = c(
    "Gex_nUMI",
    "TSSEnrichment"
  ),
  testMethod = "wilcoxon"
)

saveRDS(
  proj_anno,
  file = file.path(REBUILD_ROOT, "proj_anno_QC_doublet_Harmony.rds"),
  compress = FALSE
)

broad_markers <- list(

  T = c(
    "CD3D", "CD3E", "CD3G",
    "TRAC", "CD247", "LTB"
  ),

  B = c(
    "CD79A", "CD79B", "MS4A1",
    "CD74", "CD37", "HLA-DRA",
    # Plasma也归入B大类
    "MZB1", "TNFRSF17", "JCHAIN"
  ),

  NK = c(
    "NKG7", "GNLY", "KLRD1",
    "FCGR3A", "PRF1", "GZMB"
  ),

  Myeloid = c(
    "LYZ", "LST1", "TYROBP",
    "FCER1G", "CTSS", "CTSD",
    "LILRB1", "S100A8", "S100A9"
  )
)

broad_marker_genes <- unique(
  unlist(broad_markers)
)

markersRNA_broad <- getMarkerFeatures(
  ArchRProj = proj_anno,
  useMatrix = "GeneExpressionMatrix",
  groupBy = "Clusters_Harmony",
  bias = c("Gex_nUMI"),
  testMethod = "wilcoxon",
  maxCells = 1000,
  threads = getArchRThreads()
)

get_marker_ids <- function(se_marker, genes) {

  rd <- as.data.frame(
    SummarizedExperiment::rowData(se_marker)
  )

  candidate_cols <- intersect(
    c(
      "name",
      "symbol",
      "gene_name",
      "gene",
      "Gene"
    ),
    colnames(rd)
  )

  if (length(candidate_cols) > 0) {
    symbols <- as.character(
      rd[[candidate_cols[1]]]
    )
  } else {
    symbols <- rownames(se_marker)
  }

  idx <- match(
    toupper(genes),
    toupper(symbols)
  )

  ids <- symbols[idx[!is.na(idx)]]

  unique(ids)
}

broad_marker_ids <- get_marker_ids(
  markersRNA_broad,
  broad_marker_genes
)

length(broad_marker_ids)

heatmap_broad <- plotMarkerHeatmap(
  seMarker = markersRNA_broad,

  subsetMarkers = broad_marker_ids,

  labelMarkers = broad_marker_ids,

  labelRows = TRUE,

  log2Norm = TRUE,

  scaleRows = TRUE,

  clusterCols = FALSE,

  transpose = TRUE,

  limits = c(-2, 2)
)

plotPDF(
  heatmap_broad,
  name = "Broad_celltype_RNA_marker_heatmap",
  ArchRProj = proj_anno,
  addDOC = FALSE,
  width = 12,
  height = 8
)

broad_umap_genes <- c(
  "CD3D",
  "MS4A1",
  "NKG7",
  "LST1",
  "MZB1"
)

broad_gene_plots <- lapply(
  broad_umap_genes,
  function(gene_i) {

    plotEmbedding(
      ArchRProj = proj_anno,
      embedding = "UMAP_Combined_Harmony",
      colorBy = "GeneExpressionMatrix",
      name = gene_i,
      log2Norm = TRUE,
      imputeWeights = NULL,
      quantCut = c(0.01, 0.99),
      size = 0.1,
      plotAs = "points"
    )
  }
)

names(broad_gene_plots) <- broad_umap_genes

plotPDF(
  plotList = broad_gene_plots,
  name = "Broad_celltype_marker_UMAP",
  ArchRProj = proj_anno,
  addDOC = FALSE,
  width = 5,
  height = 5
)

saveRDS(proj_anno, file.path(REBUILD_ROOT, "stage1_project.rds"), compress = FALSE)
writeLines(capture.output(sessionInfo()), file.path(REBUILD_ROOT, "sessionInfo.txt"))
message("阶段1完成；填写 broad_annotation.csv 后运行 STAGE=2")
if (!file.exists(file.path(REBUILD_ROOT, "broad_annotation.csv"))) {
  write.csv(data.frame(Cluster = sort(unique(as.character(proj_anno$Clusters_Harmony))),
    BroadCellType = ""), file.path(REBUILD_ROOT, "broad_annotation.csv"), row.names = FALSE)
}
}
if (STAGE == 2L) {
proj_anno <- readRDS(file.path(REBUILD_ROOT, "stage1_project.rds"))
proj_anno <- annotation_from_csv(proj_anno, "Clusters_Harmony",
  file.path(REBUILD_ROOT, "broad_annotation.csv"), "BroadCellType")
# 允许 T、NK 分开注释；T_NK 需在下游排除 NK。
library(ArchR)
library(dplyr)
library(ggplot2)

table(
  proj_anno$BroadCellType,
  useNA = "ifany"
)

TB_cells <- getCellNames(proj_anno)[
  as.character(proj_anno$BroadCellType) %in%
    c("B", "T", "T_NK")
]

length(TB_cells)

proj_TB <- subsetCells(
  ArchRProj = proj_anno,
  cellNames = TB_cells
)

table(proj_TB$BroadCellType)
table(proj_TB$Sample)


n_TB <- length(
  getCellNames(proj_TB)
)

sample_cells_TB <- min(
  20000,
  n_TB
)

proj_TB <- addIterativeLSI(
  ArchRProj = proj_TB,
  useMatrix = "TileMatrix",
  depthCol = "nFrags",
  name = "TB_LSI_ATAC",
  iterations = 2,

  clusterParams = list(
    resolution = 0.2,
    sampleCells = sample_cells_TB,
    n.start = 10
  ),

  varFeatures = 25000,
  dimsToUse = 1:30,
  firstSelection = "top",
  binarize = TRUE,
  saveIterations = FALSE,
  force = TRUE
)

proj_TB <- addIterativeLSI(
  ArchRProj = proj_TB,
  useMatrix = "GeneExpressionMatrix",
  depthCol = "Gex_nUMI",
  name = "TB_LSI_RNA",
  iterations = 2,

  clusterParams = list(
    resolution = 0.2,
    sampleCells = sample_cells_TB,
    n.start = 10
  ),

  varFeatures = 3000,
  dimsToUse = 1:30,
  firstSelection = "variable",
  binarize = FALSE,
  saveIterations = FALSE,
  force = TRUE
)

proj_TB <- addHarmony(
  ArchRProj = proj_TB,
  reducedDims = "TB_LSI_ATAC",
  dimsToUse = 1:30,
  name = "TB_Harmony_ATAC",
  groupBy = "Donor",
  corCutOff = 0.75,
  force = TRUE
)

proj_TB <- addHarmony(
  ArchRProj = proj_TB,
  reducedDims = "TB_LSI_RNA",
  dimsToUse = 1:30,
  name = "TB_Harmony_RNA",
  groupBy = "Donor",
  corCutOff = 0.75,
  force = TRUE
)

proj_TB <- addCombinedDims(
  ArchRProj = proj_TB,
  reducedDims = c(
    "TB_Harmony_ATAC",
    "TB_Harmony_RNA"
  ),
  name = "TB_Harmony_Combined"
)

names(proj_TB@reducedDims)

dim(
  getReducedDims(
    ArchRProj = proj_TB,
    reducedDims = "TB_Harmony_Combined"
  )
)

proj_TB <- addUMAP(
  ArchRProj = proj_TB,
  reducedDims = "TB_Harmony_Combined",
  name = "TB_UMAP",
  nNeighbors = 40,
  minDist = 0.4,
  metric = "cosine",
  seed = 123,
  force = TRUE
)

proj_TB <- addClusters(
  input = proj_TB,
  reducedDims = "TB_Harmony_Combined",
  name = "TB_Clusters",
  resolution = 0.8,
  method = "Seurat",
  dimsToUse = 1:30,
  seed = 123,
  force = TRUE
)

table(proj_TB$TB_Clusters)

p_TB_cluster <- plotEmbedding(
  ArchRProj = proj_TB,
  embedding = "TB_UMAP",
  colorBy = "cellColData",
  name = "TB_Clusters",
  size = 0.1,
  labelMeans = TRUE,
  labelAsFactors = FALSE
)

p_TB_broad <- plotEmbedding(
  ArchRProj = proj_TB,
  embedding = "TB_UMAP",
  colorBy = "cellColData",
  name = "BroadCellType",
  size = 0.1,
  labelMeans = TRUE,
  labelAsFactors = FALSE
)

p_TB_donor <- plotEmbedding(
  ArchRProj = proj_TB,
  embedding = "TB_UMAP",
  colorBy = "cellColData",
  name = "Donor",
  size = 0.1,
  labelMeans = FALSE
)

plotPDF(
  plotList = list(
    p_TB_cluster,
    p_TB_broad,
    p_TB_donor
  ),
  name = "TB_joint_subclustering_UMAP",
  ArchRProj = proj_TB,
  addDOC = FALSE,
  width = 7,
  height = 6
)

markers_TB <- getMarkerFeatures(proj_TB, useMatrix = "GeneExpressionMatrix",
  groupBy = "TB_Clusters", bias = "Gex_nUMI", testMethod = "wilcoxon", maxCells = 1000)
saveRDS(markers_TB, file.path(TB_OUTDIR, "TB_RNA_markers.rds"))
marker_list <- getMarkers(markers_TB, cutOff = "FDR <= 0.05 & Log2FC >= 0.25")
for (nm in names(marker_list)) write.csv(as.data.frame(marker_list[[nm]]),
  file.path(TB_OUTDIR, paste0("markers_", nm, ".csv")), row.names = FALSE)
saveRDS(proj_TB, file.path(TB_OUTDIR, "stage2_project.rds"), compress = FALSE)
saveRDS(proj_anno, file.path(REBUILD_ROOT, "stage2_parent.rds"), compress = FALSE)
if (!file.exists(file.path(TB_OUTDIR, "TB_annotation.csv"))) {
  write.csv(data.frame(Cluster = sort(unique(as.character(proj_TB$TB_Clusters))), TB_Subtype = ""),
    file.path(TB_OUTDIR, "TB_annotation.csv"), row.names = FALSE)
}
message("阶段2完成；填写 TB_annotation.csv，NK 使用 NK 标签，然后运行 STAGE=3")
}
if (STAGE == 3L) {
proj_TB <- readRDS(file.path(TB_OUTDIR, "stage2_project.rds"))
proj_anno <- readRDS(file.path(REBUILD_ROOT, "stage2_parent.rds"))
proj_TB <- annotation_from_csv(proj_TB, "TB_Clusters", file.path(TB_OUTDIR, "TB_annotation.csv"), "TB_Subtype")
p_TB_subtype <- plotEmbedding(
  ArchRProj = proj_TB,
  embedding = "TB_UMAP",
  colorBy = "cellColData",
  name = "TB_Subtype",
  size = 0.1,
  labelMeans = TRUE,
  labelAsFactors = FALSE
)

p_TB_subtype

plotPDF(
  plotList = list(
    p_TB_subtype
  ),
  name = "TB_joint_subtype_annotation_UMAP",
  ArchRProj = proj_TB,
  addDOC = FALSE,
  width = 10,
  height = 7
)
final_celltype <- as.character(
  proj_anno$BroadCellType
)

names(final_celltype) <- getCellNames(
  proj_anno
)

TB_labels <- setNames(
  as.character(
    proj_TB$TB_Subtype
  ),
  getCellNames(proj_TB)
)

final_celltype[
  names(TB_labels)
] <- TB_labels

proj_anno <- addCellColData(
  ArchRProj = proj_anno,
  data = unname(final_celltype),
  cells = getCellNames(proj_anno),
  name = "FinalCellType",
  force = TRUE
)


saveRDS(
  proj_TB,
  file.path(
    TB_OUTDIR,
    "proj_TB_joint_subclustered.rds"
  ),
  compress = FALSE
)

saveRDS(
  proj_anno,
  file.path(
    TB_OUTDIR,
    "proj_anno_with_TB_subtypes.rds"
  ),
  compress = FALSE
)


exclude_labels <- c("NK", "Myeloid", "Unassigned", "Doublet", "Other")
keep_tb <- !as.character(proj_TB$TB_Subtype) %in% exclude_labels
if (!any(keep_tb)) stop("没有可用于下游的 T/B 细胞")
proj_TB_only <- subsetCells(proj_TB, cellNames = getCellNames(proj_TB)[keep_tb])

proj_TB_down <- subsetArchRProject(proj_TB_only,
  cells = getCellNames(proj_TB_only), outputDirectory = file.path(TB_OUTDIR, "DownstreamProject"),
  dropCells = TRUE, force = FALSE)
write.csv(as.data.frame(table(Donor = proj_TB_down$Donor,
  Timepoint = proj_TB_down$Timepoint, Subtype = proj_TB_down$TB_Subtype)),
  file.path(TB_OUTDIR, "donor_time_subtype_counts.csv"), row.names = FALSE)
coverage_counts <- table(proj_TB_down$TB_Subtype, proj_TB_down$Donor)
write.csv(as.data.frame.matrix(coverage_counts), file.path(TB_OUTDIR, "coverage_cells_per_donor.csv"))
if (any(rowSums(coverage_counts >= 40) < 2)) {
  stop("部分亚型少于2名供者各有40个细胞；请检查 coverage_cells_per_donor.csv，并调整分组或阈值")
}
proj_TB_down <- addGroupCoverages(
  ArchRProj = proj_TB_down,
  groupBy = "TB_Subtype",
  sampleLabels = "Donor",
  useLabels = TRUE,
  minCells = 40,
  maxCells = 500,
  minReplicates = 2,
  maxReplicates = 5,
  sampleRatio = 0.8,
  force = TRUE,
  threads = getArchRThreads()
)

pathToMacs2 <- findMacs2()

if (
  is.null(pathToMacs2) ||
  length(pathToMacs2) == 0 ||
  !file.exists(pathToMacs2)
) {
  pathToMacs2 <- "/home/wangrong/.local/bin/macs2"
}

stopifnot(file.exists(pathToMacs2))

pathToMacs2

proj_TB_down <- addReproduciblePeakSet(
  ArchRProj = proj_TB_down,
  groupBy = "TB_Subtype",
  peakMethod = "Macs2",
  pathToMacs2 = pathToMacs2,
  reproducibility = "2",
  peaksPerCell = 500,
  maxPeaks = 150000,
  minCells = 25,
  excludeChr = c("chrM", "chrY"),
  force = TRUE
)
proj_TB_down <- addPeakMatrix(
  ArchRProj = proj_TB_down,
  force = TRUE,
  threads = getArchRThreads()
)

proj_TB_down <- addPeak2GeneLinks(
  ArchRProj = proj_TB_down,
  reducedDims = "TB_Harmony_Combined",
  useMatrix = "GeneExpressionMatrix",
  dimsToUse = 1:30,
  corCutOff = 0.75,
  k = 100,
  knnIteration = 500,
  overlapCutoff = 0.8,
  maxDist = 250000,
  scaleTo = 10000,
  log2Norm = TRUE,
  addEmpiricalPval = FALSE,
  seed = 123,
  threads = max(
    1,
    floor(getArchRThreads() / 2)
  )
)

p2g_loops_TB <- getPeak2GeneLinks(
  ArchRProj = proj_TB_down,
  corCutOff = 0.45,
  FDRCutOff = 1e-4,
  varCutOffATAC = 0.25,
  varCutOffRNA = 0.25,
  resolution = 1,
  returnLoops = TRUE
)

length(p2g_loops_TB)

saveRDS(
  p2g_loops_TB,
  file.path(TB_OUTDIR, "TB_peak_to_gene_links_loops.rds")
)

write.csv(
  as.data.frame(getPeak2GeneLinks(proj_TB_down, corCutOff = 0.45, FDRCutOff = 1e-4,
    varCutOffATAC = 0.25, varCutOffRNA = 0.25, returnLoops = FALSE)),
  file.path(TB_OUTDIR, "TB_peak_to_gene_links.csv"),
  row.names = FALSE
)

p_p2g_heatmap <- plotPeak2GeneHeatmap(
  ArchRProj = proj_TB_down,
  corCutOff = 0.45,
  FDRCutOff = 1e-4,
  varCutOffATAC = 0.25,
  varCutOffRNA = 0.25,
  k = 25,
  nPlot = 25000,
  groupBy = "TB_Subtype"
)

plotPDF(
  p_p2g_heatmap,
  name = "TB_peak_to_gene_heatmap",
  ArchRProj = proj_TB_down,
  addDOC = FALSE,
  width = 10,
  height = 8
)

proj_TB_down <- addMotifAnnotations(
  ArchRProj = proj_TB_down,
  motifSet = "cisbp",
  annoName = "Motif",
  species = "Homo sapiens",
  force = TRUE
)

proj_TB_down <- addBgdPeaks(
  ArchRProj = proj_TB_down,
  method = "chromVAR",
  force = TRUE
)

proj_TB_down <- addDeviationsMatrix(
  ArchRProj = proj_TB_down,
  peakAnnotation = "Motif",
  matrixName = "MotifMatrix",
  out = c("z", "deviations"),
  force = TRUE,
  threads = getArchRThreads()
)

markerMotifs_TB <- getMarkerFeatures(
  ArchRProj = proj_TB_down,
  useMatrix = "MotifMatrix",
  groupBy = "TB_Subtype",
  useSeqnames = "z",
  bias = c(
    "TSSEnrichment",
    "log10(nFrags)"
  ),
  testMethod = "wilcoxon",
  maxCells = 1000,
  threads = getArchRThreads()
)

motif_heatmap_TB <- plotMarkerHeatmap(
  seMarker = markerMotifs_TB,
  cutOff = "FDR <= 0.05 & MeanDiff >= 0.5",
  transpose = TRUE,
  nLabel = 10,
  log2Norm = FALSE
)

plotPDF(
  motif_heatmap_TB,
  name = "TB_subtype_chromVAR_motif_heatmap",
  ArchRProj = proj_TB_down,
  addDOC = FALSE,
  width = 12,
  height = 8
)

saveRDS(
  markerMotifs_TB,
  file.path(TB_OUTDIR, "TB_marker_motif_activity.rds")
)

saveRDS(
  proj_TB_only,
  file.path(
    TB_OUTDIR,
    "proj_TB_only_annotated.rds"
  ),
  compress = FALSE
)

saveRDS(
  proj_TB_down,
  file.path(
    TB_OUTDIR,
    "proj_TB_callpeak_P2G_motif_complete.rds"
  ),
  compress = FALSE
)
writeLines(capture.output(sessionInfo()), file.path(TB_OUTDIR, "sessionInfo.txt"))
}

p_umap <- plotEmbedding(
    ArchRProj = TB_proj,
     colorBy = "cellColData",
     name = "TB_Subtype",
     embedding = "TB_UMAP_TBonly",
     pal = new_palette,
     size = 0.3)

print(p_umap)
save(p_umap,file="/public/users/wangrong/home/wangrong/results/T_B/TB_final_downstream/obj/Fig1A.rds")