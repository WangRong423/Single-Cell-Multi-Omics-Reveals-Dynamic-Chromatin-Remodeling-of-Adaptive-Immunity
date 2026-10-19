# Reanalyse RNA with summed RAW UMI counts, donor fixed effects and edgeR QL.
# Run this first; then run 02_RNA_prime_booster_matched.R.
BASE_TB_DIR <- "/public/users/wangrong/home/wangrong/results/T_B/TB_final_downstream"
FULL_PROJECT_RDS <- "" # run through 00_run_TB_proj.R
RAW_COUNTS_RDS <- if(exists("TB_RAW_COUNTS_RDS",envir=.GlobalEnv)) TB_RAW_COUNTS_RDS else ""
# Alternatively keep ORIGINAL seRNA (SummarizedExperiment of raw counts) in R session.
# Matrix barcodes MUST match full ArchR cell names; do not invent/rewrite cell IDs.
OUTDIR <- file.path(BASE_TB_DIR,"Figure2_DEG","TB_proj_reanalysis")
MIN_CELLS <- 10L
MIN_PAIRED_DONORS <- 3L
FDR_CUT <- 0.10
ABS_LOGFC_CUT <- 0.5
FILTER_MIN_COUNT <- 5
FILTER_MIN_TOTAL_COUNT <- 10
NORMALIZATION <- "TMM"

required <- c("ArchR","Matrix","SummarizedExperiment","edgeR","statmod",
              "dplyr","tidyr","tibble","ggplot2")
absent <- required[!vapply(required,requireNamespace,logical(1),quietly=TRUE)]
if(length(absent)) stop("Missing packages: ",paste(absent,collapse=", "))
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(tibble);library(ggplot2)})
dir.create(OUTDIR,recursive=TRUE,showWarnings=FALSE)
write_table <- function(x,name) write.csv(x,file.path(OUTDIR,name),row.names=FALSE,na="NA")
days <- paste0("day",c(0,1,3,6,14,30,31,33,36,44))
types <- c("CD8_NaiveT","CD4_NaiveT","CD4_Tcm","Treg","Th2","CD4_CD8_DP",
           "Th1","gdT","CD4_Teff","CD8_Tem","Naive_B","Activated_B","Resting_B","PlasmaB")
pr <- if(exists("TB_proj_analysis",envir=.GlobalEnv)) TB_proj_analysis else
  if(exists("TB_proj",envir=.GlobalEnv)) TB_proj else NULL
if(is.null(pr)) stop("Load TB_proj and run 00_run_TB_proj.R first.")
md <- as.data.frame(ArchR::getCellColData(pr))
if(!"TB_Subtype" %in% names(md) && "TB_Subtype_All" %in% names(md)) md$TB_Subtype <- md$TB_Subtype_All
if(!all(c("TB_Subtype","Sample","Donor","Timepoint") %in% names(md)))
  stop("Full project must contain TB_Subtype, Sample, Donor, Timepoint. Verify these labels explicitly.")
cellmeta <- tibble(Cell=rownames(md),Sample=as.character(md$Sample),Donor=as.character(md$Donor),
  Timepoint=as.character(md$Timepoint),TB_Subtype=as.character(md$TB_Subtype)) |>
  mutate(TB_Subtype=recode(TB_Subtype,Activaed_B="Activated_B"))
if(anyNA(cellmeta) || anyDuplicated(cellmeta$Cell)) stop("Invalid cell metadata.")
if(!all(days %in% cellmeta$Timepoint)) stop("Project does not include all ten timepoints.")
if(any(!cellmeta$Timepoint %in% days)) stop("Timepoint labels must be day0,...,day44.")
samples <- cellmeta |> distinct(Sample,Donor,Timepoint)
if(anyDuplicated(samples[c("Donor","Timepoint")]))
  stop("Multiple samples per donor/timepoint; merge technical replicates explicitly first.")

raw <- NULL; source_label <- ""
if (nzchar(RAW_COUNTS_RDS)) {
  raw <- readRDS(RAW_COUNTS_RDS)
  source_label <- RAW_COUNTS_RDS
} else if (exists("seRNA", envir = .GlobalEnv)) {
  raw <- get("seRNA", envir = .GlobalEnv)
  source_label <- "Original session seRNA"
} else {
  stop("Provide RAW_COUNTS_RDS or original raw seRNA. ArchR normalized GeneExpressionMatrix is not used.")
}
annotation <- NULL
if(inherits(raw,"SummarizedExperiment")) {
  an <- SummarizedExperiment::assayNames(raw)
  selected <- if("counts" %in% an) "counts" else if(length(an)==1) an[1] else NA_character_
  if(is.na(selected)) stop("Select the raw-count assay explicitly.")
  x <- SummarizedExperiment::assay(raw,selected)
  rd <- as.data.frame(SummarizedExperiment::rowData(raw))
  ids <- rownames(x)
  if(is.null(ids) || anyDuplicated(ids) || all(grepl("^[0-9]+$",ids))) {
    cand <- intersect(c("gene_id","id","name"),names(rd))
    valid <- cand[vapply(cand,function(z) !anyNA(rd[[z]]) && !anyDuplicated(rd[[z]]),logical(1))]
    if(!length(valid)) stop("No unique gene identifiers; supply a raw count matrix with valid rownames.")
    ids <- as.character(rd[[valid[1]]])
  }
  rownames(x) <- ids
  annotation <- data.frame(GeneID=ids,rd,check.names=FALSE)
} else if(is.list(raw) && "counts" %in% names(raw)) {
  x <- raw$counts
} else x <- raw
if(!is.matrix(x) && !inherits(x,"Matrix")) stop("Raw input must be a gene-by-cell matrix.")
if(is.null(rownames(x)) || is.null(colnames(x)) || anyDuplicated(rownames(x)) || anyDuplicated(colnames(x)))
  stop("Unique gene IDs and full cell barcodes required; duplicated symbols must not be silently merged.")
v <- if(inherits(x,"sparseMatrix")) x@x else as.numeric(x)
if(any(!is.finite(v)) || any(v<0) || any(abs(v-round(v))>1e-8))
  stop("Input is not untransformed integer counts. Return to original UMI counts; do not round normalized values.")
rm(v,raw)
x <- methods::as(Matrix::Matrix(x,sparse=TRUE),"dgCMatrix")
# Missing RNA cells are not treated as zero-expression cells.
overlap <- intersect(cellmeta$Cell,colnames(x))
if(!length(overlap)) stop("No matching cell barcodes. Resolve sample/barcode mapping explicitly.")
cm <- cellmeta |> filter(TB_Subtype %in% types) |>
  mutate(HasRNA=Cell %in% colnames(x))
rna_coverage <- cm |> group_by(Sample,Donor,Timepoint,TB_Subtype) |>
  summarise(nAnnotatedCells=n(),nCells=sum(HasRNA),.groups="drop")
coverage <- tidyr::crossing(samples,TB_Subtype=types) |>
  left_join(rna_coverage,by=c("Sample","Donor","Timepoint","TB_Subtype")) |>
  mutate(nAnnotatedCells=replace_na(nAnnotatedCells,0L),nCells=replace_na(nCells,0L),
    Eligible=nCells>=MIN_CELLS,
    Status=case_when(nAnnotatedCells==0~"No subtype cells",nCells==0~"No matching RNA measurements",
                    nCells<MIN_CELLS~"Below RNA cell threshold",TRUE~"Cell-count eligible"))
write_table(coverage,"S01_full_ten_timepoint_RNA_coverage.csv")
cm <- cm |> filter(HasRNA)
pb <- cm |> count(Sample,Donor,Timepoint,TB_Subtype,name="nCells") |>
  arrange(TB_Subtype,Donor,Timepoint) |> mutate(PBGroup=sprintf("RNAPB%05d",row_number()))
cm <- cm |> left_join(pb,by=c("Sample","Donor","Timepoint","TB_Subtype"))
x <- x[,cm$Cell,drop=FALSE]
membership <- Matrix::sparseMatrix(i=seq_len(nrow(cm)),j=match(cm$PBGroup,pb$PBGroup),
  x=1,dims=c(nrow(cm),nrow(pb)),dimnames=list(cm$Cell,pb$PBGroup))
pb_counts <- x %*% membership
pb$LibrarySize <- as.numeric(Matrix::colSums(pb_counts))
pb$Eligible <- pb$nCells>=MIN_CELLS & pb$LibrarySize>0
write_table(pb,"S02_RNA_pseudobulk_metadata.csv")
if(!is.null(annotation)) write_table(annotation,"S03_gene_annotation.csv")
keep <- pb$Eligible
pb <- pb[keep,];pb_counts <- pb_counts[,pb$PBGroup,drop=FALSE]
saveRDS(list(counts=pb_counts,metadata=pb,source=source_label,threshold=MIN_CELLS),
        file.path(OUTDIR,"RNA_raw_pseudobulk_ten_timepoints.rds"))
rm(x,membership);invisible(gc())

# Main model: stage-wise donor fixed effects; report complete pairs separately.
stage_times <- list(Prime=days[1:6],Booster=days[6:10])
summary_rows <- list(); sig_rows <- list(); n <- 0L
for(s in types) for(st in names(stage_times)) {
  tt <- stage_times[[st]]; baseline <- tt[1]
  m <- pb |> filter(TB_Subtype==s,Timepoint %in% tt) |>
    group_by(Donor) |> filter(n_distinct(Timepoint)>=2) |> ungroup() |>
    arrange(Donor,match(Timepoint,tt))
  records <- lapply(tt[-1],function(target) {
    b <- m$Donor[m$Timepoint==baseline];t <- m$Donor[m$Timepoint==target]
    pairs <- intersect(b,t)
    tibble(TB_Subtype=s,Stage=st,Baseline=baseline,Target=target,
      nBaselineDonors=length(unique(b)),nTargetDonors=length(unique(t)),
      nPairedDonors=length(pairs),PairedDonors=paste(sort(pairs),collapse=";"),
      nSamplesInStageModel=nrow(m),nDonorsInStageModel=n_distinct(m$Donor),
      Status=if(length(pairs)>=MIN_PAIRED_DONORS) "Pending" else "Insufficient paired donors",
      nTestedGenes=NA_integer_,nUp=NA_integer_,nDown=NA_integer_,nDEG=NA_integer_)
  })
  rr <- bind_rows(records)
  fitobj <- NULL
  if(any(rr$Status=="Pending")) fitobj <- tryCatch({
    m$Donor <- factor(m$Donor);m$Timepoint <- factor(m$Timepoint,levels=tt[tt %in% m$Timepoint])
    design <- model.matrix(~0+Timepoint+Donor,m)
    if(qr(design)$rank<ncol(design) || nrow(design)<=ncol(design)) stop("Invalid design / no residual degrees of freedom")
    y <- edgeR::DGEList(counts=pb_counts[,m$PBGroup,drop=FALSE])
    kg <- edgeR::filterByExpr(y,design=design,min.count=FILTER_MIN_COUNT,min.total.count=FILTER_MIN_TOTAL_COUNT)
    if(sum(kg)<100) stop("Fewer than 100 genes passed expression filtering")
    y <- y[kg,,keep.lib.sizes=FALSE];y <- edgeR::calcNormFactors(y,method=NORMALIZATION)
    y <- edgeR::estimateDisp(y,design,robust=TRUE)
    list(fit=edgeR::glmQLFit(y,design,robust=TRUE),design=design)
  },error=function(e)e)
  if (inherits(fitobj, "error")) {
    rr$Status[rr$Status == "Pending"] <- conditionMessage(fitobj)
  } else if (!is.null(fitobj)) {
    for (j in which(rr$Status == "Pending")) {
    cv <- setNames(rep(0,ncol(fitobj$design)),colnames(fitobj$design))
    cv[paste0("Timepoint",rr$Target[j])] <- 1;cv[paste0("Timepoint",baseline)] <- -1
    r <- edgeR::topTags(edgeR::glmQLFTest(fitobj$fit,contrast=cv),n=Inf,sort.by="none")$table |>
      rownames_to_column("GeneID") |> mutate(TB_Subtype=s,Stage=st,Baseline=baseline,Target=rr$Target[j],
        Significant=FDR<=FDR_CUT & abs(logFC)>=ABS_LOGFC_CUT)
    filename <- paste(s,st,rr$Target[j],sep="_")
    write.csv(r,gzfile(file.path(OUTDIR,paste0(filename,"_all_genes.csv.gz"))),row.names=FALSE)
    sig_rows[[filename]] <- r |> filter(Significant)
    rr$nTestedGenes[j] <- nrow(r);rr$nUp[j] <- sum(r$Significant & r$logFC>0)
    rr$nDown[j] <- sum(r$Significant & r$logFC<0);rr$nDEG[j] <- sum(r$Significant);rr$Status[j] <- "OK"
    sens <- tidyr::crossing(FDRThreshold=c(.05,.10),AbsLogFCThreshold=c(.10,.25,.50)) |>
      rowwise() |> mutate(nDEG=sum(r$FDR<=FDRThreshold & abs(r$logFC)>=AbsLogFCThreshold)) |> ungroup()
    write_table(sens,paste0(filename,"_threshold_sensitivity.csv"))
    }
  }
  n <- n+1L;summary_rows[[n]] <- rr
}
su <- bind_rows(summary_rows)
write_table(su,"S04_DEG_counts_and_effective_donors.csv")
write.csv(bind_rows(sig_rows),gzfile(file.path(OUTDIR,"S05_all_significant_DEGs.csv.gz")),row.names=FALSE)
pd <- su |> select(TB_Subtype,Target,Status,nUp,nDown) |>
  pivot_longer(c(nUp,nDown),names_to="Direction",values_to="nDEG") |>
  mutate(TB_Subtype=factor(TB_Subtype,rev(types)),Target=factor(Target,days[-1]))
p <- ggplot(pd,aes(Target,TB_Subtype))+geom_tile(fill="grey97",color="white")+
  geom_point(data=filter(pd,Status=="OK"),aes(size=nDEG,color=Direction))+
  geom_text(data=filter(pd,Status!="OK"),label="NT",size=2.5)+facet_wrap(~Direction,nrow=1)+
  scale_size_area(max_size=8,name="DEG count")+scale_color_manual(values=c(nUp="#C7473D",nDown="#3F6FA3"),guide="none")+
  scale_x_discrete(drop=FALSE)+labs(x=NULL,y=NULL,title="Vaccination-associated differential expression",
    subtitle="Count pseudobulk edgeR QL; FDR <= 0.10 and |log2FC| >= 0.10; NT: not tested")+
  theme_minimal()+theme(panel.grid=element_blank(),axis.text.x=element_text(angle=45,hjust=1))
ggsave(file.path(OUTDIR,"Figure2A_count_based_DEGs.pdf"),p,width=14,height=8)
ggsave(file.path(OUTDIR,"Figure2A_count_based_DEGs.png"),p,width=14,height=8,dpi=300)
writeLines(c(
 "Summed raw measured UMI counts are modelled; no averages of normalized expression are used.",
 paste("Raw source:",source_label),
 "Integer validation is necessary, not sufficient proof of raw-data provenance; verify source counts.",
 "Main DEG threshold: BH FDR <= 0.10 AND |estimated log2FC| >= 0.10.",
 "The QL test tests zero change; logFC threshold is an additional estimate filter, not a threshold-null test.",
 "BH correction is within each subtype/contrast; there is no whole-study FDR guarantee.",
 "FDR 0.10 is an explicitly less stringent choice; FDR 0.05 sensitivity results are provided.",
 "A |log2FC| threshold of 0.10 is small; it does not by itself establish biological importance.",
 "Stage models may use extra eligible observations beyond each complete pair; pair and model counts are reported.",
 "Separate DEG counts cannot demonstrate a broader primary response. Run matched comparisons next.",
 "D30 versus D0 describes the pre-booster state difference, not a comparable post-booster interval.",
 "Coverage reflects RNA barcodes represented in the raw matrix, not just ATAC cell counts.",
 "Missing donor/timepoint samples cannot be reconstructed from cell metadata alone; verify sample manifest."
),file.path(OUTDIR,"analysis_notes.txt"))
writeLines(capture.output(sessionInfo()),file.path(OUTDIR,"sessionInfo.txt"))
message("RNA count reanalysis outputs: ",OUTDIR)


######################################################################revision
# ============================================================
# 1. 加载包与设置目录
# ============================================================
required_packages <- c("dplyr", "tidyr", "ggplot2")

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages)) {
  stop(
    "请先安装缺少的 R 包：",
    paste(missing_packages, collapse = ", ")
  )
}

library(dplyr)
library(tidyr)
library(ggplot2)

INPUT_DIR <- paste0(
  "/public/users/wangrong/home/wangrong/results/T_B/",
  "TB_final_downstream/Figure2_DEG/anno_raw_pseudobulk"
)

FDR_CUT <- 0.10
ABS_LOGFC_CUT <- 0.50

NEW_OUTDIR <- file.path(
  INPUT_DIR,
  sprintf(
    "Rescreen_FDR%.2f_absLog2FC%.2f",
    FDR_CUT,
    ABS_LOGFC_CUT
  )
)

if (!dir.exists(INPUT_DIR)) {
  stop("输入目录不存在：", INPUT_DIR)
}

dir.create(
  NEW_OUTDIR,
  recursive = TRUE,
  showWarnings = FALSE
)

threshold_text <- sprintf(
  "FDR <= %.2f and |log2FC| >= %.2f",
  FDR_CUT,
  ABS_LOGFC_CUT
)

# ============================================================
# 2. 文件读写函数
# ============================================================
read_gz_csv <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con))

  read.csv(
    con,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

write_gz_csv <- function(x, path) {
  con <- gzfile(path, open = "wt")
  on.exit(close(con))

  write.csv(
    x,
    con,
    row.names = FALSE,
    na = "NA"
  )
}

write_new_table <- function(x, filename) {
  write.csv(
    x,
    file.path(NEW_OUTDIR, filename),
    row.names = FALSE,
    na = "NA"
  )
}

# ============================================================
# 3. 读取原始汇总表
# ============================================================
summary_file <- file.path(
  INPUT_DIR,
  "S03_main_DEG_summary.csv"
)

if (!file.exists(summary_file)) {
  stop("未找到原汇总表：", summary_file)
}

old_summary <- read.csv(
  summary_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_summary_cols <- c(
  "TB_Subtype", "Stage", "Baseline", "Target",
  "nPairedDonors", "Donors", "nTestedGenes",
  "nUp", "nDown", "Status"
)

missing_cols <- setdiff(
  required_summary_cols,
  names(old_summary)
)

if (length(missing_cols)) {
  stop(
    "汇总表缺少以下列：",
    paste(missing_cols, collapse = ", ")
  )
}

key_cols <- c("TB_Subtype", "Baseline", "Target")

if (anyNA(old_summary[, c(key_cols, "Status")])) {
  stop("汇总表的比较信息或 Status 存在缺失值。")
}

if (anyDuplicated(old_summary[, key_cols])) {
  stop("汇总表存在重复比较。")
}

# 原脚本的文件名仅由亚群和目标时间组成
file_keys <- paste(
  old_summary$TB_Subtype,
  old_summary$Target,
  sep = "_"
)

if (anyDuplicated(file_keys)) {
  stop("不同记录对应同一个结果文件名，请检查汇总表。")
}

# 保留原分析亚群名称及顺序
types_plot <- unique(old_summary$TB_Subtype)

days <- paste0(
  "day",
  c(0, 1, 3, 6, 14, 30, 31, 33, 36, 44)
)

if (!all(old_summary$Target %in% days[-1])) {
  stop("汇总表包含非预期的目标时间点。")
}

# 核对原分析比较基线
expected_baseline <- ifelse(
  old_summary$Target %in% c(
    "day1", "day3", "day6", "day14", "day30"
  ),
  "day0",
  "day30"
)

if (any(old_summary$Baseline != expected_baseline)) {
  stop("汇总表中的比较基线与预期不一致，请先核对。")
}

# ============================================================
# 4. 初始化新汇总表
# ============================================================
summary_new <- old_summary

summary_new$nUp_previous <- old_summary$nUp
summary_new$nDown_previous <- old_summary$nDown
summary_new$nDEG_previous <- (
  old_summary$nUp + old_summary$nDown
)

summary_new$nUp <- NA_integer_
summary_new$nDown <- NA_integer_
summary_new$nDEG <- NA_integer_
summary_new$nValidResults <- NA_integer_
summary_new$nUnavailableResults <- NA_integer_

summary_new$FDRThreshold <- FDR_CUT
summary_new$AbsLog2FCThreshold <- ABS_LOGFC_CUT

sig_rows <- list()
sensitivity_rows <- list()
source_rows <- list()

message(
  "原汇总表中成功完成检验的比较数：",
  sum(old_summary$Status == "OK")
)

# ============================================================
# 5. 读取全基因结果并重新筛选
# ============================================================
for (i in seq_len(nrow(old_summary))) {

  rec <- old_summary[i, , drop = FALSE]

  # 原来未检验的记录保留，不改为零 DEG
  if (rec$Status != "OK") next

  key <- paste(
    rec$TB_Subtype,
    rec$Target,
    sep = "_"
  )

  input_file <- file.path(
    INPUT_DIR,
    paste0(key, "_all_genes.csv.gz")
  )

  if (!file.exists(input_file)) {
    stop(
      "汇总表显示成功，但缺少全基因结果：",
      input_file
    )
  }

  r <- read_gz_csv(input_file)

  required_result_cols <- c(
    "GeneID", "logFC", "FDR",
    "TB_Subtype", "Baseline", "Target"
  )

  missing_cols <- setdiff(
    required_result_cols,
    names(r)
  )

  if (length(missing_cols)) {
    stop(
      basename(input_file),
      " 缺少列：",
      paste(missing_cols, collapse = ", ")
    )
  }

  if (!nrow(r)) {
    stop("全基因结果为空：", basename(input_file))
  }

  if (anyNA(r$GeneID) ||
      anyDuplicated(r$GeneID)) {
    stop(
      "GeneID 缺失或重复：",
      basename(input_file)
    )
  }

  if (!is.numeric(r$FDR) ||
      !is.numeric(r$logFC)) {
    stop(
      "FDR 或 logFC 不是数值型：",
      basename(input_file)
    )
  }

  # 核对文件内比较信息
  for (cc in key_cols) {
    values <- unique(as.character(r[[cc]]))

    if (length(values) != 1L ||
        anyNA(values) ||
        values != as.character(rec[[cc]])) {
      stop(
        basename(input_file),
        " 的 ", cc, " 与原汇总表不一致。"
      )
    }
  }

  if (!is.na(rec$nTestedGenes) &&
      nrow(r) != rec$nTestedGenes) {
    stop(
      "全基因行数与原 nTestedGenes 不一致：",
      basename(input_file)
    )
  }

  valid <- is.finite(r$FDR) & is.finite(r$logFC)

  finite_fdr <- r$FDR[is.finite(r$FDR)]
  if (any(finite_fdr < 0 | finite_fdr > 1)) {
    stop(
      "FDR 超出 0–1：",
      basename(input_file)
    )
  }

  # 保存旧显著性标记
  if ("Significant" %in% names(r)) {
    r$Significant_previous <- r$Significant
  }

  # 使用原全基因 FDR，不在筛选后重新做 BH
  r$Significant <- (
    valid &
      r$FDR <= FDR_CUT &
      abs(r$logFC) >= ABS_LOGFC_CUT
  )

  r$Direction <- "Not significant"
  r$Direction[!valid] <- "Unavailable"

  r$Direction[
    r$Significant & r$logFC > 0
  ] <- "Higher"

  r$Direction[
    r$Significant & r$logFC < 0
  ] <- "Lower"

  r$Stage <- rec$Stage

  r_sig <- r[
    r$Significant,
    , drop = FALSE
  ]

  r_sig <- r_sig[
    order(r_sig$FDR, -abs(r_sig$logFC)),
    , drop = FALSE
  ]

  # 保存更新后的全基因标记
  write_gz_csv(
    r,
    file.path(
      NEW_OUTDIR,
      paste0(key, "_all_genes.csv.gz")
    )
  )

  # 保存显著及上下调列表，零 DEG 时仍输出表头
  write_new_table(
    r_sig,
    paste0(key, "_significant.csv")
  )

  write_new_table(
    r_sig[
      r_sig$Direction == "Higher",
      , drop = FALSE
    ],
    paste0(key, "_Higher.csv")
  )

  write_new_table(
    r_sig[
      r_sig$Direction == "Lower",
      , drop = FALSE
    ],
    paste0(key, "_Lower.csv")
  )

  sig_rows[[key]] <- r_sig

  summary_new$nTestedGenes[i] <- nrow(r)
  summary_new$nValidResults[i] <- sum(valid)
  summary_new$nUnavailableResults[i] <- sum(!valid)
  summary_new$nUp[i] <- sum(r$Direction == "Higher")
  summary_new$nDown[i] <- sum(r$Direction == "Lower")
  summary_new$nDEG[i] <- sum(r$Significant)

  # 阈值敏感性汇总
  sens <- tidyr::crossing(
    FDRThreshold = sort(unique(
      c(0.05, 0.10, FDR_CUT)
    )),
    AbsLogFCThreshold = sort(unique(
      c(0.10, 0.25, 0.50, ABS_LOGFC_CUT)
    ))
  )

  sens$nDEG <- vapply(
    seq_len(nrow(sens)),
    function(k) {
      as.integer(sum(
        valid &
          r$FDR <= sens$FDRThreshold[k] &
          abs(r$logFC) >= sens$AbsLogFCThreshold[k],
        na.rm = TRUE
      ))
    },
    integer(1)
  )

  sens$TB_Subtype <- rec$TB_Subtype
  sens$Stage <- rec$Stage
  sens$Baseline <- rec$Baseline
  sens$Target <- rec$Target
  sens$nTestedGenes <- nrow(r)

  sensitivity_rows[[key]] <- sens

  source_rows[[key]] <- data.frame(
    TB_Subtype = rec$TB_Subtype,
    Baseline = rec$Baseline,
    Target = rec$Target,
    SourceFile = normalizePath(input_file),
    stringsAsFactors = FALSE
  )

  message(
    key,
    ": Higher = ", summary_new$nUp[i],
    "; Lower = ", summary_new$nDown[i],
    "; Total = ", summary_new$nDEG[i]
  )
}

if (!length(source_rows)) {
  stop("没有可重新筛选的成功比较。")
}

# ============================================================
# 6. 汇总与保存结果
# ============================================================
all_significant <- bind_rows(sig_rows)
sensitivity <- bind_rows(sensitivity_rows)
source_manifest <- bind_rows(source_rows)

summary_new$Delta_nDEG <- (
  summary_new$nDEG -
    summary_new$nDEG_previous
)

write_new_table(
  summary_new,
  "S03_main_DEG_summary.csv"
)

write_new_table(
  summary_new %>%
    filter(
      Baseline == "day0",
      Target == "day30"
    ),
  "D30_vs_D0_summary.csv"
)

write_new_table(
  sensitivity,
  "S04_FDR_logFC_sensitivity.csv"
)

write_gz_csv(
  all_significant,
  file.path(
    NEW_OUTDIR,
    "S06_all_significant_DEGs.csv.gz"
  )
)

write_gz_csv(
  all_significant %>%
    filter(
      Baseline == "day0",
      Target == "day30"
    ),
  file.path(
    NEW_OUTDIR,
    "D30_vs_D0_significant_DEGs.csv.gz"
  )
)

write_new_table(
  source_manifest,
  "Source_result_files.csv"
)

# 保留原模型 QC 的副本，注明来源
qc_file <- file.path(
  INPUT_DIR,
  "S05_model_libraries.csv"
)

if (file.exists(qc_file)) {
  file.copy(
    qc_file,
    file.path(
      NEW_OUTDIR,
      "S05_model_libraries_original.csv"
    ),
    overwrite = TRUE
  )
}

print(
  as_tibble(summary_new) %>%
    select(
      TB_Subtype, Baseline, Target,
      nPairedDonors,
      nDEG_previous,
      nUp, nDown, nDEG,
      Delta_nDEG, Status
    ),
  n = Inf
)

# ============================================================
# 7. 绘图数据
# ============================================================
pd <- summary_new %>%
  pivot_longer(
    cols = c(nUp, nDown),
    names_to = "Direction",
    values_to = "Count"
  ) %>%
  mutate(
    Cell = factor(
      TB_Subtype,
      levels = rev(types_plot)
    ),
    Time = factor(
      Target,
      levels = days[-1]
    ),
    Direction = factor(
      Direction,
      levels = c("nUp", "nDown"),
      labels = c("Higher", "Lower")
    ),
    Label = case_when(
      Status != "OK" ~ "NT",
      Count > 0 ~ as.character(Count),
      TRUE ~ ""
    )
  )

# ============================================================
# 8. 上下调 DEG 数量气泡图
# ============================================================
p <- ggplot(
  pd,
  aes(Time, Cell)
) +
  geom_tile(
    fill = "#F5F5F5",
    colour = "white"
  ) +
  geom_point(
    data = pd %>%
      filter(
        Status == "OK",
        Count > 0
      ),
    aes(
      size = Count,
      colour = Direction
    ),
    alpha = 0.5
  ) +
  geom_text(
    aes(label = Label),
    size = 2.4
  ) +
  facet_wrap(
    ~ Direction,
    nrow = 1
  ) +
  scale_size_area(
    max_size = 12,
    name = "DEG count"
  ) +
  scale_colour_manual(
    values = c(
      Higher = "#D56857",
      Lower = "#5279A3"
    ),
    guide = "none"
  ) +
  scale_x_discrete(
    labels = function(x) sub("^day", "D", x),
    drop = FALSE
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      colour = "black"
    ),
    axis.text.y = element_text(
      colour = "black"
    ),
    plot.caption = element_text(
      hjust = 0
    )
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "Raw-count pseudobulk differential expression",
    subtitle = paste0(
      "Paired donors; ",
      threshold_text
    ),
    caption = paste0(
      "NT: not tested/unavailable; blank tested cells: 0 DEGs.\n",
      "D1, D3, D6, D14 and D30 vs D0; ",
      "D31, D33, D36 and D44 vs D30."
    )
  )

print(p)

plot_height <- max(
  6,
  length(types_plot) * 0.32 + 2
)

ggsave(
  filename = file.path(
    NEW_OUTDIR,
    "Figure2A_raw_count_DEGs.pdf"
  ),
  plot = p,
  width = 11,
  height = plot_height,
  bg = "white"
)

ggsave(
  filename = file.path(
    NEW_OUTDIR,
    "Figure2A_raw_count_DEGs.png"
  ),
  plot = p,
  width = 11,
  height = plot_height,
  dpi = 600,
  bg = "white"
)

# ============================================================
# 9. 单独展示 D30 vs D0 的 DEG 数量
# ============================================================
pd_d30 <- pd %>%
  filter(
    Baseline == "day0",
    Target == "day30"
  ) %>%
  mutate(
    Time = factor(Target, levels = "day30")
  )

if (nrow(pd_d30)) {

  p_d30 <- ggplot(
    pd_d30,
    aes(Time, Cell)
  ) +
    geom_tile(
      fill = "#F5F5F5",
      colour = "white"
    ) +
    geom_point(
      data = pd_d30 %>%
        filter(Status == "OK", Count > 0),
      aes(
        size = Count,
        colour = Direction
      ),
      alpha = 0.5
    ) +
    geom_text(
      aes(label = Label),
      size = 3
    ) +
    facet_wrap(~ Direction, nrow = 1) +
    scale_size_area(
      max_size = 12,
      name = "DEG count"
    ) +
    scale_colour_manual(
      values = c(
        Higher = "#D56857",
        Lower = "#5279A3"
      ),
      guide = "none"
    ) +
    scale_x_discrete(labels = c(day30 = "D30 vs D0")) +
    labs(
      x = NULL,
      y = NULL,
      title = "Pre-booster versus pre-vaccination state",
      subtitle = threshold_text,
      caption = "NT: not tested/unavailable; blank tested cells: 0 DEGs."
    ) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(colour = "black"),
      plot.caption = element_text(hjust = 0)
    )

  print(p_d30)

  ggsave(
    file.path(
      NEW_OUTDIR,
      "D30_vs_D0_DEG_counts.pdf"
    ),
    p_d30,
    width = 7,
    height = plot_height,
    bg = "white"
  )

  ggsave(
    file.path(
      NEW_OUTDIR,
      "D30_vs_D0_DEG_counts.png"
    ),
    p_d30,
    width = 7,
    height = plot_height,
    dpi = 300,
    bg = "white"
  )

} else {
  p_d30 <- NULL
}

# ============================================================
# 10. 保存完整重新筛选对象
# ============================================================
saveRDS(
  list(
    source_directory = INPUT_DIR,
    source_files = source_manifest,
    thresholds = list(
      FDR = FDR_CUT,
      abs_log2FC = ABS_LOGFC_CUT,
      rule = "FDR <= cutoff AND abs(logFC) >= cutoff"
    ),
    summary = summary_new,
    significant_DEGs = all_significant,
    sensitivity = sensitivity,
    plot_data = pd,
    plot = p,
    D30_vs_D0_plot = p_d30
  ),
  file.path(
    NEW_OUTDIR,
    "DEG_rescreened_complete_results.rds"
  )
)

writeLines(
  c(
    paste0("Main DEG threshold: ", threshold_text, "."),
    "Saved all-gene edgeR results were re-filtered; models were not refitted.",
    "Original complete donor pairs, expression filtering and TMM normalization were retained.",
    "Original log2FC, P values and BH FDR values were unchanged.",
    "BH correction remains within each subtype/contrast.",
    "The QL test tests zero change; the log2FC cutoff is an additional estimate filter.",
    "Previously untested comparisons remain unavailable, not zero-DEG comparisons.",
    "D30 versus D0 characterizes the pre-booster state; absence of DEGs does not establish equivalence.",
    paste0("Input directory: ", INPUT_DIR)
  ),
  file.path(
    NEW_OUTDIR,
    "analysis_notes.txt"
  )
)

writeLines(
  capture.output(sessionInfo()),
  file.path(
    NEW_OUTDIR,
    "sessionInfo.txt"
  )
)

message(
  "重新筛选完成。输出目录：\n",
  NEW_OUTDIR
)


# 从上到下的顺序
subtype_order <- c(
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
  "gdT",
  "Treg",
  "CD4_CD8_DP"
)

# 检查是否遗漏亚群
observed <- unique(as.character(pd$Cell))
unknown <- setdiff(observed, subtype_order)

if (length(unknown)) {
  stop("以下亚群未设置顺序：", paste(unknown, collapse = ", "))
}

subtype_order <- subtype_order[subtype_order %in% observed]

# ggplot 的 Y 轴从下往上排列，因此使用 rev()
pd$Cell <- factor(
  as.character(pd$Cell),
  levels = rev(subtype_order)
)

# 显式设置坐标轴顺序，同时作用于所有图层
p_deg <- p_deg +
  scale_y_discrete(
    limits = rev(subtype_order),
    drop = FALSE
  )

print(p_deg)

# 保存更新后的图和数据
plot_height <- max(6, length(subtype_order) * 0.35 + 2)

ggsave(
  file.path(RESULT_DIR, "Figure2A_DEGs_updated_labels.pdf"),
  p_deg,
  width = 13,
  height = plot_height,
  bg = "white"
)

ggsave(
  file.path(RESULT_DIR, "Figure2A_DEGs_updated_labels.png"),
  p_deg,
  width = 13,
  height = plot_height,
  dpi = 600,
  bg = "white"
)

write.csv(
  pd,
  file.path(RESULT_DIR, "Figure2A_DEGs_updated_labels_plot_data.csv"),
  row.names = FALSE
)

saveRDS(
  list(
    plot = p_deg,
    plot_data = pd,
    subtype_order = subtype_order,
    label_mapping = label_mapping
  ),
  file.path(RESULT_DIR, "Figure2A_DEGs_updated_labels_plot.rds")
)