#!/usr/bin/env Rscript

# ============================================================
# Extract donor x timepoint ATAC trajectories for representative
# B-cell DAR–P2G loci from the separate ArchR subtype projects.
# ============================================================
# Debug fixes retained:
#   1) Use separate naive_b_proj / activated_b_proj projects.
#   2) Build donor-timepoint pseudobulk groups from Sample names.
#   3) Use getPeakSet(proj) for peak coordinates instead of relying
#      on rowRanges(getGroupSE(...)), which was NULL in the local setup.
#   4) Match candidate DAR to the ArchR peak by genomic overlap.
# ============================================================

suppressPackageStartupMessages({
  library(ArchR)
  library(GenomicRanges)
  library(SummarizedExperiment)
  library(data.table)
})

BASE <- "/public/users/wangrong/data/public_HiC/analysis/GSE84022"
DIRECT <- file.path(BASE, "direct_contact")

NAIVE_PROJ_DIR <- file.path(BASE, "ArchR_Naive_B")
ACT_PROJ_DIR   <- file.path(BASE, "ArchR_Activated_B")

CANDIDATE_FILE <- file.path(DIRECT, "Bcell_representative_locus_per_subtype.tsv")
OUT_LONG <- file.path(DIRECT, "Bcell_top_candidates_ATAC_donor_timepoint.tsv")
OUT_SUM  <- file.path(DIRECT, "Bcell_top_candidates_ATAC_timepoint_summary.tsv")

stopifnot(file.exists(CANDIDATE_FILE))
stopifnot(dir.exists(NAIVE_PROJ_DIR), dir.exists(ACT_PROJ_DIR))

# Sample suffix mapping established for the 40 donor x timepoint samples.
suffix_to_day <- c(
  "1"="day0", "2"="day1", "3"="day3", "4"="day6", "5"="day14",
  "6"="day30", "7"="day31", "8"="day33", "9"="day36", "10"="day44"
)
day_levels <- c("day0","day1","day3","day6","day14","day30","day31","day33","day36","day44")
day_num <- c(day0=0, day1=1, day3=3, day6=6, day14=14, day30=30, day31=31, day33=33, day36=36, day44=44)

parse_sample <- function(sample_vec) {
  donor <- sub("-.*$", "", sample_vec)
  suffix <- sub("^.*-", "", sample_vec)
  day <- unname(suffix_to_day[suffix])
  if (anyNA(day)) {
    bad <- unique(sample_vec[is.na(day)])
    stop("Unrecognized Sample names: ", paste(bad, collapse=", "))
  }
  data.frame(Donor=donor, PlotDay=day, stringsAsFactors=FALSE)
}

load_project <- function(path) {
  message("Loading ArchR project: ", path)
  loadArchRProject(path=path, showLogo=FALSE)
}

prepare_project <- function(proj) {
  sample_vec <- as.character(proj$Sample)
  mm <- parse_sample(sample_vec)
  proj$HiC_Donor <- mm$Donor
  proj$HiC_PlotDay <- factor(mm$PlotDay, levels=day_levels)
  proj$HiC_PBGroup <- paste(mm$Donor, mm$PlotDay, sep="__")
  proj
}

find_peak_index <- function(proj, chr, start, end) {
  peak_set <- getPeakSet(proj)
  query <- GRanges(seqnames=chr, ranges=IRanges(start=as.integer(start), end=as.integer(end)))
  hits <- findOverlaps(query, peak_set, ignore.strand=TRUE)
  if (length(hits) == 0L) {
    stop(sprintf("No ArchR peak overlaps %s:%s-%s", chr, start, end))
  }
  idx <- unique(subjectHits(hits))
  if (length(idx) > 1L) {
    # Choose the overlapping peak whose center is closest to DAR center.
    dar_mid <- (as.numeric(start) + as.numeric(end))/2
    mids <- (start(peak_set[idx]) + end(peak_set[idx]))/2
    idx <- idx[which.min(abs(mids - dar_mid))]
  }
  as.integer(idx)
}

extract_candidate <- function(proj, cand_row) {
  peak_idx <- find_peak_index(
    proj,
    as.character(cand_row$DAR_chr),
    cand_row$DAR_start,
    cand_row$DAR_end
  )

  peak_set <- getPeakSet(proj)
  matched_peak <- peak_set[peak_idx]

  gse <- getGroupSE(
    ArchRProj=proj,
    useMatrix="PeakMatrix",
    groupBy="HiC_PBGroup",
    divideN=FALSE,
    scaleTo=NULL
  )

  mat <- assay(gse)
  if (nrow(mat) != length(peak_set)) {
    stop("PeakMatrix row count does not match getPeakSet(proj).")
  }

  counts <- as.numeric(mat[peak_idx, ])
  lib <- colSums(mat)
  cpm <- ifelse(lib > 0, counts/lib*1e6, NA_real_)
  log2cpm <- log2(cpm + 1)

  groups <- colnames(mat)
  # ArchR may sanitize group names; support both '__' and '.' separators.
  donor <- sub("(__|\\.).*$", "", groups)
  day <- sub("^.*(__|\\.)", "", groups)

  # If day extraction failed, recover from known day labels within group string.
  recover_day <- function(x) {
    hit <- day_levels[vapply(day_levels, function(d) grepl(d, x, fixed=TRUE), logical(1))]
    if (length(hit)) hit[1] else NA_character_
  }
  bad <- !(day %in% day_levels)
  if (any(bad)) day[bad] <- vapply(groups[bad], recover_day, character(1))

  data.frame(
    Candidate_ID=as.character(cand_row$Candidate_ID),
    B_subtype=as.character(cand_row$B_subtype),
    DAR_chr=as.character(cand_row$DAR_chr),
    DAR_start=as.integer(cand_row$DAR_start),
    DAR_end=as.integer(cand_row$DAR_end),
    ArchR_peak_chr=as.character(seqnames(matched_peak)),
    ArchR_peak_start=start(matched_peak),
    ArchR_peak_end=end(matched_peak),
    Donor=donor,
    PlotDay=day,
    Day=unname(day_num[day]),
    PeakCounts=counts,
    LibraryPeakCounts=lib,
    CPM=cpm,
    log2CPM=log2cpm,
    stringsAsFactors=FALSE
  )
}

cand <- fread(CANDIDATE_FILE, data.table=FALSE)
needed <- c("Candidate_ID","B_subtype","DAR_chr","DAR_start","DAR_end")
miss <- setdiff(needed, colnames(cand))
if (length(miss)) stop("Missing candidate columns: ", paste(miss, collapse=", "))

naive_b_proj <- prepare_project(load_project(NAIVE_PROJ_DIR))
activated_b_proj <- prepare_project(load_project(ACT_PROJ_DIR))

out <- list()
for (i in seq_len(nrow(cand))) {
  r <- cand[i, , drop=FALSE]
  subtype <- as.character(r$B_subtype)
  proj <- switch(
    subtype,
    Naive_B=naive_b_proj,
    Activated_B=activated_b_proj,
    stop("Unsupported B_subtype: ", subtype)
  )
  out[[length(out)+1L]] <- extract_candidate(proj, r)
}

long <- rbindlist(out, fill=TRUE)
long$PlotDay <- factor(long$PlotDay, levels=day_levels)
setorder(long, B_subtype, Candidate_ID, Donor, Day)
fwrite(long, OUT_LONG, sep="\t")

sumdt <- long[, .(
  n_donors=uniqueN(Donor[is.finite(log2CPM)]),
  mean_log2CPM=mean(log2CPM, na.rm=TRUE),
  sd_log2CPM=sd(log2CPM, na.rm=TRUE),
  sem_log2CPM=sd(log2CPM, na.rm=TRUE)/sqrt(sum(is.finite(log2CPM)))
), by=.(Candidate_ID, B_subtype, PlotDay, Day)]
setorder(sumdt, B_subtype, Candidate_ID, Day)
fwrite(sumdt, OUT_SUM, sep="\t")

message("Saved:\n  ", OUT_LONG, "\n  ", OUT_SUM)
