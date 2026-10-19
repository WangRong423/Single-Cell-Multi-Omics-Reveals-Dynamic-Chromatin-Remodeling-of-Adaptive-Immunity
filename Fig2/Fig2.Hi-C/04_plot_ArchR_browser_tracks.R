#!/usr/bin/env Rscript

# ============================================================
# ArchR browser tracks for representative T/B loci
# ============================================================
# Debug fixes retained:
#   * Do not pass unsupported maxCells/highlight arguments.
#   * Inspect formals(plotBrowserTrack) and pass only supported args.
#   * Do not use plotPDF(plotList=...) for nested gtable/list output.
#   * Recursively collect drawables and save with grid.draw().
# ============================================================

suppressPackageStartupMessages({
  library(ArchR)
  library(GenomicRanges)
  library(data.table)
  library(grid)
  library(ggplot2)
  library(gtable)
})

BASE <- "/public/users/wangrong/data/public_HiC/analysis"
BDIR <- file.path(BASE, "GSE84022")
DIRECT <- file.path(BDIR, "direct_contact")
OUTDIR <- file.path(BASE, "unified_3D_genome", "ArchR_tracks")
dir.create(OUTDIR, recursive=TRUE, showWarnings=FALSE)

NAIVE_PROJ_DIR <- file.path(BDIR, "ArchR_Naive_B")
ACT_PROJ_DIR   <- file.path(BDIR, "ArchR_Activated_B")
# Set this environment variable if the CD4 Tcm project path differs:
CD4_PROJ_DIR <- Sys.getenv("CD4_TCM_PROJ_DIR", unset="")
REP_FILE <- file.path(DIRECT, "Bcell_representative_locus_per_subtype.tsv")

# fixed longitudinal palette
day_levels <- c("day0","day1","day3","day6","day14","day30","day31","day33","day36","day44")
day_cols <- c(
  day0="#E41A1C", day1="#304B9B", day3="#287A47", day6="#8C4A98", day14="#F39C12",
  day30="#D5BF00", day31="#758FC1", day33="#B05091", day36="#BE9265", day44="#59B5C4"
)
suffix_to_day <- c("1"="day0","2"="day1","3"="day3","4"="day6","5"="day14","6"="day30","7"="day31","8"="day33","9"="day36","10"="day44")

prepare_project <- function(proj) {
  smp <- as.character(proj$Sample)
  suffix <- sub("^.*-", "", smp)
  dd <- unname(suffix_to_day[suffix])
  if (anyNA(dd)) stop("Unknown sample suffixes: ", paste(unique(smp[is.na(dd)]), collapse=", "))
  proj$PlotDay <- factor(dd, levels=day_levels)
  proj
}

safe_browser <- function(proj, region, features, title=NULL) {
  args <- list(
    ArchRProj=proj,
    groupBy="PlotDay",
    useGroups=day_levels,
    region=region,
    features=features,
    plotSummary=c("bulkTrack","featureTrack","geneTrack"),
    useMatrix="TileMatrix",
    log2Norm=FALSE,
    tileSize=250,
    minCells=25,
    normMethod="ReadsInTSS",
    pal=day_cols
  )
  ff <- names(formals(ArchR::plotBrowserTrack))
  args <- args[names(args) %in% ff]
  p <- do.call(ArchR::plotBrowserTrack, args)
  if (!is.null(title)) attr(p, "HiC_title") <- title
  p
}

collect_drawables <- function(x) {
  out <- list()
  recurse <- function(z) {
    if (inherits(z, "ggplot")) {
      out[[length(out)+1L]] <<- ggplotGrob(z)
    } else if (inherits(z, c("grob","gtable","gTree"))) {
      out[[length(out)+1L]] <<- z
    } else if (is.list(z)) {
      for (zz in z) recurse(zz)
    }
  }
  recurse(x)
  out
}

save_archr_pdf <- function(p, file, width=11, height=7) {
  grobs <- collect_drawables(p)
  if (!length(grobs)) stop("No drawable grob/ggplot was returned by plotBrowserTrack().")
  pdf(file, width=width, height=height, onefile=TRUE, useDingbats=FALSE)
  on.exit(dev.off(), add=TRUE)
  for (g in grobs) {
    grid.newpage()
    grid.draw(g)
  }
}

make_region <- function(chr, p1, p2, min_width=170000, min_flank=60000) {
  mid <- as.integer((p1+p2)/2)
  width <- max(min_width, abs(p2-p1) + 2*min_flank)
  st <- max(1L, mid - as.integer(width/2))
  en <- mid + as.integer(width/2)
  GRanges(chr, IRanges(st,en))
}

plot_B <- function(subtype, proj_dir, outname) {
  stopifnot(dir.exists(proj_dir), file.exists(REP_FILE))
  proj <- prepare_project(loadArchRProject(path=proj_dir, showLogo=FALSE))
  dat <- fread(REP_FILE, data.table=FALSE)
  r <- dat[dat$B_subtype == subtype, , drop=FALSE]
  if (!nrow(r)) stop("No representative locus for ", subtype)
  r <- r[1,]

  chr <- as.character(r$DAR_chr)
  dar_mid <- as.integer((r$DAR_start+r$DAR_end)/2)
  target_col <- intersect(c("P2G_target_pos","Target_pos_hg38","target_pos_hg38"), colnames(r))
  if (!length(target_col)) stop("No hg38 P2G target coordinate column in representative table.")
  target <- as.integer(r[[target_col[1]]])

  region <- make_region(chr, dar_mid, target)
  features <- GRangesList(
    DAR=GRanges(chr, IRanges(as.integer(r$DAR_start), as.integer(r$DAR_end))),
    `P2G target`=GRanges(chr, IRanges(target, target))
  )
  p <- safe_browser(proj, region, features, subtype)
  save_archr_pdf(p, file.path(OUTDIR, outname), width=11, height=7)
}

plot_B("Naive_B", NAIVE_PROJ_DIR, "B_NaiveB_ArchR_browser.pdf")
plot_B("Activated_B", ACT_PROJ_DIR, "C_ActivatedB_ArchR_browser.pdf")

# Optional CD4 Tcm locus from the debugged analysis.
if (nzchar(CD4_PROJ_DIR) && dir.exists(CD4_PROJ_DIR)) {
  proj <- prepare_project(loadArchRProject(path=CD4_PROJ_DIR, showLogo=FALSE))
  chr <- "chr10"
  dar <- GRanges(chr, IRanges(69045494,69045994))
  target <- 69088102L
  region <- make_region(chr, 69045744L, target)
  features <- GRangesList(DAR=dar, `P2G target`=GRanges(chr, IRanges(target,target)))
  p <- safe_browser(proj, region, features, "CD4 Tcm")
  save_archr_pdf(p, file.path(OUTDIR, "A_CD4Tcm_ArchR_browser.pdf"), width=11, height=7)
} else {
  message("CD4_TCM_PROJ_DIR not supplied; skipped CD4 Tcm browser track.")
}

message("ArchR browser PDFs written to: ", OUTDIR)
