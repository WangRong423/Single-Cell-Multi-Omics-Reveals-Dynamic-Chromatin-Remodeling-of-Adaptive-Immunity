# ============================================================
# Propeller longitudinal cell-proportion analysis
# Relative log2 fold change versus Day 0
#
# Input object:
#   immune_merged
#
# Required metadata:
#   sample
#   celltype_merged
#
# Main outputs:
#   1. Propeller transformed-proportion statistics
#   2. Donor-paired proportion log2FC relative to Day 0
#   3. Heatmap + lineage annotation + cell-count panel
#   4. p_propeller_final_log2fc saved as an RDS object
# ============================================================


# ------------------------------------------------------------
# 0. Load packages
# ------------------------------------------------------------

required_packages <- c(
  "Seurat",
  "speckle",
  "limma",
  "dplyr",
  "tidyr",
  "tibble",
  "ggplot2",
  "patchwork",
  "scales",
  "grid"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "The following packages are missing: ",
    paste(missing_packages, collapse = ", "),
    "\nInstall CRAN packages with install.packages() and Bioconductor packages ",
    "speckle/limma with BiocManager::install()."
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(speckle)
  library(limma)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(grid)
})


# ------------------------------------------------------------
# 1. General settings
# ------------------------------------------------------------

if (!exists("immune_merged")) {
  stop("Object 'immune_merged' was not found in the R environment.")
}

required_meta_columns <- c(
  "sample",
  "celltype_merged"
)

missing_meta_columns <- setdiff(
  required_meta_columns,
  colnames(immune_merged@meta.data)
)

if (length(missing_meta_columns) > 0) {
  stop(
    "The following metadata columns are missing from immune_merged: ",
    paste(missing_meta_columns, collapse = ", ")
  )
}

outdir <- "/public/users/wangrong/home/wangrong/results/T_B/merged_annotation/Propeller_cell_proportions"

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

time_map <- c(
  "1"  = "day0",
  "2"  = "day1",
  "3"  = "day3",
  "4"  = "day6",
  "5"  = "day14",
  "6"  = "day30",
  "7"  = "day31",
  "8"  = "day33",
  "9"  = "day36",
  "10" = "day44"
)

time_levels <- c(
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


# ------------------------------------------------------------
# 2. Rebuild donor and timepoint from sample
#    Example:
#      M1-1  -> donor M1, visit 1, day0
#      M1-10 -> donor M1, visit 10, day44
# ------------------------------------------------------------

sample_vec <- as.character(immune_merged$sample)

bad_samples <- unique(
  sample_vec[
    is.na(sample_vec) |
      !grepl("^.+-[0-9]+$", sample_vec)
  ]
)

if (length(bad_samples) > 0) {
  stop(
    "The following sample names do not match the expected format donor-visit: ",
    paste(bad_samples, collapse = ", ")
  )
}

immune_merged$donor_propeller <- sub(
  "-[0-9]+$",
  "",
  sample_vec
)

immune_merged$visit_number <- sub(
  "^.*-([0-9]+)$",
  "\\1",
  sample_vec
)

immune_merged$timepoint <- unname(
  time_map[
    as.character(immune_merged$visit_number)
  ]
)

immune_merged$timepoint <- factor(
  immune_merged$timepoint,
  levels = time_levels
)

immune_merged$donor_propeller <- factor(
  immune_merged$donor_propeller
)

if (anyNA(immune_merged$timepoint)) {
  unmapped_visits <- unique(
    as.character(immune_merged$visit_number)[
      is.na(immune_merged$timepoint)
    ]
  )

  stop(
    "The following visit numbers could not be mapped to timepoints: ",
    paste(unmapped_visits, collapse = ", ")
  )
}

message("\nCell counts by timepoint:")
print(
  table(
    immune_merged$timepoint,
    useNA = "ifany"
  )
)

message("\nDonor-by-timepoint cell counts:")
print(
  table(
    immune_merged$donor_propeller,
    immune_merged$timepoint
  )
)


# ------------------------------------------------------------
# 3. Build major lineage annotation
# ------------------------------------------------------------

immune_merged$major_group <- dplyr::case_when(
  grepl(
    "^B_",
    as.character(immune_merged$celltype_merged)
  ) ~ "B cells",

  grepl(
    "^CD4\\.",
    as.character(immune_merged$celltype_merged)
  ) ~ "CD4 T cells",

  grepl(
    "^CD8\\.",
    as.character(immune_merged$celltype_merged)
  ) |
    as.character(immune_merged$celltype_merged) == "MAIT" ~
    "CD8 T cells",

  TRUE ~ NA_character_
)

immune_merged$major_group <- factor(
  immune_merged$major_group,
  levels = c(
    "B cells",
    "CD4 T cells",
    "CD8 T cells"
  )
)

if (anyNA(immune_merged$major_group)) {
  unknown_celltypes <- unique(
    as.character(immune_merged$celltype_merged)[
      is.na(immune_merged$major_group)
    ]
  )

  stop(
    "The following cell types could not be assigned to a major lineage: ",
    paste(unknown_celltypes, collapse = ", ")
  )
}


# ------------------------------------------------------------
# 4. Build Propeller metadata
# ------------------------------------------------------------

meta_prop <- immune_merged@meta.data %>%
  rownames_to_column("cell") %>%
  transmute(
    cell = cell,
    celltype = as.character(celltype_merged),
    major_group = as.character(major_group),
    sample = as.character(sample),
    donor = as.character(donor_propeller),
    time = as.character(timepoint)
  ) %>%
  filter(
    !is.na(celltype),
    !is.na(major_group),
    !is.na(sample),
    !is.na(donor),
    !is.na(time),
    time %in% time_levels
  ) %>%
  mutate(
    donor = factor(donor),
    time = factor(
      time,
      levels = time_levels
    )
  )

message("\nPropeller metadata: cell counts by timepoint")
print(
  table(
    meta_prop$time,
    useNA = "ifany"
  )
)

message("\nPropeller metadata: donor-by-timepoint cell counts")
print(
  table(
    meta_prop$donor,
    meta_prop$time
  )
)

message("\nPropeller metadata: lineage-by-timepoint cell counts")
print(
  table(
    meta_prop$major_group,
    meta_prop$time
  )
)


# Ensure one sample corresponds to one donor and one timepoint
sample_mapping_check <- meta_prop %>%
  distinct(
    sample,
    donor,
    time
  ) %>%
  count(
    sample,
    name = "n_mapping"
  ) %>%
  filter(
    n_mapping != 1
  )

if (nrow(sample_mapping_check) > 0) {
  print(sample_mapping_check)

  stop(
    "At least one sample corresponds to more than one donor/timepoint."
  )
}


# ------------------------------------------------------------
# 5. Run Propeller with arcsine-square-root transformation
# ------------------------------------------------------------

propeller_list <- speckle::getTransformedProps(
  clusters = meta_prop$celltype,
  sample = meta_prop$sample,
  transform = "asin"
)

required_propeller_objects <- c(
  "Counts",
  "TransformedProps",
  "Proportions"
)

missing_propeller_objects <- setdiff(
  required_propeller_objects,
  names(propeller_list)
)

if (length(missing_propeller_objects) > 0) {
  stop(
    "The Propeller output is missing: ",
    paste(missing_propeller_objects, collapse = ", ")
  )
}


# ------------------------------------------------------------
# 6. Convert Propeller outputs to plain numeric matrices
#    This avoids table-class conversion errors involving
#    sample/clusters/Freq columns.
# ------------------------------------------------------------

to_plain_matrix <- function(x, object_name = "object") {

  if (length(dim(x)) != 2) {
    stop(object_name, " is not a two-dimensional object.")
  }

  if (is.null(rownames(x)) || is.null(colnames(x))) {
    stop(object_name, " is missing row names or column names.")
  }

  out <- matrix(
    data = as.numeric(x),
    nrow = nrow(x),
    ncol = ncol(x),
    dimnames = dimnames(x)
  )

  storage.mode(out) <- "double"

  if (anyDuplicated(rownames(out)) > 0) {
    stop(object_name, " contains duplicated cell-type names.")
  }

  if (anyDuplicated(colnames(out)) > 0) {
    stop(object_name, " contains duplicated sample names.")
  }

  out
}

tp_mat <- to_plain_matrix(
  propeller_list$TransformedProps,
  "TransformedProps"
)

counts_mat <- to_plain_matrix(
  propeller_list$Counts,
  "Counts"
)

props_mat <- to_plain_matrix(
  propeller_list$Proportions,
  "Proportions"
)

stopifnot(
  identical(dim(tp_mat), dim(counts_mat)),
  identical(dim(tp_mat), dim(props_mat)),
  identical(rownames(tp_mat), rownames(counts_mat)),
  identical(rownames(tp_mat), rownames(props_mat)),
  identical(colnames(tp_mat), colnames(counts_mat)),
  identical(colnames(tp_mat), colnames(props_mat))
)


# ------------------------------------------------------------
# 7. Build and align sample information
# ------------------------------------------------------------

sample_info <- meta_prop %>%
  distinct(
    sample,
    donor,
    time
  ) %>%
  mutate(
    sample = as.character(sample),
    donor = factor(as.character(donor)),
    time = factor(
      as.character(time),
      levels = time_levels
    )
  )

missing_sample_info <- setdiff(
  colnames(tp_mat),
  sample_info$sample
)

if (length(missing_sample_info) > 0) {
  stop(
    "The following Propeller samples were not found in sample_info: ",
    paste(missing_sample_info, collapse = ", ")
  )
}

sample_info <- sample_info[
  match(
    colnames(tp_mat),
    sample_info$sample
  ),
  ,
  drop = FALSE
]

if (anyNA(sample_info$sample)) {
  stop("Failed to align sample_info with Propeller matrix columns.")
}

stopifnot(
  identical(
    colnames(tp_mat),
    sample_info$sample
  )
)

rownames(sample_info) <- sample_info$sample

message("\nSample-level donor-by-timepoint table:")
print(
  table(
    sample_info$donor,
    sample_info$time
  )
)


# Check that each donor has exactly one Day 0 sample
day0_check <- sample_info %>%
  filter(
    time == "day0"
  ) %>%
  count(
    donor,
    name = "n_day0"
  )

all_donors <- unique(
  as.character(sample_info$donor)
)

missing_day0_donors <- setdiff(
  all_donors,
  as.character(day0_check$donor)
)

if (length(missing_day0_donors) > 0) {
  stop(
    "The following donors do not have a Day 0 sample: ",
    paste(missing_day0_donors, collapse = ", ")
  )
}

if (any(day0_check$n_day0 != 1)) {
  print(day0_check)

  stop(
    "At least one donor has more than one Day 0 sample."
  )
}


# ------------------------------------------------------------
# 8. Long-format transformed proportions
# ------------------------------------------------------------

transformed_long <- as.data.frame(
  as.table(tp_mat),
  stringsAsFactors = FALSE
)

colnames(transformed_long) <- c(
  "celltype",
  "sample",
  "transformed_prop"
)

transformed_long <- transformed_long %>%
  mutate(
    celltype = as.character(celltype),
    sample = as.character(sample),
    transformed_prop = as.numeric(transformed_prop)
  ) %>%
  left_join(
    sample_info %>%
      transmute(
        sample = as.character(sample),
        donor = as.character(donor),
        time = as.character(time)
      ),
    by = "sample"
  ) %>%
  mutate(
    donor = factor(donor),
    time = factor(
      time,
      levels = time_levels
    )
  )

unmatched_transformed_samples <- transformed_long %>%
  filter(
    is.na(donor) |
      is.na(time)
  ) %>%
  distinct(sample)

if (nrow(unmatched_transformed_samples) > 0) {
  print(unmatched_transformed_samples)

  stop(
    "Some transformed-proportion samples were not matched to donor/time."
  )
}

stopifnot(
  nrow(transformed_long) ==
    nrow(tp_mat) * ncol(tp_mat)
)


# ------------------------------------------------------------
# 9. Donor-paired transformed-proportion changes versus Day 0
# ------------------------------------------------------------

baseline_transformed_df <- transformed_long %>%
  filter(
    time == "day0"
  ) %>%
  group_by(
    donor,
    celltype
  ) %>%
  summarise(
    day0_transformed = mean(
      transformed_prop,
      na.rm = TRUE
    ),
    n_day0_samples = n_distinct(sample),
    .groups = "drop"
  )

delta_long <- transformed_long %>%
  left_join(
    baseline_transformed_df %>%
      select(
        donor,
        celltype,
        day0_transformed
      ),
    by = c(
      "donor",
      "celltype"
    )
  ) %>%
  mutate(
    delta_vs_day0 =
      transformed_prop -
      day0_transformed
  )

missing_transformed_baseline <- delta_long %>%
  filter(
    is.na(day0_transformed)
  ) %>%
  distinct(
    donor,
    celltype
  )

if (nrow(missing_transformed_baseline) > 0) {
  print(missing_transformed_baseline)

  stop(
    "Some donor-celltype combinations are missing a Day 0 transformed baseline."
  )
}

message("\nDay 0 transformed-delta check:")
print(
  delta_long %>%
    filter(
      time == "day0"
    ) %>%
    summarise(
      min_delta = min(
        delta_vs_day0,
        na.rm = TRUE
      ),
      max_delta = max(
        delta_vs_day0,
        na.rm = TRUE
      ),
      mean_abs_delta = mean(
        abs(delta_vs_day0),
        na.rm = TRUE
      )
    )
)


# ------------------------------------------------------------
# 10. Paired limma model on Propeller-transformed proportions
# ------------------------------------------------------------

present_times <- time_levels[
  time_levels %in%
    unique(as.character(sample_info$time))
]

sample_info$time <- factor(
  as.character(sample_info$time),
  levels = present_times
)

sample_info$donor <- droplevels(
  factor(sample_info$donor)
)

design <- model.matrix(
  ~ 0 + time + donor,
  data = sample_info
)

colnames(design) <- sub(
  "^time",
  "",
  colnames(design)
)

design_rank <- qr(design)$rank

if (design_rank < ncol(design)) {
  stop(
    "The limma design matrix is not full rank. Rank = ",
    design_rank,
    "; columns = ",
    ncol(design),
    "."
  )
}

comparison_times <- setdiff(
  present_times,
  "day0"
)

contrast_expressions <- paste0(
  comparison_times,
  "-day0"
)

contrast_matrix <- limma::makeContrasts(
  contrasts = contrast_expressions,
  levels = design
)

colnames(contrast_matrix) <- comparison_times

tp_for_model <- tp_mat[
  ,
  sample_info$sample,
  drop = FALSE
]

stopifnot(
  identical(
    colnames(tp_for_model),
    sample_info$sample
  )
)

fit <- limma::lmFit(
  object = tp_for_model,
  design = design
)

fit2 <- limma::contrasts.fit(
  fit = fit,
  contrasts = contrast_matrix
)

fit2 <- limma::eBayes(
  fit2,
  robust = TRUE,
  trend = FALSE
)

stats_list <- lapply(
  seq_along(comparison_times),
  function(i) {

    result_i <- limma::topTable(
      fit2,
      coef = i,
      number = Inf,
      sort.by = "none",
      adjust.method = "none"
    )

    result_i %>%
      rownames_to_column("celltype") %>%
      transmute(
        celltype = as.character(celltype),
        time = comparison_times[i],
        transformed_difference = logFC,
        AveExpr = AveExpr,
        t = t,
        P.Value = P.Value
      )
  }
)

stats_df <- bind_rows(stats_list) %>%
  mutate(
    FDR_global = p.adjust(
      P.Value,
      method = "BH"
    )
  ) %>%
  group_by(time) %>%
  mutate(
    FDR_within_time = p.adjust(
      P.Value,
      method = "BH"
    )
  ) %>%
  ungroup()


# ------------------------------------------------------------
# 11. Long-format raw proportions
# ------------------------------------------------------------

props_long <- as.data.frame(
  as.table(props_mat),
  stringsAsFactors = FALSE
)

colnames(props_long) <- c(
  "celltype",
  "sample",
  "proportion"
)

props_long <- props_long %>%
  mutate(
    celltype = as.character(celltype),
    sample = as.character(sample),
    proportion = as.numeric(proportion)
  ) %>%
  left_join(
    sample_info %>%
      transmute(
        sample = as.character(sample),
        donor = as.character(donor),
        time = as.character(time)
      ),
    by = "sample"
  ) %>%
  mutate(
    donor = factor(donor),
    time = factor(
      time,
      levels = time_levels
    )
  )

unmatched_proportion_samples <- props_long %>%
  filter(
    is.na(donor) |
      is.na(time)
  ) %>%
  distinct(sample)

if (nrow(unmatched_proportion_samples) > 0) {
  print(unmatched_proportion_samples)

  stop(
    "Some raw-proportion samples were not matched to donor/time."
  )
}


# ------------------------------------------------------------
# 12. Donor-paired log2FC of proportions versus Day 0
# ------------------------------------------------------------

# Fixed pseudocount prevents log2(0) for rare populations.
prop_pseudocount <- 1e-6

baseline_prop_df <- props_long %>%
  filter(
    time == "day0"
  ) %>%
  group_by(
    donor,
    celltype
  ) %>%
  summarise(
    day0_proportion = mean(
      proportion,
      na.rm = TRUE
    ),
    n_day0_samples = n_distinct(sample),
    .groups = "drop"
  )

log2fc_long <- props_long %>%
  left_join(
    baseline_prop_df %>%
      select(
        donor,
        celltype,
        day0_proportion
      ),
    by = c(
      "donor",
      "celltype"
    )
  ) %>%
  mutate(
    log2FC_vs_day0 = log2(
      (proportion + prop_pseudocount) /
        (day0_proportion + prop_pseudocount)
    )
  )

missing_prop_baseline <- log2fc_long %>%
  filter(
    is.na(day0_proportion)
  ) %>%
  distinct(
    donor,
    celltype
  )

if (nrow(missing_prop_baseline) > 0) {
  print(missing_prop_baseline)

  stop(
    "Some donor-celltype combinations are missing a Day 0 proportion baseline."
  )
}

message("\nDay 0 log2FC check:")
print(
  log2fc_long %>%
    filter(
      time == "day0"
    ) %>%
    summarise(
      min_log2FC = min(
        log2FC_vs_day0,
        na.rm = TRUE
      ),
      max_log2FC = max(
        log2FC_vs_day0,
        na.rm = TRUE
      ),
      mean_abs_log2FC = mean(
        abs(log2FC_vs_day0),
        na.rm = TRUE
      )
    )
)

heatmap_log2fc_df <- log2fc_long %>%
  group_by(
    celltype,
    time
  ) %>%
  summarise(
    mean_log2FC = mean(
      log2FC_vs_day0,
      na.rm = TRUE
    ),
    sd_log2FC = sd(
      log2FC_vs_day0,
      na.rm = TRUE
    ),
    n_donors = n_distinct(
      donor[
        !is.na(log2FC_vs_day0)
      ]
    ),
    se_log2FC = ifelse(
      n_donors > 1,
      sd_log2FC / sqrt(n_donors),
      NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    time = as.character(time)
  ) %>%
  left_join(
    stats_df %>%
      select(
        celltype,
        time,
        transformed_difference,
        P.Value,
        FDR_within_time,
        FDR_global
      ),
    by = c(
      "celltype",
      "time"
    )
  ) %>%
  mutate(
    significance = case_when(
      time == "day0" ~ "",
      !is.na(FDR_within_time) &
        FDR_within_time < 0.001 ~ "***",
      !is.na(FDR_within_time) &
        FDR_within_time < 0.01 ~ "**",
      !is.na(FDR_within_time) &
        FDR_within_time < 0.05 ~ "*",
      TRUE ~ ""
    )
  )


# ------------------------------------------------------------
# 13. Cell-type order and UMAP-matched colors
# ------------------------------------------------------------

celltype_order <- c(
  # B cells
  "B_Memory",
  "B_Memory_IKZF3⁺",
  "B_Memory_JAM3⁺",
  "B_Memory_PDE4D⁺",
  "B_Memory_SOX5⁺",
  "B_Naive",
  "B_Naive_ATP1B3⁺",
  "B_Naive_CARMIL1⁺",
  "B_Naive_STX7⁺",
  "B_Naive_TCF4⁺",
  "B_Plasma",

  # CD4 T cells
  "CD4.Naive ARFGEF2+",
  "CD4.Naive GRAP2+",
  "CD4.Naive RPS13+",
  "CD4.Tcm",
  "CD4.Th1",
  "CD4.Th17 TAB2+",
  "CD4.Th17 TGFBR3+",
  "CD4.Treg",

  # CD8 T cells
  "CD8.Naive IL6ST+",
  "CD8.Naive JUNB+",
  "CD8.Tcm CTBP2+",
  "CD8.Tcm PYHIN1+",
  "CD8.Tem",
  "CD8.Temra",
  "CD8.Tm IFNG+",
  "CD8.Tm STAT4+",
  "MAIT"
)

subtype_cols <- c(
  # B cells
  "B_Memory"            = "#4E79A7",
  "B_Memory_IKZF3⁺"     = "#E15759",
  "B_Memory_JAM3⁺"      = "#59A14F",
  "B_Memory_PDE4D⁺"     = "#F28E2B",
  "B_Memory_SOX5⁺"      = "#B07AA1",
  "B_Naive"             = "#76B7B2",
  "B_Naive_ATP1B3⁺"     = "#EDC948",
  "B_Naive_CARMIL1⁺"    = "#FF9DA7",
  "B_Naive_STX7⁺"       = "#9C755F",
  "B_Naive_TCF4⁺"       = "#BAB0AC",
  "B_Plasma"            = "#D62728",

  # CD4 T cells
  "CD4.Naive ARFGEF2+"  = "#1F77B4",
  "CD4.Naive GRAP2+"    = "#2CA02C",
  "CD4.Naive RPS13+"    = "#FF7F0E",
  "CD4.Tcm"             = "#9467BD",
  "CD4.Th1"             = "#8C564B",
  "CD4.Th17 TAB2+"      = "#E377C2",
  "CD4.Th17 TGFBR3+"    = "#BCBD22",
  "CD4.Treg"            = "#17BECF",

  # CD8 T cells
  "CD8.Naive IL6ST+"    = "#AEC7E8",
  "CD8.Naive JUNB+"     = "#98DF8A",
  "CD8.Tcm CTBP2+"      = "#FFBB78",
  "CD8.Tcm PYHIN1+"     = "#C5B0D5",
  "CD8.Tem"             = "#C49C94",
  "CD8.Temra"           = "#F7B6D2",
  "CD8.Tm IFNG+"        = "#DBDB8D",
  "CD8.Tm STAT4+"       = "#9EDAE5",
  "MAIT"                = "#7F7F7F"
)

celltype_order <- celltype_order[
  celltype_order %in%
    unique(heatmap_log2fc_df$celltype)
]

missing_colors <- setdiff(
  celltype_order,
  names(subtype_cols)
)

if (length(missing_colors) > 0) {
  stop(
    "The following cell types are missing colors: ",
    paste(missing_colors, collapse = ", ")
  )
}


# ------------------------------------------------------------
# 14. Major-lineage annotation
# ------------------------------------------------------------

assign_major_group <- function(x) {

  dplyr::case_when(
    grepl("^B_", x) ~ "B cells",
    grepl("^CD4\\.", x) ~ "CD4 T cells",
    grepl("^CD8\\.", x) |
      x == "MAIT" ~ "CD8 T cells",
    TRUE ~ "Other"
  )
}

row_annotation_df <- tibble(
  celltype = celltype_order
) %>%
  mutate(
    major_group = assign_major_group(celltype),
    order_id = seq_along(celltype_order)
  )

major_bar_cols <- c(
  "B cells" = "#4472C4",
  "CD4 T cells" = "#70AD47",
  "CD8 T cells" = "#ED7D31"
)


# ------------------------------------------------------------
# 15. Total cell counts
# ------------------------------------------------------------

count_df <- tibble(
  celltype = rownames(counts_mat),
  total_cells = rowSums(
    counts_mat,
    na.rm = TRUE
  )
) %>%
  mutate(
    count_label = case_when(
      total_cells >= 1000000 ~ paste0(
        round(total_cells / 1000000, 1),
        "M"
      ),

      total_cells >= 1000 ~ paste0(
        round(total_cells / 1000, 1),
        "k"
      ),

      TRUE ~ as.character(total_cells)
    )
  ) %>%
  filter(
    celltype %in% celltype_order
  )

# Verify counts against Seurat metadata
seurat_count_df <- as.data.frame(
  table(
    immune_merged$celltype_merged
  )
) %>%
  transmute(
    celltype = as.character(Var1),
    seurat_total_cells = as.numeric(Freq)
  )

count_check <- count_df %>%
  left_join(
    seurat_count_df,
    by = "celltype"
  ) %>%
  mutate(
    difference =
      total_cells -
      seurat_total_cells
  )

if (!all(count_check$difference == 0)) {
  warning(
    "Propeller counts and Seurat metadata counts are not identical for all cell types."
  )
}


# ------------------------------------------------------------
# 16. Prepare plotting data
# ------------------------------------------------------------

heatmap_log2fc_df <- heatmap_log2fc_df %>%
  filter(
    celltype %in% celltype_order,
    time %in% present_times
  ) %>%
  mutate(
    celltype = factor(
      celltype,
      levels = rev(celltype_order)
    ),
    time = factor(
      time,
      levels = present_times
    )
  )

count_df <- count_df %>%
  mutate(
    celltype = factor(
      celltype,
      levels = rev(celltype_order)
    )
  )

row_annotation_df <- row_annotation_df %>%
  mutate(
    celltype = factor(
      celltype,
      levels = rev(celltype_order)
    ),
    major_group = factor(
      major_group,
      levels = c(
        "B cells",
        "CD4 T cells",
        "CD8 T cells"
      )
    )
  )

group_label_df <- row_annotation_df %>%
  group_by(
    major_group
  ) %>%
  slice_min(
    abs(
      order_id -
        median(order_id)
    ),
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup()

time_labels <- setNames(
  paste0(
    "D",
    sub(
      "^day",
      "",
      present_times
    )
  ),
  present_times
)

log2fc_limit <- 1.5


# ------------------------------------------------------------
# 17. Plot major-lineage annotation bar
# ------------------------------------------------------------

p_major_annotation <- ggplot(
  row_annotation_df,
  aes(
    x = 1,
    y = celltype,
    fill = major_group
  )
) +
  geom_tile(
    width = 0.85,
    height = 0.97,
    color = "white",
    linewidth = 0.2
  ) +
  geom_text(
    data = group_label_df,
    aes(
      x = 1,
      y = celltype,
      label = major_group
    ),
    inherit.aes = FALSE,
    angle = 90,
    fontface = "bold",
    size = 2.6,
    color = "white"
  ) +
  scale_fill_manual(
    values = major_bar_cols,
    guide = "none"
  ) +
  scale_y_discrete(
    limits = rev(celltype_order),
    drop = FALSE,
    expand = c(0, 0)
  ) +
  theme_void() +
  theme(
    plot.margin = margin(
      t = 32,
      r = 0,
      b = 24,
      l = 0
    )
  )


# ------------------------------------------------------------
# 18. Plot log2FC heatmap
# ------------------------------------------------------------

p_propeller_heatmap_log2fc <- ggplot(
  heatmap_log2fc_df,
  aes(
    x = time,
    y = celltype,
    fill = mean_log2FC
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.4
  ) +

  # Separate prime and booster phases between Day 30 and Day 31
  geom_vline(
    xintercept = 6.5,
    color = "grey45",
    linewidth = 0.45
  ) +

  geom_text(
    aes(
      label = significance
    ),
    size = 2.8,
    fontface = "bold",
    color = "black"
  ) +

  scale_fill_gradient2(
    low = "#365B8C",
    mid = "#F7F7F7",
    high = "#D95F59",
    midpoint = 0,
    limits = c(
      -log2fc_limit,
      log2fc_limit
    ),
    breaks = c(
      -1.5,
      0,
      1.5
    ),
    labels = c(
      "-1.5",
      "0",
      "1.5"
    ),
    oob = scales::squish,
    name = expression(
      log[2] * "FC vs Day 0"
    )
  ) +

  scale_x_discrete(
    labels = time_labels,
    drop = FALSE,
    expand = c(0, 0)
  ) +

  scale_y_discrete(
    limits = rev(celltype_order),
    drop = FALSE,
    expand = c(0, 0)
  ) +

  guides(
    fill = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = grid::unit(
        28,
        "mm"
      ),
      barheight = grid::unit(
        3.2,
        "mm"
      ),
      ticks = TRUE,
      frame.colour = "black",
      ticks.colour = "black"
    )
  ) +

  labs(
    x = NULL,
    y = NULL
  ) +

  theme_classic(
    base_size = 9
  ) +

  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),

    axis.text.y = element_text(
      size = 7.2,
      color = "black",
      margin = margin(
        r = 4
      )
    ),

    axis.text.x = element_text(
      size = 8,
      color = "black"
    ),

    legend.position = "top",
    legend.justification = "left",

    legend.title = element_text(
      size = 8,
      margin = margin(
        b = 2
      )
    ),

    legend.text = element_text(
      size = 7
    ),

    panel.grid = element_blank(),

    plot.margin = margin(
      t = 2,
      r = 2,
      b = 24,
      l = 2
    )
  )


# ------------------------------------------------------------
# 19. Plot total cell counts
# ------------------------------------------------------------

p_cell_counts <- ggplot(
  count_df,
  aes(
    x = total_cells,
    y = celltype,
    fill = celltype
  )
) +
  geom_col(
    width = 0.73,
    color = NA
  ) +

  geom_text(
    aes(
      label = count_label
    ),
    hjust = -0.10,
    size = 2.5,
    color = "black"
  ) +

  scale_fill_manual(
    values = subtype_cols,
    guide = "none"
  ) +

  scale_y_discrete(
    limits = rev(celltype_order),
    drop = FALSE,
    expand = c(0, 0)
  ) +

  scale_x_sqrt(
    labels = scales::label_number(
      scale_cut = scales::cut_short_scale()
    ),
    expand = expansion(
      mult = c(
        0,
        0.30
      )
    )
  ) +

  labs(
    x = "Cell count",
    y = NULL
  ) +

  theme_classic(
    base_size = 8
  ) +

  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),

    axis.line.x = element_line(
      linewidth = 0.35
    ),

    axis.text.x = element_text(
      size = 6.8
    ),

    axis.title.x = element_text(
      size = 7.5,
      margin = margin(
        t = 4
      )
    ),

    plot.margin = margin(
      t = 32,
      r = 14,
      b = 6,
      l = 2
    )
  )


# ------------------------------------------------------------
# 20. Assemble final figure
# ------------------------------------------------------------

p_propeller_final_log2fc <-
  p_major_annotation +
  p_propeller_heatmap_log2fc +
  p_cell_counts +
  plot_layout(
    widths = c(
      0.075,
      1,
      0.28
    )
  )

print(
  p_propeller_final_log2fc
)


# ------------------------------------------------------------
# 21. Save figures
# ------------------------------------------------------------

ggsave(
  filename = file.path(
    outdir,
    "Propeller_log2FC_vs_day0_heatmap.pdf"
  ),
  plot = p_propeller_final_log2fc,
  width = 185,
  height = 120,
  units = "mm",
  limitsize = FALSE
)

ggsave(
  filename = file.path(
    outdir,
    "Propeller_log2FC_vs_day0_heatmap.tiff"
  ),
  plot = p_propeller_final_log2fc,
  width = 185,
  height = 120,
  units = "mm",
  dpi = 500,
  compression = "lzw",
  limitsize = FALSE
)


# ------------------------------------------------------------
# 22. Save p_propeller_final_log2fc as RDS
# ------------------------------------------------------------

saveRDS(
  p_propeller_final_log2fc,
  file = file.path(
    outdir,
    "p_propeller_final_log2fc.rds"
  ),
  compress = TRUE
)


# Example for reloading and displaying the saved plot:
#
# p_propeller_final_log2fc <- readRDS(
#   file.path(
#     outdir,
#     "p_propeller_final_log2fc.rds"
#   )
# )
#
# print(p_propeller_final_log2fc)


# ------------------------------------------------------------
# 23. Save result tables
# ------------------------------------------------------------

write.csv(
  transformed_long,
  file.path(
    outdir,
    "01_Propeller_transformed_proportions.csv"
  ),
  row.names = FALSE
)

write.csv(
  delta_long,
  file.path(
    outdir,
    "02_Propeller_donor_transformed_delta_vs_day0.csv"
  ),
  row.names = FALSE
)

write.csv(
  stats_df,
  file.path(
    outdir,
    "03_Propeller_limma_statistics_vs_day0.csv"
  ),
  row.names = FALSE
)

write.csv(
  props_long,
  file.path(
    outdir,
    "04_Propeller_raw_proportions.csv"
  ),
  row.names = FALSE
)

write.csv(
  log2fc_long,
  file.path(
    outdir,
    "05_Propeller_donor_log2FC_vs_day0.csv"
  ),
  row.names = FALSE
)

write.csv(
  heatmap_log2fc_df,
  file.path(
    outdir,
    "06_Propeller_mean_log2FC_heatmap_data.csv"
  ),
  row.names = FALSE
)

write.csv(
  count_df,
  file.path(
    outdir,
    "07_Celltype_total_counts.csv"
  ),
  row.names = FALSE
)

write.csv(
  count_check,
  file.path(
    outdir,
    "08_Celltype_count_QC.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 24. Save complete analysis results
# ------------------------------------------------------------

saveRDS(
  list(
    meta_prop = meta_prop,
    sample_info = sample_info,
    transformed_matrix = tp_mat,
    proportions_matrix = props_mat,
    counts_matrix = counts_mat,
    transformed_long = transformed_long,
    transformed_delta_long = delta_long,
    limma_statistics = stats_df,
    proportions_long = props_long,
    donor_log2FC_long = log2fc_long,
    heatmap_log2FC_data = heatmap_log2fc_df,
    count_data = count_df,
    count_QC = count_check,
    subtype_colors = subtype_cols,
    major_lineage_colors = major_bar_cols,
    final_plot = p_propeller_final_log2fc
  ),
  file = file.path(
    outdir,
    "Propeller_log2FC_complete_results.rds"
  ),
  compress = TRUE
)

message(
  "\nAnalysis complete.\nOutput directory: ",
  outdir,
  "\nSaved plot RDS: ",
  file.path(
    outdir,
    "p_propeller_final_log2fc.rds"
  )
)
