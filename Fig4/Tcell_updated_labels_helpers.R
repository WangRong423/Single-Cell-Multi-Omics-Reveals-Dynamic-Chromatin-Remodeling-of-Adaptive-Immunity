# Shared annotation mapping. Original project files are not overwritten.
UPDATED_BASE <- "/public/users/wangrong/home/wangrong/results/T_B/TB_final_downstream"
UPDATED_T_TYPES <- c("CD4_NaiveT", "CD4_Naive_CCR4", "CD8_Memory_GZMK",
                     "Treg", "CD4_Memory_ZEB2", "CD4_Memory", "CD8_NaiveT",
                     "CD8_Memory", "gdT", "CD4_CD8_DP")
UPDATED_B_TYPES <- c("Naive_B", "Intermediate_B", "Memory_B", "PlasmaB")
SUBTYPE_RENAME_MAP <- c(
  "CD4_Tcm" = "CD4_Naive_CCR4",
  "CD4_NaiveT2" = "CD4_Naive_CCR4",
  "CD4_Teff" = "CD8_Memory_GZMK",
  "CD8_Teff" = "CD8_Memory_GZMK",
  "Th1" = "CD4_Memory_ZEB2",
  "Th2" = "CD4_Memory",
  "CD8_Tem" = "CD8_Memory",
  "Activated_B" = "Intermediate_B",
  "Activaed_B" = "Intermediate_B",
  "Resting_B" = "Memory_B"
)
rename_tb_labels <- function(x, audit_file = NULL) {
  old <- as.character(x)
  new <- old
  i <- match(old, names(SUBTYPE_RENAME_MAP))
  new[!is.na(i)] <- unname(SUBTYPE_RENAME_MAP[i[!is.na(i)]])
  map <- unique(data.frame(Original=old, Updated=new))
  known <- map[!is.na(map$Updated) & map$Updated %in% c(UPDATED_T_TYPES, UPDATED_B_TYPES), , drop=FALSE]
  if (anyDuplicated(known$Updated)) {
    print(known)
    stop("Multiple existing groups map to one updated subtype. Verify cell membership before merging.")
  }
  unknown <- setdiff(unique(new[!is.na(new) & nzchar(new)]), c(UPDATED_T_TYPES, UPDATED_B_TYPES))
  if (length(unknown)) warning("Labels outside the T/B analysis: ", paste(unknown, collapse=", "))
  if (!is.null(audit_file)) {
    dir.create(dirname(audit_file), recursive=TRUE, showWarnings=FALSE)
    write.csv(map, audit_file, row.names=FALSE)
  }
  new
}
select_tb_project <- function() {
  for (nm in c("TB_proj", "proj_TB_deviations_complete", "proj_TB_fig3", "proj_all")) {
    if (exists(nm, envir=.GlobalEnv, inherits=FALSE)) {
      p <- get(nm, envir=.GlobalEnv)
      if (inherits(p, "ArchRProject")) {
        message("Using ArchR project: ", nm)
        return(p)
      }
    }
  }
  stop("Load TB_proj into the R session first, then source this script.")
}
