# 00_paths.R
# Every script in this folder sources this file. Run the scripts from the imss_paper folder,
# for example: Rscript 3_code/scripts/03_main_results.R
#
# RAW points to the original project data. The data stay outside this package on purpose:
# the raw files contain personal information and must not be uploaded anywhere.

RAW     <- "C:/Users/pablo/OneDrive/Documentos/pablo/RStudio/code/2023/imss/imss/1_data/bases/base_paper/"
OUT_TAB <- "3_code/output/tables/"
OUT_FIG <- "3_code/output/figures/"
OUT_CSV <- "3_code/output/csv/"
# de-identified analysis file (one row per randomised person), used when RAW is not available
PUBLIC_DATA <- "3_code/data/analysis_data.csv"
for (p in c(OUT_TAB, OUT_FIG, OUT_CSV)) dir.create(p, showWarnings = FALSE, recursive = TRUE)

# Table and figure images carry no title and no notes (Pablo, 30 October 2026): the manuscript gives
# each one its caption and, where needed, a short note as text. The title and notes the scripts write
# are kept next to each image as <image>_note.txt. Set to TRUE to draw them inside the images again.
IMAGE_TEXT <- FALSE
save_note <- function(file, title, notes) {
  writeLines(c(title, "", notes), sub("\\.png$", "_note.txt", file))
}
