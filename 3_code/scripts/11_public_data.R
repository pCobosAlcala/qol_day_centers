# 11_public_data.R
#
# Writes the de-identified analysis file for the public repository: 3_code/data/analysis_data.csv.
# It is build_analysis_data() (one row per randomised person, 120 rows, no dates) with
#   - the study folio replaced by a random number from 1 to 120 (the key is not saved anywhere);
#   - age at baseline rounded to whole years.
# Nothing else is changed: invitation, attendance, days, weeks and days per week by wave, WHOQOL-BREF and
# domain scores by wave (WHO 0 to 100 scale), sex, indigenous status, household size, education level,
# main caregiver and baseline Katz and Lawton scores. Names, birth dates and survey or entry dates are
# never in it. Scripts 03, 05, 06 and 10 run from this file when the raw data are not available.
#
#   Rscript 3_code/scripts/11_public_data.R          (run from the imss_paper folder, with the raw data)

rm(list = ls())
suppressPackageStartupMessages(library(tidyverse))
source("3_code/scripts/00_paths.R")
source("3_code/scripts/00_data.R")

dat <- build_analysis_data(quiet = TRUE)
stopifnot(nrow(dat) == 120, !any(grepl("date", names(dat))))

set.seed(20261030)
pub <- dat %>%
  mutate(id = sample(nrow(dat)), age0 = round(age0)) %>%
  arrange(id)

dir.create(dirname(PUBLIC_DATA), showWarnings = FALSE, recursive = TRUE)
write_csv(pub, PUBLIC_DATA)
cat(sprintf("written %s: %d people, %d columns\n", PUBLIC_DATA, nrow(pub), ncol(pub)))
cat("columns:", paste(names(pub), collapse = ", "), "\n")
