# How do long-term care services impact older adults' quality of life?

Data and code for the article "How do long-term care services impact older adults' quality of life? A quasi-experimental study exploring effects over time and across domains", by Pablo D. Alcala and Guillermo M. Cejudo.

The study evaluates the first day care centre for dependent older adults of the Mexican Institute for Social Security (IMSS), in Mexico City, using the randomised invitations IMSS made to the 120 older adults who expressed interest (eight rounds in 2023) and three survey waves.

## Contents

| Folder | Contents |
|---|---|
| `3_code/scripts/` | R scripts, numbered in the order they run |
| `3_code/data/analysis_data.csv` | De-identified analysis file: one row per randomised person (120) |
| `3_code/output/` | Tables and figures of the article (`tables/`, `figures/`), their full notes (`*_note.txt`) and the numbers behind them (`csv/`) |
| `submission/` | Data and script of the first submission, kept for the record |

## The analysis file

`analysis_data.csv` has one row per person: a random id (1 to 120), the randomised invitation (`z`), the WHOQOL-BREF score and its five domains at baseline, midline and endline (`_y0`, `_y1`, `_y2`; WHO 0 to 100 scale, 26 items), attendance by each follow-up survey (`d1`, `d2`), days, weeks and days per week attended (`days`, `weeks`, `dpw`), and baseline characteristics: age in whole years, sex, indigenous status, household size, education level (1 to 10), whether the main caregiver is a son or daughter, and the Katz and Lawton scores. It contains no names and no dates. It is written by `11_public_data.R` from the raw files, which contain personal information and are not public. The use of the data follows a protocol approved by the IMSS National Scientific Research Committee (PI-2023-3458).

## Running the code

Run from the repository folder, for example `Rscript 3_code/scripts/03_main_results.R`. Packages: tidyverse, AER, sandwich, lmtest, fixest, lmerTest, gridExtra, gtable, ragg, readxl, lubridate and stringi.

Without the raw files, these scripts run from the analysis file:

| Script | Produces |
|---|---|
| `03_main_results.R` | Tables 3, 4 and 5 |
| `05_robustness.R` | Tables A.8 and A.9 |
| `06_heterogeneity.R` | Table A.10 |
| `10_dose_table.R` | Table A.11 |

Tables 3, 4, 5 and A.9 reproduce exactly. Tables A.8 and A.10 use age, which the public file rounds to whole years, so a few of their numbers differ slightly from the article (by at most 0.07 points in Table A.8; in Table A.10 the split by age moves a few people).

Scripts 01, 02, 04, 07, 08 and 09 need the survey dates, the attendance register or the item-level answers, so they run only with the raw data; their outputs are in `3_code/output/`.
