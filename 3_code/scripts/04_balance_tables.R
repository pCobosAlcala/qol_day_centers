# 04_balance_tables.R  (revision, October 2026; data core and 0 to 100 scale, October 2026)
#
# Who was invited, who attended and who stayed: the balance, take up and attrition tables of the
# revision, plus the numbers behind the Background sentence that describes the people who attended.
#
#   table2_balance_baseline   MAIN TEXT. Invited (91) versus never invited (29): all 120 people
#                             interviewed at baseline, each variable on everyone who has it.
#   tableA_balance_samples    ANNEX. The same comparison inside the three estimation samples:
#                             midline (102), endline (85) and all three waves (82).
#   tableA_attended           ANNEX. Among the invited, attended versus did not attend: all invited
#                             people (attendance over the whole register), and the invited people in
#                             the midline and endline samples (attendance up to that wave's survey,
#                             the definition of the main results).
#   tableA_attrition          ANNEX. Retained versus not retained: at midline (102 vs 18), at endline
#                             (85 vs 35) and over all three waves (82 vs 38).
#   tableA_attrition_rates    ANNEX. Attrition rates at midline and endline by group (attended,
#                             invited and did not attend, never invited) and by arm, Fisher exact tests.
#   background_participants   CSV only. Age, sex, caregiver and education of those who attended,
#                             checked against the figures quoted in the submitted Background.
#   text_numbers_balance      CSV only. Every number the text or the response letter may quote from
#                             this script, with a section column.
#
# Data. Everything comes from build_analysis_data() in 00_data.R (one row per randomised person,
# rebuilt from the raw files and checked there against the files of the submitted paper), except
# three things the core does not return, read here from the same raw files and checked against it:
#   (a) the caregiver's relationship as answered (pm_limpieza1.csv, x06_03_pc): needed for the
#       share cared for by a partner and for the coding of the submitted tables; the "son or
#       daughter" flag built from it must equal the core's child_carer.
#   (b) the round in which each invited person was first picked (t_pm.csv): timing check only.
#   (c) total days in the centre's attendance register over its whole period: the core counts days
#       only up to each wave's survey, which does not exist for people not surveyed at that wave.
#       Attrition by take up needs a group for everyone, so "attended" there means more than five
#       days in the whole register. The register is read exactly as in 00_data.R, and the days up
#       to each survey recomputed from it must equal the core's days1 and days2 for every person.
#
# Changes against the version of 2 October (copies of the old script and outputs are in the
# _previous subfolders, suffix _before_00data):
#   - WHOQOL-BREF and environmental rows on the WHO 0 to 100 scale (25 x (item mean minus 1)),
#     shown with one decimal; differences are 25 times the 1 to 5 ones, p values unchanged.
#   - The three invited people who entered the centre after their midline survey now have scores
#     (rebuilt from the raw items), so every row of Table 2 except Katz, Lawton and caregiver uses
#     all 120 and the midline sample is 102 (was 99). They count as retained at midline.
#   - Katz and Lawton come from the raw items (the three people above now have them); one invited
#     person lacks Katz, Lawton and the caregiver record.
#   - Attendance: Panels B and C of the take up table use the wave specific attendance of the main
#     results; Panel A, and the attrition rates by group, use the whole register (point (c)).
#     The old time invariant flag from pm_vf.csv (28 people) is no longer used.
#   - Attrition gets a midline panel and a table of rates by group and by arm. The two midline
#     definitions of the old version (not in the analysis sample, not surveyed) now coincide.
#   - The Background numbers are those of the attenders in the endline sample, and a check stops
#     the script unless they equal Panel C of the take up table.
#
# Choices kept from the referee check of October 2026:
#   - "Cared for by a son or daughter" counts caregivers who answered "Padre o madre" OR "Hija/o"
#     (the core's child_carer). The old descriptive code treats "Hija/o" as a child caregiver (the
#     older adults are about 80, so the caregiver cannot be their parent), but control.csv was built
#     BEFORE that recode, so those caregivers were "not a child" in the submitted balance tables.
#     Missing caregiver records stay missing (control.csv had coded the one missing record as 0).
#   - Tests follow the submission: Welch t-tests for continuous variables, chi-squared tests (with
#     the default continuity correction) for binary ones, shown in per cent. Fisher exact and
#     uncorrected chi-squared p values are computed for every binary row and written to the csv
#     files; the notes name any binary row where the tests fall on different sides of 0.05 or 0.10.
#   - Timing check: one invited person's baseline record is dated after randomisation round 1
#     (whose date is still to be confirmed with IMSS); Table 2 is recomputed without them.
#
# Aggregates only: the raw files contain names, birth dates and exact dates; nothing at person level
# is written, and the survey dates returned by the core stay inside this script. Run from the
# imss_paper folder:
#   "C:/Program Files/R/R-4.6.1/bin/Rscript.exe" 3_code/scripts/04_balance_tables.R

rm(list = ls())
suppressPackageStartupMessages({
  library(tidyverse)
  library(gridExtra); library(gtable); library(grid)
})
options(scipen = 999, width = 200)

source("3_code/scripts/00_paths.R")
source("3_code/scripts/00_data.R")

# ---- 1. data: one row per randomised person (120) ----------------------------------------------
# dates = TRUE: survey dates are needed for the timing check and to check the register read in
# 1c against the core. They are never written out.
dat <- build_analysis_data(dates = TRUE)

ascii  <- function(x) stringi::stri_trans_general(x, "Latin-ASCII")     # "companera", no accents
# small counts in words for the notes ("Three invited people", "ten comparisons")
num_word <- function(n) if (n <= 10) c("no", "one", "two", "three", "four", "five", "six", "seven",
                                       "eight", "nine", "ten")[n + 1] else as.character(n)
cap <- function(s) paste0(toupper(substr(s, 1, 1)), substring(s, 2))
rd_chr <- function(f) read_csv(paste0(RAW, f), col_types = cols(.default = col_character()),
                               name_repair = "unique_quiet", progress = FALSE)

# 1a. caregiver relationship as answered at baseline (point (a) of the header)
care <- rd_chr("pm_limpieza1.csv") %>%
  filter(name == "x06_03_pc") %>%
  transmute(id = x03_01, care_lab = ascii(t0)) %>%
  mutate(child_lab       = case_when(is.na(care_lab) ~ NA_integer_,
                                     care_lab %in% c("Padre o madre", "Hija/o") ~ 1L, TRUE ~ 0L),
         partner         = case_when(is.na(care_lab) ~ NA_integer_,
                                     care_lab == "Esposa/o o companera/o" ~ 1L, TRUE ~ 0L),
         # coding of control.csv and the submitted tables: "Padre o madre" only, missing = 0
         child_submitted = as.integer(care_lab %in% "Padre o madre"),
         hija            = care_lab %in% "Hija/o")
stopifnot(nrow(care) == 120, !anyDuplicated(care$id), any(care$partner %in% 1))  # label found

# 1b. round in which each invited person was first picked. In t_pm.csv, T1 to T8 equal 1 in the
# round a person is picked (NA afterwards, 0 before), as built in 2_bases_tratamiento.R.
rnd <- rd_chr("t_pm.csv") %>% select(id = x03_01, all_of(paste0("T", 1:8))) %>%
  pivot_longer(-id, names_to = "r", values_to = "picked") %>%
  filter(picked %in% "1") %>%
  group_by(id) %>% summarise(round = min(as.integer(sub("T", "", r))), .groups = "drop")
# Round 1 date: file "Seleccion aleatoria 2ene22.xlsx" (year typo for 2023). Rounds 2 to 8 came
# later in 2023. NOT YET CONFIRMED WITH IMSS: used only for the timing check of Table 2.
ROUND1_DATE <- as.Date("2023-01-02")

# 1c. days in the attendance register over its whole period (point (c) of the header). Same reading
# as 00_data.R section 5.3: one column per opening day, header = Excel date serial, "x" = attended.
reg    <- readxl::read_excel(paste0(RAW, "CD CASSAAM 16febrero24 asistencias.xlsx"), .name_repair = "minimal")
serial <- suppressWarnings(as.numeric(names(reg)))
is_day <- !is.na(serial) & serial > 40000 & serial < 50000
stopifnot(sum(is_day) == 223, sum(names(reg) == "Folio") == 1)          # 223 opening days
day_of <- as.Date(serial[is_day], origin = "1899-12-30")
marks  <- reg[, c(which(names(reg) == "Folio"), which(is_day))]
names(marks) <- c("id", as.character(day_of))
visits <- marks %>%
  pivot_longer(-id, names_to = "day", values_to = "mark", values_transform = list(mark = as.character)) %>%
  filter(tolower(mark) %in% "x") %>%
  transmute(id, day = as.Date(day))
reg_days <- visits %>%
  left_join(dat %>% select(id, date1, date2), by = "id") %>%
  group_by(id) %>%
  summarise(days_reg = n(),                                     # whole register
            chk1 = sum(day <= date1[1]), chk2 = sum(day <= date2[1]), .groups = "drop")
# register period, for the notes (month names written out so the locale does not matter)
fmt_day <- function(x) sprintf("%d %s %d", as.integer(format(x, "%d")), month.name[as.integer(format(x, "%m"))],
                               as.integer(format(x, "%Y")))
reg_period <- sprintf("%s to %s", fmt_day(min(day_of)), fmt_day(max(day_of)))

# ---- 1d. person level analysis file -----------------------------------------------------------
p <- dat %>%
  left_join(care %>% select(id, child_lab, partner, child_submitted, hija), by = "id") %>%
  left_join(rnd, by = "id") %>%
  left_join(reg_days, by = "id") %>%
  mutate(days_reg = replace_na(days_reg, 0L),
         ever     = as.integer(days_reg > 5),                   # attended more than five days, whole register
         # estimation samples of 03: baseline WHOQOL-BREF plus that wave
         in_mid   = !is.na(whoqol_y0) & !is.na(whoqol_y1),
         in_end   = !is.na(whoqol_y0) & !is.na(whoqol_y2),
         in_all3  = in_mid & in_end,
         # take up group for the attrition rates
         grp      = case_when(z == 0 ~ "never", ever == 1 ~ "attended", TRUE ~ "invited_no"),
         # timing: baseline interview dated after the invitation (only round 1 can qualify, guard 4)
         late_base = z == 1 & round %in% 1 & date0 > ROUND1_DATE)

# guard 1: the caregiver flag built from the labels is the core's child_carer
stopifnot(identical(p$child_lab, p$child_carer))
# guard 2: the register read here gives the core's days up to each survey, for every person
stopifnot(identical(as.integer(replace_na(p$chk1, 0L))[!is.na(p$days1)], as.integer(p$days1[!is.na(p$days1)])),
          identical(as.integer(replace_na(p$chk2, 0L))[!is.na(p$days2)], as.integer(p$days2[!is.na(p$days2)])))
# guard 3: the counts and identities every output below relies on
stopifnot(nrow(p) == 120, sum(p$z) == 91,
          all(p$ever[p$z == 0] == 0),                           # nobody never invited attended
          all(p$ever[p$d1 %in% 1] == 1), all(p$ever[p$d2 %in% 1] == 1),
          sum(p$in_mid) == 102, sum(p$in_end) == 85, sum(p$in_all3) == 82,
          identical(p$in_mid, !is.na(p$date1)),                 # analysis sample = surveyed, so the
          identical(p$in_end, !is.na(p$date2)),                 # two old midline definitions coincide
          identical(is.na(p$d1), !p$in_mid), identical(is.na(p$d2), !p$in_end))
# guard 4: every invited person has a selection round and nobody else does; every baseline of a
# person picked in round 2 or later precedes round 1, so only round 1 can follow a baseline
stopifnot(all(!is.na(p$round[p$z == 1])), all(is.na(p$round[p$z == 0])),
          max(p$date0[p$round %in% 2:8]) < ROUND1_DATE)

# who is in the "attended" group of the whole register but not attending in the main results.
# Against midline: began after their midline survey, or not surveyed at midline. Against endline:
# only people not surveyed at endline (guard: nobody surveyed at endline began after that survey).
n_ever       <- sum(p$ever)
n_ever_late  <- sum(p$ever == 1 & p$d1 %in% 0)                  # began after their midline survey
n_ever_nomid <- sum(p$ever == 1 & !p$in_mid)                    # not surveyed at midline
n_ever_nofu  <- sum(p$ever == 1 & !p$in_mid & !p$in_end)        # no survey after baseline
n_ever_noend <- sum(p$ever == 1 & !p$in_end)                    # not surveyed at endline
stopifnot(n_ever == sum(p$d1, na.rm = TRUE) + n_ever_late + n_ever_nomid,
          sum(p$ever == 1 & p$d2 %in% 0) == 0,
          n_ever == sum(p$d2, na.rm = TRUE) + n_ever_noend)
ever_txt <- sprintf(paste(
  "Of the %d people who attended over the whole register, %d were attending by their midline survey (%s began",
  "attending after that survey and %s %s not surveyed at midline) and %d by their endline survey (the other %s",
  "were not surveyed at endline)."),
  n_ever, sum(p$d1, na.rm = TRUE), num_word(n_ever_late), num_word(n_ever_nomid),
  if (n_ever_nomid == 1) "was" else "were", sum(p$d2, na.rm = TRUE), num_word(n_ever_noend))
cat("people:", nrow(p), "| invited:", sum(p$z),
    "| attended more than five days, whole register:", n_ever,
    "| attending by midline / endline survey:", sum(p$d1, na.rm = TRUE), "/", sum(p$d2, na.rm = TRUE),
    "| midline / endline / all three samples:", sum(p$in_mid), sum(p$in_end), sum(p$in_all3), "\n")

# ---- 2. balance helper ------------------------------------------------------------------------
# noun: how the variable is named in the sentence about missing data
vars <- tribble(
  ~v,                 ~label,                                        ~type, ~noun,
  "age0",             "Age at baseline (years)",                     "num", "age",
  "hh_size",          "People living with the older adult",           "num", "household size",
  "educ",             "Education (1 to 10)",                         "num", "education",
  "lawton0",          "Autonomy (Lawton, 0 to 8)",                   "num", "Lawton score",
  "katz0",            "Independence (Katz, 0 to 6)",                 "num", "Katz score",
  "female",           "Female (%)",                                  "bin", "sex",
  "indigenous",       "Indigenous (%)",                              "bin", "indigenous status",
  "child_carer",      "Cared for by a son or daughter (%)",          "bin", "caregiver record",
  "whoqol_y0",        "WHOQOL-BREF at baseline (0 to 100)",          "num", "baseline WHOQOL-BREF",
  "environmental_y0", "Environmental domain at baseline (0 to 100)", "num", "baseline environmental score")
K <- nrow(vars)                                          # comparisons per panel, for Bonferroni

# one row per variable: group means (per cent for binary rows), difference, test, n per group.
# g is a 0/1 column; group 1 is shown first (invited, attended, retained).
balance <- function(df, g) {
  pmap_dfr(vars %>% select(v, label, type), function(v, label, type) {
    s <- df %>% filter(!is.na(.data[[v]]), !is.na(.data[[g]]))
    a <- s[[v]][s[[g]] == 1]; b <- s[[v]][s[[g]] == 0]
    pf <- NA_real_; pu <- NA_real_                       # Fisher and uncorrected chi-squared: binary only
    if (type == "num") {
      pv <- t.test(a, b)$p.value; test <- "Welch t-test"; k <- 1
    } else {
      # 2 x 2 table built from the same observations (rows: value, columns: group), as in 01
      tb <- table(factor(c(a, b), levels = 0:1),
                  factor(rep(c(1, 0), times = c(length(a), length(b))), levels = 0:1))
      pv <- suppressWarnings(chisq.test(tb))$p.value; test <- "chi-squared"; k <- 100
      pf <- fisher.test(tb)$p.value
      pu <- suppressWarnings(chisq.test(tb, correct = FALSE))$p.value
    }
    tibble(variable = label, type = type, test = test, n1 = length(a), n0 = length(b),
           mean1 = k * mean(a), mean0 = k * mean(b), diff = k * (mean(a) - mean(b)), p = pv,
           p_fisher = pf, p_chisq_uncorrected = pu)
  })
}

# display version: one decimal for age, percentages and the 0 to 100 scores, two for the small scales
fmt_p <- function(x) ifelse(x < 0.001, "< 0.001", sprintf("%.3f", x))
fmt_balance <- function(b, lab1, lab0) {
  dg <- ifelse(b$type == "bin" | grepl("^Age|0 to 100", b$variable), 1, 2)
  f  <- function(x) sub("^-(0\\.0*)$", "\\1", sprintf(paste0("%.", dg, "f"), x))  # no "-0.00"
  tibble(Characteristic = b$variable, !!lab1 := f(b$mean1), !!lab0 := f(b$mean0),
         Difference = f(b$diff), `p value` = fmt_p(b$p), n = sprintf("%d / %d", b$n1, b$n0))
}

# lists in words for the notes (num_word and cap are defined in section 1)
and_list <- function(x, conj = "and") if (length(x) == 1) x else
  paste(paste(x[-length(x)], collapse = ", "), conj, x[length(x)])                # "a, b and c"

# Bonferroni sentence computed from the results, so the note can never contradict the table.
# Binary rows count as passing if either the chi-squared or the Fisher exact test passes.
bonf_text <- function(b, per_panel = FALSE) {
  thr <- 0.05 / K
  b <- b %>% mutate(p = pmin(p, p_fisher, na.rm = TRUE))
  pass <- b %>% filter(p < thr)
  head <- sprintf("With %s comparisons%s the Bonferroni threshold for a 5 per cent test is %.3f;",
                  num_word(K), if (per_panel) " per panel" else "", thr)
  tail <- if (nrow(pass) == 0) sprintf("no difference passes it under either test (the smallest p value is %.3f).",
                                       min(b$p)) else
    sprintf("%s pass%s it: %s.", num_word(nrow(pass)), if (nrow(pass) == 1) "es" else "",
            paste(sprintf("%s%s (p = %s)", tolower(sub(" \\(.*", "", pass$variable)),
                          if ("panel" %in% names(pass)) paste0(", ", pass$panel) else "",
                          fmt_p(pass$p)), collapse = "; "))
  paste(head, tail)
}

# binary rows where the chi-squared and Fisher exact tests fall on different sides of 0.05 or 0.10
disagree <- function(b) b %>% filter(type == "bin", (p < 0.05) != (p_fisher < 0.05) | (p < 0.10) != (p_fisher < 0.10))
fisher_text <- function(b) {
  dd <- disagree(b)
  if (nrow(dd) == 0) return("Fisher exact tests on the binary rows give the same conclusions at the 5 and 10 per cent levels.")
  paste0("Fisher exact tests change the conclusion for ",
         paste(sprintf("%s%s (chi-squared p = %s, Fisher exact p = %s)", tolower(sub(" \\(.*", "", dd$variable)),
                       if ("panel" %in% names(dd)) paste0(", ", dd$panel) else "", fmt_p(dd$p), fmt_p(dd$p_fisher)),
               collapse = "; "), ".")
}

# sentence on missing data, built from the data: who lacks what, by arm
miss_text <- function(df) {
  na_mat <- is.na(as.matrix(df[, vars$v]))
  pat <- apply(na_mat, 1, function(r) if (any(r)) and_list(vars$noun[r], "or") else NA_character_)
  tab <- tibble(z = df$z, pat = pat) %>% filter(!is.na(pat)) %>% count(z, pat)
  if (nrow(tab) == 0) return(sprintf("Every variable is observed for all %d people.", nrow(df)))
  paste0(paste(sprintf("%s %s %s no %s", cap(num_word(tab$n)), ifelse(tab$z == 1, "invited", "never invited"),
                       ifelse(tab$n == 1, "person has", "people have"), tab$pat), collapse = "; "),
         "; every other value is observed.")
}

# stack several formatted panels into one data frame, remembering where the panel labels go
stack_panels <- function(lst) {
  rows <- list(); hdr <- integer(0); k <- 0
  for (nm in names(lst)) {
    blank <- lst[[nm]][1, ]; blank[1, ] <- as.list(rep("", ncol(blank)))
    rows <- c(rows, list(blank, lst[[nm]]))
    hdr  <- c(hdr, k + 1); k <- k + 1 + nrow(lst[[nm]])
  }
  list(df = bind_rows(rows), rows = hdr, labels = names(lst))
}

# ---- 3. rendering: the table style of 03_main_results.R ---------------------------------------
# Same theme, rules, title and note wrapping as render() in 03 (copied, since sourcing 03 would rerun
# it). Two additions: panel label rows (bold, spanning the table, thin rule above), and the spanner
# row is only added when there are spanners.
render <- function(df, file, title, notes, spans = list(), panel_rows = integer(0),
                   panel_labels = character(0), left_cols = 1, note_width = 120) {
  th <- ttheme_minimal(
    base_size = 10,
    core    = list(fg_params = list(fontsize = 10, lineheight = 1.05),
                   padding = unit(c(5, 3.2), "mm")),
    colhead = list(fg_params = list(fontsize = 10, fontface = "bold", lineheight = 1.05),
                   padding = unit(c(5, 3), "mm")))
  g <- tableGrob(df, rows = NULL, theme = th)
  # left-align the label column(s)
  for (i in which(g$layout$name %in% c("core-fg", "colhead-fg") & g$layout$l %in% left_cols)) {
    g$grobs[[i]]$x <- unit(0.03, "npc"); g$grobs[[i]]$hjust <- 0
  }
  nc <- ncol(g)
  # panel labels: the data frame row is blank; draw the label across the full width instead
  for (j in seq_along(panel_rows)) {
    rr <- panel_rows[j] + 1                                          # +1 for the column-head row
    g$heights[rr] <- g$heights[rr + 1] + unit(2.5, "mm")
    g <- gtable_add_grob(g, textGrob(panel_labels[j], x = 0.03 * g$widths[1], hjust = 0,
                         y = unit(0.42, "npc"), gp = gpar(fontface = "bold", fontsize = 10)),
                         t = rr, l = 1, r = nc, clip = "off", name = paste0("panel-", j))
    if (rr > 2) g <- gtable_add_grob(g, segmentsGrob(y0 = unit(1, "npc"), y1 = unit(1, "npc"),
                                     gp = gpar(lwd = 0.6, col = "grey45")), t = rr, l = 1, r = nc)
  }
  hr <- 1                                                            # row of the column heads
  if (length(spans)) {
    g <- gtable_add_rows(g, heights = unit(7, "mm"), pos = 0); hr <- 2
    for (sp in spans) {
      g <- gtable_add_grob(g, textGrob(sp$label, gp = gpar(fontface = "bold", fontsize = 10)),
                           t = 1, l = sp$l, r = sp$r)
      g <- gtable_add_grob(g, segmentsGrob(x0 = unit(0.04, "npc"), x1 = unit(0.96, "npc"),
                           y0 = unit(0, "npc"), y1 = unit(0, "npc"), gp = gpar(lwd = 0.8)),
                           t = 1, l = sp$l, r = sp$r)
    }
  }
  rule <- function(y) segmentsGrob(y0 = unit(y, "npc"), y1 = unit(y, "npc"), gp = gpar(lwd = 1.2))
  g <- gtable_add_grob(g, rule(1), t = 1,  l = 1, r = nc)            # top rule
  g <- gtable_add_grob(g, rule(0), t = hr, l = 1, r = nc)            # under column heads
  g <- gtable_add_grob(g, rule(0), t = nrow(g), l = 1, r = nc)       # bottom rule
  # title and notes
  # title and notes inside the image only if IMAGE_TEXT (00_paths.R); otherwise saved as text
  if (IMAGE_TEXT) {
    g <- gtable_add_rows(g, heights = unit(10, "mm"), pos = 0)
    g <- gtable_add_grob(g, textGrob(title, x = unit(0, "npc"), hjust = 0,
                         gp = gpar(fontface = "bold", fontsize = 11.5)), t = 1, l = 1, r = nc)
    # wrap the notes to the table's own width (about 15.5 characters per inch at 8.5 pt)
    pdf(NULL); tw <- convertWidth(sum(g$widths), "in", valueOnly = TRUE); invisible(dev.off())
    nt <- paste(strwrap(notes, width = min(note_width, floor(tw * 15.5))), collapse = "\n")
    nlines <- length(strsplit(nt, "\n")[[1]])
    g <- gtable_add_rows(g, heights = unit(4 + 4.2 * nlines, "mm"))
    g <- gtable_add_grob(g, textGrob(nt, x = unit(0, "npc"), y = unit(1, "npc") - unit(2.5, "mm"),
                         hjust = 0, vjust = 1, gp = gpar(fontsize = 8.5, col = "grey20", lineheight = 1.1)),
                         t = nrow(g), l = 1, r = nc)
  } else save_note(file, title, notes)
  # size the device to the table
  pdf(NULL)
  wd <- convertWidth(sum(g$widths), "in", valueOnly = TRUE) + 0.5
  ht <- convertHeight(sum(g$heights), "in", valueOnly = TRUE) + 0.4
  invisible(dev.off())
  ragg::agg_png(file, width = wd, height = ht, units = "in", res = 300, background = "white")
  grid.newpage(); pushViewport(viewport(width = unit(wd - 0.5, "in"), height = unit(ht - 0.4, "in")))
  grid.draw(g); invisible(dev.off())
  cat(sprintf("written %s (%.1f x %.1f in)\n", file, wd, ht))
}

# shared note pieces (definitions, no dashes as punctuation)
note_tests <- paste(
  "Continuous rows show means and binary rows show percentages; Difference is the first group minus",
  "the second (percentage points for binary rows). p values come from Welch t-tests for continuous",
  "variables and chi-squared tests with continuity correction for binary ones. n is the number of people",
  "with data on each variable, first group / second group.")
note_defs <- paste(
  "Age is measured at the baseline interview, from birth dates. Education is coded 1 to 10 (1 none, 3 primary,",
  "4 lower secondary, 9 university degree, 10 postgraduate). Cared for by a son or daughter refers to the",
  "main caregiver, including caregivers recorded as \"Hija/o\" (son or daughter).",
  "WHOQOL-BREF is the mean of the standard 26 items; it and the environmental domain are on the WHO 0 to 100",
  "scale (25 times the item mean minus 1).")

# ---- 4. invited versus never invited: all 120 at baseline and the estimation samples ----------
# Both comparisons are computed before either table is drawn, because the Table 2 note quotes the
# size of the indigenous gap in the estimation samples.
t2 <- balance(p, "z")
smp <- list(
  list(f = "in_mid",  lab = "Panel A. Midline sample (baseline and midline observed)",  key = "midline"),
  list(f = "in_end",  lab = "Panel B. Endline sample (baseline and endline observed)",  key = "endline"),
  list(f = "in_all3", lab = "Panel C. Observed at all three waves",                    key = "all3"))
# (sample sizes are computed outside mutate(): inside it, "p" would mean the p value column)
bs <- map_dfr(smp, function(s) { ns <- sum(p[[s$f]])
  balance(p %>% filter(.data[[s$f]]), "z") %>% mutate(panel = s$lab, n_sample = ns, .before = 1) })
short_panel <- function(b) b %>% mutate(panel = sub("^(Panel [A-C])\\..*", "\\1", panel))  # "Panel A"

# Guard 5: the notes say in words that indigenous status is the largest imbalance (smallest p value)
# in all 120 and in every estimation sample, and that it is the only binary row where the choice
# between chi-squared and Fisher matters. Stop if the data ever say otherwise.
smallest <- function(b) b$variable[which.min(b$p)]
stopifnot(smallest(t2) == "Indigenous (%)",
          all(map_chr(split(bs, bs$panel), smallest) == "Indigenous (%)"),
          all(c(disagree(t2)$variable, disagree(bs)$variable) %in% "Indigenous (%)"))
ind   <- t2 %>% filter(variable == "Indigenous (%)")
ind_s <- bs %>% filter(variable == "Indigenous (%)")      # panels A, B, C in order
stopifnot(ind$diff < 0, all(ind_s$diff < 0))              # the notes say "lower among the invited"
see_cov <- "Columns (c) and (d) of Table A.8 control for it."

# ---- 4a. MAIN TEXT Table 2 -------------------------------------------------------------------
ind_txt <- paste(
  sprintf("The largest imbalance is indigenous status (%.1f per cent of the invited, %.1f per cent of the never",
          ind$mean1, ind$mean0),
  sprintf("invited; %.1f percentage points, chi-squared p = %s, Fisher exact p = %s). The gap is %.0f to %.0f points",
          abs(ind$diff), fmt_p(ind$p), fmt_p(ind$p_fisher), min(abs(ind_s$diff)), max(abs(ind_s$diff))),
  "in the estimation samples (Table A.5).", see_cov)   # the estimation samples table is A.5; A.4 is the take up table
render(fmt_balance(t2, "Invited", "Never invited"),
  file  = paste0(OUT_TAB, "table2_balance_baseline.png"),
  title = "Baseline characteristics by randomised invitation",
  notes = paste(
    sprintf("Sample: all %d people interviewed at baseline; %d were invited in one of the eight randomisation",
            nrow(p), sum(p$z)),
    sprintf("rounds and %d were never invited.", sum(p$z == 0)), note_tests, miss_text(p), note_defs,
    ind_txt, bonf_text(t2)), note_width = 118)

# Timing check (referee, October 2026). One invited person's baseline record is dated after
# randomisation round 1 (date not yet confirmed with IMSS). Recompute Table 2
# without anyone whose baseline is dated after their invitation; csv and console only.
n_late  <- sum(p$late_base)
t2_late <- balance(p %>% filter(!late_base), "z")
late_dp <- max(abs(c(t2_late$p - t2$p, t2_late$p_fisher - t2$p_fisher)), na.rm = TRUE)
write_csv(bind_rows(
    t2      %>% mutate(sample = "All interviewed at baseline", .before = 1),
    t2_late %>% mutate(sample = sprintf("Timing check: excluding %d invited with baseline dated after round 1", n_late),
                       .before = 1)) %>%
    mutate(across(where(is.double), ~ round(.x, 4))),
  paste0(OUT_CSV, "table2_balance_baseline.csv"))

# ---- 5. ANNEX: the same comparison inside the estimation samples ------------------------------
pan <- stack_panels(setNames(
  map(smp, function(s) fmt_balance(bs %>% filter(panel == s$lab), "Invited", "Never invited")),
  map_chr(smp, function(s) sprintf("%s: %d people (%d invited, %d never invited)", s$lab, sum(p[[s$f]]),
                                   sum(p[[s$f]] & p$z == 1), sum(p[[s$f]] & p$z == 0)))))
f3 <- function(x, f) and_list(f(x))                                     # "a, b and c"
render(pan$df, panel_rows = pan$rows, panel_labels = pan$labels,
  file  = paste0(OUT_TAB, "tableA_balance_samples.png"),
  title = "Balance by randomised invitation in the estimation samples and in the three wave panel",
  notes = paste(
    "Same comparison as Table 2, invited versus never invited, within the samples used for estimation.",
    "Panel A is the midline ANCOVA sample and Panel B the endline ANCOVA sample of the main results table;",
    "Panel C is the sample of the trajectories in Figure 4 (everyone observed at all three waves).",
    sprintf("Indigenous status is the largest imbalance in every panel: the share indigenous is %s percentage points",
            f3(abs(ind_s$diff), function(x) sprintf("%.1f", x))),
    sprintf("lower among the invited (chi-squared p = %s; Fisher exact p = %s).",
            f3(ind_s$p, fmt_p), f3(ind_s$p_fisher, fmt_p)), see_cov,
    note_tests, miss_text(p %>% filter(in_mid | in_end)), note_defs,
    bonf_text(short_panel(bs), per_panel = TRUE)),
  note_width = 118)
write_csv(bs %>% mutate(across(where(is.double), ~ round(.x, 4))), paste0(OUT_CSV, "tableA_balance_samples.csv"))

# ---- 6. ANNEX: attended versus did not attend, among the invited ------------------------------
# Panel A: every invited person, attendance over the whole register (point (c) of the header).
# Panels B and C: the invited people of each ANCOVA sample, attendance up to that wave's survey,
# which is the attendance instrumented in the main results.
inv <- p %>% filter(z == 1)
smp_a <- list(
  # (labels kept short so the rendered panel rows fit the table width; the note gives the definitions)
  list(df = inv %>% mutate(att = ever),
       lab = "Panel A. All invited people, whole register",   key = "all_invited"),
  list(df = inv %>% filter(in_mid) %>% mutate(att = d1),
       lab = "Panel B. Invited people in the midline sample", key = "midline"),
  list(df = inv %>% filter(in_end) %>% mutate(att = d2),
       lab = "Panel C. Invited people in the endline sample", key = "endline"))
ba <- map_dfr(smp_a, function(s) balance(s$df, "att") %>%
                mutate(panel = s$lab, n_sample = nrow(s$df), .before = 1))
pan <- stack_panels(setNames(
  map(smp_a, function(s) fmt_balance(ba %>% filter(panel == s$lab), "Attended", "Did not attend")),
  map_chr(smp_a, function(s) sprintf("%s: %d people (%d attended, %d did not)", s$lab, nrow(s$df),
                                     sum(s$df$att == 1), sum(s$df$att == 0)))))
render(pan$df, panel_rows = pan$rows, panel_labels = pan$labels,
  file  = paste0(OUT_TAB, "tableA_attended.png"),
  title = "Who took up the offer: invited people who attended and who did not",
  notes = paste(
    "Sample: people invited in one of the eight randomisation rounds. Attended means more than five days at",
    sprintf("the centre (more than one week). Panel A counts days over the centre's whole attendance register (%s),", reg_period),
    "so it classifies every invited person, including those not surveyed after baseline. Panels B and C count",
    "days up to that wave's survey, the attendance measure of the main results, among the invited people of",
    "each ANCOVA sample.", ever_txt, "Nobody who was not invited attended. The comparison describes who took",
    "up the offer (compliers) against invited people who did not (never takers).", note_tests,
    miss_text(inv), fisher_text(short_panel(ba)), note_defs, bonf_text(short_panel(ba), per_panel = TRUE)),
  note_width = 118)
write_csv(ba %>% mutate(across(where(is.double), ~ round(.x, 4))), paste0(OUT_CSV, "tableA_attended.csv"))

# ---- 7. ANNEX: attrition ------------------------------------------------------------------------
# 7a. retained versus not retained, baseline characteristics. "Not retained" rather than "lost":
# 3 of the 38 people not observed at all three waves missed only the midline interview and are in
# the endline sample.
smp_r <- list(
  list(f = "in_mid",  lab = "Panel A. Midline sample",              key = "midline"),
  list(f = "in_end",  lab = "Panel B. Endline sample",              key = "endline"),
  list(f = "in_all3", lab = "Panel C. Observed at all three waves", key = "all3"))
br <- map_dfr(smp_r, function(s) { ns <- nrow(p)
  balance(p %>% mutate(g = as.integer(.data[[s$f]])), "g") %>% mutate(panel = s$lab, n_sample = ns, .before = 1) })
n_endonly <- sum(p$in_end & !p$in_all3)                  # in the endline sample, missed midline
n_midonly <- sum(p$in_mid & !p$in_all3)                  # in the midline sample, missed endline
stopifnot(all(is.na(p$date1[p$in_end & !p$in_all3])))    # the endline-only people were not surveyed at midline

# 7b. attrition rates by take up group and by arm. "Not observed" = not surveyed at that wave, which
# is the same as being outside that wave's ANCOVA sample (guard 3).
grp_def <- list(                                         # label, people
  "Attended the centre"     = p$grp == "attended",
  "Invited, did not attend" = p$grp == "invited_no",
  "Never invited"           = p$grp == "never",
  "All invited"             = p$z == 1,
  "Everyone"                = rep(TRUE, nrow(p)))
lost_def <- list(Midline = !p$in_mid, Endline = !p$in_end)
# (counts computed before tibble(): inside it, "lost" would mean the new count column)
rates <- imap_dfr(lost_def, function(is_lost, w) imap_dfr(grp_def, function(g, gl) {
  n_g <- sum(g); n_lost <- sum(is_lost & g)
  tibble(wave = w, group = gl, n = n_g, lost = n_lost, pct = 100 * n_lost / n_g) }))
stopifnot(all(rates$n[rates$group == "Everyone"] == nrow(p)),
          sum(rates$n[rates$group %in% names(grp_def)[1:3] & rates$wave == "Midline"]) == nrow(p),
          all(rates$lost[rates$group == "Everyone"] == c(sum(!p$in_mid), sum(!p$in_end))),
          all(rates$pct < 100))
# tests of equal attrition: Fisher exact (and, for the csv, the uncorrected chi-squared)
test_def <- list(
  "Three groups"                                = c("attended", "invited_no", "never"),
  "Attended versus invited, did not attend"     = c("attended", "invited_no"),
  "Attended versus never invited"               = c("attended", "never"),
  "Invited, did not attend versus never invited" = c("invited_no", "never"),
  "By arm: invited versus never invited"        = "arm")
tests <- imap_dfr(lost_def, function(lost, w) imap_dfr(test_def, function(gs, tl) {
  if (identical(gs, "arm")) { keep <- rep(TRUE, nrow(p)); g <- p$z } else { keep <- p$grp %in% gs; g <- p$grp }
  tb <- table(factor(g[keep]), factor(lost[keep], levels = c(FALSE, TRUE)))
  tibble(wave = w, test = tl, p_fisher = fisher.test(tb)$p.value,
         p_chisq_uncorrected = suppressWarnings(chisq.test(tb, correct = FALSE))$p.value)
}))

# rendered rates table: Panel A counts and per cents, Panel B Fisher exact p values
rate_cell <- function(gl, w) with(rates %>% filter(group == gl, wave == w), sprintf("%d (%.1f%%)", lost, pct))
rate_tab <- tibble(Group = names(grp_def),
                   `People at\nbaseline` = map_chr(names(grp_def), ~ as.character(sum(grp_def[[.x]]))),
                   `Not surveyed\nat midline` = map_chr(names(grp_def), rate_cell, w = "Midline"),
                   `Not surveyed\nat endline` = map_chr(names(grp_def), rate_cell, w = "Endline"))
test_tab <- tibble(Group = names(test_def), `People at\nbaseline` = "",
                   `Not surveyed\nat midline` = fmt_p(tests$p_fisher[tests$wave == "Midline"]),
                   `Not surveyed\nat endline` = fmt_p(tests$p_fisher[tests$wave == "Endline"]))
pan_rt <- stack_panels(list(
  "Panel A. People not surveyed, number (per cent of the group)" = rate_tab,
  "Panel B. Tests of equal attrition, Fisher exact p values"     = test_tab))

# sentence on attrition by arm, quoted in the note of the retained versus not retained table (7c)
a_txt <- function(w) {
  r <- function(gl) rates %>% filter(group == gl, wave == w)
  ri <- r("All invited"); rn <- r("Never invited")
  sprintf("%d of %d invited (%.1f per cent) and %d of %d never invited (%.1f per cent; Fisher exact p = %.3f)",
          ri$lost, ri$n, ri$pct, rn$lost, rn$n, rn$pct,
          tests$p_fisher[tests$wave == w & tests$test == "By arm: invited versus never invited"])
}
arm_txt <- paste0("Attrition by arm: ", a_txt("Midline"), " were not surveyed at midline, and ", a_txt("Endline"),
                  " at endline.")

render(pan_rt$df, panel_rows = pan_rt$rows, panel_labels = pan_rt$labels,
  file  = paste0(OUT_TAB, "tableA_attrition_rates.png"),
  title = "Attrition by take up and by arm",
  notes = paste(
    sprintf("Sample: all %d people interviewed at baseline. Not surveyed means no interview at that wave, which is",
            nrow(p)),
    "the same as being outside that wave's ANCOVA sample. Attended the centre means more than five days in",
    sprintf("the centre's attendance register over its whole period (%s), so that every person has a group,", reg_period),
    "whether surveyed or not. The main results instead count days up to each wave's survey, which exists only",
    "for people surveyed at that wave.", ever_txt, "Nobody who was not invited attended.",
    "Panel B tests whether the share not surveyed is the same in the groups named in each row (Fisher exact test",
    "on the groups by surveyed table). The take up groups are defined after randomisation: to count as",
    "attending, a person had to stay alive and available, so differences between take up groups describe who",
    "stayed rather than test for attrition bias. Only the comparison by arm (invited versus never invited)",
    "tests differential attrition between the randomised groups."),
  note_width = 110)
write_csv(rates %>% mutate(pct = round(pct, 4)), paste0(OUT_CSV, "tableA_attrition_rates.csv"))
write_csv(tests %>% mutate(across(where(is.double), ~ round(.x, 4))), paste0(OUT_CSV, "tableA_attrition_tests.csv"))

# 7c. rendered balance table, retained versus not retained
pan <- stack_panels(setNames(
  map(smp_r, function(s) fmt_balance(br %>% filter(panel == s$lab), "Retained", "Not retained")),
  map_chr(smp_r, function(s) sprintf("%s: %d retained, %d not retained", s$lab, sum(p[[s$f]]), sum(!p[[s$f]])))))
render(pan$df, panel_rows = pan$rows, panel_labels = pan$labels,
  file  = paste0(OUT_TAB, "tableA_attrition.png"),
  title = "Attrition: baseline characteristics of people retained and not retained",
  notes = paste(
    sprintf("Sample: all %d people interviewed at baseline. Retained means observed at baseline and midline", nrow(p)),
    "(Panel A, the midline ANCOVA sample), at baseline and endline (Panel B, the endline ANCOVA sample) or",
    sprintf("at all three waves (Panel C). %s of the %d not retained in Panel C were observed at endline and",
            cap(num_word(n_endonly)), sum(!p$in_all3)),
    sprintf("are in the endline sample; they missed only the midline interview. %s", arm_txt),
    "Rates by take up are in Table A.7.",
    note_tests, miss_text(p), fisher_text(short_panel(br)), note_defs, bonf_text(short_panel(br), per_panel = TRUE)),
  note_width = 118)
write_csv(br %>% mutate(across(where(is.double), ~ round(.x, 4))), paste0(OUT_CSV, "tableA_attrition.csv"))

# ---- 8. Background sentence: who attended -----------------------------------------------------
# The submitted text: 81 years, 70.8 per cent women, 58.3 per cent primary education, 69.6 per cent
# cared for by children, 21.7 per cent by partners. The text now quotes the attenders in the endline
# sample (attendance by the endline survey), the "Attended" column of Panel C of tableA_attended.
# The same figures for everyone who attended over the whole register are kept for reference.
att_end <- p %>% filter(in_end, d2 == 1)
att_reg <- p %>% filter(ever == 1)
bg_one <- function(df, lab) {
  pc <- function(x) 100 * mean(x, na.rm = TRUE)
  nn <- function(x) sum(!is.na(x))
  tribble(
    ~measure,                                          ~value,                      ~denominator,             ~submitted, ~definition,
    "People",                                          nrow(df),                    nrow(df),                 NA,   "attended the centre for more than five days",
    "Mean age at baseline (years)",                    mean(df$age0),               nn(df$age0),              81,   "from birth dates at the baseline interview; the submission computed age at the script run date",
    "Women (%)",                                       pc(df$female),               nn(df$female),            70.8, "sex recorded at baseline",
    "Cared for by a son or daughter (%)",              pc(df$child_carer),          nn(df$child_carer),       69.6, "main caregiver answered Padre o madre or Hija/o; people without a caregiver record excluded",
    "Cared for by a spouse or partner (%)",            pc(df$partner),              nn(df$partner),           21.7, "main caregiver answered Esposa/o o companera/o; people without a caregiver record excluded",
    "Cared for by a son or daughter, control.csv coding (%)", pc(df$child_submitted), nrow(df),             NA,   "Padre o madre only, missing counted as no: the coding of the submitted balance tables",
    "Primary education as highest level (%)",          pc(df$educ == 3),            nn(df$educ),              58.3, "education code 3 (Primaria) exactly",
    "Primary education or less (%)",                   pc(df$educ <= 3),            nn(df$educ),              NA,   "education codes 1 to 3 (none, preschool, primary)",
    "Education, mean on the 1 to 10 scale",            mean(df$educ, na.rm = TRUE), nn(df$educ),              NA,   "1 none, 3 primary, 4 lower secondary, 9 university degree, 10 postgraduate",
    "WHOQOL-BREF at baseline (0 to 100)",              mean(df$whoqol_y0),          nn(df$whoqol_y0),         NA,   "WHO 0 to 100 scale") %>%
    mutate(sample = lab, .before = 1)
}
bg_raw <- bind_rows(bg_one(att_end, "Attended by the endline survey, endline sample (the text)"),
                    bg_one(att_reg, "Attended over the whole register, all invited (reference)"))

# Guard 6: the Background numbers are the "Attended" column of Panel C of tableA_attended, exactly
pc_att <- ba %>% filter(panel == smp_a[[3]]$lab)
bg_end <- bg_raw %>% filter(grepl("endline sample", sample))
same <- function(m, v) isTRUE(all.equal(bg_end$value[bg_end$measure == m], pc_att$mean1[pc_att$variable == v], tolerance = 1e-12))
stopifnot(bg_end$value[bg_end$measure == "People"] == pc_att$n1[pc_att$variable == "Age at baseline (years)"],
          same("Mean age at baseline (years)", "Age at baseline (years)"),
          same("Women (%)", "Female (%)"),
          same("Cared for by a son or daughter (%)", "Cared for by a son or daughter (%)"),
          same("WHOQOL-BREF at baseline (0 to 100)", "WHOQOL-BREF at baseline (0 to 100)"),
          setequal(att_end$id, inv$id[inv$in_end & inv$d2 == 1]))
bg <- bg_raw %>%
  mutate(value = round(value, 1),
         reproduced = case_when(is.na(submitted) ~ NA_character_,
                                abs(value - submitted) < 0.05 ~ "yes", TRUE ~ "no"))
write_csv(bg, paste0(OUT_CSV, "background_participants.csv"))

# ---- 9. numbers for the text and the response letter ------------------------------------------
tn <- function(section, item, value, detail) tibble(section = section, item = item, value = round(value, 4), detail = detail)
# long form of a balance table: means, difference and p values of every row
tn_bal <- function(section, b, prefix, g1, g0) {
  key <- vars$v[match(b$variable, vars$label)]
  bind_rows(
    tn(section, paste0(prefix, key, "_", g1), b$mean1, paste0(b$variable, ", n = ", b$n1)),
    tn(section, paste0(prefix, key, "_", g0), b$mean0, paste0(b$variable, ", n = ", b$n0)),
    tn(section, paste0(prefix, key, "_diff"), b$diff, paste(b$variable, "difference, first minus second group")),
    tn(section, paste0(prefix, key, "_p"), b$p, paste(b$variable, b$test, "p value as shown in the table")),
    tn(section, paste0(prefix, key, "_p_fisher"), b$p_fisher, paste(b$variable, "Fisher exact p value")) %>%
      filter(!is.na(value)))
}
# rows of a panelled balance table with p < 0.10 under either test, as the text may mention them
tn_notable <- function(section, b, keys, g1, g0) {
  map2_dfr(unique(b$panel), keys, function(pl, k) {
    bb <- b %>% filter(panel == pl, pmin(p, p_fisher, na.rm = TRUE) < 0.10)
    if (nrow(bb) == 0) return(tn(section, paste0(k, "_rows_p_below_010_n"), 0, pl))
    bind_rows(tn(section, paste0(k, "_rows_p_below_010_n"), nrow(bb), pl), tn_bal(section, bb, paste0(k, "_"), g1, g0))
  })
}

# Caregiver recode: how many caregivers it reclassifies, and what the row looked like in the 82
# people of the submitted Table 2 under each coding (the submitted coding counts a missing record as
# "not a child", as control.csv did).
cg <- function(df, v) {
  s <- df %>% filter(!is.na(.data[[v]]))
  tb <- table(factor(s[[v]], levels = 0:1), factor(s$z, levels = 0:1))
  c(inv = 100 * mean(s[[v]][s$z == 1]), not = 100 * mean(s[[v]][s$z == 0]),
    p = suppressWarnings(chisq.test(tb))$p.value)
}
cg_old82 <- cg(p %>% filter(in_all3), "child_submitted"); cg_new82 <- cg(p %>% filter(in_all3), "child_carer")

txt <- bind_rows(
  # samples
  tn("samples", "n_baseline", nrow(p), "everyone randomised, all interviewed at baseline"),
  tn("samples", "n_invited", sum(p$z), "invited in one of the eight rounds"),
  tn("samples", "n_never_invited", sum(p$z == 0), "never invited"),
  map_dfr(smp, function(s) bind_rows(
    tn("samples", paste0("n_", s$key, "_sample"), sum(p[[s$f]]), s$lab),
    tn("samples", paste0("n_", s$key, "_sample_invited"), sum(p[[s$f]] & p$z == 1), s$lab),
    tn("samples", paste0("n_", s$key, "_sample_never_invited"), sum(p[[s$f]] & p$z == 0), s$lab))),
  tn("samples", "n_endline_not_midline", n_endonly, "in the endline sample, not surveyed at midline"),
  tn("samples", "n_midline_not_endline", n_midonly, "in the midline sample, not surveyed at endline"),
  tn("samples", "n_attending_midline", sum(p$d1, na.rm = TRUE), "more than five days by the midline survey, midline sample"),
  tn("samples", "n_attending_endline", sum(p$d2, na.rm = TRUE), "more than five days by the endline survey, endline sample"),
  tn("samples", "n_attended_register", n_ever, sprintf("more than five days in the whole register (%s)", reg_period)),
  tn("samples", "n_attended_register_began_after_midline", n_ever_late, "of those: surveyed at midline with five days or fewer by then"),
  tn("samples", "n_attended_register_no_followup", n_ever_nofu, "of those: no survey after baseline"),
  tn("samples", "n_attended_register_not_at_endline", n_ever_noend, "of those: not surveyed at endline"),
  tn("samples", "n_invited_any_day_register", sum(p$days_reg > 0 & p$z == 1), "invited with at least one day in the register"),
  tn("samples", "n_missing_katz_lawton", sum(is.na(p$katz0)), "all invited; Katz and Lawton missing together"),
  tn("samples", "n_missing_caregiver", sum(is.na(p$child_carer)), "all invited"),
  # Table 2, every row
  tn_bal("table2", t2, "t2_", "invited", "not_invited"),
  # indigenous gap: Table 2 and the estimation samples
  tn("indigenous", "indigenous_invited_pct_120", ind$mean1, "all 120 at baseline"),
  tn("indigenous", "indigenous_not_invited_pct_120", ind$mean0, "all 120 at baseline"),
  tn("indigenous", "indigenous_gap_pp_120", ind$diff, "invited minus never invited, percentage points"),
  tn("indigenous", "indigenous_p_chisq_corrected_120", ind$p, "as shown in Table 2"),
  tn("indigenous", "indigenous_p_chisq_uncorrected_120", ind$p_chisq_uncorrected, "no continuity correction"),
  tn("indigenous", "indigenous_p_fisher_120", ind$p_fisher, "Fisher exact"),
  map_dfr(seq_along(smp), function(i) { k <- paste0(smp[[i]]$key, sum(p[[smp[[i]]$f]])); bind_rows(
    tn("indigenous", paste0("indigenous_invited_pct_", k), ind_s$mean1[i], ind_s$panel[i]),
    tn("indigenous", paste0("indigenous_not_invited_pct_", k), ind_s$mean0[i], ind_s$panel[i]),
    tn("indigenous", paste0("indigenous_gap_pp_", k), ind_s$diff[i], ind_s$panel[i]),
    tn("indigenous", paste0("indigenous_p_chisq_", k), ind_s$p[i], ind_s$panel[i]),
    tn("indigenous", paste0("indigenous_p_fisher_", k), ind_s$p_fisher[i], ind_s$panel[i])) }),
  # Bonferroni and smallest p values
  tn("bonferroni", "bonferroni_threshold", 0.05 / K, sprintf("%d comparisons per table or panel", K)),
  tn("bonferroni", "min_p_table2_any_test", min(pmin(t2$p, t2$p_fisher, na.rm = TRUE)), "Table 2, chi-squared or Fisher"),
  tn("bonferroni", "n_pass_table2", sum(pmin(t2$p, t2$p_fisher, na.rm = TRUE) < 0.05 / K), "rows of Table 2 below the threshold"),
  tn("bonferroni", "n_pass_balance_samples", sum(pmin(bs$p, bs$p_fisher, na.rm = TRUE) < 0.05 / K), "rows of the annex below the threshold"),
  tn("bonferroni", "n_pass_attended", sum(pmin(ba$p, ba$p_fisher, na.rm = TRUE) < 0.05 / K), "rows of the take up table below the threshold"),
  tn("bonferroni", "n_pass_attrition", sum(pmin(br$p, br$p_fisher, na.rm = TRUE) < 0.05 / K), "rows of the attrition table below the threshold"),
  # take up table: counts, and rows with p < 0.10
  map_dfr(smp_a, function(s) bind_rows(
    tn("attended", paste0("attended_", s$key, "_n_attended"), sum(s$df$att == 1), s$lab),
    tn("attended", paste0("attended_", s$key, "_n_not"), sum(s$df$att == 0), s$lab),
    tn("attended", paste0("attended_", s$key, "_min_p"),
       min(pmin(ba$p, ba$p_fisher, na.rm = TRUE)[ba$panel == s$lab]), paste(s$lab, "smallest p value, either test")))),
  tn_notable("attended", ba, paste0("attended_", map_chr(smp_a, "key")), "attended", "not_attended"),
  # attrition: rows with p < 0.10 in the retained versus not retained table
  map_dfr(smp_r, function(s) tn("attrition", paste0("retention_", s$key, "_min_p"),
                                min(pmin(br$p, br$p_fisher, na.rm = TRUE)[br$panel == s$lab]),
                                paste(s$lab, "smallest p value, either test"))),
  tn_notable("attrition", br, paste0("retention_", map_chr(smp_r, "key")), "retained", "not_retained"),
  # attrition rates by group and arm, and their tests
  rates %>% transmute(section = "attrition_rates",
                      item = paste0("attrition_", tolower(wave), "_", gsub("[^a-z]+", "_", tolower(group)), "_pct"),
                      value = round(pct, 4), detail = sprintf("%d of %d not surveyed", lost, n)),
  tests %>% transmute(section = "attrition_rates",
                      item = paste0("attrition_", tolower(wave), "_test_", gsub("[^a-z]+", "_", tolower(test)), "_p_fisher"),
                      value = round(p_fisher, 4), detail = paste(test, "Fisher exact")),
  # caregiver recode
  tn("caregiver", "caregiver_hija_recoded_n_120", sum(p$hija), "caregivers recorded as Hija/o, now counted as a child"),
  tn("caregiver", "caregiver_hija_recoded_n_82", sum(p$hija & p$in_all3), "the same among the 82 of the submitted balance tables"),
  tn("caregiver", "caregiver_child_invited_pct_82_submitted", cg_old82["inv"], "submitted coding, 82 people"),
  tn("caregiver", "caregiver_child_not_invited_pct_82_submitted", cg_old82["not"], "submitted coding, 82 people"),
  tn("caregiver", "caregiver_child_p_82_submitted", cg_old82["p"], "submitted coding, chi-squared"),
  tn("caregiver", "caregiver_child_invited_pct_82_corrected", cg_new82["inv"], "corrected coding, 82 people"),
  tn("caregiver", "caregiver_child_not_invited_pct_82_corrected", cg_new82["not"], "corrected coding, 82 people"),
  tn("caregiver", "caregiver_child_p_82_corrected", cg_new82["p"], "corrected coding, chi-squared"),
  # timing check
  tn("timing", "timing_late_baseline_n", n_late, "invited, baseline dated after round 1 (2 Jan 2023, to confirm)"),
  tn("timing", "timing_late_baseline_in_midline_sample_n", sum(p$late_base & p$in_mid), "of those"),
  tn("timing", "timing_late_baseline_in_endline_sample_n", sum(p$late_base & p$in_end), "of those"),
  tn("timing", "timing_table2_max_abs_change_p", late_dp, "Table 2 without them, chi-squared and Fisher"),
  tn("timing", "timing_table2_indigenous_p_fisher_without", t2_late$p_fisher[t2_late$variable == "Indigenous (%)"], "Table 2 without them"),
  # Background: attenders in the endline sample (unrounded; Panel C of tableA_attended)
  bg_end %>% transmute(section = "background",
                       item = paste0("background_", gsub("^_|_$", "", gsub("[^a-z0-9]+", "_", tolower(measure)))),
                       value = round(value, 4), detail = sprintf("n = %d; %s", denominator, definition)))
stopifnot(!anyDuplicated(txt$item))
write_csv(txt, paste0(OUT_CSV, "text_numbers_balance.csv"))

# ---- 10. console summary ------------------------------------------------------------------------
show <- function(b, title) {
  cat("\n===", title, "===\n")
  print(b %>% mutate(across(c(mean1, mean0, diff), ~ round(.x, 2)), across(c(p, p_fisher), ~ round(.x, 3))) %>%
          select(any_of("panel"), variable, n1, n0, mean1, mean0, diff, p, p_fisher) %>% as.data.frame(),
        row.names = FALSE)
}
show(t2, "Table 2: invited vs never invited, all at baseline")
cat(bonf_text(t2), "\n")
show(t2_late, sprintf("Timing check: Table 2 without %d invited with baseline after round 1", n_late))
show(bs, "Annex: invited vs never invited, estimation samples")
show(ba, "Annex: attended vs did not attend, among invited")
show(br, "Annex: retained vs not retained")
cat("\n=== attrition rates ===\n"); print(as.data.frame(rates %>% mutate(pct = round(pct, 1))), row.names = FALSE)
print(as.data.frame(tests %>% mutate(across(where(is.double), ~ round(.x, 3)))), row.names = FALSE)
cat("\n=== Background: who attended ===\n"); print(as.data.frame(bg %>% select(sample, measure, value, denominator, submitted, reproduced)), row.names = FALSE)
cat("\n=== text numbers ===\n"); print(as.data.frame(txt), row.names = FALSE)
cat("\nwritten to", OUT_TAB, "and", OUT_CSV, "\n")
