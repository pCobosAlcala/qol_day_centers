# 06_heterogeneity.R  (revision, October 2026; data core, 0 to 100 scale, fully interacted model)
#
# Heterogeneity of the ENDLINE effect by baseline characteristics. EXPLORATORY.
#   Reviewer 2 (R2.14) asks whether the effect differs by age; Pablo also wants sex and other
#   pre-treatment characteristics. The invitation is randomised and every characteristic is measured
#   at baseline, so splitting the sample by them is legitimate, but with 85 people (23 not invited)
#   each half has few controls and many tests are run: the table is a description, not a finding.
#
# Specification (the main table's, endline only, scores on the WHO 0 to 100 scale):
#   Within each subgroup g:  Y_i2 = a0 + a1*Y_i0 + a2*Z_i + e_i          (ANCOVA ITT, HC2 errors)
#                            IV: attendance by the endline survey (more than five days, d2)
#                            instrumented by Z, baseline value in both stages (as in 03)
#   Difference:              Y_i2 = (b0 + b1*Y_i0 + b2*Z_i) + G_i*(c0 + c1*Y_i0 + c2*Z_i) + e_i
#                            on the whole endline sample, FULLY interacted: intercept, baseline slope
#                            and invitation effect all differ by subgroup. c2 is then exactly the
#                            first subgroup's ITT minus the second's (the separate regressions are
#                            this model split in two), and its HC2 p value is the test shown. With
#                            HC2 the two halves' errors do not mix either, so the SE of c2 is
#                            sqrt(SE1^2 + SE0^2). Both facts are checked below. (The September version
#                            used a common baseline slope, so its difference did not match the gap
#                            between the two ITTs shown.) The same fully interacted IV gives the
#                            difference between the two subgroup IVs (csv only).
#   Magnitudes: ITT as per cent of each subgroup's OWN control group baseline mean first (Pablo's
#   preference), then the ITT in points and its SE, and its p value. No stars anywhere: a subgroup
#   clearing a threshold says nothing about heterogeneity, which is judged by the difference column.
#   Take-up: share of the invited in each subgroup who attended. Take-up differs across subgroups,
#   so a gap in the ITT can come from a gap in take-up rather than in the effect of attending; the
#   subgroup IV (the paper's primary estimate) is shown next to the ITT. The IV is shown whenever
#   someone invited in the subgroup attended; a robust first-stage F below 10 is flagged.
#   Multiple testing: 8 characteristics x 2 outcomes = 16 difference tests; the note reports how
#   many fall below p = 0.10, the number expected by chance (1.6), and the count after Holm.
#
# Subgroups (medians computed in the 85-person endline sample; the cut used is printed):
#   age at baseline (older = at or above the median, about 80.3), sex, indigenous,
#   education (last level completed, above the median = lower secondary or more),
#   people living with the older adult (item 05.05 "how many people live with you?", 0 = alone;
#   above the median = 3 or more), main caregiver is a son or daughter (corrected flag of 00_data.R,
#   child_carer; one missing), baseline Lawton autonomy (at or above the median; one missing),
#   baseline WHOQOL-BREF composite (at or above the median).
#   Ties: "at or above" for the near-continuous scores; "above" for the two counts where many
#   people sit exactly on the median, because "at or above" would leave only 12 or 20 people below.
#
# Data: build_analysis_data() in 00_data.R (the September version read pm_who.csv, pm.csv and
# control.csv on the 1 to 5 scale and rebuilt the caregiver flag itself; a copy of it and of its
# outputs is in the _previous subfolders, suffix _before_00data). The endline sample and every split
# are the same people as before, so ITTs, IVs, SEs and the fully interacted differences are exactly
# 25 times the old ones and p values do not change; per cents of the control baseline do change,
# because the zero of the scale moves (see 03, section 6).
#
# Outputs: 3_code/output/tables/tableA_heterogeneity.png and 3_code/output/csv/tableA_heterogeneity.csv
#   Aggregates only (subgroup counts, means, estimates); no individual rows leave this script.
# Run from the imss_paper folder:
#   "C:/Program Files/R/R-4.6.1/bin/Rscript.exe" 3_code/scripts/06_heterogeneity.R

rm(list = ls())
suppressPackageStartupMessages({
  library(tidyverse); library(AER); library(sandwich); library(lmtest)
  library(gridExtra); library(gtable); library(grid)
})
options(scipen = 999, width = 200)

source("3_code/scripts/00_paths.R")
source("3_code/scripts/00_data.R")

# special characters written as escapes so this file stays plain ASCII
NBSP <- "\u00a0"                                         # non-breaking space: keeps "p = 0.10" on one line
DAG  <- "\u2020"                                         # dagger: fewer than 8 not invited people
DDAG <- "\u2021"                                         # double dagger: weak first stage (F below 10)
nb   <- function(x) gsub(" ", NBSP, x, fixed = TRUE)     # glue a short phrase for the note wrapper

# ---- 1. data: one row per randomised person, scores on 0 to 100 ---------------------------
dat <- build_analysis_data()                             # checks against the old files run inside
outcomes <- OUTCOMES[c("whoqol", "environmental")]       # the two outcomes discussed in the text

# baseline characteristics (00_data.R, section 5.4): no dates, no names
chars <- dat %>%
  select(id, age0, female, indigenous, hh_size, educ, child_carer, lawton0, who0 = whoqol_y0)

# endline ANCOVA samples: everyone observed at baseline and endline; the same 85 for both outcomes
smp <- map(names(outcomes), ~ ancova_sample(dat, .x, 2)) %>% set_names(names(outcomes))
stopifnot(nrow(smp$whoqol) == 85, setequal(smp$whoqol$id, smp$environmental$id),
          identical(smp$whoqol$d, smp$environmental$d))
s85 <- chars %>% filter(id %in% smp$whoqol$id)

# ---- 2. subgroup indicators -------------------------------------------------------------
# Medians come from the endline estimation sample (85 people, the same for both outcomes).
med <- s85 %>% summarise(across(c(age0, educ, hh_size, lawton0, who0), ~ median(.x, na.rm = TRUE)))
cat(sprintf(paste0("\nmedians in the endline sample (n = %d): age %.2f | education code %g | people living ",
                   "with the older adult %g | Lawton %g | WHOQOL-BREF %.2f (0 to 100)\n"),
            nrow(s85), med$age0, med$educ, med$hh_size, med$lawton0, med$who0))
cat(sprintf("people living with the older adult in the 85: %d live alone, %d with 1 or 2, %d with 3 or more\n",
            sum(s85$hh_size == 0), sum(s85$hh_size %in% 1:2), sum(s85$hh_size >= 3)))
cat(sprintf("education in the 85: %d primary or less (codes 1 to 3), %d lower secondary or more (4 to 10)\n",
            sum(s85$educ <= 3), sum(s85$educ >= 4)))
cat(sprintf("main caregiver in the 85: %d a son or daughter, %d someone else, %d missing; Lawton missing: %d\n",
            sum(s85$child_carer == 1, na.rm = TRUE), sum(s85$child_carer == 0, na.rm = TRUE),
            sum(is.na(s85$child_carer)), sum(is.na(s85$lawton0))))

# The education and household labels below name the cut in words, so they hold only for these
# medians; stop if the data ever move them.
if (med$educ != 3) stop("median education is no longer code 3 (primary): rewrite the education labels")
if (med$hh_size != 2) stop("median household count is no longer 2: rewrite the household labels")
if (!any(s85$hh_size == 0)) stop("nobody lives alone any more: rewrite the household label")

# G = 1 is the subgroup listed first in the table; the difference is always G = 1 minus G = 0.
chars <- chars %>% mutate(
  g_age    = as.integer(age0 >= med$age0),               # older half
  g_sex    = as.integer(female == 1),                    # women
  g_indig  = as.integer(indigenous == 1),                # considers herself or himself indigenous
  g_edu    = as.integer(educ > med$educ),                # completed more than primary
  g_hh     = as.integer(hh_size > med$hh_size),          # lives with 3 or more people
  g_child  = child_carer,                                # main caregiver is a son or daughter (NA if missing)
  g_lawton = as.integer(lawton0 >= med$lawton0),         # more autonomous half (one NA)
  g_who    = as.integer(who0 >= med$who0))               # higher baseline quality of life

# Education codes of 00_data.R (EDUC_CODES): 1 none, 2 preschool, 3 primary, 4 lower secondary
# (secundaria), 5 technical after lower secondary, 6 upper secondary (preparatoria), 7 teacher
# training, 8 technical after upper secondary, 9 university, 10 postgraduate. With the cut at 3,
# every code above it is beyond primary, so the order among 4 to 10 does not matter here.
spec <- tribble(
  ~g,         ~title,                               ~lab1,                                                  ~lab0,
  "g_age",    "Age at baseline",                    sprintf("%.1f years or older", med$age0),               sprintf("Younger than %.1f", med$age0),
  "g_sex",    "Sex",                                "Women",                                                "Men",
  "g_indig",  "Indigenous",                         "Indigenous",                                           "Not indigenous",
  "g_edu",    "Education (last level completed)",   "Lower secondary or more",                              "Primary or less",
  "g_hh",     "People living with the older adult", "3 or more",                                            "Alone, 1 or 2",
  "g_child",  "Main caregiver",                     "A son or daughter",                                    "Someone else",
  "g_lawton", "Lawton autonomy at baseline",        sprintf("Score %g to 8 (more autonomous)", med$lawton0), sprintf("Score 0 to %g", med$lawton0 - 1),
  "g_who",    "WHOQOL-BREF at baseline",            sprintf("%.1f points or more", med$who0),               sprintf("Below %.1f points", med$who0))

# ---- 3. estimation helpers (from 03, with one guard) ------------------------------------
# HC2 as in 03. Guard: in a very small cell a point with leverage 1 makes HC2 return NaN rather
# than an error, so non-finite HC2 also falls back to HC1, with a warning; the type is recorded.
rob <- function(m) {
  V <- tryCatch(vcovHC(m, type = "HC2"), error = function(e) NULL)
  if (is.null(V) || any(!is.finite(V))) {
    warning("HC2 not finite, HC1 used"); V <- vcovHC(m, type = "HC1"); attr(V, "type") <- "HC1"
  } else attr(V, "type") <- "HC2"
  V
}
pull_est <- function(m, term) {
  V  <- rob(m); ct <- coeftest(m, vcov. = V); ci <- coefci(m, vcov. = V, level = 0.95)
  tibble(est = ct[term, 1], se = ct[term, 2], p = ct[term, 4], lo = ci[term, 1], hi = ci[term, 2],
         vcov = attr(V, "type"))
}
na_est <- tibble(est = NA_real_, se = NA_real_, p = NA_real_, lo = NA_real_, hi = NA_real_, vcov = NA_character_)

# Everything reported for one group of people (a subgroup, or the whole sample).
#   sd0: SD of the baseline outcome in the WHOLE endline sample, a common yardstick for the
#   SD version of the ITT (csv only)
cell_est <- function(sub, sd0) {
  itt <- pull_est(lm(y ~ y0 + z, sub), "z")
  fs  <- pull_est(lm(d ~ y0 + z, sub), "z")              # first stage within the group
  tu  <- mean(sub$d[sub$z == 1])                         # take-up: share of the invited who attended
  ok  <- tu > 0                                          # IV exists only if some invited person attended
  iv  <- if (ok) pull_est(ivreg(y ~ d + y0 | z + y0, data = sub), "d") else na_est
  cb  <- mean(sub$y0[sub$z == 0])                        # this group's control baseline mean
  tibble(n = nrow(sub), n_invited = sum(sub$z == 1), n_not_invited = sum(sub$z == 0),
         n_attended = sum(sub$d == 1), takeup = tu, control_base = cb,
         itt_pct = 100 * itt$est / cb, est = itt$est, se = itt$se, p = itt$p, lo = itt$lo, hi = itt$hi,
         itt_sd = itt$est / sd0, vcov = itt$vcov,
         fs_est = fs$est, fs_se = fs$se, fs_F = (fs$est / fs$se)^2, iv_estimable = ok,
         iv_pct = 100 * iv$est / cb, iv_est = iv$est, iv_se = iv$se, iv_p = iv$p,
         iv_lo = iv$lo, iv_hi = iv$hi, iv_vcov = iv$vcov)
}
FEW  <- 8                                                # flag: fewer than 8 not invited people
WEAK <- 10                                               # flag: robust first-stage F below 10

# ---- 4. estimation ----------------------------------------------------------------------
res <- list()
for (v in names(outcomes)) {
  s   <- smp[[v]] %>% left_join(chars %>% select(id, starts_with("g_")), by = "id")
  sd0 <- sd(s$y0)

  # whole endline sample: must reproduce the main table (guard below)
  res[[paste(v, "all")]] <- cell_est(s, sd0) %>%
    mutate(outcome = unname(outcomes[v]), g = "all", characteristic = "All",
           subgroup = "Whole endline sample", group = NA_integer_)

  for (k in seq_len(nrow(spec))) {
    gv <- spec$g[k]
    s2 <- s %>% filter(!is.na(.data[[gv]])) %>% mutate(G = .data[[gv]])
    # fully interacted ITT: the coefficient on z:G is the G = 1 ITT minus the G = 0 ITT
    dif <- pull_est(lm(y ~ (y0 + z) * G, s2), "z:G")
    # fully interacted IV (csv only): the coefficient on d:G is the gap between the subgroup IVs
    ivd <- tryCatch(pull_est(ivreg(y ~ (d + y0) * G | (z + y0) * G, data = s2), "d:G"),
                    error = function(e) na_est)
    one  <- cell_est(s2 %>% filter(G == 1), sd0)
    zero <- cell_est(s2 %>% filter(G == 0), sd0)
    # guards: the difference shown is exactly the gap between the two subgroup estimates, and
    # under HC2 its SE combines the two subgroup SEs (the halves do not share errors)
    stopifnot(abs(dif$est - (one$est - zero$est)) < 1e-8,
              abs(dif$se  - sqrt(one$se^2 + zero$se^2)) < 1e-8)
    if (one$iv_estimable && zero$iv_estimable)
      stopifnot(abs(ivd$est - (one$iv_est - zero$iv_est)) < 1e-8)
    for (gg in c(1, 0)) {
      res[[paste(v, gv, gg)]] <- (if (gg == 1) one else zero) %>% mutate(
        outcome = unname(outcomes[v]), g = gv, characteristic = spec$title[k],
        subgroup = if (gg == 1) spec$lab1[k] else spec$lab0[k], group = gg,
        diff_est = dif$est, diff_se = dif$se, diff_p = dif$p, diff_lo = dif$lo, diff_hi = dif$hi,
        diff_vcov = dif$vcov,
        iv_diff_est = ivd$est, iv_diff_se = ivd$se, iv_diff_p = ivd$p)
    }
  }
}
res <- bind_rows(res) %>%
  mutate(few_controls = n_not_invited < FEW,
         weak_iv      = iv_estimable & fs_F < WEAK) %>%
  relocate(outcome, characteristic, subgroup, group, n, n_invited, n_not_invited, few_controls,
           n_attended, takeup, control_base, itt_pct, est, se, p, lo, hi, itt_sd,
           iv_estimable, fs_F, weak_iv, iv_pct, iv_est, iv_se, iv_p, iv_lo, iv_hi,
           diff_est, diff_se, diff_p, diff_lo, diff_hi, iv_diff_est, iv_diff_se, iv_diff_p)

# guard 1: the whole-sample rows reproduce 03's endline results (main_results.csv, magnitudes.csv)
chk <- res %>% filter(g == "all")                        # rows in the order whoqol, environmental
stopifnot(all(chk$n == 85), all(chk$n_not_invited == 23),
          isTRUE(all.equal(round(chk$est, 2),    c(4.68, 8.94))),
          isTRUE(all.equal(round(chk$se, 2),     c(2.57, 3.18))),
          isTRUE(all.equal(round(chk$p, 3),      c(0.072, 0.006))),
          isTRUE(all.equal(round(chk$iv_est, 2), c(11.81, 22.99))),
          isTRUE(all.equal(round(chk$control_base, 1), c(53.7, 55.6))),
          isTRUE(all.equal(round(chk$itt_pct, 1), c(8.7, 16.1))))
# guard 2: n and take-up are the same in both outcomes (same people), so the table shows them once
tw2 <- res %>% filter(g != "all") %>% select(outcome, g, group, n_invited, n_not_invited, takeup) %>%
  pivot_wider(names_from = outcome, values_from = c(n_invited, n_not_invited, takeup))
stopifnot(all(tw2$`n_invited_WHOQOL-BREF` == tw2$n_invited_Environmental),
          all(tw2$`n_not_invited_WHOQOL-BREF` == tw2$n_not_invited_Environmental),
          all(tw2$`takeup_WHOQOL-BREF` == tw2$takeup_Environmental))
if (any(c(res$vcov, res$iv_vcov, res$diff_vcov) == "HC1", na.rm = TRUE))
  cat("NOTE: HC2 was not finite in at least one cell; HC1 used there (see the vcov columns)\n")

# ---- 5. multiple testing ----------------------------------------------------------------
# 8 characteristics x 2 outcomes = 16 tests of a difference in the ITT. If no true difference
# exists, about 10 per cent of them clear p < 0.10 by chance. The tests are not independent (the
# composite contains the environmental items), so the count is a rough benchmark. Holm adjusted
# p values across the 16 tests control the chance of any false positive.
difs <- res %>% filter(g != "all", group == 1) %>%
  select(outcome, characteristic, diff_est, diff_p, iv_diff_p) %>%
  mutate(diff_p_holm = p.adjust(diff_p, "holm"))
n_tests  <- nrow(difs)
expected <- 0.10 * n_tests
observed <- sum(difs$diff_p < 0.10)
holm_obs <- sum(difs$diff_p_holm < 0.10)
min_holm <- min(difs$diff_p_holm)
res <- res %>% left_join(difs %>% select(outcome, characteristic, diff_p_holm),
                         by = c("outcome", "characteristic"))
hits <- difs %>% filter(diff_p < 0.10)
stopifnot(n_tests == 16)
cat(sprintf("\ndifference tests: %d | expected below p = 0.10 by chance: %.1f | observed: %d%s | below 0.10 after Holm: %d | smallest Holm p: %.3f | smallest raw p: %.3f\n",
            n_tests, expected, observed,
            if (observed > 0) paste0(" (", paste(hits$characteristic, hits$outcome, collapse = "; "), ")") else "",
            holm_obs, min_holm, min(difs$diff_p)))
cat(sprintf("IV differences below p = 0.10: %d of %d (csv only)\n",
            sum(difs$iv_diff_p < 0.10, na.rm = TRUE), sum(!is.na(difs$iv_diff_p))))
sub_rows <- res %>% filter(g != "all")
cat(sprintf("subgroup ITTs with p < 0.10: %d of %d (shown without stars: not a test of heterogeneity)\n",
            sum(sub_rows$p < 0.10), nrow(sub_rows)))

write_csv(res, paste0(OUT_CSV, "tableA_heterogeneity.csv"))

print(res %>% mutate(across(c(takeup, control_base, itt_pct, est, se, p, fs_F, iv_est, iv_se, iv_p,
                                diff_est, diff_se, diff_p, diff_p_holm, iv_diff_est, iv_diff_p),
                              ~ round(.x, 3))) %>%
        select(outcome, characteristic, subgroup, n_invited, n_not_invited, few_controls, takeup,
               control_base, itt_pct, est, se, p, fs_F, iv_est, iv_se, iv_p, diff_est, diff_se, diff_p,
               diff_p_holm, iv_diff_est, iv_diff_p) %>%
        as.data.frame(), row.names = FALSE)

# ---- 6. rendering (render() from 03, plus bold rows and a note wrapped by measured width) --
# Changes against 03's render():
#   bold_rows  the characteristic rows, which carry the difference
#   notes      a vector of paragraphs; each is wrapped by MEASURING the printed width of the text
#              in the note's font, so lines fill the table's width and are even (strwrap counts
#              characters, which in a proportional font leaves ragged lines well short of the
#              table edge). Non-breaking spaces glue short phrases such as "p = 0.10".
wrap_note <- function(paragraphs, width_in, fontsize) {
  wd <- function(s) convertWidth(grobWidth(textGrob(gsub(NBSP, " ", s, fixed = TRUE),
                                                    gp = gpar(fontsize = fontsize))), "in", valueOnly = TRUE)
  out <- character(0)
  for (par in paragraphs) {
    words <- strsplit(gsub(" +", " ", trimws(par)), " ", fixed = TRUE)[[1]]
    cur <- ""
    for (w in words) {
      cand <- if (cur == "") w else paste(cur, w)
      if (cur != "" && wd(cand) > width_in) { out <- c(out, cur); cur <- w } else cur <- cand
    }
    out <- c(out, cur)
  }
  gsub(NBSP, " ", out, fixed = TRUE)
}

# Measurements are taken on a throwaway ragg device, the same engine that draws the PNG, so the
# wrapped lines and the height reserved for the note match what is drawn (03 sizes the note at a
# fixed 4.2 mm per line instead).
measure_on <- function() {
  f <- tempfile(fileext = ".png")
  ragg::agg_png(f, width = 30, height = 30, units = "in", res = 72)
  f
}
measure_off <- function(f) { invisible(dev.off()); unlink(f) }

render <- function(df, spans, file, title, notes, left_cols = 1, bold_rows = NULL, note_size = 8.5) {
  # horizontal cell padding 3.5 mm (5 in 03): thirteen columns would otherwise make the table
  # about 14 inches wide
  th <- ttheme_minimal(
    base_size = 10,
    core    = list(fg_params = list(fontsize = 10, lineheight = 1.05),
                   padding = unit(c(3.5, 3.2), "mm")),
    colhead = list(fg_params = list(fontsize = 10, fontface = "bold", lineheight = 1.05),
                   padding = unit(c(3.5, 3), "mm")))
  g <- tableGrob(df, rows = NULL, theme = th)
  # left-align the label column(s)
  for (i in which(g$layout$name %in% c("core-fg", "colhead-fg") & g$layout$l %in% left_cols)) {
    g$grobs[[i]]$x <- unit(0.03, "npc"); g$grobs[[i]]$hjust <- 0
  }
  # bold body rows (data row k sits in gtable row k + 1, under the column heads)
  for (i in which(g$layout$name == "core-fg" & g$layout$t %in% (bold_rows + 1))) {
    g$grobs[[i]]$gp <- gpar(fontsize = 10, fontface = "bold", lineheight = 1.05)
  }
  nc <- ncol(g)
  # spanner row
  g <- gtable_add_rows(g, heights = unit(7, "mm"), pos = 0)
  for (sp in spans) {
    g <- gtable_add_grob(g, textGrob(sp$label, gp = gpar(fontface = "bold", fontsize = 10)),
                         t = 1, l = sp$l, r = sp$r)
    g <- gtable_add_grob(g, segmentsGrob(x0 = unit(0.04, "npc"), x1 = unit(0.96, "npc"),
                         y0 = unit(0, "npc"), y1 = unit(0, "npc"), gp = gpar(lwd = 0.8)),
                         t = 1, l = sp$l, r = sp$r)
  }
  rule <- function(y) segmentsGrob(y0 = unit(y, "npc"), y1 = unit(y, "npc"), gp = gpar(lwd = 1.2))
  g <- gtable_add_grob(g, rule(1), t = 1, l = 1, r = nc)             # top rule
  g <- gtable_add_grob(g, rule(0), t = 2, l = 1, r = nc)             # under column heads
  g <- gtable_add_grob(g, rule(0), t = nrow(g), l = 1, r = nc)       # bottom rule
  # title
  # title and notes inside the image only if IMAGE_TEXT (00_paths.R); otherwise saved as text
  if (IMAGE_TEXT) {
    g <- gtable_add_rows(g, heights = unit(10, "mm"), pos = 0)
    g <- gtable_add_grob(g, textGrob(title, x = unit(0, "npc"), hjust = 0,
                         gp = gpar(fontface = "bold", fontsize = 11.5)), t = 1, l = 1, r = nc)
    # notes: wrapped to the table's own width (2 per cent kept as a margin), paragraphs kept apart;
    # the row is as tall as the drawn text plus 2.5 mm above and 3 mm below
    mf <- measure_on()
    tw <- convertWidth(sum(g$widths), "in", valueOnly = TRUE)
    nt <- wrap_note(notes, width_in = 0.98 * tw, fontsize = note_size)
    note_grob <- textGrob(paste(nt, collapse = "\n"), x = unit(0, "npc"),
                          y = unit(1, "npc") - unit(2.5, "mm"), hjust = 0, vjust = 1,
                          gp = gpar(fontsize = note_size, col = "grey20", lineheight = 1.1))
    note_mm <- convertHeight(grobHeight(note_grob), "mm", valueOnly = TRUE)
    measure_off(mf)
    g <- gtable_add_rows(g, heights = unit(note_mm + 5.5, "mm"))
    g <- gtable_add_grob(g, note_grob, t = nrow(g), l = 1, r = nc)
  } else save_note(file, title, notes)
  # size the device to the table
  mf <- measure_on()
  wd <- convertWidth(sum(g$widths), "in", valueOnly = TRUE) + 0.5
  ht <- convertHeight(sum(g$heights), "in", valueOnly = TRUE) + 0.4
  measure_off(mf)
  ragg::agg_png(file, width = wd, height = ht, units = "in", res = 300, background = "white")
  grid.newpage(); pushViewport(viewport(width = unit(wd - 0.5, "in"), height = unit(ht - 0.4, "in")))
  grid.draw(g); invisible(dev.off())
  cat(sprintf("written %s (%.1f x %.1f in)\n", file, wd, ht))
}

# ---- 7. Table A: heterogeneity ----------------------------------------------------------
# One block per characteristic: a bold row carrying the difference in the ITT and its p value,
# then the two subgroups. Per outcome: per cent of the subgroup's control baseline, ITT in points
# (SE), its p value, IV in points (SE), difference. No stars.
pc    <- function(x) sprintf("%+.1f%%", x)
pts   <- function(e, s) sprintf("%.1f (%.1f)", e, s)
fmt_p <- function(p) ifelse(p < 0.1, sprintf("%.3f", p), sprintf("%.2f", p))   # three decimals near the thresholds
ivc   <- function(r) if (!r$iv_estimable) "n/e" else paste0(pts(r$iv_est, r$iv_se), if (r$weak_iv) paste0(" ", DDAG) else "")
half  <- function(r) c(pc(r$itt_pct), pts(r$est, r$se), fmt_p(r$p), ivc(r), "")
dcell <- function(r) sprintf("%+.1f, p = %s\n[%.1f to %.1f]", r$diff_est, fmt_p(r$diff_p), r$diff_lo, r$diff_hi)

cols <- c("Subgroup", "Invited /\nnever invited", "Invited who\nattended",
          "Per cent of\ncontrol baseline", "ITT, points\n(SE)", "p", "IV, points\n(SE)", "Difference in ITT,\np [95% CI]",
          "Per cent of\ncontrol baseline ", "ITT, points \n(SE)", "p ", "IV, points \n(SE)", "Difference in ITT, \np [95% CI]")
pick <- function(o, gv, gg) res %>% filter(outcome == o, g == gv, if (is.na(gg)) is.na(group) else group == gg)
O1 <- unname(outcomes["whoqol"]); O2 <- unname(outcomes["environmental"])

row_of <- function(label, gv, gg) {
  r1 <- pick(O1, gv, gg); r2 <- pick(O2, gv, gg)
  c(label, sprintf("%d / %d%s", r1$n_invited, r1$n_not_invited, if (r1$few_controls) paste0(" ", DAG) else ""),
    sprintf("%.0f%%", 100 * r1$takeup), half(r1), half(r2))
}
rows <- list(row_of("All, endline sample", "all", NA))
bold <- integer(0)
for (k in seq_len(nrow(spec))) {
  gv <- spec$g[k]
  rows[[length(rows) + 1]] <- c(spec$title[k], "", "", "", "", "", "", dcell(pick(O1, gv, 1)),
                                "", "", "", "", dcell(pick(O2, gv, 1)))
  bold <- c(bold, length(rows))
  rows[[length(rows) + 1]] <- row_of(paste0("   ", spec$lab1[k]), gv, 1)
  rows[[length(rows) + 1]] <- row_of(paste0("   ", spec$lab0[k]), gv, 0)
}
tab <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE) %>% set_names(cols)

# pieces of the note, computed from the results
n85    <- nrow(s85)
miss   <- res %>% filter(g != "all", outcome == O1) %>%
  group_by(characteristic) %>% summarise(n = sum(n), .groups = "drop") %>% filter(n < n85)
short  <- c("Lawton autonomy at baseline" = "Lawton", "Main caregiver" = "main caregiver")
miss_txt <- if (nrow(miss) == 0) "" else paste0(
  " (", paste(unique(miss$n), collapse = " and "), " in the ",
  paste(coalesce(short[miss$characteristic], miss$characteristic), collapse = " and "), " splits",
  if (all(n85 - miss$n == 1)) ", where one answer is missing)" else ", where answers are missing)")
small  <- res %>% filter(outcome == O1, few_controls) %>% arrange(n_not_invited)
tu_all <- pick(O1, "all", NA)$takeup
tu_rng <- res %>% filter(outcome == O1, g != "all") %>% pull(takeup) %>% range()
weak   <- res %>% filter(weak_iv)
count_txt <- function(x) if (x == 0) "none does" else if (x == 1) "one does" else paste(x, "do")
holm_txt  <- if (observed == 0) {
  sprintf("so none survives a Holm adjustment for the %d tests either (smallest adjusted p = %.2f)",
          n_tests, min_holm)
} else {
  sprintf("and after a Holm adjustment for the %d tests %s (smallest adjusted p = %.2f)",
          n_tests, count_txt(holm_obs), min_holm)
}

notes <- c(
  paste0(
    "Sample: the ", n85, " people observed at baseline and at endline", miss_txt, ". Within each ",
    "subgroup the endline score is regressed on its own baseline value and the randomised invitation ",
    "(ANCOVA, intention to treat, ITT) with robust HC2 standard errors; p is the p value of that ITT. ",
    "Scores are on the WHO 0 to 100 scale. 'Per cent of control baseline' is the ITT divided by the ",
    "baseline mean of the never invited people in the same subgroup, so the two subgroups have different ",
    "bases and their per cents should not be compared with each other; compare subgroups with the ",
    "difference column, in points. 'Invited who attended' is the share ",
    "of the subgroup's invited people who attended the centre on more than five days by the endline ",
    "survey (", sprintf("%.0f", 100 * tu_all), " per cent overall, ",
    nb(sprintf("%.0f to %.0f", 100 * tu_rng[1], 100 * tu_rng[2])), " across subgroups), so a gap in the ",
    "ITT can reflect a gap in take-up rather than in the effect of attending. IV is the effect of ",
    "attending, with attendance instrumented by the invitation within the subgroup and the baseline ",
    "value in both stages",
    if (nrow(weak) > 0) paste0("; ", DDAG, " marks a robust first-stage F below ", WEAK, " (weak instrument); ",
                               "the first stage includes the outcome's own baseline value, so the F and the flag ",
                               "can differ between the two outcomes for the same subgroup")
    else paste0(" (robust first-stage F of ", sprintf("%.0f", min(res$fs_F)), " or more in every subgroup)"),
    if (any(!res$iv_estimable)) "; n/e: not estimable, nobody invited in the subgroup attended" else "",
    "."),
  paste0(
    "'Difference in ITT' is the first subgroup's ITT minus the second's, in points, with its HC2 p value and ",
    "95 per cent confidence interval in brackets, from one ",
    "regression on the whole sample in which the intercept, the baseline slope and the effect of the ",
    "invitation all differ by subgroup, so it equals the gap between the two ITTs shown. Heterogeneity is ",
    "judged by this difference column, not by whether either subgroup's own ITT is significant. The ",
    "analysis is exploratory: ", n_tests, " differences are tested (", nrow(spec), " characteristics, two ",
    "outcomes), so about ", sprintf("%.1f", expected), " would fall below ", nb("p = 0.10"),
    " by chance alone; ", count_txt(observed), ", ", gsub("p = ", nb("p = "), holm_txt, fixed = TRUE), "."),
  paste0(
    "Splits at the median of this sample: age ", nb(sprintf("%.1f years", med$age0)), " at the baseline ",
    "interview, baseline Lawton score ", nb(sprintf("%g of 8", med$lawton0)), " activities and baseline ",
    "WHOQOL-BREF ", nb(sprintf("%.1f points", med$who0)), " (at or above versus below). Education (last ",
    "level completed; lower secondary is the Mexican secundaria) and the number of people living with the ",
    "older adult are split above versus at or below the median (primary; ", med$hh_size, " people), ",
    "because many people sit exactly on it. Main caregiver: the caregiver said the older adult is their ",
    "parent, or gave the same relation in reverse ('son or daughter'). ", DAG, " fewer than ", FEW,
    " not invited people in the subgroup (",
    paste(sprintf("%d for '%s'", small$n_not_invited, tolower(small$subgroup)), collapse = ", "),
    "): those estimates rest on very few controls."))

render(tab,
  spans = list(list(label = "WHOQOL-BREF", l = 4, r = 8),
               list(label = "Environmental domain", l = 9, r = 13)),
  file  = paste0(OUT_TAB, "tableA_heterogeneity.png"),
  title = "Effect of the invitation at endline by baseline characteristics (exploratory)",
  bold_rows = bold,
  notes = notes)
