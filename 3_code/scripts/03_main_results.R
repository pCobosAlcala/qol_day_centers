# 03_main_results.R  (revision, September 2026; data core and 0 to 100 scale, October 2026)
#
# Main results under the specification agreed on 16 September 2026:
#   ANCOVA run separately per wave:   Y_iw = a0 + a1*Y_i0 + a2*Z_i + e_i
#   Sample: everyone with a baseline plus that wave (102 at midline, 85 at endline)
#   No covariates. Robust HC2 standard errors (one row per person).
#   IV leads, ITT alongside: attendance by that wave's survey (more than five days at the centre)
#   instrumented by the randomised invitation, with the baseline value in both stages.
#   WHOQOL-BREF: mean of the standard 26 items; every score on the WHO 0 to 100 scale,
#   25 x (item mean minus 1). Magnitudes: per cent of the control group's baseline mean first,
#   then points and standard deviations.
#
# Data: build_analysis_data() in 00_data.R (scores rebuilt from the raw items, wave specific
# attendance from the centre's register, all randomised people). Changes against the version
# of 23 September, which read pm_who.csv (copies of the old script and outputs are in the
# _previous subfolders, suffix _before_00data):
#   - scale: 0 to 100 instead of item means on 1 to 5, so estimates, standard errors and
#     confidence bounds are 25 times larger; p values, F statistics and SD effects do not change.
#     Per cents of the control baseline DO change, because the scale's zero moves (see section 6).
#   - midline: three invited people who entered the centre after their midline survey are now in
#     the sample (99 + 3 = 102, z = 1, attended = 0), and attendance is measured by the midline
#     survey date, which moves one person with 11 days by then from "not attending" to
#     "attending". Midline estimates therefore change.
#   - endline: same 85 people and the same attendance coding (checked in 00_data.R), so the
#     endline estimates are exactly 25 times the old ones.
#
# Outputs: three PNG tables in 3_code/output/tables/, the underlying numbers in 3_code/output/csv/
# (aggregates only). Run from the imss_paper folder:
#   "C:/Program Files/R/R-4.6.1/bin/Rscript.exe" 3_code/scripts/03_main_results.R

rm(list = ls())
suppressPackageStartupMessages({
  library(tidyverse); library(AER); library(sandwich); library(lmtest)
  library(gridExtra); library(gtable); library(grid)
})
options(scipen = 999, width = 200)

source("3_code/scripts/00_paths.R")
source("3_code/scripts/00_data.R")

# ---- 1. data: one row per randomised person, scores on 0 to 100 ---------------------------
dat <- build_analysis_data()                 # checks against the old files run inside
outcomes <- OUTCOMES                         # whoqol, general, physical, psychological, social, environmental

cat("\nattendance by wave against the old time invariant d, in each ANCOVA sample of the composite:\n")
print(as.data.frame(attr(dat, "checks")$d_vs_old), row.names = FALSE)

# ---- 2. estimation ----------------------------------------------------------------------
rob <- function(m) {                                     # HC2; HC1 only if HC2 fails, and say so
  tryCatch(vcovHC(m, type = "HC2"),
           error = function(e) { warning("HC2 failed, HC1 used: ", conditionMessage(e))
                                 vcovHC(m, type = "HC1") })
}
pull_est <- function(m, term) {
  V  <- rob(m); ct <- coeftest(m, vcov. = V); ci <- coefci(m, vcov. = V, level = 0.95)
  tibble(est = ct[term, 1], se = ct[term, 2], p = ct[term, 4], lo = ci[term, 1], hi = ci[term, 2])
}

res <- list(); fs <- list(); mag <- list()
for (v in names(outcomes)) {
  for (wv in 1:2) {
    s   <- ancova_sample(dat, v, wv)                     # id, z, d, y0, y, weeks, dpw at this wave
    itt <- lm(y ~ y0 + z, s)
    iv  <- ivreg(y ~ d + y0 | z + y0, data = s)
    f1  <- lm(d ~ y0 + z, s)
    fz  <- pull_est(f1, "z")
    res[[paste(v, wv)]] <- bind_rows(
      pull_est(itt, "z") %>% mutate(model = "ITT"),
      pull_est(iv,  "d") %>% mutate(model = "IV")) %>%
      mutate(outcome = unname(outcomes[v]), wave = wv, n = nrow(s), fsF = (fz$est / fz$se)^2)
    # magnitudes, all relative to the control group's own baseline mean in this sample;
    # SD = standard deviation of the baseline outcome in this sample (both arms)
    c0 <- s %>% filter(z == 0)
    mag[[paste(v, wv)]] <- tibble(
      outcome = unname(outcomes[v]), wave = wv, n = nrow(s),
      control_base       = mean(c0$y0),
      sd0                = sd(s$y0),
      control_change_pts = mean(c0$y) - control_base,
      control_change_pct = 100 * control_change_pts / control_base,
      control_change_sd  = control_change_pts / sd0,
      itt_pts            = unname(coef(itt)["z"]),
      itt_pct            = 100 * itt_pts / control_base,
      itt_sd             = itt_pts / sd0,
      iv_pts             = unname(coef(iv)["d"]),
      iv_pct             = 100 * iv_pts / control_base,
      iv_sd              = iv_pts / sd0)
    # first stage for three exposure measures, reported once, on the composite's sample
    if (v == "whoqol") {
      for (ex in c("d", "weeks", "dpw")) {
        m <- lm(as.formula(paste(ex, "~ y0 + z")), s); r <- pull_est(m, "z")
        fs[[paste(ex, wv)]] <- r %>% mutate(
          exposure = case_when(ex == "d"     ~ "Attended the centre (more than five days)",
                               ex == "weeks" ~ "Weeks attended",
                               TRUE          ~ "Days attended per week"),
          wave = wv, n = nrow(s), F = (r$est / r$se)^2,
          invited_mean = mean(s[[ex]][s$z == 1]), control_mean = mean(s[[ex]][s$z == 0]))
      }
    }
  }
}
res <- bind_rows(res); fs <- bind_rows(fs); mag <- bind_rows(mag)
write_csv(res, paste0(OUT_CSV, "main_results.csv"))
write_csv(fs,  paste0(OUT_CSV, "first_stage.csv"))
write_csv(mag, paste0(OUT_CSV, "magnitudes.csv"))

cat("\nn per wave:", paste(unique(res$n[res$wave == 1]), collapse = "/"), "and",
    paste(unique(res$n[res$wave == 2]), collapse = "/"), "\n\n")
print(res %>% mutate(across(c(est, se, p, lo, hi, fsF), ~ round(.x, 3))) %>%
        select(outcome, wave, model, est, se, p, lo, hi, fsF) %>% as.data.frame(), row.names = FALSE)
cat("\n")
print(fs %>% mutate(across(c(est, se, p, F, invited_mean, control_mean), ~ round(.x, 3))) %>%
        select(exposure, wave, n, est, se, p, F, invited_mean, control_mean) %>% as.data.frame(),
      row.names = FALSE)
cat("\n")
print(mag %>% mutate(across(where(is.double), ~ round(.x, 2))) %>% as.data.frame(), row.names = FALSE)

# ---- 3. rendering helpers ---------------------------------------------------------------
star <- function(p) ifelse(p < .01, "***", ifelse(p < .05, "**", ifelse(p < .1, "*", "")))
# points on the 0 to 100 scale, two decimals throughout (one decimal would print the midline
# environmental lower bound of -0.026 as "-0.0"); the 1 to 5 version used three and two
cell <- function(r) sprintf("%.2f%s\n(%.2f)\n[%.2f, %.2f]", r$est, star(r$p), r$se, r$lo, r$hi)

render <- function(df, spans, file, title, notes, left_cols = 1, note_width = 120) {
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
  # spanner row (taller when a label has two lines)
  two <- any(grepl("\n", vapply(spans, function(s) s$label, character(1)), fixed = TRUE))
  g <- gtable_add_rows(g, heights = unit(if (two) 11 else 7, "mm"), pos = 0)
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

# ---- 4. Table: main results -------------------------------------------------------------
n1 <- unique(res$n[res$wave == 1]); n2 <- unique(res$n[res$wave == 2])
stopifnot(length(n1) == 1, length(n2) == 1)              # same sample for all six outcomes
main <- map_dfr(outcomes, function(o) {
  g <- function(w, m) cell(res %>% filter(outcome == o, wave == w, model == m))
  tibble(Outcome = o,
         `Attending\n(IV)` = g(1, "IV"), `Invited\n(ITT)` = g(1, "ITT"),
         `Attending \n(IV)` = g(2, "IV"), `Invited \n(ITT)` = g(2, "ITT"))
})
fr <- function(w) { f <- res %>% filter(wave == w, model == "IV") %>% pull(fsF); sprintf("%.1f to %.1f", min(f), max(f)) }
main <- bind_rows(main,
  tibble(Outcome = "First-stage F", `Attending\n(IV)` = fr(1), `Invited\n(ITT)` = "",
         `Attending \n(IV)` = fr(2), `Invited \n(ITT)` = ""),
  tibble(Outcome = "Individuals", `Attending\n(IV)` = as.character(n1), `Invited\n(ITT)` = as.character(n1),
         `Attending \n(IV)` = as.character(n2), `Invited \n(ITT)` = as.character(n2)))
render(main,
  spans = list(list(label = "Midline\n(about 3 months after entry)", l = 2, r = 3),
               list(label = "Endline\n(about 9 months after entry)", l = 4, r = 5)),
  file  = paste0(OUT_TAB, "table_main_results.png"),
  title = "Effect of the day centre on quality of life (WHOQOL-BREF)",
  notes = paste(
    "Each cell reports the estimate in points of the WHO 0 to 100 scale, the robust (HC2) standard",
    "error in parentheses and the 95 per cent confidence interval in brackets. ANCOVA: the outcome at",
    "each wave is regressed on its own baseline value and the randomised invitation, separately for each",
    "wave. ITT is the effect of the invitation (the offer of a place). IV is the effect of attending, with",
    "attendance instrumented by the invitation and the baseline value in both stages; attending means more",
    "than five days at the centre by that wave's survey. Each wave uses everyone observed at baseline and at",
    "that wave. Months are counted from entry to the centre for those attending; those not attending were",
    "surveyed about six months and a year after baseline. First-stage F is the robust F on the excluded",
    "instrument, shown as the range across the",
    "six outcomes. WHOQOL-BREF is the mean of the standard 26 items; all scores are on the WHO 0 to 100",
    "scale (25 times the item mean minus 1). *** p < 0.01, ** p < 0.05, * p < 0.1."))

# ---- 5. Table: first stage --------------------------------------------------------------
fs_tab <- map_dfr(unique(fs$exposure), function(e) {
  g <- function(w) { r <- fs %>% filter(exposure == e, wave == w)
    sprintf("%.3f%s\n(%.3f)\nF = %.1f", r$est, star(r$p), r$se, r$F) }
  m <- function(w) { r <- fs %>% filter(exposure == e, wave == w)
    sprintf("%.2f vs %.2f", r$invited_mean, r$control_mean) }
  tibble(`Exposure` = e, `Effect of\ninvitation` = g(1), `Invited vs\nnever invited` = m(1),
         `Effect of \ninvitation` = g(2), `Invited vs \nnever invited` = m(2))
})
# attendance among those attending by endline, for the note, from the data
att2 <- ancova_sample(dat, "whoqol", 2) %>% filter(d == 1)
stopifnot(all(fs$control_mean == 0))                     # nobody uninvited attended
render(fs_tab,
  spans = list(list(label = sprintf("Midline, n = %d", n1), l = 2, r = 3),
               list(label = sprintf("Endline, n = %d", n2), l = 4, r = 5)),
  file  = paste0(OUT_TAB, "table_first_stage.png"),
  title = "First stage: what the invitation delivered",
  notes = paste(
    "Each exposure is regressed on the randomised invitation and baseline WHOQOL-BREF. Robust (HC2) standard",
    "errors in parentheses; F is the robust F on the invitation. Attending means more than five days at the",
    "centre up to that wave's survey; weeks attended are those days divided by five; days per week are days",
    "attended per seven calendar days since the first visit. Nobody who was not invited attended the centre,",
    "so non-compliance is one-sided and all three exposures are zero in that group.",
    sprintf("The %d people attending by endline had attended %.0f weeks on average, at %.1f days per week.",
            nrow(att2), mean(att2$weeks), mean(att2$dpw)),
    "*** p < 0.01, ** p < 0.05, * p < 0.1."), note_width = 110)

# ---- 6. Table: magnitudes ---------------------------------------------------------------
# Pablo's format: per cent of the control group's baseline mean first, then in parentheses the
# same number in points of the 0 to 100 scale and in baseline standard deviations.
# Per cents depend on where the scale starts: on 0 to 100 the zero is the worst possible score,
# so the control baseline (the denominator) is smaller than on the 1 to 5 item scale and the per
# cents are larger than in the September version. Points and SD do not depend on the scale.
cellm <- function(pct, pts, sdu) sprintf("%+.1f%%\n(%+.1f points,\n%+.2f SD)", pct, pts, sdu)
mag_tab <- map_dfr(outcomes, function(o) {
  a <- mag %>% filter(outcome == o, wave == 1); b <- mag %>% filter(outcome == o, wave == 2)
  tibble(Outcome = o,
         `Control\ngroup`    = cellm(a$control_change_pct, a$control_change_pts, a$control_change_sd),
         `Invited\n(ITT)`    = cellm(a$itt_pct, a$itt_pts, a$itt_sd),
         `Attending\n(IV)`   = cellm(a$iv_pct, a$iv_pts, a$iv_sd),
         `Control \ngroup`   = cellm(b$control_change_pct, b$control_change_pts, b$control_change_sd),
         `Invited \n(ITT)`   = cellm(b$itt_pct, b$itt_pts, b$itt_sd),
         `Attending \n(IV)`  = cellm(b$iv_pct, b$iv_pts, b$iv_sd))
})
w2 <- mag %>% filter(outcome == outcomes[["whoqol"]], wave == 2)  # for the reading example
render(mag_tab,
  spans = list(list(label = "Midline\n(about 3 months after entry)", l = 2, r = 4),
               list(label = "Endline\n(about 9 months after entry)", l = 5, r = 7)),
  file  = paste0(OUT_TAB, "table_magnitudes.png"),
  title = "How large are the effects",
  notes = paste(
    "Each cell gives a change as a per cent of the control group's own baseline mean on the WHO 0 to 100",
    "scale, then in parentheses the same change in points of that scale and in standard deviations of the",
    "baseline outcome in the same sample. 'Control group' is how much that group changed over the period,",
    "without the programme. 'Invited' and 'Attending' are the ITT and IV estimates of Table 3. Months are",
    "counted from entry to the centre for those attending; those not attending were surveyed about six",
    "months and a year after baseline.",
    sprintf(paste("Reading the endline row for WHOQOL-BREF: the control group's score %s by %.1f per cent",
                  "of its baseline mean of %.1f points, and the invitation raised scores by %.1f per",
                  "cent of that baseline (%.1f points)."),
            ifelse(w2$control_change_pts < 0, "fell", "rose"), abs(w2$control_change_pct),
            w2$control_base, w2$itt_pct, w2$itt_pts)), note_width = 118)
