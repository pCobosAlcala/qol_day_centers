# 08_text_numbers.R  (revision, October 2026)
#
# EVERY NUMBER THE MANUSCRIPT TEXT QUOTES THAT NO TABLE PRODUCES, in one CSV with a short label each.
# All WHOQOL-BREF scores are on the WHO 0 to 100 scale (25 x (item mean minus 1)), built in 00_data.R.
#
# Main specification (fixed, see 03_main_results.R): ANCOVA run separately per wave, the follow-up
# score on its own baseline value plus the randomised invitation z; everyone observed at baseline and
# at that wave (102 at midline, 85 at endline); robust HC2 errors; IV = attending (more than five days
# at the centre up to that wave's survey) instrumented by z, baseline in both stages; no covariates.
#
# Sections of the output (column `section`):
#   a  ANCOVA baseline coefficient, composite, each wave, and its HC2 test against 1
#      (1 is what a difference in differences or change score model imposes)
#   b  dose among those attending by endline (d = 1): endline score on baseline plus weeks attended,
#      and on baseline plus days per week (HC2); range of days attended among them
#   c  baseline Lawton of heavy versus light attenders, endline sample, two definitions
#   d  timing: months from entry to each survey (attenders) or baseline to each survey (the rest),
#      and the joint test that time to endline is unrelated to six baseline characteristics
#   e  first stage F, classical and robust (HC2), each wave
#   f  randomisation inference p values for the ITT (91 of 120 redrawn as invited, then the ITT is
#      estimated on that wave's sample), composite and environmental, both waves
#   g  multiple testing: Holm and Bonferroni adjusted p over the six outcomes, each wave, ITT and IV
#   h  suggestive check of the exclusion restriction: invited people who did not attend against the
#      never invited (SELF SELECTED comparison, not experimental)
#   i  midline robustness: dropping invitees surveyed at midline before their own invitation (a range,
#      1 to 4 people, because the dates of rounds 7 and 8 are uncertain), and the round 1 invitee whose
#      baseline is dated after round 1; both waves without the round 7 and 8 invitees; and counts that
#      document the dating anomaly of rounds 7 and 8
#   j  de Chaisemartin and Behaghel (2020) waiting list correction (weight (S minus 1)/S on invited
#      takers), composite and environmental, ITT and IV, both waves
#   k  attrition by group (attended anywhere in the register, invited not attending, never invited)
#      and by arm, each wave
#   l  context (draft 3, October 2026): substitutes for the centre among those not attending
#      (caregiver item 06.05), use of the centre from its attendance register (days of operation,
#      people per day), the endline composite with the sociodemographic covariates dropping one at
#      a time, and WHOQOL-BREF items missing in the surveys held
#
# Columns: section, item (machine name), label (what the number is), value, se, p, lo, hi, n, note.
#   lo and hi are a 95 per cent confidence interval for estimates and a min to max range for
#   descriptive rows; the label says which. In section (g) value is the estimate and p is the ADJUSTED
#   p (the unadjusted HC2 p is in the note). Aggregates only: no individual row, name or date is
#   written (survey and entry dates are used internally, as 00_data.R allows with dates = TRUE).
#
# Output: 3_code/output/csv/text_numbers.csv. Run from the imss_paper folder:
#   "C:/Program Files/R/R-4.6.1/bin/Rscript.exe" 3_code/scripts/08_text_numbers.R

rm(list = ls())
suppressPackageStartupMessages({
  library(tidyverse); library(AER); library(sandwich); library(lmtest)
})
options(scipen = 999, width = 200)

source("3_code/scripts/00_paths.R")
source("3_code/scripts/00_data.R")

# =============================================================================================
# 0. DATA AND HELPERS
# =============================================================================================
dat <- build_analysis_data(dates = TRUE)      # one row per randomised person (120); dates internal only
n_missing_items <- attr(dat, "checks")$n_missing_items   # for (l); joins below drop the attribute
stopifnot(!is.null(n_missing_items))
OUT <- c(whoqol = "WHOQOL-BREF", environmental = "Environmental")   # the two outcomes the text quotes
WAVE <- c("midline", "endline")               # wave 1, wave 2
MONTH <- 365.25 / 12                          # days per month, for the timing section

# ---- days attended in the WHOLE register (18 Jan 2023 to 16 Feb 2024) ----------------------
# 00_data.R computes the register total internally but does not export it; (b), (j) and (k) need it.
# Until 00_data.R exports days_all, read the register exactly as 00_data.R does (one column per
# opening day, header = Excel date serial, "x" or "X" = attended) and check the parse against what
# 00_data.R does export:
#   - entry_date (first day attended) for ALL 120, which checks the whole register row of everyone,
#     including people not surveyed at a wave;
#   - days1 and days2 for everyone surveyed at that wave (102 and 85 people; NA otherwise);
#   - days_all is at least days1 and days2, and zero exactly when there is no entry date.
reg <- readxl::read_excel(paste0(RAW, "CD CASSAAM 16febrero24 asistencias.xlsx"), .name_repair = "minimal")
serial <- suppressWarnings(as.numeric(names(reg)))
is_day <- !is.na(serial) & serial > 40000 & serial < 50000
stopifnot(sum(is_day) == 223, sum(names(reg) == "Folio") == 1, sum(names(reg) == "Ronda") == 1)
marks <- reg[, c(which(names(reg) == "Folio"), which(is_day))]
names(marks) <- c("id", as.character(as.Date(serial[is_day], origin = "1899-12-30")))
visits <- marks %>%
  pivot_longer(-id, names_to = "day", values_to = "mark", values_transform = list(mark = as.character)) %>%
  filter(tolower(mark) %in% "x") %>% transmute(id, day = as.Date(day))
chk_days <- dat %>% select(id, date1, date2, days1, days2, entry_date) %>%
  left_join(visits, by = "id", relationship = "one-to-many") %>%
  group_by(id, date1, date2, days1, days2, entry_date) %>%
  summarise(days_all = sum(!is.na(day)),
            entry_chk = if (all(is.na(day))) as.Date(NA) else min(day, na.rm = TRUE),
            d1_chk = sum(day <= date1[1], na.rm = TRUE), d2_chk = sum(day <= date2[1], na.rm = TRUE),
            .groups = "drop")
stopifnot(nrow(chk_days) == 120,
          identical(is.na(chk_days$entry_date), is.na(chk_days$entry_chk)),          # all 120
          all(chk_days$entry_date == chk_days$entry_chk, na.rm = TRUE),
          identical(chk_days$days_all == 0, is.na(chk_days$entry_date)),
          identical(is.na(chk_days$days1), is.na(chk_days$date1)),                  # surveyed = counted
          identical(is.na(chk_days$days2), is.na(chk_days$date2)),
          all(chk_days$days1 == chk_days$d1_chk, na.rm = TRUE),                     # 102 at midline
          all(chk_days$days2 == chk_days$d2_chk, na.rm = TRUE),                     # 85 at endline
          all(chk_days$days_all >= chk_days$days1, na.rm = TRUE),
          all(chk_days$days_all >= chk_days$days2, na.rm = TRUE))
dat <- dat %>% left_join(chk_days %>% select(id, days_all), by = "id") %>%
  mutate(ever_att = as.integer(days_all > 5))      # attended more than five days in the whole register
stopifnot(all(dat$ever_att[dat$z == 0] == 0))      # one sided non compliance, whole register
# the register's own round label (column "Ronda"), used only to document the round 7 and 8 anomaly in (i)
reg_round <- tibble(id = as.character(reg[["Folio"]]),       # base extraction: the sheet has blank
                    ronda = suppressWarnings(as.integer(as.character(reg[["Ronda"]]))))  # column names
stopifnot(!anyDuplicated(reg_round$id[!is.na(reg_round$id)]))

# ---- estimation helpers (HC2, as 03_main_results.R) ----------------------------------------
rob <- function(m) {
  tryCatch(vcovHC(m, type = "HC2"),
           error = function(e) { warning("HC2 failed, HC1 used: ", conditionMessage(e))
                                 vcovHC(m, type = "HC1") })
}
pull_est <- function(m, term) {
  V  <- rob(m); ct <- coeftest(m, vcov. = V); ci <- coefci(m, vcov. = V, level = 0.95)
  tibble(est = ct[term, 1], se = ct[term, 2], p = ct[term, 4], lo = ci[term, 1], hi = ci[term, 2],
         n = length(residuals(m)))
}

# ---- output accumulator: one row per quoted number -----------------------------------------
rows <- list()
add <- function(section, item, label, value, se = NA, p = NA, lo = NA, hi = NA, n = NA, note = "") {
  rows[[length(rows) + 1]] <<- tibble(section = section, item = item, label = label,
                                      value = as.numeric(value), se = as.numeric(se), p = as.numeric(p),
                                      lo = as.numeric(lo), hi = as.numeric(hi), n = as.numeric(n),
                                      note = note)
}
add_est <- function(section, item, label, r, note = "") {      # r = one row of pull_est()
  add(section, item, label, r$est, r$se, r$p, r$lo, r$hi, r$n, note)
}

# =============================================================================================
# (a) ANCOVA BASELINE COEFFICIENT AND ITS TEST AGAINST 1 (composite, HC2)
# =============================================================================================
# A change score (or individual fixed effects) model forces the coefficient on baseline to 1.
# The ANCOVA estimates it; a value well below 1 is why ANCOVA is the more precise estimator.
for (wv in 1:2) {
  s <- ancova_sample(dat, "whoqol", wv)
  m <- lm(y ~ y0 + z, s); V <- rob(m)
  b <- unname(coef(m)["y0"]); se <- sqrt(V["y0", "y0"])
  t1 <- (b - 1) / se; p1 <- 2 * pt(-abs(t1), df.residual(m))         # HC2 Wald test of b = 1
  ci <- b + c(-1, 1) * qt(0.975, df.residual(m)) * se
  add("a", paste0("baseline_coef_", WAVE[wv]),
      sprintf("WHOQOL-BREF, %s ANCOVA: coefficient on baseline score (HC2 SE, 95%% CI)", WAVE[wv]),
      b, se, NA, ci[1], ci[2], nrow(s))
  add("a", paste0("baseline_coef_eq1_", WAVE[wv]),
      sprintf("WHOQOL-BREF, %s: HC2 test that the baseline coefficient equals 1 (value = F, 1 df)", WAVE[wv]),
      t1^2, NA, p1, NA, NA, nrow(s), note = "p is the test of equality to 1")
}

# =============================================================================================
# (b) DOSE AMONG THOSE ATTENDING BY ENDLINE (d = 1 in the endline sample; not experimental)
# =============================================================================================
for (v in names(OUT)) {
  s <- ancova_sample(dat, v, 2) %>% filter(d == 1)
  for (ex in c("weeks", "dpw")) {
    m <- lm(as.formula(paste("y ~ y0 +", ex)), s)
    add_est("b", paste0("dose_", ex, "_", v),
            sprintf("%s endline score on baseline plus %s, among those attending by endline (points per %s; HC2)",
                    OUT[[v]], ifelse(ex == "weeks", "weeks attended", "days attended per week"),
                    ifelse(ex == "weeks", "week", "day per week")),
            pull_est(m, ex), note = "association only: who attends more is a choice")
  }
}
att2 <- dat %>% filter(!is.na(whoqol_y0), !is.na(whoqol_y2), d2 == 1)
add("b", "dose_days_by_survey_range", "Days attended up to the endline survey, those attending by endline: mean (range in lo to hi)",
    mean(att2$days2), lo = min(att2$days2), hi = max(att2$days2), n = nrow(att2))
add("b", "dose_days_register_range", "Days attended in the whole register, those attending by endline: mean (range in lo to hi)",
    mean(att2$days_all), lo = min(att2$days_all), hi = max(att2$days_all), n = nrow(att2),
    note = "register runs 18 January 2023 to 16 February 2024")
add("b", "dose_weeks_by_survey_range", "Weeks attended up to the endline survey, those attending by endline: mean (range in lo to hi)",
    mean(att2$weeks2), lo = min(att2$weeks2), hi = max(att2$weeks2), n = nrow(att2))
add("b", "dose_dpw_by_survey_range", "Days per week up to the endline survey, those attending by endline: mean (range in lo to hi)",
    mean(att2$dpw2), lo = min(att2$dpw2), hi = max(att2$dpw2), n = nrow(att2))
e85 <- dat %>% filter(!is.na(whoqol_y0), !is.na(whoqol_y2))
add("b", "brief_attenders_n", "Endline sample: invited people with 1 to 5 days by the survey (counted as not attending)",
    sum(e85$z == 1 & e85$days2 %in% 1:5), n = nrow(e85))
add("b", "gap_6_to_35_days_n", "Endline sample: people with 6 to 35 days by the survey (none expected)",
    sum(e85$days2 %in% 6:35), n = nrow(e85))

# =============================================================================================
# (c) BASELINE LAWTON (0 to 8 instrumental activities) OF HEAVY VERSUS LIGHT ATTENDERS, endline sample
# =============================================================================================
# Definition 1: invited people in the endline sample who attended at least once by the endline
#   survey (days2 >= 1, which includes those with 1 to 5 days); heavy = 20 weeks (100 days) or more
#   by the survey, light = less.
# Definition 2: those counted as attending by endline (d2 = 1, more than 5 days); heavy = more than
#   100 days by the survey, light = fewer than 100 (nobody has exactly 100, checked).
# Tests: Welch t test (as in the September note) and Wilcoxon rank sum (Lawton is a 0 to 8 count).
lawton_cmp <- function(df, tag, label) {
  df <- df %>% filter(!is.na(lawton0))
  h <- df$lawton0[df$heavy]; l <- df$lawton0[!df$heavy]
  pt_ <- t.test(h, l)$p.value
  pw_ <- suppressWarnings(wilcox.test(h, l, exact = FALSE)$p.value)
  add("c", paste0("lawton_heavy_", tag), paste0(label, ": mean baseline Lawton, heavy"), mean(h), n = length(h))
  add("c", paste0("lawton_light_", tag), paste0(label, ": mean baseline Lawton, light"), mean(l), n = length(l))
  add("c", paste0("lawton_diff_", tag), paste0(label, ": heavy minus light (p = Welch t test)"),
      mean(h) - mean(l), p = pt_, n = length(h) + length(l))
  add("c", paste0("lawton_wilcoxon_", tag), paste0(label, ": Wilcoxon rank sum p"), NA, p = pw_,
      n = length(h) + length(l))
}
d1_def <- e85 %>% filter(z == 1, days2 >= 1) %>% mutate(heavy = weeks2 >= 20)
stopifnot(!any(e85$days2 == 100 & e85$d2 == 1))
d2_def <- e85 %>% filter(d2 == 1) %>% mutate(heavy = days2 > 100)
lawton_cmp(d1_def, "20weeks", "Invited who attended at least once by endline, 20 weeks or more vs less")
lawton_cmp(d2_def, "100days", "Attending by endline (more than 5 days), more than 100 days vs fewer")
add("c", "lawton_missing_20weeks", "Definition 1 group: people without a baseline Lawton score (dropped)",
    sum(is.na(d1_def$lawton0)))

# =============================================================================================
# (d) TIMING OF THE FOLLOW UP SURVEYS
# =============================================================================================
# The follow up clock starts at ENTRY for those attending by that survey (d = 1) and at the BASELINE
# interview for everyone else. Months = days / (365.25 / 12).
tim <- function(wv) {
  dw <- paste0("d", wv); dt <- paste0("date", wv)
  dat %>% filter(!is.na(whoqol_y0), !is.na(.data[[paste0("whoqol_y", wv)]])) %>%
    mutate(att = .data[[dw]] == 1,
           m_entry = as.numeric(.data[[dt]] - entry_date) / MONTH,
           m_base  = as.numeric(.data[[dt]] - date0) / MONTH,
           m_clock = if_else(att, m_entry, m_base))
}
for (wv in 1:2) {
  x <- tim(wv); a <- x %>% filter(att); r <- x %>% filter(!att)
  add("d", paste0("months_entry_to_", WAVE[wv], "_attenders"),
      sprintf("Months from entry to the %s survey, those attending by then: mean (range in lo to hi)", WAVE[wv]),
      mean(a$m_entry), lo = min(a$m_entry), hi = max(a$m_entry), n = nrow(a))
  add("d", paste0("months_entry_to_", WAVE[wv], "_attenders_sd"),
      sprintf("Months from entry to the %s survey, those attending by then: SD", WAVE[wv]), sd(a$m_entry), n = nrow(a))
  add("d", paste0("months_base_to_", WAVE[wv], "_rest"),
      sprintf("Months from baseline to the %s survey, everyone else: mean (range in lo to hi)", WAVE[wv]),
      mean(r$m_base), lo = min(r$m_base), hi = max(r$m_base), n = nrow(r))
  add("d", paste0("months_base_to_", WAVE[wv], "_rest_sd"),
      sprintf("Months from baseline to the %s survey, everyone else: SD", WAVE[wv]), sd(r$m_base), n = nrow(r))
  add("d", paste0("months_base_to_", WAVE[wv], "_all"),
      sprintf("Months from baseline to the %s survey, whole %s sample: mean (range in lo to hi)", WAVE[wv], WAVE[wv]),
      mean(x$m_base), lo = min(x$m_base), hi = max(x$m_base), n = nrow(x))
}
x2 <- tim(2) %>% filter(att)
add("d", "months_base_to_entry_attenders", "Months from baseline to entry, those attending by endline: mean (range in lo to hi)",
    mean(as.numeric(x2$entry_date - x2$date0) / MONTH),
    lo = min(as.numeric(x2$entry_date - x2$date0) / MONTH), hi = max(as.numeric(x2$entry_date - x2$date0) / MONTH),
    n = nrow(x2))
# joint test: is time to endline related to six baseline characteristics (baseline WHOQOL-BREF, age,
# Lawton, Katz, sex, education; the set of the September check)? Classical F and HC2 Wald F.
joint_test <- function(df, yvar) {
  df <- df %>% filter(!is.na(lawton0), !is.na(katz0))
  m <- lm(as.formula(paste(yvar, "~ whoqol_y0 + age0 + lawton0 + katz0 + female + educ")), df)
  fs <- summary(m)$fstatistic
  k <- names(coef(m))[-1]; b <- coef(m)[k]; V <- rob(m)[k, k]
  Fw <- as.numeric(t(b) %*% solve(V) %*% b) / length(k)
  list(F = unname(fs[1]), p = unname(pf(fs[1], fs[2], fs[3], lower.tail = FALSE)), df1 = fs[2], df2 = fs[3],
       Fr = Fw, pr = pf(Fw, length(k), df.residual(m), lower.tail = FALSE), n = nrow(df))
}
x2all <- tim(2)
jt <- joint_test(x2all, "m_clock")
add("d", "timing_joint_F_classical", sprintf("Joint test, months to endline (entry clock for attenders, baseline clock otherwise) on six baseline characteristics: classical F (%d, %d df)", jt$df1, jt$df2),
    jt$F, p = jt$p, n = jt$n, note = "baseline WHOQOL-BREF, age, Lawton, Katz, sex, education")
add("d", "timing_joint_F_hc2", "Same joint test, HC2 Wald F",
    jt$Fr, p = jt$pr, n = jt$n, note = "baseline WHOQOL-BREF, age, Lawton, Katz, sex, education")
jb <- joint_test(x2all, "m_base")
add("d", "timing_joint_F_classical_baseclock", "Joint test with months from baseline to endline for everyone: classical F",
    jb$F, p = jb$p, n = jb$n, note = "secondary: one clock for all")
add("d", "timing_joint_F_hc2_baseclock", "Joint test with months from baseline to endline for everyone: HC2 Wald F",
    jb$Fr, p = jb$pr, n = jb$n, note = "secondary: one clock for all")

# =============================================================================================
# (e) FIRST STAGE F ON THE INVITATION, CLASSICAL AND ROBUST (HC2), EACH WAVE
# =============================================================================================
fsF <- list()
for (v in names(OUTCOMES)) for (wv in 1:2) {
  s <- ancova_sample(dat, v, wv); m <- lm(d ~ y0 + z, s)
  fsF[[paste(v, wv)]] <- tibble(outcome = v, wave = wv, n = nrow(s),
    F_classical = unname(summary(m)$coefficients["z", 3])^2,
    F_hc2 = unname(coef(m)["z"] / sqrt(rob(m)["z", "z"]))^2)
}
fsF <- bind_rows(fsF)
for (wv in 1:2) {
  f <- fsF %>% filter(wave == wv); w <- f %>% filter(outcome == "whoqol")
  add("e", paste0("fs_F_classical_", WAVE[wv]),
      sprintf("First stage F on the invitation, %s, composite sample: classical (range over six outcomes in lo to hi)", WAVE[wv]),
      w$F_classical, lo = min(f$F_classical), hi = max(f$F_classical), n = w$n)
  add("e", paste0("fs_F_hc2_", WAVE[wv]),
      sprintf("First stage F on the invitation, %s, composite sample: robust HC2 (range over six outcomes in lo to hi)", WAVE[wv]),
      w$F_hc2, lo = min(f$F_hc2), hi = max(f$F_hc2), n = w$n,
      note = "robust F exceeds classical because attendance is zero for every never invited person")
}

# =============================================================================================
# (f) RANDOMISATION INFERENCE FOR THE ITT
# =============================================================================================
# The eight rounds drew only from those not yet picked, so the design is one uniform draw of 91
# invited out of 120. Each draw re-assigns z to all 120, then the ITT ANCOVA is estimated on that
# wave's observed sample (whose number of invited then varies, as it would have). Two sided p =
# (1 + number of draws with a statistic at least as large in absolute value as the observed one) /
# (1 + number of draws): the observed assignment counts as one member of the reference set, so p is
# never 0 and the test is exact (valid at any number of draws). Shown for the coefficient
# (unstudentised) and for its HC2 t statistic (studentised). The studentised version is the one to
# quote: it stays valid when the two arms differ in outcome variance (Wu and Ding 2021, JASA), while
# the unstudentised one does not; at midline the never invited are more dispersed and the two
# versions diverge. Under the sharp null of no effect for anyone this is also the test of the IV
# estimate. Monte Carlo SE of each p = sqrt(p (1 - p) / B), in the se column.
#
# Speed: 100,000 draws x 4 outcome waves. The ITT and its HC2 SE are computed for all draws at once
# by Frisch Waugh Lovell: partial the intercept and baseline out of y and z (they do not depend on
# the draw), then b = sum(zt yt) / sum(zt^2), leverage h = h(1, y0) + zt^2 / sum(zt^2), and
# HC2 Var(b) = sum(zt^2 e^2 / (1 - h)) / sum(zt^2)^2, which is exactly the z element of the full
# HC2 sandwich. Checked against lm + sandwich at the observed assignment and against the direct
# matrix version (itt_fast) on the first 200 draws.
B_RI <- 100000
set.seed(20261002)
draws <- replicate(B_RI, { zz <- integer(nrow(dat)); zz[sample.int(nrow(dat), 91)] <- 1L; zz })  # 120 x B
itt_fast <- function(y, y0, z) {                  # ITT coefficient and its HC2 t statistic, one draw
  if (var(z) == 0) return(c(NA_real_, NA_real_))
  X <- cbind(1, y0, z); XtXi <- solve(crossprod(X)); b <- XtXi %*% crossprod(X, y)
  e <- as.vector(y - X %*% b); h <- rowSums((X %*% XtXi) * X)
  meat <- crossprod(X * (e / sqrt(1 - h))); V <- XtXi %*% meat %*% XtXi
  c(b[3], b[3] / sqrt(V[3, 3]))
}
itt_many <- function(y, y0, Z, chunk = 10000) {  # same, for every column of Z (n x B), by FWL
  W <- cbind(1, y0); P <- W %*% solve(crossprod(W), t(W))           # projection on (1, y0)
  hW <- diag(P); yt <- as.vector(y - P %*% y)
  out <- matrix(NA_real_, 2, ncol(Z))
  for (st in seq(1, ncol(Z), by = chunk)) {
    j  <- st:min(ncol(Z), st + chunk - 1)
    Zt <- Z[, j, drop = FALSE] - P %*% Z[, j, drop = FALSE]          # z partialled out
    sz <- colSums(Zt^2)
    b  <- colSums(Zt * yt) / sz
    E  <- yt - sweep(Zt, 2, b, "*")                                   # full model residuals
    H  <- hW + sweep(Zt^2, 2, sz, "/")                                # full model leverages
    se <- sqrt(colSums(Zt^2 * E^2 / (1 - H))) / sz
    bad <- sz < 1e-10                                                 # no variation in z
    out[1, j] <- ifelse(bad, NA, b); out[2, j] <- ifelse(bad, NA, b / se)
  }
  out
}
for (v in names(OUT)) for (wv in 1:2) {
  s <- ancova_sample(dat, v, wv); idx <- match(s$id, dat$id)
  obs <- itt_fast(s$y, s$y0, s$z)
  ref <- pull_est(lm(y ~ y0 + z, s), "z")                       # agrees with lm + sandwich
  stopifnot(abs(obs[1] - ref$est) < 1e-9, abs(obs[2] - ref$est / ref$se) < 1e-9,
            max(abs(itt_many(s$y, s$y0, matrix(s$z)) - obs)) < 1e-9)
  sim <- itt_many(s$y, s$y0, draws[idx, , drop = FALSE])
  chk <- sapply(1:200, function(k) itt_fast(s$y, s$y0, draws[idx, k]))
  stopifnot(max(abs(chk - sim[, 1:200]), na.rm = TRUE) < 1e-8)
  ok <- !is.na(sim[1, ]); B_ok <- sum(ok)
  p_b <- (1 + sum(abs(sim[1, ok]) >= abs(obs[1]) - 1e-12)) / (1 + B_ok)
  p_t <- (1 + sum(abs(sim[2, ok]) >= abs(obs[2]) - 1e-12)) / (1 + B_ok)
  add("f", paste0("ri_p_", WAVE[wv], "_", v),
      sprintf("%s, %s ITT: randomisation inference p, unstudentised (on the coefficient; value = observed ITT in points; se = Monte Carlo SE of p)", OUT[[v]], WAVE[wv]),
      obs[1], sqrt(p_b * (1 - p_b) / B_ok), p_b, n = nrow(s),
      note = sprintf("%d draws of 91 invited out of 120 (%d valid); p = (1 + count) / (1 + draws)", B_RI, B_ok))
  add("f", paste0("ri_p_t_", WAVE[wv], "_", v),
      sprintf("%s, %s ITT: randomisation inference p, studentised (on the HC2 t statistic; value = observed t; se = Monte Carlo SE of p). Preferred", OUT[[v]], WAVE[wv]),
      obs[2], sqrt(p_t * (1 - p_t) / B_ok), p_t, n = nrow(s),
      note = sprintf("%d draws (%d valid); p = (1 + count) / (1 + draws); HC2 p for comparison %.4f", B_RI, B_ok, ref$p))
}
rm(draws)

# =============================================================================================
# (g) MULTIPLE TESTING ACROSS THE SIX OUTCOMES (Holm and Bonferroni), each wave, ITT and IV
# =============================================================================================
mt <- list()
for (v in names(OUTCOMES)) for (wv in 1:2) {
  s <- ancova_sample(dat, v, wv)
  mt[[paste(v, wv)]] <- bind_rows(
    pull_est(lm(y ~ y0 + z, s), "z") %>% mutate(model = "ITT"),
    pull_est(ivreg(y ~ d + y0 | z + y0, data = s), "d") %>% mutate(model = "IV")) %>%
    mutate(outcome = v, wave = wv)
}
mt <- bind_rows(mt) %>% group_by(wave, model) %>%
  mutate(p_holm = p.adjust(p, "holm"), p_bonf = p.adjust(p, "bonferroni")) %>% ungroup()
# the same estimates as 03's main table, if it has been run on the current data
if (file.exists(paste0(OUT_CSV, "main_results.csv"))) {
  mr <- read_csv(paste0(OUT_CSV, "main_results.csv"), show_col_types = FALSE)
  cmp <- mt %>% mutate(outcome = unname(OUTCOMES[outcome])) %>%
    inner_join(mr, by = c("outcome", "wave", "model"), suffix = c("", "_03"))
  stopifnot(nrow(cmp) == 24, max(abs(cmp$est - cmp$est_03)) < 1e-8, max(abs(cmp$p - cmp$p_03)) < 1e-8)
}
# Layout: value = the estimate in points, se = its HC2 SE, p = the ADJUSTED p (the number to quote),
# note = the unadjusted HC2 p. So the p column always holds the p the label names.
for (i in seq_len(nrow(mt))) {
  r <- mt[i, ]
  add("g", sprintf("holm_%s_%s_%s", WAVE[r$wave], tolower(r$model), r$outcome),
      sprintf("%s, %s %s: Holm adjusted p over six outcomes (p column; value = estimate in points)", OUTCOMES[[r$outcome]], WAVE[r$wave], r$model),
      r$est, r$se, r$p_holm, n = r$n, note = sprintf("unadjusted HC2 p %.4f", r$p))
  add("g", sprintf("bonf_%s_%s_%s", WAVE[r$wave], tolower(r$model), r$outcome),
      sprintf("%s, %s %s: Bonferroni adjusted p over six outcomes (p column; value = estimate in points)", OUTCOMES[[r$outcome]], WAVE[r$wave], r$model),
      r$est, r$se, r$p_bonf, n = r$n, note = sprintf("unadjusted HC2 p %.4f", r$p))
}

# =============================================================================================
# (h) EXCLUSION RESTRICTION, SUGGESTIVE: INVITED NON ATTENDERS AGAINST THE NEVER INVITED (endline)
# =============================================================================================
# If the invitation affected quality of life only through attending, invited people who did not
# attend should look like the never invited. This group SELECTED ITSELF out of attending, so the
# comparison is not experimental and is only suggestive. Variant 2 also drops the invited people
# with 1 to 5 days (they had some exposure to the centre).
for (v in names(OUT)) {
  s <- ancova_sample(dat, v, 2) %>% left_join(dat %>% select(id, days2), by = "id")
  h1 <- s %>% filter(!(z == 1 & d == 1))                    # invited with d = 0, plus never invited
  h2 <- h1 %>% filter(!(z == 1 & days2 >= 1))               # invited who never set foot in the centre
  r1 <- pull_est(lm(y ~ y0 + z, h1), "z"); r2 <- pull_est(lm(y ~ y0 + z, h2), "z")
  add_est("h", paste0("excl_nonatt_vs_never_", v),
          sprintf("%s endline: invited not attending (n %d) minus never invited (n %d), ANCOVA, HC2", OUT[[v]],
                  sum(h1$z == 1), sum(h1$z == 0)), r1, note = "SELF SELECTED comparison: suggestive only")
  add_est("h", paste0("excl_zero_days_vs_never_", v),
          sprintf("%s endline: invited with zero days (n %d) minus never invited (n %d), ANCOVA, HC2", OUT[[v]],
                  sum(h2$z == 1), sum(h2$z == 0)), r2, note = "SELF SELECTED comparison: suggestive only")
}

# =============================================================================================
# (i) MIDLINE ROBUSTNESS TO THE TIMING OF THE INVITATION
# =============================================================================================
# Round in which each invited person was picked: t_pm.csv, T1 to T8 = 1 in that round (as in 04).
tpm <- read_csv(paste0(RAW, "t_pm.csv"), show_col_types = FALSE, name_repair = "unique_quiet") %>%
  select(id = x03_01, all_of(paste0("T", 1:8)))
rnd <- tpm %>% pivot_longer(-id, names_to = "r", values_to = "picked") %>% filter(picked %in% 1) %>%
  group_by(id) %>% summarise(round = min(as.integer(sub("T", "", r))), .groups = "drop")
stopifnot(!anyDuplicated(rnd$id), setequal(rnd$id, dat$id[dat$z == 1]),
          identical(as.integer(table(rnd$round)), c(30L, 5L, 7L, 17L, 9L, 17L, 3L, 3L)))
# ROUND DATES ARE NOT CONFIRMED WITH IMSS, and two of our sources disagree for rounds 7 and 8.
#
# Source 1, our own files (round_dates_bounds.csv in code/2023/imss/imss/2_script/paper/revision/
# output/csv; decision_log.md, section 6). The date below is the LATEST date each round can have
# taken place according to those files:
#   1 selection file "2ene22" (year typo): 2 Jan 2023;  2 in the 23 Jan follow up snapshot: 23 Jan;
#   3 "around 20 to 25 January": 25 Jan;  4 first entries 21 Feb: 20 Feb;  5 first entries 29 March:
#   28 March;  6, 7 and 8: the date the randomisation script was last saved (aleatorizacion_pc6.R
#   26 April, pc7.R 26 July, pc8.R 11 September 2023). For rounds 7 and 8 the files give no lower
#   bound at all ("late July", "early September" are the script save dates restated).
# Source 2, the attendance register: a round cannot have taken place after the first day any of its
#   invitees attended the centre. For rounds 7 and 8 this bound is months EARLIER than the script
#   save dates: every round 7 and 8 invitee who ever attended did so before the save date of their
#   round's script, and the register's own Ronda column labels two of the three t_pm.csv round 8
#   invitees as round 7 and leaves the third blank. Either rounds 7 and 8 are misdated by the script
#   dates, or invitees attended before being drawn. Counted below; to raise with IMSS.
#
# So the people surveyed at midline before their own invitation are a RANGE, not a number:
#   file dates:     surveyed before the file date of their round                          (upper count)
#   register bound: surveyed before the earlier of the file date and the first attendance
#                   of any invitee of their round                                         (lower count)
# The register bound is itself only an upper bound on the round date, so the lower count is "at most"
# (it would be lower still if the round took place earlier). Both sets are dropped in turn below.
# Neither count should be stated as fact until IMSS confirms the dates of rounds 7 and 8.
ROUND_FILE <- tribble(
  ~round, ~file_date,
  1,      as.Date("2023-01-02"),
  2,      as.Date("2023-01-23"),
  3,      as.Date("2023-01-25"),
  4,      as.Date("2023-02-20"),
  5,      as.Date("2023-03-28"),
  6,      as.Date("2023-04-26"),
  7,      as.Date("2023-07-26"),
  8,      as.Date("2023-09-11"))
# register bound: first day attended by any invitee of the round (rounds with no entrant: none)
ROUND_REG <- dat %>% filter(z == 1) %>% left_join(rnd, by = "id") %>% group_by(round) %>%
  summarise(reg_bound = if (all(is.na(entry_date))) as.Date(NA) else min(entry_date, na.rm = TRUE),
            .groups = "drop")
ROUND_DATES <- ROUND_FILE %>% left_join(ROUND_REG, by = "round") %>%
  mutate(latest_both = pmin(file_date, reg_bound, na.rm = TRUE))   # earlier of the two sources

mid <- dat %>% filter(!is.na(whoqol_y0), !is.na(whoqol_y1)) %>%
  left_join(rnd, by = "id") %>% left_join(ROUND_DATES, by = "round") %>%
  mutate(pre_file  = z == 1 & date1 < file_date,             # before the round's file date
         pre_reg   = z == 1 & date1 < latest_both,           # before the register bound as well
         late_base = z == 1 & round %in% 1 & date0 > ROUND_FILE$file_date[1])
stopifnot(nrow(mid) == 102,
          all(mid$pre_file[mid$pre_reg]),                    # the lower set sits inside the upper one
          # rounds 1 to 6 cannot matter: every invitee's midline survey is after round 6's file date
          min(mid$date1[mid$z == 1]) > max(ROUND_FILE$file_date[ROUND_FILE$round <= 6]),
          all(mid$round[mid$pre_file] %in% 7:8),
          sum(mid$late_base) == 1, all(mid$d1[mid$pre_file] == 0))
pre_file_ids <- mid$id[mid$pre_file]; pre_reg_ids <- mid$id[mid$pre_reg]; late_ids <- mid$id[mid$late_base]
by_round <- function(flag) paste(names(table(mid$round[flag])), table(mid$round[flag]), sep = " x ", collapse = ", ")
add("i", "pre_invitation_n_file_dates",
    "Invited people surveyed at midline before their round's FILE date (upper count; rounds 7 and 8 dated by script save dates)",
    length(pre_file_ids), n = sum(mid$z == 1),
    note = sprintf("by round: %s. Depends on the script save dates, which the register contradicts; do not state as fact until IMSS confirms rounds 7 and 8", by_round(mid$pre_file)))
add("i", "pre_invitation_n_register_bound",
    "Invited people surveyed at midline before the register bound of their round (at most this many; lower count)",
    length(pre_reg_ids), n = sum(mid$z == 1),
    note = sprintf("by round: %s. Register bound = the earlier of the file date and the first attendance of any invitee of the round; 0 if the round was earlier still", by_round(mid$pre_reg)))
add("i", "pre_invitation_late_entrants_n_file_dates",
    "Of the upper count, invited late entrants (entered the centre after their midline survey)",
    sum(mid$pre_file & coalesce(mid$entry_date > mid$date1, FALSE)), n = length(pre_file_ids))
add("i", "pre_invitation_late_entrants_n_register_bound",
    "Of the lower count, invited late entrants (entered the centre after their midline survey)",
    sum(mid$pre_reg & coalesce(mid$entry_date > mid$date1, FALSE)), n = length(pre_reg_ids))
add("i", "late_baseline_n", "Invited in round 1, baseline dated after round 1 (2 January 2023, to confirm); in the midline sample",
    sum(mid$late_base), note = "not in the endline sample")

# ---- the round 7 and 8 dating anomaly, in counts (all 120, not only the midline sample) ----------
r78 <- dat %>% filter(z == 1) %>% left_join(rnd, by = "id") %>% left_join(ROUND_FILE, by = "round") %>%
  filter(round %in% 7:8)
# days attended strictly before the round's script save date, per invitee (kept internal)
r78$days_before_file <- vapply(seq_len(nrow(r78)), function(i)
  sum(visits$id == r78$id[i] & visits$day < r78$file_date[i]), numeric(1))
stopifnot(nrow(r78) == 6)
add("i", "r78_entrants_before_file_date",
    "Round 7 and 8 invitees who ever attended: number whose FIRST attendance precedes their round's script save date",
    sum(!is.na(r78$entry_date) & r78$entry_date < r78$file_date), n = sum(!is.na(r78$entry_date)),
    note = sprintf("n = round 7 and 8 invitees with any attendance (%d of %d invitees)", sum(!is.na(r78$entry_date)), nrow(r78)))
add("i", "r78_more_than_5_days_before_file_date",
    "Round 7 and 8 invitees with more than five days attended before their round's script save date",
    sum(r78$days_before_file > 5), n = nrow(r78))
reg_cmp <- dat %>% filter(z == 1) %>% left_join(rnd, by = "id") %>% left_join(reg_round, by = "id")
stopifnot(nrow(reg_cmp) == 91,                                   # the note below says rounds 1 to 7 agree
          all(reg_cmp$round[is.na(reg_cmp$ronda) | reg_cmp$ronda != reg_cmp$round] == 8))
add("i", "register_round_disagrees_n",
    "Invited people whose round in the register's Ronda column differs from t_pm.csv (blank counted as differing)",
    sum(is.na(reg_cmp$ronda) | reg_cmp$ronda != reg_cmp$round), n = nrow(reg_cmp),
    note = sprintf("all are t_pm round 8 invitees (%d): %d labelled round 7 in the register, %d blank; rounds 1 to 7 agree",
                   sum(reg_cmp$round == 8), sum(reg_cmp$round == 8 & reg_cmp$ronda %in% 7),
                   sum(reg_cmp$round == 8 & is.na(reg_cmp$ronda))))

# ---- midline estimates without each set ----------------------------------------------------
excl <- list(pre_file = pre_file_ids, pre_reg = pre_reg_ids, late = late_ids,
             both_file = union(pre_file_ids, late_ids), both_reg = union(pre_reg_ids, late_ids))
excl_lab <- c(
  pre_file  = sprintf("without the %d surveyed before their round's file date, the upper count", length(pre_file_ids)),
  pre_reg   = sprintf("without the %d surveyed before their round's register bound, the lower count", length(pre_reg_ids)),
  late      = "without the round 1 invitee with a late baseline",
  both_file = sprintf("without the upper count and the late baseline (%d people)", length(excl$both_file)),
  both_reg  = sprintf("without the lower count and the late baseline (%d people)", length(excl$both_reg)))
for (v in names(OUT)) for (k in names(excl)) {
  s <- ancova_sample(dat, v, 1) %>% filter(!id %in% excl[[k]])
  add_est("i", sprintf("mid_%s_itt_%s", k, v), sprintf("%s, midline ITT %s (HC2)", OUT[[v]], excl_lab[[k]]),
          pull_est(lm(y ~ y0 + z, s), "z"))
  add_est("i", sprintf("mid_%s_iv_%s", k, v), sprintf("%s, midline IV %s (HC2)", OUT[[v]], excl_lab[[k]]),
          pull_est(ivreg(y ~ d + y0 | z + y0, data = s), "d"))
}

# ---- both waves without every round 7 and 8 invitee ------------------------------------------
# Removes any doubt about the dates of rounds 7 and 8. The comparison stays randomised: the invitees
# of rounds 1 to 6 are a uniform draw of 85 of the 120, and the never invited are a uniform draw of
# 29 of the 35 left over, so both groups are random subsets of the same 120.
r78_ids <- r78$id
for (v in names(OUT)) for (wv in 1:2) {
  s <- ancova_sample(dat, v, wv) %>% filter(!id %in% r78_ids)
  nd <- sum(ancova_sample(dat, v, wv)$id %in% r78_ids)
  add_est("i", sprintf("no_r78_%s_itt_%s", WAVE[wv], v),
          sprintf("%s, %s ITT without the round 7 and 8 invitees (%d in this sample; HC2)", OUT[[v]], WAVE[wv], nd),
          pull_est(lm(y ~ y0 + z, s), "z"), note = "invited in rounds 1 to 6 against never invited: still randomised")
  add_est("i", sprintf("no_r78_%s_iv_%s", WAVE[wv], v),
          sprintf("%s, %s IV without the round 7 and 8 invitees (%d in this sample; HC2)", OUT[[v]], WAVE[wv], nd),
          pull_est(ivreg(y ~ d + y0 | z + y0, data = s), "d"), note = "invited in rounds 1 to 6 against never invited: still randomised")
}

# =============================================================================================
# (j) DE CHAISEMARTIN AND BEHAGHEL (2020) WAITING LIST CORRECTION
# =============================================================================================
# With a waiting list, offers continue until the seats are filled, so the last person offered is a
# taker by construction and takers are over represented among the offered. Their fix (single list):
# weight every invited attender by (S minus 1)/S, everyone else by 1, where S is the number of takers,
# and run the same ITT and 2SLS with those weights. It is a weight on observations, not a rescaling of
# the outcome (the old pipeline rescaled the outcome, which is wrong). Two versions, each weighting
# exactly the set of takers its S counts:
#   S list   = everyone who attended more than five days in the whole register (the design quantity);
#              the weight goes on every invited person in the sample who is a taker in that sense
#              (ever_att = 1), whether or not they had started by that wave's survey
#   S sample = those attending by that wave's survey in that wave's estimation sample (d = 1); the
#              weight goes on them
# In both, the IV's endogenous variable is still d (attending by that survey), as in 03.
S_list <- sum(dat$ever_att)
add("j", "dcb_S_list", "Takers in the whole waiting list: attended more than five days anywhere in the register (all 120)", S_list)
for (v in names(OUT)) for (wv in 1:2) {
  s <- ancova_sample(dat, v, wv) %>% left_join(dat %>% select(id, ever_att), by = "id")
  S_smp <- sum(s$d == 1)
  stopifnot(all(s$ever_att[s$d == 1] == 1))                     # attending by a survey implies ever
  for (Sv in c("list", "sample")) {
    S     <- if (Sv == "list") S_list else S_smp
    taker <- if (Sv == "list") s$ever_att else s$d
    s$w <- 1 - s$z * taker / S                                  # (S - 1)/S on invited takers
    ri <- pull_est(lm(y ~ y0 + z, s, weights = w), "z")
    rv <- pull_est(ivreg(y ~ d + y0 | z + y0, data = s, weights = w), "d")
    nt <- sum(s$z * taker)
    add_est("j", sprintf("dcb_%s_itt_%s_S%s", WAVE[wv], v, Sv),
            sprintf("%s, %s ITT with the waiting list weight, S = %d (%s)", OUT[[v]], WAVE[wv], S, Sv), ri,
            note = sprintf("weight (S minus 1)/S on %d invited takers in this sample (%s)", nt,
                           ifelse(Sv == "list", "more than five days anywhere in the register", "more than five days by this survey")))
    add_est("j", sprintf("dcb_%s_iv_%s_S%s", WAVE[wv], v, Sv),
            sprintf("%s, %s IV with the waiting list weight, S = %d (%s)", OUT[[v]], WAVE[wv], S, Sv), rv,
            note = sprintf("weight (S minus 1)/S on %d invited takers in this sample (%s)", nt,
                           ifelse(Sv == "list", "more than five days anywhere in the register", "more than five days by this survey")))
  }
}

# =============================================================================================
# (k) ATTRITION BY GROUP AND BY ARM
# =============================================================================================
# Lost at a wave = not surveyed at that wave (everyone has a baseline, so this is exactly "not in that
# wave's ANCOVA sample"). Groups use attendance over the WHOLE register, because people lost at a
# wave have no survey date to count attendance up to. So "attended" HERE means more than five days
# anywhere in the register (18 January 2023 to 16 February 2024), NOT "more than five days by that
# wave's survey" as everywhere else in the paper: it includes 3 invitees who began after their
# midline survey and people with no follow up survey. Every label below says so.
REG_DEF <- "more than five days anywhere in the register, 18 January 2023 to 16 February 2024, not by the survey"
GRP_LAB <- c("attended" = "attended (more than five days anywhere in the register)",
             "invited not attending" = "invited not attending (five days or fewer in the whole register)",
             "never invited" = "never invited")
grp <- dat %>% mutate(group = case_when(z == 0 ~ "never invited",
                                        ever_att == 1 ~ "attended",
                                        TRUE ~ "invited not attending"),
                      arm = if_else(z == 1, "invited", "never invited"))
for (wv in 1:2) {
  g <- grp %>% mutate(lost = is.na(.data[[paste0("whoqol_y", wv)]]))
  for (gg in c("attended", "invited not attending", "never invited")) {
    x <- g %>% filter(group == gg)
    add("k", sprintf("attrition_%s_%s", WAVE[wv], gsub(" ", "_", gg)),
        sprintf("Per cent not surveyed at %s, %s (value = per cent; n = group size)", WAVE[wv], GRP_LAB[[gg]]),
        100 * mean(x$lost), n = nrow(x),
        note = sprintf("%d of %d; groups by attendance %s", sum(x$lost), nrow(x), REG_DEF))
  }
  add("k", sprintf("attrition_%s_group_fisher", WAVE[wv]),
      sprintf("Attrition at %s across the three groups (attended = more than five days anywhere in the register): Fisher exact p", WAVE[wv]), NA,
      p = fisher.test(table(g$group, g$lost))$p.value, n = nrow(g), note = paste("groups by attendance", REG_DEF))
  for (aa in c("invited", "never invited")) {
    x <- g %>% filter(arm == aa)
    add("k", sprintf("attrition_%s_arm_%s", WAVE[wv], gsub(" ", "_", aa)),
        sprintf("Per cent not surveyed at %s, %s arm (value = per cent; n = arm size)", WAVE[wv], aa),
        100 * mean(x$lost), n = nrow(x), note = sprintf("%d of %d", sum(x$lost), nrow(x)))
  }
  add("k", sprintf("attrition_%s_arm_fisher", WAVE[wv]),
      sprintf("Attrition at %s by arm: Fisher exact p", WAVE[wv]), NA,
      p = fisher.test(table(g$arm, g$lost))$p.value, n = nrow(g))
}
add("k", "attended_whole_register_n", "Invited people who attended more than five days anywhere in the register (the attended group of this section)",
    sum(grp$group == "attended"), n = sum(grp$z == 1),
    note = sprintf("of them, %d had no midline survey and %d no endline survey; %d of those surveyed at midline had not yet attended more than five days by it",
                   sum(grp$group == "attended" & is.na(grp$whoqol_y1)), sum(grp$group == "attended" & is.na(grp$whoqol_y2)),
                   sum(grp$group == "attended" & !is.na(grp$whoqol_y1) & grp$d1 %in% 0)))

# =============================================================================================
# (l) CONTEXT: SUBSTITUTES, USE OF THE CENTRE, COVARIATE SENSITIVITY, MISSING ITEMS
# =============================================================================================
# (1) Substitutes for the centre. Caregiver item 06.05 asks whether the older adult currently
#     attends a day centre ("casa de dia"), club or community space where they interact with
#     others. Source: the raw questionnaire file the data core reads (pm_limpieza1.csv, item
#     x06_05_pc, one answer per person and wave, keyed by study folio). At midline the answers
#     separate the IMSS centre from other places ("Si, el Centro de Dia del IMSS", "Si, otro"); at
#     baseline and endline they are only yes or no. Nobody who was not invited attended the IMSS
#     centre, so a yes from a never invited person refers to another place. At endline a yes from
#     an invited person who did not attend is ambiguous (it may refer to a few days at the IMSS
#     centre), so that row is descriptive only. Counts among the people surveyed at each wave.
# (2) Use of the centre: days of operation in the attendance register (one column per day) and
#     people attending per day ("x" or "X" marks, as above), over every register row.
# (3) Sensitivity of the endline composite to the sociodemographic covariates: column (c) of
#     Table A.8 (05_robustness.R) minus one covariate at a time, ITT and IV, HC2; the full set is
#     first checked against 05's tableA_covariates.csv.
# (4) WHOQOL-BREF items missing in the surveys held, as counted by 00_data.R while scoring.

# ---- (1) substitutes -----------------------------------------------------------------------
q0605 <- .read_raw_csv(RAW, "pm_limpieza1.csv") %>%
  filter(name == "x06_05_pc") %>%
  transmute(id = x03_01, a0 = .ascii(t0), a1 = .ascii(t1), a2 = .ascii(t2))
stopifnot(nrow(q0605) == 120, setequal(q0605$id, dat$id),
          all(na.omit(c(q0605$a0, q0605$a2)) %in% c("No", "Si")),
          all(na.omit(q0605$a1) %in% c("No", "Si, el Centro de Dia del IMSS", "Si, otro")))
cx <- dat %>% select(id, z, d1, d2, whoqol_y1, whoqol_y2) %>% left_join(q0605, by = "id")
stopifnot(!any(cx$a1[cx$z == 0] %in% "Si, el Centro de Dia del IMSS"))   # one sided, as everywhere
cnt <- function(x, yes) c(yes = sum(x %in% yes), answered = sum(!is.na(x)), people = length(x))
SUB_LAB <- "attends a day centre, club or community space"
m1 <- cx %>% filter(!is.na(whoqol_y1)); e2 <- cx %>% filter(!is.na(whoqol_y2))
sub_rows <- list(
  list("subst_baseline_all",   "Baseline, all 120 people",                      cnt(cx$a0, "Si")),
  list("subst_baseline_never", "Baseline, never invited",                       cnt(cx$a0[cx$z == 0], "Si")),
  list("subst_midline_never",  "Midline sample, never invited: another place ('Si, otro')",
       cnt(m1$a1[m1$z == 0], "Si, otro")),
  list("subst_midline_invited_not_attending", "Midline sample, invited not attending: another place ('Si, otro')",
       cnt(m1$a1[m1$z == 1 & m1$d1 == 0], "Si, otro")),
  list("subst_endline_never",  "Endline sample, never invited: yes (another place)",
       cnt(e2$a2[e2$z == 0], "Si")),
  list("subst_endline_invited_not_attending", "Endline sample, invited not attending: yes (place not stated; descriptive only)",
       cnt(e2$a2[e2$z == 1 & e2$d2 == 0], "Si")))
for (r in sub_rows) {
  k <- r[[3]]
  add("l", r[[1]], sprintf("Caregiver item 06.05, %s: %s (value = yes; n = answered)", r[[2]], SUB_LAB),
      k[["yes"]], n = k[["answered"]], note = sprintf("%d of the %d people in the group answered", k[["answered"]], k[["people"]]))
}

# ---- (2) use of the centre from the register -------------------------------------------------
en_date <- function(d) sprintf("%d %s %s", as.integer(format(d, "%d")), month.name[as.integer(format(d, "%m"))],
                               format(d, "%Y"))                       # English month names in any locale
daily <- marks %>%
  pivot_longer(-id, names_to = "day", values_to = "mark", values_transform = list(mark = as.character)) %>%
  group_by(day) %>% summarise(n = sum(tolower(mark) %in% "x"), .groups = "drop") %>%
  mutate(day = as.Date(day))
stopifnot(nrow(daily) == 223)
add("l", "register_days_operation", "Days of operation in the centre's attendance register (one column per day)",
    nrow(daily), n = nrow(daily),
    note = sprintf("first day %s, last day %s; days with at least one person attending: %d",
                   en_date(min(daily$day)), en_date(max(daily$day)), sum(daily$n > 0)))
add("l", "register_people_per_day", "People attending per day of operation, register: mean (range in lo to hi)",
    mean(daily$n), lo = min(daily$n), hi = max(daily$n), n = nrow(daily),
    note = sprintf("median %.0f; every register row counted, study participant or not", median(daily$n)))

# ---- (3) leave one sociodemographic covariate out, endline composite -----------------------
SOC <- c("age0", "female", "indigenous", "hh_size", "educ", "child_carer")   # column (c) of Table A.8
s_end <- ancova_sample(dat, "whoqol", 2) %>% left_join(dat %>% select(id, all_of(SOC)), by = "id")
fit_set <- function(xs) {
  s <- s_end %>% drop_na(all_of(xs)); rhs <- paste(c("y0", xs), collapse = " + ")
  bind_rows(pull_est(lm(as.formula(paste("y ~ z +", rhs)), s), "z") %>% mutate(model = "ITT"),
            pull_est(ivreg(as.formula(paste("y ~ d +", rhs, "| z +", rhs)), data = s), "d") %>%
              mutate(model = "IV"))
}
full_c <- fit_set(SOC)
cov05 <- read_csv(paste0(OUT_CSV, "tableA_covariates.csv"), show_col_types = FALSE) %>%
  filter(outcome == "WHOQOL-BREF", wave == 2, column == "c") %>% select(model, est_05 = est, p_05 = p)
chk_c <- full_c %>% inner_join(cov05, by = "model")
stopifnot(nrow(chk_c) == 2, max(abs(chk_c$est - chk_c$est_05)) < 1e-8, max(abs(chk_c$p - chk_c$p_05)) < 1e-8)
loo <- map_dfr(SOC, function(dr) fit_set(setdiff(SOC, dr)) %>% mutate(dropped = dr))
for (mm in c("ITT", "IV")) {
  x <- loo %>% filter(model == mm); f <- full_c %>% filter(model == mm)
  add("l", sprintf("loo_soc_endline_whoqol_%s", tolower(mm)),
      sprintf(paste("WHOQOL-BREF endline %s with the sociodemographic covariates of Table A.8, column (c),",
                    "dropping one covariate at a time: range of the six p values in lo to hi (value = number of models)"), mm),
      nrow(x), lo = min(x$p), hi = max(x$p), n = max(x$n),
      note = sprintf("estimates %.2f to %.2f points; with all six covariates %.2f points, p %.3f",
                     min(x$est), max(x$est), f$est, f$p))
}

# ---- (4) missing WHOQOL-BREF items -----------------------------------------------------------
n_surveys <- sum(!is.na(dat$whoqol_y0)) + sum(!is.na(dat$whoqol_y1)) + sum(!is.na(dat$whoqol_y2))
add("l", "whoqol_missing_items", "WHOQOL-BREF items missing in the surveys held (value = items missing; n = surveys)",
    n_missing_items, n = n_surveys, note = "26 items per survey; the WHO missing item rules of 00_data.R never bind")

# =============================================================================================
# WRITE AND PRINT
# =============================================================================================
tn <- bind_rows(rows) %>% mutate(across(c(value, se, p, lo, hi), ~ round(.x, 4)))
stopifnot(!anyDuplicated(tn$item))
write_csv(tn, paste0(OUT_CSV, "text_numbers.csv"), na = "")
cat(sprintf("\nwritten %stext_numbers.csv (%d rows)\n\n", OUT_CSV, nrow(tn)))
print(as.data.frame(tn %>% select(section, item, value, se, p, lo, hi, n)), row.names = FALSE)
