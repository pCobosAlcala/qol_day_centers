# 05_robustness.R  (revision, October 2026; data core and 0 to 100 scale)
#
# Two annex tables that check the main results against other modelling choices, plus the numbers
# the text and the response letter quote.
#
#   (1) Covariates (reviewer 2, comment 15: "introduce additional controls sequentially").
#       Same per wave ANCOVA as the main table, with five covariate columns:
#         (a) none, the main specification (must reproduce 03's main_results.csv exactly)
#         (b) baseline functionality: Katz and Lawton
#         (c) sociodemographics: age at baseline, sex, indigenous, household size, education,
#             main caregiver is a son or daughter (the CORRECTED caregiver flag of 00_data.R; the
#             September version used the old control.csv flag, which counted only "Padre o madre")
#         (d) both sets
#         (e) Lin (2013): invitation fully interacted with the centred baseline outcome, Katz and
#             Lawton
#       Every column keeps the outcome's own baseline value (it is an ANCOVA). Robust HC2 errors.
#       Table: IV (the primary estimate) and ITT for the composite and the Environmental domain,
#       ITT for the other four parts. The csv holds IV and ITT for all six outcomes.
#       People missing a covariate drop out of that column, so n is reported per column.
#
#   (2) Alternative specifications (reviewer 2, comments 3 and 6). Composite and Environmental:
#         (a) ANCOVA per wave, the main specification (HC2; must reproduce 03 exactly)
#         (b) ANCOVA pooling both follow-up waves: one average effect
#         (c) difference in differences by wave with individual fixed effects
#         (d) the same with a random intercept per person (ITT only), the submission's family
#         (e) AS TREATED: attended versus did not attend, the same per wave ANCOVA by OLS with HC2
#             errors (reviewer 2, comment 3). Not randomised: subject to self selection.
#         (f) csv only: the submitted model (linear time, DC = invitation x post, random intercept),
#             on the submitted data (82 people, 27 items, 1 to 5 scale) and on the current data
#             (120 people, 26 items, 0 to 100), so the response letter can quote it.
#       Errors for (b) and (c): clustered by person, CR2 (Bell and McCaffrey 2002), computed with
#       sandwich::vcovCL(type = "HC2"); with the default cadjust = TRUE sandwich cancels the
#       G/(G-1) factor for HC2, so this is the textbook CR2. Tests and intervals use a t
#       distribution with G - 1 degrees of freedom. For (c) the fixed effects are absorbed by
#       demeaning within person before the regression; the CR2 adjustment then uses the absorbed
#       design, which equals the full dummy model with a generalised inverse (Pustejovsky and
#       Tipton 2018; checked numerically on these data when this script was written).
#       Samples: (a) and (e) everyone observed at baseline and at that wave (102 and 85);
#       (b) and (c) everyone observed at baseline and at least one follow-up (people seen only at
#       baseline have no follow-up row in (b) and nothing to difference in (c)); (d) and (f) on
#       the current data: everyone with a baseline (120).
#
# Data: build_analysis_data() in 00_data.R (scores on the WHO 0 to 100 scale, wave specific
# attendance, all 120 randomised people). The September version read pm_who.csv, pm.csv and
# control.csv on the 1 to 5 scale with 99 people at midline; copies of that script and its outputs
# are in the _previous subfolders, suffix _before_00data.
#
# Outputs
#   tables: tableA_covariates.png, tableA_specifications.png                         (OUT_TAB)
#   csv   : tableA_covariates.csv, tableA_specifications.csv, text_numbers_robustness.csv (OUT_CSV)
# Aggregates only: no names, birth dates, exact dates or individual rows are written or printed.
# Run from the imss_paper folder:
#   "C:/Program Files/R/R-4.6.1/bin/Rscript.exe" 3_code/scripts/05_robustness.R

rm(list = ls())
suppressPackageStartupMessages({
  library(tidyverse); library(AER); library(sandwich); library(lmtest)
  library(gridExtra); library(gtable); library(grid)
})
# fixest (only to check the within transformation) and lmerTest (random effects) are called with
# :: so that their exported names do not mask anything here
stopifnot(requireNamespace("fixest", quietly = TRUE), requireNamespace("lmerTest", quietly = TRUE))
options(scipen = 999, width = 200)

source("3_code/scripts/00_paths.R")
source("3_code/scripts/00_data.R")

# =============================================================================================
# 1. DATA: one row per randomised person (120), scores on the WHO 0 to 100 scale
# =============================================================================================
dat <- build_analysis_data()                 # checks against the old files run inside
outcomes <- OUTCOMES                         # whoqol, general, physical, psychological, social, environmental

FUNC <- c("katz0", "lawton0")                                              # set (b)
SOC  <- c("age0", "female", "indigenous", "hh_size", "educ", "child_carer") # set (c)
covs <- dat %>% select(id, all_of(c(FUNC, SOC)))

# denominators for per cents: the control group's baseline mean and the baseline SD (both arms)
# in each per wave ANCOVA sample, exactly as in 03's magnitudes table, so that every per cent in
# this script is comparable with the main text
denom <- map_dfr(names(outcomes), function(v) map_dfr(1:2, function(wv) {
  s <- ancova_sample(dat, v, wv)
  tibble(outcome = outcomes[[v]], wave = wv,
         control_baseline = mean(s$y0[s$z == 0]), sd_baseline = sd(s$y0))
}))

# coverage of the covariates in each ANCOVA sample (aggregates only)
cat("\ncovariates missing in each ANCOVA sample of the composite:\n")
for (wv in 1:2) {
  s <- ancova_sample(dat, "whoqol", wv) %>% left_join(covs, by = "id")
  cat(sprintf("  %s (n = %d): ", c("midline", "endline")[wv], nrow(s)))
  miss <- colSums(is.na(s[, c(FUNC, SOC)]))
  cat(paste(names(miss), miss, sep = " ", collapse = ", "), "\n")
}
cat(sprintf("mean age at baseline, all 120: %.1f\n", mean(dat$age0)))

# =============================================================================================
# 2. ESTIMATION HELPERS
# =============================================================================================
# HC2, one row per person (as in 03_main_results.R): HC1 only if HC2 fails, and say so
rob <- function(m) {
  tryCatch(vcovHC(m, type = "HC2"),
           error = function(e) { warning("HC2 failed, HC1 used: ", conditionMessage(e))
                                 vcovHC(m, type = "HC1") })
}
pull_est <- function(m, term) {
  V  <- rob(m); ct <- coeftest(m, vcov. = V); ci <- coefci(m, vcov. = V, level = 0.95)
  tibble(est = ct[term, 1], se = ct[term, 2], p = ct[term, 4], lo = ci[term, 1], hi = ci[term, 2])
}

# CR2, clustered by person. For an ivreg object sandwich warns that clustered HC2 is meant for
# (generalised) linear models; it then builds the adjustment from the projected regressors, the
# standard 2SLS analogue (checked against a direct computation). Only that warning is muffled.
cr2 <- function(m, cl) {
  withCallingHandlers(vcovCL(m, cluster = cl, type = "HC2"),
    warning = function(w) {
      if (grepl("clustered HC2/HC3 are only applicable", conditionMessage(w), fixed = TRUE))
        invokeRestart("muffleWarning")
    })
}
pull_cl <- function(m, term, cl) {
  V  <- cr2(m, cl); G <- length(unique(cl))
  ct <- coeftest(m, vcov. = V, df = G - 1)
  ci <- coefci(m, parm = term, vcov. = V, level = 0.95, df = G - 1)
  tibble(est = ct[term, 1], se = ct[term, 2], p = ct[term, 4], lo = ci[term, 1], hi = ci[term, 2])
}

# =============================================================================================
# 3. COVARIATE TABLE (reviewer 2, comment 15)
# =============================================================================================
# header of each column: letter on the first line, content on the second
sets <- list(
  a = list(lab = "(a)\nNone (main)",            x = character(0)),
  b = list(lab = "(b)\nKatz and Lawton",        x = FUNC),
  c = list(lab = "(c)\nSociodemographics",      x = SOC),
  d = list(lab = "(d)\nBoth (b) and (c)",       x = c(FUNC, SOC)),
  e = list(lab = "(e)\nLin (2013), interacted", x = FUNC))   # (e) also interacts y0, below

cov_res <- list()
for (v in names(outcomes)) {
  for (wv in 1:2) {
    s0 <- ancova_sample(dat, v, wv) %>% left_join(covs, by = "id")
    for (k in names(sets)) {
      xs <- sets[[k]]$x
      s  <- s0 %>% drop_na(all_of(xs))     # baseline, this wave, and every covariate of the column
      if (k != "e") {
        rhs <- paste(c("y0", xs), collapse = " + ")
        itt <- lm(as.formula(paste("y ~ z +", rhs)), s)
        iv  <- ivreg(as.formula(paste("y ~ d +", rhs, "| z +", rhs)), data = s)
        fs  <- lm(as.formula(paste("d ~ z +", rhs)), s)
      } else {
        # Lin (2013): centre the covariates at the sample mean and interact them with the
        # invitation; the coefficient on z is then the average effect. The IV keeps the
        # z x covariate terms as included regressors and instruments attendance with z alone,
        # so it equals the Lin adjusted ITT divided by the Lin adjusted first stage.
        s <- s %>% mutate(y0c = y0 - mean(y0), katzc = katz0 - mean(katz0),
                          lawtonc = lawton0 - mean(lawton0))
        lc  <- "y0c + katzc + lawtonc + z:y0c + z:katzc + z:lawtonc"
        itt <- lm(as.formula(paste("y ~ z +", lc)), s)
        iv  <- ivreg(as.formula(paste("y ~ d +", lc, "| z +", lc)), data = s)
        fs  <- lm(as.formula(paste("d ~ z +", lc)), s)
      }
      fz <- pull_est(fs, "z")
      # the no covariate ANCOVA on the same people, to separate the covariates' effect from the
      # effect of losing people who miss a covariate
      same_itt <- pull_est(lm(y ~ z + y0, s), "z")
      same_iv  <- pull_est(ivreg(y ~ d + y0 | z + y0, data = s), "d")
      cov_res[[paste(v, wv, k)]] <- bind_rows(
        pull_est(itt, "z") %>% mutate(model = "ITT", est_nocov_same_sample = same_itt$est),
        pull_est(iv,  "d") %>% mutate(model = "IV",  est_nocov_same_sample = same_iv$est)) %>%
        mutate(outcome = outcomes[[v]], wave = wv, column = k, n = nrow(s),
               fsF = (fz$est / fz$se)^2, .before = 1)
    }
  }
}
cov_res <- bind_rows(cov_res) %>%
  mutate(column_label = gsub("\n", " ", map_chr(column, ~ sets[[.x]]$lab)), .after = column) %>%
  mutate(fsF = if_else(model == "IV", fsF, NA_real_)) %>%
  left_join(denom, by = c("outcome", "wave")) %>%
  mutate(pct_of_control_baseline = 100 * est / control_baseline, sd_units = est / sd_baseline)

# guard: column (a) must reproduce 03's main_results.csv (same samples, same models)
main_csv <- paste0(OUT_CSV, "main_results.csv")
if (!file.exists(main_csv)) stop("run 03_main_results.R first: ", main_csv, " is missing")
main <- read_csv(main_csv, show_col_types = FALSE)
cmp_main <- cov_res %>% filter(column == "a") %>%
  inner_join(main, by = c("outcome", "wave", "model"), suffix = c("", "_03"))
stopifnot(nrow(cmp_main) == 24)
gap_main <- with(cmp_main, max(abs(c(est - est_03, se - se_03, p - p_03, n - n_03))))
if (gap_main > 1e-8) stop("column (a) does not reproduce main_results.csv (largest gap ", gap_main, ")")
cat(sprintf("\ncolumn (a) reproduces main_results.csv for 6 outcomes x 2 waves x ITT and IV (largest gap %.1e)\n",
            gap_main))

write_csv(cov_res, paste0(OUT_CSV, "tableA_covariates.csv"))

cat("\n=== covariates: estimate (HC2 se) p, n ===\n")
for (m in c("IV", "ITT")) {
  cat("--", m, "\n")
  print(cov_res %>% filter(model == m) %>%
          mutate(cellv = sprintf("%.2f (%.2f) p=%.3f n=%d", est, se, p, n)) %>%
          select(outcome, wave, column, cellv) %>%
          pivot_wider(names_from = column, values_from = cellv) %>% as.data.frame(), row.names = FALSE)
}

# =============================================================================================
# 4. ALTERNATIVE SPECIFICATIONS (reviewer 2, comments 3 and 6)
# =============================================================================================
spec_out <- c(whoqol = "WHOQOL-BREF", environmental = "Environmental")
spec <- list(); submitted_current <- list(); pooled_extra <- list()

# long panel of one outcome: every wave observed, everyone with a baseline (all 120). Attendance
# is wave specific: 0 at baseline (nobody had attended), d1 at midline, d2 at endline.
long_for <- function(v) {
  dat %>%
    transmute(id, z, d_0 = 0L, d_1 = d1, d_2 = d2,
              y_0 = .data[[paste0(v, "_y0")]], y_1 = .data[[paste0(v, "_y1")]],
              y_2 = .data[[paste0(v, "_y2")]]) %>%
    pivot_longer(-c(id, z), names_to = c(".value", "wave"), names_sep = "_") %>%
    filter(!is.na(y)) %>%
    mutate(wave = as.integer(wave), mid = as.integer(wave == 1), end = as.integer(wave == 2),
           zm = z * mid, ze = z * end, dm = d * mid, de = d * end)
}

for (v in names(spec_out)) {
  o <- spec_out[[v]]

  # ---- (a) ANCOVA per wave, main (HC2), and (e) as treated OLS on the same samples ------------
  for (wv in 1:2) {
    s   <- ancova_sample(dat, v, wv)
    itt <- lm(y ~ y0 + z, s)
    iv  <- ivreg(y ~ d + y0 | z + y0, data = s)
    at  <- lm(y ~ y0 + d, s)          # as treated: attenders against everyone else, invited or not
    eff <- c("Midline", "Endline")[wv]
    spec[[paste(v, "a", wv)]] <- bind_rows(
      pull_est(itt, "z") %>% mutate(model = "ITT"), pull_est(iv, "d") %>% mutate(model = "IV")) %>%
      mutate(outcome = o, spec = "a", effect = eff, errors = "HC2",
             n_rows = nrow(s), n_people = n_distinct(s$id))
    spec[[paste(v, "e", wv)]] <- pull_est(at, "d") %>%
      mutate(model = "As treated OLS", outcome = o, spec = "e", effect = eff, errors = "HC2",
             n_rows = nrow(s), n_people = n_distinct(s$id),
             n_attending = sum(s$d), n_not_attending = sum(1 - s$d))
  }

  # ---- (b) ANCOVA pooling both follow-up waves ------------------------------------------------
  # one row per person and follow-up wave (the two per wave ANCOVA samples stacked), a wave
  # dummy lets the level differ by wave; the coefficient on the invitation (or on attendance at
  # that wave, instrumented by the invitation) is the average effect over the two waves.
  # People seen only at baseline have no row here.
  st <- bind_rows(ancova_sample(dat, v, 1) %>% mutate(endline = 0),
                  ancova_sample(dat, v, 2) %>% mutate(endline = 1))
  itt <- lm(y ~ z + y0 + endline, st)
  iv  <- ivreg(y ~ d + y0 + endline | z + y0 + endline, data = st)
  fsb <- pull_cl(lm(d ~ z + y0 + endline, st), "z", st$id)
  spec[[paste(v, "b")]] <- bind_rows(
    pull_cl(itt, "z", st$id) %>% mutate(model = "ITT"),
    pull_cl(iv,  "d", st$id) %>% mutate(model = "IV")) %>%
    mutate(outcome = o, spec = "b", effect = "Average of both waves", errors = "CR2, clustered by person",
           n_rows = nrow(st), n_people = n_distinct(st$id))
  # denominators for the pooled per cents: each person once, in the pooled sample
  ppl <- st %>% distinct(id, z, y0)
  stopifnot(!anyDuplicated(ppl$id))
  pooled_extra[[v]] <- tibble(outcome = o, control_baseline = mean(ppl$y0[ppl$z == 0]),
                              sd_baseline = sd(ppl$y0), first_stage_F = (fsb$est / fsb$se)^2,
                              first_stage = fsb$est)

  # ---- (c) difference in differences by wave, individual fixed effects -----------------------
  # all three waves for everyone with a baseline and at least one follow-up (people seen only at
  # baseline carry no information once each person has their own intercept). The fixed effects
  # are absorbed by demeaning every variable within person; the coefficients equal those of the
  # model with one dummy per person (checked against fixest below).
  lg <- long_for(v)
  fe <- lg %>% group_by(id) %>% filter(n() > 1) %>%
    mutate(across(c(y, mid, end, zm, ze, dm, de), ~ .x - mean(.x), .names = "{.col}_w")) %>%
    ungroup()
  itt <- lm(y_w ~ mid_w + end_w + zm_w + ze_w - 1, fe)
  iv  <- ivreg(y_w ~ mid_w + end_w + dm_w + de_w - 1 | mid_w + end_w + zm_w + ze_w - 1, data = fe)
  chk_fe <- fixest::feols(y ~ mid + end + zm + ze | id, fe, notes = FALSE)
  chk_iv <- fixest::feols(y ~ mid + end | id | dm + de ~ zm + ze, fe, notes = FALSE)
  stopifnot(max(abs(coef(itt)[c("zm_w", "ze_w")] - coef(chk_fe)[c("zm", "ze")])) < 1e-8,
            max(abs(coef(iv)[c("dm_w", "de_w")] - coef(chk_iv)[c("fit_dm", "fit_de")])) < 1e-8)
  spec[[paste(v, "c")]] <- bind_rows(
    pull_cl(itt, "zm_w", fe$id) %>% mutate(model = "ITT", effect = "Midline"),
    pull_cl(itt, "ze_w", fe$id) %>% mutate(model = "ITT", effect = "Endline"),
    pull_cl(iv,  "dm_w", fe$id) %>% mutate(model = "IV",  effect = "Midline"),
    pull_cl(iv,  "de_w", fe$id) %>% mutate(model = "IV",  effect = "Endline")) %>%
    mutate(outcome = o, spec = "c", errors = "CR2, clustered by person",
           n_rows = nrow(fe), n_people = n_distinct(fe$id))

  # ---- (d) the same with a random intercept per person (ITT only) ----------------------------
  # nothing absorbs the invitation main effect now, so z stays in; zm and ze are the change in
  # the gap since baseline, as in (c). Everyone with a baseline, including people seen only at
  # baseline. Model based errors, Satterthwaite degrees of freedom (lmerTest), as in the submission.
  re <- lmerTest::lmer(y ~ mid + end + z + zm + ze + (1 | id), data = lg)
  ct <- summary(re)$coefficients[c("zm", "ze"), ]
  spec[[paste(v, "d")]] <- tibble(
    est = ct[, "Estimate"], se = ct[, "Std. Error"], p = ct[, "Pr(>|t|)"],
    lo = ct[, "Estimate"] - qt(0.975, ct[, "df"]) * ct[, "Std. Error"],
    hi = ct[, "Estimate"] + qt(0.975, ct[, "df"]) * ct[, "Std. Error"],
    model = "ITT", effect = c("Midline", "Endline"),
    outcome = o, spec = "d", errors = "Model based (random intercept), Satterthwaite",
    n_rows = nrow(lg), n_people = n_distinct(lg$id))

  submitted_current[[v]] <- lg
}

# ---- (f) the submitted model: value ~ time + DC + time:DC + (1 | person) -----------------------
# time is linear (0, 1, 2); DC = invitation x post, i.e. the submission set z to zero at baseline,
# so there is no term for the baseline gap between arms. The submission reported time:DC.
fit_submitted <- function(df, outcome, data_label, scale_label) {
  df <- df %>% mutate(time = wave, DC = z * as.integer(wave > 0))
  m  <- lmerTest::lmer(y ~ time + DC + time:DC + (1 | id), data = df)
  ct <- summary(m)$coefficients
  # implied gap between arms at each wave: DC + time x (time:DC), Satterthwaite df
  L  <- rbind(mid = c(0, 0, 1, 1), end = c(0, 0, 1, 2))
  colnames(L) <- rownames(ct)
  gp <- lmerTest::contest(m, L, joint = FALSE)
  bind_rows(
    tibble(term = c("DC (invitation x post)", "time x DC (reported in the submission)"),
           est = ct[c("DC", "time:DC"), "Estimate"], se = ct[c("DC", "time:DC"), "Std. Error"],
           p = ct[c("DC", "time:DC"), "Pr(>|t|)"], dfs = ct[c("DC", "time:DC"), "df"]),
    tibble(term = c("Implied gap at midline (DC + time:DC)", "Implied gap at endline (DC + 2 time:DC)"),
           est = gp$Estimate, se = gp$`Std. Error`, p = gp$`Pr(>|t|)`, dfs = gp$df)) %>%
    mutate(lo = est - qt(0.975, dfs) * se, hi = est + qt(0.975, dfs) * se) %>% select(-dfs) %>%
    mutate(outcome = outcome, data = data_label, scale = scale_label, spec = "f", model = "ITT",
           effect = NA_character_, errors = "Model based (random intercept), Satterthwaite",
           n_rows = nrow(df), n_people = n_distinct(df$id))
}
# the submitted analysis file (82 people with all three waves, 27 item total, 1 to 5 item means)
# (in the public repository the same file is submission/paper/db.csv)
db_file <- if (file.exists(paste0(RAW, "paper/db.csv"))) paste0(RAW, "paper/db.csv") else "submission/paper/db.csv"
db <- read_csv(db_file, show_col_types = FALSE, name_repair = "unique_quiet") %>%
  select(-1) %>%
  transmute(id = x03_01, z, wave = time, y = value, var)
sub <- bind_rows(
  fit_submitted(submitted_current[["whoqol"]], "WHOQOL-BREF",
                "Current data: everyone with a baseline (120), 26 items", "WHO 0 to 100"),
  fit_submitted(submitted_current[["environmental"]], "Environmental",
                "Current data: everyone with a baseline (120), 26 items", "WHO 0 to 100"),
  fit_submitted(db %>% filter(var == "whoqol"), "WHOQOL-BREF",
                "Submitted data: 82 people with all three waves, 27 items", "1 to 5 item mean, as submitted"),
  fit_submitted(db %>% filter(var == "Enviromental"), "Environmental",
                "Submitted data: 82 people with all three waves, 27 items", "1 to 5 item mean, as submitted"))

spec_lab <- c(a = "(a) ANCOVA by wave (main)", b = "(b) ANCOVA, both waves pooled",
              c = "(c) Difference in differences, individual fixed effects",
              d = "(d) Difference in differences, random effects",
              e = "(e) As treated OLS: attended versus did not attend (not randomised, subject to self selection)",
              f = "(f) Submitted model (reference only)")
spec <- bind_rows(spec)

# per cents: per wave effects against the control baseline of that wave's ANCOVA sample (03's
# denominators), the pooled effect against the pooled sample's control baseline
pooled_extra <- bind_rows(pooled_extra)
spec <- spec %>%
  left_join(denom %>% mutate(effect = c("Midline", "Endline")[wave]) %>% select(-wave),
            by = c("outcome", "effect")) %>%
  left_join(pooled_extra %>% transmute(outcome, effect = "Average of both waves",
                                       cb_p = control_baseline, sd_p = sd_baseline),
            by = c("outcome", "effect")) %>%
  mutate(control_baseline = coalesce(control_baseline, cb_p),
         sd_baseline = coalesce(sd_baseline, sd_p)) %>% select(-cb_p, -sd_p) %>%
  mutate(pct_of_control_baseline = 100 * est / control_baseline, sd_units = est / sd_baseline,
         data = "Current data", scale = "WHO 0 to 100", term = effect)
stopifnot(!anyNA(spec$control_baseline))

spec_csv <- bind_rows(spec, sub) %>%
  mutate(specification = spec_lab[spec]) %>%
  select(outcome, spec, specification, effect, term, model, est, se, p, lo, hi, errors,
         n_rows, n_people, n_attending, n_not_attending, control_baseline, sd_baseline,
         pct_of_control_baseline, sd_units, data, scale)
write_csv(spec_csv, paste0(OUT_CSV, "tableA_specifications.csv"))

# guard: (a) reproduces 03's main table (est, se, p and n) for the two outcomes, both waves
cmp_a <- spec %>% filter(spec == "a") %>% mutate(wave = match(effect, c("Midline", "Endline"))) %>%
  inner_join(main, by = c("outcome", "wave", "model"), suffix = c("", "_03"))
stopifnot(nrow(cmp_a) == 8,
          max(abs(with(cmp_a, c(est - est_03, se - se_03, p - p_03, n_rows - n)))) < 1e-8)
cat("\nspecification (a) reproduces main_results.csv\n")

cat("\n=== alternative specifications (points, 0 to 100) ===\n")
print(spec_csv %>% mutate(across(c(est, se, p, lo, hi, pct_of_control_baseline, sd_units), ~ round(.x, 3))) %>%
        select(outcome, spec, term, model, est, se, p, lo, hi, pct = pct_of_control_baseline, sd = sd_units,
               n_rows, n_people, scale) %>%
        as.data.frame(), row.names = FALSE)
cat("\npooled IV first stage (CR2):\n")
print(pooled_extra %>% mutate(across(where(is.double), ~ round(.x, 3))) %>% as.data.frame(), row.names = FALSE)

# =============================================================================================
# 5. NUMBERS QUOTED IN THE TEXT AND THE RESPONSE LETTER
# =============================================================================================
# cells of the covariate table: IV and ITT for the composite and Environmental, ITT for the rest
in_table <- function(df) df %>% filter(model == "ITT" | outcome %in% spec_out)

# (1) pooled ANCOVA, average effect over both waves
t_pooled <- spec %>% filter(spec == "b") %>%
  left_join(pooled_extra %>% select(outcome, first_stage_F), by = "outcome") %>%
  transmute(section = "Pooled ANCOVA, average effect over both waves (CR2 errors)",
            outcome, wave = "Average of both waves", model, column = NA_character_,
            estimate = est, se, p, lo, hi, pct_of_control_baseline, sd_units, control_baseline,
            sd_baseline, first_stage_F = if_else(model == "IV", first_stage_F, NA_real_),
            n_rows, n_people)

# (2) covariate ranges across columns (a) to (e), for every cell of the table
cov_tab <- in_table(cov_res)
main_a  <- cov_tab %>% filter(column == "a") %>% select(outcome, wave, model, est_main = est, p_main = p)
cmp     <- cov_tab %>% left_join(main_a, by = c("outcome", "wave", "model"))
t_range <- cmp %>% group_by(outcome, wave, model) %>%
  summarise(est_main = est_main[1], p_main = p_main[1], est_min = min(est), est_max = max(est),
            pct_min = min(pct_of_control_baseline), pct_max = max(pct_of_control_baseline),
            p_min = min(p), p_max = max(p), max_abs_change = max(abs(est - est_main)),
            n_min = min(n), n_max = max(n), control_baseline = control_baseline[1],
            .groups = "drop") %>%
  mutate(section = "Covariates: range across columns (a) to (e)",
         wave = c("Midline", "Endline")[wave], column = "a to e",
         pct_main = 100 * est_main / control_baseline)

# (3) which cells of the table change side of p = 0.10 against column (a)
t_cross <- cmp %>% filter(column != "a", (p < 0.10) != (p_main < 0.10)) %>%
  transmute(section = "Covariates: cell crossing p = 0.10 against column (a)",
            outcome, wave = c("Midline", "Endline")[wave], model, column,
            estimate = est, se, p, est_main, p_main, pct_of_control_baseline, control_baseline,
            est_nocov_same_sample, n_people = n,
            note = sprintf("%s: p %.3f in (a), %.3f in (%s)",
                           ifelse(p < 0.10, "becomes significant at 10 per cent",
                                  "loses significance at 10 per cent"), p_main, p, column))
n_cells <- nrow(cmp %>% filter(column != "a"))
t_cross_n <- tibble(
  section = "Covariates: number of cells crossing p = 0.10 against column (a)",
  outcome = "All cells of the table", wave = "Both waves", model = "IV and ITT", column = "b to e",
  estimate = nrow(t_cross),
  note = sprintf("%d of %d cells (%d outcome x wave x estimator rows times 4 columns)",
                 nrow(t_cross), n_cells, n_cells / 4))
t_move <- tibble(
  section = "Covariates: largest absolute change against column (a), points",
  outcome = "All cells of the table", wave = "Both waves", model = c("IV", "ITT"), column = "b to e") %>%
  mutate(estimate = map_dbl(model, ~ max(abs(cmp$est - cmp$est_main)[cmp$model == .x])))

# (4) as treated against IV and ITT, per wave
t_at <- spec %>% filter(spec %in% c("a", "e")) %>%
  transmute(section = "As treated OLS against the randomised estimates (per wave ANCOVA)",
            outcome, wave = effect, model, column = NA_character_, estimate = est, se, p, lo, hi,
            pct_of_control_baseline, sd_units, control_baseline, sd_baseline, n_rows, n_people,
            note = if_else(model == "As treated OLS",
                           sprintf("%d attending versus %d not attending; not randomised", n_attending,
                                   n_not_attending), NA_character_))

txt <- bind_rows(t_pooled, t_range, t_cross, t_cross_n, t_move, t_at) %>%
  select(section, outcome, wave, model, column, estimate, se, p, lo, hi, pct_of_control_baseline,
         sd_units, control_baseline, sd_baseline, first_stage_F, est_main, p_main, pct_main,
         est_nocov_same_sample,
         est_min, est_max, pct_min, pct_max, p_min, p_max, max_abs_change, n_min, n_max,
         n_rows, n_people, note)
write_csv(txt, paste0(OUT_CSV, "text_numbers_robustness.csv"))

cat("\n=== text numbers ===\n")
cat("-- pooled ANCOVA\n")
print(t_pooled %>% mutate(across(where(is.double), ~ round(.x, 3))) %>%
        select(outcome, model, estimate, se, p, lo, hi, pct_of_control_baseline, sd_units, control_baseline,
               first_stage_F, n_rows, n_people) %>% as.data.frame(), row.names = FALSE)
cat("-- covariate ranges, columns (a) to (e)\n")
print(t_range %>% mutate(across(where(is.double), ~ round(.x, 3))) %>%
        select(outcome, wave, model, est_main, est_min, est_max, pct_main, pct_min, pct_max,
               p_main, p_min, p_max, max_abs_change, n_min, n_max) %>% as.data.frame(), row.names = FALSE)
cat("-- cells crossing p = 0.10:", nrow(t_cross), "of", n_cells, "\n")
if (nrow(t_cross)) print(t_cross %>% mutate(across(where(is.double), ~ round(.x, 3))) %>%
                           select(outcome, wave, model, column, est_main, p_main, estimate, p,
                                  est_nocov_same_sample, n_people) %>%
                           as.data.frame(), row.names = FALSE)
cat("-- largest absolute change against (a):\n")
print(t_move %>% select(model, estimate) %>% as.data.frame(), row.names = FALSE)
cat("-- as treated against IV and ITT\n")
print(t_at %>% mutate(across(where(is.double), ~ round(.x, 3))) %>%
        select(outcome, wave, model, estimate, se, p, pct_of_control_baseline, sd_units, note) %>%
        as.data.frame(), row.names = FALSE)

# =============================================================================================
# 6. RENDERING. render() is 03_main_results.R's, plus bold_rows, which sets panel header rows in
#    bold (row numbers count the data rows of df, header excluded).
# =============================================================================================
star  <- function(p) ifelse(p < .01, "***", ifelse(p < .05, "**", ifelse(p < .1, "*", "")))
fp    <- function(p) ifelse(p < 0.001, "p < 0.001", sprintf("p = %.3f", p))
# points on the 0 to 100 scale, two decimals as in the main table, then the p value
cellp <- function(r, prefix = "") sprintf("%s%.2f%s\n(%.2f)\n%s", prefix, r$est, star(r$p), r$se, fp(r$p))

render <- function(df, spans, file, title, notes, left_cols = 1, note_width = 120, bold_rows = NULL) {
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
  # bold panel rows (tableGrob row 1 is the column header, so data row r sits at t = r + 1); bold text
  # is wider than the plain text the column widths were set for, so widen its column to fit
  for (i in which(g$layout$name == "core-fg" & g$layout$t %in% (bold_rows + 1))) {
    g$grobs[[i]]$gp <- gpar(fontsize = 10, fontface = "bold", lineheight = 1.05)
    l <- g$layout$l[i]
    g$widths[l] <- max(g$widths[l], grobWidth(g$grobs[[i]]) + unit(6.4, "mm"))
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

# ---- Table A: covariates, one panel per wave ------------------------------------------------
# rows of each panel: IV first (the primary estimate) then ITT for the composite and Environmental,
# ITT for the four other parts, in the order of the main table
show <- tribble(
  ~v,              ~m,    ~lab_out,        ~lab_est,
  "whoqol",        "IV",  "WHOQOL-BREF",   "Attending (IV)",
  "whoqol",        "ITT", "",              "Invited (ITT)",
  "general",       "ITT", "General",       "Invited (ITT)",
  "physical",      "ITT", "Physical",      "Invited (ITT)",
  "psychological", "ITT", "Psychological", "Invited (ITT)",
  "social",        "ITT", "Social",        "Invited (ITT)",
  "environmental", "IV",  "Environmental", "Attending (IV)",
  "environmental", "ITT", "",              "Invited (ITT)")
stopifnot(nrow(show) == nrow(in_table(cov_res)) / 10)   # 10 = 2 waves x 5 columns
heads <- map_chr(sets, "lab")
rng <- function(x, fmt) if (abs(max(x) - min(x)) < 1e-9) sprintf(fmt, x[1]) else
  paste(sprintf(fmt, min(x)), "to", sprintf(fmt, max(x)))
cov_panel <- function(w, head) {
  body <- pmap_dfr(show, function(v, m, lab_out, lab_est) {
    cells <- map(names(sets), function(k)
      cellp(cov_res %>% filter(outcome == outcomes[[v]], wave == w, column == k, model == m)))
    tibble(Outcome = lab_out, Estimate = lab_est, !!!set_names(cells, heads))
  })
  # first stage F on the invitation, range over the two IV outcomes; individuals per column
  fsr <- map(names(sets), function(k)
    rng(cov_res %>% filter(wave == w, column == k, model == "IV", outcome %in% spec_out) %>% pull(fsF), "%.1f"))
  nr  <- map(names(sets), function(k) rng(cov_res %>% filter(wave == w, column == k) %>% pull(n), "%d"))
  blank <- set_names(as.list(rep("", length(sets))), heads)
  bind_rows(tibble(Outcome = head, Estimate = "", !!!blank), body,
            tibble(Outcome = "First-stage F (IV)", Estimate = "", !!!set_names(fsr, heads)),
            tibble(Outcome = "Individuals", Estimate = "", !!!set_names(nr, heads)))
}
cov_tab_png <- bind_rows(cov_panel(1, "Midline, about 3 months after entry"), cov_panel(2, "Endline, about 9 months after entry"))
render(cov_tab_png,
  spans = list(list(label = "Covariate adjustment (every column also controls for the outcome's baseline value)",
                    l = 3, r = 7)),
  file  = paste0(OUT_TAB, "tableA_covariates.png"),
  title = "Effect of the day centre with added baseline covariates",
  left_cols = 1:2,
  bold_rows = which(cov_tab_png$Outcome %in% c("Midline, about 3 months after entry", "Endline, about 9 months after entry")),
  notes = paste(
    "Each cell reports the estimate in points of the WHO 0 to 100 scale, the robust (HC2) standard error",
    "in parentheses and the p value. ANCOVA run separately for each wave: the outcome at that wave is",
    "regressed on its own baseline value and either the randomised invitation (ITT, the effect of the",
    "invitation) or attendance instrumented by the invitation (IV, the effect of attending; attending",
    "means more than five days at the centre by that wave's survey), plus the covariates of the column.",
    "(a) Main specification, no covariates. (b) Adds the baseline Katz (basic activities of daily living,",
    "0 to 6) and Lawton (instrumental activities, 0 to 8) scores. (c) Adds age at baseline, sex, indigenous,",
    "household size, education (ten levels, entered as a number) and whether the main caregiver is a son or",
    "daughter. (d) Adds both sets. (e) Lin (2013): the invitation is fully interacted with the baseline",
    "outcome, Katz and Lawton, each centred at its sample mean; the IV version keeps these interactions and",
    "instruments attendance with the invitation alone. Sample: everyone observed at baseline and at that",
    "wave; people missing a covariate drop out of that column (row Individuals). First-stage F is the",
    "robust F on the invitation, shown as the range over the two IV outcomes. IV estimates for the other",
    "domains are available on request. Months are counted from entry to the centre for those attending.",
    "WHOQOL-BREF is the mean of the standard 26 items; all scores",
    "are on the WHO 0 to 100 scale (25 times the item mean minus 1). *** p < 0.01, ** p < 0.05, * p < 0.1."),
  note_width = 150)

# ---- Table A: alternative specifications (composite and Environmental) ----------------------
get_sp <- function(o, k, e, m) spec %>% filter(outcome == o, spec == k, effect == e, model == m)
sp_cell <- function(o, k, e, m) {
  if (k == "e") {                                        # as treated: OLS under the attending column
    if (m == "ITT") return("Not\napplicable")
    return(cellp(get_sp(o, k, e, "As treated OLS"), prefix = "OLS "))
  }
  r <- get_sp(o, k, e, m); if (nrow(r) == 0) "Not\nestimated" else cellp(r)
}
sp_n <- function(k, e) {
  r <- spec %>% filter(spec == k, effect == e)
  ppl <- unique(r$n_people); rows <- unique(r$n_rows)
  stopifnot(length(ppl) == 1, length(rows) == 1)        # same sample for both outcomes
  if (k %in% c("a", "e")) sprintf("%d people", ppl) else sprintf("%d people\n%d person waves", ppl, rows)
}
rows_def <- tribble(
  ~k,  ~e,                      ~lab,
  "a", "Midline",               "(a) ANCOVA by wave\n(main specification)",
  "a", "Endline",               "",
  "b", "Average of both waves", "(b) ANCOVA, both\nwaves pooled",
  "c", "Midline",               "(c) Difference in\ndifferences, individual\nfixed effects",
  "c", "Endline",               "",
  "d", "Midline",               "(d) Difference in\ndifferences, random\neffects",
  "d", "Endline",               "",
  "sep", "",                    "Not randomised, subject\nto self selection",
  "e", "Midline",               "(e) As treated: attended\nversus did not attend\n(OLS ANCOVA by wave)",
  "e", "Endline",               "")
spec_tab <- pmap_dfr(rows_def, function(k, e, lab) {
  if (k == "sep") return(tibble(Specification = lab, Effect = "",
                                `Attending\n(IV)` = "", `Invited\n(ITT)` = "",
                                `Attending \n(IV)` = "", `Invited \n(ITT)` = "", Sample = ""))
  tibble(Specification = lab,
         Effect = if (e == "Average of both waves") "Average of\nboth waves" else e,
         `Attending\n(IV)`         = sp_cell("WHOQOL-BREF", k, e, "IV"),
         `Invited\n(ITT)`          = sp_cell("WHOQOL-BREF", k, e, "ITT"),
         `Attending \n(IV)`        = sp_cell("Environmental", k, e, "IV"),
         `Invited \n(ITT)`         = sp_cell("Environmental", k, e, "ITT"),
         Sample = sp_n(k, e))
})
# numbers for the note, from the estimation samples
n_bc   <- unique(spec$n_people[spec$spec %in% c("b", "c")]); stopifnot(length(n_bc) == 1)
n_d    <- unique(spec$n_people[spec$spec == "d"]);          stopifnot(length(n_d) == 1)
at_end <- spec %>% filter(spec == "e", effect == "Endline", outcome == "WHOQOL-BREF")
render(spec_tab,
  spans = list(list(label = "WHOQOL-BREF", l = 3, r = 4), list(label = "Environmental", l = 5, r = 6)),
  file  = paste0(OUT_TAB, "tableA_specifications.png"),
  title = "Effect of the day centre under alternative specifications",
  left_cols = 1:2,
  bold_rows = which(rows_def$k == "sep"),
  notes = paste(
    "Each cell reports the estimate in points of the WHO 0 to 100 scale, the standard error in parentheses",
    "and the p value. Attending (IV) is the effect of attending, with attendance (more than five days at the",
    "centre by that wave's survey) instrumented by the randomised invitation; Invited (ITT) is",
    "the effect of the invitation. (a) Main specification: the outcome at each wave regressed on its own",
    "baseline value and the invitation, separately for each wave; robust (HC2) errors; everyone observed at",
    "baseline and at that wave. (b) Both follow-up waves stacked, one row per person and wave, with a wave",
    "dummy, the baseline value and the invitation: one average effect over the two waves. (c) Difference in",
    "differences on all three waves: invitation times midline and invitation times endline, with wave",
    "dummies and individual fixed effects (for IV, attendance at each wave instrumented by invitation times",
    sprintf(paste("that wave). (b) and (c) use the %d people observed at baseline and at least once more;",
                  "people seen only at baseline drop out; errors clustered by person (CR2), p values from a t",
                  "distribution with %d degrees of freedom. (d) As (c) with a random intercept per person",
                  "instead of fixed effects; all %d people with a",
                  "baseline; model based errors, Satterthwaite p values; no IV version."),
            n_bc, n_bc - 1, n_d),
    sprintf(paste("(e) Not randomised: the ANCOVA of (a) with actual attendance in place of the invitation,",
                  "by OLS with HC2 errors. It compares the people who chose to attend (%d at endline) with",
                  "everyone else, invited or not (%d), so it is subject to self selection and is not an",
                  "estimate of the effect of the centre."),
            at_end$n_attending, at_end$n_not_attending),
    "WHOQOL-BREF is the mean of the standard 26 items; all scores are on the WHO 0 to 100 scale (25 times",
    "the item mean minus 1). *** p < 0.01, ** p < 0.05, * p < 0.1."),
  note_width = 140)

cat("\nwritten to", OUT_TAB, "and", OUT_CSV, "\n")
