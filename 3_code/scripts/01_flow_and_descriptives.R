# 01_flow_and_descriptives.R  (revision, September 2026; data core, October 2026)
#
# Two things, both on the data built by 00_data.R:
#   (1) THE NUMBERS BEHIND THE FLOW FIGURE (02_figures.R draws them; it reads the two csvs below)
#       - per randomisation round: people not yet invited (the pool), invited, took up a place
#       - who took up a place, and when relative to their midline survey
#       - the estimation samples of 03 by wave and group, and who was not surveyed
#   (2) WHOQOL-BREF on the standard 26 items against the 27 items of the submission
#       (the submitted total also had 16.03, a WHOQOL-100 question), same specification as 03.
#
# Definitions are the core's (00_data.R):
#   invited (z)          invited in any of the eight rounds (t1.csv); rounds come from t_pm.csv
#   attending at a wave  more than five days at the centre up to and including that wave's survey
#                        (d1, d2), the attendance rule of the main results
#   took up a place      the same rule over the whole register (18 January 2023 to 16 February
#                        2024): more than five days at the centre at any time
#   samples              everyone with the composite at baseline and at that wave (ancova_sample)
#
# Changes on 2 October 2026 (copy of the old script in _previous, suffix _before_00data):
#   - data from build_analysis_data() instead of pm_limpieza1.csv, pm_vf.csv, db.csv and pm.csv;
#     "took up a place" is now the core's rule (more than five days), not the old enrolment flag
#     of pm_vf.csv. 33 people took up a place, not 28: the 3 who entered only after their midline
#     survey (flag missing in pm_vf.csv) and 2 more who attended more than five days but had the
#     flag at 0. Per round: 8, 0, 3, 8, 3, 7, 2, 2 (the September decision log has 8, 0, 3, 7, 3,
#     5, 2, 0, which is the old flag).
#   - the balance tables (Table 2, A3, A4) and the first stage are no longer written here.
#     04_balance_tables.R and 03_main_results.R replaced them. The old csvs were moved to
#     output/csv/_previous with the suffix _superseded_by_04 or _superseded_by_03.
#   - the 26 versus 27 item comparison now uses the 0 to 100 scale, the per wave samples of 03
#     (102 and 85, not the 82 seen at all three waves), HC2 errors, and IV at both waves.
#
# Outputs (aggregates only; no folio, exact date or row of any person is written):
#   3_code/output/csv/flow_rounds.csv            one row per round plus the total
#   3_code/output/csv/flow_counts.csv            every count drawn in the flow figure
#   3_code/output/csv/whoqol_26_vs_27_items.csv  ITT and IV per wave under both item sets
#
# Run from the imss_paper folder:
#   "C:/Program Files/R/R-4.6.1/bin/Rscript.exe" 3_code/scripts/01_flow_and_descriptives.R

rm(list = ls())
suppressPackageStartupMessages({
  library(tidyverse); library(AER); library(sandwich); library(lmtest)
})
options(scipen = 999, width = 200)

source("3_code/scripts/00_paths.R")
source("3_code/scripts/00_data.R")

# ---- 0. the analysis data, with dates (INTERNAL: dates only leave as month ranges) -------------
dat <- build_analysis_data(dates = TRUE)

# month range in words, for aggregates only ("October 2022 to January 2023")
month_range <- function(d) {
  r <- format(range(d, na.rm = TRUE), "%B %Y")
  if (r[1] == r[2]) r[1] else paste(r[1], "to", r[2])
}

# =============================================================================================
# 1. RANDOMISATION ROUNDS (t_pm.csv, built in the old 2_bases_tratamiento.R)
# =============================================================================================
# Column Tk: 1 = drawn in round k, 0 = in the pool of round k and not drawn, NA = already invited
# in an earlier round (so not in the pool). The lower case tk columns track whether an invited
# person was still taking part; they are not needed here.
tp <- .read_raw_csv(RAW, "t_pm.csv") %>% select(id = x03_01, all_of(paste0("T", 1:8)))
rk <- tp %>%
  pivot_longer(-id, names_to = "round", values_to = "drawn") %>%
  mutate(round = as.integer(sub("T", "", round)), drawn = as.integer(drawn))
rounds <- rk %>% group_by(round) %>%
  summarise(pool = sum(!is.na(drawn)), invited = sum(drawn %in% 1L), .groups = "drop")
inv_round <- rk %>% filter(drawn %in% 1L) %>% group_by(id) %>%
  summarise(round = min(round), times = n(), .groups = "drop")

# design checks: one draw per person; the pool shrinks only by the people drawn (so the eight
# rounds are one uniform draw of 91 from 120 in stages); ever drawn is the core's z
stopifnot(nrow(tp) == 120, all(inv_round$times == 1),
          rounds$pool[1] == 120,
          all(rounds$pool[-1] == head(rounds$pool - rounds$invited, -1)),
          setequal(inv_round$id, dat$id[dat$z == 1]), sum(rounds$invited) == 91)

# =============================================================================================
# 2. ATTENDANCE OVER THE WHOLE REGISTER
# =============================================================================================
# The core returns attendance up to each survey (days1, days2), not over the whole register.
# The register is read here exactly as in 00_data.R (section 5.3): one column per opening day,
# header = Excel date serial, "x" or "X" = attended. As a check, the days counted here up to each
# survey date must equal the core's days1 and days2 for every person surveyed.
reg <- readxl::read_excel(paste0(RAW, "CD CASSAAM 16febrero24 asistencias.xlsx"),
                          .name_repair = "minimal")
serial <- suppressWarnings(as.numeric(names(reg)))
is_day <- !is.na(serial) & serial > 40000 & serial < 50000
stopifnot(sum(is_day) == 223, sum(names(reg) == "Folio") == 1)
marks <- reg[, c(which(names(reg) == "Folio"), which(is_day))]
names(marks) <- c("id", as.character(as.Date(serial[is_day], origin = "1899-12-30")))
visits <- marks %>%
  pivot_longer(-id, names_to = "day", values_to = "mark",
               values_transform = list(mark = as.character)) %>%
  filter(tolower(mark) %in% "x") %>%
  transmute(id, day = as.Date(day))

p <- dat %>%
  select(id, z, d1, d2, days1, days2, date0, date1, date2, entry_date,
         whoqol_y0, whoqol_y1, whoqol_y2) %>%
  left_join(inv_round %>% select(id, round), by = "id") %>%
  left_join(visits %>% count(id, name = "days_all"), by = "id") %>%
  mutate(days_all = replace_na(days_all, 0L))

# consistency with the core: days up to each survey, recounted from the same register
chk_days <- p %>%
  left_join(visits, by = "id", relationship = "one-to-many") %>%
  group_by(id) %>%
  summarise(k1 = sum(day <= date1[1], na.rm = TRUE), k2 = sum(day <= date2[1], na.rm = TRUE),
            days1 = days1[1], days2 = days2[1], .groups = "drop")
stopifnot(all(chk_days$k1[!is.na(chk_days$days1)] == chk_days$days1[!is.na(chk_days$days1)]),
          all(chk_days$k2[!is.na(chk_days$days2)] == chk_days$days2[!is.na(chk_days$days2)]))

p <- p %>%
  mutate(
    took   = days_all > 5,                                # took up a place: the core's rule, ever
    in_mid = !is.na(whoqol_y0) & !is.na(whoqol_y1),       # midline estimation sample
    in_end = !is.na(whoqol_y0) & !is.na(whoqol_y2),       # endline estimation sample
    # among those who took up a place: when, relative to their own midline survey
    took_when = case_when(
      !took                       ~ NA_character_,
      d1 %in% 1L                  ~ "attending by midline",        # more than five days by then
      !is.na(date1) & d1 %in% 0L  ~ "entered after midline",       # surveyed first, entered later
      is.na(date1) & is.na(date2) ~ "no follow up survey",
      TRUE                        ~ "other"))

# every person who took up a place falls in exactly one of the three cases above, and nobody
# never invited took up a place (one sided non compliance)
stopifnot(!any(p$took_when %in% "other"), sum(p$took & p$z == 0) == 0)
# the people attending at endline are all among those attending by midline
stopifnot(all(p$took_when[p$d2 %in% 1L] == "attending by midline"))
# those counted as entering after the midline survey had no day at all before it
late <- p %>% filter(took_when %in% "entered after midline")
stopifnot(all(late$days1 == 0), all(late$entry_date > late$date1))

rounds <- rounds %>%
  left_join(p %>% filter(took) %>% count(round, name = "took_up_place"), by = "round") %>%
  mutate(took_up_place = replace_na(took_up_place, 0L))
stopifnot(sum(rounds$took_up_place) == sum(p$took))

# ---- reconciliation with the old enrolment flag (pm_vf.csv), printed only ----------------------
# The September decision log took "entered" from pm_vf.csv (variable t). Printed so the change is
# on record; that file is not used for anything else.
vf  <- .read_raw_csv(RAW, "pm_vf.csv") %>% distinct(id = x03_01, t_old = as.integer(t))
rec <- p %>% left_join(vf, by = "id")
cat("\nPer round, people who took up a place: core rule (more than five days) against the old flag\n")
print(rounds %>%
        left_join(rec %>% filter(t_old %in% 1L) %>% count(round, name = "old_flag"), by = "round") %>%
        mutate(old_flag = replace_na(old_flag, 0L)) %>% as.data.frame(), row.names = FALSE)
cat(sprintf(paste0("Took up a place: %d under the core rule. Of them, %d also under the old flag, %d had the ",
                   "old flag missing (entered after their midline survey), %d had the old flag at 0.\n",
                   "Old flag 1 but not taking up a place under the core rule: %d.\n"),
            sum(rec$took), sum(rec$took & rec$t_old %in% 1L), sum(rec$took & is.na(rec$t_old)),
            sum(rec$took & rec$t_old %in% 0L), sum(!rec$took & rec$t_old %in% 1L)))

# =============================================================================================
# 3. FLOW COUNTS
# =============================================================================================
grp_of <- function(z, d) case_when(z == 0 ~ "Never invited", d == 1 ~ "Attending",
                                   TRUE ~ "Invited, not attending")
smp <- bind_rows(
  p %>% filter(in_mid) %>% transmute(wave = "midline", grp = grp_of(z, d1)),
  p %>% filter(in_end) %>% transmute(wave = "endline", grp = grp_of(z, d2))) %>%
  count(wave, grp)
n_smp <- function(w, g) smp$n[smp$wave == w & smp$grp == g]

# the samples are those of 03 (ancova_sample on the composite)
stopifnot(sum(smp$n[smp$wave == "midline"]) == nrow(ancova_sample(dat, "whoqol", 1)),
          sum(smp$n[smp$wave == "endline"]) == nrow(ancova_sample(dat, "whoqol", 2)))

took_mid <- p %>% filter(took_when %in% "attending by midline")
flow <- tribble(
  ~item,                                ~n,                                             ~months,
  "expressed_interest_baseline",        nrow(p),                                        month_range(p$date0),
  "invited",                            sum(p$z == 1),                                  NA,
  "never_invited",                      sum(p$z == 0),                                  NA,
  "took_up_place",                      sum(p$took),                                    month_range(p$entry_date[p$took]),
  "did_not_take_up_place",              sum(p$z == 1 & !p$took),                        NA,
  "took_up_attending_by_midline",       nrow(took_mid),                                 NA,
  "took_up_attending_by_midline_surveyed_endline", sum(!is.na(took_mid$date2)),         NA,
  "took_up_entered_after_midline",      nrow(late),                                     NA,
  "took_up_no_follow_up_survey",        sum(p$took_when %in% "no follow up survey"),    NA,
  "midline_attending",                  n_smp("midline", "Attending"),                  NA,
  "midline_invited_not_attending",      n_smp("midline", "Invited, not attending"),     NA,
  "midline_never_invited",              n_smp("midline", "Never invited"),              NA,
  "midline_total",                      sum(p$in_mid),                                  month_range(p$date1),
  "endline_attending",                  n_smp("endline", "Attending"),                  NA,
  "endline_invited_not_attending",      n_smp("endline", "Invited, not attending"),     NA,
  "endline_never_invited",              n_smp("endline", "Never invited"),              NA,
  "endline_total",                      sum(p$in_end),                                  month_range(p$date2),
  "midline_not_surveyed_invited",       sum(p$z == 1 & !p$in_mid),                      NA,
  "midline_not_surveyed_never_invited", sum(p$z == 0 & !p$in_mid),                      NA,
  "endline_not_surveyed_invited",       sum(p$z == 1 & !p$in_end),                      NA,
  "endline_not_surveyed_never_invited", sum(p$z == 0 & !p$in_end),                      NA,
  "all_three_waves",                    sum(p$in_mid & p$in_end),                       NA)
fv <- setNames(flow$n, flow$item)

# arithmetic of the figure: every box adds up
stopifnot(
  fv["invited"] + fv["never_invited"] == fv["expressed_interest_baseline"],
  fv["took_up_place"] + fv["did_not_take_up_place"] == fv["invited"],
  fv["took_up_attending_by_midline"] + fv["took_up_entered_after_midline"] +
    fv["took_up_no_follow_up_survey"] == fv["took_up_place"],
  fv["midline_attending"] == fv["took_up_attending_by_midline"],
  fv["endline_attending"] == fv["took_up_attending_by_midline_surveyed_endline"],
  fv["midline_attending"] + fv["midline_invited_not_attending"] + fv["midline_never_invited"] == fv["midline_total"],
  fv["endline_attending"] + fv["endline_invited_not_attending"] + fv["endline_never_invited"] == fv["endline_total"],
  fv["midline_total"] == 102, fv["endline_total"] == 85,
  fv["midline_attending"] + fv["midline_invited_not_attending"] + fv["midline_not_surveyed_invited"] == fv["invited"],
  fv["endline_attending"] + fv["endline_invited_not_attending"] + fv["endline_not_surveyed_invited"] == fv["invited"],
  fv["midline_never_invited"] + fv["midline_not_surveyed_never_invited"] == fv["never_invited"],
  fv["endline_never_invited"] + fv["endline_not_surveyed_never_invited"] == fv["never_invited"])

rounds_out <- bind_rows(
  rounds %>% mutate(round = as.character(round)),
  tibble(round = "Total", pool = NA_integer_, invited = sum(rounds$invited),
         took_up_place = sum(rounds$took_up_place)))
write_csv(rounds_out, paste0(OUT_CSV, "flow_rounds.csv"), na = "")
write_csv(flow, paste0(OUT_CSV, "flow_counts.csv"), na = "")

cat("\n=== flow: rounds ===\n"); print(as.data.frame(rounds_out), row.names = FALSE)
cat("\n=== flow: counts ===\n"); print(as.data.frame(flow), row.names = FALSE)

# =============================================================================================
# 4. WHOQOL-BREF ON 26 ITEMS (STANDARD) AGAINST THE 27 ITEMS OF THE SUBMISSION
# =============================================================================================
# The submitted composite was the mean of the 26 BREF items plus 16.03 ("do you get the kind of
# support you need from others?"), not reversed, as rebuilt and checked in 00_data.R (5.6 a).
# On the WHO 0 to 100 scale the 27 item score is 25 x ((26 x item mean of the 26 + answer to
# 16.03) / 27 minus 1). Specification exactly as in 03: ANCOVA per wave, HC2, IV with the
# baseline value in both stages.
q03 <- .read_raw_csv(RAW, "pm_limpieza1.csv") %>%
  filter(name == "x16_03") %>%
  select(id = x03_01, t0, t1, t2) %>%
  pivot_longer(-id, names_to = "wave", values_to = "label") %>%
  transmute(id, wave = as.integer(sub("t", "", wave)),
            q03 = unname(WHOQOL_RESPONSES[.ascii(label)])) %>%
  pivot_wider(names_from = wave, values_from = q03, names_prefix = "q03_")
w <- dat %>% select(id, z, d1, d2, starts_with("whoqol_y")) %>% left_join(q03, by = "id")
for (k in 0:2) {
  m26 <- w[[paste0("whoqol_y", k)]] / 25 + 1                          # back to the 1 to 5 item mean
  w[[paste0("whoqol27_y", k)]] <- 25 * ((26 * m26 + w[[paste0("q03_", k)]]) / 27 - 1)
}
# 16.03 was answered in every survey held, so both versions cover the same person waves
stopifnot(identical(is.na(w$whoqol_y0), is.na(w$whoqol27_y0)),
          identical(is.na(w$whoqol_y1), is.na(w$whoqol27_y1)),
          identical(is.na(w$whoqol_y2), is.na(w$whoqol27_y2)))

# correlation of the two versions over every survey held (person by wave)
long26 <- c(w$whoqol_y0, w$whoqol_y1, w$whoqol_y2)
long27 <- c(w$whoqol27_y0, w$whoqol27_y1, w$whoqol27_y2)
r2627  <- cor(long26, long27, use = "complete.obs")
n_surv <- sum(!is.na(long26))

hc2 <- function(m, term) {
  ct <- unclass(coeftest(m, vcov. = vcovHC(m, type = "HC2")))
  list(est = unname(ct[term, 1]), se = unname(ct[term, 2]), p = unname(ct[term, 4]))
}
versions <- c(whoqol = "26 items (standard, used)", whoqol27 = "27 items (submitted)")
cmp <- map_dfr(names(versions), function(v) {
  map_dfr(1:2, function(wv) {
    s <- w %>%
      transmute(z, d = .data[[paste0("d", wv)]], y0 = .data[[paste0(v, "_y0")]],
                y = .data[[paste0(v, "_y", wv)]]) %>%
      filter(!is.na(y0), !is.na(y))
    e_itt <- hc2(lm(y ~ y0 + z, s), "z")              # e_ names: tibble() below creates itt, iv
    e_iv  <- hc2(ivreg(y ~ d + y0 | z + y0, data = s), "d")
    tibble(version = unname(versions[v]), wave = c("midline", "endline")[wv], n = nrow(s),
           baseline_mean = mean(s$y0), baseline_sd = sd(s$y0),
           itt = e_itt$est, itt_se = e_itt$se, itt_p = e_itt$p,
           iv = e_iv$est, iv_se = e_iv$se, iv_p = e_iv$p,
           itt_sd = e_itt$est / sd(s$y0))
  })
}) %>%
  mutate(corr_26_27_all_surveys = r2627, surveys = n_surv)

# the 26 item rows must be the composite of the main results table (03)
mr_file <- paste0(OUT_CSV, "main_results.csv")
if (file.exists(mr_file)) {
  mr <- read_csv(mr_file, show_col_types = FALSE) %>% filter(outcome == "WHOQOL-BREF")
  a  <- cmp %>% filter(version == versions[["whoqol"]]) %>% mutate(wv = if_else(wave == "midline", 1, 2))
  for (i in seq_len(nrow(a))) {
    stopifnot(abs(a$itt[i] - mr$est[mr$model == "ITT" & mr$wave == a$wv[i]]) < 1e-9,
              abs(a$iv[i]  - mr$est[mr$model == "IV"  & mr$wave == a$wv[i]]) < 1e-9)
  }
  cat("\n26 item rows reproduce the WHOQOL-BREF rows of main_results.csv (03).\n")
}

write_csv(cmp %>% mutate(across(where(is.double), ~ round(.x, 4))),
          paste0(OUT_CSV, "whoqol_26_vs_27_items.csv"))
cat("\n=== WHOQOL-BREF, 26 versus 27 items (0 to 100 scale, HC2) ===\n")
print(cmp %>% mutate(across(where(is.double), ~ round(.x, 3))) %>% select(-surveys) %>% as.data.frame(),
      row.names = FALSE)
cat(sprintf("correlation of the two versions over the %d surveys held: %.4f\n", n_surv, r2627))
dd <- cmp %>% select(version, wave, itt, iv) %>%
  pivot_wider(names_from = version, values_from = c(itt, iv))
cat(sprintf("largest change from 27 to 26 items: ITT %.2f points, IV %.2f points\n",
            max(abs(dd[[2]] - dd[[3]])), max(abs(dd[[4]] - dd[[5]]))))

cat("\nwritten to", OUT_CSV, ": flow_rounds.csv, flow_counts.csv, whoqol_26_vs_27_items.csv\n")
