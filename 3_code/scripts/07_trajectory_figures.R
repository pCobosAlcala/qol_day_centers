# 07_trajectory_figures.R  (revision, October 2026; data core and 0 to 100 scale, 2 October 2026)
#
# The two descriptive figures of the paper, rebuilt on the analysis data of 00_data.R:
#
#   Figure A  figure_individual_changes.png
#             One vertical arrow per person, from the baseline to the endline WHOQOL-BREF composite,
#             by group. Sample: the ENDLINE estimation sample (85 people observed at baseline and
#             endline), the same people behind the endline column of the main table (03).
#             Submitted version: 1_paper/figures_submitted/fig1.jpg.
#
#   Figure B  figure_trajectories.png
#             Composite and the five parts of the WHOQOL-BREF at baseline, midline and endline,
#             by group. Thin lines are individuals; the thick line with dots is the group mean.
#             Sample: people observed at all three waves (82). Submitted version:
#             1_paper/figures_submitted/fig2.jpg (boxplots plus means); here the boxplots are
#             replaced by individual lines, as agreed for the revision.
#
# Timing, stated in both captions and computed here from the survey and entry dates:
#   baseline to endline is about a year for EVERYONE (11.6 months on average); "about nine months"
#   is entry to endline and applies only to those who attended (baseline came about three months
#   before entry). The main table's "about 3 / about 9 months" are entry based.
#
# Groups (the same three groups everywhere in the revision):
#   Never invited             z = 0          (never drawn in any of the eight rounds)
#   Invited, did not attend   z = 1, d2 = 0
#   Attended                  z = 1, d2 = 1  (more than five days at the centre by the endline survey)
# Nobody with z = 0 attended (one-sided non-compliance), so the three groups are exhaustive.
# d2 is the endline attendance flag of 00_data.R; it equals the old time invariant d for all 85.
# In the 82 of Figure B every attender had already attended more than five days by midline (d1 = 1),
# so the "Attended" group is the same whichever wave defines it (checked below).
#
# Changes against the version of 2 October that read pm_who.csv (copies of the old script, figures
# and csv files are in the _previous subfolders, suffix _before_00data):
#   - data from build_analysis_data() in 00_data.R instead of pm_who.csv and pm_limpieza1.csv;
#     the 85 and the 82 are the same people and the same groups as before (checked below)
#   - scores on the WHO 0 to 100 scale, 25 x (item mean minus 1), instead of item means on 1 to 5.
#     Points are 25 times larger; the per cents of the group's own baseline CHANGE, because the
#     scale's zero moves (the denominator shrinks), as in the magnitudes table of 03
#   - Figure B gets a note in every panel (change in the group mean, baseline to endline)
#   - captions give baseline to endline (about a year, everyone) and entry to endline (about nine
#     months, attenders) from the data
#
# Outputs: two PNGs (300 dpi) in OUT_FIG and three csv files of GROUP aggregates in OUT_CSV:
#   trajectory_group_means.csv      group means and SD by outcome and wave, 82 people, 0 to 100
#   individual_changes_endline.csv  rose / fell / unchanged and mean change by group, 85 people
#   trajectory_timing.csv           mean months between surveys and from entry, by sample and group
# No individual rows, identifiers or dates leave this script (dates are read with dates = TRUE and
# only group means of the gaps between them are kept).
# Run from the imss_paper folder:
#   "C:/Program Files/R/R-4.6.1/bin/Rscript.exe" 3_code/scripts/07_trajectory_figures.R

rm(list = ls())
suppressPackageStartupMessages({ library(tidyverse) })
options(scipen = 999, width = 200)

source("3_code/scripts/00_paths.R")
source("3_code/scripts/00_data.R")

# ---- 1. data: one row per randomised person, scores on 0 to 100 ---------------------------
# dates = TRUE adds the survey dates (date0, date1, date2) and the first day attended
# (entry_date). They are used only to compute group means of the gaps between surveys.
dat <- build_analysis_data(dates = TRUE)

MONTH <- 30.44                                            # average month length in days

# the three groups, in the order of the submitted figures (control first, attenders last);
# defined by endline attendance (d2), see header
glev <- c("Never invited", "Invited, did not attend", "Attended")
grp_of <- function(z, d) factor(case_when(z == 0 ~ glev[1], d == 0 ~ glev[2], d == 1 ~ glev[3]),
                                levels = glev)

# minus sign for printed numbers: the true minus (U+2212), so nothing on a figure reads as a dash
minus <- function(x) gsub("-", "\u2212", x, fixed = TRUE)
# signed number with one decimal for the panel notes: "+4.2", minus 5.4, and a plain "0.0" when
# the change rounds to zero (a sign on zero would suggest a direction that is not there)
sgn <- function(x, digits = 1) {
  out <- sprintf(paste0("%+.", digits, "f"), x)
  out[round(x, digits) == 0] <- sprintf(paste0("%.", digits, "f"), 0)
  minus(out)
}

# ---- 2. Figure A: individual change, baseline to endline (endline sample, 85) -------------
# ancova_sample() gives exactly the endline ANCOVA sample of 03: everyone with the composite at
# baseline and endline, with d = attendance by the endline survey (d2).
chg <- ancova_sample(dat, "whoqol", 2) %>%
  left_join(dat %>% select(id, date0, date2, entry_date), by = "id") %>%
  mutate(group = grp_of(z, d),
         # Round the change before classifying it, as a safeguard. On 00_data.R a person whose
         # answers did not change comes out as an exact zero, but if the composite arithmetic ever
         # changes, a floating point residue of about 1e-14 could make them count as "Score fell"
         # and be drawn as a zero length arrow. The smallest genuine change is one
         # answer step on one item, 25 / 26 = 0.96 points, so rounding to 8 decimals cannot touch it.
         change    = round(y - y0, 8),
         direction = factor(case_when(change > 0 ~ "Score rose", change < 0 ~ "Score fell",
                                      TRUE ~ "No change"),
                            levels = c("Score rose", "Score fell", "No change")),
         months_b2e     = as.numeric(date2 - date0) / MONTH,               # baseline to endline
         months_entry2e = if_else(d == 1, as.numeric(date2 - entry_date) / MONTH, NA_real_))
n_end <- nrow(chg)
stopifnot(n_end == 85, !any(is.na(chg$group)),
          min(abs(chg$change[chg$change != 0])) > 0.9,     # no residue left among the "changes"
          !any(is.na(chg$months_b2e)), all(!is.na(chg$months_entry2e[chg$d == 1])))

# timing for the caption: "about a year" for everyone, "about nine months after entry" for attenders
gap_b2e     <- mean(chg$months_b2e)
gap_entry2e <- mean(chg$months_entry2e, na.rm = TRUE)
stopifnot(abs(gap_b2e - 12) < 1, abs(gap_entry2e - 9) < 0.75)
cat(sprintf("baseline to endline (85): mean %.1f months, range %.1f to %.1f\n",
            gap_b2e, min(chg$months_b2e), max(chg$months_b2e)))
cat(sprintf("entry to endline (%d attenders): mean %.1f months, SD %.1f, range %.1f to %.1f\n",
            sum(chg$d), gap_entry2e, sd(chg$months_entry2e, na.rm = TRUE),
            min(chg$months_entry2e, na.rm = TRUE), max(chg$months_entry2e, na.rm = TRUE)))

# Horizontal position: people spread evenly across each panel, ordered by baseline score (even
# spacing keeps arrows from overlapping and makes the figure reproducible; the submitted figure
# used random jitter).
chg <- chg %>%
  group_by(group) %>%
  arrange(y0, .by_group = TRUE) %>%
  mutate(xpos = (row_number() - 0.5) / n()) %>%
  ungroup()

# Panel summaries for the strip and a two line note inside each panel. Rose, fell and unchanged
# use the same rounded change as the arrows, so the counts always match what is drawn.
# Per cents first (Pablo's preference): change in the group mean as a per cent of the group's own
# baseline mean, then the same change in points of the 0 to 100 scale in parentheses.
sumA <- chg %>%
  group_by(group) %>%
  summarise(n = n(), rose = sum(change > 0), fell = sum(change < 0), unchanged = sum(change == 0),
            mean_baseline = mean(y0), mean_endline = mean(y), mean_change = mean(y - y0),
            # size of the changes among those whose score fell or rose (the text says the declines
            # of attenders who fell were smaller than in the other groups)
            mean_change_fell = mean(change[change < 0]), mean_change_rose = mean(change[change > 0]),
            pct_change = 100 * mean(y - y0) / mean(y0),
            mean_months_baseline_to_endline = mean(months_b2e),
            mean_months_entry_to_endline = if (all(d == 1)) mean(months_entry2e) else NA_real_,
            .groups = "drop") %>%
  mutate(strip = sprintf("%s\n(n = %d)", group, n),
         note  = sprintf("%d rose, %d fell, %d unchanged\nmean change %s%% (%s points)",
                         rose, fell, unchanged, sgn(pct_change), sgn(mean_change)))
stopifnot(all(sumA$rose + sumA$fell + sumA$unchanged == sumA$n))
strip_A <- setNames(sumA$strip, as.character(sumA$group))

capA <- paste0(
  "Each arrow is one person, from the baseline score to the endline score. Baseline to endline is about ",
  "a year for everyone (", sprintf("%.1f", gap_b2e), " months on average); for those who attended, ",
  "endline came about nine months after entry (", sprintf("%.1f", gap_entry2e), " months on average). ",
  "Sample: the ", n_end, " people observed at baseline and endline (the endline estimation sample). ",
  "Attended means more than five days at the centre by the endline survey. People are ordered by ",
  "baseline score within each panel; a grey dot is a person whose score did not change. WHOQOL-BREF ",
  "composite of the standard 26 items on the WHO 0 to 100 scale (25 times the item mean minus 1). ",
  "Panel notes: change in the group mean as a per cent of the group's own baseline mean, and in ",
  "points of the 0 to 100 scale in parentheses.")

figA <- ggplot(chg) +
  geom_segment(data = chg %>% filter(direction != "No change"),
               aes(x = xpos, xend = xpos, y = y0, yend = y, colour = direction, linetype = direction),
               arrow = arrow(length = unit(0.12, "cm")), alpha = 0.75, linewidth = 0.7) +
  geom_point(data = chg %>% filter(direction == "No change"),
             aes(x = xpos, y = y0), colour = "grey40", size = 1.2) +
  # note at the foot of each panel, only if IMAGE_TEXT (00_paths.R); the lowest score is about 20
  (if (IMAGE_TEXT) geom_text(data = sumA, aes(x = 0.5, y = 6.5, label = note),
                             size = 3.0, colour = "grey25", lineheight = 1.0)) +
  scale_colour_manual(values = c("Score rose" = "#2C7FB8", "Score fell" = "#D7301F"),
                      breaks = c("Score rose", "Score fell"), name = NULL) +
  scale_linetype_manual(values = c("Score rose" = "solid", "Score fell" = "dashed"),
                        breaks = c("Score rose", "Score fell"), name = NULL) +
  facet_wrap(~ group, nrow = 1, labeller = as_labeller(strip_A)) +
  scale_y_continuous(breaks = seq(0, 100, 25)) +
  scale_x_continuous(expand = expansion(add = 0.04)) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(x = NULL, y = "WHOQOL-BREF (0 to 100), baseline to endline",
       caption = if (IMAGE_TEXT) str_wrap(capA, width = 128)) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "top", legend.direction = "horizontal", legend.text = element_text(size = 11),
        legend.key.width = unit(1.2, "cm"),
        strip.text = element_text(size = 11, face = "bold"),
        strip.background = element_rect(fill = "grey95", colour = NA),
        panel.border = element_rect(colour = "grey80", fill = NA, linewidth = 0.7),
        panel.spacing = unit(1.0, "lines"),
        plot.caption = element_text(hjust = 0, size = 8.5, colour = "grey30", lineheight = 1.1),
        plot.caption.position = "plot") +
  guides(colour = guide_legend(override.aes = list(linewidth = 1, alpha = 1)))

stopifnot(min(c(chg$y0, chg$y)) > 15)                     # the panel note area stays free of data
ggsave(paste0(OUT_FIG, "figure_individual_changes.png"), figA, width = 7.5, height = 5.8,
       dpi = 300, bg = "white", device = ragg::agg_png)
if (!IMAGE_TEXT) save_note(paste0(OUT_FIG, "figure_individual_changes.png"), "", c(capA, sumA$note))
cat("figure_individual_changes.png written\n")

# ---- 3. Figure B: trajectories, people observed at all three waves (82) -------------------
# A person is in the panel if the composite exists at all three waves; the five parts are then
# present too (the composite requires them), so every row of the figure uses the same 82 people.
ycols <- paste0(rep(names(OUTCOMES), each = 3), "_y", 0:2)
p82 <- dat %>%
  filter(!is.na(whoqol_y0), !is.na(whoqol_y1), !is.na(whoqol_y2)) %>%
  mutate(group = grp_of(z, d2))
n_all <- nrow(p82)
stopifnot(n_all == 82, !anyNA(p82[, ycols]), !any(is.na(p82$group)),
          identical(p82$d1, p82$d2))                      # every attender had attended by midline

trj <- p82 %>%
  select(id, group, all_of(ycols)) %>%
  pivot_longer(all_of(ycols), names_to = c("outcome", "wave"), names_pattern = "(.*)_y(\\d)",
               values_to = "value") %>%
  mutate(outcome = factor(OUTCOMES[outcome], levels = OUTCOMES),
         wave    = factor(as.integer(wave), levels = 0:2, labels = c("Baseline", "Midline", "Endline")))

# group means by outcome and wave, with the change from the group's own baseline mean
means <- trj %>%
  group_by(outcome, group, wave) %>%
  summarise(n = n(), mean = mean(value), sd = sd(value), .groups = "drop") %>%
  group_by(outcome, group) %>%
  mutate(change_from_baseline_pts = mean - mean[wave == "Baseline"],
         change_from_baseline_pct = 100 * change_from_baseline_pts / mean[wave == "Baseline"]) %>%
  ungroup()

# panel note: baseline to endline change in the group mean, per cent of the group's own baseline
# first, points in parentheses
notesB <- means %>%
  filter(wave == "Endline") %>%
  mutate(note = sprintf("%s%% (%s points)", sgn(change_from_baseline_pct),
                        sgn(change_from_baseline_pts)))

# timing for the caption, in the 82: months from baseline to each later survey for everyone, and
# from entry to each later survey for attenders (all of whom had entered by midline)
tim82 <- p82 %>%
  mutate(b2m = as.numeric(date1 - date0) / MONTH, b2e = as.numeric(date2 - date0) / MONTH,
         e2m = as.numeric(date1 - entry_date) / MONTH, e2e = as.numeric(date2 - entry_date) / MONTH)
stopifnot(all(tim82$e2m[tim82$d2 == 1] > 0))              # entry came before midline for attenders
b2m82 <- mean(tim82$b2m); b2e82 <- mean(tim82$b2e)
e2m82 <- mean(tim82$e2m[tim82$d2 == 1]); e2e82 <- mean(tim82$e2e[tim82$d2 == 1])
stopifnot(abs(b2e82 - 12) < 1, abs(e2e82 - 9) < 0.75)

# column strips carry each group's size
nB <- p82 %>% count(group)
strip_B <- setNames(sprintf("%s\n(n = %d)", nB$group, nB$n), as.character(nB$group))

capB <- paste0(
  "Thin lines are individuals; the thick line with dots is the group mean. Sample: the ", n_all,
  " people observed at all three waves. Baseline to endline is about a year for everyone (",
  sprintf("%.1f", b2e82), " months on average) and baseline to midline about ", sprintf("%.0f", b2m82),
  " months, in overlapping calendar windows. For those who attended, midline came about ",
  sprintf("%.1f", e2m82), " months and endline about nine months (", sprintf("%.1f", e2e82),
  " months on average) after entry. Attended means more than five days at the centre by the endline survey; ",
  "all of them had passed that mark by midline. Panel notes: change in the group mean from baseline to ",
  "endline, as a per cent of the group's own baseline mean, and in points in parentheses. Scores on the ",
  "WHO 0 to 100 scale (25 times the item mean minus 1); WHOQOL-BREF is the composite of the standard ",
  "26 items.")

figB <- ggplot(trj, aes(x = wave, y = value)) +
  geom_line(aes(group = id), linewidth = 0.3, alpha = 0.18, colour = "grey25") +
  geom_line(data = means, aes(y = mean, group = 1), linewidth = 1.4, colour = "black") +
  geom_point(data = means, aes(y = mean), size = 2.6, colour = "black") +
  # the note sits in a band below 0, which no score can reach (only if IMAGE_TEXT, 00_paths.R)
  (if (IMAGE_TEXT) geom_text(data = notesB, aes(x = 2, y = -11, label = note),
                             size = 2.7, colour = "grey20")) +
  facet_grid(rows = vars(outcome), cols = vars(group), switch = "y",
             labeller = labeller(group = strip_B)) +
  scale_y_continuous(breaks = c(0, 50, 100)) +
  scale_x_discrete(expand = expansion(add = 0.35)) +
  coord_cartesian(ylim = c(if (IMAGE_TEXT) -19 else -3, 100)) +
  labs(x = NULL, y = NULL, caption = if (IMAGE_TEXT) str_wrap(capB, width = 142)) +
  theme_minimal(base_size = 11) +
  theme(panel.spacing = unit(0.5, "lines"),
        panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0, hjust = 1, size = 10.5),
        strip.text.x = element_text(size = 10.5, face = "bold", margin = margin(b = 4)),
        panel.border = element_rect(colour = "grey85", fill = NA, linewidth = 0.5),
        plot.caption = element_text(hjust = 0, size = 8.5, colour = "grey30", lineheight = 1.1),
        plot.caption.position = "plot",
        plot.margin = margin(5, 8, 5, 5))

ggsave(paste0(OUT_FIG, "figure_trajectories.png"), figB, width = 8, height = 11,
       dpi = 300, bg = "white", device = ragg::agg_png)
if (!IMAGE_TEXT) save_note(paste0(OUT_FIG, "figure_trajectories.png"), "", c(capB, notesB$note))
cat("figure_trajectories.png written\n")

# ---- 4. numbers for the text (group aggregates only) --------------------------------------
means_out <- means %>%
  mutate(sample = sprintf("observed at all three waves (n = %d)", n_all),
         scale  = "WHO 0 to 100",
         across(c(mean, sd, change_from_baseline_pts, change_from_baseline_pct), ~ round(.x, 3))) %>%
  arrange(outcome, group, wave)
write_csv(means_out, paste0(OUT_CSV, "trajectory_group_means.csv"))

chg_out <- sumA %>% select(-strip, -note) %>%
  mutate(sample = sprintf("observed at baseline and endline (n = %d)", n_end),
         scale  = "WHO 0 to 100",
         across(c(mean_baseline, mean_endline, mean_change, pct_change, mean_change_fell, mean_change_rose),
                ~ round(.x, 3)),
         across(starts_with("mean_months"), ~ round(.x, 1)))
write_csv(chg_out, paste0(OUT_CSV, "individual_changes_endline.csv"))

# mean months between surveys, and from entry for attenders, by sample and group
timing <- bind_rows(
  chg %>% transmute(sample = sprintf("baseline and endline (n = %d)", n_end), group, d,
                    b2m = NA_real_, b2e = months_b2e, e2m = NA_real_, e2e = months_entry2e),
  tim82 %>% transmute(sample = sprintf("all three waves (n = %d)", n_all), group, d = d2,
                      b2m, b2e, e2m = if_else(d2 == 1, e2m, NA_real_), e2e = if_else(d2 == 1, e2e, NA_real_))) %>%
  group_by(sample, group) %>%
  summarise(n = n(),
            months_baseline_to_midline = mean(b2m), months_baseline_to_endline = mean(b2e),
            months_entry_to_midline = mean(e2m), months_entry_to_endline = mean(e2e),
            .groups = "drop") %>%
  mutate(across(starts_with("months"), ~ round(.x, 1)))
write_csv(timing, paste0(OUT_CSV, "trajectory_timing.csv"))

# ---- 5. self-checks against the version built on pm_who.csv -------------------------------
# (a) rose / fell / unchanged: the rescaling multiplies every change by 25, so the counts of the
#     version verified by hand on 2 October 2026 must hold exactly
ref_dir <- tribble(~group, ~rose, ~fell, ~unchanged,
  "Never invited", 8L, 14L, 1L, "Invited, did not attend", 13L, 23L, 2L, "Attended", 12L, 11L, 1L)
cmp_dir <- chg_out %>% mutate(group = as.character(group)) %>%
  inner_join(ref_dir, by = "group", suffix = c("", "_ref"))
stopifnot(nrow(cmp_dir) == 3, all(cmp_dir$rose == cmp_dir$rose_ref),
          all(cmp_dir$fell == cmp_dir$fell_ref), all(cmp_dir$unchanged == cmp_dir$unchanged_ref),
          all(cmp_dir$n == c(23, 38, 24)))
# (b) composite group means of the 82, back on the 1 to 5 scale (mean / 25 + 1), must equal the
#     values reported for the 26 item composite before the switch (same people, same groups)
ref <- tribble(~group, ~wave, ~ref,
  "Attended", "Baseline", 3.208, "Attended", "Midline", 3.106, "Attended", "Endline", 3.301,
  "Never invited", "Baseline", 3.101, "Never invited", "Midline", 2.978, "Never invited", "Endline", 2.918,
  "Invited, did not attend", "Baseline", 3.309, "Invited, did not attend", "Midline", 3.229,
  "Invited, did not attend", "Endline", 3.172)
cmp <- means %>% filter(outcome == "WHOQOL-BREF") %>%
  mutate(group = as.character(group), wave = as.character(wave)) %>%
  inner_join(ref, by = c("group", "wave"))
stopifnot(nrow(cmp) == 9, all(abs(round(cmp$mean / 25 + 1, 3) - cmp$ref) < 1e-9))
cat("\ncheck passed: counts and composite group means match the version built on pm_who.csv\n")

# ---- 6. console summary -------------------------------------------------------------------
cat("\n=== Figure A: individual change baseline to endline (85), 0 to 100 ===\n")
print(as.data.frame(chg_out %>% select(-sample, -scale)), row.names = FALSE)
cat("\n=== Figure B: group means by wave (82), 0 to 100 ===\n")
print(means_out %>% select(outcome, group, wave, n, mean, sd, change_from_baseline_pts,
                           change_from_baseline_pct) %>% as.data.frame(), row.names = FALSE)
cat("\n=== timing (months) ===\n")
print(as.data.frame(timing), row.names = FALSE)
cat("\nwritten to", OUT_FIG, "and", OUT_CSV, "\n")
