# 02_figures.R  (revision, September 2026; data core, October 2026)
#
# Figure A: participant flow, from the 120 who expressed interest, through the eight randomisation
#           rounds, to the estimation samples of the main results (102 at midline, 85 at endline)
# Figure B: timing of the surveys, as a SCHEMATIC (5 October 2026): illustrative positions only,
#           no real person and no computed individual timing; the numbers it stands for are
#           recomputed from the data and printed to the console (group means only)
# Saved as PNG for the coauthor meeting. Publication versions can be re exported as TIFF.
#
# Change on 5 October 2026 (copies in the _previous subfolders, suffix _before_schematic): Figure B
# redrawn as a schematic at Pablo's request; the flow figure is unchanged.
#
# Definitions are the core's (00_data.R), and the flow numbers are computed in
# 01_flow_and_descriptives.R (flow_rounds.csv, flow_counts.csv), so RUN 01 FIRST. This script
# only draws them, after checking them against the core.
#   attending at a wave  more than five days at the centre up to that wave's survey (d1, d2)
#   took up a place      more than five days at the centre at any time in the register
#   samples              baseline and that wave observed (ancova_sample): 102 midline, 85 endline
#
# Changes on 2 October 2026 (copies of the old script and figures in the _previous subfolders,
# suffix _before_00data):
#   - flow: per round counts added (pool, invited, took up a place); "took up a place" follows the
#     core's rule, 33 people, not the 28 of the old enrolment flag; the three who entered only after
#     their midline survey are now in the analysis (counted as not attending at midline) instead of
#     "excluded"; the estimation samples are 102 and 85 by group; the box with the people observed
#     at all three waves is gone; dates are given to the month only.
#   - timeline: same anonymised design, but built from the core (survey dates, entry date = first
#     day in the register, which equals the old date_primer_dia for every attender) and on the
#     estimation samples of each wave instead of the 82 people observed at all three waves.
#
# Run from the imss_paper folder:
#   "C:/Program Files/R/R-4.6.1/bin/Rscript.exe" 3_code/scripts/02_figures.R

rm(list = ls())
suppressPackageStartupMessages({library(tidyverse); library(patchwork)})
options(scipen = 999, width = 200)

source("3_code/scripts/00_paths.R")
source("3_code/scripts/00_data.R")

# INTERNAL: dates of identifiable people are used for means and month ranges, never written
dat <- build_analysis_data(dates = TRUE)

# =============================================================================================
# Figure A. Participant flow
# =============================================================================================
f_rounds <- paste0(OUT_CSV, "flow_rounds.csv"); f_counts <- paste0(OUT_CSV, "flow_counts.csv")
if (!file.exists(f_rounds) || !file.exists(f_counts))
  stop("flow_rounds.csv or flow_counts.csv missing: run 01_flow_and_descriptives.R first")
rounds <- read_csv(f_rounds, show_col_types = FALSE, col_types = cols(round = "c"))
counts <- read_csv(f_counts, show_col_types = FALSE)
n  <- setNames(counts$n, counts$item)                       # n["took_up_place"] etc.
mo <- setNames(counts$months, counts$item)                  # month ranges, words

# ---- the csvs must agree with the core as it stands now ----------------------------------------
grp_of <- function(z, d) case_when(z == 0 ~ "Never invited", d == 1 ~ "Attending",
                                   TRUE ~ "Invited, not attending")
core_smp <- bind_rows(
  ancova_sample(dat, "whoqol", 1) %>% transmute(wave = "midline", grp = grp_of(z, d)),
  ancova_sample(dat, "whoqol", 2) %>% transmute(wave = "endline", grp = grp_of(z, d))) %>%
  count(wave, grp)
cs <- function(w, g) core_smp$n[core_smp$wave == w & core_smp$grp == g]
stopifnot(n["expressed_interest_baseline"] == nrow(dat), n["invited"] == sum(dat$z),
          n["never_invited"] == sum(dat$z == 0),
          n["midline_attending"] == cs("midline", "Attending"),
          n["midline_invited_not_attending"] == cs("midline", "Invited, not attending"),
          n["midline_never_invited"] == cs("midline", "Never invited"),
          n["endline_attending"] == cs("endline", "Attending"),
          n["endline_invited_not_attending"] == cs("endline", "Invited, not attending"),
          n["endline_never_invited"] == cs("endline", "Never invited"),
          sum(as.numeric(rounds$invited[rounds$round != "Total"])) == n["invited"],
          sum(as.numeric(rounds$took_up_place[rounds$round != "Total"])) == n["took_up_place"])

# Round dates are not in the data. Our files date round 1 on 2 January 2023 and rounds 7 and 8 by the
# save dates of their selection scripts (late July and early September 2023), but the attendance
# register shows invitees of rounds 7 and 8 attending before those dates, and it labels two round 8
# invitees as round 7 (08_text_numbers.R, section i). We do not ask IMSS to confirm them (decided with
# Pablo, 5 October 2026), so the heading gives only the year.
rounds_txt <- "during 2023"

# ---- building blocks: boxes (rectangle plus centred text) and plain text cells -----------------
box <- function(x, y, w, h, label, fill = "white") {
  tibble(xmin = x - w / 2, xmax = x + w / 2, ymin = y - h / 2, ymax = y + h / 2,
         x = x, y = y, label = label, fill = fill)
}
cell <- function(x, y, label, hjust = 0.5, fontface = "plain") {
  tibble(x = x, y = y, label = label, hjust = hjust, fontface = fontface)
}

# ---- layout (x from 0 to 10.4, y from 0.4 to 14) -----------------------------------------------
# Row 1 top box; row 2 the table of rounds; row 3 the two arms; row 4 took up a place or not;
# row 5 when those who took up a place started, against their own midline survey; row 6 the
# estimation samples by wave and group.
X_INV <- 3.0; X_NEV <- 7.9                                   # the two arms
X_TOOK <- 2.2; X_NOT <- 5.9                                  # took up a place / did not
X_A <- 1.15; X_B <- 3.0; X_C <- 4.75                         # the three cases of those who took up a place
X_COL <- c(att = 4.1, inv = X_NOT, nev = X_NEV, tot = 9.45)  # columns of the sample table

b <- bind_rows(
  box(5.2, 13.45, 6.4, 0.9,
      sprintf("%d older adults expressed interest\nBaseline survey completed, %s",
              n[["expressed_interest_baseline"]], mo[["expressed_interest_baseline"]])),
  box(5.2, 11.27, 9.9, 2.3, ""),                                                       # rounds table
  box(X_INV, 8.9, 3.6, 0.8, sprintf("%d invited\nin one of the eight rounds", n[["invited"]])),
  box(X_NEV, 8.9, 2.6, 0.8, sprintf("%d never invited\n(control group)", n[["never_invited"]])),
  box(X_TOOK, 7.3, 3.9, 1.1,
      sprintf("%d took up a place\n(first entries %s)",          # no day threshold in the figure (Pablo, 30 October)
              n[["took_up_place"]], sub(" 2023 to ", " to ", mo[["took_up_place"]]))),
  box(X_NOT, 7.3, 2.6, 1.1, sprintf("%d did not take\nup a place", n[["did_not_take_up_place"]])),
  box(X_A, 5.45, 1.8, 1.45,
      sprintf("%d attending by\ntheir midline survey\n(%d of them surveyed\nagain at endline)",
              n[["took_up_attending_by_midline"]], n[["took_up_attending_by_midline_surveyed_endline"]])),
  box(X_B, 5.45, 1.75, 1.45,
      sprintf("%d entered only\nafter their\nmidline survey", n[["took_up_entered_after_midline"]])),
  box(X_C, 5.45, 1.65, 1.45,
      sprintf("%d had no\nfollow up survey", n[["took_up_no_follow_up_survey"]]), "grey95"),
  box(5.2, 2.125, 9.9, 2.35, "", "grey97")                                               # samples table
)

# ---- table of rounds: pool, invited, took up a place ---------------------------------------------
rx  <- 3.75 + (0:7) * 0.7                                    # x of rounds 1 to 8; total at 9.45
rr  <- rounds %>% filter(round != "Total") %>% mutate(round = as.integer(round)) %>% arrange(round)
tot <- rounds %>% filter(round == "Total")
ry  <- c(head = 11.75, pool = 11.3, inv = 10.85, took = 10.4)
t_rounds <- bind_rows(
  cell(5.2, 12.17, sprintf("Eight randomisation rounds %s", rounds_txt), fontface = "bold"),
  cell(0.5, ry[["head"]], "Round", hjust = 0, fontface = "bold"),
  cell(rx, ry[["head"]], as.character(rr$round), fontface = "bold"),
  cell(9.45, ry[["head"]], "Total", fontface = "bold"),
  cell(0.5, ry[["pool"]], "Not yet invited (pool)", hjust = 0),
  cell(rx, ry[["pool"]], as.character(rr$pool)),
  cell(0.5, ry[["inv"]], "Invited", hjust = 0),
  cell(rx, ry[["inv"]], as.character(rr$invited)),
  cell(9.45, ry[["inv"]], as.character(tot$invited)),
  cell(0.5, ry[["took"]], "Took up a place", hjust = 0),
  cell(rx, ry[["took"]], as.character(rr$took_up_place)),
  cell(9.45, ry[["took"]], as.character(tot$took_up_place)))

# ---- table of estimation samples ----------------------------------------------------------------
sy <- c(title = 3.05, head = 2.6, mid = 2.05, end = 1.45)
t_smp <- bind_rows(
  cell(5.2, sy[["title"]], "Estimation samples: people observed at baseline and at that wave",
       fontface = "bold"),
  cell(X_COL[["att"]], sy[["head"]], "Attending", fontface = "bold"),
  cell(X_COL[["inv"]], sy[["head"]], "Invited,\nnot attending", fontface = "bold"),
  cell(X_COL[["nev"]], sy[["head"]], "Never invited", fontface = "bold"),
  cell(X_COL[["tot"]], sy[["head"]], "Total", fontface = "bold"),
  cell(0.5, sy[["mid"]], sprintf("Midline survey\n%s", sub(" 2023 to ", " to ", mo[["midline_total"]])), hjust = 0),
  cell(X_COL, sy[["mid"]], as.character(c(n[["midline_attending"]], n[["midline_invited_not_attending"]],
                                          n[["midline_never_invited"]], n[["midline_total"]]))),
  cell(0.5, sy[["end"]], sprintf("Endline survey\n%s", mo[["endline_total"]]), hjust = 0),
  cell(X_COL, sy[["end"]], as.character(c(n[["endline_attending"]], n[["endline_invited_not_attending"]],
                                          n[["endline_never_invited"]], n[["endline_total"]]))))
# no notes inside the figure (Pablo, 30 October 2026); the people not surveyed are in the text and Table A.7
cells <- bind_rows(t_rounds, t_smp)

# ---- connectors (plain segments, as in the earlier versions) -------------------------------------
seg <- tribble(
  ~x,     ~xend,  ~y,    ~yend,
  5.2,    5.2,    13.0,  12.425,                             # top box to rounds table
  5.2,    5.2,    10.12, 9.6,                                # rounds table to the split
  X_INV,  X_NEV,  9.6,   9.6,
  X_INV,  X_INV,  9.6,   9.3,
  X_NEV,  X_NEV,  9.6,   9.3,
  X_INV,  X_INV,  8.5,   8.15,                               # invited to took up / did not
  X_TOOK, X_NOT,  8.15,  8.15,
  X_TOOK, X_TOOK, 8.15,  7.85,
  X_NOT,  X_NOT,  8.15,  7.85,
  X_TOOK, X_TOOK, 6.75,  6.45,                               # took up a place to its three cases
  X_A,    X_C,    6.45,  6.45,
  X_A,    X_A,    6.45,  6.175,
  X_B,    X_B,    6.45,  6.175,
  X_C,    X_C,    6.45,  6.175,
  X_A,    X_A,    4.725, 3.85,                               # attending by midline: Attending column
  X_A,    X_COL[["att"]], 3.85, 3.85,
  X_COL[["att"]], X_COL[["att"]], 3.85, 3.3,
  X_B,    X_B,    4.725, 4.25,                               # entered after midline: not attending
  X_B,    X_NOT,  4.25,  4.25,
  X_NOT,  X_NOT,  6.75,  3.3,                                # did not take up a place
  X_NEV,  X_NEV,  8.5,   3.3)                                # never invited
# the "no follow up survey" box has no connector below it: those two are in no sample

figA <- ggplot() +
  geom_segment(data = seg, aes(x = x, xend = xend, y = y, yend = yend),
               linewidth = 0.4, colour = "grey30") +
  geom_rect(data = b, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
            colour = "grey20", linewidth = 0.4) +
  geom_text(data = b %>% filter(label != ""), aes(x = x, y = y, label = label),
            size = 3.1, lineheight = 1.05) +
  # rule under the header row of each table
  annotate("segment", x = 0.45, xend = 9.95, y = ry[["head"]] - 0.2, yend = ry[["head"]] - 0.2,
           linewidth = 0.3, colour = "grey50") +
  annotate("segment", x = 0.45, xend = 9.95, y = sy[["head"]] - 0.3, yend = sy[["head"]] - 0.3,
           linewidth = 0.3, colour = "grey50") +
  geom_text(data = cells, aes(x = x, y = y, label = label, hjust = hjust, fontface = fontface),
            size = 3.1, lineheight = 1.05) +
  scale_fill_identity() +
  coord_cartesian(xlim = c(0.1, 10.3), ylim = c(0.9, 13.95), expand = FALSE) +
  theme_void() +
  theme(plot.margin = margin(10, 10, 10, 10))
ggsave(paste0(OUT_FIG, "figure_flow.png"), figA, width = 9.5, height = 10.75, dpi = 300, bg = "white",
       device = ragg::agg_png)
cat("figure_flow.png written\n")

# =============================================================================================
# Figure B. Timing of the surveys, as a schematic
# =============================================================================================
# Redrawn on 5 October 2026 at Pablo's request. The earlier version (two clocks, group means with
# standard deviations, two real example participants) is in the _previous subfolders, suffix
# _before_schematic. Nothing below is a real person or a computed timing: the positions only
# illustrate the design, which is what the figure has to convey:
#   - the baseline survey ran from October 2022 to January 2023, close to the opening of the
#     centre (month 0 = January 2023);
#   - for those attending, the follow ups came about 3 and 9 months after their own ENTRY; two
#     example rows (entry in January and in February 2023) show that the clock moves with entry;
#   - those not attending (invited who declined, and never invited) were surveyed in the same
#     calendar windows as those attending.
# The axis shows month NUMBERS only. Month names would make the schematic windows read as the real
# ones, which were wider (midline April to September 2023, endline October 2023 to February 2024).
# The facts the schematic stands for are recomputed from the data at the end of this section and
# printed to the console only (group means; no individual dates).

ev_lev  <- c("Baseline survey", "Entry to the day centre", "Midline survey", "Endline survey")
row_lev <- c("Attending: entered in\nJanuary 2023 (example)",
             "Attending: entered in\nFebruary 2023 (example)",
             "Not attending (invited who\ndeclined, and never invited)")
lab3 <- "about 3 months\nafter entry"
lab9 <- "about 9 months\nafter entry"
labw <- "same calendar\nwindow"

# ---- schematic positions: months since the centre opened; lab = text above the marker ---------
# Row r is drawn at height y = 4 - r (row 1 on top).
sch <- tribble(
  ~r, ~ev,                         ~x,   ~lab,
  1,  "Baseline survey",           -0.3, NA,     # just before the opening (all rows alike)
  1,  "Entry to the day centre",    0,   NA,
  1,  "Midline survey",             3,   lab3,   # entry + 3
  1,  "Endline survey",             9,   lab9,   # entry + 9
  2,  "Baseline survey",           -0.3, NA,
  2,  "Entry to the day centre",    1,   NA,     # entered one month later ...
  2,  "Midline survey",             4,   lab3,   # ... so both follow ups move one month later
  2,  "Endline survey",            10,   lab9,
  3,  "Baseline survey",           -0.3, NA,     # no entry; surveyed in the same windows
  3,  "Midline survey",             3.5, labw,
  3,  "Endline survey",             9.5, labw) %>%
  mutate(y = 4 - r, ev = factor(ev, levels = ev_lev))

# ---- survey windows: each covers its months on the axis plus or minus 0.4 ----------------------
# Baseline: months -3 to 0 (October 2022 to January 2023). Midline and endline: the months of the
# two attending examples (3 and 4; 9 and 10), which the "not attending" row shares.
# Baseline band as narrow as the others, just before the opening (Pablo, 30 October 2026); the real
# baseline ran from October 2022 to January 2023, as the text says.
win_b <- tibble(xmin = c(-1.2, 2.6, 8.6), xmax = c(0.4, 4.4, 10.4),
                lab  = c("Baseline\nsurveys", "Midline\nsurveys", "Endline\nsurveys"))
xb <- sch$x[sch$ev == "Baseline survey"]
stopifnot(all(xb > -1.2 & xb < 0.4),
          all(between(sch$x[sch$ev == "Midline survey"], 2.6, 4.4)),
          all(between(sch$x[sch$ev == "Endline survey"], 8.6, 10.4)))

Y_TOP <- 3.62                                                 # baseline of the window labels
figB <- ggplot(sch, aes(x = x, y = y)) +
  # survey windows, translucent so the monthly grid shows through
  geom_rect(data = win_b, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "grey50", alpha = 0.13) +
  geom_text(data = win_b, aes(x = (xmin + xmax) / 2, y = Y_TOP, label = lab), inherit.aes = FALSE,
            vjust = 0, size = 2.8, lineheight = 0.95, colour = "grey25") +
  # opening of the centre
  # (the dashed line stops below the window labels so it does not cross them)
  annotate("segment", x = 0, xend = 0, y = -Inf, yend = Y_TOP - 0.08, linetype = "dashed",
           colour = "grey45", linewidth = 0.4) +
  annotate("text", x = 0.15, y = 0.52, label = "Centre opens", hjust = 0, size = 2.8,
           colour = "grey25") +
  # each row: its path, its markers, and the labels above the follow ups
  geom_line(aes(group = y), linewidth = 0.5, colour = "grey60") +
  geom_point(aes(shape = ev, fill = ev), size = 3.1, colour = "grey15") +
  geom_text(data = sch %>% filter(!is.na(lab)), aes(y = y + 0.17, label = lab),
            vjust = 0, size = 2.8, lineheight = 0.95, colour = "grey20") +
  scale_shape_manual(values = c(21, 24, 22, 23), name = NULL, drop = FALSE) +
  scale_fill_manual(values = c("white", "grey20", "grey70", "grey45"), name = NULL, drop = FALSE) +
  # one gridline per month; true minus sign on negative months, so no label reads as a dash
  scale_x_continuous("Months since the centre opened in January 2023", breaks = -1:12,
                     labels = function(b) gsub("-", "−", b, fixed = TRUE)) +
  scale_y_continuous(NULL, breaks = 3:1, labels = row_lev) +
  coord_cartesian(xlim = c(-1.4, 12.6), ylim = c(0.38, 4.3), expand = FALSE) +
  theme_minimal(base_size = 9.5) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.3),
        axis.text.y = element_text(colour = "grey15", lineheight = 0.95),
        legend.position = "top", legend.margin = margin(0, 0, 0, 0),
        plot.caption = element_text(hjust = 0, size = 7.5, colour = "grey30"),
        plot.caption.position = "plot",
        plot.margin = margin(4, 10, 4, 4))
ggsave(paste0(OUT_FIG, "figure_timeline.png"), figB, width = 6.3, height = 3.9, dpi = 300, bg = "white",
       device = ragg::agg_png)
cat("figure_timeline.png written (schematic)\n\n")

# ---- console check: the facts the schematic stands for (group means; nothing is plotted) -------
# Estimation samples of each wave (as in 03). Attending = more than five days by that wave's survey.
# Opening = first day anyone attended (earliest entry date in the register).
mo_ <- function(a, b) as.numeric(b - a) / 30.44              # months between two dates
s1 <- ancova_sample(dat, "whoqol", 1) %>% select(id, d) %>% mutate(wave = 1)
s2 <- ancova_sample(dat, "whoqol", 2) %>% select(id, d) %>% mutate(wave = 2)
tl <- bind_rows(s1, s2) %>%
  left_join(dat %>% select(id, z, date0, date1, date2, entry_date), by = "id") %>%
  mutate(date_w = if_else(wave == 1, date1, date2), grp = grp_of(z, d))
stopifnot(all(!is.na(tl$entry_date[tl$d == 1])))             # every attender has an entry date
open_date <- min(dat$entry_date, na.rm = TRUE)
cat("Opening (first day anyone attended):", format(open_date, "%B %Y"), "\n")
cat("Baseline survey window:",
    paste(format(range(dat$date0, na.rm = TRUE), "%B %Y"), collapse = " to "), "\n\n")
chk <- tl %>%
  mutate(wave = if_else(wave == 1, "midline", "endline")) %>%
  group_by(wave, grp) %>%
  summarise(n = n(),
            months_from_opening = round(mean(mo_(open_date, date_w)), 1),
            months_from_entry   = round(mean(mo_(entry_date, date_w)), 1),   # attending only
            .groups = "drop") %>%
  mutate(months_from_entry = if_else(grp == "Attending", months_from_entry, NA_real_)) %>%
  arrange(desc(wave), grp)
print(as.data.frame(chk), row.names = FALSE)
