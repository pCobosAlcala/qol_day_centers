# 10_dose_table.R  (revision, October 2026)
#
# Annex table for the "Length and frequency of use" subsection: the endline sample (85) split by
# what people did, not by what was randomised. Descriptive only: attendance and its intensity were
# chosen, so these groups are not comparable the way invited and never invited are.
#
#   Never invited | Invited, did not attend | Attended 100 days or fewer | Attended more than 100 days
#   (attending = more than five days at the centre by the endline survey; days counted up to that
#    survey from the attendance register; the 100 day cut is the one 08_text_numbers.R uses for the
#    Lawton comparison quoted in the text)
#
# For each group: n; days, weeks and days per week attended by endline; baseline autonomy (Lawton);
# WHOQOL-BREF at baseline and endline (WHO 0 to 100 scale); change in per cent of the group's own
# baseline and in points; and the endline difference from the never invited controlling for the
# baseline score (ANCOVA, robust HC2 errors), which is descriptive.
#
# Output: 3_code/output/tables/tableA_dose.png and 3_code/output/csv/tableA_dose.csv (aggregates).
# Run from the imss_paper folder.

rm(list = ls())
suppressPackageStartupMessages({
  library(tidyverse); library(sandwich); library(lmtest)
  library(gridExtra); library(gtable); library(grid)
})
options(scipen = 999, width = 200)
source("3_code/scripts/00_paths.R")
source("3_code/scripts/00_data.R")

# ---- 1. the endline sample and the four groups -------------------------------------------------
dat <- build_analysis_data()
e <- dat %>%
  filter(!is.na(whoqol_y0), !is.na(whoqol_y2)) %>%
  mutate(group = case_when(
    z == 0                    ~ "Never invited",
    d2 == 0                   ~ "Invited, did not attend",
    d2 == 1 & days2 <= 100    ~ "Attended 100 days or fewer",
    d2 == 1 & days2 >  100    ~ "Attended more than 100 days"),
    group = factor(group, levels = c("Never invited", "Invited, did not attend",
                                     "Attended 100 days or fewer", "Attended more than 100 days")))
stopifnot(nrow(e) == 85, !any(is.na(e$group)), !any(e$days2 == 100 & e$d2 == 1, na.rm = TRUE))

# ---- 2. group summaries ------------------------------------------------------------------------
summ <- e %>% group_by(group) %>% summarise(
  n = n(),
  days = mean(days2), days_min = min(days2), days_max = max(days2),
  weeks = mean(weeks2), dpw = mean(dpw2),
  lawton = mean(lawton0, na.rm = TRUE),
  y0 = mean(whoqol_y0), y2 = mean(whoqol_y2), .groups = "drop") %>%
  mutate(change_pts = y2 - y0, change_pct = 100 * change_pts / y0)

# descriptive ANCOVA: endline on baseline plus group dummies, never invited as reference
m  <- lm(whoqol_y2 ~ whoqol_y0 + group, data = e)
ct <- coeftest(m, vcov. = vcovHC(m, type = "HC2"))
ctrl_base <- summ$y0[summ$group == "Never invited"]
diffs <- tibble(group = levels(e$group)[-1],
                diff = ct[paste0("group", levels(e$group)[-1]), 1],
                se   = ct[paste0("group", levels(e$group)[-1]), 2],
                p    = ct[paste0("group", levels(e$group)[-1]), 4]) %>%
  mutate(diff_pct = 100 * diff / ctrl_base)
out <- summ %>% mutate(group = as.character(group)) %>% left_join(diffs, by = "group")
write_csv(out, paste0(OUT_CSV, "tableA_dose.csv"))

# ---- 3. the PNG, in the style of 03 ------------------------------------------------------------
render <- function(df, file, title, notes, note_width = 120) {
  th <- ttheme_minimal(
    base_size = 10,
    core    = list(fg_params = list(fontsize = 10, lineheight = 1.05), padding = unit(c(5, 3.2), "mm")),
    colhead = list(fg_params = list(fontsize = 10, fontface = "bold", lineheight = 1.05),
                   padding = unit(c(5, 3), "mm")))
  g <- tableGrob(df, rows = NULL, theme = th)
  for (i in which(g$layout$name %in% c("core-fg", "colhead-fg") & g$layout$l == 1)) {
    g$grobs[[i]]$x <- unit(0.03, "npc"); g$grobs[[i]]$hjust <- 0
  }
  nc <- ncol(g)
  rule <- function(y) segmentsGrob(y0 = unit(y, "npc"), y1 = unit(y, "npc"), gp = gpar(lwd = 1.2))
  g <- gtable_add_grob(g, rule(1), t = 1, l = 1, r = nc)
  g <- gtable_add_grob(g, rule(0), t = 1, l = 1, r = nc)
  g <- gtable_add_grob(g, rule(0), t = nrow(g), l = 1, r = nc)
  # title and notes inside the image only if IMAGE_TEXT (00_paths.R); otherwise saved as text
  if (IMAGE_TEXT) {
    g <- gtable_add_rows(g, heights = unit(10, "mm"), pos = 0)
    g <- gtable_add_grob(g, textGrob(title, x = unit(0, "npc"), hjust = 0,
                         gp = gpar(fontface = "bold", fontsize = 11.5)), t = 1, l = 1, r = nc)
    pdf(NULL); tw <- convertWidth(sum(g$widths), "in", valueOnly = TRUE); invisible(dev.off())
    nt <- paste(strwrap(notes, width = min(note_width, floor(tw * 15.5))), collapse = "\n")
    nlines <- length(strsplit(nt, "\n")[[1]])
    g <- gtable_add_rows(g, heights = unit(4 + 4.2 * nlines, "mm"))
    g <- gtable_add_grob(g, textGrob(nt, x = unit(0, "npc"), y = unit(1, "npc") - unit(2.5, "mm"),
                         hjust = 0, vjust = 1, gp = gpar(fontsize = 8.5, col = "grey20", lineheight = 1.1)),
                         t = nrow(g), l = 1, r = nc)
  } else save_note(file, title, notes)
  pdf(NULL)
  wd <- convertWidth(sum(g$widths), "in", valueOnly = TRUE) + 0.5
  ht <- convertHeight(sum(g$heights), "in", valueOnly = TRUE) + 0.4
  invisible(dev.off())
  ragg::agg_png(file, width = wd, height = ht, units = "in", res = 300, background = "white")
  grid.newpage(); pushViewport(viewport(width = unit(wd - 0.5, "in"), height = unit(ht - 0.4, "in")))
  grid.draw(g); invisible(dev.off())
  cat(sprintf("written %s (%.1f x %.1f in)\n", file, wd, ht))
}

mns <- function(x) ifelse(x < 0, paste0("−", formatC(abs(x), format = "f", digits = 1)),
                          formatC(x, format = "f", digits = 1))           # true minus sign, never a dash
show <- out %>% transmute(
  Group = group,
  `n` = as.character(n),
  # people not attending are shown as such (the text never details their few visits, Pablo's rule)
  `Days attended\nby endline` = case_when(group == "Never invited" ~ "0",
                                          group == "Invited, did not attend" ~ "Not attending",
                                          TRUE ~ sprintf("%.0f\n(%d to %d)", days, days_min, days_max)),
  `Days per\nweek` = case_when(group == "Never invited" ~ "0",
                               group == "Invited, did not attend" ~ "Not attending",
                               TRUE ~ sprintf("%.1f", dpw)),
  `Lawton at\nbaseline (0 to 8)` = sprintf("%.1f", lawton),
  `WHOQOL-BREF\nbaseline` = sprintf("%.1f", y0),
  `WHOQOL-BREF\nendline` = sprintf("%.1f", y2),
  `Change, per cent\nof own baseline\n(points)` = sprintf("%s%%\n(%s)", mns(change_pct), mns(change_pts)),
  `Difference from never\ninvited at endline,\nper cent (points), p` =
    ifelse(is.na(diff), "Reference",
           sprintf("%s%%\n(%s), p = %s", mns(diff_pct), mns(diff),
                   ifelse(p < 0.01, sprintf("%.3f", p), sprintf("%.2f", p)))))

render(show, paste0(OUT_TAB, "tableA_dose.png"),
  "Attendance and quality of life at endline, by what people did (descriptive)",
  paste("Endline sample (85 people). Attending means more than five days at the centre by the endline",
        "survey; days are counted up to that survey from the attendance register. WHOQOL-BREF on the WHO",
        "0 to 100 scale. The last column comes from a regression of the endline score on the baseline score",
        "and the four groups (robust HC2 errors), in per cent of the never invited group's baseline score.",
        "Only the invitation was randomised: attendance and its intensity were chosen, so the groups are",
        "not comparable and these differences are descriptive, not effects."))
print(out, width = 200)
