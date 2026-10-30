# 09_environment_items.R  (revision, October 2026)
#
# WHICH ITEMS MOVE THE ENVIRONMENTAL DOMAIN? (Reviewer 2)
#
# The environmental domain is the only part of the WHOQOL-BREF with a clear effect at endline.
# Reviewer 2 asks why. This script runs the main specification of 03_main_results.R on each of
# the domain's eight items, one at a time:
#   ANCOVA run separately per wave:   Y_iw = a0 + a1*Y_i0 + a2*Z_i + e_i
#   ITT = effect of the invitation; IV = effect of attending (more than five days at the centre by
#   that wave's survey) instrumented by the invitation, with the item's baseline value in both
#   stages. Robust HC2 errors, one row per person, no covariates. Each item on the WHO 0 to 100
#   scale, 25 x (answer minus 1), so a single item takes the values 0, 25, 50, 75 or 100.
#   Holm adjusted p values over the eight items, separately for each wave and estimator.
#
# The eight environmental items (WHO number, question of the Mexican form, as mapped in 00_data.R
# and checked against the question titles of cuestionarios_CD_imss.xlsx in October 2026):
#    8  16.09  how safe do you feel in your daily life
#    9  16.10  how healthy is your physical environment (pollution, noise, traffic, climate)
#   12  16.13  have you enough money to meet your needs
#   13  16.14  how available is the information you need in your daily life
#   14  16.15  to what extent do you have the opportunity for leisure (recreational) activities
#   23  16.25  how satisfied are you with the conditions of your living place
#   24  16.26  how satisfied are you with your access to health services
#   25  16.27  how satisfied are you with your transport
# None is reversed in the WHO scoring.
#
# Data: build_analysis_data() in 00_data.R gives the domain scores, the invitation z and the wave
# attendance d1, d2, but not single items. The item answers are read here from the same raw file
# (pm_limpieza1.csv) with the core's own item map (WHOQOL_ITEMS) and answer labels
# (WHOQOL_RESPONSES). CHECK: the mean of the eight item scores reproduces the core's environmental
# domain for every person and wave, and the domain row of the table reproduces main_results.csv.
#
# Outputs (aggregates only): 3_code/output/tables/tableA_environment_items.png and
# 3_code/output/csv/tableA_environment_items.csv. Run from the imss_paper folder:
#   "C:/Program Files/R/R-4.6.1/bin/Rscript.exe" 3_code/scripts/09_environment_items.R

rm(list = ls())
suppressPackageStartupMessages({
  library(tidyverse); library(AER); library(sandwich); library(lmtest)
  library(gridExtra); library(gtable); library(grid)
})
options(scipen = 999, width = 200)

source("3_code/scripts/00_paths.R")
source("3_code/scripts/00_data.R")

# ---- 1. data -------------------------------------------------------------------------------
dat <- build_analysis_data()                 # one row per randomised person; checks run inside

# the eight environmental items, taken from the core's map, with short labels for the table
env <- WHOQOL_ITEMS %>% filter(part == "environmental") %>% arrange(who)
labels <- c(`8`  = "Safety in daily life",
            `9`  = "Healthy physical environment",
            `12` = "Enough money for needs",
            `13` = "Information for daily life",
            `14` = "Opportunity for leisure",
            `23` = "Home (living place conditions)",
            `24` = "Access to health services",
            `25` = "Transport")
stopifnot(nrow(env) == 8, identical(as.character(env$who), names(labels)), !any(env$reverse))
env <- env %>% mutate(label    = unname(labels[as.character(who)]),
                      question = sub("^x16_", "16.", src))          # x16_09 -> 16.09

# item answers, long: one row per person, item and wave (the same file and reader as the core)
lim <- .read_raw_csv(RAW, "pm_limpieza1.csv") %>%
  select(id = x03_01, name, t0, t1, t2) %>%
  filter(name %in% env$src)
stopifnot(!anyDuplicated(lim[, c("id", "name")]))
ans <- lim %>%
  pivot_longer(c(t0, t1, t2), names_to = "wave", values_to = "label") %>%
  transmute(id, src = name, wave = as.integer(sub("t", "", wave)), label = .ascii(label))
unknown <- setdiff(unique(na.omit(ans$label)), names(WHOQOL_RESPONSES))
if (length(unknown)) stop("unmapped answer label(s): ", paste(unknown, collapse = " | "))

# keep the surveys that took place (the person-waves with a domain score in the core) and put
# every item on the WHO 0 to 100 scale
dom <- dat %>%
  select(id, starts_with("environmental_y")) %>%
  pivot_longer(-id, names_to = "wave", values_to = "domain") %>%
  mutate(wave = as.integer(sub("environmental_y", "", wave))) %>%
  filter(!is.na(domain))
items <- dom %>%
  tidyr::crossing(src = env$src) %>%
  left_join(ans, by = c("id", "wave", "src")) %>%
  mutate(answer = unname(WHOQOL_RESPONSES[label]),           # 1 to 5
         score  = 25 * (answer - 1))                         # 0 to 100

# CHECK 1: no missing item in a survey held, and the item mean is the core's domain score
stopifnot(!anyNA(items$score))
chk <- items %>% group_by(id, wave) %>%
  summarise(item_mean = mean(score), domain = domain[1], .groups = "drop")
stopifnot(max(abs(chk$item_mean - chk$domain)) < 1e-9)
cat(sprintf("\nCHECK 1 passed: the mean of the eight items reproduces the environmental domain for all %d person-waves (largest gap %.1e).\n",
            nrow(chk), max(abs(chk$item_mean - chk$domain))))

# wide by item and wave, joined to the invitation and the wave attendance
wide <- items %>%
  mutate(col = paste0(src, "_y", wave)) %>%
  select(id, col, score) %>%
  pivot_wider(names_from = col, values_from = score) %>%
  right_join(dat %>% select(id, z, d1, d2), by = "id")
stopifnot(nrow(wide) == 120)

# ---- 2. estimation (same functions as 03_main_results.R) -------------------------------------
rob <- function(m) {                                     # HC2; HC1 only if HC2 fails, and say so
  tryCatch(vcovHC(m, type = "HC2"),
           error = function(e) { warning("HC2 failed, HC1 used: ", conditionMessage(e))
                                 vcovHC(m, type = "HC1") })
}
pull_est <- function(m, term) {
  V  <- rob(m); ct <- coeftest(m, vcov. = V); ci <- coefci(m, vcov. = V, level = 0.95)
  tibble(est = ct[term, 1], se = ct[term, 2], p = ct[term, 4], lo = ci[term, 1], hi = ci[term, 2])
}

# ITT, IV, first stage and magnitudes for one ANCOVA sample s (columns z, d, y0, y)
fit_one <- function(s) {
  itt <- lm(y ~ y0 + z, s)
  iv  <- ivreg(y ~ d + y0 | z + y0, data = s)
  fz  <- pull_est(lm(d ~ y0 + z, s), "z")
  c0  <- s %>% filter(z == 0)
  base <- mean(c0$y0); sd0 <- sd(s$y0)                   # control baseline mean; baseline SD, both arms
  bind_rows(pull_est(itt, "z") %>% mutate(model = "ITT"),
            pull_est(iv,  "d") %>% mutate(model = "IV")) %>%
    mutate(n = nrow(s), fsF = (fz$est / fz$se)^2,
           control_base       = base,
           sd0                = sd0,
           control_change_pts = mean(c0$y) - base,
           control_change_pct = 100 * control_change_pts / base,
           control_change_sd  = control_change_pts / sd0,
           pct                = 100 * est / base,          # estimate as per cent of control baseline
           sd_units           = est / sd0)                 # estimate in baseline SD
}

res <- list()
for (k in seq_len(nrow(env))) {
  for (wv in 1:2) {
    y0 <- paste0(env$src[k], "_y0"); yw <- paste0(env$src[k], "_y", wv)
    s  <- wide %>% filter(!is.na(.data[[y0]]), !is.na(.data[[yw]])) %>%
      transmute(id, z, d = .data[[paste0("d", wv)]], y0 = .data[[y0]], y = .data[[yw]])
    res[[paste(k, wv)]] <- fit_one(s) %>%
      mutate(who_item = env$who[k], question = env$question[k], item = env$label[k], wave = wv)
  }
}
# the domain itself, for reference (the mean of the eight items; not part of the Holm family)
for (wv in 1:2) {
  res[[paste("domain", wv)]] <- fit_one(ancova_sample(dat, "environmental", wv)) %>%
    mutate(who_item = NA_real_, question = NA_character_,
           item = "Environmental domain (mean of the eight items)", wave = wv)
}
res <- bind_rows(res)

# Holm adjustment over the eight items, separately for each wave and estimator
res <- res %>%
  group_by(wave, model, is_item = !is.na(who_item)) %>%
  mutate(p_holm = if (is_item[1]) p.adjust(p, method = "holm") else NA_real_) %>%
  ungroup() %>% select(-is_item)
stopifnot(all(count(res %>% filter(!is.na(who_item)), wave, model)$n == 8))

# samples: with no missing item, every item uses the domain's sample (102 midline, 85 endline)
stopifnot(all(res$n[res$wave == 1] == 102), all(res$n[res$wave == 2] == 85))

# CHECK 2: the domain row reproduces 03's main results (if 03 has been run)
mr_file <- paste0(OUT_CSV, "main_results.csv")
if (file.exists(mr_file)) {
  mr <- read_csv(mr_file, show_col_types = FALSE) %>%
    filter(outcome == OUTCOMES[["environmental"]]) %>%
    select(wave, model, est03 = est, se03 = se, p03 = p)
  cmp <- res %>% filter(is.na(who_item)) %>% inner_join(mr, by = c("wave", "model"))
  stopifnot(nrow(cmp) == 4,
            max(abs(c(cmp$est - cmp$est03, cmp$se - cmp$se03, cmp$p - cmp$p03))) < 1e-9)
  cat("CHECK 2 passed: the domain row reproduces main_results.csv of 03.\n")
} else {
  cat("CHECK 2 skipped: main_results.csv not found (run 03 first).\n")
}

# ---- 3. save and print -----------------------------------------------------------------------
out <- res %>%
  select(item, who_item, question, wave, model, n, est, se, p, p_holm, lo, hi, fsF,
         control_base, sd0, control_change_pts, control_change_pct, control_change_sd,
         pct, sd_units) %>%
  arrange(wave, model, is.na(who_item), who_item)
write_csv(out, paste0(OUT_CSV, "tableA_environment_items.csv"))
cat("written", paste0(OUT_CSV, "tableA_environment_items.csv"), "\n\n")

print(out %>% mutate(across(c(est, se, p, p_holm, lo, hi, fsF, control_base, sd0, control_change_pts,
                              control_change_pct, pct, sd_units), ~ round(.x, 3))) %>%
        select(item, wave, model, est, se, p, p_holm, lo, hi, fsF, control_base, sd0,
               control_change_pts, control_change_pct, pct, sd_units) %>% as.data.frame(),
      row.names = FALSE)

# how far the eight item estimates are from the domain estimate when averaged (they need not
# match exactly: each item regression controls for the item's own baseline, not the domain's)
avg <- res %>% group_by(wave, model) %>%
  summarise(mean_of_items = mean(est[!is.na(who_item)]), domain = est[is.na(who_item)], .groups = "drop")
cat("\nmean of the eight item estimates against the domain estimate:\n")
print(avg %>% mutate(across(where(is.double), ~ round(.x, 3))) %>% as.data.frame(), row.names = FALSE)

# ---- 4. rendering helper (the PNG table style of 03_main_results.R) --------------------------
star <- function(p) ifelse(p < .01, "***", ifelse(p < .05, "**", ifelse(p < .1, "*", "")))

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

# ---- 5. Table A: the eight items -------------------------------------------------------------
# cell: estimate in points with stars (unadjusted p), HC2 SE, unadjusted p, Holm p (items only)
fp <- function(p) ifelse(p < 0.001, "< 0.001", sprintf("= %.3f", p))   # never print "p = 0.000"
cell <- function(r) {
  top <- sprintf("%.2f%s\n(%.2f)\np %s", r$est, star(r$p), r$se, fp(r$p))
  if (is.na(r$p_holm)) top else sprintf("%s\nHolm p %s", top, fp(r$p_holm))
}
row_label <- function(k) {
  if (is.na(k$who_item)) "Environmental domain\n(mean of the eight items)"
  else sprintf("%s\nWHO item %d, question %s", k$item, k$who_item, k$question)
}
keys <- res %>% distinct(item, who_item, question) %>% arrange(is.na(who_item), who_item)
tab <- map_dfr(seq_len(nrow(keys)), function(i) {
  k <- keys[i, ]
  g <- function(w, m) cell(res %>% filter(item == k$item, wave == w, model == m))
  b <- function(w) sprintf("%.1f", res %>% filter(item == k$item, wave == w, model == "ITT") %>% pull(control_base))
  tibble(Item = row_label(k),
         `Control\nbaseline` = b(1), `Attending\n(IV)` = g(1, "IV"), `Invited\n(ITT)` = g(1, "ITT"),
         `Control \nbaseline` = b(2), `Attending \n(IV)` = g(2, "IV"), `Invited \n(ITT)` = g(2, "ITT"))
})
fr <- function(w) { f <- res %>% filter(wave == w, model == "IV") %>% pull(fsF); sprintf("%.1f to %.1f", min(f), max(f)) }
n1 <- unique(res$n[res$wave == 1]); n2 <- unique(res$n[res$wave == 2])
tab <- bind_rows(tab,
  tibble(Item = "First-stage F", `Control\nbaseline` = "", `Attending\n(IV)` = fr(1), `Invited\n(ITT)` = "",
         `Control \nbaseline` = "", `Attending \n(IV)` = fr(2), `Invited \n(ITT)` = ""),
  tibble(Item = "Individuals", `Control\nbaseline` = "", `Attending\n(IV)` = as.character(n1),
         `Invited\n(ITT)` = as.character(n1), `Control \nbaseline` = "",
         `Attending \n(IV)` = as.character(n2), `Invited \n(ITT)` = as.character(n2)))

render(tab,
  spans = list(list(label = "Midline\n(about 3 months after entry)", l = 2, r = 4),
               list(label = "Endline\n(about 9 months after entry)", l = 5, r = 7)),
  file  = paste0(OUT_TAB, "tableA_environment_items.png"),
  title = "Effects on the eight items of the environmental domain (WHOQOL-BREF)",
  notes = paste(
    "Item level estimates are exploratory: each item is a single answer on a five point scale, and eight",
    "items are tested at each wave. Each cell reports the estimate in points of the WHO 0 to 100 scale",
    "(25 times the answer minus 1, so an item takes the values 0, 25, 50, 75 or 100), the robust (HC2)",
    "standard error in parentheses, the unadjusted p value and the Holm adjusted p value. The Holm adjustment",
    "is over the eight items, separately for each wave and for ITT and IV. Each item is estimated as in the",
    "main results (Table 3): ANCOVA of the item at each wave on its own baseline value and the randomised",
    "invitation, separately for each wave. ITT is the effect of the invitation. IV is the effect of attending,",
    "instrumented by the invitation with the baseline value in both stages; attending means more than five",
    "days at the centre by that wave's survey (months are counted from entry for those attending). Control",
    "baseline is the item's mean at baseline among those never invited, in that wave's sample. The last row",
    "is the domain estimate of Table 3",
    "(the mean of the eight items), which is not part of the adjustment. First-stage F is the robust F on",
    "the excluded instrument, shown as the range across the nine regressions. Items are numbered as in the",
    "WHO form and in section 16 of the study questionnaire. Stars use the unadjusted p value:",
    "*** p < 0.01, ** p < 0.05, * p < 0.1."), note_width = 125)
