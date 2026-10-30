# Libraries and setup ----
pacman::p_load(tidyverse, AER, broom, broom.mixed, lme4, lmerTest,
               kableExtra, parallel)
options(scipen = 999)
rm(list = ls())

# Viz----
data <- read_csv("C://Users/pablo/OneDrive/Documentos/pablo/RStudio/code/2023/imss/imss/1_data/bases/base_paper/gerontologist/db.csv")

db = data %>%
  mutate(var = case_when(var == "whoqol" ~ "WHOQOL-BREF",
                         T ~ var)) %>%
  group_by(x03_01) %>%
  mutate(value_ = mean(mean_previos, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    z = case_when(
      z == 0 ~ 1,
      z == 1 & value_ == 0 ~ 2,
      z == 1 & value_ > 0  ~ 3,
      TRUE ~ NA_real_
    ),
    z = factor(
      z,
      levels = c(1, 2, 3),
      labels = c("Waiting\nlist", "Invited\nNot attended", "Attended")
    ),
    time = factor(
      time,
      levels = c(0, 1, 2),
      labels = c("Baseline", "Short-\nterm", "Long-\nterm")
    ),
    var = factor(
      var,
      levels = c(
        "WHOQOL-BREF",
        "General",
        "Physical",
        "Psychological",
        "Social",
        "Enviromental"
      )
    )
  ) %>%
  filter(!is.na(z), !is.na(time))

# Aggregated means
agg <- db %>%
  group_by(var, z, time) %>%
  summarise(
    val = mean(value, na.rm = TRUE),
    n   = sum(!is.na(value)),
    sd  = sd(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(se = sd / sqrt(n))

# Plot
p = ggplot(db, aes(x = time, y = value)) +
  geom_boxplot(width = 0.55, linewidth = 0.7,
               alpha = 0.10, outlier.alpha = 0) +
  geom_line(data = agg, aes(y = val, group = 1),
            linewidth = 1.6, color = "black") +
  geom_point(data = agg, aes(y = val),
             size = 3.2, color = "black") +
  facet_grid(
    rows = vars(var),
    cols = vars(z),
    switch = "y"
  ) +
  theme_minimal() +
  theme(
    panel.spacing = unit(0.2, "lines"),
    strip.text.y.left = element_text(angle = 0),
    strip.text.x = element_text(margin = margin(b = 4)),
    plot.margin = margin(5, 5, 5, 5),
    axis.title.x = element_blank()
  ) +
  labs(y = NULL)

p

# Save
ggsave(
  "3_viz/paper/psi4.jpg",
  plot = p,
  width = 8,
  height = 10,
  dpi = 1000
)

ggsave(
  "3_viz/paper/figure1.tif",
  plot = p,
  width = 8,
  height = 10,
  dpi = 100,
  device = "tiff"
)


# Econometric analysis----
# Data
data <- data %>%
  mutate(across(c(d, z), ~ case_when(time == 0 ~ 0, TRUE ~ .x)),
         z = as.factor(z),
         d = as.factor(d))

# Variables
deps <- c("Enviromental", "Psychological", "Social", "Physical", "General", "whoqol")
inds <- c("z", "value_previos", "value_months", "mean_previos")
ctrls <- c("x04_01", "x05_01", "x05_02", "x05_05", "x05_07_", "x06_03_pc_")
ivs <- c("d", "value_previos", "value_months", "mean_previos")

# IV model
results_iv <- list()
for (dv in deps) {
  df <- data %>% filter(var == dv)
  for (iv in ivs) {
    instr <- "z"  # Instrument
    # IV model without sociodemographic controls
    f_iv1 <- as.formula(paste("value ~ time +", iv, "+ time:", iv, "+ factor(x03_01) | ",
                              "time +", instr, "+ time:", instr, "+ factor(x03_01)"))
    m_iv1 <- ivreg(f_iv1, data = df)
    t_iv1 <- broom::tidy(m_iv1) %>%
      mutate(dv = dv, iv = iv, controls = "No", model = "IV")

    # IV model with sociodemographic controls
    f_iv2 <- as.formula(paste("value ~ time +", iv, "+ time:", iv, "+",
                              paste(ctrls, collapse = "+"), "+ factor(x03_01) | ",
                              "time +", instr, "+ time:", instr, "+",
                              paste(ctrls, collapse = "+"), "+ factor(x03_01)"))
    m_iv2 <- ivreg(f_iv2, data = df)
    t_iv2 <- broom::tidy(m_iv2) %>%
      mutate(dv = dv, iv = iv, controls = "Yes", model = "IV")

    results_iv[[paste(dv, iv, "iv", sep = "_")]] <- bind_rows(t_iv1, t_iv2)
  }
}

results_df_iv <- bind_rows(results_iv)

# OLS
fit_model <- function(formula, data, binary) {
  if (binary) {
    model <- glmer(formula, data = data, family = binomial())
  } else {
    model <- lmerTest::lmer(formula, data = data)
  }
  broom.mixed::tidy(model)
}

results <- list()
for (dv in deps) {
  df <- data %>% filter(var == dv)
  binary <- dv %in% c("x19_02", "x21_01")
  for (iv in inds) {

    # Random effects, no controls
    f1 <- as.formula(paste("value ~ time +", iv, "+ time:", iv, "+ (1 | x03_01)"))
    r1 <- fit_model(f1, df, binary) %>%
      mutate(dv = dv, iv = iv, interaction = "Yes", controls = "No", RE = "Yes")

    # Random effects with controls
    f2 <- as.formula(paste("value ~ time +", iv, "+ time:", iv, "+",
                           paste(ctrls, collapse = "+"), "+ (1 | x03_01)"))
    r2 <- fit_model(f2, df, binary) %>%
      mutate(dv = dv, iv = iv, interaction = "Yes", controls = "Yes", RE = "Yes")

    results[[paste(dv, iv, sep = "_")]] <- bind_rows(r1, r2)

  }
}

results_df <- bind_rows(results) %>%
  select(dv, iv, interaction, controls, RE, term, estimate, std.error, statistic, p.value)

# =============================================================================
# TABLES
# =============================================================================

add_stars <- function(p) {
  case_when(p < 0.01 ~ "***", p < 0.05 ~ "**", p < 0.1 ~ "*", TRUE ~ "")
}

fmt_cell <- function(estimate, std_error, p_value) {
  paste0(sprintf("%.3f", estimate), add_stars(p_value),
         " (", sprintf("%.3f", std_error), ")")
}

# "Enviromental" matches the typo in the data
dv_levels <- c("whoqol", "General", "Physical", "Psychological", "Social", "Enviromental")
dv_labels <- c("WHOQOL-BREF", "General", "Physical", "Psychological", "Social", "Environmental")

# N and control mean (baseline, time == 0) per outcome
stats_dv <- data %>%
  group_by(var) %>%
  summarise(
    N              = n(),
    `Control mean` = round(mean(value[time == 0], na.rm = TRUE), 2),
    .groups = "drop"
  ) %>%
  rename(dv = var)

# --- Table 1: OLS (mixed effects, with controls) ----------------------------
# z and d are factors here -> interaction terms are "time:z1" and "time:d1"

ols_wide <- results_df %>%
  filter(
    controls == "No",
    RE       == "Yes",
    term %in% c("time:z1", "time:value_previos", "time:mean_previos")
  ) %>%
  mutate(
    col = case_when(
      term == "time:z1"            ~ "Attendance",
      term == "time:value_previos" ~ "Length of use",
      term == "time:mean_previos"  ~ "Frequency"
    ),
    cell    = fmt_cell(estimate, std.error, p.value),
    Outcome = as.character(factor(dv, levels = dv_levels, labels = dv_labels))
  ) %>%
  select(dv, Outcome, col, cell) %>%
  pivot_wider(names_from = col, values_from = cell) %>%
  left_join(stats_dv, by = "dv") %>%
  arrange(match(Outcome, dv_labels)) %>%
  select(Outcome, N, `Control mean`, Attendance, `Length of use`, Frequency)

ols_wide %>%
  kable(format    = "pipe",
        align     = c("l", "r", "r", "c", "c", "c"),
        col.names = c("Outcome", "N", "Control mean",
                      "Invitation", "Length of use", "Frequency"),
        caption   = "Table 1. Effects on quality of life — OLS (mixed effects, with controls)") %>%
  kable_styling(full_width = FALSE) %>%
  add_header_above(c(" " = 3, "Coefficient (time \u00d7 exposure)" = 3)) %>%
  footnote(general       = "*** p<0.01, ** p<0.05, * p<0.1. Standard errors in parentheses. Coefficient = interaction term time \u00d7 exposure. Mixed effects: random intercept by individual.",
           general_title = "Note: ")

# --- Table 2: IV (instrument z, with controls) ------------------------------

iv_wide <- results_df_iv %>%
  filter(
    controls == "Yes",
    term %in% c("time:d1", "time:value_previos", "time:mean_previos")
  ) %>%
  mutate(
    col = case_when(
      term == "time:d1"            ~ "Attendance",
      term == "time:value_previos" ~ "Length of use",
      term == "time:mean_previos"  ~ "Frequency"
    ),
    cell    = fmt_cell(estimate, std.error, p.value),
    Outcome = as.character(factor(dv, levels = dv_levels, labels = dv_labels))
  ) %>%
  select(dv, Outcome, col, cell) %>%
  pivot_wider(names_from = col, values_from = cell) %>%
  left_join(stats_dv, by = "dv") %>%
  arrange(match(Outcome, dv_labels)) %>%
  select(Outcome, N, `Control mean`, Attendance, `Length of use`, Frequency)

iv_wide %>%
  kable(format    = "pipe",
        align     = c("l", "r", "r", "c", "c", "c"),
        col.names = c("Outcome", "N", "Control mean",
                      "Attendance", "Length of use", "Frequency"),
        caption   = "Table 2. Effects on quality of life — IV (instrument: z, with controls)") %>%
  kable_styling(full_width = FALSE) %>%
  add_header_above(c(" " = 3, "Coefficient (time \u00d7 exposure)" = 3)) %>%
  footnote(general       = "*** p<0.01, ** p<0.05, * p<0.1. Standard errors in parentheses. Instrument: invitation letter (z) for attendance (d).",
           general_title = "Note: ")


data %>%
  filter(var == "whoqol", time == 0) %>%
  ggplot(aes(x = factor(z), y = value)) +

  annotate("rect",
           xmin = -Inf, xmax = Inf,
           ymin = 2.643, ymax = 3,
           fill = "grey70", alpha = 0.25) +

  geom_jitter(width = 0.001, height = 0, alpha = 0.7, size = 1.6) +

  coord_cartesian(ylim = c(1, 5)) +
  scale_y_continuous(breaks = 1:5) +
  scale_x_discrete(expand = expansion(mult = c(0.02, 0.02))) +

  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank()
  ) +
  labs(x = NULL,
       y = "WHOQOL-BREF (baseline)")

#
data %>%
  filter(var == "whoqol", time %in% c(0,2)) %>%
  ggplot(aes(x = 1, y = value,
             group = x03_01)) +
  annotate("rect",
           xmin = -Inf, xmax = Inf,
           ymin = 2.643, ymax = 3,
           fill = "grey70", alpha = 0.25) +

  geom_jitter(width = 0.002, height = 0, alpha = 0.7, size = 1.6) +
  coord_cartesian(ylim = c(1, 5)) +
  scale_y_continuous(breaks = 1:5) +
  scale_x_discrete(expand = expansion(mult = c(0.02, 0.02))) +

  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank()
  ) +
  labs(x = NULL,
       y = "WHOQOL-BREF (baseline)")




  #
p = data %>%
  filter(var == "whoqol", time %in% c(0, 2)) %>%
  mutate(
    z_cat = case_when(
      z == 0 ~ 1,
      z == 1 & value_previos == 0 ~ 2,
      z == 1 & value_previos > 0  ~ 3,
      TRUE ~ NA_real_
    )
  ) %>%
  group_by(x03_01) %>%
  mutate(group = max(z_cat, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    group = factor(
      group,
      levels = c(1, 2, 3),
      labels = c("Waiting\nlist",
                 "Invited\nNot attended",
                 "Attended")
    )
  ) %>%
  select(x03_01, time, value, group) %>%
  pivot_wider(names_from = time, values_from = value, names_prefix = "t_") %>%
  distinct(x03_01, group, t_0, t_2) %>%
  mutate(
    x_jit = 1 + runif(n(), -0.03, 0.03),
    direction = case_when(
      t_2 > t_0 ~ "Up",
      t_2 < t_0 ~ "Down",
      TRUE      ~ "No change"
    ),
    # order for legend
    direction = factor(direction, levels = c("Up", "No change", "Down"))
  ) %>%
  filter(direction != "No change") %>%
  ggplot() +
  geom_segment(
    aes(x = x_jit, xend = x_jit, y = t_0, yend = t_2, color = direction, linetype = direction),
    arrow = arrow(length = unit(0.12, "cm")),
    alpha = 0.7,
    linewidth = 0.8
  ) +
  scale_color_manual(
    values = c(
      "Up" = "#2C7FB8",
      # "No change" = "grey60",
      "Down" = "#D7301F"
    ),
    breaks = c("Up",  "Down")
  ) +
  facet_wrap(~ group, nrow = 1) +
  scale_y_continuous(limits = c(1, 5), breaks = 1:5) +
  coord_cartesian(xlim = c(0.97, 1.03)) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    # legend on top
    legend.position = "top",
    legend.direction = "horizontal",
    legend.justification = "center",

    # make facet titles bigger
    strip.text = element_text(size = 12, face = "bold"),
    strip.background = element_rect(fill = "grey95", color = NA),

    # clearly separate facets but keep minimal
    panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.7),
    panel.spacing = unit(1.0, "lines"),

    panel.grid.minor = element_blank()
  ) +
  labs(x = NULL, y = "WHOQOL-BREF") +
  guides(
    color = guide_legend(
      override.aes = list(
        linewidth = 1,   # grosor de la flecha en la leyenda
        alpha = 1
      )
    )
  )


ggsave(
  "3_viz/paper/ups.jpg",
  plot = p,
  width = 7,
  height = 5,
  dpi = 1000
)

ggsave(
  "3_viz/paper/figure2.tif",
  plot = p,
  width = 7,
  height = 5,
  dpi = 100,
  device = "tiff"
)
