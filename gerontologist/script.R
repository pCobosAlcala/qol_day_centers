# Libraries and setup ----
pacman::p_load(tidyverse, AER, broom, broom.mixed, lme4, lmerTest,
               kableExtra, parallel)
options(scipen = 999)
rm(list = ls())

# Viz----
data <- read_csv("C://Users/pablo/OneDrive/Documentos/pablo/RStudio/2023/imss/imss/1_data/bases/base_paper/gerontologist/db.csv")

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

# IV table
table_iv = results_df_iv %>%
  filter(term %in% c("time:mean_previos", "time:value_months", "time:value_previos", "time:d1")) %>%
  mutate(sign = case_when(
    p.value < 0.01 ~ "***",
    p.value < 0.05 ~ "**",
    p.value < 0.1 ~ "*",
    TRUE ~ ""
  )) %>%
  mutate(term_clean = case_when(
    str_detect(term, "mean_previos") ~ "mean_previos",
    str_detect(term, "value_months") ~ "value_months",
    str_detect(term, "value_previos") ~ "value_previos",
    str_detect(term, "d1") ~ "d",
    TRUE ~ term
  )) %>% 
  filter(controls == "No") %>% 
  mutate(est_se = paste0(ifelse(model == "IV", "<i>", ""),
                         sprintf("%.3f", estimate), sign, " (",
                         sprintf("%.3f", std.error), ")",
                         ifelse(model == "IV", "</i>", ""))) %>%
  select(dv, controls, term_clean, est_se) %>% 
  unite(temp, controls, term_clean) %>% 
  mutate(temp = str_replace_all(temp, c("No_" = ""))) %>% 
  mutate(id = paste(dv, temp, sep = "_")) %>% 
  select(id, est_se)


# OLS
deps <- c("Enviromental", "Psychological", "Social", "Physical", "General", "whoqol")
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

# Table
results_df %>% 
  filter(controls == "No") %>% 
  filter(term %in% c("time:mean_previos", "time:value_months", "time:value_previos", "time:z1")) %>% 
  mutate(sign = case_when(
    p.value < 0.01 ~ "***",
    p.value < 0.05 ~ "**",
    p.value < 0.1 ~ "*",
    TRUE ~ ""
  )) %>% 
  mutate(term_clean = case_when(
    str_detect(term, "mean_previos") ~ "mean_previos",
    str_detect(term, "value_months") ~ "value_months",
    str_detect(term, "value_previos") ~ "value_previos",
    str_detect(term, "z1") ~ "d",
    TRUE ~ term
  )) %>% 
  mutate(est_se = paste0(sprintf("%.3f", estimate), sign, " (", sprintf("%.3f", std.error), ")")) %>%
  select(dv, term_clean, est_se) %>% 
  mutate(id = paste(dv, term_clean, sep = "_")) %>% 
  left_join(table_iv, by = "id") %>% 
  select(-id) %>% 
  mutate(est_se = paste(est_se.x, est_se.y, sep = "<br>")) %>% 
  select(-est_se.x, -est_se.y) %>% 
  
  pivot_wider(names_from = term_clean, values_from = est_se) %>% 
  mutate(dv_label = dv) %>% 
  mutate(dv = str_replace_all(dv, c("whoqol" = "WHOQOL-BREF"))) %>% 
  mutate(dv_label = str_replace_all(dv_label, c("whoqol" = "WHOQOL-BREF"))) %>% 
  mutate(dv_label = factor(dv_label, levels = c("WHOQOL-BREF", "General", "Physical", "Psychological", "Social",   "Enviromental"))) %>% 
  arrange(dv_label) %>% 
  select( -dv) %>% 
  select(dv_label, everything()) %>% 
  rename(`Dependent Variable` = dv_label) %>% 
  select(-value_months) %>% 
  kable(format = "html", escape = FALSE,
        col.names = c("DV",
                      "Invitation/attendance",
                      "Length of use", "Frequency"),
        align = "lcccc") %>% 
  kable_styling(full_width = FALSE, bootstrap_options = c("striped", "hover"))
