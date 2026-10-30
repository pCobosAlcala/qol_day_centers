# 00_data.R  (revision, October 2026)
#
# THE ANALYSIS DATA, BUILT ONCE AND IN ONE PLACE, FROM THE RAW FILES.
#
# Every analysis script sources this file right after 00_paths.R and then calls the builder:
#
#   source("3_code/scripts/00_paths.R")
#   source("3_code/scripts/00_data.R")
#   dat <- build_analysis_data()                 # one row per randomised person (120)
#   s   <- ancova_sample(dat, "whoqol", 2)       # the endline ANCOVA sample of one outcome
#
# Sourcing only defines objects (labels, item maps, functions). Nothing is read or written until
# build_analysis_data() is called, and that function never writes a file.
#
# ---------------------------------------------------------------------------------------------
# WHAT build_analysis_data() RETURNS: a tibble with one row per person interviewed at baseline,
# which is everyone randomised (120: 91 ever invited, 29 never invited).
#
#   id                     study folio (PM001 ...). Person level: never write rows out.
#   z                      1 = invited in any of the eight randomisation rounds (t1.csv), 0 = never
#   <outcome>_y0, _y1, _y2 WHOQOL-BREF at baseline, midline (about 3 months), endline (about 9),
#                          on the WHO 0 to 100 scale = 25 x (item mean minus 1). NA if not surveyed.
#                          outcomes: whoqol (composite of the standard 26 items), general (items 1
#                          and 2), physical, psychological, social, environmental (the four domains)
#   d1, d2                 1 = attended the centre on more than five days (more than one week) up to
#                          and including that wave's survey date; NA if not surveyed at that wave
#   days1, days2           days attended up to and including that wave's survey date
#   weeks1, weeks2         days / 5: weeks of a five day week (the old value_previos1, value_previos2)
#   dpw1, dpw2             days attended per week since the first visit: 7 x days / calendar days
#                          from the first visit to the survey (5 x the old mean_previos); 0 if days = 0
#   age0                   age in years at the baseline interview, from the birth date
#   female, indigenous     1 = woman; 1 = considers herself or himself indigenous
#   hh_size                people living with the older adult
#   educ                   last completed level, 1 none, 2 preschool, 3 primary, 4 lower secondary,
#                          5 technical after lower secondary, 6 upper secondary, 7 teacher training,
#                          8 technical after upper secondary, 9 university, 10 postgraduate
#   child_carer            1 = main caregiver is a son or daughter ("Padre o madre" or "Hija/o",
#                          the corrected coding of 04 and 06), NA if no caregiver record
#   katz0, lawton0         Katz (0 to 6 basic activities) and Lawton (0 to 8 instrumental
#                          activities) at baseline, number of activities done independently
#
#   With dates = TRUE it also returns date0, date1, date2 (survey dates) and entry_date (first day
#   attended). INTERNAL ONLY: survey and entry dates of identifiable people are never written out.
#
# ---------------------------------------------------------------------------------------------
# HOW EACH PIECE IS BUILT, AND HOW IT IS CHECKED (verify = TRUE, the default, stops on any failure)
#
# (1) WHOQOL-BREF from the raw items (pm_limpieza1.csv, long by item). The Mexican questionnaire
#     numbers its 28 questions 16.01 to 16.28 plus an interviewer check 16.99. The map to the WHO
#     items 1 to 26 (WHOQOL_ITEMS below) is the one the old pipeline used (pm_qol_desagregado.R,
#     which wrote pm_who.csv), checked against the question texts in cuestionarios_CD_imss.xlsx
#     (1_data/bases/base_final, subfolder "Envios" with an accented i). Not scored: 16.01 ("are you
#     currently ill?", unscored in the WHO form too), 16.03 ("do you get the kind of support you
#     need from others?", a WHOQOL-100 question that is not in the BREF, but was in the submitted
#     27 item total), 16.99 (interviewer check). Two caveats on the questionnaire itself:
#       - WHO item 2 ("how satisfied are you with your health?") is 16.17. Its header in the
#         export opens with the instruction for the satisfaction block, but in all three PM
#         sheets of cuestionarios_CD_imss.xlsx it ends with the satisfaction with health
#         question itself, so 16.17 is the right item (checked by the referee, 2 October 2026).
#       - WHO item 19 ("how satisfied are you with yourself?") is worded in 16.21 as "how
#         satisfied do you feel with your life?".
#     Scoring follows the WHO manual: items 3, 4 and 26 are reversed (6 minus the answer); a
#     missing item takes the person's own mean of the other items of its domain; a domain with
#     more than two missing items (Social: more than one) is not scored; an assessment with more
#     than 20 per cent of the 26 items missing is discarded. The WHO gives no rule for the two
#     general items; we score "general" only when both are present. The composite is the mean
#     of the 26 items (missing items imputed as above), which equals the item weighted mean of
#     the five parts. In these data no assessment that took place has a missing item, so the
#     missing item rules never bind.
#     CHECK: on the 1 to 5 scale the five parts and the old 27 item composite reproduce
#     pm_who.csv for its 117 people at every wave (tolerance 1e-9).
#
# (2) Everyone randomised. Three invited people with baseline and midline questionnaires never
#     reached pm.csv or pm_who.csv: the old 3_asistencias.R set their attendance to missing
#     because they entered the centre after their midline survey. They are included here with
#     z = 1, d1 = 0, and no endline survey. Midline sample 99 + 3 = 102; endline stays 85.
#
# (3) Attendance by wave from the centre's register (CD CASSAAM 16febrero24 asistencias.xlsx, one
#     column per opening day, header = Excel date serial, "x" or "X" = attended) and the survey
#     dates (rows "marca_" of pm_limpieza1.csv). Days are counted up to and including the survey
#     date, as the old value_previos did. Other marks are not attendance: a word spelled one
#     letter per row down the column of 18 September 2023, and "Teatro" on 15 February 2024.
#     CHECK: weeks and days per week reproduce value_previos and mean_previos of pm_who.csv for
#     the 117, z reproduces pm_who.csv, d2 equals the old time invariant d in the endline sample,
#     and nobody never invited has any attendance mark (one sided non compliance).
#
# (4) Baseline covariates recoded from the raw baseline answers with the project's own coding,
#     as in 04_balance_tables.R. Katz and Lawton follow the item coding of the old pipeline.
#     Note: that Katz coding is more lenient than the standard index (toileting or transfer with
#     assistance and occasional incontinence count as independent). Katz enters only balance
#     tables and covariate checks, never the main estimates.
#     CHECK: sex, indigenous, household size, education and the OLD caregiver flag reproduce
#     control.csv, and Katz and Lawton reproduce pm.csv, for the 117 people those files hold.
#
# Answer labels are matched after removing accents (stringi Latin-ASCII), so this file is plain
# ASCII and does not depend on the encoding of the script. Any label not listed here stops the
# build. Written in the tidyverse, like the other scripts of this folder.

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(tibble)
})

# =============================================================================================
# 1. LABELS SHARED BY EVERY SCRIPT
# =============================================================================================
# column stem in the data (names) and the label printed in tables and figures (values),
# in a neutral order: composite first, then the parts in WHOQOL order
OUTCOMES <- c(whoqol = "WHOQOL-BREF", general = "General", physical = "Physical",
              psychological = "Psychological", social = "Social", environmental = "Environmental")

# =============================================================================================
# 2. WHOQOL-BREF: ITEM MAP (Mexican questionnaire, section 16, to WHO items 1 to 26)
# =============================================================================================
WHOQOL_ITEMS <- tribble(
  ~who, ~src,     ~part,           ~reverse, ~content,
  1,    "x16_02", "general",       FALSE,    "overall quality of life",
  2,    "x16_17", "general",       FALSE,    "satisfaction with health (16.17, see header)",
  3,    "x16_04", "physical",      TRUE,     "pain prevents doing what is needed",
  4,    "x16_05", "physical",      TRUE,     "needs medical treatment to function",
  5,    "x16_06", "psychological", FALSE,    "enjoys life",
  6,    "x16_07", "psychological", FALSE,    "life is meaningful",
  7,    "x16_08", "psychological", FALSE,    "able to concentrate",
  8,    "x16_09", "environmental", FALSE,    "feels safe in daily life",
  9,    "x16_10", "environmental", FALSE,    "healthy physical environment",
  10,   "x16_11", "physical",      FALSE,    "enough energy",
  11,   "x16_12", "psychological", FALSE,    "accepts bodily appearance",
  12,   "x16_13", "environmental", FALSE,    "enough money",
  13,   "x16_14", "environmental", FALSE,    "information available",
  14,   "x16_15", "environmental", FALSE,    "opportunity for leisure",
  15,   "x16_16", "physical",      FALSE,    "able to get around",
  16,   "x16_18", "physical",      FALSE,    "satisfied with sleep",
  17,   "x16_19", "physical",      FALSE,    "satisfied with daily activities",
  18,   "x16_20", "physical",      FALSE,    "satisfied with capacity for work",
  19,   "x16_21", "psychological", FALSE,    "satisfied with self (asked as 'with your life')",
  20,   "x16_22", "social",        FALSE,    "satisfied with personal relationships",
  21,   "x16_23", "social",        FALSE,    "satisfied with sex life",
  22,   "x16_24", "social",        FALSE,    "satisfied with support from friends",
  23,   "x16_25", "environmental", FALSE,    "satisfied with living place",
  24,   "x16_26", "environmental", FALSE,    "satisfied with access to health services",
  25,   "x16_27", "environmental", FALSE,    "satisfied with transport",
  26,   "x16_28", "psychological", TRUE,     "frequency of negative feelings")

# WHO domain sizes, as a check on the map: general 2, physical 7, psychological 6, social 3,
# environmental 8 (26 items)
stopifnot(nrow(WHOQOL_ITEMS) == 26, !anyDuplicated(WHOQOL_ITEMS$src),
          identical(as.integer(table(WHOQOL_ITEMS$part)[c("general", "physical", "psychological",
                                                          "social", "environmental")]),
                    c(2L, 7L, 6L, 3L, 8L)),
          identical(sort(WHOQOL_ITEMS$who[WHOQOL_ITEMS$reverse]), c(3, 4, 26)))

# answer label (accents removed) to the 1 to 5 answer, in the item's own direction (before the
# reversal of items 3, 4 and 26). Same values as the old pipeline.
WHOQOL_RESPONSES <- c(
  # digits: most forms record the satisfaction block (items 2, 16 to 25) as 1 to 5,
  # 1 = very dissatisfied, 5 = very satisfied
  "1" = 1, "2" = 2, "3" = 3, "4" = 4, "5" = 5,
  # item 1, rating of quality of life: very poor ... very good
  "Muy mala" = 1, "Muy pobre" = 1, "Mala" = 2, "Pobre" = 2,
  "Ni mala ni buena" = 3, "Ni pobre ni buena" = 3, "Buena" = 4, "Muy buena" = 5,
  # items 3 to 15, "how much", "to what extent", "how well": not at all ... extremely, completely
  "Para nada" = 1, "Nada en absoluto" = 1, "No en absoluto" = 1, "Ninguna en absoluto" = 1,
  "Un poco" = 2, "Pocas" = 2,
  "Moderadamente" = 3, "En cantidad moderada" = 3, "Mas o menos" = 3,
  "Mucho" = 4, "Bastante" = 4, "Mayormente" = 4, "La mayor parte del tiempo" = 4,
  "En cantidad extrema" = 5, "Extremadamente" = 5, "Completamente" = 5, "Totalmente" = 5,
  # satisfaction block written as words (13 baseline forms): very dissatisfied ... very satisfied;
  # the second spelling of the middle answer is a typo of the form itself
  "Muy insatisfecha/o" = 1, "Insatisfecha/o" = 2,
  "Ni satisfecha/o ni insatisfecha/o" = 3, "Ni insatisfecha/o ni insatisfecha/o" = 3,
  "Satisfecha/o" = 4, "Muy satisfecha/o" = 5,
  # item 26, frequency of negative feelings: never ... always
  "Nunca" = 1, "Rara vez" = 2, "Con cierta frecuencia" = 3, "Con frecuencia" = 3,
  "Muy frecuentemente" = 4, "Muy seguido" = 4, "Siempre" = 5)

# =============================================================================================
# 3. KATZ AND LAWTON: answer label (accents removed) to 1 = independent, 0 = dependent.
#    Coding of the old pipeline (pm.R and pm_qol_desagregado.R, identical in this part).
#    Lawton has 8 activities, x62_01 and x62_04 to x62_09 (there is no x62_03).
# =============================================================================================
FUNCTION_ITEMS <- tribble(
  ~item,    ~label,                                                                                ~value,
  # Katz, basic activities of daily living
  "x61_01", "No recibe asistencia (puede entrar y salir de la tina u otra forma de bano)",         1,
  "x61_01", "Recibe asistencia durante el bano en una sola parte del cuerpo (ej. espalda o pierna)", 1,
  "x61_01", "Recibe asistencia durante el bano en mas de una parte del cuerpo",                    0,
  "x61_02", "Puede tomar las prendas y vestirse completamente, sin asistencia",                    1,
  "x61_02", "Puede tomar las prendas y vestirse sin asistencia excepto en abrocharse los zapatos", 1,
  "x61_02", "Recibe asistencia para tomar las prendas y vestirse",                                 0,
  "x61_03", paste("Usa el sanitario ninguna asistencia (puede utilizar algun objeto de soporte como",
                  "baston o silla de ruedas, o puede arreglar su ropa o el panal que usa, en caso de hacerlo)"), 1,
  "x61_03", paste("Recibe asistencia al ir al bano en limpiarse (puede manejar por si misma/o el",
                  "panal o comodo vaciandolo en caso de usarlo)"),                                 1,
  "x61_03", "No va al bano por si misma/o",                                                        0,
  "x61_04", paste("Se mueve dentro y fuera de la cama y silla sin ninguna asistencia (puede estar",
                  "utilizando un auxiliar de la marcha u objeto de soporte)"),                     1,
  "x61_04", "Puede moverse dentro y fuera de la cama y silla con asistencia",                      1,
  "x61_04", "No puede salir de la cama",                                                           0,
  "x61_05", "Tiene control total de esfinteres",                                                   1,
  "x61_05", "Tiene accidentes ocasionales que no afectan su vida social",                          1,
  "x61_05", "Necesita ayuda para supervision del control de esfinteres, utiliza sonda o es incontinente", 0,
  "x61_06", "Se alimenta por si sola/o sin asistencia alguna",                                     1,
  "x61_06", "Se alimenta sola/o y tiene asistencia solo para cortar la carne o untar mantequilla", 1,
  "x61_06", "Recibe asistencia en la alimentacion, o se alimenta parcial o totalmente por via enteral o parenteral", 0,
  # Lawton, instrumental activities of daily living
  "x62_01", "Utiliza el telefono por iniciativa propia",                                           1,
  "x62_01", "Es capaz de marcar bien algunos numeros familiares",                                  1,
  "x62_01", "Es capaz de contestar al telefono, pero no de marcar",                                1,
  "x62_01", "No es capaz de usar el telefono",                                                     0,
  "x62_02", "Realiza todas las compras necesarias independientemente",                             1,
  "x62_02", "Realiza independientemente pequenas compras",                                         0,
  "x62_02", "Necesita ir acompanada/o para hacer cualquier compra",                                0,
  "x62_02", "Totalmente incapaz de comprar",                                                       0,
  "x62_04", "Organiza, prepara y sirve las comidas por si sola/o adecuadamente",                   1,
  "x62_04", "Prepara adecuadamente las comidas si se le proporcionan los ingredientes",            0,
  "x62_04", "Prepara, calienta y sirve las comidas, pero no sigue una dieta adecuada",             0,
  "x62_04", "Necesita que le preparen y sirvan las comidas",                                       0,
  "x62_05", "Mantiene la casa sola/o o con ayuda ocasional (para trabajos pesados)",               1,
  "x62_05", "Realiza tareas ligeras, como lavar los platos o hacer las camas",                     1,
  "x62_05", "Realiza tareas ligeras, pero no puede mantener un adecuado nivel de limpieza",        1,
  "x62_05", "Necesita ayuda en todas las labores de la casa",                                      1,
  "x62_05", "No participa en ninguna labor de la casa",                                            0,
  "x62_06", "Lava por si sola/o toda su ropa",                                                     1,
  "x62_06", "Lava por si sola/o pequenas prendas",                                                 1,
  "x62_06", "Todo el lavado de ropa debe ser realizado por otra persona",                          0,
  "x62_07", "Viaja sola/o en transporte publico o conduce su propio coche",                        1,
  "x62_07", "Es capaz de coger un taxi, pero no usa otro medio de transporte",                     1,
  "x62_07", "Viaja en transporte publico cuando va acompanada/o por otra persona",                 1,
  "x62_07", "Solo utiliza el taxi o el automovil con ayuda de otras personas",                     0,
  "x62_07", "No viaja",                                                                            0,
  "x62_08", "Es capaz de tomar su medicacion a la hora y con la dosis correcta",                   1,
  "x62_08", "Toma su medicacion si la dosis le es preparada previamente",                          0,
  "x62_08", "No es capaz de administrarse su medicacion",                                          0,
  "x62_09", "Se encarga de sus asuntos economicos por si sola/o",                                  1,
  "x62_09", "Realiza las compras de cada dia, pero necesita ayuda en las grandes compras, bancos", 1,
  "x62_09", "No es capaz de manejar dinero",                                                       0)
stopifnot(!anyDuplicated(FUNCTION_ITEMS[, c("item", "label")]),
          length(unique(FUNCTION_ITEMS$item)) == 14)               # 6 Katz + 8 Lawton activities

# education, last completed level (accents removed), coded as in 04_balance_tables.R
EDUC_CODES <- c("Ninguno" = 1, "Preescolar o kinder" = 2, "Primaria" = 3, "Secundaria" = 4,
                "Estudios tecnicos con secundaria terminada" = 5, "Preparatoria o bachillerato" = 6,
                "Normal basica" = 7, "Estudios tecnicos con preparatoria terminada" = 8,
                "Licenciatura o ingenieria (profesional)" = 9, "Maestria o doctorado" = 10)

# =============================================================================================
# 4. SMALL HELPERS
# =============================================================================================
.ascii <- function(x) stringi::stri_trans_general(x, "Latin-ASCII")   # drops accents: "Si", "bano"

# readers: everything as text first, so no column is guessed into a date or a number
.read_raw_csv <- function(raw, f) {
  read_csv(paste0(raw, f), col_types = cols(.default = col_character()),
           name_repair = "unique_quiet", progress = FALSE)
}
.read_raw_num <- function(raw, f) {           # processed files of the old pipeline (checks only)
  suppressMessages(read_csv(paste0(raw, f), show_col_types = FALSE, name_repair = "unique_quiet",
                            progress = FALSE))
}

# =============================================================================================
# 5. THE BUILDER
# =============================================================================================
build_analysis_data <- function(dates = FALSE, verify = TRUE, quiet = FALSE, raw = RAW) {

  # ---- 5.0 public copy: without the raw files, read the de-identified analysis file written by
  #      11_public_data.R (same columns, random ids, age rounded to whole years, no dates)
  if (!dir.exists(raw)) {
    if (dates) stop("This script needs the survey dates, which are not public; it runs only with the raw data.")
    if (!quiet) message("Raw data not found: using the de-identified file ", PUBLIC_DATA)
    return(read_csv(PUBLIC_DATA, show_col_types = FALSE))
  }

  # ---- 5.1 raw questionnaire file, long by item: x03_01 (folio), name (item), t0, t1, t2 ------
  lim <- .read_raw_csv(raw, "pm_limpieza1.csv") %>%
    select(id = x03_01, name, t0, t1, t2)

  # survey dates: the form's timestamp ("marca_"), date part only. A wave without a timestamp
  # was not surveyed.
  sv <- lim %>% filter(name == "marca_") %>%
    transmute(id, date0 = as.Date(substr(t0, 1, 10)), date1 = as.Date(substr(t1, 1, 10)),
              date2 = as.Date(substr(t2, 1, 10)))
  stopifnot(nrow(sv) == 120, !anyDuplicated(sv$id), all(!is.na(sv$date0)))

  # randomised invitation: ever selected in the eight rounds (built in the old 2_bases_tratamiento.R)
  sel <- .read_raw_csv(raw, "t1.csv") %>% transmute(id = x03_01, z = as.integer(t1))
  stopifnot(setequal(sel$id, sv$id), !anyDuplicated(sel$id), all(sel$z %in% 0:1))

  # ---- 5.2 WHOQOL-BREF ------------------------------------------------------------------------
  # one row per person, wave and WHO item, for every wave the person was surveyed
  held <- sv %>% pivot_longer(-id, names_to = "wave", values_to = "date") %>%
    filter(!is.na(date)) %>% transmute(id, wave = as.integer(sub("date", "", wave)))
  ans <- lim %>% filter(name %in% WHOQOL_ITEMS$src) %>%
    pivot_longer(c(t0, t1, t2), names_to = "wave", values_to = "label") %>%
    transmute(id, src = name, wave = as.integer(sub("t", "", wave)), label = .ascii(label))
  unknown <- setdiff(unique(na.omit(ans$label)), names(WHOQOL_RESPONSES))
  if (length(unknown)) stop("unmapped WHOQOL answer label(s): ", paste(unknown, collapse = " | "))
  items <- held %>%
    tidyr::crossing(src = WHOQOL_ITEMS$src) %>%           # all 26 items for every survey held
    left_join(ans, by = c("id", "wave", "src")) %>%
    left_join(WHOQOL_ITEMS %>% select(src, who, part, reverse), by = "src") %>%
    mutate(answer = unname(WHOQOL_RESPONSES[label]),       # 1 to 5 in the item's own direction
           score  = if_else(reverse, 6 - answer, answer))  # WHO reversal of items 3, 4 and 26

  # WHO missing item rules (see header). With no missing item a part is just its item mean.
  parts <- items %>%
    group_by(id, wave) %>% mutate(miss_all = sum(is.na(score))) %>%
    group_by(id, wave, part) %>%
    summarise(k = n(), miss = sum(is.na(score)), miss_all = miss_all[1],
              mean_obs = if (all(is.na(score))) NA_real_ else mean(score, na.rm = TRUE),
              .groups = "drop") %>%
    mutate(max_miss = case_when(part == "social" ~ 1L, part == "general" ~ 0L, TRUE ~ 2L),
           value = if_else(miss <= max_miss & miss_all <= 5, mean_obs, NA_real_))  # 5 of 26 = 19%
  # composite = mean of the 26 items, each missing item imputed by its own domain mean, which is
  # the item weighted mean of the five parts; not scored if any part is not
  comp <- parts %>% group_by(id, wave) %>%
    summarise(value = if (any(is.na(value))) NA_real_ else sum(value * k) / sum(k),
              .groups = "drop") %>%
    mutate(part = "whoqol")
  who15 <- bind_rows(parts %>% select(id, wave, part, value), comp)   # 1 to 5 scale, long
  n_missing_items <- sum(is.na(items$score))                          # reported by verify

  # WHO 0 to 100 scale, wide: <outcome>_y<wave>
  who100 <- who15 %>%
    mutate(value = 25 * (value - 1), col = paste0(part, "_y", wave)) %>%
    select(id, col, value) %>%
    pivot_wider(names_from = col, values_from = value)
  ycols <- paste0(rep(names(OUTCOMES), each = 3), "_y", 0:2)
  for (cc in setdiff(ycols, names(who100))) who100[[cc]] <- NA_real_ # (a wave nobody answered)
  who100 <- who100 %>% select(id, all_of(ycols))

  # ---- 5.3 attendance register ------------------------------------------------------------------
  reg <- readxl::read_excel(paste0(raw, "CD CASSAAM 16febrero24 asistencias.xlsx"),
                            .name_repair = "minimal")
  serial <- suppressWarnings(as.numeric(names(reg)))
  is_day <- !is.na(serial) & serial > 40000 & serial < 50000           # Excel date serials
  stopifnot(sum(is_day) == 223, sum(names(reg) == "Folio") == 1)       # 223 opening days
  day_of <- as.Date(serial[is_day], origin = "1899-12-30")             # 18 Jan 2023 to 16 Feb 2024
  marks <- reg[, c(which(names(reg) == "Folio"), which(is_day))]
  names(marks) <- c("id", as.character(day_of))
  visits <- marks %>%
    pivot_longer(-id, names_to = "day", values_to = "mark",
                 values_transform = list(mark = as.character)) %>%
    filter(tolower(mark) %in% "x") %>%                                 # a day attended
    transmute(id, day = as.Date(day))
  # the register holds one folio (PM066) that is not in the survey; it has no attendance marks
  stopifnot(all(unique(visits$id) %in% sv$id))

  att <- sv %>%
    left_join(visits, by = "id", relationship = "one-to-many") %>%
    group_by(id, date1, date2) %>%
    summarise(entry_date = if (all(is.na(day))) as.Date(NA) else min(day, na.rm = TRUE),
              days_all   = sum(!is.na(day)),
              days1      = sum(day <= date1[1], na.rm = TRUE),         # up to and including
              days2      = sum(day <= date2[1], na.rm = TRUE),         # the survey date
              .groups = "drop") %>%
    mutate(
      # a wave without a survey has no "before the survey": attendance measures are NA there
      days1  = if_else(is.na(date1), NA_integer_, days1),
      days2  = if_else(is.na(date2), NA_integer_, days2),
      d1     = as.integer(days1 > 5),                                  # more than one week
      d2     = as.integer(days2 > 5),
      weeks1 = days1 / 5,                                              # old value_previos
      weeks2 = days2 / 5,
      # old mean_previos x 5: days attended per 7 calendar days since the first visit
      dpw1   = if_else(days1 == 0, 0, 7 * days1 / as.numeric(date1 - entry_date)),
      dpw2   = if_else(days2 == 0, 0, 7 * days2 / as.numeric(date2 - entry_date))) %>%
    select(id, entry_date, days_all, d1, d2, days1, days2, weeks1, weeks2, dpw1, dpw2)
  stopifnot(all(is.finite(att$dpw1[!is.na(att$dpw1)])), all(is.finite(att$dpw2[!is.na(att$dpw2)])))

  # ---- 5.4 baseline covariates ------------------------------------------------------------------
  base <- lim %>%
    filter(name %in% c("x04_01", "x05_01", "x05_02", "x05_05", "x05_07", "x06_03_pc")) %>%
    select(id, name, t0) %>%
    pivot_wider(names_from = name, values_from = t0)
  cov <- base %>%
    left_join(sv %>% select(id, date0), by = "id") %>%
    transmute(
      id,
      # birth dates come in three layouts; age at the baseline interview, as in 04
      age0 = as.numeric(date0 - as.Date(suppressWarnings(lubridate::parse_date_time(
        x04_01, orders = c("Ymd", "mdY", "dmY"), quiet = TRUE)))) / 365.25,
      female      = as.integer(.ascii(x05_01) == "Mujer"),
      indigenous  = as.integer(.ascii(x05_02) == "Si"),
      hh_size     = as.numeric(x05_05),
      educ        = unname(EDUC_CODES[.ascii(x05_07)]),
      child_carer = case_when(is.na(x06_03_pc) ~ NA_integer_,
                              x06_03_pc %in% c("Padre o madre", "Hija/o") ~ 1L, TRUE ~ 0L),
      child_old   = as.integer(x06_03_pc %in% "Padre o madre"))      # control.csv coding, check only
  if (any(is.na(cov$educ) & !is.na(base$x05_07))) stop("unmapped education label")

  # Katz and Lawton at baseline: number of activities done independently among those answered
  fn <- lim %>% filter(grepl("^x6[12]_", name), !is.na(t0)) %>%
    transmute(id, item = name, label = .ascii(t0)) %>%
    left_join(FUNCTION_ITEMS, by = c("item", "label"))
  if (any(is.na(fn$value))) stop("unmapped Katz or Lawton label")
  fn <- fn %>% group_by(id) %>%
    summarise(katz0   = if (any(grepl("^x61", item))) sum(value[grepl("^x61", item)]) else NA_real_,
              lawton0 = if (any(grepl("^x62", item))) sum(value[grepl("^x62", item)]) else NA_real_,
              .groups = "drop")

  # ---- 5.5 one row per person -------------------------------------------------------------------
  dat <- sv %>%
    left_join(sel,    by = "id") %>%
    left_join(who100, by = "id") %>%
    left_join(att,    by = "id") %>%
    left_join(cov,    by = "id") %>%
    left_join(fn,     by = "id") %>%
    arrange(id)

  # design facts every script relies on
  stopifnot(nrow(dat) == 120, sum(dat$z) == 91,
            all(dat$days_all[dat$z == 0] == 0),                      # nobody uninvited attended
            sum(!is.na(dat$whoqol_y0)) == 120,
            sum(!is.na(dat$whoqol_y0) & !is.na(dat$whoqol_y1)) == 102,
            sum(!is.na(dat$whoqol_y0) & !is.na(dat$whoqol_y2)) == 85,
            identical(!is.na(dat$whoqol_y1), !is.na(dat$date1)),       # scored = surveyed
            identical(!is.na(dat$whoqol_y2), !is.na(dat$date2)),
            all(complete.cases(dat[, c("age0", "female", "indigenous", "hh_size", "educ")])))

  # ---- 5.6 checks against the files behind the submitted paper ----------------------------------
  if (verify) {
    old_who <- .read_raw_num(raw, "pm_who.csv") %>% select(-1)
    ids117  <- unique(old_who$x03_01)
    stopifnot(length(ids117) == 117,
              setequal(setdiff(dat$id, ids117), c("PM018", "PM049", "PM086")))

    # (a) WHOQOL on the 1 to 5 scale: five parts, and the old 27 item total (26 items + 16.03)
    q03 <- lim %>% filter(name == "x16_03") %>%
      pivot_longer(c(t0, t1, t2), names_to = "wave", values_to = "label") %>%
      transmute(id, wave = as.integer(sub("t", "", wave)), q03 = unname(WHOQOL_RESPONSES[.ascii(label)]))
    tot27 <- items %>% group_by(id, wave) %>% summarise(s26 = sum(score), .groups = "drop") %>%
      left_join(q03, by = c("id", "wave")) %>%
      transmute(id, wave, part = "whoqol27", value = (s26 + q03) / 27)
    old_names <- c(General = "general", Physical = "physical", Psychological = "psychological",
                   Social = "social", Enviromental = "environmental", whoqol = "whoqol27")
    cmp_who <- old_who %>% filter(var %in% names(old_names)) %>%
      transmute(id = x03_01, wave = as.integer(time), part = unname(old_names[var]), old = value) %>%
      full_join(bind_rows(who15 %>% filter(part != "whoqol"), tot27) %>% filter(id %in% ids117),
                by = c("id", "wave", "part"))
    stopifnot(!anyNA(cmp_who$old), !anyNA(cmp_who$value))           # same person-waves scored
    chk_who <- cmp_who %>% group_by(part, wave) %>%
      summarise(n = n(), max_abs_diff = max(abs(old - value)), .groups = "drop")
    if (max(chk_who$max_abs_diff) > 1e-9) stop("rebuilt WHOQOL does not reproduce pm_who.csv")

    # (b) attendance: weeks and days per week against value_previos and mean_previos, and z
    old_att <- old_who %>% distinct(id = x03_01, wave = as.integer(time), z_old = z, d_old = d,
                                    value_previos, mean_previos) %>% filter(wave > 0)
    stopifnot(!anyDuplicated(old_att[, c("id", "wave")]))
    cmp_att <- old_att %>%
      left_join(dat %>% select(id, z, weeks1, weeks2, dpw1, dpw2), by = "id") %>%
      mutate(weeks = if_else(wave == 1, weeks1, weeks2), dpw = if_else(wave == 1, dpw1, dpw2))
    chk_att <- cmp_att %>% filter(!is.na(weeks)) %>% group_by(wave) %>%
      summarise(n = n(), weeks_max_abs_diff = max(abs(weeks - value_previos)),
                share_max_abs_diff = max(abs(dpw / 5 - mean_previos)), .groups = "drop")
    if (max(chk_att$weeks_max_abs_diff, chk_att$share_max_abs_diff) > 1e-9 ||
        any(cmp_att$value_previos[is.na(cmp_att$weeks)] != 0) ||      # old file: 0 if not surveyed
        any(cmp_att$z != cmp_att$z_old)) stop("rebuilt attendance or z does not reproduce pm_who.csv")

    # (c) wave attendance d1, d2 against the old time invariant d, in each ANCOVA sample of the
    #     composite (old d is missing for the three late entrants, who are not in pm_who.csv)
    d_old <- old_att %>% distinct(id, d_old)
    chk_d <- bind_rows(
      dat %>% filter(!is.na(whoqol_y0), !is.na(whoqol_y1)) %>% transmute(id, wave = 1L, d_new = d1, days = days1),
      dat %>% filter(!is.na(whoqol_y0), !is.na(whoqol_y2)) %>% transmute(id, wave = 2L, d_new = d2, days = days2)) %>%
      left_join(d_old, by = "id") %>%
      group_by(wave, d_old, d_new) %>%
      summarise(n = n(), min_days = min(days), max_days = max(days), .groups = "drop")
    if (any(with(chk_d, wave == 2 & (is.na(d_old) | d_old != d_new))))
      stop("d2 differs from the old d in the endline sample")

    # (d) covariates against control.csv and pm.csv for the 117 people they hold
    ctr <- .read_raw_num(raw, "control.csv") %>% select(-1) %>% rename(id = x03_01) %>%
      inner_join(dat, by = "id")
    stopifnot(nrow(ctr) == 117)
    if (!all(ctr$female == 1 - ctr$x05_01) || !all(ctr$indigenous == ctr$x05_02) ||
        !all(ctr$hh_size == ctr$x05_05) || !all(ctr$educ == ctr$x05_07_) ||
        !all(ctr$child_old == ctr$x06_03_pc_)) stop("covariates do not reproduce control.csv")
    kl_old <- .read_raw_num(raw, "pm.csv") %>% select(-1) %>%
      filter(var %in% c("katz", "lawton"), time == 0) %>%
      select(id = x03_01, var, value) %>% pivot_wider(names_from = var, values_from = value) %>%
      right_join(dat %>% filter(id %in% ids117) %>% select(id, katz0, lawton0), by = "id")
    if (!identical(is.na(kl_old$katz), is.na(kl_old$katz0)) ||
        !identical(is.na(kl_old$lawton), is.na(kl_old$lawton0)) ||
        any(kl_old$katz != kl_old$katz0, na.rm = TRUE) ||
        any(kl_old$lawton != kl_old$lawton0, na.rm = TRUE)) stop("Katz or Lawton do not reproduce pm.csv")

    checks <- list(whoqol = chk_who, attendance = chk_att, d_vs_old = chk_d,
                   n_missing_items = n_missing_items)
    if (!quiet) {
      cat(sprintf(paste0("00_data: 120 people (91 invited); samples %d midline, %d endline. Checks ",
                         "passed: WHOQOL, attendance, z and covariates reproduce pm_who.csv, pm.csv ",
                         "and control.csv for the 117 (largest WHOQOL gap %.1e); %d missing items in ",
                         "surveys held.\n"),
                  sum(!is.na(dat$whoqol_y0) & !is.na(dat$whoqol_y1)),
                  sum(!is.na(dat$whoqol_y0) & !is.na(dat$whoqol_y2)),
                  max(chk_who$max_abs_diff), n_missing_items))
    }
    attr(dat, "checks") <- checks
  }

  # ---- 5.7 return: dates only on request, and the check-only column dropped ---------------------
  keep_private <- if (dates) c("date0", "date1", "date2", "entry_date") else character(0)
  out <- dat %>%
    select(id, z, all_of(ycols), d1, d2, days1, days2, weeks1, weeks2, dpw1, dpw2,
           age0, female, indigenous, hh_size, educ, child_carer, katz0, lawton0,
           all_of(keep_private))
  attr(out, "checks") <- attr(dat, "checks")
  out
}

# =============================================================================================
# 6. ONE ANCOVA SAMPLE
# =============================================================================================
# Everyone observed at baseline and at `wave` (1 midline, 2 endline) for one outcome, with the
# columns the per wave ANCOVA and IV need:
#   id, z, d (attended more than five days by that survey), y0 (baseline), y (that wave),
#   weeks, dpw (weeks attended and days per week by that survey)
ancova_sample <- function(dat, outcome, wave) {
  stopifnot(outcome %in% names(OUTCOMES), wave %in% 1:2)
  y0 <- paste0(outcome, "_y0"); yw <- paste0(outcome, "_y", wave)
  dat %>%
    filter(!is.na(.data[[y0]]), !is.na(.data[[yw]])) %>%
    transmute(id, z, d = .data[[paste0("d", wave)]], y0 = .data[[y0]], y = .data[[yw]],
              weeks = .data[[paste0("weeks", wave)]], dpw = .data[[paste0("dpw", wave)]])
}
