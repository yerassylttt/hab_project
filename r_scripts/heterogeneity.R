# --- Packages ---------------------------------------------------------------
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(stringi)
  library(metafor)
})

# --- Paths ------------------------------------------------------------------
infile     <- "../data/processed/final_meta_table.csv"
out_dir    <- "../results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- Helpers ----------------------------------------------------------------
fix_num <- function(x) {
  x <- as.character(x)
  x <- stringi::stri_replace_all_regex(x, "\\p{Zs}+", "")  # remove Unicode spaces
  x <- gsub(",", ".", x, fixed = TRUE)                     # decimal comma -> dot
  suppressWarnings(as.numeric(x))
}

map_matrix <- function(x) {
  key <- tolower(str_trim(x))
  dplyr::case_when(
    str_detect(key, "^(gi)\\b|gastro|stomach|colon") ~ "GI contents",
    str_detect(key, "fec|faec|scat|mecon")           ~ "Feces",
    str_detect(key, "urine")                         ~ "Urine",
    str_detect(key, "serum|plasma|blood")            ~ "Serum/Plasma",
    str_detect(key, "liver")                         ~ "Liver",
    str_detect(key, "kidney")                        ~ "Kidney",
    str_detect(key, "bile")                          ~ "Bile",
    str_detect(key, "milk")                          ~ "Milk",
    str_detect(key, "amniot|allanto")                ~ "Fetal fluids",
    str_detect(key, "other|unknown|misc")            ~ "Other",
    TRUE                                            ~ "Other"
  )
}

map_region <- function(x) {
  key <- tolower(str_trim(x))
  dplyr::case_when(
    str_detect(key, "\\bcalifornia\\b|central california|san[- ]?miguel") ~ "California",
    str_detect(key, "bering.*chukchi|\\bchukchi\\b|\\bbering sea\\b")     ~ "Arctic (NE Pacific)",
    str_detect(key, "northern alaska|\\balaska\\b")                        ~ "Arctic (NE Pacific)",
    str_detect(key, "\\bscotland\\b")                                      ~ "NE Atlantic",
    str_detect(key, "kattegat|kattegart|denmark")                          ~ "NE Atlantic",
    str_detect(key, "st\\.?\\s*lawrence|quebec|quebeck")                   ~ "NW Atlantic",
    str_detect(key, "western sahara|cap\\s*blanc")                         ~ "West Africa",
    TRUE                                                                   ~ "Other/Unknown"
  )
}

bt_overall <- function(fit) {
  if (is.null(fit)) return(list(pooled = NA_real_, ci.lb = NA_real_, ci.ub = NA_real_))
  p <- predict(fit, transf = exp)
  list(pooled = as.numeric(p$pred), ci.lb = as.numeric(p$ci.lb), ci.ub = as.numeric(p$ci.ub))
}

get_I2_mv <- function(fit) {
  if (is.null(fit)) return(NA_real_)
  out <- try(i2(fit), silent = TRUE)
  if (!inherits(out, "try-error") && !is.null(out$I2)) return(as.numeric(out$I2))
  tau2 <- suppressWarnings(fit$sigma2[1]); vbar <- mean(fit$vi)
  as.numeric(tau2 / (tau2 + vbar))
}

safe_tau2_ci <- function(fit) {
  out <- c(NA_real_, NA_real_)
  ci_all <- try(confint(fit), silent = TRUE)
  if (!inherits(ci_all, "try-error")) {
    rand <- try(ci_all[["random"]], silent = TRUE)
    if (!inherits(rand, "try-error") && !is.null(rand)) {
      if (is.data.frame(rand) && all(c("ci.lb", "ci.ub") %in% names(rand))) {
        out <- c(as.numeric(rand$ci.lb[1]), as.numeric(rand$ci.ub[1]))
      } else if (is.matrix(rand) && all(c("ci.lb", "ci.ub") %in% colnames(rand))) {
        out <- c(as.numeric(rand[1, "ci.lb"]), as.numeric(rand[1, "ci.ub"]))
      }
    }
  }
  out
}

fit_one <- function(df) {
  if (nrow(df) < 2) return(NULL)
  try(rma.mv(yi, vi, random = ~ 1 | study, data = df, method = "REML"), silent = TRUE)
}

# --- Load & Clean -----------------------------------------------------------
message("Reading: ", infile)
dat_raw <- read_csv(infile, show_col_types = FALSE)

dat <- dat_raw %>%
  mutate(
    min  = fix_num(min),
    max  = fix_num(max),
    mean = fix_num(mean),
    sd   = fix_num(sd),
    n    = suppressWarnings(as.integer(n))
  ) %>%
  filter(!is.na(mean), !is.na(sd), !is.na(n), n >= 3)

dat_log <- dat %>%
  mutate(
    mean_log   = log(mean),
    min_log    = ifelse(min > 0, log(min), NA_real_),
    max_log    = ifelse(max > 0, log(max), NA_real_),
    sd_log     = sd / mean,                 # coefficient of variation
    matrix_grp = map_matrix(matrix)
  )

# --- BY TOXIN ---------------------------------------------------------------
es_tox <- dat_log %>%
  transmute(
    toxin = as.factor(toxin),
    study = as.factor(study),
    yi    = mean_log,
    vi    = (sd^2) / (n * (mean^2)),
    n     = n
  ) %>%
  filter(is.finite(yi), is.finite(vi), vi > 0)

levels_tox <- levels(droplevels(es_tox$toxin))
fits <- lapply(levels_tox, function(tx) fit_one(dplyr::filter(es_tox, toxin == tx)))
names(fits) <- levels_tox

het_row <- function(fit, tx) {
  if (is.null(fit) || inherits(fit, "try-error")) {
    return(dplyr::tibble(
      toxin = tx, tau2 = NA_real_, tau2_ci_lb = NA_real_, tau2_ci_ub = NA_real_,
      I2 = NA_real_, H2 = NA_real_, k = 0,
      pooled = NA_real_, pooled_ci_lb = NA_real_, pooled_ci_ub = NA_real_
    ))
  }
  ci_tau <- safe_tau2_ci(fit)
  I2     <- get_I2_mv(fit)
  bt     <- bt_overall(fit)
  dplyr::tibble(
    toxin        = tx,
    tau2         = as.numeric(fit$sigma2[1]),
    tau2_ci_lb   = ci_tau[1],
    tau2_ci_ub   = ci_tau[2],
    I2           = I2,
    H2           = ifelse(is.na(I2), NA_real_, 1/(1 - I2)),
    k            = fit$k,
    pooled       = bt$pooled,
    pooled_ci_lb = bt$ci.lb,
    pooled_ci_ub = bt$ci.ub
  )
}

heterogeneity_by_toxin <- dplyr::bind_rows(mapply(het_row, fits, names(fits), SIMPLIFY = FALSE))
readr::write_csv(heterogeneity_by_toxin, file.path(out_dir, "heterogeneity_by_toxin.csv"), na = "")
message("Wrote: ", file.path(out_dir, "heterogeneity_by_toxin.csv"))

# --- BY TOXIN × MATRIX ------------------------------------------------------
es_mat <- dat_log %>%
  transmute(
    toxin  = as.factor(toxin),
    study  = as.factor(study),
    matrix = as.factor(matrix_grp),
    yi     = mean_log,
    vi     = (sd^2) / (n * (mean^2))
  ) %>%
  filter(is.finite(yi), is.finite(vi), vi > 0)

heterogeneity_by_matrix <- es_mat %>%
  group_by(toxin, matrix) %>%
  group_modify(~{
    fit <- fit_one(.x)
    if (is.null(fit) || inherits(fit, "try-error")) {
      dplyr::tibble(
        k = nrow(.x),
        n_studies = dplyr::n_distinct(.x$study),
        tau2 = NA_real_, tau2_ci_lb = NA_real_, tau2_ci_ub = NA_real_,
        I2 = NA_real_, H2 = NA_real_,
        pooled = NA_real_, pooled_ci_lb = NA_real_, pooled_ci_ub = NA_real_
      )
    } else {
      ci_tau <- safe_tau2_ci(fit)
      I2     <- get_I2_mv(fit)
      bt     <- bt_overall(fit)
      dplyr::tibble(
        k            = fit$k,
        n_studies    = dplyr::n_distinct(.x$study),
        tau2         = as.numeric(fit$sigma2[1]),
        tau2_ci_lb   = ci_tau[1],
        tau2_ci_ub   = ci_tau[2],
        I2           = I2,
        H2           = ifelse(is.na(I2), NA_real_, 1/(1 - I2)),
        pooled       = bt$pooled,
        pooled_ci_lb = bt$ci.lb,
        pooled_ci_ub = bt$ci.ub
      )
    }
  }) %>%
  ungroup() %>%
  arrange(toxin, matrix)

readr::write_csv(heterogeneity_by_matrix, file.path(out_dir, "heterogeneity_by_matrix.csv"), na = "")
message("Wrote: ", file.path(out_dir, "heterogeneity_by_matrix.csv"))

# --- BY TOXIN × REGION ------------------------------------------------------
# map region once on the long table to allow inspection if needed
dat_reg <- dat_log %>% mutate(region_grp = map_region(location))

# (Optional) Inspect mapping
# dat_reg %>% distinct(location, region_grp) %>% arrange(region_grp, location)

es_reg <- dat_reg %>%
  transmute(
    toxin   = as.factor(toxin),
    study   = as.factor(study),
    region  = as.factor(region_grp),
    yi      = mean_log,
    vi      = (sd^2) / (n * (mean^2))
  ) %>%
  filter(is.finite(yi), is.finite(vi), vi > 0)

heterogeneity_by_region <- es_reg %>%
  group_by(toxin, region) %>%
  group_modify(~{
    fit <- fit_one(.x)
    if (is.null(fit) || inherits(fit, "try-error")) {
      dplyr::tibble(
        k = nrow(.x),
        n_studies = dplyr::n_distinct(.x$study),
        tau2 = NA_real_, tau2_ci_lb = NA_real_, tau2_ci_ub = NA_real_,
        I2 = NA_real_, H2 = NA_real_,
        pooled = NA_real_, pooled_ci_lb = NA_real_, pooled_ci_ub = NA_real_
      )
    } else {
      ci_tau <- safe_tau2_ci(fit)
      I2     <- get_I2_mv(fit)
      bt     <- bt_overall(fit)
      dplyr::tibble(
        k            = fit$k,
        n_studies    = dplyr::n_distinct(.x$study),
        tau2         = as.numeric(fit$sigma2[1]),
        tau2_ci_lb   = ci_tau[1],
        tau2_ci_ub   = ci_tau[2],
        I2           = I2,
        H2           = ifelse(is.na(I2), NA_real_, 1/(1 - I2)),
        pooled       = bt$pooled,
        pooled_ci_lb = bt$ci.lb,
        pooled_ci_ub = bt$ci.ub
      )
    }
  }) %>%
  ungroup() %>%
  arrange(toxin, region)

readr::write_csv(heterogeneity_by_region, file.path(out_dir, "heterogeneity_by_region.csv"), na = "")
message("Wrote: ", file.path(out_dir, "heterogeneity_by_region.csv"))

message("✅ Done.")
