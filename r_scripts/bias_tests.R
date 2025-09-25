# =========================
# Packages
# =========================
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(stringi)
  library(metafor)
  library(ggplot2)  # only for device checks; base plotting used for funnel()
})

# =========================
# Paths
# =========================
infile   <- "../data/processed/final_meta_table.csv"
out_dir  <- "../results"
plot_dir <- file.path(out_dir, "funnel_plots_svg")

dir.create(out_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# =========================
# Helpers
# =========================
fix_num <- function(x) {
  x <- as.character(x)
  x <- stringi::stri_replace_all_regex(x, "\\p{Zs}+", "")  # remove Unicode spaces
  x <- gsub(",", ".", x, fixed = TRUE)                     # decimal comma -> dot
  suppressWarnings(as.numeric(x))
}

map_matrix <- function(x){
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

map_region <- function(x){
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

# Egger’s test for grouped data frames
egger_by_group <- function(es_df, group_vars = NULL) {
  stopifnot(all(c("toxin", "yi", "vi") %in% names(es_df)))
  es_df %>%
    group_by(across(all_of(c("toxin", group_vars)))) %>%
    group_modify(~{
      if (nrow(.x) < 3) {
        return(tibble(
          k = nrow(.x),
          egger_intercept = NA_real_,
          egger_se        = NA_real_,
          egger_z         = NA_real_,
          egger_p         = NA_real_
        ))
      }
      fit <- try(rma(yi, vi, data = .x, method = "REML"), silent = TRUE)
      if (inherits(fit, "try-error")) {
        return(tibble(
          k = nrow(.x),
          egger_intercept = NA_real_,
          egger_se        = NA_real_,
          egger_z         = NA_real_,
          egger_p         = NA_real_
        ))
      }
      test <- try(regtest(fit, model = "rma"), silent = TRUE)
      if (inherits(test, "try-error")) {
        return(tibble(
          k = nrow(.x),
          egger_intercept = NA_real_,
          egger_se        = NA_real_,
          egger_z         = NA_real_,
          egger_p         = NA_real_
        ))
      }
      tibble(
        k               = nrow(.x),
        egger_intercept = unname(test$beta[1]),
        egger_se        = unname(test$se[1]),
        egger_z         = unname(test$zval[1]),
        egger_p         = unname(test$pval[1])
      )
    }) %>%
    ungroup() %>%
    arrange(across(all_of(c("toxin", group_vars))))
}

# Funnel plot for one toxin; saves SVG file if save_dir provided
plot_funnel_for_toxin <- function(df, toxin_name, save_dir = NULL) {
  subdat <- df %>% filter(toxin == toxin_name)
  if (nrow(subdat) < 3) {
    message("Skipping ", toxin_name, " — fewer than 3 rows")
    return(invisible(NULL))
  }
  fit <- try(rma(yi, vi, data = subdat, method = "REML"), silent = TRUE)
  if (inherits(fit, "try-error")) {
    message("Model failed for ", toxin_name)
    return(invisible(NULL))
  }
  
  # open device if saving
  if (!is.null(save_dir)) {
    svg(filename = file.path(save_dir, paste0("funnel_", toxin_name, ".svg")),
        width = 6, height = 6)
    on.exit(dev.off(), add = TRUE)
  }
  
  # add Egger p to title if available
  egger_p <- tryCatch(regtest(fit, model = "rma")$pval[1], error = function(e) NA_real_)
  ttl <- if (is.finite(egger_p)) {
    sprintf("Funnel plot for %s (Egger p = %.3f)", toxin_name, egger_p)
  } else {
    sprintf("Funnel plot for %s", toxin_name)
  }
  
  funnel(fit,
         xlab   = "Effect size (log scale)",
         ylab   = "Standard Error",
         refline = fit$b,
         shade  = c("white", "gray95", "gray90"),
         level  = c(90, 95, 99),
         main   = ttl)
  abline(v = fit$b, col = "red", lwd = 2, lty = 2)
  invisible(NULL)
}

# =========================
# Load & prepare effect sizes
# =========================
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
    min_log    = ifelse(is.na(min) | min <= 0, NA_real_, log(min)),
    max_log    = ifelse(is.na(max) | max <= 0, NA_real_, log(max)),
    sd_log     = sd / mean,
    matrix_grp = map_matrix(matrix),
    region_grp = map_region(location)
  )

# Base effect sizes (by study, toxin)
es_base <- dat_log %>%
  transmute(
    toxin = as.factor(toxin),
    study = as.factor(study),
    yi    = mean_log,
    vi    = (sd^2) / (n * (mean^2))
  ) %>%
  filter(is.finite(yi), is.finite(vi), vi > 0)

# Matrix-annotated
es_matrix <- dat_log %>%
  transmute(
    toxin  = as.factor(toxin),
    study  = as.factor(study),
    matrix = as.factor(matrix_grp),
    yi     = mean_log,
    vi     = (sd^2) / (n * (mean^2))
  ) %>%
  filter(is.finite(yi), is.finite(vi), vi > 0)

# Region-annotated
es_region <- dat_log %>%
  transmute(
    toxin  = as.factor(toxin),
    study  = as.factor(study),
    region = as.factor(region_grp),
    yi     = mean_log,
    vi     = (sd^2) / (n * (mean^2))
  ) %>%
  filter(is.finite(yi), is.finite(vi), vi > 0)

# =========================
# Egger’s tests (CSV outputs)
# =========================
egger_by_toxin  <- egger_by_group(es_base)
egger_by_matrix <- egger_by_group(es_matrix, group_vars = "matrix")
egger_by_region <- egger_by_group(es_region, group_vars = "region")

write_csv(egger_by_toxin,  file.path(out_dir, "eggers_by_toxin.csv"),  na = "")
write_csv(egger_by_matrix, file.path(out_dir, "eggers_by_matrix.csv"), na = "")
write_csv(egger_by_region, file.path(out_dir, "eggers_by_region.csv"), na = "")

message("Wrote: ",
        file.path(out_dir, "eggers_by_toxin.csv"), ", ",
        file.path(out_dir, "eggers_by_matrix.csv"), ", ",
        file.path(out_dir, "eggers_by_region.csv"))

# =========================
# Funnel plots by toxin (SVGs)
# =========================
toxins <- levels(droplevels(es_base$toxin))
invisible(lapply(toxins, function(tx) plot_funnel_for_toxin(es_base, tx, save_dir = plot_dir)))

message("✅ Funnel plots written to: ", normalizePath(plot_dir))
message("✅ Done.")
