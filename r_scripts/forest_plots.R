# =========================
# Packages
# =========================
suppressPackageStartupMessages({
  library(tidyverse)
  library(glue)
  library(ggtext)
  library(patchwork)
})

# =========================
# Paths
# =========================
infile  <- "../data/processed/final_meta_table.csv"
out_dir <- "../results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# =========================
# Helpers
# =========================
num_only <- function(x){
  x <- as.character(x)
  x <- str_replace_all(x, ",", ".")
  as.numeric(str_replace_all(x, "[^0-9.]", ""))
}

map_species <- function(x){
  key <- str_to_lower(str_trim(x))
  case_when(
    str_detect(key, "calif")            ~ "CSL",
    str_detect(key, "steller")          ~ "SSL",
    str_detect(key, "harbo[u]?r seal")  ~ "HS",
    str_detect(key, "fur seal")         ~ "FS",
    str_detect(key, "monk seal")        ~ "MS",
    str_detect(key, "ringed seal")      ~ "RgS",
    str_detect(key, "ribbon seal")      ~ "RbS",
    str_detect(key, "spotted seal")     ~ "SS",
    str_detect(key, "bearded seal")     ~ "BS",
    TRUE ~ str_to_upper(str_extract_all(key, "\\b\\w") %>% sapply(paste0, collapse = ""))
  )
}

map_matrix <- function(x){
  key <- str_to_lower(str_trim(x))
  case_when(
    key %in% c("urine","maternal urine","fetal urine") ~ ifelse(key=="maternal urine","MU", ifelse(key=="fetal urine","FU","Ur")),
    key %in% c("feces","faeces","stool")               ~ "Fe",
    key == "serum"                                     ~ "Se",
    key == "plasma"                                    ~ "Pl",
    key == "blood"                                     ~ "Bd",
    key == "blubber"                                   ~ "Bl",
    key == "kidney"                                    ~ "Ki",
    key == "liver"                                     ~ "Li",
    key == "brain"                                     ~ "Br",
    key == "muscle"                                    ~ "Mu",
    key == "gi"                                        ~ "GI",
    str_detect(key, "stomach content")                 ~ "StC",
    str_detect(key, "colon content")                   ~ "ColC",
    key == "amniotic fluid"                            ~ "Am",
    key == "fetal gastric fluid"                       ~ "FG",
    key == "aqueous humor"                             ~ "Aq",
    key == "abdominal fluid"                           ~ "Ab",
    key == "milk"                                      ~ "Mi",
    TRUE                                               ~ str_to_title(str_sub(key, 1, 2))
  )
}

norm_toxin <- function(x){
  key <- str_to_lower(str_trim(x))
  case_when(
    str_detect(key, "domoic|\\bda\\b")            ~ "DA",
    str_detect(key, "saxitox|\\bstx\\b|paralytic")~ "STX",
    TRUE ~ toupper(key)
  )
}

pretty_10pow_labels <- function(ks){
  sapply(ks, function(k){
    v <- 10^k
    format(v, big.mark = ",", trim = TRUE, scientific = FALSE)
  })
}

# cleaned “one file” parser; toxin_tag is "DA" or "STX"
clean_one <- function(df, toxin_tag){
  df %>%
    rename(
      Publication     = any_of("study"),
      Year            = any_of("year"),
      Toxin           = any_of("toxin"),
      Species         = any_of("species"),
      `Sample source` = any_of(c("Sample source","matrix","Sample")),
      Method          = any_of(c("Method","Detection method","method")),
      Min             = any_of("min"),
      Max             = any_of("max"),
      Stranded        = any_of("stranded")
    ) %>%
    mutate(
      toxin_norm = norm_toxin(Toxin),
      Method     = dplyr::recode(Method, "PSP mouse bioassay" = "PSP MBA")
    ) %>%
    # keep only the requested toxin
    filter(toxin_norm == toxin_tag) %>%
    # drop explicit non-stranding rows (robust “no” checks)
    mutate(Stranded = as.character(Stranded),
           stranded_key = str_to_lower(str_trim(coalesce(Stranded, "")))) %>%
    filter(!str_detect(stranded_key, "^(no|none|0)$")) %>%
    mutate(
      event_year = na_if(as.character(Year), ""),
      sp         = map_species(Species),
      mat        = map_matrix(`Sample source`),
      min_raw    = as.character(Min),
      max_raw    = as.character(Max),
      min_lt     = str_detect(min_raw, "^\\s*<"),
      max_lt     = str_detect(max_raw, "^\\s*<"),
      min_ngg    = num_only(min_raw),
      max_ngg    = num_only(max_raw)
    ) %>%
    filter(!(min_lt & max_lt)) %>%                   # drop ND-only rows
    filter(is.finite(max_ngg)) %>%                   # need numeric Max
    filter(is.na(min_ngg) | max_ngg >= min_ngg) %>% # sanity
    transmute(Publication, event_year, sp, mat,
              method = Method, min_ngg, max_ngg)
}

# Color dictionaries
sp_cols <- c(
  CSL = "#1b9e77", SSL = "#7570b3", HS  = "#d95f02", BS  = "#1f78b4",
  RgS = "#e7298a", SS  = "#a6761d", GS  = "#66a61e", MS  = "#e6ab02", FS  = "#666666"
)

method_colors <- c(
  "ELISA"    = "#1F78B4",
  "DBA"      = "#33A02C",
  "LC-MS/MS" = "#E31A1C",
  "PSP MBA"  = "#6A3D9A"
)

# =========================
# Plot builder for a given toxin (“DA” / “STX”)
# =========================
build_range_forest <- function(in_csv, toxin_tag){
  raw <- readr::read_csv(in_csv, show_col_types = FALSE)
  df_all <- clean_one(raw, toxin_tag)
  
  # Collapse duplicates
  df_plot <- df_all %>%
    group_by(Publication, event_year, sp, mat, method) %>%
    summarise(
      min_ngg = if (all(is.na(min_ngg))) NA_real_ else min(min_ngg, na.rm = TRUE),
      max_ngg = max(max_ngg, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(is.finite(max_ngg))
  
  # Method ordering: ELISA first
  method_order <- c("ELISA", setdiff(sort(unique(df_plot$method)), "ELISA"))
  df_plot <- df_plot %>% mutate(method = factor(method, levels = method_order, ordered = TRUE))
  
  # Labels
  df_plot <- df_plot %>%
    mutate(
      label_plain = glue("{Publication} | {sp}"),
      sp_color    = ifelse(!is.na(sp_cols[sp]), sp_cols[sp], "#333333"),
      label_md    = glue("{Publication} | ",
                         "<span style='color:{sp_color}; font-weight:600'>{sp}</span>")
    ) %>%
    arrange(method, mat, Publication, sp, min_ngg, max_ngg) %>%
    group_by(label_plain) %>%
    mutate(label_id = row_number()) %>%
    ungroup() %>%
    mutate(label_key = paste(label_plain, label_id, sep = "__"))
  
  y_levels <- rev(df_plot$label_key)
  df_plot <- df_plot %>%
    mutate(
      yf    = factor(label_key, levels = y_levels),
      y_idx = as.numeric(yf)
    )
  lab_map <- setNames(df_plot$label_md, df_plot$yf)
  n_rows  <- nlevels(df_plot$yf)
  
  # Panels: left labels
  p_labels <- ggplot(df_plot, aes(y = yf, x = 0)) +
    geom_point(alpha = 0) +
    scale_y_discrete(limits = levels(df_plot$yf), labels = lab_map) +
    scale_x_continuous(expand = c(0,0)) +
    theme_void() +
    theme(axis.text.y = ggtext::element_markdown(size = 9, hjust = 1),
          plot.margin = margin(r = 0))
  
  # Matrix strip
  p_matrix <- ggplot(df_plot, aes(y = as.numeric(yf), x = 1, fill = mat)) +
    geom_tile(width = 1, height = 1) +
    geom_text(aes(label = mat), color = "white", fontface = "bold", size = 3) +
    scale_y_continuous(limits = c(0.5, n_rows + 0.5), expand = c(0, 0)) +
    theme_void() +
    theme(legend.position = "none", plot.margin = margin(l = 0, r = 0))
  
  # Method positions and strip
  method_pos <- df_plot %>%
    group_by(method) %>%
    summarise(y_mid = (min(y_idx) + max(y_idx))/2, .groups = "drop")
  
  p_method <- ggplot(df_plot, aes(y = yf, x = 1, fill = method)) +
    geom_tile(width = 0.32, height = 1, show.legend = FALSE) +
    scale_fill_manual(values = method_colors) +
    geom_text(data = method_pos,
              aes(y = y_mid, x = 1, label = method),
              inherit.aes = FALSE,
              hjust = 0.5, vjust = 0.5,
              color = "white", fontface = "bold", size = 3.2) +
    scale_y_discrete(limits = levels(df_plot$yf), expand = c(0, 0)) +
    theme_void() +
    theme(
      legend.position = "none",
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0),
      plot.title.position = "plot",
      plot.title = element_text(margin = margin(b = 2))
    )
  
  # Transform to log10 once; compute dynamic x range
  df_plot_tx <- df_plot %>%
    mutate(
      x_min = ifelse(is.finite(min_ngg) & min_ngg > 0, log10(min_ngg), NA_real_),
      x_max = ifelse(is.finite(max_ngg) & max_ngg > 0, log10(max_ngg), NA_real_)
    )
  
  x_lo <- floor(min(df_plot_tx$x_min, df_plot_tx$x_max, na.rm = TRUE))
  x_hi <- ceiling(max(df_plot_tx$x_min, df_plot_tx$x_max, na.rm = TRUE))
  x_lo <- max(0, x_lo)   # don’t go below 1 on the linear axis
  x_hi <- min(6, x_hi)   # cap at 10^6 as before
  breaks <- seq(x_lo, x_hi, by = 1)
  labels <- pretty_10pow_labels(breaks)
  
  # Method background bands
  method_bands_tx <- df_plot_tx %>%
    group_by(method) %>%
    summarise(
      ymin = min(as.numeric(yf)) - 0.5,
      ymax = max(as.numeric(yf)) + 0.5,
      .groups = "drop"
    ) %>%
    mutate(xmin = x_lo, xmax = x_hi)
  
  title_txt <- glue("{toxin_tag} concentration ranges by study / matrix / method")
  
  p_main <- ggplot(df_plot_tx, aes(y = yf)) +
    geom_rect(data = method_bands_tx,
              aes(ymin = ymin, ymax = ymax, xmin = -Inf, xmax = Inf, fill = method),
              inherit.aes = FALSE, alpha = 0.20, colour = NA) +
    scale_fill_manual(values = method_colors) +
    guides(fill = "none") +
    geom_segment(aes(x = x_min, xend = x_max, yend = yf),
                 linewidth = 0.6, colour = "grey60", na.rm = TRUE) +
    geom_point(aes(x = x_min), size = 2.6, colour = "blue", na.rm = TRUE) +
    geom_point(aes(x = x_max), size = 2.6, colour = "red",  na.rm = TRUE) +
    scale_y_discrete(limits = levels(df_plot$yf), labels = NULL) +
    scale_x_continuous(
      name   = "Concentration (ng/g or ml, log10)",
      breaks = breaks,
      labels = labels,
      limits = c(x_lo, x_hi),
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    theme_minimal(base_size = 12) +
    theme(axis.title.y = element_blank(),
          axis.text.y  = element_blank(),
          axis.ticks.y = element_blank(),
          plot.margin  = margin(l = 0, r = 0)) +
    labs(title = title_txt)
  
  # Assemble
  (p_labels | p_matrix | p_method | p_main) +
    plot_layout(widths = c(0.01, 0.08, 0.08, 1))
}

# =========================
# Build & Save: DA and STX
# =========================
p_DA  <- build_range_forest(infile, "DA")
p_STX <- build_range_forest(infile, "STX")

ggsave(filename = file.path(out_dir, "forest_DA_methods.svg"),
       plot = p_DA, device = "svg", width = 12, height = 8)

ggsave(filename = file.path(out_dir, "forest_STX_methods.svg"),
       plot = p_STX, device = "svg", width = 12, height = 8)

message("✅ Wrote: ",
        file.path(out_dir, "forest_DA_methods.svg"), " and ",
        file.path(out_dir, "forest_STX_methods.svg"))
