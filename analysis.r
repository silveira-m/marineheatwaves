## ============================================================
## Marine Heatwave (MHW) Time Series Analysis — Portugal (OISST, 1981-2025)
## Standalone R script — analyses ALL 6 variables in the long CSV:
##   activity_degC_days_m2, number_of_events, sum_area_km2,
##   mean_intensity_degC, mean_duration_days, mean_area_km2
## Each variable gets its own results folder "<variable>_results/"
## with the same STL/GESD, annual, breakpoint, seasonal, monthly,
## and DFA outputs as the original single-variable script.
##
## FIX (this version): only the three SUM-type metrics
## (activity_degC_days_m2, number_of_events, sum_area_km2) are zero-filled
## for months with no MHW event. The three MEAN-type metrics
## (mean_intensity_degC, mean_duration_days, mean_area_km2) are left as NA
## in those months, since they are conditional averages ("average
## intensity/duration/area of events that occurred") and not totals -- a
## month with no event has no such average to report, so it must be
## excluded (na.rm = TRUE) rather than treated as a value of 0. Previously
## all six variables were zero-filled, which biased every annual/seasonal/
## monthly/STL summary of the three mean-type metrics toward 0 in years or
## seasons with few event-months (e.g. mean_duration_days could appear
## below its true floor of min_duration_mhw).
## ============================================================

## ---- Setup ----
required_pkgs <- c(
  "tidyverse", "cowplot", "anomalize", "strucchange", "MARSS",
  "broom", "Kendall", "trend", "zoo", "ggh4x"
)
missing_pkgs <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(missing_pkgs) > 0) {
  install.packages(missing_pkgs, repos = "https://cloud.r-project.org")
}
invisible(lapply(required_pkgs, library, character.only = TRUE))

# Replacement for archived ggseas::tsdf()
tsdf <- function(timeseries, colname = "x") {
  out <- as.data.frame(timeseries)
  if (is.null(colnames(timeseries))) names(out) <- "y"
  out[[colname]] <- as.numeric(stats::time(timeseries))
  out[c(colname, setdiff(names(out), colname))]
}

# For the three MEAN-type metrics, a period (month or year) with zero MHW
# events is a genuine NA ("no events -> no average to report"), not a 0 --
# see the note above data_wide. Most of the pipeline (mean(), lm(), and
# Kendall::MannKendall()) already tolerates that NA correctly via internal
# na.rm/na.omit handling. A few functions do NOT tolerate any NA at all --
# stl(), trend::sens.slope(), and strucchange::breakpoints() -- and will
# error out (na.fail) if given a series with gaps. For those specific
# calls only, fill_gaps() linearly interpolates internal NA runs and
# carries the nearest valid value out to the series edges, purely so the
# function has a complete series to run on. This never touches the
# NA-preserving data written to CSV or drawn in geom_line() plots -- those
# correctly show the true gap.
fill_gaps <- function(x) {
  if (!anyNA(x)) return(x)
  x %>%
    zoo::na.approx(na.rm = FALSE) %>%
    zoo::na.locf(na.rm = FALSE) %>%
    zoo::na.locf(fromLast = TRUE, na.rm = FALSE)
}

theme_data <- theme_minimal(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey92", color = NA),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )
theme_set(theme_data)

palette_zone <- c(NW = "#1b9e77", SW = "#d95f02", S = "#7570b3")

linear_label <- function(model) {
  co <- coef(model)
  s  <- summary(model)
  sprintf("slope = %s, R^2 = %s, p = %s",
          signif(co[2], 3), signif(s$r.squared, 3), signif(s$coefficients[2, 4], 3))
}

## ---- Load data (once per bathymetry file) ----
# >>> EDIT THIS: list every bathymetry CSV you want processed in this run.
# Each file gets its own results folder (named after the bathymetry depth
# detected in its filename) under output_base_dir. Add or remove paths as
# needed -- everything below runs once per entry in this vector.
data_files <- c(
  "/Users/miguelsilveira/Documents/GitHub/marineheatwaves/DATA/OUTPUT/MHW_summary_table_long_oisst.19812025.Portugal_1000m_bathymetry.csv"
  # , "/path/to/MHW_summary_table_long_oisst.19812025.Portugal_200m_bathymetry.csv"
  # , "/path/to/MHW_summary_table_long_oisst.19812025.Portugal_500m_bathymetry.csv"
)
stopifnot(all(file.exists(data_files)))

# Only the three SUM-type event-summary metrics are zero-filled for months
# with no MHW event: a month with no event genuinely contributed 0 total
# activity, 0 events, and 0 km^2 of area, so 0 is the correct value there.
#
# The three MEAN-type metrics (mean_intensity_degC, mean_duration_days,
# mean_area_km2) are deliberately left as NA in months with no event: they
# are conditional averages ("average intensity/duration/area of the events
# that occurred"), not totals, so a month with no event has no such value
# to report. na.rm = TRUE in every downstream summarise()/mean() call then
# correctly excludes those months instead of treating them as 0.
sum_type_vars  <- c("activity_degC_days_m2", "number_of_events", "sum_area_km2")
mean_type_vars <- c("mean_intensity_degC", "mean_duration_days", "mean_area_km2")

# >>> EDIT THIS to control where every results_dir/ folder is created <<<
# All output (CSVs and PNGs) is written under:
#   <output_base_dir>/<bathymetry_tag>/<variable>_results/
# One subfolder per bathymetry file in data_files above. Absolute path so
# output always lands in the same place regardless of the R session's
# working directory (getwd()).
output_base_dir <- "/Users/miguelsilveira/Documents/GitHub/marineheatwaves/TSA results"
if (!dir.exists(output_base_dir)) dir.create(output_base_dir, recursive = TRUE)
message("Output base directory: ", normalizePath(output_base_dir))

# Loads and prepares one bathymetry CSV. Returns a list with the prepared
# data_wide tibble and the bathymetry_tag detected from its filename.
load_bathymetry_file <- function(data_file) {
  raw <- read_csv(data_file, show_col_types = FALSE)
  stopifnot(all(c("Year", "Month", "Zone", "Parameter", "Value") %in% names(raw)))

  # Pull the bathymetry depth out of the filename (e.g. "...1000m_bathymetry.csv"
  # -> "1000m") so results from different bathymetry cutoffs land in separate
  # folders instead of overwriting each other. Falls back to "unknown_depth"
  # if the filename doesn't follow the expected "<depth>m_bathymetry" pattern.
  bathymetry_tag <- str_extract(data_file, "\\d+m(?=_bathymetry)")
  if (is.na(bathymetry_tag)) bathymetry_tag <- "unknown_depth"
  message("Bathymetry tag detected: ", bathymetry_tag, "  (", basename(data_file), ")")

  raw <- raw %>%
    mutate(
      Zone = factor(Zone, levels = c("NW", "SW", "S")),
      date = as.Date(sprintf("%d-%02d-01", Year, Month)),
      month_name = factor(month.abb[Month], levels = month.abb),
      season = factor(
        case_when(
          Month %in% c(12, 1, 2)  ~ "DJF",
          Month %in% c(3, 4, 5)   ~ "MAM",
          Month %in% c(6, 7, 8)   ~ "JJA",
          Month %in% c(9, 10, 11) ~ "SON"
        ),
        levels = c("DJF", "MAM", "JJA", "SON")
      )
    )

  data_wide <- raw %>%
    select(Year, Month, date, month_name, season, Zone, Parameter, Value) %>%
    pivot_wider(names_from = Parameter, values_from = Value) %>%
    arrange(Zone, date) %>%
    mutate(
      across(all_of(sum_type_vars), ~ replace_na(.x, 0))
      # mean_type_vars are intentionally NOT zero-filled here -- see note above.
    )

  list(data_wide = data_wide, bathymetry_tag = bathymetry_tag)
}

## ---- Variable configuration ----
## agg_fun: how to aggregate the monthly value into an annual/seasonal/
##   monthly-climatology figure. "sum" = totals over the period (activity,
##   event counts, total area). "mean" = averages over the period (mean
##   intensity, mean duration, mean area) — summing these wouldn't be
##   physically meaningful.
all_variables <- c(
  "activity_degC_days_m2", "number_of_events", "sum_area_km2",
  "mean_intensity_degC", "mean_duration_days", "mean_area_km2"
)

agg_fun_for <- function(var) {
  if (var %in% sum_type_vars) "sum" else "mean"
}

y_label_for <- function(var) {
  switch(var,
    "activity_degC_days_m2" = expression(degree*C~days~m^-2),
    "number_of_events"      = "Number of events",
    "sum_area_km2"          = expression(Total~area~(km^2)),
    "mean_intensity_degC"   = expression(Mean~intensity~(degree*C)),
    "mean_duration_days"    = "Mean duration (days)",
    "mean_area_km2"         = expression(Mean~area~(km^2))
  )
}

display_name_for <- function(var) {
  switch(var,
    "activity_degC_days_m2" = "MHW Activity",
    "number_of_events"      = "MHW Event Count",
    "sum_area_km2"          = "MHW Total Area",
    "mean_intensity_degC"   = "MHW Mean Intensity",
    "mean_duration_days"    = "MHW Mean Duration",
    "mean_area_km2"         = "MHW Mean Area"
  )
}

t_window <- 133
alpha    <- 0.05

## ---- Full analysis pipeline for one variable ----
analyse_variable <- function(var) {

  message("== Analysing: ", var, " ==")

  agg_fun     <- agg_fun_for(var)
  y_lab       <- y_label_for(var)
  disp_name   <- display_name_for(var)
  results_dir <- file.path(output_base_dir, bathymetry_tag, paste0(var, "_results"))
  if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

  # Generic working frame: rename target column to "value" so the rest of
  # the pipeline is identical regardless of which variable is being run.
  dw <- data_wide %>% rename(value = all_of(var))

  ## ---- STL decomposition + GESD anomaly detection ----
  ## NOTE: for the three mean-type variables, `value` may now contain NA
  ## in months with no event (see fix above). ts() will carry these NAs
  ## through, and stl() requires a series with no NAs. We linearly
  ## interpolate purely for the purposes of STL/anomaly detection (a
  ## smooth trend/seasonal decomposition needs a continuous series) while
  ## all other outputs (annual, seasonal, monthly summaries) continue to
  ## use na.rm = TRUE on the original NA-preserving data.
  run_stl_for_zone <- function(zone) {
    df_zone <- dw %>% filter(Zone == zone) %>% arrange(date)
    ts_zone <- ts(fill_gaps(df_zone$value),
                  start = c(df_zone$Year[1], df_zone$Month[1]), frequency = 12)

    ts_zone %>%
      stl(s.window = "periodic", t.window = t_window, robust = TRUE) %>%
      .$time.series %>%
      tsdf() %>%
      as_tibble() %>%
      anomalize(remainder, alpha = alpha, method = "gesd") %>%
      mutate(observed = trend + remainder + seasonal, Zone = zone)
  }

  data_stl <- map_dfr(levels(dw$Zone), run_stl_for_zone) %>%
    mutate(Zone = factor(Zone, levels = levels(dw$Zone)))

  write_csv(data_stl, file.path(results_dir, paste0(var, "_stl_anomalies.csv")))

  plot_stl <- data_stl %>%
    pivot_longer(cols = c(observed, trend, seasonal, remainder),
                 names_to = "component", values_to = "val") %>%
    mutate(component = factor(component, levels = c("observed", "trend", "seasonal", "remainder")))

  obs_panel <- ggplot(filter(plot_stl, component == "observed"), aes(x, val, color = Zone)) +
    geom_line() +
    geom_point(data = filter(plot_stl, component == "observed", anomaly == "Yes"),
               shape = 21, size = 2.5, fill = "white") +
    facet_grid2(component ~ Zone, scales = "free_y", independent = "y") +
    scale_color_manual(values = palette_zone) +
    labs(title = paste0("Observed ", disp_name, " & anomalies"), x = NULL, y = y_lab) +
    theme(legend.position = "none")

  decomp_panel <- ggplot(filter(plot_stl, component != "observed"), aes(x, val, color = Zone)) +
    geom_line() +
    geom_point(data = filter(plot_stl, component == "remainder", anomaly == "Yes"), size = 1.5) +
    facet_grid2(component ~ Zone, scales = "free_y", independent = "y") +
    scale_color_manual(values = palette_zone) +
    labs(x = NULL, y = NULL) +
    theme(legend.position = "none", strip.text.x = element_blank())

  stl_combo <- plot_grid(obs_panel, decomp_panel, ncol = 1, rel_heights = c(1.2, 2))
  print(stl_combo)
  ggsave(file.path(results_dir, paste0(var, "_stl.png")), stl_combo,
         width = 260, height = 220, units = "mm", dpi = 400)

  ## ---- Annual means/totals, trends & breakpoints ----
  data_annual <- dw %>%
    group_by(Zone, Year) %>%
    summarise(
      value_annual = if (agg_fun == "sum") sum(value, na.rm = TRUE) else mean(value, na.rm = TRUE),
      .groups = "drop"
    )

  write_csv(data_annual, file.path(results_dir, paste0(var, "_annual_summary.csv")))

  annual_trend_tests <- data_annual %>%
    group_by(Zone) %>%
    summarise(
      # Years the mean-type metrics had zero events all year (true NA,
      # gap-filled below only for the tests that require a complete series).
      n_years_no_events = sum(is.na(value_annual)),
      lm_slope  = coef(lm(value_annual ~ Year))[2],           # lm() drops NA rows on its own
      lm_p      = summary(lm(value_annual ~ Year))$coefficients[2, 4],
      lm_r2     = summary(lm(value_annual ~ Year))$r.squared,
      mk_tau    = Kendall::MannKendall(value_annual)$tau,     # MannKendall() tolerates NA on its own
      mk_p      = Kendall::MannKendall(value_annual)$sl,
      sen_slope = trend::sens.slope(fill_gaps(value_annual))$estimates,  # sens.slope() has no NA tolerance
      .groups = "drop"
    )
  write_csv(annual_trend_tests, file.path(results_dir, paste0(var, "_annual_trend_tests.csv")))
  print(annual_trend_tests)

  annual_plot <- ggplot(data_annual, aes(Year, value_annual, color = Zone, fill = Zone)) +
    geom_smooth(method = "lm", alpha = 0.15) +
    geom_line() +
    geom_point(size = 1.2) +
    facet_grid(Zone ~ ., scales = "free_y") +
    scale_color_manual(values = palette_zone) +
    scale_fill_manual(values = palette_zone) +
    labs(title = paste0("Annual ", ifelse(agg_fun == "sum", "total ", "mean "), disp_name),
         x = NULL, y = y_lab) +
    theme(legend.position = "none")
  print(annual_plot)
  ggsave(file.path(results_dir, paste0(var, "_annual.png")), annual_plot,
         width = 150, height = 170, units = "mm", dpi = 400)

  find_breaks <- function(zone) {
    df_zone <- data_annual %>% filter(Zone == zone) %>% arrange(Year)
    # breakpoints() has no NA tolerance (like sens.slope above), so feed it
    # the gap-filled series; the years/values plotted below still come from
    # the real, NA-preserving df_zone$value_annual.
    ts_zone <- ts(fill_gaps(df_zone$value_annual), start = min(df_zone$Year), frequency = 1)
    bp <- tryCatch(strucchange::breakpoints(ts_zone ~ 1), error = function(e) NULL)

    out <- tibble(Year = df_zone$Year, value = df_zone$value_annual, Zone = zone, breakpoint = FALSE)
    if (!is.null(bp) && length(bp$breakpoints) > 0 && !anyNA(bp$breakpoints)) {
      bp_years <- df_zone$Year[bp$breakpoints]
      out$breakpoint[out$Year %in% bp_years] <- TRUE
    }
    out
  }

  data_bp <- map_dfr(levels(dw$Zone), find_breaks) %>%
    mutate(Zone = factor(Zone, levels = levels(dw$Zone)))

  bp_plot <- ggplot(data_bp, aes(Year, value)) +
    geom_line(color = "grey60") +
    geom_vline(data = filter(data_bp, breakpoint), aes(xintercept = Year),
               color = "firebrick", linetype = "dashed") +
    geom_text(data = filter(data_bp, breakpoint), aes(x = Year, label = Year),
              y = Inf, angle = 90, vjust = -0.4, hjust = 1.1, size = 3, color = "firebrick") +
    facet_grid(Zone ~ ., scales = "free_y") +
    labs(title = paste0("Trend changes / breakpoints — ", disp_name), x = NULL, y = y_lab)
  print(bp_plot)
  ggsave(file.path(results_dir, paste0(var, "_annual_breakpoints.png")), bp_plot,
         width = 150, height = 170, units = "mm", dpi = 400)

  ## ---- Seasonal means (climatological — always mean, regardless of agg_fun) ----
  data_seasonal <- dw %>%
    group_by(Zone, season, Year) %>%
    summarise(value_mean = mean(value, na.rm = TRUE), .groups = "drop")

  write_csv(
    data_seasonal %>% pivot_wider(names_from = Zone, values_from = value_mean),
    file.path(results_dir, paste0(var, "_seasonal_summary.csv"))
  )

  seasonal_trend_tests <- data_seasonal %>%
    group_by(Zone, season) %>%
    summarise(
      lm_slope = coef(lm(value_mean ~ Year))[2],
      lm_p     = summary(lm(value_mean ~ Year))$coefficients[2, 4],
      .groups = "drop"
    )
  write_csv(seasonal_trend_tests, file.path(results_dir, paste0(var, "_seasonal_trend_tests.csv")))

  seasonal_box <- ggplot(dw, aes(x = season, y = value, color = Zone)) +
    geom_boxplot() +
    facet_grid(Zone ~ .) +
    scale_color_manual(values = palette_zone) +
    labs(title = paste0("Seasonal distribution of monthly ", disp_name), x = NULL, y = y_lab) +
    theme(legend.position = "none")

  seasonal_trend_plot <- ggplot(data_seasonal, aes(Year, value_mean, color = Zone, fill = Zone)) +
    geom_point(alpha = 0.3) +
    geom_smooth(method = "lm", se = TRUE) +
    facet_grid(Zone ~ season, scales = "free_y") +
    scale_color_manual(values = palette_zone) +
    scale_fill_manual(values = palette_zone) +
    labs(title = paste0("Seasonal trends — ", disp_name), x = NULL, y = y_lab) +
    theme(legend.position = "none")

  print(seasonal_box)
  print(seasonal_trend_plot)
  ggsave(file.path(results_dir, paste0(var, "_seasonal_box.png")), seasonal_box,
         width = 183, height = 150, units = "mm", dpi = 400)
  ggsave(file.path(results_dir, paste0(var, "_seasonal_trend.png")), seasonal_trend_plot,
         width = 183, height = 150, units = "mm", dpi = 400)

  ## ---- Monthly means / climatology ----
  data_monthly <- dw %>%
    group_by(Zone, month_name, Month) %>%
    summarise(value_mean = mean(value, na.rm = TRUE), .groups = "drop")

  monthly_trend_tests <- dw %>%
    group_by(Zone, month_name) %>%
    summarise(
      lm_slope = coef(lm(value ~ Year))[2],
      lm_p     = summary(lm(value ~ Year))$coefficients[2, 4],
      .groups = "drop"
    )
  write_csv(monthly_trend_tests, file.path(results_dir, paste0(var, "_monthly_trend_tests.csv")))

  monthly_climatology <- ggplot(data_monthly, aes(month_name, value_mean, group = Zone, color = Zone)) +
    geom_line() +
    geom_point() +
    scale_color_manual(values = palette_zone) +
    labs(title = paste0("Monthly climatology — ", disp_name), x = NULL, y = y_lab, color = "Zone") +
    theme(axis.text.x = element_text(angle = 90))

  monthly_box <- ggplot(dw, aes(month_name, value, color = Zone)) +
    geom_boxplot(outlier.size = 0.7) +
    facet_grid(Zone ~ .) +
    scale_color_manual(values = palette_zone) +
    labs(title = paste0("Monthly means — ", disp_name), x = NULL, y = y_lab) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90))

  monthly_combo <- plot_grid(monthly_box, monthly_climatology, rel_widths = c(2, 1.3), nrow = 1)
  print(monthly_combo)
  ggsave(file.path(results_dir, paste0(var, "_monthly.png")), monthly_combo,
         width = 220, height = 120, units = "mm", dpi = 400)

  ## ---- Dynamic Factor Analysis (common trends across zones) ----
  dfa_wide <- data_annual %>%
    select(Zone, Year, value_annual) %>%
    pivot_wider(names_from = Zone, values_from = value_annual) %>%
    arrange(Year)

  dfa_mat   <- dfa_wide %>% select(-Year) %>% as.matrix() %>% t()
  # Guard against zero-variance rows (e.g. a zone with no events at all),
  # which would make z-scoring blow up to NaN/Inf.
  row_sd <- apply(dfa_mat, 1, sd, na.rm = TRUE)
  if (any(row_sd == 0 | is.na(row_sd))) {
    message("  Skipping DFA for ", var, ": one or more zones have zero variance.")
    return(invisible(NULL))
  }
  dfa_mat_z <- t(scale(t(dfa_mat)))
  rownames(dfa_mat_z) <- colnames(dfa_wide)[-1]
  colnames(dfa_mat_z) <- dfa_wide$Year

  fit_dfa <- function(n_trends) {
    MARSS(dfa_mat_z, model = list(m = n_trends, R = "diagonal and equal"),
          form = "dfa", z.score = FALSE, silent = TRUE, control = list(maxit = 2000))
  }

  dfa_models <- map(1:min(2, nrow(dfa_mat_z) - 1), fit_dfa)
  aicc_vals  <- map_dbl(dfa_models, ~ .x$AICc)
  best_dfa   <- dfa_models[[which.min(aicc_vals)]]

  print(tibble(n_trends = seq_along(aicc_vals), AICc = aicc_vals))

  trends_est <- as.data.frame(t(best_dfa$states)) %>%
    mutate(Year = dfa_wide$Year) %>%
    pivot_longer(-Year, names_to = "trend", values_to = "trend_value")

  dfa_plot <- ggplot(trends_est, aes(Year, trend_value, color = trend)) +
    geom_line(linewidth = 1) +
    labs(title = paste0(disp_name, " — MARSS/DFA common trend(s), best model m = ", best_dfa$call$model$m),
         x = NULL, y = "Trend (standardized units)")
  print(dfa_plot)
  ggsave(file.path(results_dir, paste0(var, "_dfa.png")), dfa_plot,
         width = 183, height = 120, units = "mm", dpi = 400)

  loadings <- as.data.frame(coef(best_dfa, type = "matrix")$Z)
  rownames(loadings) <- rownames(dfa_mat_z)
  write_csv(loadings %>% rownames_to_column("Zone"), file.path(results_dir, paste0(var, "_dfa_loadings.csv")))
  print(loadings)

  invisible(NULL)
}

## ---- Run for every bathymetry file, and every variable within it ----
for (data_file in data_files) {
  loaded <- load_bathymetry_file(data_file)
  data_wide      <- loaded$data_wide       # read by analyse_variable() below
  bathymetry_tag <- loaded$bathymetry_tag  # read by analyse_variable() below

  for (v in all_variables) {
    analyse_variable(v)
  }
}

## ---- Done ----
sessionInfo()