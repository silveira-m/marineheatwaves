## ============================================================
## Marine Heatwave (MHW) Time Series Analysis — Portugal (OISST, 1981-2025)
## Standalone R script (no external files needed besides the CSV below)
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

## ---- Load data ----
# >>> EDIT THIS if the CSV is not in the same folder <<<
data_file <- "/Users/miguelsilveira/Documents/GitHub/marineheatwaves/DATA/OUTPUT/MHW_summary_table_long_oisst.19812025.Portugal_800m_bathymetry.csv"
stopifnot(file.exists(data_file))

variable    <- "MHW800m"
results_dir <- paste0(variable, "_results")
if (!dir.exists(results_dir)) dir.create(results_dir)

raw <- read_csv(data_file, show_col_types = FALSE)
stopifnot(all(c("Year", "Month", "Zone", "Parameter", "Value") %in% names(raw)))

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
    number_of_events = replace_na(number_of_events, 0),
    sum_area_km2      = replace_na(sum_area_km2, 0)
  )

## ---- STL decomposition + GESD anomaly detection ----
t_window <- 133
alpha    <- 0.05

run_stl_for_zone <- function(zone) {
  df_zone <- data_wide %>% filter(Zone == zone) %>% arrange(date)
  ts_zone <- ts(df_zone$activity_degC_days_m2,
                start = c(df_zone$Year[1], df_zone$Month[1]), frequency = 12)
  
  ts_zone %>%
    stl(s.window = "periodic", t.window = t_window, robust = TRUE) %>%
    .$time.series %>%
    tsdf() %>%
    as_tibble() %>%
    anomalize(remainder, alpha = alpha, method = "gesd") %>%
    mutate(observed = trend + remainder + seasonal, Zone = zone)
}

data_stl <- map_dfr(levels(data_wide$Zone), run_stl_for_zone) %>%
  mutate(Zone = factor(Zone, levels = levels(data_wide$Zone)))

write_csv(data_stl, file.path(results_dir, paste0(variable, "_stl_anomalies.csv")))

plot_stl <- data_stl %>%
  pivot_longer(cols = c(observed, trend, seasonal, remainder),
               names_to = "component", values_to = "value") %>%
  mutate(component = factor(component, levels = c("observed", "trend", "seasonal", "remainder")))

obs_panel <- ggplot(filter(plot_stl, component == "observed"), aes(x, value, color = Zone)) +
  geom_line() +
  geom_point(data = filter(plot_stl, component == "observed", anomaly == "Yes"),
             shape = 21, size = 2.5, fill = "white") +
  facet_grid2(component ~ Zone, scales = "free_y", independent = "y") +
  scale_color_manual(values = palette_zone) +
  labs(title = "Observed MHW activity & anomalies", x = NULL, y = expression(degree*C~days~m^-2)) +
  theme(legend.position = "none")

decomp_panel <- ggplot(filter(plot_stl, component != "observed"), aes(x, value, color = Zone)) +
  geom_line() +
  geom_point(data = filter(plot_stl, component == "remainder", anomaly == "Yes"), size = 1.5) +
  facet_grid2(component ~ Zone, scales = "free_y", independent = "y") +
  scale_color_manual(values = palette_zone) +
  labs(x = NULL, y = NULL) +
  theme(legend.position = "none", strip.text.x = element_blank())

stl_combo <- plot_grid(obs_panel, decomp_panel, ncol = 1, rel_heights = c(1.2, 2))
print(stl_combo)
ggsave(file.path(results_dir, paste0(variable, "_stl.png")), stl_combo,
       width = 260, height = 220, units = "mm", dpi = 400)

## ---- Annual means, trends & breakpoints ----
data_annual <- data_wide %>%
  group_by(Zone, Year) %>%
  summarise(
    activity_total       = sum(activity_degC_days_m2, na.rm = TRUE),
    events_total          = sum(number_of_events, na.rm = TRUE),
    area_total_km2        = sum(sum_area_km2, na.rm = TRUE),
    intensity_mean_degC   = mean(mean_intensity_degC, na.rm = TRUE),
    duration_mean_days    = mean(mean_duration_days, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(across(c(intensity_mean_degC, duration_mean_days), ~ ifelse(is.nan(.x), NA, .x)))

write_csv(data_annual, file.path(results_dir, paste0(variable, "_annual_summary.csv")))

annual_trend_tests <- data_annual %>%
  group_by(Zone) %>%
  summarise(
    lm_slope  = coef(lm(activity_total ~ Year))[2],
    lm_p      = summary(lm(activity_total ~ Year))$coefficients[2, 4],
    lm_r2     = summary(lm(activity_total ~ Year))$r.squared,
    mk_tau    = Kendall::MannKendall(activity_total)$tau,
    mk_p      = Kendall::MannKendall(activity_total)$sl,
    sen_slope = trend::sens.slope(activity_total)$estimates,
    .groups = "drop"
  )
write_csv(annual_trend_tests, file.path(results_dir, paste0(variable, "_annual_trend_tests.csv")))
print(annual_trend_tests)

annual_plot <- ggplot(data_annual, aes(Year, activity_total, color = Zone, fill = Zone)) +
  geom_smooth(method = "lm", alpha = 0.15) +
  geom_line() +
  geom_point(size = 1.2) +
  facet_grid(Zone ~ ., scales = "free_y") +
  scale_color_manual(values = palette_zone) +
  scale_fill_manual(values = palette_zone) +
  labs(title = "Annual total MHW activity", x = NULL, y = expression(degree*C~days~m^-2)) +
  theme(legend.position = "none")
print(annual_plot)
ggsave(file.path(results_dir, paste0(variable, "_annual.png")), annual_plot,
       width = 150, height = 170, units = "mm", dpi = 400)

find_breaks <- function(zone) {
  df_zone <- data_annual %>% filter(Zone == zone) %>% arrange(Year)
  ts_zone <- ts(df_zone$activity_total, start = min(df_zone$Year), frequency = 1)
  bp <- tryCatch(strucchange::breakpoints(ts_zone ~ 1), error = function(e) NULL)
  
  out <- tibble(Year = df_zone$Year, value = df_zone$activity_total, Zone = zone, breakpoint = FALSE)
  if (!is.null(bp) && length(bp$breakpoints) > 0 && !anyNA(bp$breakpoints)) {
    bp_years <- df_zone$Year[bp$breakpoints]
    out$breakpoint[out$Year %in% bp_years] <- TRUE
  }
  out
}

data_bp <- map_dfr(levels(data_wide$Zone), find_breaks) %>%
  mutate(Zone = factor(Zone, levels = levels(data_wide$Zone)))

bp_plot <- ggplot(data_bp, aes(Year, value)) +
  geom_line(color = "grey60") +
  geom_vline(data = filter(data_bp, breakpoint), aes(xintercept = Year),
             color = "firebrick", linetype = "dashed") +
  geom_text(data = filter(data_bp, breakpoint), aes(x = Year, y = max(data_bp$value), label = Year),
            angle = 90, vjust = -0.4, hjust = 1, size = 3, color = "firebrick") +
  facet_grid(Zone ~ ., scales = "free_y") +
  labs(title = "Trend changes / breakpoints in annual MHW activity", x = NULL,
       y = expression(degree*C~days~m^-2))
print(bp_plot)
ggsave(file.path(results_dir, paste0(variable, "_annual_breakpoints.png")), bp_plot,
       width = 150, height = 170, units = "mm", dpi = 400)

## ---- Seasonal means ----
data_seasonal <- data_wide %>%
  group_by(Zone, season, Year) %>%
  summarise(activity_mean = mean(activity_degC_days_m2, na.rm = TRUE), .groups = "drop")

write_csv(
  data_seasonal %>% pivot_wider(names_from = Zone, values_from = activity_mean),
  file.path(results_dir, paste0(variable, "_seasonal_summary.csv"))
)

seasonal_trend_tests <- data_seasonal %>%
  group_by(Zone, season) %>%
  summarise(
    lm_slope = coef(lm(activity_mean ~ Year))[2],
    lm_p     = summary(lm(activity_mean ~ Year))$coefficients[2, 4],
    .groups = "drop"
  )
write_csv(seasonal_trend_tests, file.path(results_dir, paste0(variable, "_seasonal_trend_tests.csv")))

seasonal_box <- ggplot(data_wide, aes(x = season, y = activity_degC_days_m2, color = Zone)) +
  geom_boxplot() +
  facet_grid(Zone ~ .) +
  scale_color_manual(values = palette_zone) +
  labs(title = "Seasonal distribution of monthly MHW activity", x = NULL,
       y = expression(degree*C~days~m^-2)) +
  theme(legend.position = "none")

seasonal_trend_plot <- ggplot(data_seasonal, aes(Year, activity_mean, color = Zone, fill = Zone)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", se = TRUE) +
  facet_grid(Zone ~ season, scales = "free_y") +
  scale_color_manual(values = palette_zone) +
  scale_fill_manual(values = palette_zone) +
  labs(title = "Seasonal trends", x = NULL, y = expression(degree*C~days~m^-2)) +
  theme(legend.position = "none")

print(seasonal_box)
print(seasonal_trend_plot)
ggsave(file.path(results_dir, paste0(variable, "_seasonal_box.png")), seasonal_box,
       width = 183, height = 150, units = "mm", dpi = 400)
ggsave(file.path(results_dir, paste0(variable, "_seasonal_trend.png")), seasonal_trend_plot,
       width = 183, height = 150, units = "mm", dpi = 400)

## ---- Monthly means / climatology ----
data_monthly <- data_wide %>%
  group_by(Zone, month_name, Month) %>%
  summarise(activity_mean = mean(activity_degC_days_m2, na.rm = TRUE), .groups = "drop")

monthly_trend_tests <- data_wide %>%
  group_by(Zone, month_name) %>%
  summarise(
    lm_slope = coef(lm(activity_degC_days_m2 ~ Year))[2],
    lm_p     = summary(lm(activity_degC_days_m2 ~ Year))$coefficients[2, 4],
    .groups = "drop"
  )
write_csv(monthly_trend_tests, file.path(results_dir, paste0(variable, "_monthly_trend_tests.csv")))

monthly_climatology <- ggplot(data_monthly, aes(month_name, activity_mean, group = Zone, color = Zone)) +
  geom_line() +
  geom_point() +
  scale_color_manual(values = palette_zone) +
  labs(title = "Monthly climatology of MHW activity", x = NULL,
       y = expression(degree*C~days~m^-2), color = "Zone") +
  theme(axis.text.x = element_text(angle = 90))

monthly_box <- ggplot(data_wide, aes(month_name, activity_degC_days_m2, color = Zone)) +
  geom_boxplot(outlier.size = 0.7) +
  facet_grid(Zone ~ .) +
  scale_color_manual(values = palette_zone) +
  labs(title = "Monthly means", x = NULL, y = expression(degree*C~days~m^-2)) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90))

monthly_combo <- plot_grid(monthly_box, monthly_climatology, rel_widths = c(2, 1.3), nrow = 1)
print(monthly_combo)
ggsave(file.path(results_dir, paste0(variable, "_monthly.png")), monthly_combo,
       width = 220, height = 120, units = "mm", dpi = 400)

## ---- Dynamic Factor Analysis (common trends across zones) ----
dfa_wide <- data_annual %>%
  select(Zone, Year, activity_total) %>%
  pivot_wider(names_from = Zone, values_from = activity_total) %>%
  arrange(Year)

dfa_mat   <- dfa_wide %>% select(-Year) %>% as.matrix() %>% t()
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
  pivot_longer(-Year, names_to = "trend", values_to = "value")

dfa_plot <- ggplot(trends_est, aes(Year, value, color = trend)) +
  geom_line(linewidth = 1) +
  labs(title = paste0("MARSS/DFA common trend(s) — best model: m = ", best_dfa$call$model$m),
       x = NULL, y = "Trend (standardized units)")
print(dfa_plot)
ggsave(file.path(results_dir, paste0(variable, "_dfa.png")), dfa_plot,
       width = 183, height = 120, units = "mm", dpi = 400)

loadings <- as.data.frame(coef(best_dfa, type = "matrix")$Z)
rownames(loadings) <- rownames(dfa_mat_z)
print(loadings)

## ---- Done ----
sessionInfo()
