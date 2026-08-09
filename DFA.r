#===================================================
# DFA FOR MHW DATA USING MARSS
# Option A: choose initial filter by Zone or Bathymetry
# Supports multiple values for the filter
# Supports annual, monthly, and seasonal aggregation
#
# PATCHED VERSION:
#  - Fixed TimeIndex misalignment bug in plot_series_vs_fitted()
#  - Added degenerate-variance detection for selected models
#  - Wrapped all plotting calls in tryCatch so one bad model/parameter
#    can no longer halt the entire loop
#  - Default MARSS plot() restricted to plot.type panels that don't
#    call qqnorm() on possibly-empty/NA residual vectors
#===================================================

library(readr)
library(dplyr)
library(tidyr)
library(MARSS)
library(ggplot2)
library(cowplot) 
library(purrr)
library(stringr)
library(broom)

#---------------------------------------------------
# 1. USER SETTINGS
#---------------------------------------------------

file_path <- "/Users/miguelsilveira/Documents/GitHub/marineheatwaves/DATA/OUTPUT/allBat_Concatenated_SummaryTable_long.csv"

# Choose one:
filter_type  <- "Zone"              # "Zone" or "Bathymetry"
filter_value <- c("NW", "SW", "S")       # one or more values
# filter_value <- "NW"              # single zone example
# filter_value <- c(200, 400)       # example for Bathymetry, if filter_type = "Bathymetry"

# Available parameters in your file:
# "activity_degC_days_m2"
# "mean_area_km2"
# "mean_duration_days"
# "mean_intensity_degC"
# "number_of_events"
# "sum_area_km2"

#parameters_to_run <- NULL
#parameters_to_run <- c("activity_degC_days_m2")
# parameters_to_run <- c("mean_intensity_degC", "mean_duration_days")
#parameters_to_run <- c("activity_degC_days_m2", "mean_area_km2", "mean_duration_days",
#                       "mean_intensity_degC", "number_of_events", "sum_area_km2")

# Selected set: the four independent MHW "primitive" metrics -- frequency,
# duration, spatial extent, and intensity -- deliberately excluding
# activity_degC_days_m2 (~ intensity x duration) and sum_area_km2
# (~ mean_area_km2 x number_of_events), since those two are algebraically
# derived from metrics already in this set and would otherwise make any
# shared common trend partly circular rather than a genuine finding.
#parameters_to_run <- c("number_of_events", "mean_duration_days",
#                        "mean_area_km2", "mean_intensity_degC")

# Run DFA for ALL parameters. Section 8 (cross-parameter trend correlation)
# runs on whatever ends up in all_trends_long regardless of how many
# parameters this includes -- so this includes activity_degC_days_m2 and
# sum_area_km2 alongside their component metrics. When reviewing
# ..._cross_parameter_trend_correlations.csv, keep in mind that the pairs
# activity<->intensity, activity<->duration, sum_area<->mean_area, and
# sum_area<->number_of_events are expected to correlate strongly simply
# because of how those quantities are computed -- that's not new
# information, unlike a strong correlation between two metrics that aren't
# algebraically related to each other.
parameters_to_run <- NULL

# Time aggregation
time_agg <- "annual"   # "annual", "monthly", or "seasonal"

# Seasonal definition:
# Winter = 1,2,3 ; Spring = 4,5,6 ; Summer = 7,8,9 ; Autumn = 10,11,12

# Model settings
max_common_trends <- 3
R_structures <- c("diagonal and equal",
                  "diagonal and unequal",
                  "equalvarcov")

zscore_series <- TRUE
demean_series <- TRUE

# MARSS's plot() checks interactive() (whether this is an interactive R
# session) -- not whether the active graphics device is actually an
# on-screen device -- to decide whether to pause with "Hit <Return> to see
# next plot" between panels. Since this script is often run from an
# interactive console/RStudio, that pause fires even while writing to a
# png() file device, silently blocking execution until Return is pressed.
# Disabling it here ensures the script runs unattended regardless of how
# it's launched.
devAskNewPage(FALSE)

# Tolerance below which an estimated R or Q variance is flagged as
# "degenerate" (near-zero variance components break standardized
# residual diagnostics, e.g. qqnorm on the model's default plot).
degenerate_tol <- 1e-4

# If TRUE: when the AICc-best model is flagged degenerate, fall back to the
# best non-degenerate model instead (if one exists among the fitted
# candidates), rather than reporting/plotting an unstable fit. The AICc-best
# model is still recorded in the console/summary either way for reference.
prefer_nondegenerate <- TRUE

# Significance level used for factor-loading and common-trend confidence
# intervals (e.g. 0.05 = 95% CIs).
alpha_sig <- 0.05

results_dir <- "DFA_results"
if(!dir.exists(results_dir)) dir.create(results_dir)

#---------------------------------------------------
# 2. LOAD AND CLEAN DATA
#---------------------------------------------------

dat0 <- read_csv(file_path, guess_max = 100000)

dat0 <- dat0 %>%
  mutate(
    Zone = as.character(Zone),
    Parameter = as.character(Parameter),
    Bathymetry = as.numeric(Bathymetry),
    Year = as.integer(Year),
    Month = as.integer(Month),
    Value = as.numeric(Value)
  ) %>%
  mutate(
    Season = case_when(
      Month %in% c(1, 2, 3)    ~ "Winter",
      Month %in% c(4, 5, 6)    ~ "Spring",
      Month %in% c(7, 8, 9)    ~ "Summer",
      Month %in% c(10, 11, 12) ~ "Autumn",
      TRUE ~ NA_character_
    )
  ) %>%
  mutate(
    Season = factor(Season, levels = c("Winter", "Spring", "Summer", "Autumn"))
  )

dat0 <- dat0 %>%
  mutate(Value = ifelse(is.infinite(Value), NA, Value)) %>%
  mutate(Value = ifelse(abs(Value) > 1e10, NA, Value))

#---------------------------------------------------
# 3. APPLY INITIAL FILTER (supports multiple values)
#---------------------------------------------------

if(filter_type == "Zone"){
  
  dat1 <- dat0 %>% filter(Zone %in% filter_value)
  filter_label <- paste0("Zone_", paste(filter_value, collapse = "-"))
  
  # Each series will combine Zone + Bathymetry, so multiple zones don't get mixed
  dat1 <- dat1 %>%
    mutate(SeriesID = paste0("Zone_", Zone, "_Bath_", Bathymetry))
  
} else if(filter_type == "Bathymetry"){
  
  dat1 <- dat0 %>% filter(Bathymetry %in% as.numeric(filter_value))
  filter_label <- paste0("Bathymetry_", paste(filter_value, collapse = "-"))
  
  # Each series will combine Zone + Bathymetry, so multiple depths don't get mixed
  dat1 <- dat1 %>%
    mutate(SeriesID = paste0("Zone_", Zone, "_Bath_", Bathymetry))
  
} else {
  stop("filter_type must be 'Zone' or 'Bathymetry'")
}

if(nrow(dat1) == 0) stop("No data after filtering.")

if(is.null(parameters_to_run)){
  parameters_to_run <- sort(unique(dat1$Parameter))
}

#---------------------------------------------------
# 4. HELPER FUNCTIONS
#---------------------------------------------------

make_time_id <- function(df, time_agg = "annual"){
  if(time_agg == "annual"){
    df %>%
      mutate(Time = Year,
             TimeLabel = as.character(Year))
  } else if(time_agg == "monthly"){
    df %>%
      mutate(Time = Year + (Month - 1) / 12,
             TimeLabel = paste0(Year, "_", sprintf("%02d", Month)))
  } else if(time_agg == "seasonal"){
    season_index <- c(Winter = 1, Spring = 2, Summer = 3, Autumn = 4)
    df %>%
      mutate(
        SeasonNum = unname(season_index[as.character(Season)]),
        Time = Year + (SeasonNum - 1) / 4,
        TimeLabel = paste0(Year, "_", as.character(Season))
      )
  } else {
    stop("time_agg must be 'annual', 'monthly', or 'seasonal'")
  }
}

aggregate_time_series <- function(df, time_agg = "annual"){
  df2 <- make_time_id(df, time_agg)
  
  if(time_agg == "annual"){
    out <- df2 %>%
      group_by(Parameter, SeriesID, Year, Time, TimeLabel) %>%
      summarise(Value = mean(Value, na.rm = TRUE), .groups = "drop")
  } else if(time_agg == "monthly"){
    out <- df2 %>%
      group_by(Parameter, SeriesID, Year, Month, Time, TimeLabel) %>%
      summarise(Value = mean(Value, na.rm = TRUE), .groups = "drop")
  } else if(time_agg == "seasonal"){
    out <- df2 %>%
      group_by(Parameter, SeriesID, Year, Season, Time, TimeLabel) %>%
      summarise(Value = mean(Value, na.rm = TRUE), .groups = "drop")
  }
  
  out
}

build_marss_matrix <- function(df_param){
  time_order <- df_param %>%
    distinct(Time, TimeLabel) %>%
    arrange(Time)
  
  wide <- df_param %>%
    select(SeriesID, TimeLabel, Value) %>%
    pivot_wider(names_from = TimeLabel, values_from = Value)
  
  series_names <- wide$SeriesID
  wide$SeriesID <- NULL
  
  mat <- as.matrix(wide)
  rownames(mat) <- series_names
  storage.mode(mat) <- "numeric"
  
  ordered_cols <- time_order$TimeLabel
  mat <- mat[, ordered_cols, drop = FALSE]
  
  return(mat)
}

safe_aicc <- function(fit){
  tryCatch(MARSSaic(fit)$AICc, error = function(e) NA_real_)
}

# --- NEW: degeneracy check -------------------------------------------------
# Flags models where an estimated R or Q variance component is near zero.
# These "degenerate" fits are the usual cause of qqnorm() (and other
# standardized-residual diagnostics) failing with "y is empty or has only
# NAs" in the default MARSS plot() call, since standardized residuals are
# undefined (0/0) when the corresponding variance collapses to ~0.
#
# Also identifies WHICH series (for R) or WHICH common trend (for Q) is
# carrying the near-zero variance, using the row/col names MARSS attaches
# to the parameter matrices (these should match your SeriesID values for R,
# since that's what was used as the rownames of the input matrix `y`). This
# turns "the model is degenerate" into "series X is being fit almost
# perfectly (near-zero observation error), which is what's degenerate" --
# so you don't have to go hunting through the model_selection CSV to find
# which series/trend to investigate or drop.
check_degenerate <- function(fit, tol = 1e-4){
  out <- list(degenerate = FALSE, min_R = NA_real_, min_Q = NA_real_,
              degenerate_R_series = NA_character_,
              degenerate_Q_trend = NA_character_)
  
  tryCatch({
    par_mat <- coef(fit, type = "matrix")
    
    if(!is.null(par_mat$R)){
      Rdiag <- diag(par_mat$R)
      out$min_R <- suppressWarnings(min(Rdiag, na.rm = TRUE))
      
      if(!is.na(out$min_R) && out$min_R < tol){
        rn <- rownames(par_mat$R)
        if(is.null(rn)) rn <- paste0("series_", seq_along(Rdiag))
        out$degenerate_R_series <- paste(rn[which(Rdiag < tol)], collapse = "; ")
      }
    }
    if(!is.null(par_mat$Q)){
      Qdiag <- diag(par_mat$Q)
      out$min_Q <- suppressWarnings(min(Qdiag, na.rm = TRUE))
      
      if(!is.na(out$min_Q) && out$min_Q < tol){
        qn <- rownames(par_mat$Q)
        if(is.null(qn)) qn <- paste0("trend_", seq_along(Qdiag))
        out$degenerate_Q_trend <- paste(qn[which(Qdiag < tol)], collapse = "; ")
      }
    }
    
    out$degenerate <- (!is.na(out$min_R) && out$min_R < tol) ||
                       (!is.na(out$min_Q) && out$min_Q < tol)
  }, error = function(e) NULL)
  
  out
}
# ----------------------------------------------------------------------------

# --- NEW: AICc-based "significance" of model comparison --------------------
# MARSS model selection isn't a nested-model design in general (different m
# AND different R structures aren't all nested in one another), so a formal
# likelihood-ratio test isn't generally applicable across the whole
# candidate set. The standard practice (Burnham & Anderson 2002) is instead
# to classify support for each model *relative to the best model* using
# delta AICc thresholds. delta = 0 for the best model itself.
classify_aicc_support <- function(delta_aicc){
  cut(delta_aicc,
      breaks = c(-Inf, 2, 4, 7, 10, Inf),
      labels = c("Strong support", "Moderate support", "Weak support",
                 "Very weak support", "No support"),
      right = FALSE)
}
# ----------------------------------------------------------------------------

# --- NEW: significance of factor loadings (Z) -------------------------------
# Uses MARSSparamCIs() to get standard errors / confidence intervals for the
# estimated parameters, then extracts the Z (loadings) block. A loading is
# treated as "significant" when its CI excludes zero. Where a standard error
# is available, a two-sided z-test p-value is also computed.
# Returns NULL (rather than erroring) if MARSSparamCIs fails -- this can
# happen for models with a near-singular Hessian, which is common for
# degenerate or otherwise unstable fits.
get_loading_significance <- function(fit, alpha = 0.05){
  fit_ci <- tryCatch(MARSSparamCIs(fit, alpha = alpha), error = function(e) NULL)
  if(is.null(fit_ci)) return(NULL)
  
  Z_est <- tryCatch(coef(fit, type = "matrix")$Z, error = function(e) NULL)
  if(is.null(Z_est)) return(NULL)
  
  reshape_like <- function(x){
    if(is.null(x)) return(NULL)
    if(identical(dim(x), dim(Z_est))) return(x)
    tryCatch(matrix(x, nrow = nrow(Z_est), ncol = ncol(Z_est),
                     dimnames = dimnames(Z_est)),
             error = function(e) NULL)
  }
  
  Z_up <- reshape_like(tryCatch(coef(fit_ci, type = "matrix", what = "par.upCI")$Z,
                                 error = function(e) NULL))
  Z_lo <- reshape_like(tryCatch(coef(fit_ci, type = "matrix", what = "par.lowCI")$Z,
                                 error = function(e) NULL))
  Z_se <- reshape_like(tryCatch(coef(fit_ci, type = "matrix", what = "par.se")$Z,
                                 error = function(e) NULL))
  
  if(is.null(Z_up) || is.null(Z_lo)) return(NULL)
  
  significant <- (Z_lo > 0) | (Z_up < 0)   # CI excludes zero
  
  p_value <- NULL
  if(!is.null(Z_se)){
    z_stat <- Z_est / Z_se
    p_value <- 2 * pnorm(-abs(z_stat))
  }
  
  list(estimate = Z_est, se = Z_se, ci_lo = Z_lo, ci_up = Z_up,
       p_value = p_value, significant = significant)
}
# ----------------------------------------------------------------------------

# --- NEW: significance of common trends over time ---------------------------
# Extracts the smoother state covariance (VtT, an m x m x T array) from
# MARSSkfss() to build a pointwise confidence band around each estimated
# trend. A trend is flagged "significant" at a given time point when that
# band excludes zero at that point. Returns NULL if the smoother output
# isn't available (e.g. an unstable/degenerate fit).
get_trend_ci <- function(fit, alpha = 0.05){
  ks <- tryCatch(MARSSkfss(fit), error = function(e) NULL)
  if(is.null(ks) || is.null(ks$xtT) || is.null(ks$VtT)) return(NULL)
  
  xtT <- ks$xtT              # m x TT
  VtT <- ks$VtT               # m x m x TT
  m  <- nrow(xtT)
  TT <- ncol(xtT)
  
  z_crit <- qnorm(1 - alpha / 2)
  se <- matrix(NA_real_, nrow = m, ncol = TT)
  
  for(t in seq_len(TT)){
    Vt <- VtT[, , t]
    if(is.null(dim(Vt))) Vt <- matrix(Vt, nrow = m, ncol = m)  # m = 1 edge case
    se[, t] <- sqrt(pmax(diag(Vt), 0))
  }
  
  ci_lo <- xtT - z_crit * se
  ci_up <- xtT + z_crit * se
  significant <- (ci_lo > 0) | (ci_up < 0)   # CI excludes zero
  
  list(se = se, ci_lo = ci_lo, ci_up = ci_up, significant = significant)
}
# ----------------------------------------------------------------------------

get_fitted_trends <- function(fit){
  ks <- MARSSkfss(fit)
  ks$xtT
}

get_loadings <- function(fit){
  coef(fit, type = "matrix")$Z
}

plot_loadings <- function(Z, file_out, parameter_name, sig_stars = NULL){
  dfz <- as.data.frame(as.table(Z))
  names(dfz) <- c("Series", "Trend", "Loading")
  dfz$Stars <- if(!is.null(sig_stars)) as.vector(sig_stars) else ""
  
  p <- ggplot(dfz, aes(x = Trend, y = Series, fill = Loading)) +
    geom_tile() +
    geom_text(aes(label = Stars), color = "black", size = 4, fontface = "bold") +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b") +
    labs(title = paste("Factor loadings -", parameter_name),
         subtitle = "* p<0.05, ** p<0.01, *** p<0.001 (or * if only CI, no SE, available)",
         x = "Common trend", y = "Series") +
    theme_minimal(base_size = 11)
  
  ggsave(file_out, p, width = 8, height = 5, dpi = 300)
  p
}

plot_common_trends <- function(trends, time_labels, file_out, parameter_name,
                                ci_lo = NULL, ci_up = NULL, significant = NULL){
  trend_names <- rownames(trends)
  if(is.null(trend_names)) trend_names <- paste0("Trend", seq_len(nrow(trends)))
  
  dft <- as.data.frame(t(trends))
  names(dft) <- trend_names
  dft$TimeLabel <- time_labels
  dfl <- pivot_longer(dft, cols = -TimeLabel, names_to = "Trend", values_to = "Value")
  
  # FIX: the original version used `seq_along(TimeLabel)` directly inside
  # aes(x = ...). At that point TimeLabel is the pivoted-long column, whose
  # length is ntime * ntrend (repeated per trend), not ntime -- so
  # seq_along() returned 1:(ntime*ntrend) instead of a repeating 1:ntime
  # index. pivot_longer() here iterates TimeLabel-major, Trend-minor (same
  # ordering reasoning as the fix in plot_series_vs_fitted), so the correct
  # index is built explicitly below instead of relying on seq_along(TimeLabel).
  dfl$TimeIndex <- rep(seq_along(time_labels), each = length(trend_names))
  
  have_ci <- !is.null(ci_lo) && !is.null(ci_up)
  if(have_ci){
    lo_df <- as.data.frame(t(ci_lo)); names(lo_df) <- trend_names
    lo_df$TimeLabel <- time_labels
    lo_long <- pivot_longer(lo_df, cols = -TimeLabel, names_to = "Trend", values_to = "CI_low")
    
    up_df <- as.data.frame(t(ci_up)); names(up_df) <- trend_names
    up_df$TimeLabel <- time_labels
    up_long <- pivot_longer(up_df, cols = -TimeLabel, names_to = "Trend", values_to = "CI_high")
    
    dfl <- dfl %>%
      left_join(lo_long, by = c("TimeLabel", "Trend")) %>%
      left_join(up_long, by = c("TimeLabel", "Trend"))
  }
  
  have_sig <- !is.null(significant)
  if(have_sig){
    sig_df <- as.data.frame(t(significant)); names(sig_df) <- trend_names
    sig_df$TimeLabel <- time_labels
    sig_long <- pivot_longer(sig_df, cols = -TimeLabel, names_to = "Trend", values_to = "Significant")
    dfl <- dfl %>% left_join(sig_long, by = c("TimeLabel", "Trend"))
  }
  
  p <- ggplot(dfl, aes(x = TimeIndex, y = Value, color = Trend, group = Trend))
  
  if(have_ci){
    p <- p + geom_ribbon(aes(ymin = CI_low, ymax = CI_high, fill = Trend),
                          alpha = 0.15, color = NA)
  }
  
  p <- p +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    geom_line(linewidth = 0.8)
  
  # Explicit markers on top of the line at every time point where the trend's
  # CI excludes zero, in addition to the shaded ribbon -- makes significant
  # points identifiable even where the ribbon is hard to read (e.g. printed
  # in black & white, or where trends overlap).
  if(have_sig){
    # NOTE: filter on the logical column directly (not isTRUE(), which is
    # not vectorized -- see the isTRUE() bug fixed elsewhere in this script).
    # `%in% TRUE` also safely drops any NA rows instead of erroring on them.
    p <- p + geom_point(data = dfl %>% filter(Significant %in% TRUE),
                         size = 2.2, shape = 21, fill = "white", stroke = 1.1)
  }
  
  p <- p +
    scale_x_continuous(breaks = seq_along(time_labels), labels = time_labels) +
    labs(title = paste("Estimated common trends -", parameter_name),
         subtitle = if(have_ci) "Shaded band = 95% CI; open circles mark points where the trend is significantly non-zero" else NULL,
         x = "Time", y = "Trend value") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  
  ggsave(file_out, p, width = 10, height = 5, dpi = 300)
  p
}

plot_series_vs_fitted <- function(y, fit, time_labels, file_out, parameter_name){
  ks <- MARSSkfss(fit)
  yhat <- ks$ytT
  
  obs <- as.data.frame(t(y))
  fitd <- as.data.frame(t(yhat))
  
  obs$TimeLabel <- time_labels
  fitd$TimeLabel <- time_labels
  
  obs_l <- pivot_longer(obs, cols = -TimeLabel, names_to = "Series", values_to = "Observed")
  fit_l <- pivot_longer(fitd, cols = -TimeLabel, names_to = "Series", values_to = "Fitted")
  
  dd <- left_join(obs_l, fit_l, by = c("TimeLabel", "Series"))
  
  # FIX: pivot_longer() on `obs`/`fitd` iterates row-by-row (TimeLabel is the
  # outer loop) and expands each row across the Series columns (inner loop).
  # So dd is ordered TimeLabel-major, Series-minor: (t1,s1),(t1,s2)...(t1,sK),
  # (t2,s1)... The original `rep(seq_along(time_labels), times = nrow(y))`
  # produced a Series-major sequence instead, which mismatched every row's
  # true time index. `each =` reproduces the correct TimeLabel-major order.
  dd$TimeIndex <- rep(seq_along(time_labels), each = nrow(y))
  
  p <- ggplot(dd, aes(TimeIndex)) +
    geom_line(aes(y = Observed), color = "grey60") +
    geom_line(aes(y = Fitted), color = "#d95f02") +
    facet_wrap(~Series, scales = "free_y") +
    scale_x_continuous(breaks = seq_along(time_labels), labels = time_labels) +
    labs(title = paste("Observed vs fitted -", parameter_name),
         x = "Time", y = "Value") +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  
  ggsave(file_out, p, width = 12, height = 7, dpi = 300)
  p
}

plot_canonical_correlations <- function(y, trends, file_out, parameter_name, alpha = 0.05){
  cors <- cor(t(y), t(trends), use = "pairwise.complete.obs")
  
  n_series <- nrow(y)
  n_trends <- nrow(trends)
  pvals <- matrix(NA_real_, nrow = n_series, ncol = n_trends, dimnames = dimnames(cors))
  
  # Standard t-test for a Pearson correlation, using the pairwise sample
  # size for each (series, trend) pair specifically -- since y can contain
  # NAs, different pairs may be based on different numbers of usable time
  # points, so n isn't a single shared value across the whole matrix.
  for(i in seq_len(n_series)){
    for(j in seq_len(n_trends)){
      ok <- !is.na(y[i, ]) & !is.na(trends[j, ])
      n_pair <- sum(ok)
      r <- cors[i, j]
      if(n_pair > 2 && !is.na(r) && abs(r) < 1){
        tstat <- r * sqrt(n_pair - 2) / sqrt(1 - r^2)
        pvals[i, j] <- 2 * pt(-abs(tstat), df = n_pair - 2)
      }
    }
  }
  
  stars <- matrix("", nrow = n_series, ncol = n_trends, dimnames = dimnames(cors))
  stars[!is.na(pvals) & pvals < 0.001] <- "***"
  stars[!is.na(pvals) & pvals >= 0.001 & pvals < 0.01]  <- "**"
  stars[!is.na(pvals) & pvals >= 0.01  & pvals < 0.05]  <- "*"
  
  dfc <- as.data.frame(as.table(cors))
  names(dfc) <- c("Series", "Trend", "Correlation")
  dfc$Label <- paste0(round(dfc$Correlation, 2), as.vector(stars))
  
  p <- ggplot(dfc, aes(x = Trend, y = Series, fill = Correlation)) +
    geom_tile() +
    geom_text(aes(label = Label), size = 3) +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", limits = c(-1,1)) +
    labs(title = paste("Correlations: time series vs common trends -", parameter_name),
         subtitle = "* p<0.05, ** p<0.01, *** p<0.001",
         x = "Common trend", y = "Time series") +
    theme_minimal(base_size = 11)
  
  ggsave(file_out, p, width = 8, height = 5, dpi = 300)
  
  out <- as.data.frame(as.table(cors))
  names(out) <- c("Series", "Trend", "Correlation")
  out$p_value <- as.vector(pvals)
  out$significant <- !is.na(out$p_value) & out$p_value < alpha
  out
}

# --- NEW: visualize AICc-based model support --------------------------------
# Bar chart of delta AICc per candidate model, colored by the Burnham &
# Anderson support category from classify_aicc_support(). Puts the
# significance-of-model-comparison information (currently only in the
# console message and model_selection CSV) into a plot as well.
plot_model_comparison <- function(model_sel, file_out, parameter_name){
  support_levels <- c("Strong support", "Moderate support", "Weak support",
                       "Very weak support", "No support")
  
  df <- model_sel %>%
    mutate(model_label = paste0("m=", m, ", ", R_structure),
           model_label = factor(model_label, levels = rev(model_label[order(delta_AICc)])),
           AICc_support = factor(AICc_support, levels = support_levels))
  
  p <- ggplot(df, aes(x = model_label, y = delta_AICc, fill = AICc_support)) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(values = c("Strong support" = "#1a9850",
                                  "Moderate support" = "#91cf60",
                                  "Weak support" = "#fee08b",
                                  "Very weak support" = "#fc8d59",
                                  "No support" = "#d73027"),
                       drop = FALSE) +
    labs(title = paste("Model comparison (delta AICc) -", parameter_name),
         subtitle = "Lower delta AICc = better supported, relative to the AICc-best model",
         x = NULL, y = "Delta AICc", fill = "Support") +
    theme_minimal(base_size = 11)
  
  ggsave(file_out, p, width = 8, height = 5, dpi = 300)
  p
}

fit_dfa_models <- function(y, max_common_trends, R_structures,
                           zscore_series = TRUE, demean_series = TRUE,
                           degenerate_tol = 1e-4){
  
  n_ts <- nrow(y)
  max_m <- min(max_common_trends, n_ts - 1)
  if(max_m < 1) stop("Need at least 2 series for DFA.")
  if(max_m < max_common_trends){
    cat("Note: requested max_common_trends =", max_common_trends,
        "reduced to", max_m, "given", n_ts, "usable series.\n")
  }
  
  model_table <- list()
  fit_list <- list()
  k <- 1
  
  for(m in 1:max_m){
    for(Rstr in R_structures){
      model_list <- list(m = m, R = Rstr)
      
      fit_try <- tryCatch(
        MARSS(
          y,
          model = model_list,
          form = "dfa",
          z.score = zscore_series,
          demean = demean_series,
          method = "kem",
          control = list(maxit = 5000, allow.degen = TRUE),
          silent = TRUE
        ),
        error = function(e) NULL
      )
      
      if(!is.null(fit_try)){
        aicc <- safe_aicc(fit_try)
        logLik_val <- tryCatch(logLik(fit_try), error = function(e) NA_real_)
        converged <- tryCatch(fit_try$convergence == 0, error = function(e) NA)
        degen_check <- check_degenerate(fit_try, tol = degenerate_tol)
        
        model_table[[k]] <- data.frame(
          m = m,
          R_structure = Rstr,
          AICc = aicc,
          logLik = as.numeric(logLik_val),
          converged = converged,
          degenerate = degen_check$degenerate,
          min_R_variance = degen_check$min_R,
          min_Q_variance = degen_check$min_Q,
          degenerate_R_series = degen_check$degenerate_R_series,
          degenerate_Q_trend = degen_check$degenerate_Q_trend,
          stringsAsFactors = FALSE
        )
        
        fit_list[[paste0("m", m, "_", make.names(Rstr))]] <- fit_try
        k <- k + 1
      }
    }
  }
  
  model_df <- bind_rows(model_table) %>% arrange(AICc)
  list(summary = model_df, fits = fit_list)
}

#---------------------------------------------------
# 5. AGGREGATE
#---------------------------------------------------

dat2 <- aggregate_time_series(dat1, time_agg = time_agg)

#---------------------------------------------------
# 6. RUN DFA PARAMETER BY PARAMETER
#---------------------------------------------------

all_model_summaries <- list()

# Stores each parameter's extracted common trends (long format, one column
# per trend named "<parameter>_<TrendName>", plus TimeLabel), so that after
# the main loop finishes, section 8 can correlate trends ACROSS parameters
# without fitting one large joint DFA over all parameters at once (see
# conversation notes on why that's risky: definitional relationships between
# MHW metrics, distributional mismatch, and dimensionality/degeneracy).
all_trends_long <- list()

for(param_i in parameters_to_run){
  
  cat("\n==============================\n")
  cat("Running DFA for:", param_i, "\n")
  cat("==============================\n")
  
  dp <- dat2 %>% filter(Parameter == param_i)
  
  y <- build_marss_matrix(dp)
  time_labels <- colnames(y)
  
  keep_rows <- apply(y, 1, function(x) sum(!is.na(x)) >= 5)
  y <- y[keep_rows, , drop = FALSE]
  
  if(nrow(y) < 2){
    cat("Skipping", param_i, "- fewer than 2 usable series.\n")
    next
  }
  
  param_clean <- str_replace_all(param_i, "[^A-Za-z0-9_]", "_")
  out_dir_param <- file.path(results_dir, paste0(filter_label, "_", time_agg, "_", param_clean))
  if(!dir.exists(out_dir_param)) dir.create(out_dir_param, recursive = TRUE)
  
  write.csv(y, file.path(out_dir_param, paste0(param_clean, "_matrix_used.csv")))
  
  fits_out <- fit_dfa_models(
    y = y,
    max_common_trends = max_common_trends,
    R_structures = R_structures,
    zscore_series = zscore_series,
    demean_series = demean_series,
    degenerate_tol = degenerate_tol
  )
  
  model_sel <- fits_out$summary
  
  if(nrow(model_sel) > 0){
    # delta AICc and Burnham & Anderson support classification, relative to
    # the AICc-best model in this candidate set (delta = 0 for that model).
    model_sel <- model_sel %>%
      mutate(delta_AICc = AICc - min(AICc, na.rm = TRUE),
             AICc_support = as.character(classify_aicc_support(delta_AICc)))
    
    if(nrow(model_sel) >= 2){
      runner_up <- model_sel %>% slice(2)
      cat("Model comparison:", param_i, "- best: m =", model_sel$m[1], ",",
          model_sel$R_structure[1], "(AICc =", round(model_sel$AICc[1], 2), ")",
          "vs runner-up: m =", runner_up$m, ",", runner_up$R_structure,
          "(AICc =", round(runner_up$AICc, 2), ", delta AICc =",
          round(runner_up$delta_AICc, 2), "->", runner_up$AICc_support, ")\n")
    }
    # NOTE: delta_AICc / AICc_support are still computed and written to
    # ..._model_selection.csv above -- only the bar-chart plot of them was
    # removed, per request to keep plot-level significance markers limited
    # to the loadings heatmap and the canonical correlations plot.
  }
  
  all_model_summaries[[param_i]] <- model_sel
  write.csv(model_sel,
            file.path(out_dir_param, paste0(param_clean, "_model_selection.csv")),
            row.names = FALSE)
  
  if(nrow(model_sel) == 0){
    cat("No valid models for", param_i, "\n")
    next
  }
  
  best <- model_sel %>% slice(1)
  
  if(isTRUE(best$degenerate)){
    cat("WARNING:", param_i, "- AICc-best model (m =", best$m, ",", best$R_structure,
        ") has a near-zero variance estimate (min R =", round(best$min_R_variance, 6),
        ", min Q =", round(best$min_Q_variance, 6), ").\n")
    if(!is.na(best$degenerate_R_series)){
      cat("  -> Series with near-zero OBSERVATION variance (R):",
          best$degenerate_R_series,
          "-- this series is being fit almost perfectly by the trend(s).",
          "Consider inspecting, aggregating, or dropping it.\n")
    }
    if(!is.na(best$degenerate_Q_trend)){
      cat("  -> Common trend with near-zero PROCESS variance (Q):",
          best$degenerate_Q_trend,
          "-- this trend is essentially deterministic/flat.",
          "Consider reducing max_common_trends.\n")
    }
    
    if(isTRUE(prefer_nondegenerate)){
      # NOTE: isTRUE() is not vectorized -- it only evaluates a single
      # logical value. Using filter(!isTRUE(degenerate)) here was a bug:
      # applied to the whole `degenerate` column it silently returned a
      # single scalar FALSE for the entire vector, making `!isTRUE(...)`
      # a no-op TRUE that dplyr recycled across every row -- so nondegen
      # was actually identical to model_sel (unfiltered), and slice(1)
      # kept re-selecting the same degenerate AICc-best row. `!degenerate`
      # is the correct, properly vectorized row-wise filter.
      nondegen <- model_sel %>% filter(!is.na(degenerate) & !degenerate)
      
      if(nrow(nondegen) > 0){
        fallback <- nondegen %>% slice(1)
        delta_aicc <- fallback$AICc - best$AICc
        cat("  -> Falling back to best NON-degenerate model instead: m =", fallback$m,
            ",", fallback$R_structure, "(AICc =", round(fallback$AICc, 2),
            ", delta AICc vs. degenerate best =", round(delta_aicc, 2), ")\n")
        best <- fallback
      } else {
        cat("  -> No non-degenerate model available for", param_i,
            "- proceeding with the degenerate AICc-best model.",
            "Diagnostic plots relying on standardized residuals may be skipped or unreliable.\n")
      }
    } else {
      cat("  -> prefer_nondegenerate is FALSE; proceeding with the degenerate model as-is.\n")
    }
  }
  
  best_name <- paste0("m", best$m, "_", make.names(best$R_structure))
  best_fit <- fits_out$fits[[best_name]]
  
  saveRDS(best_fit, file.path(out_dir_param, paste0(param_clean, "_best_dfa_model.rds")))
  
  sink(file.path(out_dir_param, paste0(param_clean, "_best_model_summary.txt")))
  cat("Parameter:", param_i, "\n")
  cat("Filter:", filter_label, "\n")
  cat("Time aggregation:", time_agg, "\n")
  cat("Selected model (m =", best$m, ",", best$R_structure, ") - degenerate:",
      isTRUE(best$degenerate), "\n")
  if(isTRUE(best$degenerate)){
    if(!is.na(best$degenerate_R_series)){
      cat("Series with near-zero observation variance (R):", best$degenerate_R_series, "\n")
    }
    if(!is.na(best$degenerate_Q_trend)){
      cat("Trend with near-zero process variance (Q):", best$degenerate_Q_trend, "\n")
    }
  }
  cat("\n")
  print(best)
  cat("\nMODEL OBJECT:\n")
  print(best_fit)
  cat("\nTIDY:\n")
  print(tryCatch(broom::tidy(best_fit), error = function(e) "broom::tidy not available"))
  cat("\nGLANCE:\n")
  print(tryCatch(broom::glance(best_fit), error = function(e) "broom::glance not available"))
  sink()
  
  # --- Default MARSS diagnostic plots ---------------------------------------
  # Each plot.type is rendered to its OWN file in its own tryCatch, rather
  # than one multi-panel plot() call to a single png(). Two reasons:
  #  1) A single png() device only ever keeps the LAST page drawn to it --
  #     multiple plot.type panels sent to one file would silently overwrite
  #     each other, so one file per panel is required to keep them all.
  #  2) A degenerate model (near-zero R or Q variance) can make a specific
  #     panel's internal Cholesky-based CI computation fail (e.g.
  #     "state.resids.xtT", or anything prefixed "std."/"qqplot.") without
  #     affecting the other panels. Isolating each panel in its own
  #     tryCatch means one failing panel no longer costs you the panels
  #     that would otherwise have rendered fine.
  default_plot_types <- c("fitted.ytT", "xtT", "model.resids.ytT", "state.resids.xtT")
  
  for(pt in default_plot_types){
    tryCatch({
      png(file.path(out_dir_param, paste0(param_clean, "_MARSS_", pt, ".png")),
          width = 1600, height = 1200, res = 200)
      plot(best_fit, plot.type = pt)
      dev.off()
    }, error = function(e){
      # Defensive cleanup: close any graphics device left open by a plot()
      # call that errored out partway through, so a stray open device
      # can't corrupt the next panel's (or next parameter's) png() output.
      while(!is.null(dev.list())) dev.off()
      cat("Skipping MARSS plot type '", pt, "' for ", param_i,
          " - error: ", conditionMessage(e), "\n", sep = "")
    })
  }
  
  trends <- get_fitted_trends(best_fit)
  Z <- get_loadings(best_fit)
  
  write.csv(Z, file.path(out_dir_param, paste0(param_clean, "_loadings_Z.csv")))
  write.csv(trends, file.path(out_dir_param, paste0(param_clean, "_common_trends.csv")))
  
  # --- Collect trends for cross-parameter correlation analysis (section 8) ---
  trend_names_out <- rownames(trends)
  if(is.null(trend_names_out)) trend_names_out <- paste0("Trend", seq_len(nrow(trends)))
  
  trend_export <- as.data.frame(t(trends))
  names(trend_export) <- paste0(param_clean, "_", trend_names_out)
  trend_export$TimeLabel <- time_labels
  
  all_trends_long[[param_i]] <- trend_export
  
  # --- Significance: factor loadings ----------------------------------------
  loading_sig <- tryCatch(get_loading_significance(best_fit, alpha = alpha_sig),
                           error = function(e) NULL)
  sig_stars <- NULL
  
  if(!is.null(loading_sig)){
    sig_stars <- matrix("", nrow = nrow(Z), ncol = ncol(Z), dimnames = dimnames(Z))
    if(!is.null(loading_sig$p_value)){
      sig_stars[loading_sig$p_value < 0.001] <- "***"
      sig_stars[loading_sig$p_value >= 0.001 & loading_sig$p_value < 0.01]  <- "**"
      sig_stars[loading_sig$p_value >= 0.01  & loading_sig$p_value < 0.05]  <- "*"
    } else {
      sig_stars[loading_sig$significant] <- "*"
    }
    
    loading_sig_df <- as.data.frame(as.table(Z))
    names(loading_sig_df) <- c("Series", "Trend", "Estimate")
    loading_sig_df$CI_low  <- as.vector(loading_sig$ci_lo)
    loading_sig_df$CI_high <- as.vector(loading_sig$ci_up)
    if(!is.null(loading_sig$p_value)) loading_sig_df$p_value <- as.vector(loading_sig$p_value)
    loading_sig_df$significant <- as.vector(loading_sig$significant)
    loading_sig_df$stars <- as.vector(sig_stars)
    
    write.csv(loading_sig_df,
              file.path(out_dir_param, paste0(param_clean, "_loadings_significance.csv")),
              row.names = FALSE)
  } else {
    cat("Could not compute loading confidence intervals for", param_i,
        "(MARSSparamCIs failed or returned an unexpected structure) -",
        "skipping loadings significance table.\n")
  }
  
  # --- Significance: common trends over time --------------------------------
  trend_ci <- tryCatch(get_trend_ci(best_fit, alpha = alpha_sig), error = function(e) NULL)
  
  if(!is.null(trend_ci)){
    trend_names <- rownames(trends)
    if(is.null(trend_names)) trend_names <- paste0("Trend", seq_len(nrow(trends)))
    
    to_long <- function(mat, value_name){
      df <- as.data.frame(t(mat)); names(df) <- trend_names
      df$TimeLabel <- time_labels
      pivot_longer(df, cols = -TimeLabel, names_to = "Trend", values_to = value_name)
    }
    
    trend_sig_out <- to_long(trends, "Estimate") %>%
      left_join(to_long(trend_ci$se, "SE"), by = c("TimeLabel", "Trend")) %>%
      left_join(to_long(trend_ci$ci_lo, "CI_low"), by = c("TimeLabel", "Trend")) %>%
      left_join(to_long(trend_ci$ci_up, "CI_high"), by = c("TimeLabel", "Trend")) %>%
      left_join(to_long(trend_ci$significant, "significant"), by = c("TimeLabel", "Trend"))
    
    write.csv(trend_sig_out,
              file.path(out_dir_param, paste0(param_clean, "_common_trends_significance.csv")),
              row.names = FALSE)
  } else {
    cat("Could not compute trend confidence intervals for", param_i,
        "- skipping trend significance table.\n")
  }
  
  tryCatch({
    plot_loadings(
      Z = Z,
      file_out = file.path(out_dir_param, paste0(param_clean, "_loadings_heatmap.png")),
      parameter_name = param_i,
      sig_stars = sig_stars
    )
  }, error = function(e){
    cat("Skipping loadings plot for", param_i, "- error:", conditionMessage(e), "\n")
  })
  
  tryCatch({
    plot_common_trends(
      trends = trends,
      time_labels = time_labels,
      file_out = file.path(out_dir_param, paste0(param_clean, "_common_trends_plot.png")),
      parameter_name = param_i
      # NOTE: ci_lo / ci_up / significant intentionally not passed here --
      # per request, significance markers are shown only on the loadings
      # heatmap and canonical correlations plot, not on this one. The
      # underlying numbers are still written to
      # ..._common_trends_significance.csv above if you want them.
    )
  }, error = function(e){
    cat("Skipping common trends plot for", param_i, "- error:", conditionMessage(e), "\n")
  })
  
  tryCatch({
    plot_series_vs_fitted(
      y = y,
      fit = best_fit,
      time_labels = time_labels,
      file_out = file.path(out_dir_param, paste0(param_clean, "_observed_vs_fitted.png")),
      parameter_name = param_i
    )
  }, error = function(e){
    cat("Skipping observed-vs-fitted plot for", param_i, "- error:", conditionMessage(e), "\n")
  })
  
  cors <- NULL
  tryCatch({
    cors <- plot_canonical_correlations(
      y = y,
      trends = trends,
      file_out = file.path(out_dir_param, paste0(param_clean, "_series_trends_correlations.png")),
      parameter_name = param_i,
      alpha = alpha_sig
    )
  }, error = function(e){
    cat("Skipping canonical correlations plot for", param_i, "- error:", conditionMessage(e), "\n")
  })
  
  if(!is.null(cors)){
    write.csv(cors,
              file.path(out_dir_param, paste0(param_clean, "_series_trends_correlations.csv")),
              row.names = FALSE)
  }
}

#---------------------------------------------------
# 7. SAVE GLOBAL MODEL SELECTION TABLE
#---------------------------------------------------

if(length(all_model_summaries) > 0){
  global_model_table <- bind_rows(all_model_summaries, .id = "Parameter")
  write.csv(global_model_table,
            file.path(results_dir, paste0(filter_label, "_", time_agg, "_all_parameters_model_selection.csv")),
            row.names = FALSE)
}

cat("\nAll DFA analyses completed.\n")

#---------------------------------------------------
# 8. CROSS-PARAMETER TREND CORRELATIONS
#---------------------------------------------------
# Correlates the extracted common trends AGAINST EACH OTHER, across
# parameters, instead of fitting one large joint DFA over all parameters at
# once. This answers "do the common trends found for different metrics
# (activity, duration, area, etc.) move together over time?" while avoiding
# the risks of a single joint model: (a) some MHW metrics are definitionally
# related (e.g. activity ~ intensity x duration), so a joint DFA could just
# recover that identity rather than a genuine shared driver; (b) combining
# ~15 series per parameter across 6 parameters would push series count into
# a range where degenerate/non-converging fits (already seen with far fewer
# series) become far more likely; (c) different parameter types (counts vs.
# continuous, different skew) don't share one Gaussian state-space model as
# comfortably as series of the same type do.
#
# DFA above still runs for ALL parameters in parameters_to_run. This section
# is restricted separately to the four independent "primitive" metrics --
# frequency, duration, spatial extent, and intensity -- excluding
# activity_degC_days_m2 (~ intensity x duration) and sum_area_km2
# (~ mean_area_km2 x number_of_events), since those two are algebraically
# derived from metrics already in this set and would make any correlation
# here partly circular rather than a genuine finding. Change
# cross_parameter_include below if you want a different subset (or set it
# to NULL to include every parameter that was run).
cross_parameter_include <- c("number_of_events", "mean_duration_days",
                              "mean_area_km2", "mean_intensity_degC")

trends_for_cross <- if(is.null(cross_parameter_include)){
  all_trends_long
} else {
  all_trends_long[names(all_trends_long) %in% cross_parameter_include]
}

if(length(trends_for_cross) >= 2){
  
  # Merge all parameters' trends on TimeLabel. Different parameters can, in
  # principle, have slightly different available time points; unmatched
  # points become NA after the join and are simply excluded pairwise when
  # computing each correlation (same approach used elsewhere in this script).
  trend_wide <- Reduce(function(a, b) full_join(a, b, by = "TimeLabel"), trends_for_cross)
  
  trend_cols <- setdiff(names(trend_wide), "TimeLabel")
  trend_mat <- t(as.matrix(trend_wide[, trend_cols, drop = FALSE]))
  colnames(trend_mat) <- trend_wide$TimeLabel
  storage.mode(trend_mat) <- "numeric"
  
  n_trends_total <- nrow(trend_mat)
  cross_cors  <- matrix(NA_real_, n_trends_total, n_trends_total,
                         dimnames = list(trend_cols, trend_cols))
  cross_pvals <- cross_cors
  
  for(i in seq_len(n_trends_total)){
    for(j in seq_len(n_trends_total)){
      xi <- trend_mat[i, ]; xj <- trend_mat[j, ]
      ok <- !is.na(xi) & !is.na(xj)
      n_pair <- sum(ok)
      if(n_pair > 2){
        r <- suppressWarnings(cor(xi, xj, use = "pairwise.complete.obs"))
        cross_cors[i, j] <- r
        if(!is.na(r) && abs(r) < 1){
          tstat <- r * sqrt(n_pair - 2) / sqrt(1 - r^2)
          cross_pvals[i, j] <- 2 * pt(-abs(tstat), df = n_pair - 2)
        }
      }
    }
  }
  
  # Long-format table: one row per trend pair, with the source parameter
  # split out of each trend's name and a same_parameter flag, so
  # within-parameter pairs (a normal by-product of DFA -- trends within one
  # model are fit to be close to orthogonal) can be filtered out from the
  # genuinely cross-parameter relationships you actually asked about.
  cross_df <- as.data.frame(as.table(cross_cors))
  names(cross_df) <- c("Trend_A", "Trend_B", "Correlation")
  cross_df$p_value <- as.vector(cross_pvals)
  cross_df$significant <- !is.na(cross_df$p_value) & cross_df$p_value < alpha_sig
  
  extract_param <- function(trend_col){
    sub("_(X?[0-9]+|Trend[0-9]+)$", "", trend_col)
  }
  cross_df$Parameter_A <- extract_param(as.character(cross_df$Trend_A))
  cross_df$Parameter_B <- extract_param(as.character(cross_df$Trend_B))
  cross_df$same_parameter <- cross_df$Parameter_A == cross_df$Parameter_B
  
  write.csv(cross_df,
            file.path(results_dir, paste0(filter_label, "_", time_agg,
                                           "_cross_parameter_trend_correlations.csv")),
            row.names = FALSE)
  
  tryCatch({
    stars_mat <- matrix("", n_trends_total, n_trends_total, dimnames = dimnames(cross_cors))
    stars_mat[!is.na(cross_pvals) & cross_pvals < 0.001] <- "***"
    stars_mat[!is.na(cross_pvals) & cross_pvals >= 0.001 & cross_pvals < 0.01]  <- "**"
    stars_mat[!is.na(cross_pvals) & cross_pvals >= 0.01  & cross_pvals < 0.05]  <- "*"
    
    plot_df <- as.data.frame(as.table(cross_cors))
    names(plot_df) <- c("Trend_A", "Trend_B", "Correlation")
    plot_df$Label <- paste0(round(plot_df$Correlation, 2), as.vector(stars_mat))
    
    p <- ggplot(plot_df, aes(x = Trend_A, y = Trend_B, fill = Correlation)) +
      geom_tile() +
      geom_text(aes(label = Label), size = 2.6) +
      scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", limits = c(-1, 1)) +
      labs(title = paste("Cross-parameter common trend correlations -", filter_label, time_agg),
           subtitle = "* p<0.05, ** p<0.01, *** p<0.001. Within-parameter trend pairs included for reference.",
           x = NULL, y = NULL) +
      theme_minimal(base_size = 10) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    
    ggsave(file.path(results_dir, paste0(filter_label, "_", time_agg,
                                          "_cross_parameter_trend_correlations.png")),
           p, width = max(8, n_trends_total * 0.6), height = max(6, n_trends_total * 0.5), dpi = 300)
  }, error = function(e){
    cat("Skipping cross-parameter trend correlation plot - error:", conditionMessage(e), "\n")
  })
  
  top_cross <- cross_df %>%
    filter(!same_parameter, significant, as.character(Trend_A) < as.character(Trend_B)) %>%
    arrange(desc(abs(Correlation))) %>%
    head(10)
  
  if(nrow(top_cross) > 0){
    cat("\nStrongest significant CROSS-PARAMETER trend correlations:\n")
    print(top_cross[, c("Trend_A", "Trend_B", "Correlation", "p_value")])
  } else {
    cat("\nNo significant cross-parameter trend correlations found (alpha =", alpha_sig, ").\n")
  }
  
} else {
  cat("\nCross-parameter trend correlation skipped - fewer than 2 parameters produced a valid DFA model.\n")
}