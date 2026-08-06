#===================================================
# DFA FOR MHW DATA USING MARSS
# Option A: choose initial filter by Zone or Bathymetry
# Supports multiple values for the filter
# Supports annual, monthly, and seasonal aggregation
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
filter_value <- c("NW", "SW")       # one or more values
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
parameters_to_run <- c("activity_degC_days_m2")
# parameters_to_run <- c("mean_intensity_degC", "mean_duration_days")
#parameters_to_run <- c("activity_degC_days_m2", "mean_area_km2", "mean_duration_days",
#                       "mean_intensity_degC", "number_of_events", "sum_area_km2")

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

get_fitted_trends <- function(fit){
  ks <- MARSSkfss(fit)
  ks$xtT
}

get_loadings <- function(fit){
  coef(fit, type = "matrix")$Z
}

plot_loadings <- function(Z, file_out, parameter_name){
  dfz <- as.data.frame(as.table(Z))
  names(dfz) <- c("Series", "Trend", "Loading")
  
  p <- ggplot(dfz, aes(x = Trend, y = Series, fill = Loading)) +
    geom_tile() +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b") +
    labs(title = paste("Factor loadings -", parameter_name),
         x = "Common trend", y = "Series") +
    theme_minimal(base_size = 11)
  
  ggsave(file_out, p, width = 8, height = 5, dpi = 300)
  p
}

plot_common_trends <- function(trends, time_labels, file_out, parameter_name){
  dft <- as.data.frame(t(trends))
  dft$TimeLabel <- time_labels
  dfl <- pivot_longer(dft, cols = -TimeLabel, names_to = "Trend", values_to = "Value")
  
  p <- ggplot(dfl, aes(x = seq_along(TimeLabel), y = Value, color = Trend, group = Trend)) +
    geom_line(linewidth = 0.8) +
    scale_x_continuous(breaks = seq_along(time_labels), labels = time_labels) +
    labs(title = paste("Estimated common trends -", parameter_name),
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
  dd$TimeIndex <- rep(seq_along(time_labels), times = nrow(y))
  
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

plot_canonical_correlations <- function(y, trends, file_out, parameter_name){
  cors <- cor(t(y), t(trends), use = "pairwise.complete.obs")
  
  dfc <- as.data.frame(as.table(cors))
  names(dfc) <- c("Series", "Trend", "Correlation")
  
  p <- ggplot(dfc, aes(x = Trend, y = Series, fill = Correlation)) +
    geom_tile() +
    geom_text(aes(label = round(Correlation, 2)), size = 3) +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", limits = c(-1,1)) +
    labs(title = paste("Correlations: time series vs common trends -", parameter_name),
         x = "Common trend", y = "Time series") +
    theme_minimal(base_size = 11)
  
  ggsave(file_out, p, width = 8, height = 5, dpi = 300)
  cors
}

fit_dfa_models <- function(y, max_common_trends, R_structures,
                           zscore_series = TRUE, demean_series = TRUE){
  
  n_ts <- nrow(y)
  max_m <- min(max_common_trends, n_ts - 1)
  if(max_m < 1) stop("Need at least 2 series for DFA.")
  
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
        
        model_table[[k]] <- data.frame(
          m = m,
          R_structure = Rstr,
          AICc = aicc,
          logLik = as.numeric(logLik_val),
          converged = converged,
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
    demean_series = demean_series
  )
  
  model_sel <- fits_out$summary
  all_model_summaries[[param_i]] <- model_sel
  write.csv(model_sel,
            file.path(out_dir_param, paste0(param_clean, "_model_selection.csv")),
            row.names = FALSE)
  
  if(nrow(model_sel) == 0){
    cat("No valid models for", param_i, "\n")
    next
  }
  
  best <- model_sel %>% slice(1)
  best_name <- paste0("m", best$m, "_", make.names(best$R_structure))
  best_fit <- fits_out$fits[[best_name]]
  
  saveRDS(best_fit, file.path(out_dir_param, paste0(param_clean, "_best_dfa_model.rds")))
  
  sink(file.path(out_dir_param, paste0(param_clean, "_best_model_summary.txt")))
  cat("Parameter:", param_i, "\n")
  cat("Filter:", filter_label, "\n")
  cat("Time aggregation:", time_agg, "\n\n")
  print(best)
  cat("\nMODEL OBJECT:\n")
  print(best_fit)
  cat("\nTIDY:\n")
  print(tryCatch(broom::tidy(best_fit), error = function(e) "broom::tidy not available"))
  cat("\nGLANCE:\n")
  print(tryCatch(broom::glance(best_fit), error = function(e) "broom::glance not available"))
  sink()
  
  png(file.path(out_dir_param, paste0(param_clean, "_MARSS_default_plot.png")),
      width = 1600, height = 1200, res = 200)
  plot(best_fit)
  dev.off()
  
  trends <- get_fitted_trends(best_fit)
  Z <- get_loadings(best_fit)
  
  write.csv(Z, file.path(out_dir_param, paste0(param_clean, "_loadings_Z.csv")))
  write.csv(trends, file.path(out_dir_param, paste0(param_clean, "_common_trends.csv")))
  
  plot_loadings(
    Z = Z,
    file_out = file.path(out_dir_param, paste0(param_clean, "_loadings_heatmap.png")),
    parameter_name = param_i
  )
  
  plot_common_trends(
    trends = trends,
    time_labels = time_labels,
    file_out = file.path(out_dir_param, paste0(param_clean, "_common_trends_plot.png")),
    parameter_name = param_i
  )
  
  plot_series_vs_fitted(
    y = y,
    fit = best_fit,
    time_labels = time_labels,
    file_out = file.path(out_dir_param, paste0(param_clean, "_observed_vs_fitted.png")),
    parameter_name = param_i
  )
  
  cors <- plot_canonical_correlations(
    y = y,
    trends = trends,
    file_out = file.path(out_dir_param, paste0(param_clean, "_series_trends_correlations.png")),
    parameter_name = param_i
  )
  
  write.csv(cors,
            file.path(out_dir_param, paste0(param_clean, "_series_trends_correlations.csv")))
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