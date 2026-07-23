# Minimal climb power-duration modelling script
# Input columns expected: date, time, power, HR, VAM, name

library(readxl)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(grid)

# -----------------------------------------------------------------------------
# User inputs
# -----------------------------------------------------------------------------
file_path <- "climb_results.xlsx"

mode <- "watts_for_time"   # choose: "watts_for_time" or "time_for_watts"

user_time   <- "34:00"     # used when mode == "watts_for_time"
user_watts  <- 371         # used when mode == "time_for_watts"
user_weight <- 80          # user-defined body weight in kg

known_weights <- c("Lennart Herrlich" = 80)

plot_min_sec <- 30 * 60
plot_max_sec <- 90 * 60

# Model selection:
# "auto" chooses the best model by weighted LOOCV RMSE.
# Or manually set one of:
# "log_linear", "log_quadratic", "reciprocal",
# "reciprocal_quadratic", "combined_log_reciprocal", "power_law_like"
selected_model_name <- "auto"

# Weighting for model selection:
# Higher values emphasize fast times more strongly.
# 0.0 = ordinary LOOCV RMSE
# 1.0 = moderate fast-time emphasis
# 2.0 = strong fast-time emphasis
fast_time_weight_strength <- 2

# -----------------------------------------------------------------------------
# Shared plotting theme
# -----------------------------------------------------------------------------
discrete_palette_fallback <- function(n) {
  if (n <= 8) {
    return(RColorBrewer::brewer.pal(n, "Dark2"))
  }
  
  return(grDevices::hcl.colors(n, palette = "Dark 3"))
}

options(
  ggplot2.discrete.fill = function(...) {
    ggplot2::discrete_scale(
      aesthetics = "fill",
      palette = discrete_palette_fallback,
      ...
    )
  },
  ggplot2.discrete.colour = function(...) {
    ggplot2::discrete_scale(
      aesthetics = "colour",
      palette = discrete_palette_fallback,
      ...
    )
  }
)

grey_theme <- theme_gray(base_size = 9, base_family = "Arial") +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 9, margin = margin(b = 10)),
    plot.caption = element_text(size = 7, color = "grey30"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "grey20"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid.major = element_line(color = "white", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    legend.key.size = unit(0.4, "cm")
  )

theme_set(grey_theme)

dark2_curve_colour <- RColorBrewer::brewer.pal(8, "Dark2")[1]


# More y-axis breaks
axis_break_count <- 12

y_breaks_pretty <- function(x) {
  pretty(x, n = axis_break_count)
}

# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------
parse_time_to_sec <- function(x) {
  parts <- strsplit(as.character(x), ":")
  
  sapply(parts, function(p) {
    p <- as.numeric(p)
    
    if (length(p) == 2) return(p[1] * 60 + p[2])
    if (length(p) == 3) return(p[1] * 3600 + p[2] * 60 + p[3])
    
    NA_real_
  })
}

sec_to_time <- function(sec) {
  sec <- round(sec)
  
  h <- sec %/% 3600
  m <- (sec %% 3600) %/% 60
  s <- sec %% 60
  
  if (h > 0) {
    sprintf("%d:%02d:%02d", h, m, s)
  } else {
    sprintf("%d:%02d", m, s)
  }
}

clean_number <- function(x) {
  as.numeric(gsub("[^0-9.]", "", as.character(x)))
}

weighted_rmse <- function(actual, predicted, weights) {
  sqrt(weighted.mean((actual - predicted)^2, weights, na.rm = TRUE))
}

ordinary_rmse <- function(actual, predicted) {
  sqrt(mean((actual - predicted)^2, na.rm = TRUE))
}

# -----------------------------------------------------------------------------
# 1. Import and clean file
# -----------------------------------------------------------------------------
df <- read_excel(file_path) %>%
  rename_with(tolower) %>%
  mutate(
    time_sec = parse_time_to_sec(time),
    power = clean_number(power),
    hr = clean_number(hr),
    vam = clean_number(vam),
    weight = known_weights[name],
    power_to_weight = power / weight
  ) %>%
  filter(
    !is.na(time_sec),
    !is.na(power),
    !is.na(power_to_weight),
    time_sec >= plot_min_sec,
    time_sec <= plot_max_sec
  ) %>%
  arrange(time_sec)

if (nrow(df) < 6) {
  stop("Not enough valid rows between 30 and 90 minutes to compare models robustly.")
}

# Fast-time emphasis: shorter durations receive larger selection weights.
df <- df %>%
  mutate(
    selection_weight = (max(time_sec) / time_sec) ^ fast_time_weight_strength
  )

# -----------------------------------------------------------------------------
# 2. Candidate model definitions
# -----------------------------------------------------------------------------
candidate_models <- list(
  log_linear = power_to_weight ~ log(time_sec),
  
  log_quadratic = power_to_weight ~ log(time_sec) + I(log(time_sec)^2),
  
  reciprocal = power_to_weight ~ I(1 / time_sec),
  
  reciprocal_quadratic = power_to_weight ~ I(1 / time_sec) + I(1 / time_sec^2),
  
  combined_log_reciprocal = power_to_weight ~ log(time_sec) + I(1 / time_sec),
  
  power_law_like = power_to_weight ~ log(time_sec) + I(sqrt(1 / time_sec))
)

# -----------------------------------------------------------------------------
# 3. Model-selection phase
# -----------------------------------------------------------------------------
evaluate_model <- function(model_name, model_formula, data) {
  full_fit <- lm(model_formula, data = data)
  
  full_pred <- predict(full_fit, newdata = data)
  
  loocv_pred <- rep(NA_real_, nrow(data))
  
  for (i in seq_len(nrow(data))) {
    train_data <- data[-i, , drop = FALSE]
    test_data <- data[i, , drop = FALSE]
    
    loocv_fit <- tryCatch(
      lm(model_formula, data = train_data),
      error = function(e) NULL
    )
    
    if (!is.null(loocv_fit)) {
      loocv_pred[i] <- predict(loocv_fit, newdata = test_data)
    }
  }
  
  data.frame(
    model = model_name,
    n_parameters = length(coef(full_fit)),
    aic = AIC(full_fit),
    bic = BIC(full_fit),
    rmse_in_sample_wkg = ordinary_rmse(data$power_to_weight, full_pred),
    loocv_rmse_wkg = ordinary_rmse(data$power_to_weight, loocv_pred),
    weighted_loocv_rmse_wkg = weighted_rmse(
      actual = data$power_to_weight,
      predicted = loocv_pred,
      weights = data$selection_weight
    ),
    fast_edge_error_wkg = mean(
      full_pred[data$time_sec <= quantile(data$time_sec, 0.25)] -
        data$power_to_weight[data$time_sec <= quantile(data$time_sec, 0.25)],
      na.rm = TRUE
    )
  )
}

model_selection <- bind_rows(
  lapply(
    names(candidate_models),
    function(model_name) {
      evaluate_model(
        model_name = model_name,
        model_formula = candidate_models[[model_name]],
        data = df
      )
    }
  )
) %>%
  arrange(weighted_loocv_rmse_wkg)

write.csv(
  model_selection,
  "model_selection_results.csv",
  row.names = FALSE
)

if (selected_model_name == "auto") {
  selected_model_name <- model_selection$model[1]
}

if (!selected_model_name %in% names(candidate_models)) {
  stop("selected_model_name must be 'auto' or one of the candidate model names.")
}

fit <- lm(candidate_models[[selected_model_name]], data = df)

message("Selected model: ", selected_model_name)
message("Model-selection table saved to: model_selection_results.csv")

# -----------------------------------------------------------------------------
# 4. Generate fitted curve from 30 to 90 minutes
# -----------------------------------------------------------------------------
curve_grid <- tibble(
  time_sec = seq(plot_min_sec, plot_max_sec, length.out = 500)
)

curve_pred <- predict(
  fit,
  newdata = curve_grid,
  interval = "confidence"
) %>%
  as.data.frame()

curve_plot_data <- bind_cols(curve_grid, curve_pred)

# -----------------------------------------------------------------------------
# 5. Plot all candidate curves for model comparison
# -----------------------------------------------------------------------------
all_model_curves <- bind_rows(
  lapply(names(candidate_models), function(model_name) {
    model_fit <- lm(candidate_models[[model_name]], data = df)
    
    tibble(
      time_sec = curve_grid$time_sec,
      power_to_weight = predict(model_fit, newdata = curve_grid),
      model = model_name
    )
  })
)

p_model_comparison <- ggplot(df, aes(x = time_sec / 60, y = power_to_weight)) +
  geom_point(size = 1.8, alpha = 0.8) +
  geom_line(
    data = all_model_curves,
    aes(
      x = time_sec / 60,
      y = power_to_weight,
      colour = model
    ),
    inherit.aes = FALSE,
    linewidth = 0.75
  ) +
  scale_x_continuous(
    breaks = seq(30, 90, by = 2.5),
    labels = function(x) paste0(x, " min")
  ) +
  scale_y_continuous(
    breaks = y_breaks_pretty,
    labels = function(x) paste0(round(x, 2), " W/kg")
  ) +
  coord_cartesian(xlim = c(30, 90)) +
  labs(
    title = "Candidate model comparison",
    subtitle = paste0(
      "Selected model: ", selected_model_name,
      " | Selection metric: weighted LOOCV RMSE"
    ),
    x = "Time",
    y = "Power-to-weight",
    colour = "Model",
    caption = "Fast-time observations receive higher weight during model selection."
  )

ggsave(
  filename = "model_comparison_curves.png",
  plot = p_model_comparison,
  width = 9,
  height = 5,
  dpi = 300
)

# -----------------------------------------------------------------------------
# 6. Plot selected fitted curve with uncertainty ribbon
# -----------------------------------------------------------------------------
p_curve <- ggplot(df, aes(x = time_sec / 60, y = power_to_weight)) +
  geom_point(size = 1.8, alpha = 0.8) +
  geom_ribbon(
    data = curve_plot_data,
    aes(
      x = time_sec / 60,
      ymin = lwr,
      ymax = upr
    ),
    inherit.aes = FALSE,
    alpha = 0.25,
    fill = dark2_curve_colour
  ) +
  geom_line(
    data = curve_plot_data,
    aes(x = time_sec / 60, y = fit),
    inherit.aes = FALSE,
    linewidth = 0.9,
    colour = dark2_curve_colour
  ) +
  scale_x_continuous(
    breaks = seq(30, 90, by = 2.5),
    labels = function(x) paste0(x, " min")
  ) +
  scale_y_continuous(
    breaks = y_breaks_pretty,
    labels = function(x) paste0(round(x, 2), " W/kg")
  ) +
  coord_cartesian(xlim = c(30, 90)) +
  labs(
    title = "Selected power-to-weight curve fit",
    subtitle = paste0(
      "Model: ", selected_model_name,
      " | 95% confidence interval"
    ),
    x = "Time",
    y = "Power-to-weight",
    caption = "The selected model minimizes weighted LOOCV RMSE with extra emphasis on faster times."
  )

ggsave(
  filename = "selected_power_to_weight_curve_fit.png",
  plot = p_curve,
  width = 8,
  height = 5,
  dpi = 300
)

# -----------------------------------------------------------------------------
# 7. Helper for inverse prediction: estimate time from target W/kg
# -----------------------------------------------------------------------------
predict_wkg_at_time <- function(time_sec_value) {
  predict(
    fit,
    newdata = tibble(time_sec = time_sec_value)
  )[[1]]
}

predict_time_from_wkg <- function(target_wkg) {
  f_root <- function(t) {
    predict_wkg_at_time(t) - target_wkg
  }
  
  lower_value <- f_root(plot_min_sec)
  upper_value <- f_root(plot_max_sec)
  
  if (lower_value * upper_value > 0) {
    warning(
      "Target watts/kg is outside the modelled 30–90 min range. ",
      "Returning nearest boundary estimate."
    )
    
    boundary_times <- c(plot_min_sec, plot_max_sec)
    boundary_errors <- abs(sapply(boundary_times, f_root))
    
    return(boundary_times[which.min(boundary_errors)])
  }
  
  uniroot(
    f = f_root,
    interval = c(plot_min_sec, plot_max_sec)
  )$root
}

# -----------------------------------------------------------------------------
# 8. Predict watts for target time OR predict time for target watts
# -----------------------------------------------------------------------------
if (mode == "watts_for_time") {
  
  target_time_sec <- parse_time_to_sec(user_time)
  
  if (target_time_sec < plot_min_sec || target_time_sec > plot_max_sec) {
    warning("User-defined time is outside the plotted 30–90 min range.")
  }
  
  pred <- predict(
    fit,
    newdata = tibble(time_sec = target_time_sec),
    interval = "confidence"
  ) %>%
    as.data.frame()
  
  predicted_watts <- pred$fit * user_weight
  lower_watts <- pred$lwr * user_weight
  upper_watts <- pred$upr * user_weight
  
  message(
    "Predicted power for ", user_time, " at ", user_weight, " kg: ",
    round(predicted_watts), " W",
    " [", round(lower_watts), "–", round(upper_watts), " W]"
  )
  
  watts_curve_plot_data <- curve_plot_data %>%
    mutate(
      fit_watts = fit * user_weight,
      lwr_watts = lwr * user_weight,
      upr_watts = upr * user_weight
    )
  
  p_prediction <- ggplot() +
    geom_point(
      data = df,
      aes(x = time_sec / 60, y = power_to_weight * user_weight),
      size = 1.8,
      alpha = 0.8
    ) +
    geom_ribbon(
      data = watts_curve_plot_data,
      aes(
        x = time_sec / 60,
        ymin = lwr_watts,
        ymax = upr_watts
      ),
      alpha = 0.25,
      fill = dark2_curve_colour
    ) +
    geom_line(
      data = watts_curve_plot_data,
      aes(x = time_sec / 60, y = fit_watts),
      linewidth = 0.9,
      colour = dark2_curve_colour
    ) +
    geom_point(
      aes(x = target_time_sec / 60, y = predicted_watts),
      size = 2.7,
      colour = dark2_curve_colour
    ) +
    annotate(
      "text",
      x = target_time_sec / 60,
      y = predicted_watts,
      label = paste0(
        round(predicted_watts), " W at ", user_weight, " kg\n",
        "95% CI: ", round(lower_watts), "–", round(upper_watts), " W"
      ),
      vjust = -0.8,
      size = 3
    ) +
    scale_x_continuous(
      breaks = seq(30, 90, by = 2.5),
      labels = function(x) paste0(x, " min")
    ) +
    scale_y_continuous(
      breaks = y_breaks_pretty,
      labels = function(x) paste0(round(x), " W")
    ) +
    coord_cartesian(xlim = c(30, 90)) +
    labs(
      title = "Predicted watts for target time",
      subtitle = paste0(
        "Target time: ", user_time,
        " | Weight: ", user_weight, " kg",
        " | Model: ", selected_model_name
      ),
      x = "Time",
      y = "Power",
      caption = "Prediction curve is scaled from W/kg to watts using the user-defined body weight."
    )
  
  ggsave(
    filename = "prediction_watts_for_time.png",
    plot = p_prediction,
    width = 8,
    height = 5,
    dpi = 300
  )
  
} else if (mode == "time_for_watts") {
  
  target_power_to_weight <- user_watts / user_weight
  
  predicted_time_sec <- predict_time_from_wkg(target_power_to_weight)
  
  message(
    "Predicted time for ", user_watts, " W at ", user_weight, " kg: ",
    sec_to_time(predicted_time_sec)
  )
  
  watts_curve_plot_data <- curve_plot_data %>%
    mutate(
      fit_watts = fit * user_weight,
      lwr_watts = lwr * user_weight,
      upr_watts = upr * user_weight
    )
  
  p_prediction <- ggplot() +
    geom_point(
      data = df,
      aes(x = time_sec / 60, y = power_to_weight * user_weight),
      size = 1.8,
      alpha = 0.8
    ) +
    geom_ribbon(
      data = watts_curve_plot_data,
      aes(
        x = time_sec / 60,
        ymin = lwr_watts,
        ymax = upr_watts
      ),
      alpha = 0.25,
      fill = dark2_curve_colour
    ) +
    geom_line(
      data = watts_curve_plot_data,
      aes(x = time_sec / 60, y = fit_watts),
      linewidth = 0.9,
      colour = dark2_curve_colour
    ) +
    geom_point(
      aes(x = predicted_time_sec / 60, y = user_watts),
      size = 2.7,
      colour = dark2_curve_colour
    ) +
    annotate(
      "text",
      x = predicted_time_sec / 60,
      y = user_watts,
      label = paste0(
        sec_to_time(predicted_time_sec), "\n",
        user_watts, " W at ", user_weight, " kg"
      ),
      vjust = -0.8,
      size = 3
    ) +
    scale_x_continuous(
      breaks = seq(30, 90, by = 2.5),
      labels = function(x) paste0(x, " min")
    ) +
    scale_y_continuous(
      breaks = y_breaks_pretty,
      labels = function(x) paste0(round(x), " W")
    ) +
    coord_cartesian(xlim = c(30, 90)) +
    labs(
      title = "Predicted time for target watts",
      subtitle = paste0(
        "Target power: ", user_watts,
        " W | Weight: ", user_weight, " kg",
        " | Model: ", selected_model_name
      ),
      x = "Time",
      y = "Power",
      caption = "Inverse prediction is solved numerically within the 30–90 min model range."
    )
  
  ggsave(
    filename = "prediction_time_for_watts.png",
    plot = p_prediction,
    width = 8,
    height = 5,
    dpi = 300
  )
  
} else {
  stop("mode must be either 'watts_for_time' or 'time_for_watts'")
}