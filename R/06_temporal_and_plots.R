# Plotting helpers and temporal outputs.

interview_theme <- function() {
  ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 18, color = "#172026"),
      plot.subtitle = ggplot2::element_text(size = 11, color = "#425466"),
      plot.caption = ggplot2::element_text(size = 9, color = "#6b7280"),
      axis.title = ggplot2::element_text(face = "bold", color = "#172026"),
      axis.text = ggplot2::element_text(color = "#334155"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(color = "#d8dee8", linewidth = 0.35),
      plot.background = ggplot2::element_rect(fill = "#f7f4ee", color = NA),
      panel.background = ggplot2::element_rect(fill = "#fcfbf7", color = NA),
      legend.background = ggplot2::element_rect(fill = "#fcfbf7", color = NA),
      legend.key = ggplot2::element_rect(fill = "#fcfbf7", color = NA)
    )
}


interview_pretty_label <- function(x) {
  out <- gsub("_", " ", x)
  out <- gsub("\\bno2\\b", "NO2", out, ignore.case = TRUE)
  out <- gsub("\\bndvi\\b", "NDVI", out, ignore.case = TRUE)
  out <- gsub("\\bgndvi\\b", "GNDVI", out, ignore.case = TRUE)
  out <- gsub("\\bndwi\\b", "NDWI", out, ignore.case = TRUE)
  out <- gsub("\\bevi\\b", "EVI", out, ignore.case = TRUE)
  out <- gsub("\\bsavi\\b", "SAVI", out, ignore.case = TRUE)
  out <- gsub("\\bmean\\b", "Mean", out, ignore.case = TRUE)
  out <- gsub("\\blatest\\b", "Latest", out, ignore.case = TRUE)
  out <- gsub("\\bfitted\\b", "Fitted", out, ignore.case = TRUE)
  out <- gsub("\\bresidual\\b", "Residual", out, ignore.case = TRUE)
  trimws(out)
}


interview_value_palette <- function() {
  c("#fff5f0", "#fcbba1", "#fb6a4a", "#cb181d", "#67000d")
}


interview_residual_palette <- function() {
  c("#2166ac", "#f7f7f7", "#b2182b")
}


interview_index_palette <- function(index_name) {
  index_name <- tolower(index_name)
  if (index_name %in% c("ndvi", "evi", "gndvi", "savi")) {
    return(c("#6b4f2a", "#c2b280", "#dfe8b2", "#8dc26f", "#1f6b3a"))
  }
  if (index_name %in% c("ndwi")) {
    return(c("#5c4d3a", "#dfe8f2", "#7bb6d9", "#2c7fb8", "#0b4f6c"))
  }
  c("#f7f7f7", "#cccccc", "#4d4d4d")
}


interview_quantile_factor <- function(x, n = 5, prefix = "Q") {
  probs <- seq(0, 1, length.out = n + 1)
  breaks <- unique(stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE, type = 7))
  if (length(breaks) < 3) {
    return(factor(rep("Single range", length(x))))
  }
  cut(
    x,
    breaks = breaks,
    include.lowest = TRUE,
    ordered_result = TRUE,
    labels = paste0(prefix, seq_len(length(breaks) - 1))
  )
}


interview_prepare_plot_sf <- function(sf_obj) {
  sf::st_transform(sf_obj, 4326)
}


interview_plot_sensor_surface <- function(surface_raster, points_sf, value_col, output_file, config) {
  vals <- terra::as.data.frame(surface_raster, xy = TRUE, na.rm = TRUE)
  names(vals)[3] <- "value"
  pts <- interview_prepare_plot_sf(points_sf)
  p <- ggplot2::ggplot() +
    ggplot2::geom_raster(data = vals, ggplot2::aes(x = x, y = y, fill = value)) +
    ggplot2::geom_contour(data = vals, ggplot2::aes(x = x, y = y, z = value), bins = 10, color = "#ffffff", alpha = 0.45, linewidth = 0.2) +
    ggplot2::geom_sf(data = pts, ggplot2::aes(color = .data[[value_col]]), size = 1.5, alpha = 0.85) +
    ggplot2::scale_fill_gradientn(colors = interview_value_palette(), name = interview_pretty_label(value_col)) +
    ggplot2::scale_color_gradientn(colors = interview_value_palette(), name = interview_pretty_label(value_col)) +
    ggplot2::coord_sf(expand = FALSE) +
    interview_theme() +
    ggplot2::labs(
      title = paste(interview_pretty_label(value_col), "interpolated surface"),
      subtitle = paste0("IDW on a ", config$surfaces$resolution, " m grid"),
      caption = "Exploratory interpolation from monitoring sites."
    )
  ggplot2::ggsave(output_file, p, width = config$plots$width, height = config$plots$height, dpi = config$plots$dpi)
}


interview_plot_scatter <- function(dat, predictor_var, response_var, output_file, config) {
  fit <- stats::lm(stats::as.formula(paste(response_var, "~", predictor_var)), data = dat)
  glance <- broom::glance(fit)
  p <- ggplot2::ggplot(dat, ggplot2::aes(x = .data[[predictor_var]], y = .data[[response_var]])) +
    ggplot2::geom_point(color = "#7f1d1d", fill = "#fca5a5", shape = 21, stroke = 0.3, size = 2.2, alpha = 0.85) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "#991b1b", fill = "#fecaca") +
    interview_theme() +
    ggplot2::labs(
      title = paste(interview_pretty_label(response_var), "vs", interview_pretty_label(predictor_var)),
      subtitle = paste0("Linear screening relationship | R^2 = ", format(round(glance$r.squared, 3), nsmall = 3)),
      x = interview_pretty_label(predictor_var),
      y = interview_pretty_label(response_var)
    )
  ggplot2::ggsave(output_file, p, width = config$plots$width, height = config$plots$height, dpi = config$plots$dpi)
}


interview_plot_point_map <- function(sf_obj, value_col, output_file, title, diverging = FALSE, midpoint = 0, config, legend_label = NULL) {
  pts <- interview_prepare_plot_sf(sf_obj)
  legend_label <- legend_label %||% interview_pretty_label(value_col)
  p <- ggplot2::ggplot(pts) +
    ggplot2::geom_sf(ggplot2::aes(color = .data[[value_col]]), size = 2.2, alpha = 0.9) +
    ggplot2::coord_sf(expand = FALSE) +
    interview_theme() +
    ggplot2::labs(title = title, color = legend_label, x = NULL, y = NULL)

  if (diverging) {
    p <- p + ggplot2::scale_color_gradient2(
      low = interview_residual_palette()[1],
      mid = interview_residual_palette()[2],
      high = interview_residual_palette()[3],
      midpoint = midpoint,
      name = legend_label
    )
  } else {
    p <- p + ggplot2::scale_color_gradientn(colors = interview_value_palette(), name = legend_label)
  }

  ggplot2::ggsave(output_file, p, width = config$plots$width, height = config$plots$height, dpi = config$plots$dpi)
}


interview_plot_quantile_map <- function(sf_obj, value_col, output_file, title, config, legend_label = NULL) {
  pts <- interview_prepare_plot_sf(sf_obj)
  legend_label <- legend_label %||% interview_pretty_label(value_col)
  pts$quantile_class <- interview_quantile_factor(pts[[value_col]], n = 5, prefix = "Q")
  palette <- c("#eff3ff", "#bdd7e7", "#6baed6", "#3182bd", "#08519c")
  p <- ggplot2::ggplot(pts) +
    ggplot2::geom_sf(ggplot2::aes(color = quantile_class), size = 2.2, alpha = 0.9) +
    ggplot2::coord_sf(expand = FALSE) +
    interview_theme() +
    ggplot2::labs(title = title, color = legend_label, x = NULL, y = NULL) +
    ggplot2::scale_color_manual(values = palette[seq_len(length(levels(pts$quantile_class)))], drop = FALSE)
  ggplot2::ggsave(output_file, p, width = config$plots$width, height = config$plots$height, dpi = config$plots$dpi)
}


interview_plot_index_with_points <- function(index_raster_path, points_sf, point_value_col, output_file, index_name, config) {
  r <- terra::rast(index_raster_path)
  df <- terra::as.data.frame(r, xy = TRUE, na.rm = TRUE)
  names(df)[3] <- "index_value"
  pts <- sf::st_transform(points_sf, sf::st_crs(terra::crs(r)))
  p <- ggplot2::ggplot() +
    ggplot2::geom_raster(data = df, ggplot2::aes(x = x, y = y, fill = index_value), alpha = 0.95) +
    ggplot2::geom_sf(data = pts, ggplot2::aes(color = .data[[point_value_col]]), size = 1.8, alpha = 0.88) +
    ggplot2::coord_sf(expand = FALSE) +
    interview_theme() +
    ggplot2::labs(
      title = paste(interview_pretty_label(index_name), "with monitoring sites"),
      subtitle = paste("Background:", interview_pretty_label(index_name), "| Points:", interview_pretty_label(point_value_col)),
      fill = interview_pretty_label(index_name),
      color = interview_pretty_label(point_value_col),
      x = NULL,
      y = NULL
    ) +
    ggplot2::scale_fill_gradientn(colors = interview_index_palette(index_name)) +
    ggplot2::scale_color_gradientn(colors = interview_value_palette())
  ggplot2::ggsave(output_file, p, width = config$plots$width, height = config$plots$height, dpi = config$plots$dpi)
}


interview_plot_residual_vs_fitted <- function(dat, output_file, config) {
  p <- ggplot2::ggplot(dat, ggplot2::aes(x = fitted_value, y = residual_value)) +
    ggplot2::geom_point(color = "#1d4ed8", size = 2.1, alpha = 0.85) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "#b91c1c") +
    interview_theme() +
    ggplot2::labs(
      title = "Residuals vs fitted values",
      subtitle = "Residuals should center around zero without strong structure",
      x = "Fitted value",
      y = "Residual"
    )
  ggplot2::ggsave(output_file, p, width = config$plots$width, height = config$plots$height, dpi = config$plots$dpi)
}


interview_plot_time_series <- function(df, time_col, value_col, output_file, config) {
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[time_col]], y = .data[[value_col]])) +
    ggplot2::geom_line(color = "#991b1b", linewidth = 0.8) +
    ggplot2::geom_point(color = "#7f1d1d", size = 2) +
    interview_theme() +
    ggplot2::labs(
      title = paste(gsub("_", " ", value_col), "through time"),
      x = "Time",
      y = interview_pretty_label(value_col)
    )
  ggplot2::ggsave(output_file, p, width = config$plots$width, height = config$plots$height, dpi = config$plots$dpi)
}


interview_plot_gwr_recommendation <- function(decision_df, output_file, config) {
  palette <- c(
    "Suggested" = "#1f9d55",
    "Use cautiously" = "#d97706",
    "Not suggested" = "#b91c1c"
  )
  p <- ggplot2::ggplot(decision_df, ggplot2::aes(x = metric, y = 1, fill = status)) +
    ggplot2::geom_tile(color = "#ffffff", linewidth = 0.7, height = 0.65) +
    ggplot2::geom_text(ggplot2::aes(label = short_note), color = "#172026", size = 3.7) +
    ggplot2::scale_fill_manual(values = palette, drop = FALSE) +
    interview_theme() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      legend.position = "none"
    ) +
    ggplot2::labs(
      title = "GWR recommendation summary",
      subtitle = "Quick check on whether a local model is worth pursuing"
    )
  ggplot2::ggsave(output_file, p, width = 11, height = 3.8, dpi = config$plots$dpi)
}


run_temporal_analysis <- function(point_integration_sf, output_dir, config) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)

  dat <- sf::st_drop_geometry(point_integration_sf)
  time_cols <- grep("_scene_datetime$", names(dat), value = TRUE)
  if (length(time_cols) == 0) {
    writeLines("No temporal raster scene columns were available.", con = file.path(output_dir, "README.md"))
    return(NULL)
  }

  rows <- list()
  for (time_col in time_cols) {
    index_name <- sub("_scene_datetime$", "", time_col)
    value_col <- paste0(index_name, "_value")
    if (!value_col %in% names(dat)) next
    tmp <- dat[, c(time_col, value_col), drop = FALSE]
    names(tmp) <- c("scene_datetime", "value")
    tmp$index_name <- index_name
    rows[[length(rows) + 1]] <- tmp
  }

  out <- do.call(rbind, rows)
  out <- out[!is.na(out$scene_datetime) & is.finite(out$value), , drop = FALSE]
  if (nrow(out) == 0) {
    writeLines("Temporal analysis was skipped because no raster scene timestamps were available.", con = file.path(output_dir, "README.md"))
    return(NULL)
  }
  summary_df <- out |>
    dplyr::group_by(index_name, scene_datetime) |>
    dplyr::summarise(mean_value = mean(value, na.rm = TRUE), median_value = stats::median(value, na.rm = TRUE), .groups = "drop")

  utils::write.csv(summary_df, file.path(output_dir, "tables", "temporal_summary.csv"), row.names = FALSE)
  unique_indices <- unique(summary_df$index_name)
  for (index_name in unique_indices) {
    tmp <- summary_df[summary_df$index_name == index_name, , drop = FALSE]
    if (nrow(tmp) < 2) next
    interview_plot_time_series(tmp, "scene_datetime", "mean_value", file.path(output_dir, "figures", paste0(index_name, "_temporal_profile.png")), config)
  }
  summary_df
}
