# Regression, Moran's I, and GWR diagnostics.

run_exploratory_analysis <- function(point_integration_sf, output_dir, config) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_dir, "vectors"), recursive = TRUE, showWarnings = FALSE)

  response_var <- config$exploratory$response_var
  predictor_var <- config$exploratory$predictor_var
  dat <- point_integration_sf
  keep <- is.finite(dat[[response_var]]) & is.finite(dat[[predictor_var]])
  dat <- dat[keep, , drop = FALSE]
  if (nrow(dat) < 10) {
    stop("Not enough complete cases for exploratory analysis.")
  }

  model_df <- sf::st_drop_geometry(dat)
  fit <- stats::lm(stats::as.formula(paste(response_var, "~", predictor_var)), data = model_df)
  model_df$fitted_value <- stats::fitted(fit)
  model_df$residual_value <- stats::residuals(fit)
  dat$fitted_value <- model_df$fitted_value
  dat$residual_value <- model_df$residual_value

  regression_summary <- broom::glance(fit)
  regression_terms <- broom::tidy(fit)
  utils::write.csv(regression_summary, file.path(output_dir, "tables", "regression_summary.csv"), row.names = FALSE)
  utils::write.csv(regression_terms, file.path(output_dir, "tables", "regression_terms.csv"), row.names = FALSE)

  interview_plot_scatter(model_df, predictor_var, response_var, file.path(output_dir, "figures", "scatter_relationship.png"), config)
  interview_plot_residual_vs_fitted(model_df, file.path(output_dir, "figures", "residual_vs_fitted.png"), config)
  interview_plot_point_map(dat, response_var, file.path(output_dir, "figures", "observed_map.png"), paste("Observed", interview_pretty_label(response_var), "at monitoring sites"), FALSE, 0, config, legend_label = interview_pretty_label(response_var))
  interview_plot_point_map(dat, "fitted_value", file.path(output_dir, "figures", "fitted_map.png"), paste("Fitted", interview_pretty_label(response_var), "at monitoring sites"), FALSE, 0, config, legend_label = paste("Fitted", interview_pretty_label(response_var)))
  interview_plot_point_map(dat, "residual_value", file.path(output_dir, "figures", "residual_map.png"), "Regression residuals at monitoring sites", TRUE, 0, config, legend_label = paste("Residual", interview_pretty_label(response_var)))

  coords <- sf::st_coordinates(sf::st_transform(dat, sf::st_crs(terra::crs(terra::vect(dat)))))
  knn <- spdep::knearneigh(coords, k = min(config$exploratory$k_neighbors, nrow(dat) - 1))
  nb <- spdep::knn2nb(knn)
  lw <- spdep::nb2listw(nb, style = "W")

  moran_global <- spdep::moran.test(model_df$residual_value, lw)
  moran_perm <- spdep::moran.mc(model_df$residual_value, lw, nsim = config$exploratory$permutations)
  moran_df <- data.frame(
    statistic = unname(moran_global$statistic),
    p_value = moran_global$p.value,
    simulated_p_value = moran_perm$p.value,
    stringsAsFactors = FALSE
  )
  utils::write.csv(moran_df, file.path(output_dir, "tables", "morans_i_residuals.csv"), row.names = FALSE)

  gwr_outputs <- list()
  if (isTRUE(config$exploratory$run_gwr) && requireNamespace("GWmodel", quietly = TRUE)) {
    sp_data <- methods::as(dat, "Spatial")
    gwr_formula <- stats::as.formula(paste(response_var, "~", predictor_var))
    bw <- GWmodel::bw.gwr(
      formula = gwr_formula,
      data = sp_data,
      approach = "AICc",
      kernel = "bisquare",
      adaptive = TRUE
    )
    gwr_fit <- GWmodel::gwr.basic(
      formula = gwr_formula,
      data = sp_data,
      bw = bw,
      kernel = "bisquare",
      adaptive = TRUE
    )

    gwr_sf <- sf::st_as_sf(gwr_fit$SDF)
    sf::st_write(gwr_sf, file.path(output_dir, "vectors", "gwr_outputs.gpkg"), delete_dsn = TRUE, quiet = TRUE)

    coef_col <- predictor_var
    if (coef_col %in% names(gwr_sf)) {
      interview_plot_point_map(gwr_sf, coef_col, file.path(output_dir, "figures", "gwr_local_coefficient_map.png"), paste("GWR local coefficient for", interview_pretty_label(predictor_var)), TRUE, 0, config, legend_label = paste("GWR coefficient:", interview_pretty_label(predictor_var)))
    }
    if ("Local_R2" %in% names(gwr_sf)) {
      interview_plot_quantile_map(gwr_sf, "Local_R2", file.path(output_dir, "figures", "gwr_local_r2_map.png"), "GWR local R2 quantile classes", config, legend_label = "Local R2 quantiles")
    }

    gwr_table <- as.data.frame(gwr_sf)
    utils::write.csv(gwr_table, file.path(output_dir, "tables", "gwr_outputs.csv"), row.names = FALSE)

    gwr_resid_col <- grep("residual", names(gwr_sf), ignore.case = TRUE, value = TRUE)
    if (length(gwr_resid_col) > 0) {
      gwr_resid <- gwr_sf[[gwr_resid_col[1]]]
      gwr_moran <- spdep::moran.test(gwr_resid, lw)
      gwr_outputs$moran <- data.frame(
        statistic = unname(gwr_moran$statistic),
        p_value = gwr_moran$p.value,
        stringsAsFactors = FALSE
      )
      utils::write.csv(gwr_outputs$moran, file.path(output_dir, "tables", "morans_i_gwr_residuals.csv"), row.names = FALSE)
    }

    local_r2_iqr <- if ("Local_R2" %in% names(gwr_sf)) {
      stats::IQR(gwr_sf$Local_R2, na.rm = TRUE)
    } else {
      NA_real_
    }

    gwr_decision <- data.frame(
      metric = c("Residual spatial pattern", "After GWR residual pattern", "Local model variability", "Overall recommendation"),
      status = c(
        if (moran_df$simulated_p_value[1] < 0.05) "Suggested" else "Not suggested",
        if (!is.null(gwr_outputs$moran) && gwr_outputs$moran$p_value[1] < 0.05) "Use cautiously" else "Suggested",
        if (!is.na(local_r2_iqr) && local_r2_iqr > 0.08) "Suggested" else "Use cautiously",
        if (moran_df$simulated_p_value[1] < 0.05 && !is.null(gwr_outputs$moran) && gwr_outputs$moran$p_value[1] < 0.05) "Use cautiously" else if (moran_df$simulated_p_value[1] < 0.05) "Suggested" else "Not suggested"
      ),
      short_note = c(
        if (moran_df$simulated_p_value[1] < 0.05) "Residuals cluster" else "Weak clustering",
        if (!is.null(gwr_outputs$moran) && gwr_outputs$moran$p_value[1] < 0.05) "Still clustered" else "Mostly reduced",
        if (!is.na(local_r2_iqr) && local_r2_iqr > 0.08) "Spatial variation present" else "Limited local contrast",
        if (moran_df$simulated_p_value[1] < 0.05 && !is.null(gwr_outputs$moran) && gwr_outputs$moran$p_value[1] < 0.05) "Use with caution" else if (moran_df$simulated_p_value[1] < 0.05) "Reasonable next step" else "Stay global"
      ),
      stringsAsFactors = FALSE
    )
    utils::write.csv(gwr_decision, file.path(output_dir, "tables", "gwr_recommendation.csv"), row.names = FALSE)
    interview_plot_gwr_recommendation(gwr_decision, file.path(output_dir, "figures", "gwr_recommendation.png"), config)
  }

  list(
    regression_summary = regression_summary,
    regression_terms = regression_terms,
    moran = moran_df,
    gwr = gwr_outputs
  )
}
