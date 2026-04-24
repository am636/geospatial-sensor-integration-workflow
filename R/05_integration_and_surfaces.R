# Point and polygon integration plus sensor surfaces.

interview_extract_value_from_raster <- function(raster_path, points_sf, method = "bilinear") {
  x <- terra::rast(raster_path)
  pts <- sf::st_transform(points_sf, sf::st_crs(terra::crs(x)))
  out <- terra::extract(x, terra::vect(pts), method = method, bind = FALSE)
  out[, 2]
}


interview_match_scene_times <- function(point_times, scene_times) {
  if (length(scene_times) == 0) {
    return(integer(length(point_times)))
  }
  if (length(scene_times) == 1) {
    return(rep(1L, length(point_times)))
  }
  if (all(is.na(scene_times)) || all(is.na(point_times))) {
    return(rep(1L, length(point_times)))
  }
  vapply(point_times, function(pt) {
    if (is.na(pt)) {
      return(1L)
    }
    which.min(abs(as.numeric(difftime(scene_times, pt, units = "secs"))))
  }, integer(1))
}


run_point_raster_integration <- function(sensor_sites_sf, index_inventory, output_dir, config) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_dir, "vectors"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_dir, "figures"), recursive = TRUE, showWarnings = FALSE)

  dat <- sensor_sites_sf
  selected_indices <- unique(index_inventory$index_name)

  for (index_name in selected_indices) {
    subset_idx <- index_inventory[index_inventory$index_name == index_name, , drop = FALSE]
    subset_idx <- subset_idx[order(subset_idx$scene_datetime), , drop = FALSE]
    scene_times <- subset_idx$scene_datetime
    sensor_times <- dat$timestamp_utc
    matched <- interview_match_scene_times(sensor_times, scene_times)
    dat[[paste0(index_name, "_scene_datetime")]] <- scene_times[matched]
    dat[[paste0(index_name, "_value")]] <- NA_real_

    for (i in seq_len(nrow(subset_idx))) {
      keep <- matched == i
      if (!any(keep)) next
      dat[[paste0(index_name, "_value")]][keep] <- interview_extract_value_from_raster(
        raster_path = subset_idx$raster_path[i],
        points_sf = dat[keep, , drop = FALSE],
        method = "bilinear"
      )
    }
  }

  sf::st_write(dat, file.path(output_dir, "vectors", "point_raster_integration.gpkg"), delete_dsn = TRUE, quiet = TRUE)
  utils::write.csv(sf::st_drop_geometry(dat), file.path(output_dir, "tables", "point_raster_integration.csv"), row.names = FALSE)

  point_label <- config$exploratory$response_var
  if (!point_label %in% names(dat)) {
    point_label <- intersect(c("mean_no2", "latest_no2"), names(dat))[1]
  }
  if (!is.na(point_label) && length(point_label) == 1) {
    for (index_name in selected_indices) {
      idx_rows <- index_inventory[index_inventory$index_name == index_name, , drop = FALSE]
      if (nrow(idx_rows) == 0) next
      interview_plot_index_with_points(
        index_raster_path = idx_rows$raster_path[1],
        points_sf = dat,
        point_value_col = point_label,
        output_file = file.path(output_dir, "figures", paste0(index_name, "_with_sensor_points.png")),
        index_name = index_name,
        config = config
      )
    }
  }

  readme_lines <- c(
    "# Point Raster Integration",
    "",
    "- This stage extracts all selected Sentinel-2 indices to the monitoring sites.",
    "- Sensor CRS is resolved upstream and the reprojected sensor layer is reused here.",
    "- If raster scene timestamps exist, extraction respects scene timing.",
    "- If raster scene timestamps are missing, extraction is treated as single-scene integration.",
    "- Historical sensor records are condensed into site summaries before extraction."
  )
  writeLines(readme_lines, con = file.path(output_dir, "README.md"))
  dat
}


run_polygon_raster_integration <- function(vector_inventory, index_inventory, output_dir, config) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_dir, "vectors"), recursive = TRUE, showWarnings = FALSE)

  if (nrow(vector_inventory) == 0) {
    writeLines("No polygon vector was available for polygon-raster integration.", con = file.path(output_dir, "README.md"))
    return(NULL)
  }

  polygon_row <- vector_inventory[grepl("POLYGON", toupper(vector_inventory$geometry_type)), , drop = FALSE]
  if (nrow(polygon_row) == 0) {
    writeLines("Vector data were present, but no polygon layer was detected.", con = file.path(output_dir, "README.md"))
    return(NULL)
  }

  polygons <- sf::st_read(polygon_row$file_path[1], quiet = TRUE)
  poly_id <- config$polygon$id_col %||% names(polygons)[1]
  results <- list()

  for (i in seq_len(nrow(index_inventory))) {
    r <- terra::rast(index_inventory$raster_path[i])
    polys_raster_crs <- sf::st_transform(polygons, sf::st_crs(terra::crs(r)))
    extracted <- terra::extract(r, terra::vect(polys_raster_crs), fun = mean, na.rm = TRUE, bind = TRUE)
    names(extracted)[ncol(extracted)] <- paste0(index_inventory$index_name[i], "_mean")
    df <- as.data.frame(extracted)
    df$scene_datetime <- index_inventory$scene_datetime[i]
    df$index_name <- index_inventory$index_name[i]
    results[[length(results) + 1]] <- df
  }

  out <- do.call(rbind, results)
  utils::write.csv(out, file.path(output_dir, "tables", "polygon_raster_integration.csv"), row.names = FALSE)
  out
}


interview_make_surface_template <- function(reference_raster_path, resolution = 100) {
  ref <- terra::rast(reference_raster_path)
  terra::rast(
    xmin = terra::xmin(ref),
    xmax = terra::xmax(ref),
    ymin = terra::ymin(ref),
    ymax = terra::ymax(ref),
    resolution = resolution,
    crs = terra::crs(ref)
  )
}


interview_idw_surface <- function(points_sf, value_col, template, idp = 2, nmax = 12) {
  pts <- points_sf[is.finite(points_sf[[value_col]]), , drop = FALSE]
  pts <- sf::st_transform(pts, sf::st_crs(terra::crs(template)))
  pts_sp <- methods::as(pts, "Spatial")

  xy <- terra::xyFromCell(template, 1:terra::ncell(template))
  grid_df <- data.frame(x = xy[, 1], y = xy[, 2], cell_id = seq_len(nrow(xy)))
  grid_sp <- sp::SpatialPixelsDataFrame(
    points = grid_df[, c("x", "y"), drop = FALSE],
    data = data.frame(cell_id = grid_df$cell_id),
    proj4string = sp::CRS(SRS_string = sf::st_crs(pts)$wkt)
  )

  form <- stats::as.formula(paste(value_col, "~ 1"))
  idw_sp <- gstat::idw(formula = form, locations = pts_sp, newdata = grid_sp, idp = idp, nmax = nmax)
  out <- terra::rast(idw_sp)
  pred_idx <- grep("pred$", names(out))
  if (length(pred_idx) == 0) pred_idx <- 1
  out <- out[[pred_idx[1]]]
  names(out) <- paste0(value_col, "_surface")
  out
}


run_sensor_surface_module <- function(sensor_sites_sf, reference_raster_path, output_dir, config) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_dir, "rasters"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)

  template <- interview_make_surface_template(reference_raster_path, resolution = config$surfaces$resolution)
  distance_raster <- terra::distance(template, terra::vect(sensor_sites_sf))
  terra::writeRaster(distance_raster, file.path(output_dir, "rasters", "distance_to_nearest_monitor_m.tif"), overwrite = TRUE)

  summary_rows <- list()
  for (value_col in config$surfaces$value_cols) {
    if (!value_col %in% names(sensor_sites_sf)) next
    surface <- interview_idw_surface(
      points_sf = sensor_sites_sf,
      value_col = value_col,
      template = template,
      idp = config$surfaces$idp,
      nmax = config$surfaces$nmax
    )
    raster_path <- file.path(output_dir, "rasters", paste0(value_col, "_idw_", config$surfaces$resolution, "m.tif"))
    terra::writeRaster(surface, raster_path, overwrite = TRUE)
    interview_plot_sensor_surface(surface, sensor_sites_sf, value_col, file.path(output_dir, "figures", paste0(value_col, "_surface_map.png")), config)
    vals <- terra::values(surface, mat = FALSE)
    vals <- vals[is.finite(vals)]
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      value_col = value_col,
      points_used = sum(is.finite(sensor_sites_sf[[value_col]])),
      raster_min = min(vals),
      raster_mean = mean(vals),
      raster_max = max(vals),
      stringsAsFactors = FALSE
    )
  }

  summary_df <- do.call(rbind, summary_rows)
  utils::write.csv(summary_df, file.path(output_dir, "tables", "sensor_surface_summary.csv"), row.names = FALSE)
  summary_df
}
