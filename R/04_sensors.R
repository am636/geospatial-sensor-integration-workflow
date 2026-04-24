# Sensor preparation and CRS inference.

interview_drop_empty_columns <- function(dat) {
  keep <- !vapply(dat, function(x) all(is.na(x) | trimws(as.character(x)) == ""), logical(1))
  dat[, keep, drop = FALSE]
}


interview_parse_sensor_time <- function(x) {
  suppressWarnings(
    lubridate::parse_date_time(
      x,
      orders = c("Ymd HMS", "Ymd HM", "Ymd", "Y-m-d H:M:S", "Y-m-d", "d/m/Y H:M:S", "d/m/Y"),
      tz = "UTC"
    )
  )
}


interview_score_crs_candidate <- function(df, x_col, y_col, candidate_crs, raster_crs, raster_extent) {
  sf_obj <- tryCatch(
    sf::st_as_sf(df, coords = c(x_col, y_col), crs = candidate_crs, remove = FALSE),
    error = function(e) NULL
  )
  if (is.null(sf_obj)) {
    return(data.frame(candidate_crs = candidate_crs, inside_fraction = 0, inside_points = 0))
  }

  sf_obj <- tryCatch(sf::st_transform(sf_obj, sf::st_crs(raster_crs)), error = function(e) NULL)
  if (is.null(sf_obj)) {
    return(data.frame(candidate_crs = candidate_crs, inside_fraction = 0, inside_points = 0))
  }

  coords <- sf::st_coordinates(sf_obj)
  inside <- coords[, 1] >= terra::xmin(raster_extent) &
    coords[, 1] <= terra::xmax(raster_extent) &
    coords[, 2] >= terra::ymin(raster_extent) &
    coords[, 2] <= terra::ymax(raster_extent)
  inside[is.na(inside)] <- FALSE

  data.frame(
    candidate_crs = candidate_crs,
    inside_fraction = mean(inside),
    inside_points = sum(inside),
    stringsAsFactors = FALSE
  )
}


interview_infer_sensor_crs <- function(df, raster_reference_path, sensors_config) {
  if (!is.null(sensors_config$manual_source_crs)) {
    return(list(source_crs = sensors_config$manual_source_crs, score_table = NULL))
  }

  ref <- terra::rast(raster_reference_path)
  candidates <- lapply(
    sensors_config$candidate_crs,
    function(crs_code) interview_score_crs_candidate(
      df = df,
      x_col = sensors_config$x_col,
      y_col = sensors_config$y_col,
      candidate_crs = crs_code,
      raster_crs = terra::crs(ref),
      raster_extent = terra::ext(ref)
    )
  )
  score_table <- do.call(rbind, candidates)
  score_table <- score_table[order(-score_table$inside_fraction, -score_table$inside_points), , drop = FALSE]
  list(source_crs = score_table$candidate_crs[1], score_table = score_table)
}


interview_prepare_sensor_table <- function(sensor_csv_path, raster_reference_path, config, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_dir, "vectors"), recursive = TRUE, showWarnings = FALSE)

  dat <- utils::read.csv(sensor_csv_path, stringsAsFactors = FALSE)
  dat <- interview_sanitize_data_frame(dat)
  dat <- interview_drop_empty_columns(dat)
  dat <- dat[is.finite(dat[[config$sensors$x_col]]) & is.finite(dat[[config$sensors$y_col]]), , drop = FALSE]

  inferred <- interview_infer_sensor_crs(dat, raster_reference_path, config$sensors)
  dat$timestamp_utc <- if (!is.null(config$sensors$time_col) && config$sensors$time_col %in% names(dat)) {
    interview_parse_sensor_time(dat[[config$sensors$time_col]])
  } else if (!is.null(config$sensors$year_col) && config$sensors$year_col %in% names(dat)) {
    as.POSIXct(paste0(dat[[config$sensors$year_col]], "-01-01 00:00:00"), tz = "UTC")
  } else {
    as.POSIXct(config$fallback_scene_datetime, tz = "UTC")
  }

  dat$sensor_id <- dat[[config$sensors$id_col]]

  source_sf <- sf::st_as_sf(
    dat,
    coords = c(config$sensors$x_col, config$sensors$y_col),
    crs = inferred$source_crs,
    remove = FALSE
  )

  raster_crs <- terra::crs(terra::rast(raster_reference_path))
  raster_epsg <- sf::st_crs(raster_crs)$epsg %||% NA_integer_
  raster_sf <- sf::st_transform(source_sf, sf::st_crs(raster_crs))
  lonlat_sf <- sf::st_transform(raster_sf, 4326)
  lonlat_coords <- sf::st_coordinates(lonlat_sf)
  raster_sf$lon <- lonlat_coords[, 1]
  raster_sf$lat <- lonlat_coords[, 2]

  alias_map <- config$sensors$value_aliases
  value_cols <- intersect(config$sensors$value_cols, names(raster_sf))
  site_fields <- c(config$sensors$id_col, config$sensors$x_col, config$sensors$y_col)
  raw_dat <- sf::st_drop_geometry(raster_sf)

  site_summary <- raw_dat |>
    dplyr::group_by(.data[[config$sensors$id_col]], .data[[config$sensors$x_col]], .data[[config$sensors$y_col]]) |>
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(value_cols),
        list(
          mean = ~ mean(.x, na.rm = TRUE),
          latest = ~ .x[which.max(timestamp_utc)]
        ),
        .names = "{.fn}_{.col}"
      ),
      latest_year = suppressWarnings(max(.data[[config$sensors$year_col]], na.rm = TRUE)),
      n_records = dplyr::n(),
      lon = dplyr::first(lon),
      lat = dplyr::first(lat),
      timestamp_utc = max(timestamp_utc, na.rm = TRUE),
      .groups = "drop"
    )

  names(site_summary)[1:3] <- c("sensor_id", config$sensors$x_col, config$sensors$y_col)
  for (src_name in names(alias_map)) {
    alias <- alias_map[[src_name]]
    old_mean <- paste0("mean_", src_name)
    old_latest <- paste0("latest_", src_name)
    if (old_mean %in% names(site_summary)) names(site_summary)[names(site_summary) == old_mean] <- paste0("mean_", alias)
    if (old_latest %in% names(site_summary)) names(site_summary)[names(site_summary) == old_latest] <- paste0("latest_", alias)
  }

  site_sf <- sf::st_as_sf(
    site_summary,
    coords = c(config$sensors$x_col, config$sensors$y_col),
    crs = inferred$source_crs,
    remove = FALSE
  )
  site_sf <- sf::st_transform(site_sf, sf::st_crs(raster_crs))

  utils::write.csv(raw_dat, file.path(output_dir, "tables", "sensor_records_clean.csv"), row.names = FALSE)
  utils::write.csv(site_summary, file.path(output_dir, "tables", "sensor_sites_summary.csv"), row.names = FALSE)
  sf::st_write(raster_sf, file.path(output_dir, "vectors", "sensor_records_raster_crs.gpkg"), delete_dsn = TRUE, quiet = TRUE)
  sf::st_write(site_sf, file.path(output_dir, "vectors", "sensor_sites_summary_raster_crs.gpkg"), delete_dsn = TRUE, quiet = TRUE)
  if (!is.null(inferred$score_table)) {
    utils::write.csv(inferred$score_table, file.path(output_dir, "tables", "sensor_crs_candidate_scores.csv"), row.names = FALSE)
  }

  readme_lines <- c(
    "# Sensors Output",
    "",
    paste0("- Source CRS used: EPSG:", inferred$source_crs),
    paste0("- Raster CRS: EPSG:", raster_epsg),
    "- `sensor_records_clean.csv` keeps the cleaned source records.",
    "- `sensor_sites_summary.csv` stores one summary row per monitoring site.",
    "- Historical records stay available in the cleaned table before site-level summarisation.",
    "- Repeated site observations are collapsed into mean and latest summaries.",
    "- The GeoPackage layers are written in the raster CRS."
  )
  writeLines(readme_lines, con = file.path(output_dir, "README.md"))

  list(
    source_crs = inferred$source_crs,
    score_table = inferred$score_table,
    sensor_records = raster_sf,
    sensor_sites = site_sf,
    sensor_records_table = raw_dat,
    sensor_sites_table = site_summary
  )
}
