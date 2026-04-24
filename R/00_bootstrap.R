# Shared helpers and defaults.

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}


interview_require_packages <- function(extra = character(0)) {
  base_pkgs <- c(
    "terra", "sf", "ggplot2", "dplyr", "tidyr", "lubridate",
    "gstat", "sp", "spdep", "GWmodel", "broom"
  )
  needed <- unique(c(base_pkgs, extra))
  missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      "Missing required package(s): ",
      paste(missing, collapse = ", "),
      ". Install them before running the workflow."
    )
  }
  invisible(TRUE)
}


interview_source_modules <- function(project_dir) {
  module_dir <- file.path(project_dir, "R")
  module_files <- sort(list.files(module_dir, pattern = "^[0-9]{2}_.*\\.R$", full.names = TRUE))
  for (file_path in module_files) {
    if (basename(file_path) == "00_bootstrap.R") {
      next
    }
    source(file_path)
  }
  invisible(module_files)
}


interview_create_dirs <- function(project_dir) {
  dirs <- c(
    file.path(project_dir, "data"),
    file.path(project_dir, "data", "inbox"),
    file.path(project_dir, "output"),
    file.path(project_dir, "output", "00_initial_inspection"),
    file.path(project_dir, "output", "01_raster_harmonization"),
    file.path(project_dir, "output", "02_sentinel2_indices"),
    file.path(project_dir, "output", "03_sensors"),
    file.path(project_dir, "output", "04_point_raster_integration"),
    file.path(project_dir, "output", "05_polygon_raster_integration"),
    file.path(project_dir, "output", "06_sensor_surfaces"),
    file.path(project_dir, "output", "07_temporal_analysis"),
    file.path(project_dir, "output", "08_exploratory_analysis")
  )
  for (dir_path in dirs) {
    dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  }
  invisible(dirs)
}


interview_write_markdown_table <- function(df) {
  if (nrow(df) == 0) {
    return("_No rows available._")
  }

  dat <- df
  for (nm in names(dat)) {
    if (inherits(dat[[nm]], c("POSIXct", "POSIXt"))) {
      dat[[nm]] <- format(dat[[nm]], "%Y-%m-%d %H:%M:%S", tz = "UTC")
    }
  }
  dat[is.na(dat)] <- ""

  lines <- c(
    paste0("| ", paste(names(dat), collapse = " | "), " |"),
    paste0("|", paste(rep("---", ncol(dat)), collapse = "|"), "|")
  )

  for (i in seq_len(nrow(dat))) {
    lines <- c(lines, paste0("| ", paste(as.character(dat[i, ]), collapse = " | "), " |"))
  }

  paste(lines, collapse = "\n")
}


interview_sanitize_text <- function(x) {
  out <- iconv(as.character(x), from = "", to = "UTF-8", sub = "")
  out[is.na(out)] <- ""
  out
}


interview_sanitize_data_frame <- function(dat) {
  for (nm in names(dat)) {
    if (is.character(dat[[nm]]) || is.factor(dat[[nm]])) {
      dat[[nm]] <- interview_sanitize_text(dat[[nm]])
    }
  }
  dat
}


default_interview_config <- function(project_dir = normalizePath(getwd(), winslash = "/", mustWork = TRUE)) {
  list(
    project_dir = project_dir,
    project_name = "Geospatial Sensor Integration Workflow",
    author_name = "",
    author_email = "",
    data_dir = file.path(project_dir, "data", "inbox"),
    output_dir = file.path(project_dir, "output"),
    selected_indices = c("ndvi", "evi", "savi", "gndvi", "ndwi"),
    selected_sentinel_layers = NULL,
    fallback_scene_datetime = NULL,
    cloud_mask = list(
      enabled = FALSE,
      keep_scl_classes = c(4, 5, 6)
    ),
    harmonization = list(
      force = FALSE,
      manual_reference = NULL,
      manual_target_crs = NULL,
      continuous_method = "bilinear",
      categorical_method = "near"
    ),
    sensors = list(
      enabled = TRUE,
      x_col = "x",
      y_col = "y",
      time_col = NULL,
      year_col = "Year",
      id_col = "Location",
      value_cols = c("Concentration_ugm3"),
      value_aliases = c(Concentration_ugm3 = "no2"),
      manual_source_crs = NULL,
      candidate_crs = c(3857, 27700, 4326, 32630)
    ),
    polygon = list(
      enabled = TRUE,
      id_col = NULL,
      keep_attr_cols = NULL,
      stats = c("mean", "median")
    ),
    surfaces = list(
      enabled = FALSE,
      value_cols = c("mean_no2", "latest_no2"),
      resolution = 100,
      idp = 2,
      nmax = 12
    ),
    temporal = list(
      enabled = TRUE,
      max_series = 12
    ),
    exploratory = list(
      enabled = TRUE,
      response_var = "mean_no2",
      predictor_var = "ndvi_value",
      k_neighbors = 8,
      permutations = 499,
      run_gwr = TRUE
    ),
    plots = list(
      width = 10,
      height = 7,
      dpi = 220
    )
  )
}
