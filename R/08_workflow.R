# Workflow orchestration.

run_interview_workflow <- function(config) {
  interview_require_packages()
  interview_create_dirs(config$project_dir)

  inspect_dir <- file.path(config$output_dir, "00_initial_inspection")
  harmonize_dir <- file.path(config$output_dir, "01_raster_harmonization")
  index_dir <- file.path(config$output_dir, "02_sentinel2_indices")
  sensors_dir <- file.path(config$output_dir, "03_sensors")
  point_dir <- file.path(config$output_dir, "04_point_raster_integration")
  polygon_dir <- file.path(config$output_dir, "05_polygon_raster_integration")
  surface_dir <- file.path(config$output_dir, "06_sensor_surfaces")
  temporal_dir <- file.path(config$output_dir, "07_temporal_analysis")
  exploratory_dir <- file.path(config$output_dir, "08_exploratory_analysis")

  inspected <- run_auto_ingest_inspection(
    data_dir = config$data_dir,
    output_dir = inspect_dir,
    fallback_datetime = config$fallback_scene_datetime
  )

  if (nrow(inspected$raster_inventory) == 0) {
    stop("No raster files were found in data/inbox.")
  }

  sensor_table_path <- interview_select_sensor_table(inspected$table_inventory, config$sensors)
  if (is.null(sensor_table_path)) {
    stop("No suitable sensor table was detected in data/inbox.")
  }

  harmonized <- run_conditional_harmonization(inspected$raster_inventory, harmonize_dir, config$harmonization)
  index_inventory <- run_sentinel2_index_module(harmonized$inventory, config, index_dir)

  sensor_outputs <- interview_prepare_sensor_table(
    sensor_csv_path = sensor_table_path,
    raster_reference_path = harmonized$reference_path,
    config = config,
    output_dir = sensors_dir
  )

  point_outputs <- run_point_raster_integration(sensor_outputs$sensor_sites, index_inventory, point_dir, config)
  polygon_outputs <- if (isTRUE(config$polygon$enabled)) {
    run_polygon_raster_integration(inspected$vector_inventory, index_inventory, polygon_dir, config)
  } else {
    NULL
  }

  surface_outputs <- if (isTRUE(config$surfaces$enabled)) {
    run_sensor_surface_module(sensor_outputs$sensor_sites, harmonized$reference_path, surface_dir, config)
  } else {
    NULL
  }

  temporal_outputs <- if (isTRUE(config$temporal$enabled)) {
    run_temporal_analysis(point_outputs, temporal_dir, config)
  } else {
    NULL
  }

  predictor_name <- paste0(config$selected_indices[1], "_value")
  if (predictor_name %in% names(point_outputs)) {
    config$exploratory$predictor_var <- predictor_name
  }

  exploratory_outputs <- if (isTRUE(config$exploratory$enabled)) {
    run_exploratory_analysis(point_outputs, exploratory_dir, config)
  } else {
    NULL
  }

  run_summary <- data.frame(
    project_name = config$project_name,
    author_name = config$author_name,
    author_email = config$author_email,
    sensor_table = basename(sensor_table_path),
    raster_count = nrow(inspected$raster_inventory),
    vector_count = nrow(inspected$vector_inventory),
    scene_count = length(unique(index_inventory$scene_key)),
    indices_run = paste(unique(index_inventory$index_name), collapse = ", "),
    scene_times_found = sum(!is.na(index_inventory$scene_datetime)),
    harmonization_run = harmonized$harmonization_run,
    stringsAsFactors = FALSE
  )
  utils::write.csv(run_summary, file.path(config$output_dir, "run_summary.csv"), row.names = FALSE)

  readme_lines <- c(
    paste0("# ", config$project_name)
  )
  if (nzchar(config$author_name %||% "")) {
    readme_lines <- c(readme_lines, "", paste0("Author: ", config$author_name))
  }
  if (nzchar(config$author_email %||% "")) {
    readme_lines <- c(readme_lines, paste0("Email: ", config$author_email))
  }
  readme_lines <- c(
    readme_lines,
    "",
    "## Output folders",
    "",
    "- `00_initial_inspection`: input discovery and metadata checks",
    "- `01_raster_harmonization`: raster harmonization when required",
    "- `02_sentinel2_indices`: Sentinel-2 index rasters and maps",
    "- `03_sensors`: cleaned and reprojected sensor products",
    "- `04_point_raster_integration`: extracted raster values, linked attributes, and index-plus-point maps",
    "- `05_polygon_raster_integration`: polygon summaries with attribute carry-through",
    "- `06_sensor_surfaces`: optional interpolated sensor surfaces",
    "- `07_temporal_analysis`: multi-scene temporal summaries and time synchronisation",
    "- `08_exploratory_analysis`: regression, Moran's I, GWR diagnostics, and recommendation outputs"
  )
  writeLines(readme_lines, con = file.path(config$output_dir, "README.md"))

  list(
    inspected = inspected,
    harmonized = harmonized,
    index_inventory = index_inventory,
    sensor_outputs = sensor_outputs,
    point_outputs = point_outputs,
    polygon_outputs = polygon_outputs,
    surface_outputs = surface_outputs,
    temporal_outputs = temporal_outputs,
    exploratory_outputs = exploratory_outputs,
    run_summary = run_summary
  )
}
