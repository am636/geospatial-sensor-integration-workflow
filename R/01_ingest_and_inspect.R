# File discovery, metadata inspection, and Sentinel parsing.

interview_classify_file <- function(file_path) {
  ext <- tolower(tools::file_ext(file_path))
  if (ext %in% c("tif", "tiff", "jp2")) return("raster")
  if (ext %in% c("gpkg", "shp", "geojson")) return("vector")
  if (ext %in% c("csv")) return("table")
  if (ext %in% c("xml", "aux", "dbf", "shx", "prj", "cpg", "qpj")) return("sidecar")
  "other"
}


interview_list_files <- function(data_dir) {
  files <- list.files(data_dir, recursive = TRUE, full.names = TRUE, all.files = FALSE, no.. = TRUE)
  files <- files[!dir.exists(files)]
  files <- sort(normalizePath(files, winslash = "/", mustWork = TRUE))
  data.frame(
    file_path = files,
    file_name = basename(files),
    file_type = vapply(files, interview_classify_file, character(1)),
    stringsAsFactors = FALSE
  )
}


interview_parse_band_code <- function(file_name) {
  base_name <- tools::file_path_sans_ext(basename(file_name))
  band_match <- regmatches(base_name, regexpr("B[0-9]{1,2}|SCL", base_name, ignore.case = TRUE))
  if (length(band_match) == 0 || is.na(band_match)) {
    return(NA_character_)
  }
  band_match <- toupper(band_match)
  if (band_match == "SCL") {
    return("SCL")
  }
  digits <- gsub("[^0-9]", "", band_match)
  sprintf("B%02d", as.integer(digits))
}


interview_parse_scene_datetime <- function(file_name, fallback_datetime = NULL) {
  base_name <- tools::file_path_sans_ext(basename(file_name))
  ts_match <- regmatches(base_name, regexpr("[0-9]{8}T[0-9]{6}", base_name))
  if (length(ts_match) == 1 && !is.na(ts_match) && nzchar(ts_match)) {
    return(as.POSIXct(ts_match, format = "%Y%m%dT%H%M%S", tz = "UTC"))
  }

  date_match <- regmatches(base_name, regexpr("[0-9]{8}", base_name))
  if (length(date_match) == 1 && !is.na(date_match) && nzchar(date_match)) {
    return(as.POSIXct(paste0(date_match, " 00:00:00"), format = "%Y%m%d %H:%M:%S", tz = "UTC"))
  }

  if (!is.null(fallback_datetime) && nzchar(fallback_datetime)) {
    return(as.POSIXct(fallback_datetime, tz = "UTC"))
  }

  as.POSIXct(NA, tz = "UTC")
}


interview_parse_tile_id <- function(file_name) {
  base_name <- tools::file_path_sans_ext(basename(file_name))
  tile_match <- regmatches(base_name, regexpr("T[0-9]{2}[A-Z]{3}", base_name, ignore.case = TRUE))
  if (length(tile_match) == 1 && !is.na(tile_match) && nzchar(tile_match)) {
    return(toupper(tile_match))
  }
  "UNKNOWN_TILE"
}


interview_raster_metadata <- function(file_path, fallback_datetime = NULL) {
  x <- terra::rast(file_path)
  ext <- terra::ext(x)
  scene_datetime <- interview_parse_scene_datetime(file_path, fallback_datetime = fallback_datetime)
  data.frame(
    file_name = basename(file_path),
    file_path = normalizePath(file_path, winslash = "/", mustWork = TRUE),
    band_code = interview_parse_band_code(file_path),
    tile_id = interview_parse_tile_id(file_path),
    scene_datetime = scene_datetime,
    scene_date = if (!is.na(scene_datetime)) format(scene_datetime, "%Y-%m-%d", tz = "UTC") else NA_character_,
    crs = terra::crs(x),
    resolution_x = terra::res(x)[1],
    resolution_y = terra::res(x)[2],
    xmin = ext[1],
    xmax = ext[2],
    ymin = ext[3],
    ymax = ext[4],
    rows = terra::nrow(x),
    cols = terra::ncol(x),
    layers = terra::nlyr(x),
    stringsAsFactors = FALSE
  )
}


interview_vector_metadata <- function(file_path) {
  x <- sf::st_read(file_path, quiet = TRUE)
  data.frame(
    file_name = basename(file_path),
    file_path = normalizePath(file_path, winslash = "/", mustWork = TRUE),
    geometry_type = paste(unique(as.character(sf::st_geometry_type(x))), collapse = ", "),
    crs = sf::st_crs(x)$input %||% NA_character_,
    features = nrow(x),
    stringsAsFactors = FALSE
  )
}


interview_table_metadata <- function(file_path) {
  dat <- utils::read.csv(file_path, stringsAsFactors = FALSE)
  dat <- interview_sanitize_data_frame(dat)
  empty_cols <- vapply(dat, function(x) all(is.na(x) | trimws(as.character(x)) == ""), logical(1))
  data.frame(
    file_name = basename(file_path),
    file_path = normalizePath(file_path, winslash = "/", mustWork = TRUE),
    rows = nrow(dat),
    cols = ncol(dat),
    column_names = paste(names(dat), collapse = ", "),
    empty_columns = paste(names(dat)[empty_cols], collapse = ", "),
    stringsAsFactors = FALSE
  )
}


interview_select_sensor_table <- function(table_inventory, sensors_config) {
  if (nrow(table_inventory) == 0) {
    return(NULL)
  }

  scored <- lapply(seq_len(nrow(table_inventory)), function(i) {
    dat <- utils::read.csv(table_inventory$file_path[i], stringsAsFactors = FALSE)
    dat <- interview_sanitize_data_frame(dat)
    names_lc <- tolower(names(dat))
    score <- 0
    score <- score + ifelse(tolower(sensors_config$x_col) %in% names_lc, 3, 0)
    score <- score + ifelse(tolower(sensors_config$y_col) %in% names_lc, 3, 0)
    score <- score + ifelse(tolower(sensors_config$id_col) %in% names_lc, 2, 0)
    score <- score + sum(tolower(names(sensors_config$value_aliases)) %in% names_lc)
    data.frame(
      file_path = table_inventory$file_path[i],
      file_name = table_inventory$file_name[i],
      score = score,
      stringsAsFactors = FALSE
    )
  })

  scored <- do.call(rbind, scored)
  scored <- scored[order(-scored$score, scored$file_name), , drop = FALSE]
  scored$file_path[1]
}


run_auto_ingest_inspection <- function(data_dir, output_dir, fallback_datetime = NULL) {
  interview_require_packages()
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  inventory <- interview_list_files(data_dir)
  raster_files <- inventory[inventory$file_type == "raster", , drop = FALSE]
  vector_files <- inventory[inventory$file_type == "vector", , drop = FALSE]
  table_files <- inventory[inventory$file_type == "table", , drop = FALSE]

  raster_inventory <- if (nrow(raster_files) > 0) {
    do.call(rbind, lapply(raster_files$file_path, interview_raster_metadata, fallback_datetime = fallback_datetime))
  } else {
    data.frame()
  }

  vector_inventory <- if (nrow(vector_files) > 0) {
    do.call(rbind, lapply(vector_files$file_path, interview_vector_metadata))
  } else {
    data.frame()
  }

  table_inventory <- if (nrow(table_files) > 0) {
    do.call(rbind, lapply(table_files$file_path, interview_table_metadata))
  } else {
    data.frame()
  }

  utils::write.csv(inventory, file.path(output_dir, "file_inventory.csv"), row.names = FALSE)
  utils::write.csv(raster_inventory, file.path(output_dir, "raster_inventory.csv"), row.names = FALSE)
  utils::write.csv(vector_inventory, file.path(output_dir, "vector_inventory.csv"), row.names = FALSE)
  utils::write.csv(table_inventory, file.path(output_dir, "table_inventory.csv"), row.names = FALSE)

  report_lines <- c(
    "# Inspection Summary",
    "",
    paste0("- Raster files: ", nrow(raster_files)),
    paste0("- Vector files: ", nrow(vector_files)),
    paste0("- Table files: ", nrow(table_files)),
    ""
  )

  if (nrow(raster_inventory) > 0) {
    report_lines <- c(report_lines, "## Raster Inventory", "", interview_write_markdown_table(raster_inventory))
  }
  if (nrow(vector_inventory) > 0) {
    report_lines <- c(report_lines, "", "## Vector Inventory", "", interview_write_markdown_table(vector_inventory))
  }
  if (nrow(table_inventory) > 0) {
    report_lines <- c(report_lines, "", "## Table Inventory", "", interview_write_markdown_table(table_inventory))
  }
  writeLines(report_lines, con = file.path(output_dir, "inspection_summary.md"))

  list(
    file_inventory = inventory,
    raster_inventory = raster_inventory,
    vector_inventory = vector_inventory,
    table_inventory = table_inventory
  )
}
