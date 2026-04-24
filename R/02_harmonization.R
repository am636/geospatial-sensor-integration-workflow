# Conditional harmonization helpers.

interview_choose_reference_raster <- function(raster_inventory, manual_reference = NULL) {
  if (!is.null(manual_reference) && manual_reference %in% raster_inventory$file_name) {
    return(raster_inventory$file_path[match(manual_reference, raster_inventory$file_name)])
  }

  preferred <- raster_inventory[raster_inventory$band_code == "B04", , drop = FALSE]
  if (nrow(preferred) > 0) {
    preferred <- preferred[order(preferred$scene_datetime, preferred$file_name), , drop = FALSE]
    return(preferred$file_path[1])
  }

  raster_inventory$file_path[1]
}


interview_rasters_need_harmonization <- function(raster_paths, reference_path, manual_target_crs = NULL) {
  ref <- terra::rast(reference_path)
  for (path in raster_paths) {
    x <- terra::rast(path)
    same_geom <- tryCatch(terra::compareGeom(ref, x, stopOnError = FALSE), error = function(e) FALSE)
    same_crs <- identical(terra::crs(ref), terra::crs(x))
    target_ok <- is.null(manual_target_crs) || identical(terra::crs(ref), manual_target_crs) || identical(terra::crs(x), manual_target_crs)
    if (!isTRUE(same_geom) || !isTRUE(same_crs) || !isTRUE(target_ok)) {
      return(TRUE)
    }
  }
  FALSE
}


interview_harmonize_single_raster <- function(file_path, reference_raster, output_path, method = "bilinear", target_crs = NULL) {
  x <- terra::rast(file_path)

  if (!is.null(target_crs) && !identical(terra::crs(x), target_crs)) {
    x <- terra::project(x, target_crs, method = method)
  }

  same_geom <- tryCatch(terra::compareGeom(reference_raster, x, stopOnError = FALSE), error = function(e) FALSE)
  if (!isTRUE(same_geom)) {
    if (!identical(terra::crs(reference_raster), terra::crs(x))) {
      x <- terra::project(x, terra::crs(reference_raster), method = method)
    }
    x <- terra::resample(x, reference_raster, method = method)
  }

  terra::writeRaster(x, output_path, overwrite = TRUE)
  normalizePath(output_path, winslash = "/", mustWork = TRUE)
}


run_conditional_harmonization <- function(raster_inventory, output_dir, harmonization_config) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_dir, "rasters"), recursive = TRUE, showWarnings = FALSE)

  reference_path <- interview_choose_reference_raster(
    raster_inventory = raster_inventory,
    manual_reference = harmonization_config$manual_reference
  )
  reference_raster <- terra::rast(reference_path)

  needs_harmonization <- isTRUE(harmonization_config$force) ||
    interview_rasters_need_harmonization(
      raster_paths = raster_inventory$file_path,
      reference_path = reference_path,
      manual_target_crs = harmonization_config$manual_target_crs
    )

  out_inventory <- raster_inventory
  out_inventory$reference_file <- basename(reference_path)
  out_inventory$analysis_raster_path <- out_inventory$file_path
  out_inventory$harmonized <- FALSE

  if (needs_harmonization) {
    for (i in seq_len(nrow(out_inventory))) {
      method <- if (toupper(out_inventory$band_code[i]) == "SCL") {
        harmonization_config$categorical_method
      } else {
        harmonization_config$continuous_method
      }
      out_path <- file.path(output_dir, "rasters", basename(out_inventory$file_path[i]))
      out_inventory$analysis_raster_path[i] <- interview_harmonize_single_raster(
        file_path = out_inventory$file_path[i],
        reference_raster = reference_raster,
        output_path = out_path,
        method = method,
        target_crs = harmonization_config$manual_target_crs
      )
      out_inventory$harmonized[i] <- TRUE
    }
  }

  utils::write.csv(out_inventory, file.path(output_dir, "harmonization_inventory.csv"), row.names = FALSE)

  summary_df <- data.frame(
    reference_file = basename(reference_path),
    harmonization_run = needs_harmonization,
    raster_count = nrow(out_inventory),
    stringsAsFactors = FALSE
  )
  utils::write.csv(summary_df, file.path(output_dir, "harmonization_summary.csv"), row.names = FALSE)

  list(
    inventory = out_inventory,
    reference_path = normalizePath(reference_path, winslash = "/", mustWork = TRUE),
    harmonization_run = needs_harmonization
  )
}
