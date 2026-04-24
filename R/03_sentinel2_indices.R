# Sentinel-2 index logic.

interview_safe_ratio <- function(num, den) {
  out <- num / den
  out[den == 0] <- NA
  out[!is.finite(out)] <- NA
  out
}


interview_supported_indices <- function() {
  list(
    ndvi = list(
      label = "NDVI",
      required_bands = c("B08", "B04"),
      formula = function(b) interview_safe_ratio(b$B08 - b$B04, b$B08 + b$B04),
      palette = c("#7f3b08", "#f7f7f7", "#1a9850")
    ),
    evi = list(
      label = "EVI",
      required_bands = c("B08", "B04", "B02"),
      formula = function(b) 2.5 * ((b$B08 - b$B04) / (b$B08 + 6 * b$B04 - 7.5 * b$B02 + 1)),
      palette = c("#8c510a", "#f6e8c3", "#01665e")
    ),
    savi = list(
      label = "SAVI",
      required_bands = c("B08", "B04"),
      formula = function(b) 1.5 * ((b$B08 - b$B04) / (b$B08 + b$B04 + 0.5)),
      palette = c("#8c510a", "#f6e8c3", "#01665e")
    ),
    gndvi = list(
      label = "GNDVI",
      required_bands = c("B08", "B03"),
      formula = function(b) interview_safe_ratio(b$B08 - b$B03, b$B08 + b$B03),
      palette = c("#7f3b08", "#f7f7f7", "#1a9850")
    ),
    ndwi = list(
      label = "NDWI",
      required_bands = c("B03", "B08"),
      formula = function(b) interview_safe_ratio(b$B03 - b$B08, b$B03 + b$B08),
      palette = c("#543005", "#f5f5f5", "#0571b0")
    )
  )
}


interview_build_scene_table <- function(raster_inventory) {
  x <- raster_inventory
  x$scene_datetime <- as.POSIXct(x$scene_datetime, tz = "UTC")
  x$scene_key <- paste(
    x$tile_id %||% "UNKNOWN_TILE",
    ifelse(is.na(x$scene_datetime), "UNKNOWN_TIME", format(x$scene_datetime, "%Y%m%dT%H%M%S", tz = "UTC")),
    sep = "__"
  )
  x[order(x$scene_datetime, x$band_code, x$file_name), , drop = FALSE]
}


interview_apply_cloud_mask <- function(index_raster, scene_rows, cloud_mask_config) {
  if (!isTRUE(cloud_mask_config$enabled)) {
    return(index_raster)
  }
  scl_row <- scene_rows[scene_rows$band_code == "SCL", , drop = FALSE]
  if (nrow(scl_row) == 0) {
    return(index_raster)
  }
  scl <- terra::rast(scl_row$analysis_raster_path[1])
  keep <- cloud_mask_config$keep_scl_classes %||% c(4, 5, 6)
  mask <- scl
  vals <- terra::values(mask, mat = FALSE)
  vals[!vals %in% keep] <- NA
  terra::values(mask) <- vals
  terra::mask(index_raster, mask)
}


interview_compute_index_maps <- function(index_inventory, output_dir, plot_config) {
  dir.create(file.path(output_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
  supported <- interview_supported_indices()
  for (i in seq_len(nrow(index_inventory))) {
    r <- terra::rast(index_inventory$raster_path[i])
    vals <- terra::as.data.frame(r, xy = TRUE, na.rm = TRUE)
    names(vals)[3] <- "value"
    spec <- supported[[index_inventory$index_name[i]]]
    p <- ggplot2::ggplot(vals, ggplot2::aes(x = x, y = y, fill = value)) +
      ggplot2::geom_raster() +
      ggplot2::scale_fill_gradientn(colors = spec$palette, name = spec$label) +
      ggplot2::coord_equal() +
      interview_theme() +
      ggplot2::labs(
        title = paste0(spec$label, " map"),
        subtitle = paste0("Scene: ", index_inventory$scene_label[i])
      )
    ggplot2::ggsave(
      filename = file.path(output_dir, "figures", paste0(index_inventory$index_name[i], "_", index_inventory$scene_tag[i], "_map.png")),
      plot = p,
      width = plot_config$width,
      height = plot_config$height,
      dpi = plot_config$dpi
    )
  }
}


run_sentinel2_index_module <- function(prepared_inventory, config, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_dir, "rasters"), recursive = TRUE, showWarnings = FALSE)
  scenes <- interview_build_scene_table(prepared_inventory)
  if (!is.null(config$selected_sentinel_layers) && length(config$selected_sentinel_layers) > 0) {
    required_bands <- unique(unlist(lapply(interview_supported_indices()[config$selected_indices], `[[`, "required_bands")))
    keep_bands <- unique(c(config$selected_sentinel_layers, required_bands))
    scenes <- scenes[scenes$band_code %in% keep_bands, , drop = FALSE]
  }
  supported <- interview_supported_indices()
  selected <- intersect(config$selected_indices, names(supported))
  if (length(selected) == 0) {
    stop("No supported indices selected.")
  }

  index_rows <- list()
  scene_keys <- unique(scenes$scene_key)
  for (scene_key in scene_keys) {
    scene_rows <- scenes[scenes$scene_key == scene_key, , drop = FALSE]
    band_layers <- setNames(lapply(scene_rows$analysis_raster_path, terra::rast), scene_rows$band_code)
    scene_datetime <- scene_rows$scene_datetime[1]
    scene_tag <- if (is.na(scene_datetime)) "unknown_time" else format(scene_datetime, "%Y%m%dT%H%M%S", tz = "UTC")
    scene_label <- if (is.na(scene_datetime)) "Unknown time" else format(scene_datetime, "%Y-%m-%d %H:%M:%S", tz = "UTC")

    for (index_name in selected) {
      spec <- supported[[index_name]]
      if (!all(spec$required_bands %in% names(band_layers))) {
        next
      }
      bands <- band_layers[spec$required_bands]
      index_raster <- spec$formula(bands)
      names(index_raster) <- index_name
      index_raster <- interview_apply_cloud_mask(index_raster, scene_rows, config$cloud_mask)

      out_name <- paste0(index_name, "_", scene_tag, ".tif")
      out_path <- file.path(output_dir, "rasters", out_name)
      terra::writeRaster(index_raster, out_path, overwrite = TRUE)

      index_rows[[length(index_rows) + 1]] <- data.frame(
        index_name = index_name,
        index_label = spec$label,
        scene_key = scene_key,
        scene_datetime = scene_datetime,
        scene_date = as.Date(scene_datetime),
        scene_tag = scene_tag,
        scene_label = scene_label,
        raster_path = normalizePath(out_path, winslash = "/", mustWork = TRUE),
        tile_id = scene_rows$tile_id[1],
        stringsAsFactors = FALSE
      )
    }
  }

  index_inventory <- do.call(rbind, index_rows)
  utils::write.csv(index_inventory, file.path(output_dir, "index_inventory.csv"), row.names = FALSE)
  interview_compute_index_maps(index_inventory, output_dir, config$plots)
  index_inventory
}
