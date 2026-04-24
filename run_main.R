# Main entry point.

resolve_project_dir <- function() {
  script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  candidate_dirs <- character(0)

  if (length(script_arg) > 0) {
    script_path <- sub("^--file=", "", script_arg[1])
    candidate_dirs <- c(candidate_dirs, dirname(normalizePath(script_path, winslash = "/", mustWork = TRUE)))
  }

  wd <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  candidate_dirs <- c(candidate_dirs, wd)
  candidate_dirs <- unique(candidate_dirs)

  for (dir_path in candidate_dirs) {
    if (file.exists(file.path(dir_path, "R", "00_bootstrap.R"))) {
      return(normalizePath(dir_path, winslash = "/", mustWork = TRUE))
    }
  }

  stop("Could not locate the project folder containing R/00_bootstrap.R")
}

PROJECT_DIR <- resolve_project_dir()
source(file.path(PROJECT_DIR, "R", "00_bootstrap.R"))
interview_source_modules(PROJECT_DIR)

config <- default_interview_config(PROJECT_DIR)

# Small user-facing configuration block.
config$author_name <- ""
config$author_email <- ""
config$selected_indices <- c("ndvi", "evi", "savi", "gndvi", "ndwi")
config$fallback_scene_datetime <- NULL
config$cloud_mask$enabled <- FALSE
config$sensors$manual_source_crs <- NULL
config$harmonization$force <- FALSE
config$surfaces$enabled <- FALSE
config$exploratory$run_gwr <- TRUE

results <- run_interview_workflow(config)
print(results$run_summary)
