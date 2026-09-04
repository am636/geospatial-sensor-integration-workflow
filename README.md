# Geospatial sensor and Sentinel-2 integration in R

This repository contains an R workflow I use to explore how field or sensor measurements can be linked with Sentinel-2 imagery and related spatial layers.

The workflow covers:

- input and metadata checks;
- raster harmonisation when grids differ;
- Sentinel-2 index calculation (NDVI, EVI, SAVI, GNDVI and NDWI);
- cleaning and spatially referencing sensor tables;
- extraction of raster values at monitoring locations;
- optional polygon summaries and interpolation;
- simple temporal and exploratory spatial analysis.

The repository is intended as a reusable research utility rather than a finished study. Input data are not included because the workflow was written to work with local project datasets.

## Structure

- `run_main.R` — entry point and user settings
- `R/00_bootstrap.R` — shared setup and helper functions
- `R/01_ingest_and_inspect.R` — input discovery and metadata checks
- `R/02_harmonization.R` — raster alignment
- `R/03_sentinel2_indices.R` — Sentinel-2 indices
- `R/04_sensors.R` — sensor-table preparation and CRS handling
- `R/05_integration_and_surfaces.R` — point/polygon integration and optional surfaces
- `R/06_temporal_and_plots.R` — temporal summaries and plotting
- `R/07_exploratory_analysis.R` — regression, Moran's I and optional GWR checks
- `R/08_workflow.R` — workflow sequence

## Running

1. Put local inputs in `data/inbox/`.
2. Open `run_main.R`.
3. Set the small configuration block for the dataset being used.
4. Run `run_main.R`.

The default run is conservative: cloud masking and interpolation are off, rasters are only harmonised when needed, and missing scene dates remain missing rather than being inferred.

## Notes

This repository demonstrates the data-integration logic and analytical structure. Any ecological or environmental interpretation depends on the quality, scale and provenance of the input data used in a particular study.

**Author:** Ali Moayedi  
University of St Andrews
