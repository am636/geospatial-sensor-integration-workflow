# Geospatial Sensor Integration Workflow

Author: Ali Moayedi · Contact: am656@st-andrews.ac.uk

R workflow for integrating Sentinel-2 raster data with point and polygon datasets, deriving linked index products such as `NDVI`, `EVI`, `SAVI`, `GNDVI`, and `NDWI`, and building connected analysis tables for raster, point, polygon, and multi-scene time-aware workflows.

## Project layout

- `run_main.R` - main entry point
- `R/` - workflow modules
- `data/inbox/` - raw inputs
- `output/` - generated outputs

## Modules

- `00_bootstrap.R` - setup, defaults, shared helpers
- `01_ingest_and_inspect.R` - file discovery and metadata checks
- `02_harmonization.R` - raster harmonization
- `03_sentinel2_indices.R` - multi-scene Sentinel-2 index generation
- `04_sensors.R` - sensor cleaning, CRS inference, and attribute handling
- `05_integration_and_surfaces.R` - point and polygon integration
- `06_temporal_and_plots.R` - time handling, synchronisation, and plotting helpers
- `07_exploratory_analysis.R` - regression, Moran's I, and GWR diagnostics
- `08_workflow.R` - workflow orchestration

## Running the project

1. Place your input files in `data/inbox/`
2. Open `run_main.R`
3. Adjust the small config block at the top if needed
4. Run the script

Default behaviour:

- cloud masking off
- raster harmonization only when required
- interpolation surfaces off
- time-aware scene handling for multi-date Sentinel inputs
- scene-to-sensor time synchronisation where timestamps are available
- point and polygon attribute tables carried into the integrated outputs

If raster filenames do not include scene dates, the workflow keeps time fields as unknown rather than fabricating them.
