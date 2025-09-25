# HAB Project – Meta-analysis & Environmental Drivers

This repository contains data, scripts, and notebooks for analyzing **harmful algal blooms (HABs)** and their association with **marine mammal mortalities**, with a focus on the Caspian Sea.  
The project combines **meta-analysis of toxin concentrations** with **environmental drivers** (temperature, wind, chlorophyll, river discharge).

---

## Repository Structure

```
data/
└── processed/
    ├── final_meta_table.csv       # Master table of toxin concentration studies
    ├── mortalities.csv            # Pinniped mortality records
    ├── volga_discharge.xlsx       # Volga River discharge data
    ├── wind_caspian_2010_24.csv   # Wind forcing time series
    ├── chla_caspian_2021Q2.nc     # MODIS chlorophyll-a (example)
    └── era5_caspian_2025_04_06.nc # ERA5 reanalysis (SST, winds, etc.)

notebooks/
├── caspian_mortalities.ipynb      # Mortality events analysis
├── chla_analysis.ipynb            # Chlorophyll-a satellite data
├── proportions.ipynb              # Species/matrix proportions
└── temp_analysis.ipynb            # Temperature anomalies and trends

r_scripts/
├── heterogeneity.R                # Meta-analysis heterogeneity statistics
├── forest_plots.R                 # Forest/range plots for DA and STX
└── bias_tests.R                   # Egger’s test + funnel plots

results/
└── funnel_plots_svg/              # Funnel plot outputs (SVG)
```

---

## Workflows

### 1. Meta-analysis (R)
- **`heterogeneity.R`**  
  Computes pooled effect sizes, τ², I², and prediction intervals.  
  Outputs:
  - `heterogeneity_by_toxin.csv`
  - `heterogeneity_by_matrix.csv`
  - `heterogeneity_by_region.csv`

- **`forest_plots.R`**  
  Produces range-style forest plots of toxin concentrations.  
  Outputs:
  - `forest_DA_methods.svg`
  - `forest_STX_methods.svg`

- **`bias_tests.R`**  
  Runs Egger’s regression test for small-study effects.  
  Generates funnel plots for each toxin.  
  Outputs:
  - `eggers_by_toxin.csv`
  - `eggers_by_matrix.csv`
  - `eggers_by_region.csv`
  - `funnel_<TOXIN>.svg` (in `results/funnel_plots_svg/`)

### 2. Environmental analysis (Python, Jupyter notebooks)
- **`chla_analysis.ipynb`** – Chlorophyll-a variability from MODIS L3  
- **`temp_analysis.ipynb`** – SST anomalies, ΔSST offshore–coastal  
- **`caspian_mortalities.ipynb`** – Aligns mortality records with environmental drivers  
- **`proportions.ipynb`** – Visualization of sample/species proportions  

---

## Requirements

### R
- `tidyverse`
- `stringi`
- `stringr`
- `metafor`
- `patchwork`
- `ggtext`

Install in R:

```r
install.packages(c("tidyverse", "stringi", "stringr", "metafor", "patchwork", "ggtext"))
```

### Python
- `numpy`
- `pandas`
- `matplotlib`
- `xarray`
- `netCDF4`
- `cartopy` (optional, for mapping)

Install via conda:

```bash
conda install -c conda-forge numpy pandas matplotlib xarray netCDF4 cartopy
```

---

## Usage

### Running R scripts

From inside `r_scripts/`:

```bash
Rscript heterogeneity.R
Rscript forest_plots.R
Rscript bias_tests.R
```

Results will be written into `../results/`.

### Running notebooks

From the repo root:

```bash
jupyter lab
```

and open any notebook under `notebooks/`.

---

## Outputs

- **CSV summaries** (meta-analysis results, Egger’s tests)  
- **Plots** (SVG/PNG): forest plots, funnel plots  
- **Notebooks** with exploratory analysis  

---
