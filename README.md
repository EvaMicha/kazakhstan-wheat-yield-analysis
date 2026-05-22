# Geospatial Analysis of Wheat Yield Potential — Kazakhstan

**Author:** Yeszhanov Makhambet  
**Course:** Theory and Applications of Statistics  
**Type:** Academic Project — independently researched, coded, and analyzed  
**Language:** R  
**License:** MIT

---

## Overview

This project analyzes the spatial distribution of wheat agro-climatic yield potential across 220 administrative districts (ADM2) of Kazakhstan. Wheat is the primary cash crop for Kazakhstani rural farmers, making land productivity a key indicator of rural agricultural income potential.

The analysis combines three real-world datasets, produces statistical summaries, exploratory visualizations, spatial maps, a fitted variogram, and an Ordinary Kriging interpolation surface.

---

## Hypothesis

Northern districts — the steppe belt — show the highest wheat yield potential, reflecting Kazakhstan's known agro-ecological gradient. Higher yield potential corresponds to conditions that support rural agricultural income.

**Result:** Partially confirmed. The northern steppe is a major wheat belt as expected. However, East Kazakhstan's Altai mountain foothills showed the highest yields in the country, driven by orographic precipitation rather than flat terrain — an unexpected finding that contradicts the initial assumption that elevation negatively affects yield.

---

## Repository Contents

| File | Description |
|---|---|
| `Statistics_Quiz_1.R` | Complete R script — all analysis, plots, maps, variogram, and Kriging |
| `whea200b_yld.tif` | Wheat yield raster from FAO GAEZ v4 (see download instructions below) |
| `README.md` | This file |

---

## Data Sources

### `whea200b_yld.tif` — FAO GAEZ v4 (manual download required)
| Field | Value |
|---|---|
| Full name | Global Agro-Ecological Zones, version 4 |
| Source | https://gaez.fao.org → Theme 3 → Data Viewer |
| Variable | Agro-climatic potential yield — Wheat |
| Period | Baseline year ~2000 (`200b` = year 2000 baseline) |
| Water supply | Rain-fed (no irrigation assumed) |
| Input level | High (optimal fertilizer and pest management) |
| Unit | kg/ha (kilograms per hectare) |
| Resolution | ~9.25 km per pixel (5 arc-minutes) |
| Coverage | Global — clipped to Kazakhstan in the script |
| CRS | WGS84 (EPSG:4326) |

> **Note:** This is theoretical yield ceiling given climate conditions — not actual recorded farmer yields.

**How to download:**
1. Go to https://gaez.fao.org
2. Open Theme 3 — Agro-climatic Potential Yield → Data Viewer
3. Filter: Crop = Wheat, Water Supply = Rain-fed, Input Level = High
4. Click any pixel on the map → download link appears in the popup
5. Place `whea200b_yld.tif` in the same folder as the script

### GADM — Administrative Boundaries (auto-downloaded)
- Kazakhstan ADM1: 17 Oblasts
- Kazakhstan ADM2: 220 Districts
- Downloaded automatically via `geodata::gadm("KAZ", level = 1)` and `level = 2`

### SRTM — Elevation (auto-downloaded)
- 30 arc-second resolution (~1 km per pixel)
- Downloaded automatically via `geodata::elevation_30s("KAZ")`

---

## Requirements

### R Version
R 4.0 or higher

### System Libraries (Linux only)
Must be installed before R packages:

```bash
# Arch / CachyOS
yay -S udunits
sudo pacman -S gdal geos proj

# Ubuntu / Debian
sudo apt-get install libudunits2-dev libgdal-dev libgeos-dev libproj-dev
```

### R Packages
Installed automatically by the script via `pacman::p_load()`:

| Package | Purpose |
|---|---|
| `tidyverse` | Data wrangling and ggplot2 visualization |
| `sf` | Vector spatial data — polygons, points |
| `terra` | Raster data processing |
| `geodata` | Download GADM boundaries and SRTM elevation |
| `gstat` | Variogram fitting and Kriging interpolation |
| `stars` | Raster prediction grid for Kriging |
| `viridis` | Color scales for maps and plots |
| `scales` | Number formatting (comma, percent) |
| `skimr` | Detailed statistical summaries |
| `janitor` | Column name standardization |
| `patchwork` | Combining multiple ggplot2 plots |
| `ggnewscale` | Dual color scales in combined maps |

---

## How to Run

1. Download `whea200b_yld.tif` (see instructions above)
2. Place `whea200b_yld.tif` and `Statistics_Quiz_1.R` in the same folder
3. Open `Statistics_Quiz_1.R` and update the `setwd()` path at the top to match your folder
4. Run the script from top to bottom — **do not skip or rearrange sections**, each part depends on the previous
5. On first run, GADM boundaries and elevation data download automatically — internet connection required
6. Kriging in Part 8 takes approximately 1-2 minutes

---

## Script Structure

| Part | Content | Output |
|---|---|---|
| 1 | Load all three datasets | Console messages |
| 2 | CRS inspection and transformation to WGS84 | Console output |
| 3 | Clip rasters to Kazakhstan, extract values at district centroids | `kaz_data` dataframe (92 rows) |
| 4 | Statistical summary — mean, SD, CI, oblast breakdown | Console tables |
| 5 | Non-spatial graphics | 5 plots |
| 6 | Spatial maps | 3 maps |
| 7 | Experimental variogram + fitted exponential model | 1 variogram plot |
| 8 | Ordinary Kriging — prediction surface + uncertainty map | 2 raster plots |

---

## Plots Produced

| # | Type | Title |
|---|---|---|
| 1 | Histogram | Distribution of Wheat Yield Potential Across Districts |
| 2 | Bar chart | Top 20 Districts by Wheat Yield Potential |
| 3 | Boxplot | Wheat Yield Distribution by Oblast |
| 4 | Scatter | Elevation vs. Wheat Yield Potential |
| 5 | Scatter | Latitude vs. Wheat Yield (North-South Gradient) |
| 6 | Raster map | Wheat Agro-climatic Potential Yield — Kazakhstan |
| 7 | Raster map | Elevation — Kazakhstan |
| 8 | Combined map | Terrain & Wheat Yield Potential — Kazakhstan |
| 9 | Variogram | Log Wheat Yield Potential (92 Districts) |
| 10 | Kriging map | Predicted Log Wheat Yield — Kazakhstan |
| 11 | Kriging map | Estimation Variance (Uncertainty) |

---

## Key Findings

- Only **92 of 220 districts** have wheat yield data — the remaining 128 are desert or high mountain where rain-fed wheat is not viable
- National mean yield potential: **~2,463 kg/ha** (SD ~1,000 kg/ha)
- **Unexpected result:** elevation shows a positive correlation with yield — higher elevation = higher yield — because mountain foothills in East Kazakhstan receive significantly more precipitation than the flat central desert
- **Latitude pattern is U-shaped**, not linear — lowest yields occur at mid-latitudes (47-48°N), which correspond to the Betpak-Dala desert and dried Aral Sea basin
- **Variogram range ~293 km** — districts within this distance share spatially correlated yield conditions
- **Kriging uncertainty** is highest in the western desert and southern regions where no data points exist — predictions there are unreliable

---

## Known Limitations

- Yield values are agro-climatic **potential**, not actual recorded production
- Only 92 of 220 districts contributed to the variogram — sparse data reduces interpolation reliability
- A single UTM zone (EPSG:32642) was used for Kriging — Kazakhstan spans ~3,000 km east-west, introducing minor projection distortion at the western edge
- Kriging produces physically meaningless negative log-yield predictions in uninhabited desert areas far from any data point
