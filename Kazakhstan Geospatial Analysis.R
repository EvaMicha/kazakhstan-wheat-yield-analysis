# ==============================================================================
# PROJECT:  Geospatial Analysis of Wheat Yield Potential — Kazakhstan
# AUTHOR:   Makhambet Yeszhanov
# DATE:     May 2026
# COURSE:   Theory and Applications of Statistics
# TYPE:     Academic Project — independently researched, coded, and analyzed
#
# DESCRIPTION:
#   Analyzes the spatial distribution of wheat agro-climatic yield potential
#   across 220 administrative districts (ADM2) of Kazakhstan using real spatial
#   data from FAO GAEZ v4, GADM administrative boundaries, and SRTM elevation.
#   Produces statistical summaries, exploratory visualizations, spatial maps,
#   a fitted variogram, and an Ordinary Kriging interpolation surface.
#
# HYPOTHESIS:
#   Northern districts (steppe belt) show the highest wheat yield potential,
#   reflecting Kazakhstan's known agro-ecological gradient. Higher yield
#   potential corresponds to conditions that support rural agricultural income.
#
# DATA SOURCES:
#   [1] FAO GAEZ v4 — whea200b_yld.tif (manual download required, see README)
#   [2] GADM — Kazakhstan ADM1/ADM2 boundaries (auto-downloaded in script)
#   [3] SRTM — Elevation 30 arc-second (auto-downloaded in script)
#
# REQUIREMENTS:
#   - R version 4.0 or higher
#   - whea200b_yld.tif must be placed in the same folder as this script
#   - Internet connection required for first run (GADM + SRTM download)
#   - Linux users: install system libraries before running
#       yay -S udunits && sudo pacman -S gdal geos proj   # Arch/CachyOS
#       sudo apt-get install libudunits2-dev libgdal-dev  # Ubuntu/Debian
#
# LICENSE:
#   MIT License — free to use, modify, and distribute with attribution
# ==============================================================================

if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  tidyverse,
  sf,
  terra,
  geodata,
  gstat,
  stars,
  viridis,
  scales,
  skimr,
  janitor,
  patchwork,
  ggnewscale
)

# ------------------------------------------------------------------------------
# IMPORTANT: Change this path to the folder where you placed whea200b_yld.tif
# Example Windows: setwd("C:/Users/YourName/Documents/KAZ_analysis/")
# Example Mac:     setwd("/Users/YourName/Documents/KAZ_analysis/")
# Example Linux:   setwd("/home/yourname/Documents/KAZ_analysis/")
# ------------------------------------------------------------------------------
setwd()

# ==============================================================================
# PART 1: LOAD DATA
# ==============================================================================

# 1.1 Kazakhstan Administrative Boundaries from GADM
# ADM1 = 17 Oblasts (regions) | ADM2 = 220 Districts
cat("Loading GADM boundaries...\n")
kaz_adm1 <- gadm("KAZ", level = 1, path = tempdir()) %>% st_as_sf()
kaz_adm2 <- gadm("KAZ", level = 2, path = tempdir()) %>% st_as_sf()

# 1.2 Wheat Yield Raster from GAEZ v4
# whea200b_yld.tif: Agro-climatic potential yield for wheat
# "200b" = baseline period ~2000 | Rain-fed | High input | Units: kg/ha
cat("Loading wheat yield raster...\n")
wheat_rast <- rast("whea200b_yld.tif")

# 1.3 Elevation Raster from SRTM (same function used in class for Japan)
cat("Loading elevation raster...\n")
elev_rast <- elevation_30s("KAZ", path = tempdir())


# ==============================================================================
# PART 2: CRS EXPLORATION AND TRANSFORMATION TO WGS84
# ==============================================================================

cat("\n=== CRS EXPLORATION ===\n")

# 2.1 Inspect original CRS of each dataset
cat("\n--- Wheat Raster (original CRS) ---\n")
print(crs(wheat_rast, describe = TRUE))

cat("\n--- ADM2 Boundaries (original CRS) ---\n")
print(st_crs(kaz_adm2))

cat("\n--- Elevation Raster (original CRS) ---\n")
print(crs(elev_rast, describe = TRUE))

# 2.2 What is WGS84?
# WGS84 (EPSG:4326) is the global geographic coordinate system used by GPS.
# Coordinates are in degrees of latitude and longitude.
# It is the standard reference system for global datasets.
# GAEZ and GADM data are typically already in WGS84.
# We verify and explicitly transform to be sure — this is good practice.

cat("\n--- Transforming all datasets to WGS84 (EPSG:4326) ---\n")

kaz_adm1_wgs <- st_transform(kaz_adm1, 4326)
kaz_adm2_wgs <- st_transform(kaz_adm2, 4326)
wheat_rast_wgs <- project(wheat_rast, "EPSG:4326")  # terra::project for rasters
elev_rast_wgs  <- project(elev_rast,  "EPSG:4326")

cat("All datasets now in WGS84 (EPSG:4326)\n")
cat("Note: For variograms we will later reproject to UTM Zone 42N (EPSG:32642)\n")
cat("      because variograms require distances in METERS, not degrees.\n")


# ==============================================================================
# PART 3: CLIP RASTERS TO KAZAKHSTAN AND EXTRACT DATA AT DISTRICT CENTROIDS
# ==============================================================================

# 3.1 Create a single Kazakhstan boundary polygon for masking
kaz_boundary <- st_union(kaz_adm2_wgs) %>% vect()

# 3.2 Clip global rasters to Kazakhstan extent + mask (remove outside pixels)
cat("\nClipping rasters to Kazakhstan...\n")
wheat_kaz <- crop(wheat_rast_wgs, kaz_boundary) %>% mask(kaz_boundary)
elev_kaz  <- crop(elev_rast_wgs,  kaz_boundary) %>% mask(kaz_boundary)

# 3.3 Rename raster layers for clarity
names(wheat_kaz) <- "wheat_yield"
names(elev_kaz)  <- "elevation"

# 3.4 Compute district centroid coordinates (the point that represents each district)
kaz_centroids <- kaz_adm2_wgs %>%
  st_centroid() %>%
  mutate(
    lon = st_coordinates(.)[, 1],
    lat = st_coordinates(.)[, 2]
  )

# 3.5 Extract raster values at each district centroid
# This gives each of the 220 districts a wheat yield and elevation value
wheat_vals <- terra::extract(wheat_kaz, vect(kaz_centroids))
elev_vals  <- terra::extract(elev_kaz,  vect(kaz_centroids))

# 3.6 Assemble final analysis dataframe
kaz_data <- kaz_adm2_wgs %>%
  st_drop_geometry() %>%                        # Remove geometry for tabular work
  mutate(
    wheat_yield = wheat_vals[, 2],              # Column 1 = ID, Column 2 = value
    elevation   = elev_vals[, 2],
    lon         = kaz_centroids$lon,
    lat         = kaz_centroids$lat
  ) %>%
  filter(!is.na(wheat_yield) & wheat_yield > 0) %>%
  clean_names()                                 # Standardize column names

cat(paste("\nDistricts with valid wheat yield data:", nrow(kaz_data), "\n"))


# ==============================================================================
# PART 4: STATISTICAL SUMMARY
# ==============================================================================

cat("\n=== FULL STATISTICAL SUMMARY (skimr) ===\n")
kaz_data %>%
  select(wheat_yield, elevation, lat, lon) %>%
  skim()

# Confidence interval calculation (95%)
n       <- nrow(kaz_data)
mean_y  <- mean(kaz_data$wheat_yield)
sd_y    <- sd(kaz_data$wheat_yield)
se_y    <- sd_y / sqrt(n)
ci_low  <- mean_y - 1.96 * se_y
ci_high <- mean_y + 1.96 * se_y

cat("\n=== WHEAT YIELD — KEY INDICATORS ===\n")
cat(sprintf("  N districts    : %d\n",   n))
cat(sprintf("  Mean           : %.1f kg/ha\n", mean_y))
cat(sprintf("  Median         : %.1f kg/ha\n", median(kaz_data$wheat_yield)))
cat(sprintf("  Std. Deviation : %.1f kg/ha\n", sd_y))
cat(sprintf("  Min            : %.1f kg/ha\n", min(kaz_data$wheat_yield)))
cat(sprintf("  Max            : %.1f kg/ha\n", max(kaz_data$wheat_yield)))
cat(sprintf("  95%% CI         : [%.1f , %.1f] kg/ha\n", ci_low, ci_high))

# Summary by Oblast (regional breakdown)
cat("\n=== SUMMARY BY OBLAST ===\n")
oblast_summary <- kaz_data %>%
  group_by(name_1) %>%
  summarise(
    n_districts  = n(),
    mean_yield   = round(mean(wheat_yield), 1),
    sd_yield     = round(sd(wheat_yield), 1),
    min_yield    = round(min(wheat_yield), 1),
    max_yield    = round(max(wheat_yield), 1),
    mean_elev    = round(mean(elevation, na.rm = TRUE), 0)
  ) %>%
  arrange(desc(mean_yield))

print(oblast_summary)


# ==============================================================================
# PART 5: NON-SPATIAL GRAPHICS
# ==============================================================================

# --- 5.1 DISTRIBUTION: Histogram of Wheat Yield ---
p1 <- ggplot(kaz_data, aes(x = wheat_yield)) +
  geom_histogram(fill = "wheat3", color = "black", bins = 30, alpha = 0.85) +
  geom_vline(xintercept = mean_y,          color = "red",   linetype = "dashed", size = 1) +
  geom_vline(xintercept = ci_low,          color = "blue",  linetype = "dotted", size = 0.8) +
  geom_vline(xintercept = ci_high,         color = "blue",  linetype = "dotted", size = 0.8) +
  annotate("text", x = mean_y + 50, y = 35,
           label = paste0("Mean = ", round(mean_y, 0), " kg/ha"),
           color = "red", hjust = 0, size = 3.5) +
  theme_minimal() +
  labs(title = "1. Distribution of Wheat Yield Potential Across Districts",
       subtitle = "Red = Mean | Blue dotted = 95% Confidence Interval | GAEZ v4",
       x = "Yield Potential (kg/ha)", y = "Number of Districts")

print(p1)

# --- 5.2 RANKING: Top 20 Districts by Wheat Yield ---
p2 <- kaz_data %>%
  arrange(desc(wheat_yield)) %>%
  head(20) %>%
  ggplot(aes(x = reorder(name_2, wheat_yield), y = wheat_yield, fill = name_1)) +
  geom_col(color = "black", alpha = 0.9) +
  coord_flip() +
  scale_fill_viridis_d(option = "turbo") +
  scale_y_continuous(labels = comma) +
  theme_minimal() +
  labs(title = "2. Top 20 Districts by Wheat Yield Potential",
       x = "District (ADM2)", y = "Yield Potential (kg/ha)", fill = "Oblast")

print(p2)

# --- 5.3 BOXPLOT: Yield distribution by Oblast ---
p3 <- ggplot(kaz_data, aes(x = reorder(name_1, wheat_yield, FUN = median), y = wheat_yield)) +
  geom_boxplot(fill = "lightgreen", outlier.alpha = 0.4, outlier.size = 1) +
  coord_flip() +
  scale_y_continuous(labels = comma) +
  theme_minimal() +
  labs(title = "3. Wheat Yield Distribution by Oblast",
       subtitle = "Ordered by median yield",
       x = "Oblast", y = "Yield Potential (kg/ha)")

print(p3)

# --- 5.4 SCATTER: Elevation vs. Yield ---
# The Bid-Rent equivalent for agriculture: terrain affects productivity
p4 <- ggplot(kaz_data, aes(x = elevation, y = wheat_yield, color = name_1)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", color = "black", se = TRUE, linetype = "dashed") +
  scale_color_viridis_d() +
  theme_light() +
  labs(title = "4. Elevation vs. Wheat Yield Potential",
       subtitle = "East Kazakhstan mountain foothills show high yield due to precipitation",
       x = "Elevation (m)", y = "Yield Potential (kg/ha)", color = "Oblast")

print(p4)

# --- 5.5 SCATTER: Latitude vs. Yield (North-South gradient) ---
# Kazakhstan's wheat belt is in the north (Akmola, Kostanay, North Kazakhstan)
p5 <- ggplot(kaz_data, aes(x = lat, y = wheat_yield)) +
  geom_point(alpha = 0.5, color = "darkgreen", size = 2) +
  geom_smooth(method = "loess", color = "red", se = TRUE) +
  theme_minimal() +
  labs(title = "5. Latitude vs. Wheat Yield (North-South Gradient)",
       subtitle = "Lowest yields at mid-latitudes (central desert)",
       x = "Latitude (degrees North)", y = "Yield Potential (kg/ha)")

print(p5)



# ==============================================================================
# PART 6: MAPS
# ==============================================================================

# 6.1 Convert rasters to dataframes for ggplot
# Aggregate elevation first — SRTM at full resolution has ~7.5M pixels for Kazakhstan
# fact=10 reduces it 10x → ~75,000 pixels — still visually accurate, much faster
elev_kaz_agg  <- aggregate(elev_kaz,  fact = 10, fun = mean)
wheat_kaz_agg <- aggregate(wheat_kaz, fact = 2,  fun = mean)

wheat_df <- as.data.frame(wheat_kaz_agg, xy = TRUE) %>%
  rename(wheat_yield = 3) %>%
  filter(!is.na(wheat_yield) & wheat_yield > 0)

elev_df <- as.data.frame(elev_kaz_agg, xy = TRUE) %>%
  rename(elevation = 3) %>%
  filter(!is.na(elevation))

# 6.2 MAP 1: Wheat Yield Raster + Oblast Boundaries
map1 <- ggplot() +
  geom_tile(data = wheat_df, aes(x = x, y = y, fill = wheat_yield)) +
  geom_sf(data = kaz_adm2_wgs, fill = NA, color = "gray40", size = 0.1) +
  geom_sf(data = kaz_adm1_wgs, fill = NA, color = "black",  size = 0.5) +
  scale_fill_viridis_c(option = "viridis", name = "Yield\n(kg/ha)", labels = comma) +
  theme_void() +
  labs(title = "6. Wheat Agro-climatic Potential Yield — Kazakhstan",
       subtitle = "GAEZ v4 | Rain-fed, High Input | ADM2 district borders shown")

print(map1)

# 6.3 MAP 2: Elevation
map2 <- ggplot() +
  geom_tile(data = elev_df, aes(x = x, y = y, fill = elevation)) +
  geom_sf(data = kaz_adm1_wgs, fill = NA, color = "white", size = 0.4) +
  scale_fill_gradient(low = "white", high = "black", name = "Elev (m)") +
  theme_void() +
  labs(title = "7. Elevation — Kazakhstan",
       subtitle = "SRTM 30s | Flat steppe in north, Tian Shan mountains in southeast")

print(map2)

# 6.4 MAP 3: Combined — Elevation (raster background) + Wheat Yield (district points)
# This directly mirrors the class example: Japan topology + land prices
centroids_joined <- kaz_centroids %>%
  left_join(kaz_data %>% select(GID_2 = gid_2, wheat_yield), by = "GID_2")

map3 <- ggplot() +
  geom_tile(data = elev_df, aes(x = x, y = y, fill = elevation)) +
  scale_fill_gradient(low = "white", high = "black", guide = "none") +
  new_scale_fill() +  # Requires 'ggnewscale' — see note below
  geom_sf(data = centroids_joined, aes(color = wheat_yield), size = 2, alpha = 0.85) +
  scale_color_viridis_c(option = "turbo", name = "Yield (kg/ha)", labels = comma) +
  geom_sf(data = kaz_adm1_wgs, fill = NA, color = "red", size = 0.4) +
  theme_void() +
  labs(title = "8. Terrain & Wheat Yield Potential — Kazakhstan",
       subtitle = "Background: Elevation (white=low, black=high) | Dots: District yield potential")

# Note: map3 uses two color scales. Install 'ggnewscale' if not already installed:
# install.packages("ggnewscale"); library(ggnewscale)

print(map3)



# ==============================================================================
# PART 7: VARIOGRAM ANALYSIS
# ==============================================================================

# --- What is a Variogram? ---
# A variogram measures how similar values are as a function of distance between points.
# If points far apart are as different as nearby points → no spatial structure (pure nugget).
# If nearby points are more similar and differences grow with distance → spatial autocorrelation.
# This is key for understanding whether wheat yield clusters spatially (it should).

# 7.1 Prepare spatial point data for variogram
# We MUST reproject to a metric CRS — distances must be in METERS, not degrees
# EPSG:32642 = WGS84 / UTM Zone 42N — covers central Kazakhstan
kaz_sf <- kaz_data %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  st_transform(32642)

# 7.2 Log-transform wheat yield
# Yield values are right-skewed (few very high values)
# Variogram/Kriging assumes approximate normality → log-transform corrects this
# Same approach used in class for land prices (ln_price)
kaz_sf <- kaz_sf %>%
  mutate(ln_yield = log(wheat_yield)) %>%
  filter(is.finite(ln_yield))

cat("\n=== VARIOGRAM INPUT DATA ===\n")
cat(sprintf("Points for variogram: %d\n", nrow(kaz_sf)))
cat(sprintf("Mean ln_yield: %.3f\n", mean(kaz_sf$ln_yield)))
cat(sprintf("Variance of ln_yield: %.3f\n", var(kaz_sf$ln_yield)))

# 7.3 Calculate Experimental Variogram
# 'ln_yield ~ 1' = Ordinary Kriging (constant but unknown mean)
var_exp <- variogram(ln_yield ~ 1, kaz_sf)

# 7.4 Fit Theoretical Model
# Exponential model: approaches the sill gradually, never fully reaches it.
# Chosen over Spherical because the experimental variogram was noisy (only 92 points).
# psill  = semivariance plateau — initial guess near observed variance
# range  = distance at which spatial correlation becomes negligible
# nugget = variance at distance zero (micro-scale variation)
# fit.method = 6 uses Ordinary Least Squares — more robust than default
#              weighted method when data is irregular or sparse

var_model <- fit.variogram(
  var_exp,
  vgm(psill = 20, model = "Exp", range = 300000, nugget = 1),
  fit.method = 6
)

cat("\n=== FITTED VARIOGRAM MODEL PARAMETERS ===\n")
print(var_model)
cat("\nInterpretation:\n")
cat("  Nugget : micro-scale variance (below district resolution)\n")
cat("  Sill   : total variance (nugget + psill)\n")
cat("  Range  : distance beyond which districts are spatially uncorrelated\n")

plot(var_exp, var_model,
     main  = "Variogram: Log Wheat Yield Potential (Kazakhstan, 92 Districts)",
     xlab  = "Distance h (meters)",
     ylab  = expression(paste("Semivariance ", gamma, "(h)")),
     col   = "darkgreen",
     pch   = 16)


# ==============================================================================
# PART 8: ORDINARY KRIGING INTERPOLATION
# ==============================================================================

# Kriging predicts yield values at unsampled locations using the variogram structure.
# We use it to produce a continuous surface of predicted yield across Kazakhstan.

# 8.1 Create prediction grid in UTM Zone 42N
kaz_bbox_utm <- st_bbox(st_transform(kaz_adm2_wgs, 32642))
kaz_grid     <- st_as_stars(kaz_bbox_utm, dx = 25000, dy = 25000)  # 25km grid

# 8.2 Run Ordinary Kriging
cat("\nRunning Ordinary Kriging... (this may take 1-2 minutes)\n")
kaz_kriging <- krige(
  formula   = ln_yield ~ 1,
  locations = kaz_sf,
  newdata   = kaz_grid,
  model     = var_model,
  nmax      = 30   # Use 30 nearest neighbors per prediction point (speeds up computation)
)
kaz_mask <- st_transform(st_union(kaz_adm2_wgs), 32642) %>% vect()
kaz_kriging_masked <- rast(kaz_kriging) %>% mask(kaz_mask)
# 8.3 Plot: Predicted Yield Surface
plot(kaz_kriging_masked["var1.pred"],
     main = "Ordinary Kriging: Predicted Log Wheat Yield — Kazakhstan",
     col  = viridis(100, option = "viridis"),
     axes = TRUE)

# 8.4 Plot: Estimation Variance (Uncertainty)
# Where we have no data points → high variance → low confidence in prediction
# This directly visualizes where more field surveys would be most valuable
plot(kaz_kriging_masked["var1.var"],
     main = "Ordinary Kriging: Estimation Variance (Uncertainty)",
     col  = viridis(100, option = "mako"),
     axes = TRUE)


# ==============================================================================
# SUMMARY OF FINDINGS
# ==============================================================================

# 1. DATA:
#    - 92 of 220 Kazakhstan districts had valid wheat yield data
#    - 128 districts excluded — desert south and mountain peaks where
#      rain-fed wheat cultivation is not viable (raster returns NA)
#    - Source: GAEZ v4 (FAO), rain-fed, high input, baseline ~2000
#    - Combined with SRTM elevation raster

# 2. STATISTICAL PATTERNS:
#    - National mean yield potential: ~2,463 kg/ha (SD: ~1,000 kg/ha)
#    - Wide distribution — Kazakhstan is highly heterogeneous agriculturally
#    - North Kazakhstan, Akmola, Qostanay oblasts show highest median yields
#    - Positive correlation between elevation and yield — counterintuitive,
#      explained by East Kazakhstan mountain foothills receiving higher
#      orographic precipitation than flat central desert

# 3. GEOGRAPHIC PATTERNS:
#    - Clear north-south-east structure visible in maps
#    - Northern steppe belt: consistent moderate-to-high yields
#    - East Kazakhstan Altai foothills: highest yields in the country
#    - Central and southern desert: no wheat cultivation (white areas in map)
#    - Variogram range ~293km confirms spatial clustering at regional scale

# 4. HYPOTHESIS OUTCOME:
#    - Partially confirmed: northern steppe is a major wheat belt as expected
#    - Unexpected finding: East Kazakhstan mountain foothills outperform
#      the northern steppe due to precipitation advantage, not terrain flatness
#    - Key driver is moisture availability, whether from latitude-linked
#      precipitation in the north or orographic precipitation in the east

cat("\n=== Script complete. All plots and variogram fitted. ===\n")

# ==============================================================================
# SESSION INFO
# Printed for reproducibility — shows exact R and package versions used
# ==============================================================================
sessionInfo()
