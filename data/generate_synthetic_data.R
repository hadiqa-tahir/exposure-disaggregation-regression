# =============================================================================
# generate_synthetic_data.R
#
# Creates a fully synthetic dataset for the Albania DHS methods paper pipeline.
# No original DHS survey data is used in the output.
#
# Strategy:
#   - Sample random valid pixels from the Albania NDVI raster as cluster centres
#     (guarantees coordinates fall on Albanian land)
#   - Randomly assign 1,114 synthetic individuals across 322 clusters
#   - Extract NDVI within 3,000m buffer for each cluster (terra::extract)
#   - Simulate ARI outcomes at observed prevalence (~2.1%)
#   - Output matches the column structure expected by the model pipeline
#
# Output: Datasets/df_albania_1114.parquet
#   Columns: HOUSE, ARI, NDVI, DISTANCES, X_cluster, Y_cluster, x_ndvi, y_ndvi
#
# Real DHS data cannot be shared publicly. To reproduce the paper's results,
# apply for the Albania 2017-18 DHS data at https://dhsprogram.com
# =============================================================================

library(arrow)
library(dplyr)
library(terra)
library(sf)
library(here)

set.seed(20170911)
buffer_radius <- 3000    # meters (Study uses 3000m)
# =============================================================================
# 1. Load and prepare raster
# =============================================================================

albania_ndvi     <- rast(here::here("Data", "albania_ndvi.tif"))
albania_boundary <- vect(here::here("Data", "df_albania_boundary_shp",
                                    "alb_admbnda_adm0_2019c.shp"))

# Match CRS, crop, mask
albania_boundary <- project(albania_boundary, albania_ndvi)
albania_ndvi     <- crop(albania_ndvi, albania_boundary)
albania_ndvi     <- mask(albania_ndvi, albania_boundary)
plot(albania_ndvi)

cat("Raster loaded and masked.\n")

# =============================================================================
# 2. Sample 322 cluster centroids from valid raster pixels
# =============================================================================

# Get coordinates of all non-NA pixels
all_coords <- which(!is.na(values(albania_ndvi)))
all_coords  <- xyFromCell(albania_ndvi, all_coords)   # matrix: x (easting), y (northing)

# Exclude pixels within 3,000m of the raster extent edge so every sampled
# cluster centroid has a complete buffer within the raster
ext        <- ext(albania_ndvi)
all_coords <- all_coords[
  all_coords[, "x"] >= (ext$xmin + buffer_radius) &
    all_coords[, "x"] <= (ext$xmax - buffer_radius) &
    all_coords[, "y"] >= (ext$ymin + buffer_radius) &
    all_coords[, "y"] <= (ext$ymax - buffer_radius), ]

cat(sprintf("Interior pixels available for sampling: %d\n", nrow(all_coords)))

n_clusters  <- 322
cluster_idx <- sample(nrow(all_coords), n_clusters, replace = FALSE)
cluster_xy  <- all_coords[cluster_idx, ]   # n_clusters x 2, UTM 32634

cat(sprintf("%d cluster centroids sampled from raster.\n", n_clusters))

# =============================================================================
# 3. Assign 1,114 individuals across clusters (random cluster sizes)
# =============================================================================

n_individuals <- 1114

# Random assignment — naturally produces varying cluster sizes
individual_cluster <- sample(seq_len(n_clusters), n_individuals, replace = TRUE)

individuals <- data.frame(
  HOUSE     = seq_len(n_individuals),
  DHSCLUST  = individual_cluster,
  X_cluster = cluster_xy[individual_cluster, "x"],
  Y_cluster = cluster_xy[individual_cluster, "y"]
)

cat(sprintf("Cluster size range: %d to %d\n",
            min(table(individuals$DHSCLUST)),
            max(table(individuals$DHSCLUST))))

rm(cluster_idx, individual_cluster, all_coords, cluster_xy, ext)

# =============================================================================
# 6. Simulate ARI at observed prevalence (~2.1%)
# =============================================================================

ari_draw   <- as.numeric(rbinom(n_individuals, 1, prob = 0.021))
ari <- data.frame(HOUSE = seq_len(n_individuals), ARI = as.numeric(as.logical(ari_draw)))
individuals <- individuals %>% left_join(ari, by = "HOUSE")
cat(sprintf("Simulated ARI: %d / %d (%.1f%%)\n",
            sum(ari_draw), n_individuals, 100 * mean(ari_draw)))

rm(ari, ari_draw)

# =============================================================================
# 4. Extract NDVI pixels within 3,000m buffer for each cluster
# =============================================================================

albania <- st_as_sf(individuals, coords = c("X_cluster", "Y_cluster"), crs = 32634, remove = FALSE) 
plot(albania$geometry, add = TRUE)

# Create buffer_radiusm buffer and extract NDVI data 
albania <- st_buffer(albania, buffer_radius)   
ndvi <- terra::extract(albania_ndvi, albania, xy = TRUE)
names(ndvi) <- c("HOUSE", "NDVI", "x_ndvi", "y_ndvi")

cat("Extraction complete.\n")

# remove columns 
albania <- st_drop_geometry(albania)

# merge to create the dataframe
df <- merge(albania, ndvi, by = "HOUSE")
df <- df %>% 
  filter(!is.na(NDVI))


# =============================================================================
# 5. compute distances / weights
# =============================================================================

df <- df %>%
  mutate(
    DISTANCES = sqrt((X_cluster - x_ndvi)^2 + (Y_cluster - y_ndvi)^2),
    DISTANCES = if_else(DISTANCES == 0, 10, DISTANCES),
    WEIGHT    = (1 / DISTANCES) / ave(1 / DISTANCES, HOUSE, FUN = sum)
  )

cat(sprintf("Total pixel rows: %d\n", nrow(df)))
#df <- st_as_sf(df, coords = c("X_cluster", "Y_cluster"), crs = 32634, remove = FALSE) 


# =============================================================================
# 7. Save
# =============================================================================

write_parquet(df, here::here("data", "df_1114_synthetic.parquet"))
cat("Saved.")
