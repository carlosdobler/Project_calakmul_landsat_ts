
# SCRIPT TO:
#   * Generate sample points based on a unsupervised classification of 3 dates (ndmi)
#
# *************************************************************************************************


# Load libraries
library(raster)
library(sf)
library(tidyverse)

# Input dir
in_dir <- "/home/cdobler/Documents/calakmul_landsat_analysis/"

# Load three ndmi scenes

# Nueva vida
# ndmi1 <- raster(str_c(in_dir, "04_ndmi/ndmi_1987_09_24.tif"))
# ndmi2 <- raster(str_c(in_dir, "04_ndmi/ndmi_2003_03_20.tif"))
# ndmi3 <- raster(str_c(in_dir, "04_ndmi/ndmi_2017_11_29.tif"))

# La Guadalupe
ndmi1 <- raster(str_c(in_dir, "04_ndmi/ndmi_1988_11_05.tif"))
ndmi2 <- raster(str_c(in_dir, "04_ndmi/ndmi_2003_03_20.tif"))
ndmi3 <- raster(str_c(in_dir, "04_ndmi/ndmi_2017_11_29.tif"))

# Import ejidos shape
ejido <- read_sf("/media/cdobler/Neobxbaumia/Research/Data_spatial/Shapefiles/03 SYPR shapefiles/Ejidos/Ejidos_select.shp") %>% 
  st_transform(projection(ndmi1)) %>%
  select("NOM_NA") %>% 
  slice(3) # ejido of interest    *** CHANGE ***

# Crop rasters
ndmi_stack <- stack(ndmi1, ndmi2, ndmi3) %>% 
  crop(ejido)

# Unsupervised classif
ndmi_tbl <- tabularaster::as_tibble(ndmi_stack) %>% 
  spread(key = dimindex, value = cellvalue)

set.seed(246)
clusters <- kmeans(na.omit(ndmi_tbl[ , 2:4]), iter.max = 500, centers = 10, nstart = 5, algorithm = "Lloyd")

# Raster map of clusters
clusters_map <- ndmi_stack[[1]] %>% raster()

clusters_map[] <- na.omit(ndmi_tbl) %>%
  {tibble(cellindex = .$cellindex, clus = clusters$cluster)} %>% 
  right_join(ndmi_tbl) %>% 
  arrange(cellindex) %>% 
  .$clus
  
# Polygonize
clusters_map_poly <- clusters_map %>% 
  rasterToPolygons() %>% 
  st_as_sf() %>% 
  group_by(layer) %>% 
  summarize()

# Sample points
sample_pts <- st_sample(clusters_map_poly, rep(5, 10)) %>% st_sf()
sample_pts$clus <- raster::extract(clusters_map, sample_pts)

st_write(sample_pts,
         dsn = "/media/cdobler/Neobxbaumia/Research/Project_calakmul_landsat_ts/data/training_pts_lg/pts.gpkg", # *** CHANGE ***
         driver = "GPKG")
  


