#
# SCRIPT TO:
#     * Calculate NDMIs
#
# *************************************************************************************************

# Load libraries
library(raster)
library(tidyverse)
library(furrr)
library(tictoc)
#library(conflicted)

# External directory (input/output)
ext_dir <- "/home/cdobler/Documents/calakmul_landsat_analysis/"

# Get list of scenes
list_of_files <- str_c(ext_dir, "02_masked_scenes/") %>% 
  list.files() %>% 
  grep(pattern = "aux", inv = T, value = T)


# NDMI stack
plan(multicore, workers = 6)

tic()
ndmi_stack <- future_map(seq_along(list_of_files),
                         function(s) {
                           
                           sc <- str_c(ext_dir, "02_masked_scenes/", list_of_files[s]) %>% 
                             stack(bands = c(4, 5)) #   *** CHANGE BANDS ***
                           
                           ndmi <- overlay(sc, fun = function(x, y){
                             round(1000 * (x - y)/(x + y))
                           }
                           )
                           
                           date <- list_of_files[s] %>%
                             str_sub(end = -5) %>% 
                             str_split("_", simplify = T) %>% 
                             {str_c(.[3], "_", .[4], "_", .[5])}
                           
                           writeRaster(ndmi, 
                                       str_c(ext_dir,
                                             "04_ndmi/",  # *** CHANGE NAME ***
                                             "ndmi_",  # *** CHANGE NAME ***
                                             date, 
                                             ".tif"), 
                                       datatype = "INT2S")
                           
                         }
)
toc()




