#
# SCRIPT TO:
#     * Calculate Tasseled Cap components
#
# *************************************************************************************************

# Load libraries
library(raster)
library(tidyverse)
library(furrr)
library(tictoc)
library(RStoolbox)
#library(conflicted)

# External directory (input/output)
ext_dir <- "/home/cdobler/Documents/calakmul_landsat_analysis/"

# Get list of scenes
list_of_files <- str_c(ext_dir, "02_masked_scenes/") %>% 
  list.files() %>% 
  grep(pattern = "aux", inv = T, value = T)


# TC stack
plan(multicore, workers = 6)

Sys.time()
tic()
tc_stack <- future_map(seq_along(list_of_files),
                         function(s) {
                           
                           name <- str_sub(list_of_files[s], start = 7, end = -5)
                           
                           sat <- str_split(name, "_", simplify = T)[4]
                           
                           sat <- case_when(
                             sat == "4" ~ "Landsat4TM",
                             sat == "5" ~ "Landsat5TM",
                             sat == "7" ~ "Landsat7ETM",
                             sat == "8" ~ "Landsat8OLI"
                           )
                           
                           sc <- str_c(ext_dir, "02_masked_scenes/", list_of_files[s]) %>% 
                             stack()
                           
                           tc <- tasseledCap(sc, sat = sat) %>% round()
                           
                           writeRaster(tc, 
                                       str_c(ext_dir,
                                             "04_tass_cap/",
                                             "tc_",
                                             name, 
                                             ".tif"), 
                                       datatype = "INT2S")
                           
                         }
)
toc()

