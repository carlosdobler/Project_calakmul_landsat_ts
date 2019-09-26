
# SCRIPT TO:
#   - Mask all landsat scenes
#
# *************************************************************************************************

library(tidyverse)
library(raster)

# External directory (input/output)
ext_dir <- "/home/cdobler/Documents/calakmul_landsat_analysis/"

# List of landsat stacks
list_of_files <- grep(list.files(str_c(ext_dir, "00_cropped_scenes")), 
                      pattern = "aux", 
                      inv = T, 
                      value = T)

# Loop
for(i in seq_along(list_of_files)) {
  
  name <- str_split(list_of_files[i], "[.]", simplify = T)[1]
  
  # Print progress
  print(str_c("Processing ", i, "/", length(list_of_files), " --- ", name))
  
  # Load raster
  r <- stack(str_c(ext_dir, "00_cropped_scenes/", list_of_files[i]), bands = 2:7)
  
  # Load mask
  m <- raster(str_c(ext_dir, "01_qa_masks/qa_", list_of_files[i]))
  
  # Aplly mask to raster
  r_m <- overlay(r, m, fun = function(a, b){
    a[b == 0] <- NA
    return(a)
  })
  
  # Save to disk
  writeRaster(r_m, 
              str_c(ext_dir, "02_masked_scenes/rm_", name, ".tif"), 
              options="INTERLEAVE=BAND")
  
}
