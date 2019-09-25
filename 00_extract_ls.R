
# *************************************************************************************************
# Script to:
#     * decompress Landsat images 
#     * crop to extent
#     * save as multiband tif's in disk (including qa layer)
#
# *************************************************************************************************


# Load libraries
library(tidyverse)
library(lubridate)
library(raster)


# Set directory of input and output files
input_dir <- "/media/cdobler/PSEUDBOMBAX/raw_data/calakmul_landsat_time_series_surf_ref/"
output_dir <- "/home/cdobler/Documents/calakmul_landsat_analysis/"

# Create object to crop (encompasses 4 ejidos)
e <- extent(232305, 274605, 2006535, 2085825)

# Loop to extract scenes from zip
for(i in seq_along(list.files(input_dir))){
  
  file <- list.files(input_dir)[i]
  
  # Change name: sr_year_month_day_sensor
  name <- file %>%
  {
    str_c(str_sub(., 11, 14),"_", # year
          str_sub(., 15, 16),"_", # month
          str_sub(., 17, 18),"_", # day
          str_sub(., 4, 4))       # sensor
  } %>% {str_c("sr_", .)}
  
  # Print progress
  print(str_c("Processing ", i, "/", length(list.files(input_dir)), " --- ", name))
  
  # Extract bands and QA layer
  # for Landsat 4-7:
  if(str_split(name, "_", simplify = T)[5] != "8"){
    
    # Untar only selected bands
    untar(str_c(input_dir, file), list = T) %>% 
      .[grepl("_band1|_band2|_band3|_band4|_band5|_band7|pixel", .)] %>% 
      untar(str_c(input_dir, file), 
            files = ., 
            exdir = str_c(output_dir, "temp_processing/"))
    
    # Crop and stack
    ls_stack <- list.files(str_c(output_dir, "temp_processing/"), full.names = T) %>% 
      .[!grepl("aux", .)] %>% 
      stack() %>% 
      crop(e)
    
  # for Landsat 8:
  } else {
    
    # Untar
    untar(str_c(input_dir, file), list = T) %>% 
      .[grepl("_band2|_band3|_band4|_band5|_band6|_band7|pixel", .)] %>% 
      untar(str_c(input_dir, file), 
            files = ., 
            exdir = str_c(output_dir, "temp_processing/"))
    
    # Crop and stack
    ls_stack <- list.files(str_c(output_dir, "temp_processing/"), full.names = T) %>% 
      .[!grepl("aux", .)] %>% 
      stack() %>% 
      crop(e)
    
  }
  
  # Save to disk
  writeRaster(ls_stack, 
              str_c(output_dir, "00_cropped_scenes/", name, ".tif"), 
              options="INTERLEAVE=BAND")
  
  # NOTE:
  # Saved file        Original (L 4-7)  Original (L 8)
  # ***************   ****************  **************
  # Band 1 = QA      
  # Band 2 = Blue            1                 2
  # Band 3 = Green           2                 3
  # Band 4 = Red             3                 4
  # Band 5 = NIR             4                 5
  # Band 6 = SWIR 1          5                 6
  # Band 7 = SWIR 2          7                 7
  
  # Remove raw tifs
  str_c(output_dir, "temp_processing/") %>%  
    list.files(full.names = T) %>% 
    file.remove()
  
} # end of for loop
