
# *************************************************************************************************
# Script to:
#     * Decompress Landsat images 
#     * Crop to extent
#     * Calculate tasseled caps and NDMI
#     * Save as multiband tif's in disk (including qa layer)
#
# *************************************************************************************************


# Load libraries
library(tidyverse)
library(lubridate)
library(stars)

# Set directory of input and output files
input_dir <- "/media/cdobler/PSEUDBOMBAX/raw_data/calakmul_landsat_time_series_surf_ref/"
output_dir <- "/media/cdobler/PSEUDBOMBAX/processed_data/calakmul_landsat_analysis/"

# Loop through all files
seq_along(list.files(input_dir)) %>%
  map(function(sc){
    
    print(str_c("Processing scene ", sc, " / 397..."))
    
    file <- list.files(input_dir)[sc]
    sensor <- str_sub(file, 4, 4)
    date <- str_sub(file, 11, 18)
    
    # Untar
    untar(str_c(input_dir, file), list = T) %>% 
      str_subset("band|pixel") %>%
      {if(sensor == "8") str_subset(., "band1", negate = T) else .} %>%
      untar(str_c(input_dir, file), 
            files = ., 
            exdir = str_c(output_dir, "temp_processing/"))
    
    # Import and crop
    str_c(output_dir, "temp_processing/") %>%
      list.files(full.names = T) %>%
      read_stars() %>% 
      filter(y < 2091625,
             y > 1983335,
             x < 274605,
             x > 232305) %>% 
      setNames(c("qa", "b", "g", "r", "nir", "swir1", "swir2")) %>%
      
      # Create mask
      mutate(qa = case_when(sensor != "8" & qa == 66 ~ 1L,
                            sensor == "8" & qa %in% c(322, 386, 834, 898, 1346) ~ 1L,
                            TRUE ~ NA_integer_)) -> miau
    
    # Extract band
    miau %>% 
      select(qa) -> qa
    
    # Merge into one attribute
    miau %>% 
      merge() -> bnds
    
    # Apply mask
    bnds[is.na(qa)] <- NA
    
    # Calculate additional layers
    bnds %>%
      split("X1") %>% # Unfold
      
      mutate(tc_b = round((b*0.2043) + (g*0.4158) + (r*0.5524) + (nir*0.5741) + (swir1*0.3124) + (swir2*0.2303)),
             tc_g = round((b*-0.1603) + (g*0.2819) + (r*-0.4934) + (nir*0.7940) + (swir1*-0.0002) + (swir2*-0.1446)),
             tc_w = round((b*0.0315) + (g*0.2021) + (r*0.3102) + (nir*0.1594) + (swir1*-0.6806) + (swir2*-0.6109)),
             
             di = round(tc_b - (tc_g + tc_w)),
             
             ndmi = round((nir - swir1)/(nir + swir2) * 1000)) %>% 
      
      merge() %>% 
      
      # Convert to raster (smaller file size)
      as("Raster") %>% 
      
      # Export
      raster::writeRaster(str_c(output_dir, "00_cropped_masked_bands/img0_",date, "_", sensor, ".tif"),
                          datatype = "INT2S",
                          overwrite=TRUE)
    
      # write_stars(str_c(output_dir, "00_cropped_masked_bands/img0_",date, "_", sensor, ".tif"),
      #             options="INTERLEAVE=BAND",
      #             type="Int16")
    
    # Delete temporary files
    str_c(output_dir, "temp_processing/") %>%
      list.files(full.names = T) %>%
      file.remove()
    
  })


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
  