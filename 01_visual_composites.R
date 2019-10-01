
# *************************************************************************************************
# Script to:
#     * Calculate 3-month composites (for visual classification)
#     * First method: via medians
#     * Second method: via overlay (scene with less clouds on top) 
# 
# *************************************************************************************************


# Load libraries
library(tidyverse)
library(lubridate)
library(stars)
library(tictoc)


dir <- "/media/cdobler/PSEUDBOMBAX/processed_data/calakmul_landsat_analysis/"

# Assign quatrimesters to scenes
tibble(file = list.files(str_c(dir, "00_cropped_masked_bands/")),
       
       date = file %>% 
         str_split("_", simplify = T) %>% .[,2] %>% 
         as_date(),
       
       year = year(date),
       month = month(date),
       
       comp = case_when(month %in% c(12, 1, 2) ~ 1L,
                        month %in% c(3, 4, 5) ~ 2L,
                        month %in% c(6, 7, 8) ~ 3L,
                        month %in% c(9, 10, 11) ~ 4L)) %>% 
  
  mutate(year = ifelse(month == 12, year+1, year)) %>% 
  group_by(year, comp) %>% 
  mutate(rn = group_indices()) -> date_tbl


# Method 1: pixel-wise medians ********************************************************************

# Parallelize
library(furrr)
plan(multisession, workers = 6)

Sys.time()
tic()

# Loop through quatrimester-years
date_tbl %>% 
  pull(rn) %>% 
  unique() %>% 
  
  future_map(function(i){
    
    # Vector of available files within the quatrimester-year
    date_tbl %>% 
      filter(rn == i) %>% 
      pull(file) %>% 
      str_c(dir, "00_cropped_masked_bands/", .) -> files
    
    # Date in quatrimester-year format
    date_tbl %>% 
      filter(rn == i) %>% 
      summarize() %>% 
      mutate(date = str_c(year, "_", comp)) %>% 
      pull(date) -> date
    
    # Is there only one file within the quatrimester-year? Export as-is.
    if(length(files) == 1){
      
      read_stars(files, RasterIO = list(bands = 2:7)) %>% 
        
        as("Raster") %>%
        raster::writeRaster(str_c(dir, "01_three_month_composites/comp3m_",date, ".tif"),
                            datatype = "INT2S")
      
    # Are there more files? Calculate median.
    } else {
      
      # Loop through bands
      2:7 %>% 
        map(function(b){
          
          # Loop through files within the quatrimester-year
          seq_along(files) %>% 
            map(function(f){
              
              files[f] %>% 
                read_stars(RasterIO = list(bands = b))
              
            }) %>% 
            do.call(c, .) %>% 
            merge() %>% 
            
            # Calculate median
            st_apply(c(1,2), median, na.rm = T)
          
        }) %>%
        
        do.call(c, .) %>% 
        merge() %>% 
        
        # Export
        as("Raster") %>%
        raster::writeRaster(str_c(dir, "01_three_month_composites/comp3m_",date, ".tif"),
                            datatype = "INT2S")
      
    }
    
  })

toc()


# Method 2: pixel-wise overlay ********************************************************************

# Parallelize
# plan(multisession, workers = 6)

Sys.time()
tic()

# Loop through quatrimester-years
date_tbl %>% 
  pull(rn) %>% 
  unique() %>% 
  
  walk(function(i){
    
    # Vector of available files within the quatrimester-year
    date_tbl %>% 
      filter(rn == i) %>% 
      pull(file) %>% 
      str_c(dir, "00_cropped_masked_bands/", .) -> files
    
    # Date in quatrimester-year format
    date_tbl %>% 
      filter(rn == i) %>% 
      summarize() %>% 
      mutate(date = str_c(year, "_", comp)) %>% 
      pull(date) -> date
    
    # Is there only one file within the quatrimester-year? Export as-is.
    if(length(files) == 1){
      
      read_stars(files, RasterIO = list(bands = 2:7)) %>% 
        
        as("Raster") %>%
        raster::writeRaster(str_c(dir, "01_three_month_overlays/comp3m_", date, ".tif"),
                            datatype = "INT2S")
      
    # Are there more files? Overlay.
    } else {
      
      # Obtain info on gaps and arrange scenes from less gaps to more
      files %>% 
        read_stars(RasterIO = list(bands = 1)) %>%
        as_tibble() %>%
        gather(-(1:2), key = date, value = val) %>% 
        filter(is.na(val)) %>% 
        group_by(date) %>% 
        count() %>% 
        arrange(n) %>% 
        pull(date) -> files_clouds
      
      # Loop through bands
      2:7 %>% 
        map(function(b){
          
          str_c(dir, "00_cropped_masked_bands/", files_clouds) %>% 
            read_stars(RasterIO = list(bands = b)) %>% 
            
            merge() %>% 
            as("Raster") %>% 
            
            # Merge in oder given
            raster::merge()
            
        }) %>% 
        
        raster::stack() %>% 
        
        # Export
        raster::writeRaster(str_c(dir, "01_three_month_overlays/comp3m_", date, ".tif"),
                            datatype = "INT2S")
      
    }
    
  })

toc()