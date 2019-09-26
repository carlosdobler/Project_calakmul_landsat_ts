
# *************************************************************************************************
# Script to:
#     * Make composites based on median value of tri/cuatrimester (for visual classification)
# 
# *************************************************************************************************



# 01. Prep ----------------------------------------------------------------------------------------
library(tidyverse)
library(raster)
library(lubridate)
library(furrr)
library(tictoc)

# tri or quatri?
agg <- "tri"  # *** CHANGE ***

# External directory (input/output)
ext_dir <- "/home/cdobler/Documents/calakmul_landsat_analysis/"

# List of landsat stacks
list_of_files <- grep(list.files(str_c(ext_dir, "02_masked_scenes")), 
                      pattern = "aux", 
                      inv = T, 
                      value = T)

# Table with dates (modified for months of subsequent year)
if (agg == "quatri") {
  
  date_tbl <- tibble(file = list_of_files) %>% 
    mutate(year = as.integer(str_split(file, "_", simplify = T)[,3]),
           month = str_split(file, "_", simplify = T)[,4],
           comp = case_when(month %in% c("12", "01", "02") ~ "1st",
                              month %in% c("03", "04", "05") ~ "2nd",
                              month %in% c("06", "07", "08") ~ "3rd",
                              month %in% c("09", "10", "11") ~ "4th"),
           year = ifelse(month == "12", year+1L, year))
  
} else if (agg == "tri") {
  
  date_tbl <- tibble(file = list_of_files) %>% 
    mutate(year = as.integer(str_split(file, "_", simplify = T)[,3]),
           month = str_split(file, "_", simplify = T)[,4],
           comp = case_when(month %in% c("02", "03", "04", "05") ~ "1st",
                              month %in% c("06", "07", "08", "09") ~ "2nd",
                              month %in% c("10", "11", "12", "01") ~ "3rd"),
           year = ifelse(month == "01", year-1L, year))
  
}


# 02. MEAN of all pixels in a year/cuatrimester -------------------------------------------------
plan(multicore, workers = 6)

# Loop years
Sys.time()
tic()
future_map(unique(date_tbl$year), function(y){
  
  #Loop cuatrimesters
  sort(unique(date_tbl$comp)) %>% map(function(c) {
    
    # Print progress
    # print(str_c("Processing year ", y, ", cuatri ", c))
    
    t <- date_tbl %>% 
      filter(year == y & comp == c)
    
    # If only one date for year/cuatri
    if(nrow(t) == 1){
      
      r <- t %>% 
        .$file %>% 
        str_c(ext_dir, "02_masked_scenes/", .) %>% 
        stack()
      
      writeRaster(r, 
                  str_c(ext_dir, "03_composites/", agg, "/co_med_", y, "_", c, ".tif"), 
                  options="INTERLEAVE=BAND")
      
      # If more than one date for year/cuatri  
    } else if(nrow(t) > 1){
      
      # Loop bands
      1:6 %>% map(function(b){
        
        t %>% 
          .$file %>% 
          str_c(ext_dir, "02_masked_scenes/", .) %>% 
          map(function(i) raster(i, band = b)) %>% 
          stack() %>% 
          calc(median, na.rm = T)
        
      }) %>% 
        stack() -> r
      
      writeRaster(r, 
                  str_c(ext_dir, "03_composites/", agg, "/co_med_", y, "_", c, ".tif"), 
                  options="INTERLEAVE=BAND")
      
    } #  end of else 
    
    }) # end of cuatri loop
  }) # end of year loop
toc()

