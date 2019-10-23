#
# *************************************************************************************************
# Script to:
#     * Map land cover using model
#     * SCRIPT TO BE RUN IN A VIRTUAL MACHINE OF 64 CORES
#
# *************************************************************************************************


# Load libraries
library(tidyverse)
library(lubridate)
library(stars)
library(bfast)
library(furrr)
library(tictoc)
library(imputeTS)
library(randomForest)


# Input/output dir
dir <- "/media/cdobler/PSEUDBOMBAX/processed_data/calakmul_landsat_analysis/"

# Load model
model_randfor_full <- readRDS(glue::glue(here::here("out_training//model_randfor_full.RDS")))

# Create stack of NBR layers
str_c(dir, "00_cropped_masked_bands") %>% 
  list.files(full.names = T) -> all_files_list

all_files_list %>% 
  #.[1:10] %>% # REMOVE!
  str_split("/", simplify = T) %>% 
  .[,8] %>% 
  str_split("_", simplify = T) %>% 
  .[, 2] %>% 
  as_date() -> dates

tic()
all_files_list %>% 
  #.[1:10] %>% # REMOVE!
  read_stars(RasterIO = list(bands = 14,
                             nXSize = 100, # 14 + 14 cells surrounding
                             nYSize = 100)
             ) %>% 
  merge() %>% 
  st_set_dimensions(3, values = dates) %>% 
  st_set_dimensions(names = c("x", "y", "date")) -> nbr_stack
toc() # ~ 6 min

# Tibble of coordinates
nbr_stack %>% 
  slice(date, 1) %>% 
  as_tibble() %>% 
  select(-X) %>% 
  mutate(id_cell = row_number()) -> coord_tbl


# MACRO LOOP *****
plan(multicore, workers = 60)

Sys.time()
tic("MAIN")

seq_vect <- seq(1, nrow(coord_tbl), 10000)
tbb <- vector("list", length(seq_vect))

# Outer loop: saves in list each 10,000 runs
for (ind in 1:length(seq_vect)) {
  
  # upper limit of range of cells to be processed
  s1 <- seq_vect[ind]
  
  # lower limit
  if(s1 != last(seq_vect)) {
    s2 <- s1 + 9999
  } else {
    s2 <- nrow(coord_tbl)
  }
  
  tic(str_c("*****  ", ind, " / ",  length(seq_vect), "  ---  ", s1, " - ", s2, " / ", n_cells))
  
  # *****
  tbb[[ind]] <- future_map(s1:s2, function(i) {
    
    coord_tbl %>% 
      filter(id_cell == i) -> coords
    
    st_point(c(coords$x, coords$y)) %>% 
      st_sfc() %>% 
      st_set_crs(st_crs(nbr_stack)) -> pt
    
    nbr_stack %>%
      st_crop(pt) %>% 
      as_tibble() %>% 
      rename(nbr = X) %>% 
      select(date, nbr) %>% 
      
      # Aggregate by month (maximum)
      mutate(date = str_c(year(date), "-", month(date), "-01") %>% as_date()) %>% 
      group_by(date) %>% 
      summarize(nbr = max(nbr, na.rm = T)) %>% 
      ungroup() %>% 
      mutate(nbr = ifelse(is.infinite(nbr), NA, nbr)) %>% 
      
      # Grow to all months
      right_join(tibble(date = seq(as_date("1985-01-01"), as_date("2018-12-01"), by = "1 month"))) %>%
      
      # Impute NAs (linear interpolation)
      mutate(nbr_imp = ts(nbr, start = c(1985,1), frequency = 12) %>% na_interpolation()) -> classif_tbl_1
    
    # BFAST
    classif_tbl_1 %>% 
      pull(nbr_imp) -> nbr_ts
    
    # Initial parameters
    hi <- 24/408 # 24 months min. segment
    bfast_model <- 1
    class(bfast_model) <- "try-error"
    nobp_test <- T
    
    # While loop: 
    # Prevents models without breakpoints (until a certain point: hi = 0.16) or with errors
    # i.e. conditions to break the loop: model has breakpoints or hi is higher than 0.16, plus no errors
    while ((nobp_test == T & hi <= 0.16) | class(bfast_model) == "try-error") {
      
      bfast_model <- try(bfast(nbr_ts, season = "harmonic", max.iter = 5, h = hi))
      hi <- hi + 0.02
      
      # Is the model error-free? Update "nobp_test" with whether the model has breakpoints (F) or not (T)  
      if (class(bfast_model) != "try-error") {
        nobp_test <- bfast_model$nobp$Vt
      }
    }
    
    # Does the model have no breakpoints? Assign an initial and final breakpoint
    if (bfast_model$nobp$Vt == T) {
      
      brkpts <- c(0, length(bfast_model$Yt))
      
    # Does the model have breakpoints? Assign all of them, plus an initial and a final
    } else {
      
      brkpts <- c(0, 
                  bfast_model$output[[length(bfast_model$output)]]$bp.Vt$breakpoints, 
                  length(bfast_model$Yt))
      
    }
    
    # Vector of breakpoints (one id per segment)
    brkpts_id <- 1:(length(brkpts)-1) %>% 
      map(function(x) rep(x, times = brkpts[x+1]-brkpts[x])) %>% 
      unlist()
    
    # Derive temporal metrics per segment
    classif_tbl_1 %>% 
      mutate(brkpts_id = brkpts_id,
             nbr_gapless_unseas = bfast_model$output[[length(bfast_model$output)]]$Vt, # de-seasoned
             nbr_gap_unseas = ifelse(!is.na(nbr), nbr_gapless_unseas, NA), # de-seasoned with gaps (original)
             full_max = max(nbr_gap_unseas, na.rm = T),
             full_min = min(nbr_gap_unseas, na.rm = T),
             full_range = full_max - full_min) %>% 
      group_by(brkpts_id) %>% 
      mutate(seg_mean = mean(nbr_gap_unseas, na.rm = T),
             seg_mean_rel = (seg_mean - full_min)/full_range, # relative to the range of the full ts
             seg_sd = sd(nbr_gap_unseas, na.rm = T),
             seg_sd = ifelse(is.na(seg_sd), 0, seg_sd),
             seg_q10 = quantile(nbr_gap_unseas, 0.1, na.rm = T),
             seg_q90 = quantile(nbr_gap_unseas, 0.9, na.rm = T),
             seg_min = min(nbr_gap_unseas, na.rm = T),
             seg_max = max(nbr_gap_unseas, na.rm = T),
             seg_range_q = seg_q90 - seg_q10,
             seg_range = seg_max - seg_min) %>% 
      ungroup() -> classif_tbl_2
    
    # Use random forest model to predict land cover
    classif_tbl_2 %>% 
      mutate(cover = predict(model_randfor_full, classif_tbl_2, type = "class"),
             x = coords$x,
             y = coords$y) %>% 
      select(x, y, date, cover) -> classif_tbl_3
    
    # Obtain duration and sequence of cycles
    rl <- rle(as.numeric(classif_tbl_3$cover))
    
    # Final tb
    classif_tbl_3 %>% 
      mutate(duration = map(rl$lengths, 
                            function(i) seq(1,i)) %>% unlist(),
             
             cycle = pmap(list(rl$lengths, seq_along(rl$lengths)),
                          function(a, b) rep(b, a)) %>% unlist())
    
  }) %>% bind_rows() # end of future_map
  
  toc()
} # end of for

toc()

saveRDS(tbb, glue::glue(here::here("tbl_cover.RDS")))






# Preparation ----

# Set ejido
ej <- "lg"
ej_n <- 3 # 1 = NM ; 2 = NV ; 3 = LG ; 4 = NB

# Load libraries
library(conflicted)
library(raster)
library(sf)
library(tidyverse)
library(velox)
library(tictoc)
library(imputeTS)
library(bfast)
library(furrr)
library(randomForest)


# Set directory of external files (input/output)
ext_dir <- "/home/cdobler/Documents/calakmul_landsat_analysis/"

# Import ejidos shape
ejido <- read_sf("/media/cdobler/Neobxbaumia/Research/Data_spatial/Shapefiles/03 SYPR shapefiles/Ejidos/Ejidos_select.shp") %>% 
  st_transform("+proj=utm +zone=16 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0") %>% 
  slice(ej_n) %>%  # ejido of interest
  dplyr::select(NOM_NA)

# Get list of VI scenes
scene_files <- str_c(ext_dir, "04_ndmi/") %>% 
  list.files() %>% 
  grep(pattern = "aux", inv = T, value = T)

# Create vi velox
Sys.time()
tic()
vi_stack <- str_c(ext_dir, "04_ndmi/", scene_files) %>% 
  stack() %>%
  crop(extent(ejido) + 60) %>% 
  velox()
toc() # ~ 5.5 min

# Reference raster
r <- raster(str_c(ext_dir, "04_ndmi/", scene_files[1])) %>% 
  raster() %>% 
  crop(extent(ejido) + 60)

# Load model
rand_for_model <- readRDS(glue::glue(here::here("output/model_randfor_full_{ej}.RDS")))

beepr::beep(2); Sys.sleep(0.5); beepr::beep(2)



# Loop ----
plan(multicore, workers = 6)

Sys.time()
tic("MAIN")

seq_vect <- seq(1, ncell(r), 1000)
tbb <- vector("list", length(seq_vect))

for (ind in 1:51) { #:length(seq_vect)) {
  
  s1 <- seq_vect[ind]
  
  if(s1 != last(seq_vect)) {
    s2 <- s1 + 999
  } else {
    s2 <- ncell(r)
  }
  
  tic(str_c("*****  ", ind, " / ",  length(seq_vect), "  ---  ", s1, " - ", s2, " / ", ncell(r)))
  
  # *****
  tbb[[ind]] <- future_map(s1:s2, function(i) {
    
    # Get location
    pt <- st_point(c(xFromCell(r, i), yFromCell(r, i))) %>% st_geometry()
    
    # Extract vi values
    vi_date_tbl <- vi_stack$extract_points(pt) %>%
      t() %>%
      as_tibble() %>%
      rename(vi = V1) %>%
      
      # Obtain dates
      mutate(date = scene_files,
             year = as.integer(str_sub(date, 6, 9)),
             month = str_sub(date, 11, 12),
             qtr = case_when(month %in% c("12", "01", "02") ~ 1L,
                             month %in% c("03", "04", "05") ~ 2L,
                             month %in% c("06", "07", "08") ~ 3L,
                             month %in% c("09", "10", "11") ~ 4L),
             year = ifelse(month == "12", year+1L, year)) %>%
      dplyr::select(year, qtr, vi, -date, -month) %>%
      
      # Aggregate per quarter
      group_by(year, qtr) %>%
      summarize(vi = median(vi, na.rm = T)) %>% # median? max?
      ungroup() %>%
      mutate(date = str_c(year, "_", qtr)) %>%
      dplyr::select(-year, -qtr) %>%
      
      # Grow to include all quarters
      right_join(tibble(date = str_c(rep(1985:2018, each = 4), "_", c(1, 2, 3, 4))), by = "date") %>%
      mutate(cellr = as.integer(i))
    
    # bfast model
    vi_ts <- vi_date_tbl %>%
      .$vi %>%
      ts(start = 1985, frequency = 4) %>%
      na.interpolation()
    
    hi <- 0.06
    bfast_model <- 1
    class(bfast_model) <- "try-error"
    nobp_test <- T
    
    # While to prevent models without breakpoints (until a certain point: hi < 0.2) + without errors
    while ((nobp_test == T & hi < 0.2) | class(bfast_model) == "try-error") {
      
      
      bfast_model <- try(bfast(vi_ts, season = "harmonic", max.iter = 5, h = hi))
      hi <- hi + 0.02
      
      if (class(bfast_model) != "try-error") {
        nobp_test <- bfast_model$nobp$Vt
      }
      
    }
    
    # breakpoints
    if (bfast_model$nobp$Vt == T) {
      
      brkpts <- c(0, length(bfast_model$Yt))
      
    } else {
      
      brkpts <- c(0,
                  bfast_model$output[[length(bfast_model$output)]]$bp.Vt$breakpoints,
                  length(bfast_model$Yt))
      
    }
    
    # Vector of breakpoints (one number per segment)
    brkpts_vector <- 1:(length(brkpts)-1) %>%
      map(function(x) {rep(x, times = brkpts[x+1]-brkpts[x])}) %>%
      unlist()
    
    # Aggregate segments and calculate variables to classify
    expl_tbl <- vi_date_tbl %>%
      mutate(vi_interp_unseas = bfast_model$output[[length(bfast_model$output)]]$Vt,
             vi_gap_unseas = ifelse(!is.na(vi), vi_interp_unseas, NA),
             brkpts_vector = brkpts_vector,
             full_max = max(vi_gap_unseas, na.rm = T),
             full_min = min(vi_gap_unseas, na.rm = T),
             full_range = full_max - full_min) %>%
      group_by(brkpts_vector) %>%
      mutate(seg_mean = mean(vi_gap_unseas, na.rm = T),
             seg_sd = sd(vi_gap_unseas, na.rm = T),
             seg_q10 = quantile(vi_gap_unseas, 0.1, na.rm = T),
             seg_q90 = quantile(vi_gap_unseas, 0.9, na.rm = T),
             seg_min = min(vi_gap_unseas, na.rm = T),
             seg_max = max(vi_gap_unseas, na.rm = T),
             seg_rel = (seg_mean - full_min)/full_range) %>%  # relative to the range of the full ts
      ungroup() %>%
      mutate(seg_sd = ifelse(is.na(seg_sd), 0, seg_sd))
    
    # Predict based on model
    cover_type <- tibble(cover = predict(rand_for_model, expl_tbl, type = "class"),
                         date = str_c(rep(1985:2018, each = 4), "_", c(1, 2, 3, 4)))
    
    rl <- rle(as.numeric(cover_type$cover))
    
    # Final tb
    cover_type %>% 
      mutate(duration = map(rl$lengths, function(i) seq(1,i)) %>% unlist(),
             
             rl = pmap(
               list(rl$lengths,
                    1:length(rl$lengths)),
               function(a, b) rep(b, a)) %>% unlist(),
             
             cellr = i
      )
    
  }) %>% bind_rows() # end of future_map
  
  toc()
} # end of for

toc()
beepr::beep(2); Sys.sleep(0.2); beepr::beep(2); Sys.sleep(0.2); beepr::beep(2)

saveRDS(tbb, glue::glue(here::here("output/tbl_cover_{ej}.RDS")))

# tbb %>% bind_rows() %>% .[!complete.cases(.),] %>% .$cellr %>% unique()

