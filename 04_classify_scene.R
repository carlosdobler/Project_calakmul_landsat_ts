#
# *************************************************************************************************
# Script to:
#     * Map land cover using two-step model
#     * SCRIPT TO BE RUN IN A GCP VIRTUAL MACHINE OF 64 CORES
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
dir <- "/media/cdobler/PSEUDBOMBAX/processed_data/calakmul_landsat_analysis/" # $$$
dir <- "/mnt/disks/input_data/"

# Load models
model_randfor_change_full <-here::here("out_training/model_randfor_change_full.RDS") %>% readRDS() # $$$
model_randfor_full <- readRDS(glue::glue(here::here("out_training/model_randfor_full.RDS"))) # $$$

model_randfor_change_full <- str_c(dir, "model_randfor_change_full.RDS") %>% readRDS()
model_randfor_full <- str_c(dir, "model_randfor_full.RDS") %>% readRDS()

# Obtain files names
str_c(dir, "00_cropped_masked_bands") %>%
  list.files(full.names = T) %>% 
  .[str_detect(., "aux", negate = T)] -> all_files_list

# Extract dates
all_files_list %>% # $$$
  str_split("/", simplify = T) %>%
  .[,8] %>%
  str_split("_", simplify = T) %>%
  .[, 2] %>%
  as_date() -> dates

all_files_list %>%
  str_split("/", simplify = T) %>%
  .[,6] %>%
  str_split("_", simplify = T) %>%
  .[, 2] %>%
  as_date() -> dates

# Number of rows
all_files_list %>% 
  .[1] %>% 
  read_stars(RasterIO = list(bands = 14)) %>% 
  st_bbox() %>% 
  {.[4] - .[2]} %>% 
  {./30} %>% 
  as.vector() -> y_pixels

# How many rows per block?
y_block <- 20

# Sequence to index start of each block
seq(1, y_pixels, y_block) -> y_seq

# Initialize list to save
# vector("list", length(y_seq)) -> tbb
read_rds(here::here("out_training/tbb_tmp.RDS")) -> tbb # $$$
read_rds("/mnt/disks/input_data/tbb_tmp.RDS") -> tbb

# Future options
options(future.globals.maxSize= 10000*1024^2)
plan(multicore, workers = 63)

options(future.fork.enable = TRUE) # $$$
plan(multicore, workers = 7) # $$$


# Loop along blocks

tic("MAIN")

# for (ind in seq_along(y_seq)){
for (ind in 55:length(y_seq)){
  
  tic(str_c(" --- Block ", ind, " of ", length(y_seq),  " processed --- "))
  
  # Import block (stack of nbr)
  # first block? No nYOff argument
  if(ind == 1) {
    
    all_files_list %>%
      read_stars(RasterIO = list(bands = c(14),
                                 nYSize = y_block)) %>%
      merge() %>%
      st_set_dimensions(3, values = dates) %>%
      st_set_dimensions(names = c("x", "y", "date")) -> nbr_stack
    
    # last block? Adjust nYSize
  } else if (ind == length(y_seq)){
    
    all_files_list %>%
      read_stars(RasterIO = list(bands = c(14),
                                 nYOff = y_seq[ind],
                                 nYSize = y_pixels - y_seq[ind] + 1)) %>%
      merge() %>%
      st_set_dimensions(3, values = dates) %>%
      st_set_dimensions(names = c("x", "y", "date")) -> nbr_stack
    
    # all other blocks
  } else {
    
    all_files_list %>%
      read_stars(RasterIO = list(bands = c(14),
                                 nYOff = y_seq[ind],
                                 nYSize = y_block)) %>%
      merge() %>%
      st_set_dimensions(3, values = dates) %>%
      st_set_dimensions(names = c("x", "y", "date")) -> nbr_stack
    
  }
  
  # Table of coordinates (of block)
  nbr_stack %>% 
    slice(date, 1) %>% 
    as_tibble() %>% 
    select(-X) %>% 
    mutate(id = row_number()) -> coord_tbl
  
  
  # LOOP ALONG CELLS 1: differentiate cells that did not change in all the t-s
  
  coord_tbl %>% 
    pull(id) %>%
    future_map(function(i) {
    # map(function(i) {
      
      # print(i)
      
      # Create a point (cell center)
      coord_tbl %>% 
        filter(id == i) -> pt_tbl
      
      st_point(c(pt_tbl$x, pt_tbl$y)) %>%
        st_sfc() %>%
        st_set_crs(st_crs(nbr_stack)) -> pt
      
      # Extract nbr values
      nbr_stack %>%
        st_crop(pt) %>%
        as_tibble() %>%
        rename(nbr = X) %>%
        select(date, nbr) %>% 
      
      # Classify whether its cover changed or not
        summarize(ts_mean = mean(nbr, na.rm = T),
                  ts_sd = sd(nbr, na.rm = T),
                  ts_min = min(nbr, na.rm = T),
                  ts_q01 = quantile(nbr, 0.01, na.rm = T),
                  ts_q05 = quantile(nbr, 0.05, na.rm = T),
                  ts_q10 = quantile(nbr, 0.1, na.rm = T),
                  ts_q90 = quantile(nbr, 0.9, na.rm = T),
                  ts_q95 = quantile(nbr, 0.95, na.rm = T),
                  ts_q99 = quantile(nbr, 0.99, na.rm = T)) %>% 
        
        predict(model_randfor_change_full, ., type = "class") -> change_cover
      
      # Did it remain as Forest all the time-series? Generate all-Forest tbl
      if (change_cover == "N") {
        
        tibble(x = pt_tbl$x,
               y = pt_tbl$y,
               year = 1985:2018,
               cover = factor("F", levels = c("F", "NF")),
               duration = seq(1, length(1985:2018)),
               cycle = 1)
        
      }
    }) %>% bind_rows() -> no_change_tbl
  
  print(" --- Done LOOP 1 --- ")
  
  
  # LOOP ALONG CELLS 2: classify cells that changed
  
  coord_tbl %>% 
    anti_join(no_change_tbl, by = c("x" = "x", "y" = "y")) %>% 
    mutate(id = row_number()) -> change_coord_tbl
  
  change_coord_tbl %>% 
    pull(id) %>% 
    future_map(function(i) {
    # map(function(i) {
      
      # print(str_c(i, " / ", nrow(change_coord_tbl)))
      
      # Create a point (cell center)
      change_coord_tbl %>% 
        filter(id == i) -> pt_tbl
      
      st_point(c(pt_tbl$x, pt_tbl$y)) %>%
        st_sfc() %>%
        st_set_crs(st_crs(nbr_stack)) -> pt
      
      # Extract nbr values
      nbr_stack %>%
        st_crop(pt) %>%
        as_tibble() %>%
        rename(nbr = X) %>%
        select(date, nbr) %>% 
        
        # Aggregate to monthly data
        mutate(date = str_c(year(date), "-", month(date), "-01") %>% as_date()) %>%
        group_by(date) %>%
        {suppressWarnings(summarize(., nbr = max(nbr, na.rm = T)))} %>% 
        ungroup() %>%
        mutate(nbr = ifelse(is.infinite(nbr), NA, nbr)) %>%
        
        # Grow to all months
        right_join(tibble(date = seq(as_date("1985-01-01"), as_date("2018-12-01"), by = "1 month")),
                   by = "date") %>%
        
        # Impute NAs (linear interpolation)
        mutate(nbr_imp = ts(nbr, start = c(1985,1), frequency = 12) %>% na_interpolation()) -> classif_tbl_1
      
      
      # BFAST
      # nbr in "ts" format
      classif_tbl_1 %>%
        pull(nbr_imp) -> nbr_ts
      
      # Initial parameters
      hi <- round(30/length(nbr_ts), digits = 2) # 24 months min. segment
      bfast_model <- 1
      class(bfast_model) <- "try-error"
      nobp_test <- T
      
      # While loop:
      # Prevents models without breakpoints (until a certain point: hi = 0.16) or with errors
      # i.e. conditions to break the loop: model has breakpoints or hi is higher than 0.16, plus no errors
      while ((nobp_test == T & hi <= 0.16) | class(bfast_model) == "try-error") {
        
        bfast_model <- try(bfast(nbr_ts, season = "harmonic", max.iter = 3, h = hi))
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
      
      # Calculate temporal metrics per segment
      classif_tbl_1 %>%
        mutate(brkpts_id = brkpts_id,
               nbr_gapless_unseas = bfast_model$output[[length(bfast_model$output)]]$Vt, # de-seasoned
               nbr_gap_unseas = ifelse(!is.na(nbr), nbr_gapless_unseas, NA), # de-seasoned with gaps (original)
               full_max = max(nbr_gap_unseas, na.rm = T),
               full_min = min(nbr_gap_unseas, na.rm = T),
               full_range = full_max - full_min) %>%
        group_by(brkpts_id) %>%
        {suppressWarnings(mutate(.,
                                 seg_mean = mean(nbr_gap_unseas, na.rm = T),
                                 seg_mean_rel = (seg_mean - full_min)/full_range, # relative to the range of the full ts
                                 seg_sd = sd(nbr_gap_unseas, na.rm = T),
                                 seg_q10 = quantile(nbr_gap_unseas, 0.1, na.rm = T),
                                 seg_q90 = quantile(nbr_gap_unseas, 0.9, na.rm = T),
                                 seg_min = min(nbr_gap_unseas, na.rm = T),
                                 seg_max = max(nbr_gap_unseas, na.rm = T),
                                 seg_range_q = seg_q90 - seg_q10,
                                 seg_range = seg_max - seg_min))} %>%
        ungroup() %>% 
        mutate_at(vars(seg_mean:seg_range), list(~ifelse(is.na(.) | is.infinite(.), 0, .))) -> classif_tbl_2
      
      # Use random forest model to predict land cover
      classif_tbl_2 %>%
        mutate(cover = predict(model_randfor_full, classif_tbl_2, type = "class"),
               x = pt_tbl$x,
               y = pt_tbl$y) %>%
        select(x, y, date, cover) %>% 
        
        mutate(year = year(date),
               cover_int = case_when(cover == "NF" ~ 2L,
                                     cover == "F" ~ 1L)) %>% 
        group_by(x, y, year) %>% 
        summarize(cover_int = max(cover_int)) %>%
        ungroup() %>% 
        mutate(cover = case_when(cover_int == 2 ~ "NF",
                                 cover_int == 1 ~ "F") %>% factor(levels = c("F", "NF"))) %>% 
        select(-cover_int) -> classif_tbl_3
      
      # Obtain duration and sequence of cycles
      rl <- rle(as.numeric(classif_tbl_3$cover))
      
      # Final tb
      classif_tbl_3 %>%
        mutate(duration = map(rl$lengths,
                              function(i) seq(1,i)) %>% unlist(),
               
               cycle = pmap(list(rl$lengths, seq_along(rl$lengths)),
                            function(a, b) rep(b, a)) %>% unlist())
        
    }) %>% bind_rows() -> yes_change_tbl
  
  print(" --- Done LOOP 2 --- ")
  
  bind_rows(no_change_tbl, yes_change_tbl) -> tbb[[ind]]

  toc() # Block
  
  saveRDS(tbb, "/mnt/disks/input_data/tbb_tmp.RDS")

} # end of for

saveRDS(tbb, "/mnt/disks/input_data/tbb_f.RDS")
toc() # MAIN








# MACRO LOOP *****
options(future.globals.maxSize= 6000*1024^2)
options(future.fork.enable = TRUE)
plan(multicore, workers = 6)

seq_vect <- seq(1, nrow(coord_tbl), 5000)
tbb <- vector("list", length(seq_vect))

Sys.time()
tic("MAIN")

# Outer loop: saves in list each 5,000 runs
for (ind in 1:length(seq_vect)) {

  # upper limit of range of cells to be processed
  s1 <- seq_vect[ind]

  # lower limit
  if(s1 != last(seq_vect)) {
    s2 <- s1 + 4999
  } else {
    s2 <- nrow(coord_tbl)
  }

  tic(str_c("*****  ", ind, " / ",  length(seq_vect), "  ---  ", s1, " - ", s2, " / ", nrow(coord_tbl)))

  # *****
  tbb[[ind]] <- future_map(s1:s2, function(i) {

    coord_tbl %>%
      filter(id_cell == i) -> coords

    st_point(c(coords$x, coords$y)) %>%
      st_sfc() %>%
      st_set_crs(st_crs(nbr_stack)) -> pt

    tic()
    nbr_stack %>%
      raster::extract(as(pt, "Spatial"), df = T) -> t
    toc()
      st_crop(pt) %>%
      as_tibble() %>%
      rename(nbr = X) %>%
      select(date, nbr) %>%

      # Aggregate by month (maximum)
      # Block
      # mutate(date = str_c(year(date), "-", month(date), "-01") %>% as_date()) %>%
      # group_by(date) %>%
      # summarize(nbr = max(nbr, na.rm = T)) %>%
      # ungroup() %>%
      # mutate(nbr = ifelse(is.infinite(nbr), NA, nbr)) %>%
      # 
      # # Grow to all months
      # right_join(tibble(date = seq(as_date("1985-01-01"), as_date("2018-12-01"), by = "1 month"))) %>%
      # 
      # # Impute NAs (linear interpolation)
      # mutate(nbr_imp = ts(nbr, start = c(1985,1), frequency = 12) %>% na_interpolation()) -> classif_tbl_1
    
    # Instead:
    mutate(year = year(date),
           month = month(date),
           quarter = case_when(month %in% c(12, 1, 2) ~ 1L,
                               month %in% c(3, 4, 5) ~ 2L,
                               month %in% c(6, 7, 8) ~ 3L,
                               month %in% c(9, 10, 11) ~ 4L),
           year = ifelse(month == 12, year+1, year),
           month = case_when(quarter == 1 ~ 1,
                             quarter == 2 ~ 4,
                             quarter == 3 ~ 7,
                             quarter == 4 ~ 10),
           date = str_c(year, "-", month, "-01") %>% as_date()) %>% 
      
      group_by(date, year, quarter) %>%
      summarize(nbr = mean(nbr, na.rm = T)) %>%
      ungroup() %>%
      mutate(nbr = ifelse(is.nan(nbr), NA, nbr)) %>%
      
      right_join(tibble(date = seq(as_date("1985-01-01"), as_date("2018-12-01"), by = "3 months")), by = "date") %>% 
      
      mutate(nbr_imp = ts(nbr, start = c(1985,1), frequency = 4) %>% na_interpolation()) -> classif_tbl_1
    

    # BFAST
    classif_tbl_1 %>%
      pull(nbr_imp) -> nbr_ts

    # Initial parameters
    #hi <- 24/408 # 24 months min. segment
    hi <- round(10/length(nbr_ts), digits = 2) # 10 quarters
    bfast_model <- 1
    class(bfast_model) <- "try-error"
    nobp_test <- T

    # While loop:
    # Prevents models without breakpoints (until a certain point: hi = 0.16) or with errors
    # i.e. conditions to break the loop: model has breakpoints or hi is higher than 0.16, plus no errors
    while ((nobp_test == T & hi <= 0.16) | class(bfast_model) == "try-error") {

      bfast_model <- try(bfast(nbr_ts, season = "harmonic", max.iter = 3, h = hi))
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
             #seg_sd = ifelse(is.na(seg_sd), 0, seg_sd),
             seg_q10 = quantile(nbr_gap_unseas, 0.1, na.rm = T),
             seg_q90 = quantile(nbr_gap_unseas, 0.9, na.rm = T),
             seg_min = min(nbr_gap_unseas, na.rm = T),
             seg_max = max(nbr_gap_unseas, na.rm = T),
             seg_range_q = seg_q90 - seg_q10,
             seg_range = seg_max - seg_min) %>%
      ungroup() %>% 
      mutate_at(vars(seg_mean:seg_range), list(~ifelse(is.na(.) | is.infinite(.), 0, .))) -> classif_tbl_2

    # Use random forest model to predict land cover
    classif_tbl_2 %>%
      mutate(cover = predict(model_randfor_full, classif_tbl_2, type = "class"),
             x = coords$x,
             y = coords$y) %>%
      select(x, y, date, cover) %>% 
      
      mutate(year = year(date),
             cover_int = case_when(cover == "NF" ~ 2L,
                                   cover == "F" ~ 1L)) %>% 
      group_by(x, y, year) %>% 
      summarize(cover_int = max(cover_int)) %>%
      ungroup() %>% 
      mutate(cover = case_when(cover_int == 2 ~ "NF",
                               cover_int == 1 ~ "F") %>% factor()) %>% 
      select(-cover_int) -> classif_tbl_3

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
  
  saveRDS(tbb, "/mnt/disks/input_data/tbb_nb_l.RDS")
} # end of for

toc()

#saveRDS(tbb, glue::glue(here::here("tbl_cover.RDS")))
saveRDS(tbb, "/mnt/disks/input_data/tbb_nb.RDS")

# element 138 done

