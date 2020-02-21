
# *************************************************************************************************
# Script to:
#     * Join training pts table (quarters) with monthly scenes
#     * Segment time seies of training pts (bfast)
#     * Run Random Forests
#     * Assess the accuracy of the model
#
# *************************************************************************************************

# Load libraries
library(tidyverse)
library(lubridate)
library(sf)
library(stars)
library(bfast)
library(furrr)
library(tictoc)
library(randomForest)
library(patchwork)

# Override future multicore
options(future.fork.enable = TRUE)

# Input/output dir
dir <- "/media/cdobler/PSEUDBOMBAX/processed_data/calakmul_landsat_analysis/"



# PART 1: Model change vs. forest persistence in time series **************************************

# Import training pts (coordinates table)
read_csv(here::here("out_training/training_pts_tbl.csv"),
         col_types = "ddii") %>% 
  select(-cluster) -> training_pts_tbl

# Import change/no-change table 
read_csv(here::here("out_training/training_pts_change.csv"), 
         na = "NA", 
         col_types = "if") %>% 
  left_join(training_pts_tbl, by = "pt") -> change_ts_tbl

# Extract vi values for all points
str_c(dir, "00_cropped_masked_bands") %>% 
  list.files(full.names = T) %>%
  map(function(i){
    
    # Obtain date
    i %>% 
      str_split("/", simplify = T) %>% 
      .[,8] %>% 
      str_split("_", simplify = T) %>% 
      .[, 2] %>% 
      as_date() -> date
    
    print(date)
    
    # Obtain nbr data
    i %>% 
      raster::stack(bands = 14) %>% 
      raster::extract(select(change_ts_tbl, x, y), df = T) %>% 
      rename("nbr" = 2) %>% 
      mutate(pt = change_ts_tbl$pt,
             date = date) %>% 
      select(-ID) %>% 
      as_tibble()
    
  }) %>% 
  bind_rows() -> nbr_training_pts_tbl

# nbr_training_pts_tbl %>% 
#   saveRDS(here::here("out_training/nbr_training_pts_tbl.RDS")) #                                   RDS
# 
readRDS(here::here("out_training/nbr_training_pts_tbl.RDS")) -> nbr_training_pts_tbl

nbr_training_pts_tbl %>% 
  group_by(pt) %>% 
  summarize(ts_mean = mean(nbr, na.rm = T),
            ts_sd = sd(nbr, na.rm = T),
            ts_min = min(nbr, na.rm = T),
            ts_q01 = quantile(nbr, 0.01, na.rm = T),
            ts_q05 = quantile(nbr, 0.05, na.rm = T),
            ts_q10 = quantile(nbr, 0.1, na.rm = T),
            ts_q90 = quantile(nbr, 0.9, na.rm = T),
            ts_q95 = quantile(nbr, 0.95, na.rm = T),
            ts_q99 = quantile(nbr, 0.99, na.rm = T)) %>% 
  mutate(change = change_ts_tbl$change) %>% 
  filter(!is.na(change)) -> ts_stats_training_pts

# Split the data (66% training and 33% validation)
set.seed(234)
split_train <- sample(nrow(ts_stats_training_pts), 0.66*nrow(ts_stats_training_pts), replace = F)
train_set <- ts_stats_training_pts[split_train,]
valid_set <- ts_stats_training_pts[-split_train,]

# Run random forest with train set
model_randfor_change <- randomForest(change ~ ts_mean + ts_sd + ts_min + ts_q01 + ts_q05 + ts_q10,
                                     data = train_set,
                                     importance = T)

# Importance of variables
varImpPlot(model_randfor_change)

# Assess accuracy
valid_pred <- predict(model_randfor_change, valid_set, type = "class")
mean(valid_pred == valid_set$change) # overall accuracy = 93%

# Error matrix
(valid_error <- tibble(predicted = valid_pred, observed = valid_set$change) %>% 
    group_by(observed, predicted) %>% 
    tally() %>%
    ungroup() %>% 
    mutate(accu = c("hits_change", "missed_change", "missed_persistence", "hits_persistence"),
           accu_l = c("a", "c", "b", "d")))

# Peirce skill score
valid_error %>% 
  arrange(accu_l) %>% 
  pull(n) %>% 
  #{(.[1]/(.[1]+.[3])) - (.[2]/(.[2]+.[4]))} # same as below
  {(.[1]*.[4] - .[2]*.[3])/((.[1]+.[3]) * (.[2]+.[4]))} # Peirce Skill Score = 87%

# Run Random Forest with all the data
model_randfor_change_full <- randomForest(change ~ ts_mean + ts_sd + ts_min + ts_q01 + ts_q05 + ts_q10,
                                          data = ts_stats_training_pts,
                                          importance = T)

# Importance of variables
varImpPlot(model_randfor_change_full)

# model_randfor_change_full %>% 
#   saveRDS(here::here("out_training/model_randfor_change_full.RDS")) #                              RDS




# PART 2: Wrangle data and build segment tables ***************************************************

# Import visual classification table
read_csv(here::here("out_training/training_pts_classif.csv"), 
         na = "NA", 
         col_types = str_flatten(rep("c", 451))) %>% 
  rename(quarter = date) -> visual_classif_tbl

# Which points were actually classified?
visual_classif_tbl %>%
  gather(-quarter, key = pt, value = cover) %>% 
  group_by(pt) %>% 
  summarize(sum_na = sum(is.na(cover)),
            total = n()) %>% 
  mutate(prop = sum_na/total) %>% 
  filter(prop < 1) %>% 
  pull(pt) -> pts_with_data

# Remove points with no data from visual classification data
visual_classif_tbl %>% 
  select(quarter, pts_with_data) -> visual_classif_tbl

# Match dates of scenes with year-quarters 
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
  
  mutate(year2 = ifelse(month == 12, year+1, year),
         quarter = str_c(year2, "_", comp)) %>% 
  
  select(year, month, quarter) -> year_quarter_tbl

# Join tables to covert from quarter-years to monthly data 
# Block:
left_join(year_quarter_tbl, visual_classif_tbl, by = "quarter") %>%
  select(-quarter) %>%
  gather(-(year:month), key = pts, val = cover) %>%
  mutate(cover = case_when(cover == "A" ~ 3L,
                           cover == "F" ~ 2L,
                           TRUE ~ 1L)) %>%
  group_by(year, month, pts) %>%
  summarize(cover = max(cover)) %>% # choose non-forest over forest over non-classified
  ungroup() %>%

  mutate(pts = parse_number(pts),
         cover = case_when(cover == 3 ~ "NF",
                           cover == 2 ~ "F",
                           TRUE ~ NA_character_) %>% factor(),
         date = str_c(year,"-", month, "-01") %>% as_date()) %>%
  select(-year, -month) -> training_tbl_1

# Instead:
# visual_classif_tbl %>%
#   gather(-quarter, key = pts, val = cover) %>%
#   mutate(cover = case_when(cover == "A" ~ "NF",
#                            cover == "F" ~ "F",
#                            TRUE ~ NA_character_) %>% factor(),
#          pts = parse_number(pts),
#          year = str_sub(quarter, end = 4),
#          quarter = str_sub(quarter, 6),
#          month = case_when(quarter == 1 ~ 1,
#                            quarter == 2 ~ 4,
#                            quarter == 3 ~ 7,
#                            quarter == 4 ~ 10),
#          date = str_c(year, "-", month, "-01") %>% as_date) %>%
#   select(date, year, quarter, cover, pts) -> training_tbl_1

# VI table (monthly)
# Block:
nbr_training_pts_tbl %>% 
  
  mutate(date = str_c(year(date), "-", month(date), "-01") %>% as_date()) %>%
  group_by(pt, date) %>%
  summarize(nbr = max(nbr, na.rm = T)) %>%
  ungroup() %>%
  mutate(nbr = ifelse(is.infinite(nbr), NA, nbr)) -> vi_tbl_1

# Instead
# nbr_training_pts_tbl %>% 
#   
#   mutate(year = year(date),
#          month = month(date),
#          quarter = case_when(month %in% c(12, 1, 2) ~ 1L,
#                              month %in% c(3, 4, 5) ~ 2L,
#                              month %in% c(6, 7, 8) ~ 3L,
#                              month %in% c(9, 10, 11) ~ 4L),
#          year = ifelse(month == 12, year+1, year),
#          month = case_when(quarter == 1 ~ 1,
#                            quarter == 2 ~ 4,
#                            quarter == 3 ~ 7,
#                            quarter == 4 ~ 10),
#          date = str_c(year, "-", month, "-01") %>% as_date()) %>%
#   group_by(pt, date, year, quarter) %>%
#   summarize(nbr = mean(nbr, na.rm = T)) %>% # max? min? mean?
#   ungroup() %>%
#   mutate(nbr = ifelse(is.infinite(nbr) | is.nan(nbr), NA, nbr)) %>%
#   select(-year, -quarter) -> vi_tbl_1

# Join tables (classified with vi)
training_pts_tbl %>% 
  filter(pt %in% parse_number(pts_with_data)) %>% 
  mutate(id = row_number()) %>% 
  
  pmap_df(function(x, y, pt, id){
    
    print(str_c("Processing point ", id))
    
    # Join classified data with vi values and expand table to include to all months
    # Block
    training_tbl_1 %>%
      filter(pts == pt) %>%
      select(-pts) -> cl
    
    vi_tbl_1 %>% 
      rename(pts = pt) %>% 
      filter(pts == pt) %>% 
      select(-pts) %>% 
      left_join(cl, by = "date") %>% 
    
      right_join(tibble(date = seq(as_date("1985-01-01"), as_date("2018-12-01"), by = "1 month")), 
                 by = "date") %>%

      # Impute NAs (linear interpolation)
      mutate(nbr_imp = ts(nbr, start = c(1985,1), frequency = 12) %>% imputeTS::na_interpolation(),
             pt = pt) %>%

      # Arrange columns
      select(pt, date, cover, nbr, nbr_imp)
    
    # Instead:
    # training_tbl_1 %>%
    #   filter(pts == pt) %>%
    #   select(-pts) -> cl
    # 
    # vi_tbl_1 %>% 
    #   rename(pts = pt) %>% 
    #   filter(pts == pt) %>% 
    #   select(-pts) %>% 
    #   left_join(cl, by = "date") %>% 
    #   mutate(nbr_imp = ts(nbr, start = c(1985, 1), frequency = 4) %>% imputeTS::na_interpolation(),
    #          pt = pt) %>%
    #   select(pt, date, year, quarter, cover, nbr, nbr_imp)
    
  }) -> training_tbl_2



# PART 2: BFAST (time series segmentation) ********************************************************

# Parallelize
plan(multicore, workers = 6)

Sys.time()
tic()

training_tbl_2 %>% 
  pull(pt) %>% 
  unique() %>% 

  future_map(function(p){
    
    # NBR as ts object
    training_tbl_2 %>% 
      filter(pt == p) %>% 
      pull(nbr_imp) %>% 
      ts(start = c(1985,1), frequency = 12) -> nbr_ts
      #ts(start = c(1985,1), frequency = 4) -> nbr_ts
    
    # Segment ts with BFAST
    
    # Initial parameters
    hi <- round(24/length(nbr_ts), digits = 2) # 24 months min. segment
    #hi <- round(10/length(nbr_ts), digits = 2) # 10 quarters
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
    training_tbl_2 %>% 
      filter(pt == p) %>% 
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
      ungroup()
    
  }) %>% 
  bind_rows() -> training_tbl_3
  
toc() # ~ 9 min.



# PART 3: Random forests **************************************************************************

# Remove months without cover data
training_tbl_3 %>% 
  filter(!is.na(cover)) -> training_tbl_4

# Split the data (75% training and 25% validation)
set.seed(234)
split_train <- sample(nrow(training_tbl_4), 0.75*nrow(training_tbl_4), replace = F)
train_set <- training_tbl_4[split_train,]
valid_set <- training_tbl_4[-split_train,]

# Run random forest with train set
model_randfor <- randomForest(cover ~ seg_mean + seg_mean_rel + seg_q10 + seg_q90 + seg_range_q + seg_range,
                              data = train_set,
                              importance = T)

# Importance of variables
varImpPlot(model_randfor)



# PART 4: Accuracy assessment *********************************************************************

# *** Fit vs training data ***
train_pred <- predict(model_randfor, train_set, type = "class")
mean(train_pred == train_set$cover) # overall accuracy = 96% (quarters: 91 max; 91 min; 93 mean)

# Error matrix
(train_error <- tibble(predicted = train_pred, observed = train_set$cover) %>% 
    group_by(observed, predicted) %>% 
    tally() %>%
    ungroup() %>% 
    mutate(accu = c("hits_forest", "missed_for", "missed_nonfor", "hits_nonfor"),
           accu_l = c("a", "c", "b", "d")))

# Peirce Skill Score
train_error %>% 
  arrange(accu_l) %>% 
  pull(n) %>% 
  {(.[1]*.[4] - .[2]*.[3])/((.[1]+.[3]) * (.[2]+.[4]))} # Peirce Skill Score = 93% (quarters: 82 max; 82 min; 86 mean)


# *** Fit vs validation (independent) data ***
valid_pred <- predict(model_randfor, valid_set, type = "class")
mean(valid_pred == valid_set$cover) # overall accuracy = 96.4% (quarters: 89 max; 87 min; 89 mean)

# Error matrix
(valid_error <- tibble(predicted = valid_pred, observed = valid_set$cover) %>% 
    group_by(observed, predicted) %>% 
    tally() %>%
    ungroup() %>% 
    mutate(accu = c("hits_forest", "missed_for", "missed_nonfor", "hits_nonfor"),
           accu_l = c("a", "c", "b", "d")))

# Peirce skill score
valid_error %>% 
  arrange(accu_l) %>% 
  pull(n) %>% 
  #{(.[1]/(.[1]+.[3])) - (.[2]/(.[2]+.[4]))} # same as below
  {(.[1]*.[4] - .[2]*.[3])/((.[1]+.[3]) * (.[2]+.[4]))} # Peirce Skill Score = 92.6% (quarters: 77% max; 72 min; 78 mean)



# PART 5: Re-run model with all the data **********************************************************

# Model with all data
model_randfor_full <- randomForest(cover ~ seg_mean + seg_mean_rel + seg_q10 + seg_q90 + seg_range_q + seg_range,
                                   data = training_tbl_4,
                                   importance = T)

# Importance of variables
varImpPlot(model_randfor_full)

saveRDS(model_randfor_full, glue::glue(here::here("out_training/model_randfor_full.RDS"))) #         RDS 
saveRDS(model_randfor_full, glue::glue(here::here("out_training/model_randfor_full_quarters.RDS")))



# PART 6: Visual inspection of results ************************************************************

training_tbl_classif <- training_tbl_3 %>% 
  mutate(predicted = predict(model_randfor_full, training_tbl_3, type = "class")) %>% 
  mutate(accuracy = case_when(
    (cover == "F" & predicted == "F") | (cover == "NF" & predicted == "NF") ~ "hit",
    (cover == "F" & predicted == "NF") | (cover == "NF" & predicted == "F") ~ "miss",
    TRUE ~ "not classified"),
    accuracy = factor(accuracy, levels = c("hit", "miss", "not classified")))

plan(multicore, workers = 6)

# Loop through year-quarters
str_c(rep(1985:2018, each = 4), "_", 1:4) %>%
  map(function(yq){
    
    str_c(dir, "01_three_month_overlays") %>% 
      list.files() %>% 
      str_sub(8, 13) -> available_yr_quarters
    
    # Is there a composite for this quarter? Generate data for chip.
    if(yq %in% available_yr_quarters){
      
      str_c(dir, "01_three_month_overlays/comp3m_", yq, ".tif") %>% 
        read_stars(RasterIO = list(bands = c(4, 5, 3))) -> r
      
      # Loop through points
      training_pts_tbl %>% 
        filter(pt %in% parse_number(pts_with_data)) %>% 
        
        future_pmap(function(x, y, pt){
          
          #print(pt)
          
          st_point(c(x, y)) %>% 
            st_sfc() %>% 
            st_set_crs(st_crs(r)) %>% 
            st_buffer(dist = 450, endCapStyle = "SQUARE") -> pol
          
          r %>% 
            st_crop(pol) %>% 
            as("Raster") %>% 
            RStoolbox::ggRGB(r = 1, g = 2, b = 3,
                             stretch = "lin",
                             ggObj = F) %>% 
            
            mutate(pt = pt,
                   yr_quarter = yq)
          
        }) %>% 
        bind_rows()
      
      # Is there no composite? Create empty chip.
    } else {
      
      training_pts_tbl %>% 
        filter(pt %in% parse_number(pts_with_data)) %>% 
        mutate(fill = "#FFFFFFFF",
               yr_quarter = yq) %>% 
        select(x, y, fill, pt, yr_quarter)
      
    }
    
  }) %>% 
  bind_rows() -> plot_tbl
  

# Loop through points
training_pts_tbl %>% 
  filter(pt %in% parse_number(pts_with_data)) %>% 
  
  pmap(function(x, y, pt){
    
    print(pt)
    
    # Create 9x9 square polygon
    st_point(c(x, y)) %>% 
      st_buffer(dist = 45, endCapStyle = "SQUARE") %>% 
      as("Spatial") %>% 
      broom::tidy() -> pix
    
    # Generate layout of chips
    plot_tbl %>% 
      rename(pts = pt) %>% 
      filter(pts == pt) %>%
      
      {
        ggplot() +
          geom_raster(data = ., aes(x = x, y = y, fill = fill)) +
          scale_fill_identity() +
          geom_polygon(data = pix, aes(x = long, y = lat), fill = NA, color = "black", size = 0.3) +
          coord_equal()+
          theme_void() +
          facet_wrap(~yr_quarter, ncol = 17)
      } -> chips
    
    # Change classification
    change_ts_tbl %>% 
      rename(pts = pt) %>% 
      filter(pts == pt) %>% 
      pull(change) %>% 
      as.character() -> change
    
    # Time-series plot
    training_tbl_classif %>% 
      rename(pts = pt) %>% 
      filter(pts == pt) %>% 
      group_by(brkpts_id) %>% 
      summarise(date = max(date)) %>% 
      bind_rows(tibble(brkpts_id = 0, date = as_date("1985-01-01")))-> brkpts_dates
    
    training_tbl_classif %>% 
      rename(pts = pt) %>% 
      filter(pts == pt) %>% 
      
      ggplot(aes(x = date)) +
      geom_point(aes(y = seg_mean, color = predicted), size = 1, shape = 15) +
      geom_vline(aes(xintercept = date), data = brkpts_dates, linetype = 3, color = "gray30") +
      geom_point(aes(y = nbr, shape = accuracy), size = 2) +
      scale_color_manual(values = c("#4DAF4A", "#FF7F00"), drop = F) +
      scale_shape_manual(values = c(19, 1, 4), drop = F) +
      scale_x_date(breaks = "3 years", minor_breaks = "1 year", date_labels = "%Y") +
      theme(panel.grid.minor.y = element_blank(),
            panel.grid.major.y = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank(),
            axis.title.x = element_blank(),
            axis.text.x = element_text(hjust = 0)) +
      labs(title = str_c("Point ", pt),
           subtitle = str_c("Change: ", change)) -> plot_ts
    
    # Patch time-series and chips
    p <- plot_ts + chips + plot_layout(ncol = 1, heights = c(1.5, 7))
    
    # Save
    suppressMessages(ggsave(here::here("chips_classif/", str_c("chips453_pt_", str_pad(pt, 3, pad = "0"), ".png")),
                            plot = p,
                            width = 14,
                            height = 9))
    
  })
