
# Script to:
#     * Segment time series (BFAST)
#     * Develop random forest model
#     * Assess its accuracy 


# 01. Preparation ---------------------------------------------------------------------------------

# Set ejido
ej <- "nm"
ej_n <- 1 # 1 = NM ; 2 = NV ; 3 = LG ; 4 = NB

# Load libraries
library(conflicted)
library(raster)
library(sf)
library(tidyverse)
library(lubridate)
library(randomForest)
library(tabularaster)
library(velox)
library(tictoc)
library(imputeTS)
library(bfast)
library(furrr)


# Set directory of external files (input/output)
ext_dir <- "/home/cdobler/Documents/calakmul_landsat_analysis/"

# Get list of VI scenes
scene_files <- str_c(ext_dir, "04_ndmi/") %>% 
  list.files() %>% 
  grep(pattern = "aux", inv = T, value = T)

# Empty raster
r <- raster(str_c(ext_dir, "04_ndmi/", scene_files[1])) %>% raster()



# 02. Prepare data --------------------------------------------------------------------------------

# Import training table (visual classification) 
train_tbl <- read_csv(glue::glue(here::here("output/table_classif_cuatri_{ej}.csv")), 
                      na = "NA") %>%
  gather(-date, key = pt, value = cat) %>% 
  mutate(date = str_sub(date, end = 6))

# Import points
train_pts <- read_sf(glue::glue(here::here("data/training_pts_{ej}/pts.gpkg"))) %>% 
  mutate(id = 1:n()) %>% 
  dplyr::filter(id %in% unique(train_tbl$pt)) # remove pts not classified

# train_pts <- read_sf(glue::glue(here::here("data/training_pts_{ej}/training_pts_3.shp"))) %>%
#   mutate(id = 1:n()) %>%
#   dplyr::filter(id %in% unique(train_tbl$pt)) # remove pts not classified

# Import vi stack
Sys.time()
tic()
vi_stack <- str_c(ext_dir, "04_ndmi/", scene_files) %>% 
  stack() %>%
  # crop(extent(train_pts) + 90)
  velox()
toc()  # ~3 minutes without cropping

# Extract vi vals for training points
Sys.time()
tic()
train_vi <- vi_stack$extract_points(train_pts) %>% 
  t() %>% 
  as_tibble() %>% 
  rename_all(list(~unique(train_tbl$pt))) %>% 
  mutate(date = scene_files,
         year = as.integer(str_sub(date, 6, 9)),
         month = str_sub(date, 11, 12),
         qtr = case_when(month %in% c("12", "01", "02") ~ 1L,
                         month %in% c("03", "04", "05") ~ 2L,
                         month %in% c("06", "07", "08") ~ 3L,
                         month %in% c("09", "10", "11") ~ 4L),
         year = ifelse(month == "12", year+1L, year)) %>% 
  dplyr::select(year, qtr, everything(), -date, -month)
toc() # ~5.3 minutes

# Aggregate by cuatri
train_vi_qtr <- train_vi %>% 
  group_by(year, qtr) %>%
  summarize_all(list(~median(., na.rm = T))) %>% # median? max?
  ungroup() %>% 
  mutate_all(list(~ifelse(is.infinite(.), NA, .))) %>% 
  mutate(date = str_c(year, "_", qtr)) %>% 
  dplyr::select(date, everything(), -year, -qtr) %>% 
  
  # Grow the tibble to include all cuatris
  right_join(tibble(date = str_c(rep(1985:2018, each = 4), "_", c(1, 2, 3, 4)))) %>% 
  gather(-1, key = pt, value = vi) %>% 
  mutate(pt = str_pad(pt, 2, "left", "0"),
         year = as.integer(str_split(date, "_", simplify = T)[, 1]),
         qtr = as.integer(str_split(date, "_", simplify = T)[, 2]),
         date_num = year + (qtr/4) - (1/4)) %>% 
  dplyr::select(-date, everything(), date)

# Impute missing observations
vi_ts <- unique(train_vi_qtr$pt) %>% 
  map(function(i){
    
    train_vi_qtr %>% 
      dplyr::filter(pt == i) %>% 
      .$vi %>% 
      ts(start = 1985, frequency = 4) %>% 
      na.interpolation()
    
  }) %>% 
  setNames(unique(train_vi_qtr$pt))

beepr::beep(2); Sys.sleep(0.2); beepr::beep(2)



# 03. BFAST ---------------------------------------------------------------------------------------
plan(multicore, workers = 6)

train_randfor <- future_map(unique(train_vi_qtr$pt), function(i){
  
  # Model
  hi <- 0.06
  bfast_model <- 1
  class(bfast_model) <- "try-error"
  nobp_test <- T

  # While: prevents models without breakpoints (until a certain point: hi < 0.2) or with errors
  while ((nobp_test == T & hi < 0.2) | class(bfast_model) == "try-error") {

    bfast_model <- try(bfast(vi_ts[[i]], season = "harmonic", max.iter = 5, h = hi))
    hi <- hi + 0.02
    
    if (class(bfast_model) != "try-error") {
      nobp_test <- bfast_model$nobp$Vt
    }
  }
  
  # Breakpoints
  if (bfast_model$nobp$Vt == T) {  # no breakpoints
    
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
  expl_tbl <- train_vi_qtr %>% 
    dplyr::filter(pt == i) %>% 
    mutate(vi_interp = bfast_model$Yt,
           vi_interp_unseas = bfast_model$output[[length(bfast_model$output)]]$Vt,
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
    ungroup()
  
  # Response variable
  resp_tbl <- train_tbl %>% 
    mutate(pt = str_pad(pt, 2, "left", "0")) %>% 
    dplyr::filter(pt == i) %>% 
    dplyr::select(-pt)
  
  # Join explanatory and response vars for random forest
  train_randfor <- expl_tbl %>% left_join(resp_tbl)
  
  # Return
  return(train_randfor)
  
}) %>% 
  bind_rows()

# Plot
brkpts_dates <- train_randfor %>% 
  dplyr::filter(pt == "20") %>% 
  group_by(brkpts_vector) %>% 
  summarize(date_num = max(date_num))

train_randfor %>% 
  dplyr::filter(pt == "20") %>% 
  ggplot(aes(x = date_num)) +
  geom_point(aes(y = seg_mean), color = "red", size = 0.3) +
  geom_point(aes(y = vi_gap_unseas)) +
  geom_ribbon(aes(ymin = seg_mean - seg_sd, ymax = seg_mean + seg_sd), fill = "red", alpha = 0.3) +
  #geom_line(aes(y = seg_mean), color = "red", linetype = 3) +
  geom_vline(aes(xintercept = date_num), data = brkpts_dates, linetype = 3) +
  scale_x_continuous(breaks = seq(1985, 2019, 3), minor_breaks = seq(1985, 2019, 1))
  


# 04. Random forest classification ----------------------------------------------------------------

# Remove obs with no classification
train_randfor_sub <- train_randfor %>% 
  mutate(seg_sd = ifelse(is.na(seg_sd), 0, seg_sd)) %>% 
  dplyr::filter(!is.na(cat), cat != "u") %>% 
  mutate(cat = as_factor(cat))

# Split the data
set.seed(123)
split_train <- sample(nrow(train_randfor_sub), 0.75*nrow(train_randfor_sub), replace = F)
train_set <- train_randfor_sub[split_train,]
valid_set <- train_randfor_sub[-split_train,]

# Run random forest with train set
model_randfor <- randomForest(cat ~ seg_mean + seg_sd + seg_rel +
                                seg_min + seg_max + seg_q10 + seg_q90, 
                              data = train_set,
                              importance = T)

# Importance of variables
varImpPlot(model_randfor)

# *****

# Assess the fit (with training data)
train_pred <- predict(model_randfor, train_set, type = "class")
mean(train_pred == train_set$cat) # overall accuracy
    # LG = 92.5 %
    # NV = 95.9 %
    # NM = 95.8 %

# Error matrix
(train_error <- tibble(predicted = train_pred, observed = train_set$cat) %>% 
  group_by(observed, predicted) %>% 
  tally() %>%
  ungroup() %>% 
  mutate(accu = c("hits_forest", "incorrect_for", "incorrect_nonfor", "hits_nonfor"),
         accu_l = c("a", "c", "b", "d")))

# Pierce skill score
train_error %>% 
  arrange(accu_l) %>% 
  .$n %>% 
  {(.[1]*.[4] - .[2]*.[3])/((.[1]+.[3]) * (.[2]+.[4]))}
    # LG = 79.3 %
    # NV = 87.3 %
    # NM = 85.5 %

# *****

# Assess the fit (with independent data)
valid_pred <- predict(model_randfor, valid_set, type = "class")
mean(valid_pred == valid_set$cat) # overall accuracy
    # LG = 91.7 %
    # NV = 92.4 %
    # NM = 93.6 %

# Error matrix
(valid_error <- tibble(predicted = valid_pred, observed = valid_set$cat) %>% 
  group_by(observed, predicted) %>% 
  tally() %>%
  ungroup() %>% 
  mutate(accu = c("hits_forest", "incorrect_for", "incorrect_nonfor", "hits_nonfor"),
         accu_l = c("a", "c", "b", "d")))

# Pierce skill score
valid_error %>% 
  arrange(accu_l) %>% 
  .$n %>% 
  #{(.[1]/(.[1]+.[3])) - (.[2]/(.[2]+.[4]))} # same as below
  {(.[1]*.[4] - .[2]*.[3])/((.[1]+.[3]) * (.[2]+.[4]))}
    # LG = 77.1 %
    # NV = 77.4 %
    # NM = 78.9 %


# *****

# Model with all data
model_randfor_full <- randomForest(cat ~ seg_mean + seg_sd + seg_rel + seg_min + seg_max + seg_q10 + seg_q90, 
                                   data = train_randfor_sub,
                                   importance = T)

saveRDS(model_randfor_full, glue::glue(here::here("output/model_randfor_full_{ej}.RDS")))


# Plot
result_model_tbl <- train_randfor %>%
  mutate(seg_sd = ifelse(is.na(seg_sd), 0, seg_sd)) %>% 
  mutate(cat_predict = predict(model_randfor_full, train_randfor, type = "class")) %>% 
  dplyr::select(pt, date_num, cat, cat_predict, seg_mean, vi_gap_unseas, brkpts_vector) %>% 
  mutate(col_code = factor(case_when((cat == cat_predict & cat == "f") ~ "correct forest",
                                     (cat != cat_predict & cat == "f") ~ "incorrect forest",
                                     (cat == cat_predict & cat == "a") ~ "correct non-forest",
                                     (cat != cat_predict & cat == "a") ~ "incorrect non-forest",
                                     TRUE ~ NA_character_), 
                           levels = c("correct forest", 
                                      "incorrect forest", 
                                      "correct non-forest",
                                      "incorrect non-forest")))

unique(result_model_tbl$pt) %>% 
  walk(function(i) {
    
    brkpts_dates <- train_randfor %>% 
      dplyr::filter(pt == i) %>% 
      group_by(brkpts_vector) %>% 
      summarize(date_num = max(date_num))
    
    result_model_tbl %>% 
      dplyr::filter(pt == i) %>% 
      
      ggplot(aes(x = date_num)) +
      geom_point(aes(y = vi_gap_unseas), color = "dark grey") +
      geom_point(aes(y = seg_mean, color = col_code, shape = col_code), size = 2) +
      scale_color_manual(name = NULL,
                         values = c("#4DAF4A", "#984EA3", "#FF7F00", "#E41A1C"),
                         breaks = c("correct forest", "incorrect forest", "correct non-forest", "incorrect non-forest"),
                         drop = F) +
      scale_shape_manual(values = c(15, 18, 15, 18), guide = "none", drop = F) +
      geom_vline(aes(xintercept = date_num), data = brkpts_dates, linetype = 3) +
      scale_x_continuous(breaks = seq(1985, 2019, 3), minor_breaks = seq(1985, 2019, 1)) +
      labs(subtitle = str_c("Point ", i),
           y = "VI",
           x = NULL) +
      ylim(c(-300, 600))
    
    ggsave(str_c(glue::glue(here::here("output/result_plots_{ej}/"), "pt_", i, ".png")),
           width = 9.5,
           height = 2.2)
    
  })



# **** END ****


