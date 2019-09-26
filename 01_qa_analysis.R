
# *************************************************************************************************
# Script to:
#     * create Landsat masks
#     * analyze distribution of valid pixels per month (random points)
#     * count valid pixels per year (whole scene)
#
# *************************************************************************************************




# 01. Preparation ---------------------------------------------------------------------------------

# Load libraries
library(tidyverse)
library(lubridate)
library(raster)


# External directory (input/output)
ext_dir <- "/home/cdobler/Documents/calakmul_landsat_analysis/"

# List of landsat stacks
list_of_files <- grep(list.files(str_c(ext_dir, "00_cropped_scenes")), 
                      pattern = "aux", 
                      inv = T, 
                      value = T)




# 02. Generate QA layers --------------------------------------------------------------------------

for(i in seq_along(list_of_files)){
  
  file <- list_of_files[i]
  name <- str_split(file, "[.]", simplify = T)[1]
  
  # Print progress
  print(str_c("Processing ", i, "/", length(list_of_files), " --- ", name))
  
  # Process mask + write to disk
  qa_mask <- str_c(ext_dir, "00_cropped_scenes/", file) %>% 
    raster(band = 1)%>%  # qa band
    {
      # 1 for clear pixels, 0 for else
      
      # Landsat 4-7:
      if(str_split(name, "_", simplify = T)[5] != "8"){
        calc(., fun = function(x){ifelse(x == 66, 1, 0)})
      
      # Landsat 8:
      } else {
        calc(., fun = function(x){ifelse(x %in% c(322, 386, 834, 898, 1346), 1, 0)})
      }
    }
  
  writeRaster(qa_mask, 
              str_c(ext_dir, "01_qa_masks/qa_", name, ".tif"), 
              datatype = "INT1U")
  
} # end of for loop

rm(name, i, list_of_files, file)


# 03. Monthly distribution of valid pixels per month ----------------------------------------------

list_of_masks <- grep(list.files(str_c(ext_dir, "01_qa_masks"),
                                 full.names = T), 
                      pattern = "aux", 
                      inv = T, 
                      value = T)

# Bounding polygon
ext_pol <- raster(list_of_masks[1]) %>% sf::st_bbox() %>% sf::st_as_sfc()

# Sample points
set.seed(123)
sample_pts <- sf::st_sample(ext_pol, 1000)

# Extract valid counts
qa_stack <- stack(list_of_masks)
tictoc::tic()
valid_extract <- raster::extract(qa_stack, as(sample_pts, "Spatial")) # ~1.8 min
tictoc::toc()

# Tidy tibble
valid <- valid_extract %>% 
  t() %>% 
  as_tibble(rownames = "date") %>%
  mutate(valid = pmap_dbl(dplyr::select(., V1:V1000), sum, na.rm = T)) %>%
  dplyr::select(date, valid) %>% 
  mutate(valid = valid/1000,  # percent of valid points per scene
         year = str_split(date, "_", simplify = T)[,3],
         month = str_split(date, "_", simplify = T)[,4],
         date = ymd(glue::glue("{year}-{month}-01")))

# Analysis *****

# What is the coverage of scenes across years
valid %>%
  mutate(valid_class = factor(case_when(valid < 0.1 ~ "< 10%",
                                        valid < 0.25 ~ "10-25%",
                                        valid < 0.5 ~ "25-50%",
                                        TRUE ~ "> 50%"), 
                              levels = c("> 50%", "25-50%","10-25%", "< 10%"))) %>% 
  
  ggplot(aes(x = date, y = valid*100)) +
  geom_point(aes(color = valid_class), size = 2, alpha = 0.7, show.legend = F) +
  scale_x_date(date_breaks = "3 years", date_minor_breaks = "1 year", date_labels = "%Y") +
  labs(y = "percent clear", x = NULL,
       subtitle = "Valid pixels per scene (397 scenes)")

ggsave(here::here("output/result_land_change_plots/valid_pixels.png"), width = 7, height = 3.2)

# Which month has more scenes with > X% coverage?
valid %>% 
  mutate('above_10%' = ifelse(valid >= 0.1, 1, 0),
         'above_25%' = ifelse(valid >= 0.25, 1, 0),
         'above_50%' = ifelse(valid >= 0.5, 1, 0)) %>% 
  group_by(month) %>% 
  summarize_at(vars('above_10%':'above_50%'), list(~sum)) %>% 
  gather(-month, key = valid, value = count) %>% 
  
  ggplot(aes(x = month, y = count, group = valid, fill = valid)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(subtitle = "Number of scenes per month") +
  guides(fill=guide_legend(title="Percent clear"))

ggsave(here::here("output/result_land_change_plots/valid_scenes_month.png"), width = 7, height = 3.2)

    # Cuatrimesters: 3-4-5 // 6-7-8 // 9-10-11 // 12-1-2  
    # Trimesters: 2-3-4-5 // 6-7-8-9 // 10-11-12-1

# How do cuatrimesters behave across years?
valid %>%
  mutate(cuatri = case_when(month %in% c("12", "01", "02") ~ "1st",
                            month %in% c("03", "04", "05") ~ "2nd",
                            month %in% c("06", "07", "08") ~ "3rd",
                            month %in% c("09", "10", "11") ~ "4th"),
         
         valid_class = factor(case_when(valid < 0.1 ~ "< 10%",
                                        valid < 0.25 ~ "10-25%",
                                        valid < 0.5 ~ "25-50%",
                                        TRUE ~ "> 50%"), 
                              levels = c("> 50%", "25-50%","10-25%", "< 10%"))
         ) %>% 
  
  
  ggplot(aes(x = as.integer(year), y = valid)) +
  #geom_boxplot() +
  geom_point(aes(color = valid_class), size = 2, alpha = 0.6) +
  facet_grid(cuatri~.) +
  scale_x_continuous(breaks = seq(1985, 2018, 3), minor_breaks = 1985:2018)

# And trimesters?
valid %>%
  mutate(cuatri = case_when(month %in% c("02", "03", "04", "05") ~ "1st",
                            month %in% c("06", "07", "08", "09") ~ "2nd",
                            month %in% c("10", "11", "12", "01") ~ "3rd"),
         
         valid_class = factor(case_when(valid < 0.1 ~ "< 10%",
                                        valid < 0.25 ~ "10-25%",
                                        valid < 0.5 ~ "25-50%",
                                        TRUE ~ "> 50%"), 
                              levels = c("> 50%", "25-50%","10-25%", "< 10%"))
  ) %>% 
  
  
  ggplot(aes(x = as.integer(year), y = valid)) +
  #geom_boxplot() +
  geom_point(aes(color = valid_class), size = 2, alpha = 0.6) +
  facet_grid(cuatri~.) +
  scale_x_continuous(breaks = seq(1985, 2018, 3), minor_breaks = 1985:2018)
  
  

# 04. Calculate annual per-pixel occurences -------------------------------------------------------

### PENDING ###

date_vector <- readRDS(here::here("output/qa_list.RDS"))[[2]]
gc()

tbl_date <- tibble(date = date_vector,
                   year = year(ymd(date)))

which(grepl(1986, tbl_date$year))
which(tbl_date$year == 1986)


#stackApply

for (i in 1985:2017) {
  
  which(tbl_date$year == i) %>% 
    subset(qa_stack, .) %>% 
    sum() %>% 
    writeRaster(str_c("/media/cdobler/PSEUDBOMBAX/tmp_processing/qa_annual_count/qa_annual_count_", i, ".tif"), 
                datatype = "INT1U")
}

seq(1985, 2017) %>% 
  map(function(x) {
    which(tbl_date$year == x) %>% 
      subset(qa_stack, .) %>% 
      sum()
      
  })














  











