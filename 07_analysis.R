#
# SCRIPT TO:
#     * Analyze land-cover durations
#     * Analyze linear trends
#
# *************************************************************************************************


# Preparation ----

# Set ejido
ej <- "lg"
ej_n <- 3 # 1 = NM ; 2 = NV ; 3 = LG ; 4 = NB
ejido_nom <- "La Guadalupe"

# Load libraries
library(conflicted)
library(sf)
library(tidyverse)
library(raster)

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

# Reference raster
r <- raster(str_c(ext_dir, "04_ndmi/", scene_files[1])) %>% 
  raster() %>% 
  crop(extent(ejido) + 60)


# ANALYSIS ----------------------------------------------------------------------------------------

tbb <- readRDS(glue::glue(here::here("output/tbl_cover_{ej}.RDS")))

tbb_cat <- tbb %>% 
  bind_rows() %>% 
  mutate(duration = ifelse(cover == "f" & rl == 1, duration + (10*4), duration),
         duration = ifelse(cover == "a" & rl == 1, duration + (3*4), duration),
         cat = case_when(
           cover == "f" & duration <   5*4 ~ 4L,
           cover == "f" & duration <  10*4 ~ 3L,
           cover == "f" & duration <  20*4 ~ 2L,
           cover == "f" & duration >= 20*4 ~ 1L,
           
           cover == "a" & duration <   3*4 ~ 5L,
           cover == "a" & duration <   6*4 ~ 6L,
           cover == "a" & duration <   9*4 ~ 7L,
           cover == "a" & duration >=  9*4 ~ 8L
         ),
         cat_name = factor(case_when(
           cat == 1 ~ "forest > 20",
           cat == 2 ~ "forest 10-20",
           cat == 3 ~ "forest 5-10",
           cat == 4 ~ "forest 1-5",
           cat == 5 ~ "non-forest 1-3",
           cat == 6 ~ "non-forest 3-6",
           cat == 7 ~ "non-forest 6-9",
           cat == 8 ~ "non-forest > 9"
         ), levels = c("forest > 20",
                       "forest 10-20",
                       "forest 5-10",
                       "forest 1-5",
                       "non-forest 1-3",
                       "non-forest 3-6",
                       "non-forest 6-9",
                       "non-forest > 9")
         )
  )

tbb_date <- tbb_cat %>% 
  mutate(year = str_split(date, "_", simplify = T)[ ,1] %>% as.integer(),
         qtr = str_split(date, "_", simplify = T)[ ,2] %>% as.integer()) %>% 
  dplyr::filter(qtr == 2)

# MAPS
tbb_date %>% 
  dplyr::filter(year %in% c(1997, 2002, 2007, 2012, 2017)) %>% 
  arrange(year, cellr) %>% 
  mutate(x = rep(raster::as.data.frame(r, xy = T)$x, 5),
         y = rep(raster::as.data.frame(r, xy = T)$y, 5)) %>% 
  
  ggplot(aes(x = x, y = y, fill = cat_name)) +
  geom_raster() +
  scale_fill_manual(values = rev(c('#d73027','#f46d43','#fdae61','#fee08b','#d9ef8b','#a6d96a','#66bd63','#1a9850')),
                    name = NULL,
                    drop = F) +
  coord_equal() +
  facet_grid(.~year) +
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        panel.background = element_blank(), legend.position = "none")

ggsave(glue::glue(here::here("output/result_land_change_plots/0_maps_{ej}.png")),
       width = 12.5, height = 3.5) 


# PLOT FOR VS NON-FOR
tbb_date %>% 
  mutate(cat_simple = ifelse(cat > 4, "non-forest", "forest")) %>% 
  group_by(year, cat_simple) %>% 
  tally() %>% 
  dplyr::filter(!is.na(cat_simple)) %>% 
  mutate(tot = sum(n),
         percent = n/tot*100) %>% 
  
  ggplot(aes(x = year, y = percent, fill = cat_simple)) +
  geom_area(color = "black", size = 0.4) +
  scale_fill_manual(values = rev(c('#f46d43','#66bd63')),
                    name = NULL) +
  scale_y_reverse() +
  coord_cartesian(ylim = c(25, 0)) +
  scale_x_continuous(breaks = seq(1985,2018,5), minor_breaks = 1985:2018) +
  labs(subtitle = ejido_nom,
       x = NULL, y = "% of landscape")

ggsave(glue::glue(here::here("output/result_land_change_plots/0_fvsnf_{ej}.png")),
       width = 7, height = 3.2)


# PLOT OF ALL YEARS
tbb_change <- tbb_date %>% 
  group_by(year, cat_name) %>% 
  tally() %>%
  dplyr::filter(year >= 1997) %>% 
  mutate(tot = sum(n),
         prop = n/tot*100) %>% 
  ungroup() %>%
  mutate(h1 = ifelse(!grepl("non", .$cat_name), 1, 2)) %>%
  group_by(year, h1) %>%
  mutate(prop_sum_1 = sum(prop)) %>%
  ungroup() %>%
  mutate(h2 = ifelse(!grepl("> 20", .$cat_name), 2, 1)) %>%
  group_by(year, h2) %>%
  mutate(prop_sum_2 = sum(prop)) %>%
  ungroup()

t_for_sl <- tbb_change %>% 
  dplyr::filter(h1 == 1) %>% 
  group_by(year) %>% 
  summarise(prop_sum = mean(prop_sum_1)) %>% 
  .$prop_sum %>% 
  Kendall::MannKendall(.) %>% 
  {.$sl[[1]]}

t_for <- tbb_change %>% 
  dplyr::filter(h1 == 1) %>% 
  group_by(year) %>% 
  summarise(prop_sum = mean(prop_sum_1)) %>% 
  #{summary(lm(prop_sum ~ year, data = .))} %>% 
  {summary(mblm::mblm(prop_sum ~ year, data = .))} %>% 
  # {tibble(chge_perc = round(coefficients(.)[2,1], 2), 
  #         pval = round(coefficients(.)[2,4], 2))} %>% 
  {tibble(chge_perc = round(coefficients(.)[2,1], 2), 
          pval = round(t_for_sl, 2))} %>% 
  mutate(chge_km = round(chge_perc / (ncell(r)*900/1000000) * 100, 2),
         pval = ifelse(pval < 0.01, "< 0.01", str_c("= ", pval)))

t_oldfor_sl <- tbb_change %>% 
  dplyr::filter(h2 == 1) %>% 
  .$prop_sum_2 %>% 
  Kendall::MannKendall(.) %>% 
  {.$sl[[1]]}

t_oldfor <- tbb_change %>% 
  dplyr::filter(h2 == 1) %>% 
  #{summary(lm(prop_sum_2 ~ year, data = .))} %>%
  {summary(mblm::mblm(prop_sum_2 ~ year, data = .))} %>%
  # {tibble(chge_perc = round(coefficients(.)[2,1], 2), 
  #         pval = round(coefficients(.)[2,4], 2))} %>% 
  {tibble(chge_perc = round(coefficients(.)[2,1], 2), 
          pval = round(t_oldfor_sl, 2))} %>% 
  mutate(chge_km = round(chge_perc / (ncell(r)*900/1000000) * 100, 2),
         pval = ifelse(pval < 0.01, "< 0.01", str_c("= ", pval)))
  

ggplot(tbb_change, aes(x = year)) +
  geom_area(aes(y = prop, fill = cat_name)) +
  geom_line(aes(y = prop_sum_1, group = h1), linetype = 2) +
  geom_line(aes(y = prop_sum_2, group = h2), linetype = 3) +
  #geom_area(color = "black", size = 0.4) +
  scale_fill_manual(values = rev(c('#d73027','#f46d43','#fdae61','#fee08b','#d9ef8b','#a6d96a','#66bd63','#1a9850')),
                    name = NULL) +
  scale_y_reverse() +
  coord_cartesian(ylim = c(35,0)) +
  scale_x_continuous(breaks = seq(1997,2017,5), minor_breaks = 1997:2017) +
  labs(subtitle = ejido_nom,
       x = NULL, y = "% of landscape",
       caption = glue::glue("
                            \nAnnual net forest change: {t_for$chge_km} km2 ({t_for$chge_perc}% (Theil-Sen); Mann-Kendall's p-val {t_for$pval})
                            Annual net old-growth for. change: {t_oldfor$chge_km} km2 ({t_oldfor$chge_perc}% (Theil-Sen); Mann-Kendall's p-val {t_oldfor$pval})"))

ggsave(glue::glue(here::here("output/result_land_change_plots/1_all_{ej}_ts.png")),
       width = 7, height = 3.2)

# PLOT OF FOREST
tbb_date %>% 
  group_by(year, cat) %>% 
  tally() %>% 
  dplyr::filter(year >= 1997, cat <= 4) %>% 
  mutate(tot = sum(n),
         prop = n/tot*100) %>% 
  
  ggplot(aes(x = year, y = prop, fill = factor(cat, levels = c(1:8)))) +
  geom_area(color = "black", size = 0.4) +
  scale_fill_manual(values = rev(c('#d9ef8b','#a6d96a','#66bd63','#1a9850')),
                    labels = c("forest > 20", "forest 10-20", "forest 5-10", "forest 1-5"),
                    name = NULL) +
  coord_cartesian(ylim = c(0,25)) +
  scale_x_continuous(breaks = seq(1997,2017,5), minor_breaks = 1997:2017) +
  labs(subtitle = ejido_nom,
       x = NULL, y = "% of total forest")

ggsave(glue::glue(here::here("output/result_land_change_plots/2_for_{ej}.png")),
       width = 7, height = 3.2)


# PLOT OF NON-FOREST
tbb_date %>% 
  group_by(year, cat) %>% 
  tally() %>% 
  dplyr::filter(year >= 1997, cat >= 5) %>% 
  mutate(tot = sum(n),
         prop = n/tot*100) %>% 
  
  ggplot(aes(x = year, y = prop, fill = factor(cat, levels = c(1:8)))) +
  geom_area(color = "black", size = 0.4) +
  scale_fill_manual(values = rev(c('#d73027','#f46d43','#fdae61','#fee08b')),
                    labels = c("non-forest 1-3", "non-forest 3-6", "non-forest 6-9", "non-forest > 9"),
                    name = NULL) +
  scale_x_continuous(breaks = seq(1997,2017,5), minor_breaks = 1997:2017) +
  labs(subtitle = ejido_nom,
       x = NULL, y = "% of total non-for")

ggsave(glue::glue(here::here("output/result_land_change_plots/3_nonfor_{ej}.png")),
       width = 7, height = 3.2)


# *****
# Create raster stacks
future_map(unique(tbb_date$year), function(y) {
  
  covr <- r
  
  covr[] <- tbb_date %>% 
    dplyr::filter(year == y) %>% 
    arrange(cellr) %>% 
    .$cat
  
  return(covr)
  
}) %>% 
  stack() %>% 
  writeRaster(glue::glue(here::here("output/result_mapstack_{ej}/cover_stack_{ej}.tif")),
              datatype = "INT1U", overwrite = T)


