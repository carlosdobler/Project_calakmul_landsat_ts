#
# *************************************************************************************************
# Script to:
#     * Analyze forest cover trends
#
# *************************************************************************************************

library(tidyverse)
library(lubridate)
library(furrr)


# Load forest cover table
here::here("out_training/tbb_nb.RDS") %>% 
  readRDS() %>% 
  bind_rows() -> tbb
  
tbb %>% 
  filter(is.na(cover)) # check (none)

# "/mnt/disks/input_data/tbb_nb.RDS" %>% readRDS() %>% bind_rows() %>% filter(!is.na(cover)) -> tbb


# Average duration of segments to assign to first cycles
tbb %>% 
  filter(cycle != 1) %>% # all but first 
  group_by(x, y, cover) %>% 
  summarize(duration = max(duration)) %>%
  ungroup() %>% 
  group_by(cover) %>% 
  summarize(duration_mean = mean(duration) %>% round(),
            duration_sd = sd(duration) %>% round()) -> stats_duration_1st_cycle

# Gnerate random durations for first cycles
tbb %>% 
  filter(cycle == 1, cover == "F") %>% 
  group_by(x, y) %>% 
  mutate(rand_dur_inic = sample(rnorm(50, mean = 25, sd = 8), 1) %>% round()) %>% 
  ungroup() %>% 
  mutate(duration = duration + rand_dur_inic) %>% 
  select(-rand_dur_inic) -> rand_dur_inic_F

tbb %>% 
  filter(cycle == 1, cover == "NF") %>% 
  group_by(x, y) %>% 
  mutate(rand_dur_inic = sample(rnorm(50, mean = 6, sd = 5), 1) %>% round()) %>% 
  ungroup() %>% 
  mutate(rand_dur_inic = ifelse(rand_dur_inic < 0, 0, rand_dur_inic), 
         duration = duration + rand_dur_inic) %>% 
  select(-rand_dur_inic) -> rand_dur_inic_NF

# Join tables
tbb %>% 
  filter(cycle != 1) %>% 
  bind_rows(rand_dur_inic_F) %>% 
  bind_rows(rand_dur_inic_NF) -> tbb

# * * * *

tbb %>% 
  mutate(cat = case_when(
    cover == "F" & duration < 5 ~ "Forest 1-5",
    cover == "F" & duration < 10 ~ "Forest 5-10",
    cover == "F" & duration < 20 ~ "Forest 10-20",
    cover == "F" & duration >= 20 ~ "Forest > 20",
    
    cover == "NF" & duration < 3 ~ "Non-forest 1-3",
    cover == "NF" & duration < 6 ~ "Non-forest 3-6",
    cover == "NF" & duration < 9 ~ "Non-forest 6-9",
    cover == "NF" & duration >= 3 ~ "Non-forest > 9") %>% 
      
      factor(levels = c("Forest > 20",
                        "Forest 10-20",
                        "Forest 5-10",
                        "Forest 1-5",
                        "Non-forest 1-3",
                        "Non-forest 3-6",
                        "Non-forest 6-9",
                        "Non-forest > 9"))
  ) -> tbb

tbb %>%
  group_by(x, y) %>% 
  summarize() %>% 
  nrow() -> n_cells

# * * * * * maps

tbb %>% 
  filter(year %in% c(1988, 1994, 2000, 2006, 2012, 2018)) %>% 
  
  ggplot(aes(x = x, y = y, fill = cat)) +
  geom_raster() +
  scale_fill_manual(values = rev(c('#d73027','#f46d43','#fdae61','#fee08b','#d9ef8b','#a6d96a','#66bd63','#1a9850')),
                    name = NULL,
                    drop = F) +
  coord_equal() +
  facet_grid(.~year)

# *****

tbb %>% 
  filter(year >= 1997) %>% 
  group_by(year, cat) %>% 
  tally() %>% 
  
  mutate(tot = sum(n),
         prop = n/tot*100) %>% 
  ungroup() %>%
  
  mutate(h1 = ifelse(str_detect(cat, "Non"), 2L, 1L)) %>%
  group_by(year, h1) %>%
  mutate(prop_sum_1 = sum(prop)) %>%
  ungroup() %>%
  
  mutate(h2 = ifelse(str_detect(cat, "> 20"), 1L, 2L)) %>%
  group_by(year, h2) %>%
  mutate(prop_sum_2 = sum(prop)) %>%
  ungroup() -> tbb_change


# Stats without controlling for old-growth:

tbb_change %>% 
  filter(h1 == 1) %>% 
  group_by(year, prop_sum_1) %>% 
  summarise() %>% 
  pull(prop_sum_1) %>% 
  Kendall::MannKendall() %>% 
  {.$sl[[1]]} %>% 
  round(2) -> mk_plevel

tbb_change %>% 
  filter(h1 == 1) %>% 
  group_by(year, prop_sum_1) %>% 
  summarise() %>% 
  {summary(mblm::mblm(prop_sum_1 ~ year, data = .))} %>% 
  {tibble(change_perc = round(coefficients(.)[2,1], 2),
          mk_plevel = mk_plevel)} %>% 
  mutate(change_km = round(change_perc/n_cells*900/1000000, digits = 3),
         change_km = ifelse(change_km < 0.01, "< 0.01", change_km),
         mk_plevel = ifelse(mk_plevel < 0.01, "< 0.01", str_c("= ", mk_plevel))) -> trend_stats


# Controlling for old-growth:

tbb_change %>% 
  filter(h2 == 1) %>% 
  pull(prop_sum_2) %>%
  Kendall::MannKendall() %>% 
  {.$sl[[1]]} %>% 
  round(2) -> mk_plevel_og

tbb_change %>% 
  filter(h2 == 1) %>% 
  {summary(mblm::mblm(prop_sum_2 ~ year, data = .))} %>% 
  {tibble(change_perc = round(coefficients(.)[2,1], 2),
          mk_plevel = mk_plevel_og)} %>% 
  mutate(change_km = round(change_perc/n_cells*900/1000000, digits = 3),
         change_km = ifelse(change_km < 0.01, "< 0.01", change_km),
         mk_plevel = ifelse(mk_plevel_og < 0.01, "< 0.01", str_c("= ", mk_plevel_og))) -> trend_stats_og


# Plot ****

tbb_change %>% 
  filter(h1 == 2) %>% 
  group_by(year, prop_sum_1) %>% 
  summarise() -> change_for

tbb_change %>% 
  filter(h2 == 2) %>% 
  group_by(year, prop_sum_2) %>% 
  summarise() -> change_for_og

ggplot(tbb_change, aes(x = year)) +
  geom_area(aes(y = prop, fill = cat)) +
  geom_line(data = change_for, aes(y = prop_sum_1), linetype = 2) +
  geom_line(data = change_for_og, aes(y = prop_sum_2), linetype = 3) +
  scale_fill_manual(values = rev(c('#d73027','#f46d43','#fdae61','#fee08b','#d9ef8b','#a6d96a','#66bd63','#1a9850')),
                    name = NULL) +
  scale_y_reverse() +
  coord_cartesian(ylim = c(20,0)) +
  scale_x_continuous(breaks = seq(1997,2017,5), minor_breaks = seq(1998, 2018)) +
  labs(x = NULL, y = "% of landscape",
       subtitle = "Nuevo Becal",
       caption = glue::glue("
                            \nAnnual net forest change: {trend_stats$change_km} km2 ({trend_stats$change_perc}% (Theil-Sen); Mann-Kendall's p-val {trend_stats$mk_plevel})
                            Annual net old-growth for. change: {trend_stats_og$change_km} km2 ({trend_stats_og$change_perc}% (Theil-Sen); Mann-Kendall's p-val {trend_stats_og$mk_plevel})"))

ggsave(glue::glue(here::here("nb.png")),
       width = 7, height = 3.2)
  



  