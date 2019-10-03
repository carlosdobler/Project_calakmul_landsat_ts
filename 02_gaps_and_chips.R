
# *************************************************************************************************
# Script to:
#     * Calculate and analyze gaps
#     * Choose three gap-less scenes and run unsupervised classification
#     * Randomly distribute training points among resulting classes
#     * Produce CHIPS for each point
# 
# *************************************************************************************************


# Load libraries
library(tidyverse)
library(lubridate)
library(stars)
library(furrr)
library(tictoc)
library(patchwork)

# Input/output dir
dir <- "/media/cdobler/PSEUDBOMBAX/processed_data/calakmul_landsat_analysis/"


# PART 1: Calculate gaps per scene ****************************************************************

tic()

# Loop through all scenes
seq_along(list.files(str_c(dir, "00_cropped_masked_bands/"))) %>% 
  map_df(function(i){
    
    # Read QA band and convert to vector
    str_c(dir, "00_cropped_masked_bands/") %>% 
      list.files(full.names = T) %>% 
      .[i] %>% 
      read_stars(RasterIO = list(bands = 1)) %>% 
      pull(1) -> qa_data
    
    # Calculate total number of cells
    qa_data %>% 
      length() -> total_cell
    
    # Calculate percent of scene with gaps
    qa_data %>% 
      is.na(.) %>% 
      sum() %>% 
      {./total_cell * 100} -> gaps 
    
    # Get file name (of scene)
    str_c(dir, "00_cropped_masked_bands/") %>% 
      list.files() %>% 
      .[i] -> file
    
    # Get date
    file %>% 
      str_split("_", simplify = T) %>% .[,2] %>% 
      as_date() -> date
    
    # Create tibble row
    bind_cols(gaps = gaps, date = date, file = file)
    
  }) -> gaps_tbl # this tbl can be used to analyze cloud cover and gaps

toc() # ~ 4 min

# analysis pending...



# PART 2: Create sample points ********************************************************************

# Select gap-less scenes
gaps_tbl %>%
  filter(gaps < 5) %>% 
  # slice(1, 6, 11, 13, 16) -> selected_scenes
  slice(1, 11, 16) %>% 
  pull(file) -> selected_scenes

# Create table with NDMI of all pixels of selected scenes
str_c(dir, "00_cropped_masked_bands/", selected_scenes) %>% 
  read_stars(RasterIO = list(bands = 12)) %>% 
  merge() %>% 
  as_tibble() %>% 
  spread(key = X1, value = X) %>% 
  filter(y < 2091625 - (16*30), # avoid edges
         y > 1983335 + (16*30),
         x < 274605 - (16*30),
         x > 232305 + (16*30)) -> ndmi_tbl

# Select 1/5 of rows
set.seed(246)
ndmi_tbl %>% 
  drop_na() %>% 
  sample_frac(0.2) -> ndmi_tbl

# Classify (unsupervised)
set.seed(246)
kmeans(ndmi_tbl[, -(1:2)],
       iter.max = 500,
       centers = 15,
       nstart = 5,
       algorithm = "Lloyd") -> clusters
  
ndmi_tbl %>% 
  mutate(cluster = clusters$cluster) -> ndmi_tbl

# Select 30 random points per cluster
set.seed(246)
ndmi_tbl %>% 
  group_by(cluster) %>% 
  sample_n(30) %>% 
  select(x, y, cluster) %>% 
  ungroup() %>% 
  mutate(pt = row_number()) -> training_pts_tbl

write_csv(training_pts_tbl,
          here::here("out_training/training_pts_tbl.csv"))



# PART 3: Create chips **************************************************************************** 

# Loop through points
training_pts_tbl %>% 
  pmap(function(x, y, pt, ...){
    
    print(glue::glue("Processing point {pt} of 450 ..."))
    
    # All quatrimester-years
    quatri_vect <- str_c("comp3m_", rep(1985:2018, each = 4), "_", 1:4, ".tif")
    
    # Loop through all quatrimester-years (composites)
    seq_along(quatri_vect) %>% 
      map_df(function(d){
        
        # Is there a composite for this run? Create a chip
        if(quatri_vect[d] %in% list.files(str_c(dir, "01_three_month_overlays/"))){
          
          str_c(dir, "01_three_month_overlays/", quatri_vect[d]) %>% 
            read_stars(RasterIO = list(bands = c(4, 5, 3), 
                                       nXOff = (x - 232305)/30 - (29/2-1), # n of cells between origin and point, minus 14 cells
                                       nYOff = (y - 2091615)/-30 - (29/2-1), 
                                       nXSize = 29, # 14 + 14 cells surrounding
                                       nYSize = 29)) -> chip
          
          
          # Create tibble of chip
          chip %>% 
            as("Raster") %>% 
            RStoolbox::ggRGB(r = 1, g = 2, b = 3,
                  stretch = "lin",
                  ggObj = F) %>% 
            
            mutate(pt = pt,
                   quatri = quatri_vect[d] %>% str_sub(8, 13))
          
        # Is there no composite for this run? Create an empty chip  
        } else {
          
          tibble(x = x,
                 y = y,
                 fill = "#FFFFFFFF",
                 pt = pt,
                 quatri = quatri_vect[d] %>% str_sub(8, 13))
          
        }
        
      }) -> plot_tbl
    
    # Create 9x9 square polygon
    st_point(c(x, y)) %>% 
      st_buffer(dist = 45, endCapStyle = "SQUARE") %>% 
      as("Spatial") %>% 
      broom::tidy() -> pix
    
    # Generate layout of chips
    ggplot() +
      geom_raster(data = plot_tbl, aes(x = x, y = y, fill = fill)) +
      scale_fill_identity() +
      geom_polygon(data = pix, aes(x = long, y = lat), fill = NA, color = "black", size = 0.3) +
      coord_equal()+
      theme_void() +
      facet_wrap(~quatri, ncol = 17) -> chips
    
    # Time-series plots
    # Loop through all raw scenes
    str_c(dir, "00_cropped_masked_bands") %>% 
      list.files(full.names = T) %>%
      map_df(function(i){
        
        # Obtain data on training point (tass cap and ndmi for a single pixel)
        read_stars(i, 
                   RasterIO = list(bands = c(8:10, 12), 
                                   nXOff = (x - 232305)/30 +0.5, 
                                   nYOff = (y - 2091615)/-30 +0.5, 
                                   nXSize = 1, 
                                   nYSize = 1)) %>%
          pull(1) -> vi
        
        # Obtain date
        i %>% 
          str_split("/", simplify = T) %>% 
          .[,8] %>% 
          str_split("_", simplify = T) %>% 
          .[, 2] %>% 
          as_date() -> date
        
        # Create tibble
        tibble(date = date, 
               bri = vi[1],
               gre = vi[2],
               wet = vi[3],
               ndmi = vi[4])
        
      }) -> vi_tbl
    
    vi_tbl %>% 
      mutate_at(vars(2:5), list(~(.-mean(., na.rm = T))/sd(., na.rm = T))) %>% # standardize
      gather(-date, key = vi, value = val) %>% 
      mutate(vi = factor(vi, levels = c("bri", "gre", "wet", "ndmi"))) %>% 
      
      # Generate time series plot
      ggplot(aes(x = date, y = val)) +
      geom_smooth(span = 0.1, linetype = 0) +
      geom_point(aes(color = val)) +
      scale_color_viridis_c() +
      facet_grid(vi~., scales = "free") +
      scale_x_date(breaks = "3 years", minor_breaks = "1 year", date_labels = "%Y") +
      theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
            axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(),
            axis.title.x = element_blank(),
            legend.position = "none") +
      labs(title = str_c("Point ", pt)) -> plot
    
    # Patch time-series and chips
    p <- plot + chips + plot_layout(ncol = 1, heights = c(1, 3))
    
    # Save
    suppressMessages(ggsave(here::here("out_training/", str_c("chips453_pt_", str_pad(pt, 3, pad = "0"), ".png")),
                            plot = p,
                            width = 14,
                            height = 10))
    
  })

# Generate table for visual classification
quatri_vect <- str_c(rep(1985:2018, each = 4), "_", 1:4)

matrix(NA, length(quatri_vect), nrow(training_pts_tbl)) %>%
  as_tibble() %>% 
  rename_all(list(~str_c("pt_", str_pad(seq_len(nrow(training_pts_tbl)), 3, pad = "0")))) %>% 
  mutate(date = quatri_vect) %>% 
  select(date, everything()) %>% 
  write_csv(here::here("out_training/training_pts_classif.csv"))