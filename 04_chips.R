
### Script to generate PNGs of time series chips for imported training points
# NOTE: code written for QUATRI
#
# *************************************************************************************************


# Load libraries
library(raster)
library(RStoolbox) # libgeos++-dev
library(tidyverse)
library(lubridate)
library(sf)
library(patchwork)
library(furrr)
library(tictoc)
library(velox)

# External directory (input/output)
ext_dir <- "/home/cdobler/Documents/calakmul_landsat_analysis/"

# Load training points
training_pts <- read_sf(here::here("data/training_pts_lg/pts.gpkg")) # *** CHANGE FILE ***

# Import ejidos shape
ejido <- read_sf("/media/cdobler/Neobxbaumia/Research/Data_spatial/Shapefiles/03 SYPR shapefiles/Ejidos/Ejidos_select.shp") %>% 
  st_transform(st_crs(training_pts)) %>%
  select("NOM_NA") %>% 
  slice(3) # ejido of interest    *** CHANGE ***

# Get list of scenes
list_of_files <- str_c(ext_dir, "03_composites/quatri/") %>%
  list.files() %>% 
  grep(pattern = "aux", inv = T, value = T)

# Reference raster
r <- raster(str_c(ext_dir, "03_composites/quatri/", list_of_files[1])) %>% raster()


# *************************************************************************************************

# Function to loop through scenes
func_chips_scenes <- function(sc, chanR = 3, chanG = 2, chanB = 1) {   # Change defaults
  
  # print(str_c("Processing scene ", sc, " / ", length(list_of_files)))
  
  rast <- stack(str_c(ext_dir, "03_composites/quatri/", list_of_files[sc])) %>% 
    crop(extent(ejido) + 60)
  
  # Get date
  date <- list_of_files[sc] %>%
    str_sub(end = -5) %>% 
    str_split("_", simplify = T) %>% 
    {str_c(.[3], "_", .[4])}
  
  yr <- date %>% str_sub(1, 4) %>% as.integer()
  
  # Generate bands quantile tibble
  lim_tbl <- map_df(seq_len(nlayers(rast)), function(x) {
    lower <- quantile(raster::values(rast[[x]]), 0.02, na.rm = T)
    upper <- quantile(raster::values(rast[[x]]), 0.98, na.rm = T)
    bind_cols(band = x, lower = lower, upper = upper)
  })
  
  # Matrix for ggRGB
  lim_matrix <- lim_tbl %>% 
    mutate(band_order = case_when(band == chanR ~ 1,
                                  band == chanG ~ 2,
                                  band == chanB ~ 3)) %>% 
    filter(!is.na(band_order)) %>% 
    arrange(band_order) %>% 
    dplyr::select(lower, upper) %>% 
    as.matrix()
  
  # Loop through points
  tbl_pts <- seq_len(nrow(training_pts)) %>% 
    map(function(pt){
      
      # Extent of chip (16 cells on each direction)
      e <- tabularaster::cellnumbers(r, training_pts[pt,]) %>% 
        .$cell_ %>% 
        rowColFromCell(r, .) %>% as_tibble() %>% 
        {c(.$row+15, .$row-15, .$col-15, .$col+15)} %>% 
        {c((xFromCol(r, .[3])-15), 
           (xFromCol(r, .[4])+15), 
           (yFromRow(r, .[1])-15), 
           (yFromRow(r, .[2])+15))} %>% 
        extent(.)
      
      # Number of valid pixels in chip
      valid <- raster::extract(rast[[1]], e) %>% 
      {!is.na(.)} %>%
        as.numeric() %>% 
        sum() %>% 
        {./961}   # 961: number of clear cells in the chip
      
      # Conditionals to process a chip
      # if (valid > 0.9 | 
      #     (yr > 1989 & yr < 1993 & valid > 0.7) |
      #     (yr > 2003 & yr < 2010 & valid > 0.7)) {
      
      if (valid >= 0.8) {
        
        ggRGB(rast, r = chanR, g = chanG, b = chanB, 
              ext = e, 
              limits = lim_matrix, 
              stretch = "lin", 
              ggObj = F) %>% 
          mutate(date = date,
                 point = pt)
        
      } # end of if
      
    }) %>% # end of map
    
    bind_rows()
  
  return(tbl_pts)
  
} # end of function


# *************************************************************************************************
# Obtain massive tibble with data to display chips (for all points)

plan(multicore, workers = 6)

Sys.time()
tic()
chips_tbl <- future_map(seq_along(list_of_files),
                        func_chips_scenes, chanR = 4, chanG = 5, chanB = 3)
toc() # ~ 3 min

chips_tbl <- chips_tbl %>% bind_rows()

# Save
saveRDS(chips_tbl, here::here("output/chips_tibble_comp_non_453_lg.RDS")) # *** CHANGE NAME ***


# *************************************************************************************************
# Generate chips for each point with ndmi dynamics added

chips_tbl <- readRDS(here::here("output/chips_tibble_comp_non_453_lg.RDS")) # *** CHANGE NAME ***

dates <- str_c(ext_dir, "04_ndmi/") %>% 
  list.files() %>% 
  grep(pattern = "aux", inv = T, value = T) %>% 
  str_sub(6, 15)

vi_stack <- str_c(ext_dir, "04_ndmi/") %>% 
  list.files(full.names = T) %>% 
  grep(pattern = "aux", inv = T, value = T) %>% 
  stack() 

Sys.time()
tic()
vi_stack <- vi_stack %>% 
  crop(extent(training_pts) + 90)
toc() # ~19.6 min

Sys.time()
tic()
vi_stack <- velox(vi_stack)
toc() # ~ 8 sec

# Loop for chips displays (one per point)
plan(multicore, workers = 6)

Sys.time()
tic()
future_map(seq_len(nrow(training_pts)),
           function(pt){
             
             print(pt)
             
             # Focus polygon (3x3) as tidy for ggplot
             pix <- tabularaster::cellnumbers(r, training_pts[pt,]) %>%
               .$cell_ %>%
               rowColFromCell(r, .) %>% as_tibble() %>%
               {c(.$row+1, .$row-1, .$col-1, .$col+1)} %>%
               {c((xFromCol(r, .[3])-15),
                  (xFromCol(r, .[4])+15),
                  (yFromRow(r, .[1])-15),
                  (yFromRow(r, .[2])+15))} %>%
               extent(.) %>%
               as(., 'SpatialPolygons') %>%
               broom::tidy()
             
             if (nrow(filter(chips_tbl, point == pt)) != 0) {
               
               # Chips
               chip <- chips_tbl %>%
                 filter(point == pt) %>%
                 ggplot(aes(x = x, y = y)) +
                 geom_raster(aes(fill = fill)) +
                 scale_fill_identity() +
                 geom_polygon(aes(x = long, y = lat),
                              data = pix,
                              fill = NA,
                              color = "black",
                              size = 0.3) +
                 coord_equal()+
                 theme_void() +
                 facet_wrap(vars(date), ncol = 14)
               
               # VI dynamics
               pt_analysis <- tabularaster::cellnumbers(r, training_pts[pt,]) %>%
                 .$cell_ %>%
                 xyFromCell(r, .) %>%
                 as_tibble() %>%
                 st_as_sf(., coords = c("x", "y"))
               
               tb_analysis <- tibble(
                 date = dates,
                 vi = vi_stack$extract_points(pt_analysis)[1, ]) %>% 
                 mutate(year = as.integer(str_split(date, "_", simplify = T)[ ,1]),
                        month = str_split(date, "_", simplify = T)[ ,2],
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
                 right_join(tibble(date = str_c(rep(1985:2018, each = 4), "_", c(1, 2, 3, 4)))) %>% 
                 mutate(year = as.integer(str_split(date, "_", simplify = T)[, 1]),
                        qtr = as.integer(str_split(date, "_", simplify = T)[, 2]),
                        date_num = year + (qtr/4) - (1/4),
                        pt = pt)
               
               gg <- ggplot(tb_analysis, aes(x = date_num, y = vi)) +
                 geom_point() +
                 scale_x_continuous(breaks = seq(1985, 2019, 3), minor_breaks = 1985:2019) +
                 #scale_x_date(breaks = "3 years", minor_breaks = "1 year", labels = scales::date_format("%Y")) +
                 ylim(-300, 600) +
                 theme(axis.title.x = element_blank()) +
                 labs(title = str_c("Point ", pt),
                      subtitle = "")
               
               # Patchwork
               p <- gg + chip + plot_layout(ncol = 1, heights = c(1, 4))
               
               # Save plot *** CHANGE NAME ***
               ggsave(here::here("output/chips_cuatri_453_lg/", str_c("chips_pt_", str_pad(pt, 2, pad = "0"), ".png")),
                      plot = p,
                      width = 16,
                      height = 10)
               
             } # end of if
             
             })
toc() # ~ 1.4 min


# *************************************************************************************************
# Generate the cvs to classify manually

ini_tb <- expand.grid(yr = 1985:2018, cuatri = c("1st", "2nd", "3rd", "4th")) %>% 
  mutate(date = str_c(yr, "_", cuatri)) %>% 
  select(date)

map_dfr(unique(chips_tbl$point), 
       function(pt){
         
         chips_tbl %>% 
           filter(point == pt) %>% 
           {tibble(date = unique(.$date),
                   point = str_pad(pt, 2, "left", "0"))}
         
       }) %>% 
  mutate(cover = "u") %>% 
  spread(key = point, value = cover) %>% 
  write_csv(here::here("output/table_classif_cuatri_lg.csv")) # *** CHANGE NAME ***
    
