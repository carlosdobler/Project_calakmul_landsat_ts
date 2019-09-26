
# Script to: 
#     * mask SR bands
#     * create stacks of single bands (to run timeSyncR)


library(tidyverse)
library(lubridate)
library(raster)


dir_qa <- "/media/cdobler/PSEUDBOMBAX/tmp_processing/qa_all_layers/"
files_qa <- grep(list.files(dir_qa), pattern = "aux", inv = T, value = T)

dir_sr <- "/media/cdobler/PSEUDBOMBAX/tmp_processing/sr_all_layers/"
files_sr <- grep(list.files(dir_sr), pattern = "aux", inv = T, value = T)


# Build stacks ------------------------------------------------------------------------------------

for (b in 3:7) {
  
  date_vector <- vector("character", length(files_sr))
  
  for (i in 1:length(files_sr)) {
    
    year <- str_split(files_sr[i], "_", simplify = T)[2]
    month <- str_split(files_sr[i], "_", simplify = T)[3]
    day <- str_split(files_sr[i], "_", simplify = T)[4]
    
    print(str_c("Processing ", i, "/", length(files_sr), " of band ", b, " --- ", year, "-", month, "-", day))
    
    qa_mask <- raster(str_c(dir_qa, files_qa[i]))
    
    band <- raster(str_c(dir_sr, files_sr[i]), band = b) %>% 
      overlay(., qa_mask, fun = function(a, b) {
        ifelse(b == 1, a, NA)
      })
    
    writeRaster(band, 
                str_c("/media/cdobler/PSEUDBOMBAX/tmp_processing/sr_band_stacks/bnd_", b, "/sr_", year, "_", month, "_", day, "_bnd_", b, ".tif"),
                datatype = "INT2S")
                
    date_vector[i] <- str_c(year, "-", month, "-", day)
    
  }
  
  # dir_band <- "/media/cdobler/PSEUDBOMBAX/tmp_processing/temp_files/"
  # file_band <- grep(list.files(dir_band, full.names = T),
  #                              pattern = "aux", inv = T, value = T)
  # 
  # band_stack <- stack(file_band)
  # 
  # writeRaster(band_stack,
  #             str_c("/media/cdobler/PSEUDBOMBAX/tmp_processing/sr_band_stacks/bnd_", b, ".tif"),
  #             datatype = "INT2S",
  #             options="INTERLEAVE=BAND")
  # 
  # "/media/cdobler/PSEUDBOMBAX/tmp_processing/temp_files/" %>% 
  #   list.files(full.names = T) %>% 
  #   file.remove()
  # 
  # rm(band_stack)
  
}

saveRDS(date_vector, "/media/cdobler/PSEUDBOMBAX/tmp_processing/sr_band_stacks/date_vector.RDS")

