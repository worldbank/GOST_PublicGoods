library(dplyr)
library(fuzzyjoin)
library(readxl)

lbr_io_pth <- 'P:/Data/PROJECTS/Health/healthsitesio/lbr_io.csv'
lbr_ml_pth <- 'P:/Data/PROJECTS/Health/from_tashrik/master lists/Liberia.xlsx'
lbr_hl_pth <- 'P:/Data/PROJECTS/Health/from_tashrik/hierarchy/raw_liberia_collapsed.csv'

lbr_io <- read.csv(lbr_io_pth)
lbr_ml <- read_xlsx(lbr_ml_pth, 1) %>%
  rename(
    hf_name = "HF Name"
  )
lbr_hl <- read.csv(lbr_hl_pth)

# HF Name

lbr.join <- stringdist_full_join(
  lbr_ml, 
  lbr_io[c("name", "geometry")], 
  by=c('hf_name'='name'),
  distance_col='dist',
  max_dist=50
) %>%
  arrange(dist) %>%
  group_by(hf_name) %>%
  # head(2) %>%
  filter(row_number()==1) %>%
  ungroup()

View(lbr.join[c("County", "District", "Location", "New HF ID", "hf_name", "name", "dist")])