## Data processing (l2fc calc)

#### Packages ####
library(tidyverse)
library(here)

#### Functions ####

# Calculates biovolume fold change between coculture and the mean of all mono samples
# in each channel.
# log2(µ_dual / µ_mono) = log2(µ_dual) - log2(µ_mono)
# i.e. high log2fc: increase in dual
calc_log2fc <- function(df, half_mono_mean = F) {
  if (half_mono_mean == T) {
    # DIVIDING MONO_MEAN BY 2!
    df <- df %>%
      group_by(ch) %>%
      filter(culture == "mono") %>%
      summarise(mono_mean = mean(Biofilm_Volume, na.rm = T) / 2, .groups = "drop") %>%
      select(mono_mean, ch) %>%
      full_join(df, by = "ch")
  } else {
    # Not dividing mono_mean by 2
    df <- df %>%
      group_by(ch) %>%
      filter(culture == "mono") %>%
      summarise(mono_mean = mean(Biofilm_Volume, na.rm = T), .groups = "drop") %>%
      select(mono_mean, ch) %>%
      full_join(df, by = "ch")
  }
  
  df <- df %>%
    group_by(iso_id, well, ch, rep) %>%
    # Transforming before averaging - i.e. mean(log2FC) NOT log2(mean(FC))
    mutate(log2fc = mean(log2((Biofilm_Volume + 1) / mono_mean)))
}

#### Global variables ####
l2fc_thresh <- 1
half_mono_mean <- T
bq_dir <- list(here("images_opera_bq", "TSB"), here("images_opera_bq", "0.1X_LBGM"))
data_in_dir <- list(here("data_bq_raw", "mkate_eps_TSB"), here("data_bq_raw", "mkate_eps_LBGM"))

# Add missing directories
#if (!dir.exists(here("figs", figfolder))) dir.create(here("figs", figfolder))

#### Read modified raw data ####
# Read metadata files. Double check that these are good!
md_files <- lapply(bq_dir, function(x) list.files(x, "_md", recursive = T))
layouts <- map2(md_files, bq_dir, function(x, dir) {
  lapply(x, function(file) {
    read_csv2(paste0(dir, "/", file),
              col_types = cols()
    )
  })
}) %>% 
  unlist(., recursive = F)

layouts <- lapply(layouts, function(x) {
  letter <- str_extract(x$well, pattern = "[[:upper:]]")
  number <- str_extract(x$well, pattern = "[[:digit:]]+") %>%
    str_pad(width = 2, pad = "0")
  comb <- str_c(letter, number)
  
  x <- x %>%
    mutate(well = comb)
  #   arrange(well)
})

# Read previously exported data files.
df_files <- lapply(data_in_dir, function(x) {
  lapply(x, function(y) list.files(y, "*.csv", recursive = T))
}) %>% 
  unlist(., recursive = F)

df_list <- map2(df_files, data_in_dir, function(files, dir) {
  lapply(files, function(file) {
    read_csv(paste0(dir, "/", file),
             col_types = cols()
    )
  })
})

#### Transformation ####
# FIXED IN LBGM DATA
# Temporary solution to changing wrong col names
# run and well_pos are wrong in loaded data
df_list[[1]] <- lapply(df_list[[1]], function(x) rename(x, rep = well_pos, well_pos = run))

# Temporary solution to pad wells (i.e. turn A1 to A01)
# This is fixed already in change_biofilm_filenames.R
df_list <- lapply(df_list, function(dfs) {
  lapply(dfs, function(x) {
    letter <- str_extract(x$well, pattern = "[[:upper:]]")
    number <- str_extract(x$well, pattern = "[[:digit:]]+") %>%
      str_pad(width = 2, pad = "0")
    comb <- str_c(letter, number)
    
    x <- x %>%
      mutate(well = comb) %>%
      arrange(well)
  })
})


# Average the biovolume per well per ch and calculate log2fc
# Doesn't require "group_by(rep)", because I am performing it on list-level (but I'm including it to keep the variable).
# i.e. calculates mono_mean per plate.

# Adding plate to data
l2fc_list <- map2(df_list, data_in_dir, function(dfs, dir) {
  map2(dfs, rep(list.dirs(dir, recursive = F, full.names = F), each = 3),
       function(.x, .y) mutate(.x, plate = .y)
  )
}) %>% 
  unlist(., recursive = F)

l2fc_list <- map2(l2fc_list, layouts, function(.x, .y) {
  .x %>%
    #    filter(well != "H02" | plate != "P8PL3") %>% #VERY MANUAL!!!!!!!!!
    group_by(well, rep, ch, plate) %>%
    summarise(Biofilm_Volume = mean(Biofilm_Volume, na.rm = T), .groups = "drop") %>%
    # Join metadata from layout to corresponding plate.
    left_join(.y, by = "well") %>%
    calc_log2fc(., half_mono_mean = T) %>%
    mutate(threshold = ifelse(log2fc > l2fc_thresh, "above", "below"))
})

df <- do.call("bind_rows", l2fc_list) %>%
  mutate(ch = ifelse(ch == 1, "mKate", "eGFP")) %>%
  ungroup()

#write_csv(df, here("data_processed/l2fc.csv"))