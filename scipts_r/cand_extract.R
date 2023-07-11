## Candidate extraction

#### Packages ####
library(tidyverse)
library(here)

#### Functions ####
# Convert capital letters to numbers
LETTER2num <- function(x) {
  stringr::str_extract(x, "[:upper:]") %>%
    utf8ToInt(.) - utf8ToInt("A") + 1L
}
# Takes a vector of microtiter wells and converts it to a number-based format
# with two leading zeros. Also appends the number of well positions specified
# and the ".flex" file suffix. The format converts from "A1" to "00100100X.flex".
well_to_flex <- function(wells_vector, nrep) {
  map(wells_vector, function(well) {
    stringr::str_c(
      stringr::str_pad(LETTER2num(well), 3, pad = "0"),
      stringr::str_pad(stringr::str_extract(well, "[:digit:]+"), 3, pad = "0")
    )
  }) %>%
    # Adds 'nrep' replicate files (i.e. the well positions)
    purrr::map(function(convert) {
      purrr::map(seq(1:nrep), function(rep) {
        stringr::str_c(
          convert,
          stringr::str_pad(rep, 3, pad = "0"),
          ".flex"
        )
      })
    })
}

# Copy raw images of interest to output directory for manual inspection
# Requires data frame with "plate", "rep", "iso_id", and "flex"-columns,
# all as chr vectors. Plate corresponds to first directory layer, rep to second.
# Flex corresponds to raw image file names with format "00X00X00X.flex"
extract_imgs <- function(input_dir, output_dir, sig_df) {
  dirs <- list.dirs(input_dir, recursive = F, full.names = F)
  if (!dir.exists(here(output_dir))) {
    dir.create(here(output_dir))
  }

  map(dirs, function(dir) {
    map(unique(sig_df$plate), function(plate) {
      if (grepl(plate, dir)) {
        dir.create(here(output_dir, dir))

        subdirs <- list.dirs(here(input_dir, dir),
          recursive = F, full.names = F
        )

        map(subdirs, function(subdir) {
          raw_dir <- list.dirs(here(input_dir, dir, subdir),
            recursive = F, full.names = F
          ) %>%
            .[. != "bins"]

          map(unique(sig_df$rep), function(rep) {
            if (stringr::str_detect(subdir, rep)) {
              output_subdirs <- unique(sig_df[sig_df$plate == plate & sig_df$rep == rep, ]$iso_id)
              map(output_subdirs, function(output_subdir) {
                if (!dir.exists(here(output_dir, dir, output_subdir))) {
                  dir.create(here(output_dir, dir, output_subdir))
                }

                file.copy(
                  paste0(here(input_dir, dir, subdir, raw_dir), "/", sig_df[sig_df$plate == plate & sig_df$rep == rep & sig_df$iso_id == output_subdir, ]$flex),
                  paste0(here(output_dir, dir, output_subdir), "/", rep, "_", sig_df[sig_df$plate == plate & sig_df$rep == rep & sig_df$iso_id == output_subdir, ]$flex)
                )
              })
            }
          })
        })
      }
    })
  })
}


#### Global variables ####
df <- read_csv(here("data_processed/l2fc.csv"), col_names = T)

# Add missing directories
#if (!dir.exists(here("figs", figfolder))) dir.create(here("figs", figfolder))

#### Processing ####
# Filter for top and bottom candidates
group_list <- df %>% 
  filter(culture == "dual", plate != "jena") %>% 
  group_by(ch, medium) %>% 
  arrange(log2fc) %>% 
  group_split()

split_list <- map(group_list, function(x) {
  x %>% 
    group_by(iso_id) %>% 
    summarise(log2fc = mean(log2fc), .groups = "drop") 
})

split_list_up <- map(split_list, function(x) filter(x, log2fc > 0))
split_list_down <- map(split_list, function(x) filter(x, log2fc < 0))

up <- semi_join(split_list_up[[1]], split_list_up[[2]], by = "iso_id") %>% 
  semi_join(split_list_up[[4]], by = "iso_id") %>% 
  semi_join(split_list_up[[3]], by = "iso_id")

down <- semi_join(split_list_down[[3]], split_list_down[[4]], by = "iso_id") %>% 
  semi_join(split_list_down[[2]], by = "iso_id") %>% 
  semi_join(split_list_down[[1]], by = "iso_id")


df_wide <- df %>% 
  mutate(group = str_c(medium, ch, sep = "_"))%>% 
  filter(culture == "dual", plate != "jena") %>% 
  group_by(iso_id, group) %>% 
  summarise(log2fc = mean(log2fc), .groups = "drop") %>% 
  select(iso_id, log2fc, group) %>% 
  spread(key = group, value = log2fc)


# Score isolates based on mean and sd
df_wide_scored <- df_wide %>% 
  mutate(mean = rowMeans(df_wide[-1]), sd = apply(df_wide[-1], 1, sd)) %>% 
  arrange(desc(mean), sd)

# Scoring could be better. Check Magnus' suggestion
df_wide_scored$iso_id <- fct_reorder(as.factor(df_wide_scored$iso_id), df_wide_scored$mean) %>% 
  fct_reorder(desc(df_wide_scored$sd))

df_wide_scored %>% 
  select(-mean, -sd) %>% 
  gather("group", "log2fc", -iso_id) %>% 
  ggplot(aes(x = group, y = iso_id, fill = log2fc)) +
  geom_tile() +
  scale_fill_continuous(type = "viridis")

ggsave(filename = "heatmap.png", 
       device = "png", 
       path = here("figs", figfolder), 
       width = 8,
       height = 30,
       units = "in", 
       dpi = 300, 
       scale = 1.5)


#Extract top and bottom 25
plate_info <- df %>% 
  select(iso_id, plate, well) %>% 
  group_by(iso_id, plate, well) %>% 
  distinct(iso_id)

top <- df_wide_scored %>% 
  slice_head(n = 31) %>% 
  semi_join(up, by = "iso_id") 

bot <- df_wide_scored %>% 
  slice_tail(n = 26) %>% 
  semi_join(down, by = "iso_id")


top %>%
  left_join(plate_info, by = c("iso_id")) %>% 
  write_csv(here("data_processed", "top25_allcond.csv"))

bot %>%
  left_join(plate_info, by = c("iso_id")) %>% 
  write_csv(bot, here("data_processed", "bot25_allcond.csv"))


#Extract raw images
# top_flex <- top %>%
#   left_join(df, by = c("iso_id")) %>%
#   select(iso_id, plate, rep, well, medium, ch, log2fc) %>%
#   ungroup() %>%
#   mutate(flex = well_to_flex(well, nrep = 4)) %>%
#   unnest(cols = flex) %>%
#   mutate_at("flex", as.character) %>%
#   mutate(rep = paste0("rep", rep))
# 
# bot_flex <- bot %>%
#   left_join(df, by = c("iso_id")) %>%
#   select(iso_id, plate, rep, well, medium, ch, log2fc) %>%
#   ungroup() %>%
#   mutate(flex = well_to_flex(well, nrep = 4)) %>%
#   unnest(cols = flex) %>%
#   mutate_at("flex", as.character) %>%
#   mutate(rep = paste0("rep", rep))
# 
# 
# raw_image_dirs <- list.dirs(here("images_opera"), recursive = F)
# 
# map(raw_image_dirs, function(x) {
#   dir_name <- str_extract(x, "(?<=images_opera/).*")
#   if (!dir.exists(here("data_processed", "top25", dir_name))) dir.create(here("data_processed", "top25", dir_name), recursive = T)
#   extract_imgs(x, here("data_processed", "top25", dir_name), top_flex)
# })
# 
# map(raw_image_dirs, function(x) {
#   dir_name <- str_extract(x, "(?<=images_opera/).*")
#   if (!dir.exists(here("data_processed", "bot25", dir_name))) dir.create(here("data_processed", "bot25", dir_name), recursive = T)
#   extract_imgs(x, here("data_processed", "bot25", dir_name), bot_flex)
# })
