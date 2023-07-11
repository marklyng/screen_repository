## Categorizing screen isolates

#### Packages ####
library(tidyverse)
library(here)

df <- read_csv(here("data_processed/l2fc.csv"), col_names = T) %>% 
  mutate(plate = str_replace(str_to_upper(plate), "_", ""))

df_wide <- df %>%  
  mutate(group = str_c(medium, rep, sep = "_")) %>% 
  filter(culture == "dual", plate != "JENA", medium == "0.1x LBGM") %>% 
  group_by(iso_id, group, ch) %>%
  select(iso_id, log2fc, ch, group) %>% 
  spread(key = group, value = log2fc) %>% 
  ungroup()


# Score isolates based on mean and sd
# df_wide_scored <- df_wide %>% 
#   mutate(mean = rowMeans(df_wide[-1:-2], na.rm = T), sd = apply(df_wide[-1:-2], 1, sd, na.rm = T))

df_wide_median <- df_wide %>% 
  mutate(median = apply(df_wide[-1:-2], 1, median, na.rm = T), 
         iqr = apply(df_wide[-1:-2], 1, IQR, na.rm = T),
         q1 = apply(df_wide[-1:-2], 1, quantile, 0.25, na.rm = T),
         q3 = apply(df_wide[-1:-2], 1, quantile, 0.75, na.rm = T)) 

# categories <- list(
#   up = filter(df_wide_scored, mean > 1, sd < 0.5) %>% 
#     arrange(desc(mean), sd),
#   neutral = rbind(filter(df_wide_scored, !between(mean, -1, 1), sd > 0.5),
#                  filter(df_wide_scored, between(mean, -1, 1))
# ),
#   down = filter(df_wide_scored, mean < -1, sd < 0.5) %>% 
#     arrange(mean, desc(sd))
# )

categories_median <- list(
  up = filter(df_wide_median, q1 - 1.5 * iqr > 1),
  neutral = distinct(
    rbind(filter(df_wide_median, between(q1 - 1.5 * iqr, -1, 1)),
          filter(df_wide_median, between(q3 + 1.5 * iqr, -1, 1)),
          filter(df_wide_median, q1 - 1.5 * iqr < 1, q3 + 1.5 * iqr > 1)
    )
  ),
  down = filter(df_wide_median, q3 + 1.5 * iqr < -1)
)

# subcategories <- c(
#   map(categories, .f = function(x) {
#     filter(x, ch == "mKate")
#   }),
#   map(categories, .f = function(x) {
#     filter(x, ch == "eGFP")
#   })
# )

subcategories_med <- c(
  map(categories_median, .f = function(x) {
    filter(x, ch == "mKate")
  }),
  map(categories_median, .f = function(x) {
    filter(x, ch == "eGFP")
  })
)


subcategories_med <- map2(subcategories_med, names(subcategories_med), function(x, name){
  mutate(x, cat = name)
}) %>% 
  do.call("rbind", .)

cat_list_med <- left_join(distinct(select(filter(df, culture == "dual", plate != "JENA"), 
                                  plate, well, iso_id)),
          select(subcategories_med, iso_id, ch, cat),
          by = "iso_id") %>% 
  pivot_wider(names_from = ch, values_from = cat)


# subcategories <- map2(subcategories, names(subcategories), function(x, name){
#   mutate(x, cat = name)
# }) %>% 
#   do.call("rbind", .)
# 
# cat_list <- left_join(distinct(select(filter(df, culture == "dual", plate != "JENA"), 
#                                           plate, well, iso_id)),
#                           select(subcategories, iso_id, ch, cat),
#                           by = "iso_id") %>% 
#   pivot_wider(names_from = ch, values_from = cat)



cat_list %>% 
  group_by(mKate) %>% 
  tally()

cat_list %>% 
  filter(iso_id == "p5_109")

cat_list %>% 
  filter(mKate == "up") %>% 
  write_csv(here("data/mkate_up.csv"))


#write_csv(cat_list, here("data/isolate_categories_median.csv"))
