# Change OPERA file names
# dir is your directory with files you want to rename, or parent directory with multiple other directories.
# suffix is file format of the files you want to rename. Use ".*" if you have different file formats but want to rename them all.
library(tidyverse)
library(here)

input = here("images_opera", "test")
suffix = "flex"
rename_opera_bins <- function(input, suffix = "tif") {
  bin_dirs <- list.dirs(input, recursive = T) %>%
    grep("bins$", ., value = TRUE)

  file_list <- lapply(bin_dirs, list.files)

  org <- mapply(function(dir, file) paste0(dir, "/", file),
    dir = bin_dirs, file = file_list,
    SIMPLIFY = F
  )

  files_idx <- lapply(file_list, function(x) stringr::str_extract_all(x, "0.[:digit:]", simplify = T)) %>%
    t()
  files_idx <- lapply(files_idx, function(x) {
    apply(x, 2, function(x) stringr::str_replace_all(x, "^0{1,}", "")) %>%
      as.data.frame() %>%
      dplyr::mutate(V1 = LETTERS[as.numeric(V1)])
  })

  # Requires a parent directory called "rep.[:digit:]". Returns "NA" as "pos", otherwise
  # TODO: Perhaps make a work-around?
  names <- mapply(function(dir, idx) {
    paste0(idx[, 1], stringr::str_pad(idx[, 2], 2, pad = "0"),
           "_rep", idx[, 3],
           "_", gsub("rep *", "pos", stringr::str_extract(dir, "rep.*[:digit:]")))
  }, dir = bin_dirs, idx = files_idx, SIMPLIFY = F)

  # Channel and time are being separated from their original file names.
  # Might cause problems in non-standard datasets.
  ch_t <- lapply(file_list, function(dir) {
    as_tibble(stringr::str_split(dir, "_", simplify = T)) %>%
      mutate(
        V2 = gsub("C", "ch", V2),
        V3 = paste0("time", rep(seq(length(dir) / 2), each = 2))
      ) %>%
      select(V2, V3) %>%
      apply(1, function(x) paste(x, collapse = "_"))
  })

  new <- mapply(function(name, ch_t) {
    paste0(name, "_", ch_t, ".", suffix)
  }, name = names, ch_t = ch_t, SIMPLIFY = F)

  mapply(function(org, new, name) file.rename(org, paste0(name, "/", new)),
    org = org, new = new, name = names(new)
  )
}
  
rename_opera_bins(here(input))


  
numbers_to_letters <- function(dir, suffix = "tif") {
  files <- list.files(path = dir, pattern = paste0("*.", suffix), full.names = F)

  files_idx <- sapply(files, function(x) str_extract_all(x, "0.[:digit:]", simplify = T)) %>%
    t()
  files_idx <- apply(files_idx, 2, function(x) str_replace_all(x, "^0{1,}", "")) %>% 
    as.data.frame() %>%
    mutate_at("V1", funs(LETTERS[as.numeric(.)]))
  
  names <- paste0(dir, files_idx[, 1], files_idx[, 2], "_rep", files_idx[, 3], ".", suffix)
  
  file.rename(paste0(dir, files), names)
}

fix_for_bq <- function(dir, suffix = "tif") {
  org <- list.files(path = dir, pattern = paste0(suffix))
  org <- stringr::str_sort(org, numeric = TRUE)
  ## Define new file names
  if(!any(grepl("pos", org))){
    new <- paste0("pos1_", org)
    } else {
      print("'pos' is already present in some filenames. Remove 'pos' from all filenames, before adding to the rest")
      new <- org
    }

  new <- as_tibble(stringr::str_split(new, "_", simplify = T))
  new <- mutate(new, V4 = gsub("C", "ch", new$V4))
  new[, 5] <- paste0("T", rep(seq(length(org)/2), each = 2), ".", suffix) 
  new <- apply(new, 1, function(x) paste(x, collapse = "_"))
  
  # Rename files
  file.rename(paste0(dir, org), paste0(dir, new))
}

fix_sp8_names <- function(dir, suffix = "tif") {
  org <- list.files(dir, pattern = paste0(suffix))
  new <- org %>% 
    stringr::str_replace_all(c("- " = "", ".tif" = "")) %>% 
    stringr::str_split(" ", simplify = T) %>% 
    as_tibble() %>% 
    select(meta = 3, ch = 4) %>% 
    mutate(ch = str_replace(ch, "C=", "ch"), 
           ch = ifelse(ch == "ch0", "ch1", "ch0"))
  
  meta <- separate(new, meta, c("scan", "well_L", "well_N", "rep"), sep = "_") %>% 
    select(-scan) %>% 
    mutate(rep = str_replace(rep, "R", "rep")) %>% 
    unite(well, c("well_L", "well_N"), sep = "") %>% 
    unite(fname, c("well", "rep", "ch"), sep = "_") %>% 
    mutate(fname = paste0("pos1_", fname, "_T", rep(seq(length(org)/2), each = 2), ".", suffix))

  new <- t(as.vector(meta))

  file.rename(paste0(dir, "/", org), paste0(dir, "/", new))
}


# fix_sp8_names(here("images_sp8", "mKateEPS_rep2_210310", "tifs"))
# numbers_to_letters("C:/Users/Mark Lyng/OneDrive - Danmarks Tekniske Universitet/PhD DTU/Projects/Multispecies interaction screen/images_opera/p5_pl3/mKateEPS_rep2_210326/opera_raw/")
# fix_for_bq("C:/Users/Mark Lyng/OneDrive - Danmarks Tekniske Universitet/PhD DTU/Projects/Multispecies interaction screen/images_opera/p5_pl3/mKateEPS_rep3_210407/tifs/")
# 
