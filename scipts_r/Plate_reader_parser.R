#' ---
#' PLATE READER FUNCTIONS
#'
#' Title: Plate reader parser
#' Author: Mark Lyng
#' Aim: Convert raw plate reader files to long-format data that can be plotted with ggplot
#' Necessities: A raw plate reader file, a layout file with information about the wells
#'             (see layout examples - can be long or wide format),
#'             information about the data format (is it in a wide or long format?
#'             Time series data from the Synergy HTX - Jerry - is always long).
#'             dplyr (library(tidyverse) or library(dplyr)).
#'             Source this script to your data analysis file with source("plate_reader_parser.R").
#' Date: 29/10-2020
#' Update 07/11-2020 - Refurbished comments to make code look nicer.
#' --

## This function will take your raw plate reader file as input and output a long-format data frame
## with structural metadata as indicated in your layout file.
parse_plate <- function(data, plate_layout, time_series = F, wide_data = T, wide_layout = T, time_end = NULL, measurement_identifier) {
  ## Check if layout has been provided in wide format. If yes, convert to long format
  if (wide_layout == T) {
    layout_names_idx <- grep(
      pattern = "[A-Za-z]{2}", # Looking for cells with at least two letters
      plate_layout[, 1]
    )
    wells <- as.character(plate_layout[grep(
      pattern = "^[A-Z]$",
      plate_layout[, 1]
    ), 1])

    ## Create data frame to store information
    ## A bit clunky, though it is able to handle all types of titer plates
    layout_long_test <- data.frame(well = c(paste0(LETTERS[rep(
      1:max(which(LETTERS[1:26] %in% wells)),
      each = ncol(plate_layout) - 1
    )],
    rep(1:(ncol(plate_layout) - 1),
      times = max(which(LETTERS[1:26] %in% wells))
    ))))

    ## The layout file block should start at the row right after the variable name
    for (i in seq_len(length(layout_names_idx))) {
      block_name <- plate_layout[layout_names_idx[i], 1]

      block_start <- layout_names_idx[i] + 1
      block_end <- layout_names_idx[i] + length(unique(wells))

      new_block <- plate_layout[(block_start):(block_end), 2:ncol(plate_layout)]
      new_block <- tidyr::pivot_longer(data = new_block, cols = names(new_block))

      layout_long_test[, block_name] <- new_block[, 2]
    }
    plate_layout <- layout_long_test
  }
  
  ## If data is end point, find measurement blocks from raw file information
  ## This only works for the tecan plate reader, so far. I don't know what raw files from the
  ## Synergy looks like yet, so I can't grep for anything.
  if (time_series == F) {
    start_time_idx <- which(data[, 1] == "Start Time:")
    end_idx <- which(data[, 1] == "End Time:")
    names_idx <- grep(pattern = "Label:.*", data[, 1])

    all_data <- c()
    for (i in seq_len(length(start_time_idx))) {
      block_name <- data[names_idx[i], 1]
      block_name <- gsub(pattern = "Label: ", replacement = "", block_name)

      block_start <- start_time_idx[i] + 4
      block_end_idx <- end_idx[i] - 5

      ## If data was provided in wide format, it must be reformatted into a long format
      if (wide_data == T) {
        new_block <- data[(block_start):(block_end_idx), 1:7]
        new_block <- tidyr::pivot_longer(data = new_block, cols = names(new_block[2:7]))
        new_block[, 2] <- rep(x = c(1:6), each = 4)
        new_block <- transform(new_block, well = paste(V1, variable, sep = ""))
        new_block <- cbind(new_block[4], new_block[3])

        ## Else if wide_data == F
      } else {
        new_block <- data[(block_start):(block_end_idx), 1:2]
        names(new_block)[1] <- "well"
        names(new_block)[2] <- "value"
      }

      ## Join parsed data to layout data to add layout variables
      new_block$value <- as.numeric(new_block$value)
      joined_block <- dplyr::full_join(plate_layout, new_block, by = "well")
      joined_block$measure <- block_name

      all_data <- rbind(all_data, joined_block)
    }
  }
  ## Else if time_series == T, adjust measurement blocks accordingly
  ## This only works for the Synergy HTX.
  ## I don't know what time series data from the tecan looks like, yet.
  else {
    names_idx <- which(grepl(pattern = measurement_identifier, data[, 1])) # Looking for identifiers like '1:...'
    all_data <- c()
    for (i in seq_len(length(names_idx))) {
      block_start <- names_idx[i] + 2
      
      if(is.null(time_end)) time_end <- "00:00:00"
      
      block_end <- 1 + names_idx[i] + min(which(grepl(x = data[block_start:nrow(data), 2], pattern = time_end))) # Changed to be more flexible       
      block_name <- sub(pattern = "Read [1-9]:", "", data[names_idx[i], 1])

      new_block <- data[(block_start):(block_end), 2:ncol(data)]
      names(new_block) <- new_block[1, ]
      new_block <- new_block[-1, ] %>% 
        select_if(function(x) !(all(is.na(x)) | all(x == ""))) 
      new_block <- new_block %>% 
        select(!starts_with("T°"))
      new_block <- tidyr::pivot_longer(data = new_block, cols = names(new_block[-1]), names_to = "well")

      ## Join parsed data to layout data to add layout variables
      joined_block <- dplyr::full_join(plate_layout, new_block, by = "well")
      joined_block$measure <- block_name
      all_data <- rbind(all_data, joined_block)
    }
  }

  if (any(is.na(all_data$value))) warning("You have missing values in your data. Are you sure you provided the right format?")
  return(all_data)
} # End function

#' ---
#' Subtract non-flu from flu
#' Aim: Subtract background fluorescence from fluorescent samples using non-fluorescent controls
#' Necessities: This function requires a tidy df which can be made with the plate_reader_parser function. The df must have columns named 'sample' and 'type'
#'             'sample' contains your sample names.
#'             'type' contains either WT (non-fluorescent) or FLU (fluorescent).
#'             You can provide additional grouping variables as a character vector in the 'groups' argument.
#' --

subtr_nonflu <- function(df, groups = c("measure", "time", "medium")) {
  ## Merge the values of the non-fluorescent control samples to the rows of the fluorescent samples
  ## and calculate the difference.
  df_norm <- df %>%
    filter(type == "WT") %>%
    group_by(sample, !!!syms({{ groups }})) %>%
    summarize(mean_ctrl = mean(value)) %>%
    full_join(., df %>%
      filter(type == "FLU")) %>%
    filter(measure != "OD", !is.na(mean_ctrl)) %>% # We only want fluorescent data and no NA-vals. "OD" potential problem?
    mutate(value = value - mean_ctrl)

  return(df_norm)
}


#' ---
#' Calculate plate reader log2fc
#' Aim: Calculate the difference in fluorescence between two related samples (ss and ds) as a log(FC)
#' Necessities: This function requires a tidy df which can be made with the plate_reader_parser function. The df must have columns named "sample" and
#'             "GFP_partner" and/or "RFP_partner".
#'             "sample" contains your sample names.
#'             "rep" contains the sample replicate.
#'             The samples to combine and compare are given as ss_GFP/ss_RFP and ds (character vectors of length 1 or more).
#'             Finally, you can provide additional grouping variables used to calculate mean values and group samples correctly when calculating log2fc.
#' Date: 06/11-2020
#' --

## Function to calculate log2(FC) between two types of fluorescent samples
pltrd_log2fc <- function(df, ss_GFP, ss_RFP, ds, groups = c("measure", "time", "medium", "rep", "type")) {
  ## Calculate log2(FC) for dual species fluorescence compared to single species
  if (has_name(df, "ss_partner_GFP")) {
    GFP <- df %>%
      filter(sample %in% ss_GFP, measure == "GFP") %>%
      full_join(., df %>%
        filter(sample %in% ds, measure == "GFP"),
      by = c({{ groups }}, {{ "ss_partner_GFP" }})
      ) %>%
      mutate(log2fc = log2(value.y / value.x)) %>%
      mutate(abs_log = log2((abs(value.y / value.x)))) # Abs value might not be the way to go. Doesn't reflect distance between vals.
  }

  ## Repeat above for RFP samples
  if (has_name(df, "ss_partner_RFP")) {
    RFP <- df %>%
      filter(sample %in% ss_RFP, measure == "RFP") %>%
      full_join(., df %>%
        filter(sample %in% ds, measure == "RFP"),
      by = c({{ groups }}, {{ "ss_partner_RFP" }})
      ) %>%
      mutate(log2fc = log2(value.y / value.x)) %>%
      mutate(abs_log = log2((abs(value.y / value.x))))
  }

  ## Bad way to test what we've done above. Combines data or just present output.
  if (all(has_name(df, c("ss_partner_RFP", "ss_partner_GFP")))) {
    df_norm_log2fc <- bind_rows(GFP, RFP) %>%
      select(-starts_with("ss"), -type, -starts_with("mean"), -starts_with("well")) %>%
      rename(ss = sample.x, ds = sample.y, value_ss = value.x, value_ds = value.y)
  } else if (has_name(df, c("ss_partner_GFP"))) {
    df_norm_log2fc <- GFP %>%
      select(-starts_with("ss"), -type, -starts_with("mean"), -starts_with("well")) %>%
      rename(ss = sample.x, ds = sample.y, value_ss = value.x, value_ds = value.y)
  } else if (has_name(df, c("ss_partner_RFP"))) {
    df_norm_log2fc <- RFP %>%
      select(-starts_with("ss"), -type, -starts_with("mean"), -starts_with("well")) %>%
      rename(ss = sample.x, ds = sample.y, value_ss = value.x, value_ds = value.y)
  }
  return(df_norm_log2fc)
}
