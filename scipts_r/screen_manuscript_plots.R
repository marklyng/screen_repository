# Packages
pacman::p_load(here,
               tidyverse,
               patchwork,
               ggforce,
               vegan, 
               phyloseq,                                                                                        # Bioconductor
               decontam,                                                                                        # Bioconductor
               DESeq2,                                                                                          # Bioconductor
               ANCOMBC,                                                                                         # Bioconductor
               multcompView,
               ggpubr,
               ggrepel,
               ggtree)


#### Global variables ####
l2fc_thresh = 1

custom_palette = RColorBrewer::brewer.pal(8, "Set2")
brew_blue = RColorBrewer::brewer.pal(3, "Blues")[2]

custom_red = "#c7546e"
custom_grey = "#B3B3B3"
cust_mag = "#E78AC3"
cust_green = "#A6D854"
text_size = 12.5

#### Functions ####
getPalette = colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))

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
  if (!dir.exists(here(output_dir))) {dir.create(here(output_dir))}
  
  map(dirs, function(dir) {
    map(unique(sig_df$plate), function(plate) {
      if (grepl(plate, dir)) {
        dir.create(here(output_dir, dir))
        
        subdirs <- list.dirs(here(input_dir, dir),
                             recursive = F, full.names = F
        )
        
        map(subdirs, function(subdir) {
          raw_dir <- list.dirs(here(input_dir, dir, subdir),
                               recursive = F, full.names = F) %>% 
            .[. != "bins"]
          
          map(unique(sig_df$rep), function(rep) {
            if (stringr::str_detect(subdir, rep)) {
              output_subdirs <- unique(sig_df[sig_df$plate == plate & sig_df$rep == rep, ]$iso_id)
              map(output_subdirs, function(output_subdir) {
                if (!dir.exists(here(output_dir, dir, output_subdir))) {dir.create(here(output_dir, dir, output_subdir))}
                
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




#### Data import and wrangle ####
## BiofilmQ data
# Metadata
layouts <- map(list.files(here("metadata"), pattern = "_screen", full.names = T), function(medium) {
  map(list.files(here(medium), pattern = "md.csv", full.names = T), function(md) {
    read_csv2(md, show_col_types = F) %>% 
      mutate(well = str_c(str_extract(well, pattern = "[[:upper:]]"),
                          str_pad(str_extract(well, pattern = "[[:digit:]]+"),
                                  width = 2, pad = "0")
      )
      )
  }) %>% set_names(list.files(here(medium)))
}) %>% 
  set_names(list.files(here("metadata"), pattern = "_screen")) %>% 
  map(~rep(.x, each = 3))

# Data
data <- map(list.dirs(here("data_bq_raw"), full.names = T, recursive = F), function(medium) {
  map(list.files(here(medium), pattern = ".csv", full.names = T, recursive = T), function(df) {
    read_csv(df, show_col_types = F) %>% 
      mutate(well = str_c(str_extract(well, pattern = "[[:upper:]]"),
                          str_pad(str_extract(well, pattern = "[[:digit:]]+"),
                                  width = 2, pad = "0"))) %>% 
      arrange(well)
  }) %>% set_names(list.files(here(medium), pattern = ".csv", full.names = F, recursive = T))
}) %>% 
  set_names(list.dirs(here("data_bq_raw"), full.names = F, recursive = F))

data$mkate_eps_TSB <- map(data$mkate_eps_TSB, function(x) dplyr::rename(x, rep = well_pos, well_pos = run))     # Columns are wrong in TSB data

df <- map2(layouts, data, function(md_medium, df_medium) {
  map2(md_medium, df_medium, function(md, df) {
    df %>% 
      dplyr::select(Ncells, Biofilm_Volume, well, well_pos, rep, ch) %>% 
      left_join(md, by = "well")
  }) %>% bind_rows()
}) %>% bind_rows() %>% 
  mutate(ch = ifelse(ch == 1, "mKate", "eGFP"),
         rep = as_factor(rep))

# Calculate log2(co/mono), where the mean of the monoculture is divided by two, to take into account
# that only one species is present
df_l2fc <- df %>% 
  group_by(ch, rep, plate, medium) %>%
  filter(culture == "mono") %>%
  summarise(mono_mean = mean(Biofilm_Volume, na.rm = T) / 2, .groups = "drop") %>%
  right_join(df, by = c("ch", "plate", "rep", "medium")) %>% 
  group_by(ch, plate, medium, well, rep, sample, iso_id, culture, temp, time, mono_mean) %>%                    # Average measurements from multiple positions in one well
  summarise(biovolume = mean(Biofilm_Volume, na.rm = T),
            ncells = mean(Ncells, na.rm = T),
            .groups = "drop") %>%                                                                               # Calculate co/mono ratio and log2(ratio) of each well replicate
  mutate(ratio = (biovolume + 1) / mono_mean,                                                                   # +1 gets rid of zeros, and thus -Inf
         log2fc = log2(ratio))

# Summarize over replicates
df_median <- df_l2fc %>% 
  group_by(ch, plate, medium, iso_id, culture) %>% 
  summarise(bv_median = median(biovolume),
            log2fc_median = median(log2fc),
            q1 = quantile(log2fc, 0.25),
            q3 = quantile(log2fc, 0.75),
            iqr = IQR(log2fc),
            se = sd(log2fc)/n()) %>%
  mutate(threshold = ifelse(q1 - 1.5*iqr > l2fc_thresh, "positive",
                            ifelse(q3 + 1.5*iqr < -l2fc_thresh, "negative",
                                   "neutral")))


## rpoD Amplicon seq data
# Metadata
bcs <- str_c("bc", str_pad(seq_along(1:31), 2, side = "left", pad = "0"))

md <- read_csv2(here("metadata", "amplicon_md.csv")) %>% 
  filter(reporter != "gfp",
         !scr_cat %in% c("control_neg", "none"))

# Read bowtie2 data
tables <- map2(list.files(here("data", "ampliconseq", "abundance_tables"), full.names = T), bcs, function(x, bc) {
  read_table(x, 
             col_names = c(bc, "species"), 
             col_types = cols("i", "c"))
})

feature_table <- tables %>% 
  purrr::reduce(full_join, by = "species") %>% 
  mutate(across(where(is.numeric), ~ replace_na(.x, 0))) %>% 
  relocate(bc01, .after = species) %>% 
  mutate(otu = str_c("otu", seq(1:nrow(.)))) %>%
  dplyr::select(species,
                otu,
                bc13:bc28) # Remove GFP categories and negative controls

# Phyloseq object
tax_table <- feature_table %>%
  column_to_rownames("otu") %>% 
  select(species) %>% 
  separate(species, into = c("Genus", "Species", "Strain"), sep = "_", remove = T) %>% 
  select(-Strain) %>% 
  mutate(Kingdom = "Bacteria", 
         Phylum	= "Proteobacteria", 
         Class =	"Gammaproteobacteria", 
         Order = "Pseudomonadales", 
         Family = "Pseudomonadaceae") %>% 
  relocate(Genus, Species, .after = Family) %>% 
  as.matrix() %>% 
  tax_table()

otu_table <- feature_table %>% 
  select(-species) %>% 
  column_to_rownames("otu") %>% 
  otu_table(taxa_are_rows = T)

sample_data <- md %>%
  column_to_rownames("bc") %>% 
  sample_data()

ps <- phyloseq(otu_table, tax_table, sample_data)

# Normalize for sample depth with DESeq2
deseq_counts <- DESeqDataSetFromMatrix(column_to_rownames(select(feature_table, -species), "otu"), 
                                       colData = md,
                                       design = ~scr_cat) 
deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)

vst_feature_table <- assay(deseq_counts_vst)

vst_otu <- otu_table(vst_feature_table, 
                       taxa_are_rows = T)

vst_ps <- phyloseq(vst_otu, tax_table, sample_data)


## antiSMASH data
bgc <- readxl::read_xlsx(here("data/antismash", "bgc_presence.xlsx"))


#### Figure 1: Most Pseudomonads antagonize Bacillus ####
df_median %>% 
  left_join(bgc, by = "iso_id") %>% 
  filter(ch == "mKate",
         culture == "dual") %>%
  ggplot(aes(log2fc_median, medium)) +
  geom_violin(position = "identity",
              fill = "grey95") +
  geom_vline(xintercept = c(-1,1),
             linetype = "dashed") +
  geom_point(aes(fill = threshold), 
             position = position_jitter(seed = 1), 
             alpha = 0.5,
             pch = 21,
             size = 3) +
  # geom_label_repel(aes(label = ifelse(tax_method == "wgs", str_c(iso_id, "\n", ps_subgroup), NA),                 # Is this mostly for us?
  #                      fill = threshold),
  #                  position = position_jitter(seed = 1),
  #                  alpha = 0.85,
  #                  label.size = 0.01,
  #                  max.overlaps = 7,
  #                  box.padding = 0.5,
  #                  force_pull = -0.2,
  #                  segment.curvature = -0.1,
  #                  segment.ncp = 3,
  #                  segment.angle = 20,
  #                  show.legend = F) +
  theme_bw(base_size = text_size) +
  labs(x = "Log2(foldchange)",
       y = "") +
  scale_fill_manual(name = "Category",
                    labels = c("positive" = "Positive",
                               "negative" = "Negative",
                               "neutral" = "Neutral"),
                    values = RColorBrewer::brewer.pal(11, "RdBu")[round(seq(from = 1, to = 11, length.out = 3))],
  ) +
  theme(legend.spacing.x = unit(1, "line"),
        legend.position = "top",
        legend.justification = "bottom")


ggsave(filename = "fig1_plots.pdf",
       path = here("figures", "figure1_screenoutcome"),
       dpi = 600,
       scale = 2,
       width = 180,
       height = 115,
       units = "mm"
)

#### Figure 2: Taxonomy of interaction categories ####

# Calculate relative counts
ps_rel <- transform_sample_counts(ps, function(x) x/sum(x)*100)

rel_abund_sum <- psmelt(ps_rel) %>%
  mutate(scr_cat = fct_relevel(scr_cat,
                                "control_pos",
                                "neutral",
                                "down",
                                "up")) %>% 
  group_by(scr_cat, Species) %>%
  mutate(median = median(Abundance),
         Species = ifelse(median < 2, "< 2%", str_c("P. ", Species))
         ) %>%
  group_by(Sample, scr_cat, Species) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

# Plot
rel_abund_sum %>%
  ggplot(aes(x = Sample, y = Abundance, fill = Species)) + 
  geom_bar(stat = "identity", 
           col = "grey30") + 
  scale_fill_manual(values = getPalette(n = length(unique(rel_abund_sum$Species)))) +
  coord_cartesian(xlim = c(0.35, 4.7),
                  ylim = c(-1, 101),
                  expand = F) +
  labs(x = "Four replicates per category", 
       y = "Relative abundance (%)") +
  facet_grid(. ~ scr_cat, 
             scales = "free_x",
             labeller = labeller(scr_cat = c("control_pos" = "Total library \nn = 719",
                                             "down" = "Negative \nn = 304",
                                             "neutral" = "Neutral \nn = 395",
                                             "up" = "Promoting \nn = 20"))) +
  theme_classic(base_size = text_size) + 
  theme(strip.background = element_blank(), 
        axis.text.x = element_blank())


ggsave(filename = "fig2_plots.pdf",
       path = here("figures", "figure2_rpodamplicon"),
       dpi = 600, # For pngs
       scale = 2,
       width = 90,
       height = 80,
       units = "mm"
)

#### Figure 3: Reference strains predictions ####
# Type strain pellicle intensities
ts_pel_int <- map(list.files(here("data/typestrains"), pattern = "pellicle", full.names = T), read_csv) %>% 
  bind_rows() %>% 
  mutate(
    Label = str_remove(Label, "type_strains_"),
    Label = str_remove(Label, "_pellicle"),
    Label = str_remove(Label, " - .*"),
    ch = ifelse(str_detect(Label, "1/3"), "mkate", "gfp"),
    type = ifelse(str_detect(Label, "bs[:alpha:]"), "co", "mono")) %>% 
  separate(Label, into = c("medium", "ps_strain"), sep = "_", extra = "drop") %>%
  select(-1) %>% 
  dplyr::rename(area_µm2 = Area,
                px_mean = Mean,
                px_sd = StdDev,
                sum_px_int = RawIntDen) %>% 
  mutate(ps_strain = str_remove(ps_strain, "-.*"),
         ps_strain = str_remove(ps_strain, ":.*"),
         ps_strain = ifelse(type == "mono", "none", str_remove(ps_strain, "bs")),
         ps_strain = fct_relevel(ps_strain, c("none",
                                              "lini",
                                              "poae",
                                              "koreensis",
                                              "capeferrum",
                                              "protegens")),
         ch = fct_relevel(ch, c("mkate",
                                "gfp")),
         medium = fct_relevel(medium, c("lbgm",
                                        "0.1xlbgm")),
         int_per_µm2 = sum_px_int / area_µm2)


# Type strain colony intensities
ts_col_int <- read_csv(here("data/typestrains", "colony_intensities.csv"),
                       col_names = T) %>% 
  mutate(
    Label = str_remove(Label, " - .*"),
    ch = ifelse(str_detect(Label, "1/3"), "mkate", "gfp"),
    type = ifelse(str_detect(Label, "co"), "co", "mono")) %>% 
  separate(Label, into = c("medium", "ps_strain"), sep = "_", extra = "drop") %>%
  select(-1) %>% 
  dplyr::rename(area_µm2 = Area,
                px_mean = Mean,
                px_sd = StdDev,
                sum_px_int = RawIntDen) %>% 
  mutate(ps_strain = ifelse(type == "mono", "none", str_remove(ps_strain, "bs")),
         ps_strain = fct_relevel(ps_strain, c("none",
                                              "lini",
                                              "poae",
                                              "koreensis",
                                              "capeferrum",
                                              "protegens")),
         ch = fct_relevel(ch, c("mkate",
                                "gfp")),
         medium = fct_relevel(medium, c("lbgm",
                                        "0.1xlbgm")),
         int_per_µm2 = sum_px_int / area_µm2)


# ANOVA and Tukey Kramer
pel_letters <- ts_pel_int %>%
  group_by(ch, medium) %>%
  nest() %>%
  mutate(tukey = map(data, ~ TukeyHSD(aov(px_sd ~ ps_strain, .)))) %>%
  dplyr::select(-data) %>%
  mutate(letters = map(tukey, function(x) {
    multcompLetters(x$ps_strain[, "p adj"], reversed = T)$Letters
  })) %>%
  mutate(letters = map(letters, enframe)) %>%
  unnest(letters) %>%
  dplyr::rename(
    ps_strain = name,
    letters = value
  ) %>% 
  dplyr::filter(ps_strain != "koreensis")

col_letters <- ts_col_int %>%
  filter(ch == "mkate") %>% 
  group_by(ch, medium) %>%
  nest() %>%
  mutate(tukey = map(data, ~ TukeyHSD(aov(area_µm2 ~ ps_strain, .)))) %>%
  dplyr::select(-data) %>%
  mutate(letters = map(tukey, function(x) {
    multcompLetters(x$ps_strain[, "p adj"], reversed = T)$Letters
  })) %>%
  mutate(letters = map(letters, enframe)) %>%
  unnest(letters) %>%
  dplyr::rename(
    ps_strain = name,
    letters = value
  ) %>% 
  dplyr::filter(ps_strain != "koreensis")


# Plot pellicle px intensity sd (correlates with wrinkles)
pel_sd_p <- ts_pel_int %>% 
  dplyr::filter(ps_strain != "koreensis") %>% 
  ggplot(aes(ps_strain, px_sd, fill = ch)) +
  stat_summary(aes(col = ch),
               geom = "linerange", 
               fun.max = max, 
               fun.min = min,
               linewidth = 1.5,
               position = position_dodge(width = 0.9)) +
  geom_point(pch = 21,
             size = 3,
             position = position_dodge(width = 0.9)) +
  geom_text(data = pel_letters,
            aes(x = ps_strain, 
                y = ifelse(ch == "mkate", 1700, 800), 
                label = letters, 
                col = ch),
            size = 5,
            position = position_dodge(width = 0.9)) +
  facet_grid(.~medium,
             labeller = labeller(medium = c("0.1xlbgm" = "0.1X LBGM\nmKate exposure time = 5.0s\nGFP exposure time = 5.0s",
                                            "lbgm" = "1x LBGM\nmKate exposure time = 1.5s\nGFP exposure time = 1.0s"))) +
  coord_cartesian(ylim = c(0, 1700)) +
  scale_fill_manual(name = "Channel",
                    labels = c("mkate" = "Bacillus",
                               "gfp" = "Peps"),
                    values = c("mkate" = cust_mag,
                               "gfp" = cust_green)) +
  scale_color_manual(name = "Channel",
                     labels = c("mkate" = "Bacillus",
                                "gfp" = "Peps"),
                     values = c("mkate" = cust_mag,
                                "gfp" = cust_green)) +
  scale_x_discrete(labels = c("none" = "DK1042",
                              "lini" = "+ lini",
                              "poae" = "+ poae",
                              "koreensis" = "+ koreensis",
                              "capeferrum" = "+ capeferrum",
                              "protegens" = "+ protegens"),
                   guide = guide_axis(angle = 30)) +
  labs(x = "",
       y = "Pixel intensity std. deviation") +
  theme_bw(base_size = text_size) +
  theme(
    strip.background = element_blank())


# Plot colony area faceted on media
col_area_p <- ts_col_int %>% 
  filter(ch == "mkate",
         ps_strain != "koreensis") %>% 
  ggplot(aes(ps_strain, area_µm2/10^6, fill = ch)) +
  stat_summary(aes(col = ch),
               geom = "linerange", 
               fun.max = max, 
               fun.min = min,
               linewidth = 1.5,
               position = position_dodge(width = 0.9)) +
  geom_point(pch = 21,
             size = 3,
             position = position_dodge(width = 0.9)) +
  geom_text(data = col_letters,
            aes(x = ps_strain, 
                y = ifelse(ch == "mkate", 95, 40), 
                label = letters),
            size = 5,
            position = position_dodge(width = 0.9)) +
  facet_grid(.~medium,
             labeller = labeller(medium = c("0.1xlbgm" = "0.1X LBGM",
                                            "lbgm" = "1x LBGM"))) +
  coord_cartesian(ylim = c(0, 95)) +
  scale_fill_manual(name = "Channel",
                    labels = c("mkate" = "Bacillus",
                               "gfp" = "Peps"),
                    values = c("mkate" = cust_mag,
                               "gfp" = cust_green)) +
  scale_color_manual(name = "Channel",
                     labels = c("mkate" = "Bacillus",
                                "gfp" = "Peps"),
                     values = c("mkate" = cust_mag,
                                "gfp" = cust_green)) +
  scale_x_discrete(labels = c("none" = "DK1042",
                              "lini" = "+ lini",
                              "poae" = "+ poae",
                              "koreensis" = "+ koreensis",
                              "capeferrum" = "+ capeferrum",
                              "protegens" = "+ protegens"),
                   guide = guide_axis(angle = 30)) +
  labs(x = "",
       y = "Bacillus colony area (mm2)") +
  theme_bw(base_size = text_size) +
  theme(
    strip.background = element_blank())



pel_sd_p + col_area_p +
  plot_layout(guides = "collect")

ggsave(filename = "fig3_plots.pdf",
       path = here("figures", "figure3_typestrains"),
       dpi = 600, # For png
       scale = 2,
       width = 180,
       height = 100,
       units = "mm"
)


#### Figure 4: The sec met potential of sequenced isolates ####
bgc <- readxl::read_xlsx(here("data/antismash", "bgc_presence.xlsx")) %>% 
  dplyr::filter(tax_method == "wgs")


bgc_tree <- bgc %>%
  mutate(label = iso_id) %>% 
  dplyr::filter(tax_method == "wgs", !is.na(pyoverdine_siderophore)) %>% 
  pivot_longer(cols = pyoverdine_siderophore:putisolvin_nrp, names_to = "variable")


bgc_annotation <- bgc %>% 
  mutate(label = iso_id) %>% 
  dplyr::filter(tax_method == "wgs", !is.na(pyoverdine_siderophore)) %>% 
  pivot_longer(cols = pellicle:colony, names_to = "variable")


bgc_matrix <- bgc %>% 
  select(pyoverdine_siderophore:putisolvin_nrp) %>% 
  mutate(across(everything(), ~ifelse(is.na(.x), "no", .x)),
         across(everything(), ~ifelse(.x == "yes", 1, 0)))  

bgc_clust <- bgc_matrix %>% 
  t() %>% 
  dist() %>% 
  hclust()

tree <- ape::read.tree(here("data/antismash", "tygs_tree.phy"))


bgc_summary <- bgc %>% 
  dplyr::filter(tax_method == "wgs", !is.na(pyoverdine_siderophore)) %>% 
  pivot_longer(cols = n_bgc:n_ripp, names_to = "variable", values_to = "n") %>% 
  pivot_longer(cols = pellicle:colony, names_to = "growth_mode", values_to = "growth") %>% 
  dplyr::select(iso_id, growth_mode, growth, variable, n)


bgc_summary_p <- bgc_summary %>% 
  ggplot(aes(variable, n, fill = growth)) +
  stat_summary(aes(col = growth),
               position = position_dodge(width = 0.9),
               geom = "linerange", 
               fun.max = max, 
               fun.min = min,
               linewidth = 1.5) +
  ggforce::geom_sina(pch = 21,
                     size = 3) +
  geom_pwc(label = "p.adj") +
  facet_grid(.~growth_mode,
             labeller = labeller(growth_mode = c("colony" = "Colony",
                                                 "pellicle" = "Pellicle"))) +
  scale_x_discrete(labels = c("n_bgc" = "BGC",
                              "n_nrps" = "NRPS",
                              "n_pks" = "PKS",
                              "n_ripp" = "RiPP")) +
  scale_fill_manual(values = c("no" = custom_red,
                               "yes" = brew_blue),
                    name = "",
                    labels = c("no" = "Inhibition",
                               "yes" = "No inhibition")) +
  scale_color_manual(values = c("no" = custom_red,
                                "yes" = brew_blue),
                     name = "",
                     labels = c("no" = "Inhibition",
                                "yes" = "No inhibition")) +
  labs(x = "BGC category",
       y = "Number of predicted BGCs") +
  theme_bw(base_size = text_size) +
  theme(
    strip.background = element_blank())


bgc_p <- bgc_tree %>% 
  ggplot(aes(x = variable, y = iso_id, fill = value)) +
  geom_tile(
    color = "white",
    lwd = 2,
    linetype = 1
  ) +
#  coord_fixed() +
  scale_fill_manual(values = c("no" = custom_red,
                               "yes" = brew_blue)) +
  theme_minimal(base_size = text_size) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 30))

bgc_ann_p <- bgc_annotation %>% 
  ggplot(aes(x = variable, y = iso_id, fill = value)) +
  geom_tile(
    color = "white",
    lwd = 2,
    linetype = 1
  ) +
  #  coord_fixed() +
  scale_fill_manual(values = c("no" = custom_red,
                               "yes" = brew_blue)) +
  theme_minimal(base_size = text_size) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 30))

wgs_tree <- tree %>% 
  ggtree(layout = "rectangular") + 
  geom_tiplab(align = T) +
  xlim(0, 0.2) +
  theme_void(base_size = text_size)

bgc_tree <- bgc_clust %>% 
  ggtree(layout = "rectangular") +
  geom_tippoint(size = 3) +
  layout_dendrogram()


heat_p <- bgc_p %>% 
  aplot::insert_left(bgc_ann_p,
                     width = 0.2) %>% 
  aplot::insert_left(wgs_tree,
                     width = 0.2) %>% 
  aplot::insert_top(bgc_tree,
                    height = 0.1)


ggsave(filename = "figure4_heatmap.pdf",
       plot = heat_p,
       path = here("figures", "figure4_isolateinteractions"),
       scale = 2,
       width = 180,
       height = 180,
       units = "mm",
       dpi = 600)

ggsave(filename = "figure4_plot.pdf",
       plot = bgc_summary_p,
       path = here("figures", "figure4_isolateinteractions"),
       scale = 2,
       width = 150,
       height = 60,
       units = "mm",
       dpi = 600)


#### Figure 5: Non-inhibitors are not located within the pellicle ####
pell_dist <- map2(list.files(here("data/p9-31_bs_pellicles_confocal"), full.names = T), list.files(here("data/p9-31_bs_pellicles_confocal")), function(file, filename) {
  read.csv(file) %>% 
    as_tibble() %>% 
    mutate(sample = str_remove(filename, ".csv")) %>% 
    separate(sample, into = c("medium", "type", "ch"), sep = "_") %>% 
    dplyr::rename(dist_µm = 'X.micron.',
                  int_mean = Mean)
}) %>% 
  bind_rows()


pell_dist %>% 
  ggplot(aes(x = dist_µm, y = int_mean, col = ch, linetype = type)) +
  geom_line(linewidth = 1) +
  facet_grid(. ~ medium,
             scales = "free_x",
             labeller = labeller(medium = c("0.1xlbgm" = "0.1X LBGM",
                                            "lbgm" = "1X LBGM"))) +
  labs(y = "Mean pixel intensity",
       x = "Image Z position (µm)") +
  coord_flip(clip = "off",
             expand = T) +
  expand_limits(y = 0) +
  scale_color_manual(name = "Channel",
                     values = c("gfp" = cust_green,
                                "mkate" = cust_mag),
                     labels = c("gfp" = "eGFP (P9_31)",
                                "mkate" = "mKate2 (DK1042)")) +
  scale_linetype_manual(name = "Culture type",
                        values = c("co" = "solid",
                                   "mono" = "dotted"),
                        labels = c("co" = "Coculture",
                                   "mono" = "Monoculture")) +
  theme_bw(base_size = text_size) +
  theme(strip.background = element_blank())

ggsave(filename = "fig5_plots.pdf",
       path = here("figures", "figure5_coculturepellicles"),
       dpi = 600, # For png
       scale = 2,
       width = 80,
       height = 85,
       units = "mm"
)


#### Figure 6: DAPG is not the only Bacillus-inhibiting molecule ####
phl <- readxl::read_xlsx(here("data/antismash", "bgc_presence.xlsx")) %>% 
  dplyr::filter(!is.na(phld_pcr))

phl


#### Figure S1: Image analysis workflow ####
# See figure directory

#### Figure S2: Screen results QC ####
# Calculate summary stats
median_lines <- df_median %>% 
  filter(ch == "mKate",
         culture == "dual") %>%
  group_by(medium) %>% 
  summarize(lines = median(se))

# Biovolume per plate in different culture types
p_s1_1 <- df_l2fc %>% 
  filter(ch == "mKate") %>% 
  ggplot(aes(plate, log2(biovolume), col = rep)) +
  geom_point(position = "jitter",
             alpha = 1/2) +
  stat_summary(geom = "crossbar",
               fun = median) +
  scale_color_brewer(palette = "Set2",
                     labels = c("Rep 1",
                                "Rep 2",
                                "Rep 3"),
                     name = "") +
  facet_grid(culture ~ medium,
             scales = "free_y",
             labeller = labeller(culture = c("blk" ="Blank control",
                                             "dual" = "Coculture",
                                             "mono" = "Monoculture"))) +
  theme_bw() +
  labs(x = "Plate",
       y = "Log2(Biovolume) (log2[µm3])") +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = text_size),
    text = element_text(size = text_size)
  ) +
  guides(color = guide_legend(override.aes = list(size = 3,
                                                  alpha = 1)))

# log2(FC) as a function of biovolume
p_s1_2 <- df_l2fc %>% 
  filter(culture == "dual",
         ch == "mKate") %>% 
  ggplot(aes(log2(biovolume), log2fc, col = plate)) +
  stat_summary(geom = "point", 
               fun = median, 
               alpha = 1/2, 
               show.legend = T) +
  scale_color_manual(values = c(RColorBrewer::brewer.pal(6, "BuGn")[4:6],
                                RColorBrewer::brewer.pal(6, "Oranges")[4:6],
                                RColorBrewer::brewer.pal(6, "BuPu")[4:6]),
                     name = "Plate") +
  facet_grid(. ~ medium) +
    theme_bw() +
  labs(x = "Log2(Biovolume) (Log2[µm3])",
       y = "Log2(FC)") +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = text_size),
    text = element_text(size = text_size)
  ) +
  guides(color = guide_legend(override.aes = list(size = 3,
                                                  alpha = 1)))

# Distribution of standard errors (Interactions in TSB have higher standard errors)
p_s1_3 <- df_median %>% 
  filter(ch == "mKate",
         culture == "dual") %>%
  ggplot(aes(se)) +
  geom_density(fill = "grey95") +
  geom_vline(data = median_lines, 
             aes(xintercept = lines),
             linetype = "dashed",
             linewidth = 1) +
  geom_text(data = median_lines,
            aes(x = lines + 0.175,
              label = str_c("Median = ", round(lines, 3))),
            y = 2.3) +
  coord_cartesian(xlim = c(-0.02, 2),
                  ylim = c(-0.05, 2.4),
                  expand = F) +
  facet_grid(medium ~ .) +
  theme_bw() +
  labs(x = "log2(FC) Standard error of the mean",
       y = "Density") +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = text_size),
    text = element_text(size = text_size)
  )

(p_s1_1 + p_s1_2) / p_s1_3


ggsave(filename = "figs2_plots.pdf",
       path = here("figures", "figure_s2_screenraw"),
       dpi = 600,
       scale = 2,
       width = 180,
       height = 115,
       units = "mm"
)

#### Figure S3: Amplicon seq PERMANOVA and PCoA ####
# Calculate distance matrix from DESeq2-normalized data
euc_dist <- dist(t(vst_feature_table))

# PERMANOVA with adonis
anova(betadisper(euc_dist, md$scr_cat)) # 0.4832 - betadispersion is similar across sample; adonis2 result is trustworthy
permanova <- adonis2(euc_dist ~ scr_cat, 
                     data = md, 
                     permutations = 999, 
                     method = "bray")

# Plot PCoA (of normalized data) with phyloseq
vst_pcoa <- ordinate(vst_ps, method = "MDS", distance = "euclidean")

phyloseq::plot_ordination(vst_ps, vst_pcoa) + 
  geom_point(aes(fill = scr_cat),
             pch = 21,
             size = 3) +
  stat_ellipse(aes(fill = scr_cat,
                   col = scr_cat), 
               geom = "polygon", 
               level = 0.95, 
               alpha = 0.1, 
               linewidth = 1, 
               show.legend = F) +
  annotate(geom = "text",
           x = 7.5,
           y = 8,
           label = str_c("PERMANOVA p < ", permanova$`Pr(>F)`[1])) +
  coord_fixed(sqrt(vst_pcoa$values$Eigenvalues[2] / vst_pcoa$values$Eigenvalues[1])) + 
  scale_fill_brewer(name = "Screen category", 
                    labels = c("Total library", 
                               "Inhibiting", 
                               "Neutral", 
                               "Promoting"),
                    palette = "Set2") +
  scale_color_brewer(name = "Screen category",
                     labels = c("Total library",
                                "Inhibiting",
                                "Neutral",
                                "Promoting"),
                     palette = "Set2") +
  labs(x = "PCoA 1 [41.4%]", 
       y = "PCoA 2 [20.7%]") +
  theme_bw(base_size = text_size) 

ggsave(filename = "figs3_plots.pdf",
       path = here("figures", "figure_s3_pcoa"),
       dpi = 600,
       scale = 2,
       width = 88,
       height = 50,
       units = "mm"
)


#### Table 1: Distribution of screen interactions ####
df_median %>%
  dplyr::filter(ch == "mKate", culture == "dual") %>% 
  group_by(threshold, medium) %>% 
  tally() %>% 
  write_csv(here("figures", "table1_categorydist", "table.csv"))


#### Table 2: rpoD differential abundance ####
diff_abund <- ancombc2(data = ps,
                       fix_formula = "scr_cat",
                       p_adj_method = "fdr",
                       group = "scr_cat",
                       pairwise = T,
                       alpha = 0.05,
                       iter_control = list(tol = 0.01, max_iter = 100, verbose = T),
                       verbose = T)

# Pairwise comparison with total library as reference group
diff_abund_all_otu <- out$res %>% 
  dplyr::rename(otu = "taxon") %>% 
  left_join(dplyr::select(feature_table, species, otu), by = "otu") %>% 
  as_tibble() %>%  
  dplyr::select(species, 
                starts_with("lfc"), 
                starts_with("q"),
                starts_with("diff"),
                -matches("Intercept"))

diff_abund_sig <- diff_abund_all_otu %>% 
  filter(if_any(starts_with("diff"), ~ . == "TRUE"))

write_csv(diff_abund_sig,
          file = here("data/ampliconseq", "diff_abund_sig.csv"))

# Pairwise comparisons between categories
diff_abund_all_otu_pairwise <- out$res_pair %>% 
  dplyr::rename(otu = "taxon") %>% 
  left_join(dplyr::select(feature_table, species, otu), by = "otu") %>% 
  as_tibble() %>%  
  dplyr::select(species, 
                starts_with("lfc"), 
                starts_with("q"),
                starts_with("diff"))

diff_abund_sig_pairwise <- diff_abund_all_otu %>% 
  filter(if_any(starts_with("diff"), ~ . == "TRUE"))

write_csv(diff_abund_sig_pairwise,
          file = here("data/ampliconseq", "diff_abund_sig_pairwise.csv"))


#### Dataset S1: Screen categorizations with raw values ####
category_dataset <- df_median %>% 
  ungroup() %>% 
  dplyr::filter(ch == "mKate") %>% 
  dplyr::select(-ch) %>% 
  dplyr::rename("category" = "threshold")

write_csv(category_dataset, 
          file = here("data", "dataset_s1.csv"), 
          col_names = T)
