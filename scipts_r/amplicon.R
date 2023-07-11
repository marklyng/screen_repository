### Ampliconseq - Screen metataxonomy
pacman::p_load(here, tidyverse, ggfortify, cluster, vegan, phyloseq)

library("decontam") # Bioconductor
library("ANCOMBC") # Bioconductor
library(DESeq2)


getPalette = colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))

#### Read-in and wrangle ####
bcs <- str_c("bc", str_pad(seq_along(1:31), 2, side = "left", pad = "0"))

md <- read_csv2(here("metadata", "amplicon_md.csv"))
md_no_ctrl <- md %>% 
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
  mutate(otu = str_c("otu", seq(1:nrow(.))))

feature_table_no_ctrl <- feature_table %>%
  dplyr::select(species,
                otu,
                bc13:bc28) # Remove GFP categories and negative controls

df_wide <- feature_table_no_ctrl %>% 
  dplyr::select(-otu) %>% 
  gather(key = key, value = value, 2:ncol(.)) %>% 
  spread(key = names(.[1]), value = "value") %>% 
  dplyr::rename(bc = key) %>% 
  left_join(md_no_ctrl, by = "bc")

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


# Phyloseq object without GFP and negative controls
tax_table_no_ctrl <- feature_table_no_ctrl %>%
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

otu_table_no_ctrl <- feature_table_no_ctrl %>% 
  select(-species) %>% 
  column_to_rownames("otu") %>% 
  otu_table(taxa_are_rows = T)

sample_data_no_ctrl <- md_no_ctrl %>%
  column_to_rownames("bc") %>% 
  sample_data()


ps_no_ctrl <- phyloseq(otu_table_no_ctrl, tax_table_no_ctrl, sample_data_no_ctrl)


#### Decontam ####
# https://benjjneb.github.io/decontam/vignettes/decontam_intro.html
df <- as.data.frame(sample_data(ps)) # Put sample_data into a ggplot-friendly data.frame

df$library_size <- sample_sums(ps)

df <- df[order(df$library_size), ]

df$index <- seq(nrow(df))

ggplot(df, 
       aes(index, library_size, color = scr_cat)) + 
  geom_point()




# Prevalence method (with negative controls)
sample_data(ps)$is_neg <- sample_data(ps)$scr_cat == "control_neg"

contamdf_prev <- isContaminant(ps, method = "prevalence", neg = "is_neg")

table(contamdf_prev$contaminant) # No contaminants identified with this method. Test with frequency method instead.



#### Normalize for sample depth with DESeq2 ####
deseq_counts <- DESeqDataSetFromMatrix(column_to_rownames(select(feature_table_no_ctrl, -species), "otu"), 
                                       colData = md_no_ctrl,
                                       design = ~scr_cat) 
# we have to include the "colData" and "design" arguments because they are 
# required, as they are needed for further downstream processing by DESeq2, 
# but for our purposes of simply transforming the data right now, they don't 
# matter
deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)

# and here is pulling out our transformed table
vst_feature_table <- assay(deseq_counts_vst)

# making our phyloseq object with transformed table
vst_ft_ps <- otu_table(vst_feature_table, taxa_are_rows = T)

vst_ps <- phyloseq(vst_ft_ps, tax_table_no_ctrl, sample_data_no_ctrl)


# PERMANOVA with adonis
anova(betadisper(euc_dist, md_no_ctrl$scr_cat)) # 0.4832 - betadispersion is similar across sample; adonis2 result is trustworthy

permanova <- adonis2(log_trans ~ scr_cat, 
                     data = md_no_ctrl, 
                     permutations = 999, 
                     method = "bray")


# generating and visualizing the PCoA with phyloseq
vst_pcoa <- ordinate(vst_ps, method = "MDS", distance = "euclidean")
eigen_vals <- vst_pcoa$values$Eigenvalues # allows us to scale the axes according to their magnitude of separating apart the samples

phyloseq::plot_ordination(vst_ps, vst_pcoa) + 
  geom_point(aes(fill = scr_cat),
             pch = 21,
             size = 3) +
  stat_ellipse(aes(fill = scr_cat,
                   col = scr_cat), 
               geom = "polygon", 
               level = 0.95, 
               alpha = 0.1, 
               size = 1, 
               show.legend = F) +
  annotate(geom = "text",
           x = 7.5,
           y = 8,
           label = "PERMANOVA p < 0.001") +
  coord_fixed(sqrt(eigen_vals[2] / eigen_vals[1])) + 
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
  theme_bw() + 
  theme(text = element_text(size = 12.5))


#### Normalize to 100k reads and take log() ####
ps_norm <- psmelt(ps) %>% 
  filter(!Sample %in% c("bc29", "bc30", "bc31")) %>% 
  group_by(Sample) %>% 
  mutate(log_abund_norm = log(Abundance/(sum(Abundance)/100000))) %>% 
  mutate(log_abund_norm = ifelse(log_abund_norm == "-Inf", 0, log_abund_norm)) %>% 
  ungroup()

# two-tailed t-test pairwise between total library (x) and each reporter/scr_cat-combination (y)
ps_pval <- ps_norm %>% 
  filter(scr_cat == "control_pos") %>% 
  select(OTU, rep, log_abund_norm) %>% 
  rename(log_norm_ctrl = log_abund_norm) %>% 
  full_join(filter(ps_norm, scr_cat != "control_pos"), by = c("OTU", "rep")) %>% 
  group_by(OTU, reporter, scr_cat, Species) %>% 
  summarize(tval = t.test(x = log_norm_ctrl, y = log_abund_norm)$statistic,
         pval = t.test(x = log_norm_ctrl, y = log_abund_norm)$p.value) %>% 
  ungroup() %>% 
  mutate(padj = p.adjust(pval, method = "BH"))
  

ps_pval %>% 
  filter(padj < 0.05) %>% 
  group_by(reporter, scr_cat) %>% 
  group_split()


#### Plot relative abundance ####
# Calculate relative counts
ps_rel <- transform_sample_counts(ps_no_ctrl, function(x) x/sum(x)*100)

rel_abund_sum <- psmelt(ps_rel) %>%
  group_by(scr_cat, Species) %>%
  mutate(median = median(Abundance),
         Species = ifelse(median < 2, "< 2%", str_c("P. ", Species))) %>%
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
  theme_classic(base_size = 12.5) + 
  theme(strip.background = element_blank(), 
        axis.text.x=element_blank())


ggsave(filename = "ampl_relabund_mkateonly.pdf",
       path = here("figs"),
       device = "pdf",scale = 1.5, dpi = 600)



#### Differential abundance analysis ####

diff_abund <- ancombc2(data = ps_no_ctrl,
                       fix_formula = "scr_cat",
                       p_adj_method = "fdr",
                       group = "scr_cat",
                       global = T,
                       pairwise = T,
                       alpha = 0.05,
                       iter_control = list(tol = 0.01, max_iter = 100, verbose = T),
                       verbose = T)


out$res_global %>% 
  dplyr::rename(otu = "taxon") %>% 
  left_join(dplyr::select(feature_table, species, otu), by = "otu") %>% 
  as_tibble() %>% 
  dplyr::filter(diff_abn == T)

out$res_pair %>% 
  dplyr::rename(otu = "taxon") %>% 
  left_join(dplyr::select(feature_table, species, otu), by = "otu") %>% 
  as_tibble() %>% 
  dplyr::filter(diff_scr_catup_scr_catdown == T) %>% 
  dplyr::select(species)

out$res %>% 
  dplyr::rename(otu = "taxon") %>% 
  left_join(dplyr::select(feature_table, species, otu), by = "otu") %>% 
  as_tibble() %>% 
  dplyr::filter(diff_scr_catneutral == T) %>% 
  dplyr::select(ends_with("cat_up"), species)

# select the bottom 20 with lowest p values
res.or_p <- rownames(res$q_val[,"scr_catup"])[base::order(res$q_val[,"scr_catup"])]
taxa_sig <- res.or_p[1:20]
ps_rel_sig <- prune_taxa(taxa_sig, ps_rel)

