library(here)
library(sangeranalyseR)
library(tidyverse)

# Put your raw .ab1 files in ./seq/raw
# Put your md.csv file in ./seq/raw

path <- here("seq", "tax_rpoD_seq", "scv_15_136_178")

# Tried to make a metadata file for their CSV-method, but failed.
# Still need the md to rename the files
md <- list.files(path = here(path, "raw"), pattern = "*.ab1") %>%  
  as_tibble() %>% 
  # Find the defining barcode digits
  mutate(barcode = str_extract(value, "(?<=[A-Z]{3})[:digit:]*(?=_)")) %>%
  # Join with md.csv (only needs contig and barcode)
  left_join(read_csv2(here(path, "md.csv"),
                      col_types = "ccc"),
            by = "barcode") %>% 
  # File name should be: "contig_direction.ab1"
  mutate(contig = str_c("P8-", contig),   # This is just for me, you could just make the csv-file properly in the first place
         new_file = str_c(contig, direction, sep = "_"),
         new_file = str_c(new_file, "ab1", sep = "."))


if(!dir.exists(here(path, "raw_renamed"))) dir.create(here(path, "raw_renamed"))

#Only do this when you're sure you have the right "new names"
file.copy(from = str_c(here(path, "raw"), "/", pull(md, value)),
            to = str_c(here(path, "raw_renamed"), "/", pull(md, new_file)))

my_aligned_contigs <- SangerAlignment(ABIF_Directory = here(path, "raw_renamed"),
                                      REGEX_SuffixForward = "_F.ab1$",
                                      REGEX_SuffixReverse = "_R.ab1$")

my_aligned_contigs@objectResults@readResultTable

#launchApp(my_aligned_contigs) # This thing throws an error

generateReport(my_aligned_contigs, 
               outputDir = here(path))

# Write fasta file that you can blast
writeFasta(my_aligned_contigs,
           outputDir = path)
