# pulling in samples


# install.packages("tidyverse")
library(tidyverse)
library(scales)

# list all files in a path
file_list <- list.files(path="/Volumes/vandenbrinklab/Susan/Human GVHD/TCRseq files", full.names = T)

unified_all <- file_list %>%
  map_dfr(~ read_tsv(., col_types = cols_only(sample_name = col_character(),
                                          total_templates = col_double(),
                                          productive_templates = col_double(),
                                          total_rearrangements = col_double(),
                                          productive_rearrangements = col_double(),
                                          productive_clonality = col_double(),
                                          productive_entropy = col_double(),
                                          sample_clonality = col_double(),
                                          sample_entropy = col_double(),
                                          sample_amount_ng = col_double(),
                                          max_productive_frequency = col_double(),
                                          rearrangement = col_character(),
                                          amino_acid = col_character(),
                                          frame_type = col_character(),
                                          frequency = col_double(), templates = col_double(),
                                          productive_frequency = col_double(),
                                          cdr3_length = col_double(),
                                          v_family = col_character(),
                                          v_gene = col_character(),
                                          v_allele = col_character(),
                                          d_family = col_character(),
                                          d_gene = col_character(),
                                          d_allele = col_character(),
                                          j_family = col_character(),
                                          j_gene = col_character(),
                                          j_allele = col_character(),
                                          v_resolved = col_character(),
                                          d_resolved = col_character(),
                                          j_resolved = col_character()
  )))

glimpse(unified)

unified_all %>%
  write_csv('unified_all.csv')


