# Parsing TCR data downloaded from Adaptive Biotech

Download the data as follows:
- navigate to <https://clients.adaptivebiotech.com/pub/dewolf-2022-gvhd>
- click explore this project
- log in (making an account if needed)
- click "Open in Analyses"
- click "Analyses on the nav pane
- select all
- click export -> export Sample (NOT  V2!)
- wait patiently
- continue waiting; resist the urge to try the process again
- eventually move downloaded zip to this repo, unpack, and move the tsv files to `data/raw_TCR_data/`.



The following code was used to generate the merged file. The resulting file `TCR/tcr_unified.csv` is used for later analyses including the GLIPH2 analysis. 

```{r 00-preprocess-TCR-2 }
library(tidyverse)
knitr::opts_chunk$set(echo = TRUE)

tcr_data_dir <- "data/raw_TCR_data/"

dir.create("data/TCR", showWarnings = FALSE)

tcr_raw_output <- file.path("data/TCR", "tcr_unified.csv")

if (!file.exists(tcr_raw_output)) {
  raw_files <- sort(dir(tcr_data_dir, pattern = "*.tsv", full.names = TRUE))
  if (length(raw_files) == 0) {
    stop(paste0("No files found within ", tcr_data_dir))
  }
  tcr_raw_original <- map(raw_files, .f = function(x) {
    # x = "sampleExport.2022-06-01_05-18-41//Balb-c_1_blood.tsv"
    read.csv(x, sep = "\t") %>% mutate(sample_name = gsub("\\.tsv", "", basename(x)))
  }) %>%
    bind_rows() %>%
    mutate(
      # sample_name = gsub("T_cells", "Tcells", sample_name),
      org = case_when(
        grepl("^Pt", sample_name) ~ "Human",
        grepl("^Balb-c", sample_name) ~ "Mouse",
        grepl("^Recipient_", sample_name) ~ "Mouse",
        grepl("^Donor", sample_name) ~ "Mouse",
        TRUE ~ "FixMe"
      ),
      samp_type = case_when(
        grepl("^Donor", sample_name) ~ "Donor",
        grepl("^Pt", sample_name) ~ "Patient",
        grepl("^Balb-c", sample_name) ~ "Healthy",
        grepl("^Recipient_", sample_name) ~ "Recipient",
        TRUE ~ "FixMe"
      ),
      indv_type = case_when(
        grepl("[Dd]onor", sample_name) ~ "Donor",
        grepl("^Pt", sample_name) ~ "Recipient",
        grepl("^Balb-c", sample_name) ~ "Recipient",
        grepl("^Recipient_", sample_name) ~ "Recipient",
        TRUE ~ "FixMe"
      ),
      gvhd = case_when(
        grepl("[Dd]onor", sample_name) ~ "noGvHD",
        grepl("^Pt[A-G]", sample_name) ~ "GvHD",
        grepl("^Recipient_", sample_name) ~ "GvHD",
        grepl("^Balb", sample_name) ~ "noGvHD",
        grepl("^Pt[H-J]", sample_name) ~ "noGvHD",
      ),
      indv_id = case_when(
        grepl("^Balb-c", sample_name) ~ gsub(".*?_(\\d*?)_.*", "\\1", sample_name),
        grepl("^Pt", sample_name) ~ gsub("Pt(.+?)_.*", "\\1", sample_name),
        grepl("^Donor", sample_name) ~ gsub(".*?_(\\d*?)_.*", "\\1", sample_name),
        grepl("Recipient_tx", sample_name) ~ gsub("Recipient_tx(\\d*?)_(\\d+)_.*", "\\2", sample_name),
        TRUE ~ "FixMe"
      ),
      tissue = case_when(
        indv_type == "Donor" & org == "Mouse" ~ "Tcells",
        indv_type == "Donor" & org == "Human" ~ "graft_blood",
        grepl("^Balb-c", sample_name) ~ gsub(".*_\\d_(.*)$", "\\1", sample_name),
        grepl("^Pt", sample_name) ~ gsub("Pt._(.*)$", "\\1", gsub("_TCRB", "", sample_name)),
        grepl("Recipient_tx", sample_name) ~ gsub(".*_D\\d+_(.*)$", "\\1", sample_name),
        TRUE ~ "FixMe"
      ),
      exp_rep = ifelse(grepl("Recipient_tx", sample_name), gsub("Recipient_tx(\\d*?)_(\\d+)_.*", "\\1", sample_name), NA),
      collection_day = ifelse(grepl("Recipient_tx", sample_name), gsub("Recipient_tx.*D(\\d*?)_.*", "\\1", sample_name), NA)
    ) %>%
    mutate(
      indv_label = paste(samp_type, indv_id),
      indv_label = ifelse(!is.na(exp_rep), paste0(indv_label, ", rep", exp_rep), indv_label),
      indv_label = ifelse(!is.na(collection_day), paste0(indv_label, ", day ", collection_day), indv_label),
      indv_tissue_label = paste(samp_type, indv_id, "-", gsub("_", " ", tissue)),
      sample_amount_ng = round(sample_amount_ng, 2)
    )


  write.csv(tcr_raw_original %>% arrange(sample_name, rearrangement), file = tcr_raw_output, row.names = FALSE)
}


```


 Note the rounding done to some fields, which was done to overcome some float precision issues in R. For instance, when done without the rounding, results like the following can be problematic:

```
read.csv("data/TCR/tcr_unified.csv") %>% filter(sample_name == "Balb-c_1_spleen") %>% pull(sample_amount_ng) %>% unique()
     [1] 1195.458 1195.458
```


