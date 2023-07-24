
# PREPARE INPUTS:

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

path <- "~/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out/"

count_f <- list.files(path = path, pattern = "Counts.txt", full.names = T)

MTD_f <- list.files(path = path, pattern = "METADATA.tsv", full.names = T)

recode_to <- c(`Low` = "Experimental")

library(tidyverse)

MTD <- read_tsv(MTD_f) %>%
  dplyr::mutate(pH = dplyr::recode_factor(pH, !!!recode_to)) %>%
  dplyr::mutate(pH = factor(pH, levels = c("Control", "Experimental"))) %>%
  mutate(LIBRARY_ID = gsub("^HR", "", LIBRARY_ID)) 

COUNTS <- read_tsv(count_f)

colNames <- names(read_tsv(count_f))

names(COUNTS)[names(COUNTS) %in% colNames] <- gsub("HR|.clean.newid.subset", "", colNames)

rowNames <- COUNTS$Name

COUNTS <- COUNTS %>% select_if(is.double) %>% as(., "matrix")

rownames(COUNTS) <- rowNames

head(COUNTS)

out <- list(COUNTS, MTD)

write_rds(out, file = paste0(path, "Counts_and_mtd.rds"))
