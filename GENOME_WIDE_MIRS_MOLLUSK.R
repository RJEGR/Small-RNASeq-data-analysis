
# Ricardo Gomez-Reyes
# Huang et al (2021) had (in silico) scan potential mollusks miRs using a genome-wide approach: 

# 1. Predicted mature miRs were allowed to have only 0 - 4 mismatch in sequences with know mature miRNA
# 2. The mismatched nucleotides were not permitted in know miRNA seed region (2 - 8 bp)
# 3. miRNA precursor can fold into an appropriate hairpin secondary structure that contains mature miRNA sequence within one arm of the hairpin and has the smallest possible folding energy.

# Here we generate precursor and mature fasta format for downstream analysis:

# Table S1. Genome-wide scanning miRNAs in mollusks
# Table S2. Target prediction of miRNAs in C. gigas and P. fucata.

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

path <- "~/Documents/MIRNA_HALIOTIS/GENOME_WIDE_DISCOVERY_AND_FUNCTIONAL_MIRS_IN_MOLLUSCS/"

xls_f <- list.files(path = path, pattern = "^Table S1", full.names = T)

library(tidyverse)
library(readxl)

# To load all sheets in a workbook, use lapply()


# The list of sheet names is especially useful when you want to iterate them:

# excel_sheets(xls_f)

my_custom_name_repair <- function(nms) tolower(gsub(" ", "_", nms))

workbook <- lapply(excel_sheets(xls_f), 
  read_excel, 
  path = xls_f, 
  range = cell_cols("B:D"), 
  .name_repair = my_custom_name_repair)

# name workbook list by sheet names

names(workbook) <- excel_sheets(xls_f)


workbook <-
  mapply(`[<-`, workbook, 'sp', value = names(workbook), SIMPLIFY = FALSE)

workbook <- do.call(rbind, workbook)

workbook <- workbook %>% distinct()

# hist(nchar(workbook$mature_sequence))
# hist(nchar(workbook$mirna_precursor))

workbook <- workbook %>% filter(nchar(mature_sequence) >= 18)

# =================
# Save fasta file
# ================
# help in https://astrobiomike.github.io/amplicon/dada2_workflow_ex

miRs_header <- workbook %>% unite("headers", c("mirna", "sp"), sep = "|") %>% pull(headers)
head(miRs_header <- paste(">", miRs_header, sep=""))

mature_seqs <- workbook %>% pull(mature_sequence)
hairpin_seqs <- workbook %>% pull(mirna_precursor)

mature_fasta <- c(rbind(miRs_header, mature_seqs))
hairpin_fasta <- c(rbind(miRs_header, hairpin_seqs))

write(mature_fasta, file = paste0(path, "/mature.fasta"))
write(hairpin_fasta, file = paste0(path, "/precursor.fasta"))


# exit