
# Ricardo Gomez-Reyes
# Huang et al (2021) obtained (in silico) scan potential mollusks miRs using a genome-wide approach: 

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

MTD <- list.files(path = path, pattern = "METADATA.xls", full.names = T)

MTD <- readxl::read_excel(MTD)


xls_f <- list.files(path = path, pattern = "^Table S1", full.names = T)

library(tidyverse)
library(readxl)

# To load all sheets in a workbook, use lapply()


# The list of sheet names is especially useful when you want to iterate them:

# excel_sheets(xls_f)

my_custom_name_repair <- function(nms) tolower(gsub(" ", "_", nms))

.workbook <- lapply(excel_sheets(xls_f), 
  read_excel, 
  path = xls_f, 
  range = cell_cols("B:D"), 
  .name_repair = my_custom_name_repair)

# name workbook list by sheet names

names(.workbook) <- excel_sheets(xls_f)

workbook <-
  mapply(`[<-`, .workbook, 'sp', value = names(.workbook), SIMPLIFY = FALSE)

workbook <- do.call(rbind, workbook)

workbook <- workbook %>% distinct()

# hist(nchar(workbook$mature_sequence))
# hist(nchar(workbook$mirna_precursor))

workbook <- workbook %>% filter(nchar(mature_sequence) >= 18)

workbook <- workbook %>% left_join(MTD, by = "sp")

llist <- function(x) {
  x <- paste(x, sep = '|', collapse = '|')
  x <- list(x)
  x <- unlist(x)
}

# DEBIDO A SECUENCIAS DUPLICADAS, AGRUPAREMOS IDS POR SECUENCIAS
workbook %>% distinct(mature_sequence)

workbook %>% arrange(mirna)

# workbook <- workbook %>% 
#   arrange(mirna) %>%
#   distinct(mature_sequence, sp,.keep_all = T)

# summarise(across(type, .fns = paste_headers), .groups = "drop_last"
  
paste_sp <- function(x) { 
  x <- x[!is.na(x)] 
  x <- unique(x)
  # n <- length(x)
  
  x <- paste(x, sep = '|', collapse = '|') 
  # x <- paste(n, x, sep = '|', collapse = '|') 
  }

count_sp <- function(x) { 
  x <- x[!is.na(x)] 
  x <- unique(x)
  n <- length(x)}

# workbook %>% filter(mature_sequence == "TGGACGGAGAACTGATAAGGG") %>% view()

# %>%

.workbook <- workbook %>% 
  group_by(mature_sequence) %>%
  summarise(
    n = count_sp(sp),
    across(mirna, .fns = paste_sp),
    across(sp, .fns = paste_sp), .groups = "drop_last") %>% 
  arrange(desc(mature_sequence)) %>%
  # mutate(header = ifelse(n > 5, paste0("Highly_conserved|", n), sp)) %>%
  mutate(header = paste("Mollusk_miR",1:nrow(.), sep = "_")) 

# .workbook %>% right_join(workbook, by = "mature_sequence") 

# 
# workbook <- workbook %>% 
#   # distinct(mature_sequence, sp,.keep_all = T) %>%
#   group_by(mature_sequence) %>%
#   mutate(n_sp = length(sp), 
#     across(sp, .fns = llist)) %>% 
#   ungroup() %>%
#   distinct(mature_sequence, .keep_all = T) %>% 
#   # summarise(n_sp = length(sp), 
#     # across(sp, .fns = llist)) %>% 
#   mutate(n_sp = paste0(n_sp, "_sp")) %>%
#   unite("headers", c("mirna", "n_sp"), sep = "_") 
#   

# workbook <- workbook %>%
  # mutate(n = mature_sequence) %>% unite("headers", c("headers", "n"), sep = "_") 


.workbook %>% arrange(header) %>% distinct(header) 


write_rds(.workbook, file = paste0(path, "/molluscs_mature.rds"))

#  workbook %>% unite("headers", c("mirna", "sp"), sep = "|")

# =================
# Save fasta file
# ================
# help in https://astrobiomike.github.io/amplicon/dada2_workflow_ex

miRs_header <- .workbook %>% pull(header)

head(miRs_header <- paste(">", miRs_header, sep=""))

mature_seqs <- .workbook %>% pull(mature_sequence)

mature_fasta <- c(rbind(miRs_header, mature_seqs))

write(mature_fasta, file = paste0(path, "/molluscs_mature.fa"))

# hairpin_seqs <- .workbook %>% pull(mirna_precursor)
# hairpin_fasta <- c(rbind(miRs_header, hairpin_seqs))
# write(hairpin_fasta, file = paste0(path, "/molluscs_precursor.fa"))


# exit