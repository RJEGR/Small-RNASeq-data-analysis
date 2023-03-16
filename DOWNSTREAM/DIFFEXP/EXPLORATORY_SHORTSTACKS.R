# SHORSTACKS

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

# path <- "~/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230308_test/"
path <- "~/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230314_test/"

list.files(path = path, pattern = "txt")

count_f <- list.files(path = path, pattern = "Counts.txt", full.names = T)

res_f <- list.files(path = path, pattern = "Results.txt", full.names = T)

# Count <- read_tsv(count_f)

# A tab-delimited text file giving key information for all small RNA clusters. Columns:

Results <- read_tsv(res_f)

Results %>% arrange(MajorRNA) %>% head() %>% view()

# MIRNA column: Did the locus pass all criteria to be called a MIRNA locus? If so, 'Y'. If not, 'N'.

Results %>% count(MIRNA)

Results %>% filter(MIRNA == "Y" & is.na(KnownRNAs) == T) %>% view()

Results %>% filter(MIRNA == "Y") %>% drop_na(KnownRNAs) %>% view()


# Using follow columns to generate fasta headers:
# Locus: Coordinates of the locus in Chrom:Start-Stop format, one-based, inclusive.
# Name: Name for the Locus. De novo loci are named like Cluster_1, etc. Loci provided by the user with option --locifile preserve the user-given names.

# Length: Length of the locus (base-pairs)

# MajorRNA: Sequence of the single most abundant RNA sequence at the locus.

names(Results)


fasta_prep <- Results %>%
  select()
  unite("headers", c("Name", "#Locus"), sep = "|") 

fasta_prep %>% arrange(headers) %>% distinct(headers) 


# Results %>% head() %>% view()
