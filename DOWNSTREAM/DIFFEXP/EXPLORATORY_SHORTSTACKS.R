# SHORSTACKS

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

# path <- "~/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230308_test/"
# path <- "~/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230314_test/"
path <- "~/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out/"

list.files(path = path, pattern = "txt")

# count_f <- list.files(path = path, pattern = "Counts.txt", full.names = T)

res_f <- list.files(path = path, pattern = "Results.txt", full.names = T)

# Count <- read_tsv(count_f)

# A tab-delimited text file giving key information for all small RNA clusters. Columns:

Results <- read_tsv(res_f)

# Results %>% arrange(MajorRNA) %>% sample_n(100) %>% view()

# MIRNA column:  ------
# Did the locus pass all criteria to be called a MIRNA locus? If so, 'Y'. If not, 'N'.

Results %>% dplyr::count(MIRNA)

# Q: How many clusters not identify in mirgendb and pirdb pass all criteria to be called a MIRNA locus?

Results %>% filter(is.na(KnownRNAs)) %>% dplyr::count(MIRNA)

Results %>% drop_na(KnownRNAs) %>% dplyr::count(MIRNA)


# 
# Results %>% filter(MIRNA == "N" & is.na(KnownRNAs) == F) %>% view()
# 
# Results %>% filter(MIRNA == "Y" & is.na(KnownRNAs) == T) %>% view()
# 
# Results %>% filter(MIRNA == "Y") %>% drop_na(KnownRNAs) %>% view()

# PIRNAS -----

Results %>% filter(grepl("piR-",KnownRNAs)) %>% count(MIRNA)
# Results %>% filter(grepl("piR-",KnownRNAs) & MIRNA == "Y" ) %>% view()


# Using follow columns to generate fasta headers: (not requiered for shortstacks 4.0)
# Locus: Coordinates of the locus in Chrom:Start-Stop format, one-based, inclusive.
# Name: Name for the Locus. De novo loci are named like Cluster_1, etc. Loci provided by the user with option --locifile preserve the user-given names.

# Length: Length of the locus (base-pairs)

# MajorRNA: Sequence of the single most abundant RNA sequence at the locus.

names(Results)


fasta_prep <- Results %>%
  # select()
  unite("headers", c("Name", "Locus"), sep = "|") %>%
  arrange(headers)

fasta_prep %>% distinct(headers) 

fasta_prep %>% head() %>% view()

srna_seqs <- fasta_prep %>% pull(MajorRNA)
srna_headers <- fasta_prep %>% pull(headers)

srna_fasta <- c(rbind(srna_headers, srna_seqs))

write(srna_fasta, file= paste0(path, "MajorRNA.fasta"))
# Results %>% head() %>% view()

Results %>% 
  distinct(MajorRNA, .keep_all = T) %>%
  mutate(color = "Novel") %>%
  mutate(color = ifelse(grepl("piR",KnownRNAs), "piR", color)) %>%
  mutate(color = ifelse(grepl("Mir",KnownRNAs), "miR", color)) %>% 
  # count(color)
  ggplot(aes(Length, Reads, color = color)) + geom_point() + facet_grid(MIRNA ~ color)

# Polarity distribution of transcript-mapped sRNAs ----
# 21: Number of 21 nucleotide reads aligned to the locus.
# i: Number of i to dicermax nucleotide reads aligned to the locus 


# Results %>% 
#   select(Name, Strand, any_of(as.character(21:30))) %>%
#   pivot_longer(cols = any_of(as.character(21:30)), names_to = "Length", values_to = "n") %>%
#   mutate(Length = ifelse(Length < 22, "<22", ifelse(Length > 24, ">24", Length))) %>%
#   group_by(Length, Strand) %>% 
#   summarise(n = sum(n)) %>%
#   ungroup() %>%
#   mutate(frac = n/sum(n)) %>% 
#   ggplot(aes(x = Length, y = frac, fill = Strand)) +
#   geom_col()

Results %>%
  mutate(Length = ifelse(Length < 20, "<20", ifelse(Length > 24, ">24", Length))) %>%
  count(Length, Strand) %>% mutate(frac = n/sum(n)) %>% 
  ggplot(aes(x = Length, y = frac, fill = Strand)) +
  geom_col(position = "dodge2")


# Radar charts show the fractions of sRNAs in each of the 21 possible phasing registers; the registers highlighted in magenta are those predicted by the miRNA target sites.
  
  
# DicerCall ----

# If >= 80% of all aligned reads are within the boundaries of --dicermin and --dicermax, than the DicerCall gives the size of most abundant small RNA size. If < 80% of the aligned reads are in the --dicermin and --dicermax boundaries, DicerCall is set to 'N'. Loci with a DicerCall of 'N' are unlikely to be small RNAs related to the Dicer-Like/Argonaute system of gene regulation.
Results %>% count(DicerCall)

