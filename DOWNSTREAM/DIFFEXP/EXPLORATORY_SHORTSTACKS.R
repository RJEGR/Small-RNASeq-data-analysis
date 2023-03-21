# SHORSTACKS

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

# path <- "~/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230308_test/"
# path <- "~/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230314_test/"
path <- "~/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out//"

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
Results %>% drop_na(KnownRNAs) %>% dplyr::count(MIRNA)

Results %>% filter(MIRNA == "N" & is.na(KnownRNAs) == T) %>% view()

Results %>% filter(MIRNA == "Y") %>% drop_na(KnownRNAs) %>% view()

# PIRNAS -----

Results %>% filter(grepl("piR-",KnownRNAs)) %>% count(MIRNA)
Results %>% filter(grepl("piR-",KnownRNAs) & MIRNA == "Y" ) %>% view()


# Using follow columns to generate fasta headers: (not requiered for shortstacks 4.0)
# Locus: Coordinates of the locus in Chrom:Start-Stop format, one-based, inclusive.
# Name: Name for the Locus. De novo loci are named like Cluster_1, etc. Loci provided by the user with option --locifile preserve the user-given names.

# Length: Length of the locus (base-pairs)

# MajorRNA: Sequence of the single most abundant RNA sequence at the locus.

names(Results)


fasta_prep <- Results %>%
  # select()
  unite("headers", c("Name", "#Locus"), sep = "|") 

fasta_prep %>% arrange(headers) %>% distinct(headers) 


# Results %>% head() %>% view()



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
  geom_col()


# Radar charts show the fractions of sRNAs in each of the 21 possible phasing registers; the registers highlighted in magenta are those predicted by the miRNA target sites.
  
  
# DicerCall ----

# If >= 80% of all aligned reads are within the boundaries of --dicermin and --dicermax, than the DicerCall gives the size of most abundant small RNA size. If < 80% of the aligned reads are in the --dicermin and --dicermax boundaries, DicerCall is set to 'N'. Loci with a DicerCall of 'N' are unlikely to be small RNAs related to the Dicer-Like/Argonaute system of gene regulation.
Results %>% count(DicerCall)

# alignment details ----
# mapping_type
# U: Uniquely mapped (not a multimapper).
# P: Multimapper placed using the method set by option --mmap.
# R: Multimapper placed at random.
# H: Very highly multi-mapped read (>=50 hits).
# N: Unmapped reads.


recode_to <- c(`U` = "Uniquely ", `P`= "Multimapper (mmap)",`R` = "Random mmap", `H` = "Highly mmap (>=50 hits) ", `N` = "Unmapped")


path <- "~/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out//"

f <- list.files(path = path, pattern = "alignment_details.tsv", full.names = T)
alignment_details <- read_tsv(f)

readfile <- gsub(".clean.newid.subset.bam", "", basename(alignment_details$readfile))

alignment_details$readfile <- readfile


read_lengthL <- c("<21","21","22","23","24",">24")

# mapping_type read_length   count
ylab <- "Frac"

alignment_details %>%
  group_by(mapping_type, read_length) %>% 
  summarise(n = sum(count)) %>%
  group_by(read_length) %>%
  mutate(frac = n/sum(n)) %>%
  dplyr::mutate(mapping_type = dplyr::recode_factor(mapping_type, !!!recode_to)) %>%
  dplyr::mutate(read_length = factor(read_length, levels = read_lengthL)) %>%
  ggplot(aes(x = read_length, y = frac, fill = mapping_type)) +
  geom_col(width = 0.85) + # position="dodge"
  scale_y_continuous(ylab, labels = scales::percent) +
  guides(fill = guide_legend(title = "")) +
  see::scale_fill_pizza(reverse = T) +
  theme_bw(base_family = "GillSans", base_size = 11) +
  theme(strip.background = element_rect(fill = 'white', color = 'white'),
    legend.position = 'top',
    panel.border = element_blank(),
    panel.grid.minor = element_blank()) -> ps

# https://easystats.github.io/see/articles/seecolorscales.html#overview-of-palette-colors

ggsave(ps, filename = 'ALIGNMENT_DETAILS.png', path = path, width = 6.7, height = 4, device = png)
