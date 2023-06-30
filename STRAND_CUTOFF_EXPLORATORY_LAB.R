
# STRAND CUTOFF ANALYSIS


rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

DB %>% 
  # count(Strand) %>%
  ggplot(aes(FracTop, fill = Strand)) +
  geom_histogram() + facet_wrap(~ Strand, scales = "free_y")

# default: 0.8. Loci with >80% reads on the top genomic strand are '+' stranded, loci with <20% reads on the top genomic strand are '-' stranded, and all others are unstranded '.'

library(tidyverse)

wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out"

head(DB <- read_tsv(paste0(wd, "/RNA_LOCATION_DB.tsv")))


DB %>% arrange(desc(FracTop))
# head() %>% View()

# ANALYSIS OF CANONICAL OR ISOFORM (NEILSEN ET AL 2012)
DB %>% 
  drop_na(KnownRNAs) %>%
  ggplot(aes(MajorRNAReads/Reads)) +
  geom_histogram()


SRNAS <- read_rds(paste0(wd, "/KNOWN_CLUSTERS_MIRS_PIRS.rds"))

str(which_pirs <- SRNAS %>% filter(grepl("piR", KnownRNAs)) %>% distinct(Name) %>% pull())
str(which_mirs <- SRNAS %>% filter(!grepl("piR", KnownRNAs)) %>% distinct(Name) %>% pull())

# REPLACE 

DB %>% 
  # mutate(RNAs = NA) %>%
  # mutate(RNAs = ifelse(grepl("Mir", KnownRNAs), "miRs", ifelse(grepl("piR", KnownRNAs), "piR", RNAs))) %>% 
  mutate(KnownRNAs = ifelse(is.na(KnownRNAs), Locus_type, KnownRNAs)) %>%
  group_by(KnownRNAs) %>%
  count(type, biotype, sort = T) %>% view()
  
  
