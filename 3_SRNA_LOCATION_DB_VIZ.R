
# RICARDO GOMEZ-REYES  
# AFTER RUN 2_SRNA_LOCATION_DB_PREP.R

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)


# ACORDING TO ISAAC MARTINEZ UGALDE (CEI ABREU STUDENT)
# GENOMIC ANNOTATION IS CURATED. 
# OVERLAPPING FEATURES ARE HIERARCHICALY CATEGORIZED FOR EACH GENOMIC REGION.
# THEN, ARE USED TO IDENTIFY MIRNA/SRNAS SOURCE
# 1) USING SETDIFF (FROM GENOMICRANGES; Lawrence et al., 2013) TO OBTAIN THE NON-INTERSECTED REGION AND THE INTERSECTED
# 1.2) GENOMIC CATEGORIES ARE USED: SRNA FAMILIES, EXON, INTRONS AND NOVEL REPEAT ELEMENTS.
# 2) THE LEVEL OF SRNA PRODUCTION IN EACH GENOMIC REGION IS OBTAINED USING FIND OVERLAPS (OR subsetover)
# USING THE RESIZE FUNCTION, WITH PARAMETER FIX="CENTER" TO COUNT ACCORDING THE OVERLAP OF THE CENTRAL NUCLEOTIDE OF EACH READ


library(GenomicRanges)

wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out"

head(DB <- read_tsv(paste0(wd, "/RNA_LOCATION_DB.tsv")))

# VIZ NUMBER OF INTRONIC, EXONIC, ETC. 
# FIRST, CURATE THE HIERARCHICAL FEATURE:

out <- DB %>% group_by(SRNAtype, biotype, DicerCall) %>% 
  summarise(Reads = sum(Reads), n = n()) %>%
  group_by(SRNAtype) %>%
  mutate(pct = Reads/sum(Reads))

view(DB %>% head())

out %>%
  # filter(SRNAtype == "miR") %>%
  ggplot(aes(x = DicerCall, y = Reads, fill = type)) +
  # facet_grid(~ SRNAtype, scales = "free", space = "free") +
  geom_col()

# ANALYSIS OF CANONICAL OR ISOFORM (NEILSEN ET AL 2012) ?

DB %>% 
  drop_na(KnownRNAs) %>%
  ggplot(aes(MajorRNAReads/Reads)) +
  geom_histogram()

#

DB %>% 
  # drop_na(KnownRNAs) %>%
  ggplot(aes(FracTop, MajorRNAReads/Reads, color = Strand)) +
  facet_grid(~ SRNAtype) +
  geom_point()


# ANY CUTOFF?

DB %>% count(SRNAtype)

DB %>% 
  filter()
  group_by(SRNAtype) %>%
  count(type, sort = T) %>%
  ggplot()
  
  

