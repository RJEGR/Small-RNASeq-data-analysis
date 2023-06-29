
# STRAND CUTOFF ANALYSIS

DB %>% 
  # count(Strand) %>%
  ggplot(aes(FracTop, fill = Strand)) +
  geom_histogram() + facet_wrap(~ Strand, scales = "free_y")

# default: 0.8. Loci with >80% reads on the top genomic strand are '+' stranded, loci with <20% reads on the top genomic strand are '-' stranded, and all others are unstranded '.'

wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out"

head(DB <- read_tsv(paste0(wd, "/RNA_LOCATION_DB.tsv")))


# ANALYSIS OF CANONICAL OR ISOFORM (NEILSEN ET AL 2012)
DB %>% 
  drop_na(KnownRNAs) %>%
  ggplot(aes(MajorRNAReads/Reads)) +
  geom_histogram()
