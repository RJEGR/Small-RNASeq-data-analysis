
# STRAND CUTOFF ANALYSIS

DB %>% 
  # count(Strand) %>%
  ggplot(aes(FracTop, fill = Strand)) +
  geom_histogram() + facet_wrap(~ Strand, scales = "free_y")

# default: 0.8. Loci with >80% reads on the top genomic strand are '+' stranded, loci with <20% reads on the top genomic strand are '-' stranded, and all others are unstranded '.'
