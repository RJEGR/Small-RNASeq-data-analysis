# TEST OVERLAPS IRANGE

# GENOMIC SOURCE ===
# /Users/cigom/Documents/MIRNA_HALIOTIS/RNA_SEQ_ANALYSIS/REPEAT_MASKER_OUT
# multi_genome.newid.fa.out.gff

# /Users/cigom/Documents/MIRNA_HALIOTIS/ENSEMBLE
# Haliotis_rufescens_gca023055435v1rs.xgHalRufe1.0.p.56.gtf
# multi_genome.newid.gtf


# TRANSCRIPTOMIC SOURCE ===

# /Users/cigom/Documents/MIRNA_HALIOTIS/RNA_SEQ_ANALYSIS/ASSEMBLY/STRINGTIE/REFBASED_MODE/ANNOTATION
# transcripts.fa.transdecoder.gff3

# /Users/cigom/Documents/MIRNA_HALIOTIS/RNA_SEQ_ANALYSIS/ASSEMBLY/STRINGTIE/REFBASED_MODE
# transcripts.gtf

# SRNA SOURCE ===

# path <- '/Users/cigom/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out'
# Results.gff3	knownRNAs.gff3

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(GenomicRanges)
library(tidyverse)


# RANGE 1 ====
srna_f <- list.files(path = path, pattern = "Results.gff3", full.names = T)

srna_df <- read_tsv(srna_f, col_names = F)

srna_elements <- srna_df %>% distinct(X3) %>% pull()

start_ <- srna_df$X4

end_ <- srna_df$X5

ranges_ <- IRanges(start = start_, end = end_)

seqnames_ <- srna_df$X1

score_ <- srna_df$X6

strand_ <- srna_df$X7

strand_ <- gsub("[.]", "*", strand_)

Rle(seqnames_)

gr <- GRanges(Rle(seqnames_), ranges = ranges_, strand = strand_)

# a DataFrame object containing the metadata columns. Columns cannot be named "seqnames", "ranges", "strand", "seqlevels", "seqlengths", "isCircular", "start", "end", "width", or "element".

# GRanges(seqnames = "chrZ", IRanges(start=c(5,10),end=c(35,45)),
  # strand="+", seqlengths=c(chrZ=100L))

miRNAs_loc <- c("MIRNA_hairpin", "mature_miRNA", "miRNA-star")

table(srna_df$X3)

srna_df %>% group_by(X3, X7) %>% 
  # tally(X6, sort = T) %>%
  group_by(X3, X7) %>% 
  summarise(n = n(), TotalReads = sum(X6)) %>%
  mutate(col = X3) %>% 
  mutate(col = ifelse(col %in% miRNAs_loc, "miRNA_locus", col)) %>%
  mutate(col = ifelse(grepl("siRNA", col), "siRNA_locus", col)) %>%
  ungroup() %>%
  mutate(pct = n / sum(n)) %>%
  ggplot() +
  facet_wrap(~ col, scales = "free") +
  geom_col(aes(y = n, x = X3, fill = X7), 
    position = position_dodge(width = 0.5), width = 0.4) +
  theme_bw(base_family = "GillSans", base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, 
    hjust = 1, vjust = 1, size = 14))


# THIS STEP IS TO FILTER ACCORDING TO SRNA UPPRESSED BY CONDITION

srna_df %>% select(X9) %>% 
  separate(col =X9, sep = ";", into = c("Name", "DicerCall", "MIRNA")) %>%
  count(MIRNA)


# RANGE 2 =====

path <- "/Users/cigom/Documents/MIRNA_HALIOTIS/RNA_SEQ_ANALYSIS/REPEAT_MASKER_OUT"

mask_f <- list.files(path = path, pattern = "multi_genome.newid.fa.out.gff", full.names = T)

mask_df <- read_tsv(mask_f, col_names = F, comment = "#")

table(mask_df$X4)

start_ <- mask_df$X4

end_ <- mask_df$X5

score_ <- mask_df$X6

ranges_ <- IRanges(start = start_, end = end_)

seqnames_ <- mask_df$X1

strand_ <- mask_df$X7

strand_ <- gsub("[.]", "*", strand_)

gr2 <- GRanges(Rle(seqnames_), ranges = ranges_, strand = strand_)

# SEARCH OVERLAPS

# gr %over% gr2
# gr1[gr1 %over% gr2]

fo <- findOverlaps(gr, gr2, minoverlap = 1)

# Hits object with 2656 hits and 0 metadata columns:

queryHits(fo) # position vector where coordinates from query is found

subjectHits(fo) # position vector where coordinates from subject is found

# ranges from the query for which we found a hit in the subject

index = queryHits(fo)

ovlps <- gr[index,]

srna_df[index,] %>% 
  group_by(X3, X7) %>% 
  summarise(n = n(), TotalReads = sum(X6)) %>%
  mutate(col = X3) %>% 
  mutate(col = ifelse(col %in% miRNAs_loc, "miRNA_locus", col)) %>%
  mutate(col = ifelse(grepl("siRNA", col), "siRNA_locus", col)) %>%
  ungroup() %>%
  mutate(pct = n / sum(n)) %>%
  ggplot() +
  geom_col(aes(y = n, x = col, fill = X7), 
    position = position_dodge(width = 0.5), width = 0.4) +
  theme_bw(base_family = "GillSans", base_size = 12)

# 

coverage(gr)
# maxPos <- which.max(ctcfCoverage10)
# 
# > roi <- resize(IRanges(maxPos, width = 1), 5000, “center”)
# 
# > roiCoverage <- ctcfCoverage$chr10[roi]

# RANGE 3 ====

path <- "/Users/cigom/Documents/MIRNA_HALIOTIS/RNA_SEQ_ANALYSIS/ASSEMBLY/STRINGTIE/REFBASED_MODE/"


orf_f <- list.files(path = paste0(path, "ANNOTATION"), 
  pattern = "transcripts.fa.transdecoder.gff3", full.names = T)

orf_df <- read_tsv(df_f, col_names = F, comment = "#")

gtf_f <- list.files(path = path, pattern = "^transcripts.gtf", full.names = T)

gtf_df <- read_tsv(gtf_f, col_names = F, comment = "#")

srna_df %>% select(X9) %>% 
  separate(col =X9, sep = ";", into = c("Name", "DicerCall", "MIRNA")) %>%
  count(MIRNA)


orf_df %>% distinct(X1) %>% filter(!grepl("LOC", X1))

orf_df %>% head() %>% view()
gtf_df %>% filter(grepl("LOC", X9)) %>% head() %>% view()

start_ <- df$X4

end_ <- df$X5

score_ <- mask_df$X6

ranges_ <- IRanges(start = start_, end = end_)

seqnames_ <- mask_df$X1

strand_ <- mask_df$X7

strand_ <- gsub("[.]", "*", strand_)

gr3 <- GRanges(Rle(seqnames_), ranges = ranges_, strand = strand_)


# WHICH MITOCONDRIAL
# JALGQA010000616.1 <- MITOCHONDR

srna_df %>% filter(X1 %in% "JALGQA010000616.1") %>%
  group_by(X3) %>%
  summarise(n = n(), TotalReads = sum(X6))

mask_df %>%
  filter(X1 %in% "JALGQA010000616.1") %>%
  group_by(X3) %>%
  summarise(n = n(), TotalReads = sum(X6))
