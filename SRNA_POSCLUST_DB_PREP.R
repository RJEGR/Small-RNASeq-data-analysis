
# RICARDO GOMEZ-REYES
# PART OF STEP: SRNA POS-CLUSTERING AND NOVEL CRITERIA TO IDENTIFY FUNCTIONAL SIRNAS/PIRNAS:

# PREPARE INFO REGARGIND:
# 1) SAMPLE FREQUENCY (I.E. PREVALENCE, SPREAD OR CO-OCURRANCE OF SRNAS IN THE SAMPLES)
# 2) SRNA ABUNDANCES (ROWSUM OF SRNAS. IT WOULD OR WOULDNT BE CORRELATED TO SAMPLE FREQ.)
# 3) SRNA-SEQUENCE GROUPING (I.E. BASED ON SEQ. SIMILARITY GROUP IF SEQS > 70 % Sequence Identity)

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

path <- "~/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out/"

out <- read_rds(paste0(path, "Counts_and_mtd.rds"))

dim(COUNTS <- out[[1]])

dim(MTD <- out[[2]])

# head(sort_m <- sort(rowSums(COUNTS), decreasing = T))

# any(names(sort_m) %in% rownames(COUNTS))

# head(COUNTS[match(rownames(COUNTS), names(sort_m) ),])

# ORDER SRNA MATRIX BY DECREASING COUNT (Freq. and total read count):

head(SORTED_MAT <- COUNTS[order(rowSums(COUNTS),decreasing=T),])

READS_AND_FREQ <- COUNTS %>%
  as_tibble(rownames = "Name") %>%
  pivot_longer(-Name, values_to = "Reads") %>%
  mutate(REP = substr(name, 1,nchar(name)-1)) %>%
  filter(Reads > 0) %>%
  group_by(Name, REP) %>%
  summarise(Reads = sum(Reads), Freq = n())

# READS_AND_FREQ %>% group_by(REP) %>% cor_test(Reads, Freq)
READS_AND_FREQ %>% select(-Freq) %>% pivot_wider(names_from = REP, values_from = Reads, values_fill = 0)
READS_AND_FREQ %>% select(-Reads) %>% pivot_wider(names_from = REP, values_from = Freq, values_fill = 0)


head(DB <- read_tsv(paste0(path, "/RNA_LOCATION_DB.tsv")))

DB <- READS_AND_FREQ %>% summarise(Reads = sum(Reads), Freq = sum(Freq)) %>%
  left_join(DB) %>% rename("SampleFreq" = "Freq") %>%
  arrange(match(Name, rownames(SORTED_MAT)))

DB  %>% ggplot(aes(x = Reads, y = SampleFreq, color = SRNAtype)) + 
  stat_ecdf()

s1 <- DNAString("AGTATAGATGATAGAT")
s2 <- DNAString("AGTAGATAGATGGATGATAGATA")

palign1 <- pairwiseAlignment(s1, s2)
palign1
pid(palign1)

# palign2 <-
#   pairwiseAlignment(s1, s2,
#     substitutionMatrix =
#       nucleotideSubstitutionMatrix(match = 2, mismatch = 10, baseOnly = TRUE))
# palign2
# pid(palign2, type = "PID4")

# TRY SEQUENCE SIMILARITY ANALYSIS AND 

seqs <- DB %>% distinct(MajorRNA) %>% pull() %>% head(100)

bioseq::rna(seqs)

# print(seqs <- Biostrings::RNAStringSet(seqs))

# Biostrings::pairwiseAlignment(seqs, seqs)

# alignedDNA <- DECIPHER::AlignSeqs(RNAStringSet(seqs))

bioseq::seq_consensus(bioseq::rna(seqs))



library(msa)

# method=c("ClustalW", "ClustalOmega", "Muscle")

align <- msa::msa(RNAStringSet(seqs), method = "ClustalW")

align <- msaConvert(align)$seq

print(bioseq::rna(align))

bioseq::seq_cluster(bioseq::rna(align))
