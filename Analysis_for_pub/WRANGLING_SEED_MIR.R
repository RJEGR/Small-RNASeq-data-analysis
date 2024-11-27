
# What is the number of seed families targeting same gene?
# Find the microRNA-seed family per target, with a microRNA degree > 1
# Then calculate proper density of target site by:
# Density = S / D, where S: Deegre of microRNA, and D: N seed constituted in the microRNA-seed family per target

# LOAD RELATIONAL-DB
# Use/create kmer method to find seed 



library(tidyverse)

dir <- "~/Documents/MIRNA_HALIOTIS/MIRNA_PUB_2024/"

print(DB <- read_tsv(file.path(dir, "LONGER_RELATIONAL_DB.tsv")))


SEEDDB <- DB %>% 
  mutate(Seed = substr(MajorRNA, 2,8)) %>%
  group_by(biotype_best_rank) %>%
  distinct(gene_id, Seed) %>%
  group_by(gene_id) %>%
  summarise(across(Seed, .fns = paste_col), Seed_degree = n()) %>% 
  arrange(desc(Seed_degree))


DENSITYDB <- DB %>% 
  group_by(biotype_best_rank) %>%
  distinct(STRINGID, gene_id, MajorRNA) %>% 
  group_by(STRINGID, gene_id) %>%
  summarise(across(MajorRNA, .fns = paste_col), miR_degree = n()) %>% 
  arrange(desc(miR_degree)) %>%
  left_join(SEEDDB) %>% ungroup()


DENSITYDB <- DENSITYDB %>%
  mutate(STRINGID = ifelse(is.na(STRINGID), gene_id, STRINGID)) %>%
  mutate(Diversity = 1 -(Seed_degree/miR_degree)) %>%
  arrange(desc(Seed_degree)) %>%
  mutate(STRINGID = factor(STRINGID, levels = unique(STRINGID))) 

lo = 0# floor(min(DENSITYDB$y))
up = 1 #ceiling(max(DENSITYDB$y))
mid = (lo + up)/2

DENSITYDB %>%
  ggplot(aes(x = STRINGID, y = Seed_degree, fill = Diversity)) + 
  # facet_grid(biotype_best_rank ~.) +
  geom_col() +
  theme_classic(base_family = "GillSans", base_size = 12) +
  labs(y = "Density (seed_degree)", x = "Gene target") +
  theme(legend.position = 'top', 
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 7)) +
  scale_fill_gradient2(
    low = "white", high = "gray5", mid = "gray50",
    # low = "blue", high = "red", mid = "white", 
    na.value = "white", midpoint = mid, limit = c(lo, up),
    name = "Diversity") 

Average_degree <- mean(DENSITYDB$Seed_degree)

# Test using logo motif -----

paste_col <- function(x) { 
  x <- x[!is.na(x)] 
  x <- unique(sort(x))
  x <- paste(x, sep = ';', collapse = ';')
}

TARGETDB <- DB %>%
  mutate(MajorRNA = substr(MajorRNA, 2,13)) %>%
  distinct(STRINGID, gene_id, MajorRNA) %>% 
  group_by(STRINGID, gene_id) %>%
  summarise(across(MajorRNA, .fns = paste_col), n = n()) %>% 
  arrange(desc(n)) %>%
  mutate(STRINGID = ifelse(is.na(STRINGID), gene_id, STRINGID)) %>%
  mutate(STRINGID = paste0(STRINGID, " | ", n)) %>%
  filter(n > Average_degree)

QUERY <- split(strsplit(TARGETDB$MajorRNA, ";") , TARGETDB$STRINGID)
str(QUERY <- lapply(QUERY, unlist))

library(ggseqlogo)

ggseqlogo(QUERY[5], nrow = 5, method = "probability")

LOGO <- geom_logo(QUERY[5], p = F, method = "probability")

lo = floor(min(LOGO$y))
up = ceiling(max(LOGO$y))
mid = (lo + up)/2

LOGO %>% group_by(position, order) %>% summarise(sum(y))

LOGO %>% 
  # mutate()
  ggplot(aes(y = y, x = position, fill = letter)) + 
  geom_col() +
  facet_grid(seq_group ~.)
  scale_fill_gradient2(
    low = "white", high = "gray5", mid = "gray50",
    # low = "blue", high = "red", mid = "white", 
    na.value = "white", midpoint = mid, limit = c(lo, up),
    name = NULL) 

ggseqlogo(QUERY[1:10], nrow = 10) #+
  # annotate('segment', x = 2, xend=8, y=0.6, yend=0.6, size=2) + 
  # annotate('text', x=6, y=0.5, label='Text annotation') 


# test seed analysis using multiple-sequence-alignment:  ----
# using ClustalW or muscle, (Time-computer consuming)
# As Sempere, Wheeler , et al. 
# using Clustal /MAFFT 6.85 /MUSCLE. AND 
# manually refined the alignments with RALEE/Gblocks/

library(msa)

seqs <- DB %>% distinct(MajorRNA) %>% pull()

#  "ClustalOmega"

msa_align <- msa::msa(Biostrings::RNAStringSet(seqs), 
  method = "ClustalOmega")


# msa_align <- msa::msa(Biostrings::RNAStringSet(seqs), 
#   method = "ClustalW", gapOpening = 12, 
#   gapExtension = 1)

RNAStringSet(msa_align)

.align <- msaConvert(msa_align)$seq

.align <- gsub("U","T", msaConvert(msa_align)$seq)

library(bioseq)

bioseq::seq_cluster(bioseq::dna(.align))

table(SEQ_CLUSTERS <- bioseq::seq_cluster(bioseq::dna(.align), threshold = 0.85, method = "single"))


# Testing using kd model (complicate)
library(scanMiR)
seed <- "AGCAUUAA"  # Seed sequence for hsa-miR-155-5p
data("SampleTranscript")  # Load sample transcript
matches <- findSeedMatches(SampleTranscript, seed, verbose = FALSE)

# The KdModel class contains the information concerning the sequence (12-mer) affinity of a given miRNA, and is meant to compress and make easily manipulable the dissociation constants (Kd) predictions from McGeary, Lin et al. (2019).

SEEDKD <- DB %>%
  mutate(Seed = substr(MajorRNA, 2,8)) %>%
  mutate(Seed = paste0(Seed, "AAAA")) %>%
  distinct(MajorRNA, Seed) %>%
  group_by(Seed) %>%
  summarise(across(MajorRNA, .fns = paste_col), n = n()) %>%
  arrange(desc(n))

SEEDKD <- split(strsplit(SEEDKD$MajorRNA, ";") , SEEDKD$Seed)
str(SEEDKD <- lapply(SEEDKD, unlist))


data(SampleKdModel)

assignKd <- function(x) {
  
  SeedName <- names(unlist(x))
  
  # cat("Running",SeedName)
  
  x <- gsub("U","T", x)
  
  assignKdType(x, SampleKdModel) # %>% as_tibble() %>% mutate(Seed = SeedName)
  
  
}



DF <- lapply(SEEDKD, assignKd)

do.call(rbind, DF) %>% as_tibble(rownames = "Seed") %>% arrange(desc(log_kd))

assignKdType(c("CTAGCATTAAGT","ACGTACGTACGT"), SampleKdModel)

# Creating kd model

kd <- dummyKdData()
head(kd)

mod3 <- getKdModel(kd=kd, mirseq="TTAATGCTAATCGTGATAGGGGTT", name = "my-miRNA")

getKdModel2 <- function(seq) {
  
  seq <- gsub("U","T", seq)
  
  cat(seq)
  
  mod <- getKdModel(kd=kd, mirseq=seq, name = seq, conservation = 0)
  
}

lapply(SEEDKD[10], getKdModel2)

