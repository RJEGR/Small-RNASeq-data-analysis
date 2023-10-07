# SEED SEARCH of 147 mirs, AND UPGRADE THE MAIN DB,  
# FILTER TO ONE SEED PER FAMILIES (EX. MIR-92-3p W/ AUUGCAC, UUGCACU SEED TO ONLY ONE) 
rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

path <- "~/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out/"

head(.DB <- read_tsv(paste0(path, "/RNA_LOCATION_DB.tsv")))

MIRGENEDB <- read_tsv(paste0(path, "SRNA2MIRGENEDB.tsv")) %>%
  tidyr::unite("Family", c("Family", "arm"), sep = "-")


SEED <- MIRGENEDB %>% distinct(Seed) %>% pull() %>% sort()

str(seqs.name <- SEED)

SEED <- DNAStringSet(RNAStringSet(SEED))

names(SEED) <- seqs.name


# 2) 

.DB %>% dplyr::count(SRNAtype)

DB <- .DB %>% filter(SRNAtype == "miR")

str(SEQS <- DB %>% arrange(MajorRNA) %>% pull(MajorRNA)) # %>% head(100)

str(seqs.name <- DB %>% arrange(MajorRNA) %>% select(Name) %>% pull(Name))

SEQS <- DNAStringSet(RNAStringSet(seqs))

names(SEQS) <- seqs.name

# MpairwiseAlignment(f = SEED, q = SEQS[1], type = "match") %>% as_tibble()

DF <- list()

for (i in 1:length(SEQS)) {
  
  q <- SEQS[i]
  
  DF[[i]] <- MpairwiseAlignment(f = SEED, q = q, "match") %>% as_tibble()
  
}

do.call(rbind, DF) -> DF

DF %>%
  # group_by(subject_name) %>% 
  # sample_n(1) %>%
  ggplot(aes(score, normalized_score, color = Identity)) + geom_point()


# USING SEED LENGTH = 7, 
DF <- DF %>% filter(score == 7) %>% rename("subject_name" = "Name")

DF %>% distinct(subject_name) # ONLY 87 MIRS ()

DF %>% count(subject_name, sort = T)

DF %>% arrange(pattern_name) %>% distinct(pattern_name)

DNASEQS <- as.data.frame(SEQS) %>% as_tibble(rownames = "Name")

DF %>% 
  left_join(DNASEQS, by = "Name") %>% 
  left_join(distinct(MIRGENEDB, Family, Seed), by = c("pattern_name"="Seed")) %>%
  view()
