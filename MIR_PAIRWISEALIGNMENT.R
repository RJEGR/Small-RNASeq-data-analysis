
# RICARDO GOMEZ-REYES
# # SEARCH HOMOLOBY W/ MIRBASE
# USING MOLLUSK DB ALSO AS INPUT


rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

path <- "~/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out/"

print(.DB <- read_tsv(paste0(path, "/RNA_LOCATION_DB.tsv")))

path <- "~/Documents/MIRNA_HALIOTIS/GENOME_WIDE_DISCOVERY_AND_FUNCTIONAL_MIRS_IN_MOLLUSCS/"

# STEP 0 ====

library(Biostrings)

dnaref <- readDNAStringSet(paste0(path, "/molluscs_mature.fa"))

dnarefdb <- read_rds(paste0(path, "/molluscs_mature.rds")) %>% dplyr::rename("pattern_name" = "header")

# INSTEAD, USE MIRBASE

mirbasef <- "miRBase-mature.fa"

destfile <- paste0(path, "/", mirbasef)

# download.file("https://mirbase.org/download/mature.fa", destfile = destfile)

dnaref <- DNAStringSet(readRNAStringSet(destfile))

# dnaref[grepl("let-7", names(dnaref))]

# query <- as(data.frame(row.names = query$Name, query$MajorRNA), "matrix")

query <- .DB %>% filter(SRNAtype == "miR" & is.na(KnownRNAs)) %>% 
  select(Name, MajorRNA) %>% arrange(MajorRNA) %>%
  pull(MajorRNA, name = Name)

query <- DNAStringSet(RNAStringSet(query))

# STEP 1 =====

MpairwiseAlignment <- function(f,q, type = c("match", "align")) { 
  
  # From: James F. Reid" <james.reid@ifom-ieo-campus.it
  # https://rdrr.io/bioc/microRNA/src/inst/scripts/sequence-similarity.R
  
  type <- match.arg(type)
 
  sub.mat <-
    switch(type,
      match = {
        ans <- matrix(-10000L, 4, 4,
          dimnames = list(c("A", "T", "C", "G"),
            c("A", "T", "C", "G")))
        diag(ans) <- 1L
        ans
      },
      align = {
        ans <- matrix(-1L, 4, 4,
          dimnames = list(c("A", "T", "C", "G"),
            c("A", "T", "C", "G")))
        diag(ans) <- 1L
        ans
      })
  
  gap <- 
    switch(type,
      match = -10000L,
      align = -1L)
  
  MA <-
    pairwiseAlignment(pattern = f, subject = q,
      type = "local",
      substitutionMatrix = sub.mat,
      gapOpening = gap,
      gapExtension = -10000L,
      scoreOnly = FALSE)
  
  # similar to a blastn for search mirs en genomes: (10.1101/gr.193367.115)
  # BLASTN with the following settings: − word_size 4 -reward 5 -penalty − 4 -gapopen 8 -gapextend 6. 

  # a gap opening penalty of 5 and a gap extension penalty of 2. Moreover, to generate a substitution matrix for sequence alignment, we set the match score to 1 and the mismatch score to -1. (doi.org/10.1371/journal.pcbi.1006931)
  
  
  # MA <- pairwiseAlignment(pattern = f, subject = q, type = "local", gapOpening = 5, gapExtension = 2,
  #   substitutionMatrix =
  #     nucleotideSubstitutionMatrix(match = 1, mismatch = -1, type = "DNA"))
  
  
  
  dna_as_ch <- function(x) { as.character(x)}
  
  # PID1: 100 * (identical positions) / (aligned positions + internal gap positions)
  
  # Normalize score as doi.org/10.1371/journal.pcbi.1006931
  #  The sequence similarity score obtained for each miRNA pair was further normalized to the range [0, 1] by the following equation:
  
  s <- c(sample(seq(0,100),10))
  
  (s - min(s)) / (max(s) - min(s))
  
  normalized_score <- function(x) { (x - min(x)) / (max(x) - min(x)) }
  
  OUT <- data.frame(pattern = dna_as_ch(pattern(MA)),  
    subject = dna_as_ch(subject(MA)),
    score = as.numeric(dna_as_ch(score(MA))),
    normalized_score = normalized_score(as.numeric(dna_as_ch(score(MA)))),
    Identity = pid(MA, type = "PID1"),
    # PID2 =  pid(MA, type = "PID2"),
    pattern_name = names(f),
    subject_name = names(q))
  
  
  OUT %>% arrange(desc(score))
  
}

# MpairwiseAlignment(f,q)

# pattern may be multiple sequences"

DF <- MpairwiseAlignment(f = dnaref, q = query[1], "match") %>% as_tibble()

DF %>% view()

# DF %>% left_join(dnarefdb) %>% view()

hist(DF$normalized_score)

DF <- list()

for (i in 1:length(query)) {

  q <- query[i]
    
  DF[[i]] <- MpairwiseAlignment(f = dnaref, q = q, "match") %>% as_tibble()
  
}

do.call(rbind, DF) -> DF

write_rds(DF, file = paste0(path, "/MIRS_pairwiseAlignment.rds"))

# STEP 2 =====

DF <- read_rds(paste0(path, "/MIRS_pairwiseAlignment.rds"))

# HOW TO CHOOSE THE BEST HIT?

DF <- DF %>% filter(normalized_score == 1) 

DF <- DF %>% group_by(subject_name) %>% filter(score == max(score))

DF %>% group_by(subject_name) %>% slice_head(n =1) # ?

str(mature_id <- sapply(strsplit(DF$pattern_name, " "), '[', 1))

str(g_name <- sapply(strsplit(DF$pattern_name, " "), '[', 3))

str(s_name <- sapply(strsplit(DF$pattern_name, " "), '[', 4))

DF$KnownRNAs <- mature_id

DF$Family <- substr(mature_id, 5,nchar(mature_id)-1)

DF$sp <- paste0(g_name, " ", s_name) 

# view(DF)

#   
# plot(DF$score, DF$normalized_score)

# RESOLVE SINGLE MATCHING

SINGLE_NAMES <- DF %>% count(subject_name) %>% filter(n == 1) %>% pull(subject_name)

DF1 <- DF %>% 
  filter(subject_name %in% SINGLE_NAMES) %>%
  select(subject_name,KnownRNAs, sp) 

view(DF1)

# RESOLVE MULTIPLE MATCHING

str(MULTIPLE_NAMES <- DF %>% count(subject_name) %>% filter(n != 1) %>% pull(subject_name))

DF %>% count(subject_name, sort = T)

DF %>% 
  filter(subject_name %in% MULTIPLE_NAMES) %>%
  group_by(subject_name) %>%
  summarise(n = n(), across(Family, .fns = paste_go), .groups = "drop_last") %>%
  view()

data.frame(pattern_seq = as.character(dnaref)) %>% 
  as_tibble(rownames = "pattern_name") %>%
  right_join(DF) %>% 
  filter(subject_name %in% SINGLE_NAMES) %>% 
  view()
