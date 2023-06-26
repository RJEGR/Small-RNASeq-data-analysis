
rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

path <- "~/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out/"

list.files(path = path, pattern = "txt")

res_f <- list.files(path = path, pattern = "Results.txt", full.names = T)

RESULTS <- read_tsv(res_f)

# Q: Number of novel and know miRs annotated

RES <- RESULTS %>% 
  select(Name, Chrom, KnownRNAs, MIRNA, Reads, UniqueReads, MajorRNA) %>%
  mutate(KnownRNAs = ifelse(is.na(KnownRNAs) & MIRNA == "Y", Name, KnownRNAs)) %>%
  # mutate(MIRNA = ifelse(!is.na(KnownRNAs), "Y", MIRNA)) %>%
  drop_na(KnownRNAs) %>%
  mutate(Type = ifelse(grepl("Cluster_", KnownRNAs), "Novel", "known"))

# count number of true novel + homology-identify sRNA locus

RES %>% count(Type)

SEQS <- RESULTS %>% select(Name, MajorRNA)

# Q: HOW MANY KNOW PI/MIRNAS WERE IDENTIFY ====

bind_srnas <- function(x) {
  x <- paste(x, sep = '-', collapse = '-')
}

SPLIT <- RESULTS %>%
  select(Name, KnownRNAs, MIRNA) %>%
  drop_na(KnownRNAs) %>% 
  mutate(KnownRNAs = strsplit(KnownRNAs, ";")) %>%
  unnest(KnownRNAs) %>% 
  mutate(biotype = ifelse(grepl("^piR-", KnownRNAs), "piRNA", "Mir"))

SPLIT %>% 
  filter(MIRNA == "Y") %>%   # filter(MIRNA == "N") %>% # SWICHT TO "Y" 
  distinct(Name, biotype) %>% 
  group_by(Name) %>%
  summarise(across(biotype, .fns = bind_srnas), .groups = "drop_last") %>%
  count(biotype)

# # SPLIT NOVEL/KNOWN MIRNAS ====

# WHICH CLUSTERS ARE MIRS EITHER: 
# HOMOLOGY-BASED:SHORTSCACKS PREDICTED AND HOMOLOGY-BASED TYPES ====

paste_headers <- function(x) {
  x <- paste(x, sep = '|', collapse = '|')
  # x <- list(x)
  # x <- unlist(x)
}

MIRS <- SPLIT %>% 
  # filter(MIRNA == "N") %>% # <- OMIT TO SAVE BOTH TYPES
  filter(biotype == "Mir") %>% 
  group_by(Name) %>%
  summarise(across(KnownRNAs, .fns = paste_headers), n = n(), .groups = "drop_last") 

# INCLUDE THE PIR-HOMOLOGY-BASED:SHORTSCACKS PREDICTED ==== 

MIRS <- SPLIT %>% 
  filter(MIRNA == "Y" & biotype == "piRNA") %>%
  group_by(Name) %>%
  summarise(across(KnownRNAs, .fns = paste_headers), n = n(), .groups = "drop_last") %>%
  filter(!Name %in% MIRS$Name) %>%
  rbind(MIRS) 

.MIRS <- MIRS %>% left_join(SEQS, by = "Name")

MIRS <- MIRS %>% 
  mutate(KnownRNAs = ifelse(n > 5, paste0("Highly_conserved|", n), KnownRNAs)) %>%
  select(-n)

# PROCESS FASTA FORMAT ====

MIRS <- MIRS %>%
  left_join(SEQS, by = "Name") %>%
  mutate(Name = paste0(">", Name)) %>%
  unite("Name", Name:KnownRNAs, sep = " ") 


# INCLUDE NOVEL SHORTSCACKS PREDICTED

NOVEL_MIRS <- RESULTS %>% 
  filter(is.na(KnownRNAs)) %>% 
  filter(MIRNA == "Y") %>% 
  mutate(Name = paste0(">", Name)) %>%
  mutate(Name = paste(Name, Start, End, Strand, sep = ":")) %>%
  select(Name, MajorRNA) 

# 

fasta_prep <- NOVEL_MIRS %>% rbind(MIRS)

NOVEL_MIRS %>% rbind(MIRS) %>% distinct(MajorRNA)

RESULTS %>% 
  # filter(is.na(KnownRNAs)) %>% 
  filter(MIRNA == "Y") %>% distinct(MajorRNA)

# MajorRNA: Sequence of the single most abundant RNA sequence at the locus

# AUNQUE HAY SECUENCIAS UNICAS, RECORDAR QUE EL MajorRNA ES LA SECUENCIA REPRESENTATIVA, SUGIRIENDO ISOFORMAS

# WRITE MIRNA FASTA ==== 

seqs <- fasta_prep %>% pull(MajorRNA)

headers <- fasta_prep %>% pull(Name) 

fasta <- c(rbind(headers, seqs))

write(fasta, file= paste0(path, "KNOWN_AND_NOVEL_MIRS_MajorRNA.fasta"))

# PIRNAS TO FURTHER ANALYSIS W/ PIRSCAN ======

SPLIT %>% 
  filter(MIRNA == "N") %>%  
  distinct(Name, biotype) %>% 
  group_by(Name) %>%
  summarise(across(biotype, .fns = bind_srnas), .groups = "drop_last") %>%
  count(biotype)

PIRNAS <- SPLIT %>% 
  filter(MIRNA == "N") %>% # <- OMIT TO SAVE BOTH TYPES
  filter(biotype == "piRNA") %>% 
  filter(!Name %in% .MIRS$Name) # EXCLUDE MIR:PIRS

# Because only 164 unique sequence from 831 CLUSTERS:

PIRNAS %>% left_join(SEQS, by = "Name") %>% distinct(MajorRNA)

.PIRNAS <- PIRNAS %>% 
  group_by(Name) %>%
  summarise(
    across(KnownRNAs, .fns = paste_headers),n = n(), .groups = "drop_last") %>%
  left_join(SEQS, by = "Name")

# USE CLUSTER NAME INSTEAD OF PIRNA NAME:

PIRNAS <- PIRNAS %>% 
  left_join(SEQS, by = "Name") %>% 
  group_by(MajorRNA) %>%
  distinct(Name) %>%
  summarise(
    # across(KnownRNAs, .fns = paste_headers),
    across(Name, .fns = paste_headers), n = n(), .groups = "drop_last")


# PROCESS FASTA FORMAT ====

PIRNAS <- PIRNAS %>% 
  mutate(Name = ifelse(n > 5, paste0("MultiLocus|", n), Name))

PIRNAS <- PIRNAS %>%  mutate(Name = paste0(">", Name)) 



# WRITE PIRNA FASTA ==== 

seqs <- PIRNAS %>% pull(MajorRNA)

headers <- PIRNAS %>% pull(Name) 

fasta <- c(rbind(headers, seqs))

write(fasta, file= paste0(path, "KNOWN_PIRS_MajorRNA.fasta"))


write_rds(rbind(.PIRNAS, .MIRS), file = paste0(path, "KNOWN_CLUSTERS_MIRS_PIRS.rds"))

# Q: Which other loci were for siRNA locus ====

.PIRNAS 

.MIRS

pattern <- "Results.gff3"

f <- list.files(path = path, pattern = pattern, full.names = T)

assembly <- rtracklayer::import(f)

# Assembly track is different in nrows because Parent clusters

df <- assembly %>% as_tibble() %>% 
  select(ID, width, type) %>%
  right_join(RES, by = c("ID"="Name"))

df %>% 
  group_by(type, MIRNA, Type) %>% 
  summarise(n = n(), Reads = sum(Reads)) %>% 
  arrange(desc(n)) %>% view()

# The majority of homology-based piRNAs loci are Unknown Unknown_sRNA_locus, followed by 27 to 29 siRNA_locus and 

df <- assembly %>% as_tibble() %>% 
  # select(ID, width, type) %>%
  anti_join(RES, by = c("ID"="Name"))

df %>% select(Parent) %>% arrange(desc(Parent)) %>% head() %>% pull()

# score: number of sRNA-seq aligned reads at that locus.

df %>% group_by(type, MIRNA) %>% 
  summarise(n = n(), Reads = sum(score)) %>% 
  arrange(desc(n)) %>% view()

# For a detailed analysis of the challenges of calling phasing of siRNA clusters in genome-wide analyses, see Polydore et al. (2018).

# OMIT  =====

# FURTHER STATS ====

# Number of sRNAs?

# Algunas secuencias de piRs son anotadas como miRs, ejemplo: 
# piR-bgl-168596;Esc-Mir-2-o20-v2_3p > Cluster_44687

RES <- RES %>%
  mutate(sRNA = "Novel") %>%
  mutate(sRNA = ifelse(grepl("piR-",KnownRNAs), "piR", sRNA)) %>%
  mutate(sRNA = ifelse(grepl("Mir-",KnownRNAs), "miR", sRNA))

RES %>% count(sRNA)
 
RES %>% 
  mutate(x = log10(Reads), y = log10(UniqueReads)) %>%
  ggplot(aes(x, y, color = sRNA, shape = Type)) +
  geom_point(size = 5, alpha = 0.7) +
  labs(caption = "Diff piR pattern expression from novel and miRs")

# separate know MIRS and Novel from piRs

str(query.seqs <- RES %>% 
    # filter(sRNA != "piR") %>% 
    distinct(Name) %>% 
    pull(Name))


# RES <- RESULTS %>% filter(Name %in% query.seqs)

which_cols <- c("Short", 21:30, "Long")

out <- RESULTS %>% 
  pivot_longer(cols = all_of(which_cols), names_to = "length") %>%
  select(Name, Chrom, length) %>% 
  filter(Name %in% query.seqs) %>%
  left_join(RES) 

out %>%
  group_by(length, sRNA) %>% tally() %>%
  ggplot(aes(x = length, y = n, fill = sRNA)) +
  geom_col(width = 0.85,position = "dodge2")


RESULTS %>% 
  filter(Name %in% query.seqs) %>%
  distinct(MajorRNA) %>% 
  arrange(MajorRNA) %>% 
  mutate(n = nchar(MajorRNA)) %>% 
  ggplot(aes(n)) + geom_histogram()

# revisar la anotacion de los 978 sRNAs (i.e 281 secuencias unicos) 
# y cotejar su anotacion

df <- RES %>% 
  drop_na(KnownRNAs) %>%
  mutate(KnownRNAs = strsplit(as.character(KnownRNAs), ";")) %>% 
  unnest(KnownRNAs) %>% distinct(KnownRNAs ,.keep_all = T)
  # filter(Name %in% query.seqs) 

df %>% count(sRNA)


# ENCONTRAR TABLA DE ESPECIE Y FUNCION CONOCIDA DE KNOWN MIRNAS: https://mirgenedb.org/

# sRNA2annot <- split(df$KnownRNAs, df$Name)



# MULTI-MAP DETAILS

# Reads: Number of aligned sRNA-seq reads that overlap this locus.
# UniqueReads: Number of uniquely aligned (e.g. not multi-mapping) reads that overlap this locus.

# MajorRNA: Sequence of the single most abundant RNA sequence at the locus.
# MajorRNAReads: Number of reads for the MajorRNA aligned to this locus.

RESULTS %>% filter(Length > 200) %>% view()

RESULTS %>%
  select(Name, Length, Reads, UniqueReads, MajorRNAReads) %>%
  pivot_longer(cols = c(Reads, UniqueReads, MajorRNAReads)) %>%
  ggplot(aes(Length, log10(value), color=name)) + geom_point()
  # ggplot(aes(log10(value), fill = name)) +  geom_histogram()
