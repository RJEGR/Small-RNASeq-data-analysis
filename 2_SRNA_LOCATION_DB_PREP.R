
# RICARDO GOMEZ-REYES  
# AFTER RUN 1_SRNA_LOCATION_DB_PREP.R
# 

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

library(GenomicRanges)

wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out"

head(.out <- read_rds(paste0(wd, "/SRNA_LOCATION_OUT.rds")))

.out <- do.call(rbind, .out) %>% as_tibble()

.out %>% count(type)

.out %>% count(Name)


paste_headers <- function(x) { 
  x <- x[!is.na(x)] 
  n <- length(x)
  x <- unique(x)
  
  x <- paste(x, sep = '|', collapse = '|') 
  x <- paste(n, x, sep = '|', collapse = '|') }

# .out %>% filter(Name == "Cluster_11")

out <- .out %>%
  head(100) %>%
  group_by(Name, biotype) %>%
  summarise(across(type, .fns = paste_headers), .groups = "drop_last") %>% 
  ungroup()

out %>% drop_na(gene_id)

# MERGE W/ SHORTSTACKS OUTPUT =====

# SRNA-SEQ TRANSCRIPTOME
wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out"

pattern <- "Results.gff3"

f <- list.files(path = wd, pattern = pattern, full.names = T)

assembly <- rtracklayer::import(f)

# SPLIT TRUE MIRS/PIRS

SRNAS <- read_rds(paste0(wd, "/KNOWN_CLUSTERS_MIRS_PIRS.rds"))

str(which_pirs <- SRNAS %>% filter(grepl("piR", KnownRNAs)) %>% distinct(Name) %>% pull())
str(which_mirs <- SRNAS %>% filter(!grepl("piR", KnownRNAs)) %>% distinct(Name) %>% pull())

# LOAD MODULES ?? ====

bwnet <- readRDS(paste0(wd, "/2023-06-26/bwnet.rds"))

bwmodules = WGCNA::labels2colors(bwnet$colors)

names(bwmodules) <- names(bwnet$colors)

assembly$Module <- NA

x <- names(bwmodules) # head(names(bwmodules))

lab <- bwmodules # head(bwmodules)

assembly$Module <- NA

assembly$Module <- bwmodules[match(assembly$ID, names(bwmodules))]

# table(bwmodules)

# PREPARE DB FOR VIZ ====
# 1)
# TO KEEP TRACK LOCUS TYPE PER CLUSTER
# Note: be aware score (Reads) cols for *.mature/*.star locus are documented in assembly  (262 rows)

# assembly %>% as_tibble() %>% filter(grepl(".mature|.star", ID)) %>% select()
  
.DB1 <- assembly %>% as_tibble() %>% select(ID, type) %>% rename("Locus_type"="type", "Name"="ID")

# 2)  BIND SEVERAL OUTPUT  ====
# TO KEEP TRACK READS COUNT PER CLUSTER

f <- list.files(path = wd, pattern = "Results.txt", full.names = T)

read_tsv(f) %>% arrange(desc(FracTop))
  head() %>% View()

DB <- read_tsv(f) %>% 
  select(Locus, Name, Chrom, Start, End, Strand, DicerCall, Reads, UniqueReads, FracTop, 
    MajorRNA, MajorRNAReads, MIRNA) %>%
  left_join(.DB1, by = "Name")


# 3) BIND KNOWN RNA LABEL =====

.DB <- SRNAS %>% select(-MajorRNA) %>% rename("NKnownRNAs"="n")

DB <- DB %>% left_join(.DB, by = "Name")



# 4) BIND LOCATION DB ====

# HAY QUE REVISAR QUE EFECTO TIENE STRAND EN EL ANALISIS DE OVERLAPS, 
# PUES NO COINCIDEN LAS COORDENADAS GENOMICAS (.out) CON LAS DE LA BASE DE DATOS (DB)

# DB %>% drop_na(KnownRNAs) %>% count(Locus_type, Strand, sort = T) %>% view()

bind_types <- function(x) { x <- paste(unique(x[!is.na(x)]), sep = '|', collapse = "|") }

bind_ids <- function(x) { 
  
  x <- x[!is.na(x)]
  
  x <- unique(x)
  x <- paste(x, sep = ';', collapse = ';') }


.DB <- .out %>%
  group_by(Name) %>%
  summarise(
    across(type, .fns = bind_types), 
    across(biotype, .fns = bind_types),
    across(gene_id, .fns = bind_ids),
    across(transcript_id, .fns = bind_ids),
    .groups = "drop_last") %>% 
  ungroup() 

DB <- DB %>% left_join(.DB, by = "Name")

write_tsv(DB, file = paste0(wd, "/RNA_LOCATION_DB.tsv"))

# head(read_tsv(paste0(wd, "/RNA_LOCATION_DB.tsv")))


