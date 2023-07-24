# RICARDO GOMEZ-REYES
# 1) READ GFF FROM GENOME 
# SELECT PROTEIN_CODING COORDS (UTR IS EMBED HERE)
# CREATE A DF OBJECT OF NAME gene_features W/ COLS:
# "seqnames", "gene_coords" "gene_id", "description" "type", "biotype"


# 2) READ GTF FROM GENOME 
# SELECT PROTEIN_CODING COORDS (ACCORDING TO GFF)
# # CREATE A DF OBJECT W/ COLS:
# # "seqnames", "gene_coords" "gene_id", transcript_id,  "type", "transcript_biotype"
# USE .cds.all.fa BLASTX (diamont) RESOURCE TO RETRIVE ALL GENE ANNOTATION  
# SAVE OUT AS list

# gene_features (31171 genes) contains description related to gene_id (used as targets for srna regulatory analysis)
# transcript_features (55609 isoforms) contains key values for gene_id and transcript_id (used for trinotate)

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/ENSEMBLE/"


# GET gene features from locus ====

pattern <- "multi_genome.newid.gff3$"

genome <- list.files(path = wd, pattern = pattern, full.names = T)

.gr <- rtracklayer::import(genome)

# .gr %>% as_tibble() %>% dplyr::count(type)
# .gr %>% as_tibble() %>% drop_na(description) %>% dplyr::count(biotype)

# ...
# protein_coding 31171 # <---
# ...
#

gr <- .gr[which(.gr$type == "gene")]

gr <- gr[which(gr$biotype == "protein_coding")]

# Sanity check
# 31 171 of protein_coding with gene/protein description:

gr %>% as_tibble() %>% drop_na(description) %>% dplyr::count(biotype) 


which_cols <- c("seqnames", "gene_coords", "gene_id","description", "type", "biotype")

gene_features <- gr %>% as_tibble() %>% 
  drop_na(description) %>%
  mutate(gene_coords = paste(start, end, strand, sep = ":")) %>%
  select(any_of((which_cols)))

# Sanity check

gene_features %>% distinct(gene_id) %>% nrow() # 31 171


# 2) ====

# print(read_rds(paste0(wd, "SRNA_FUNCTION_PREDICTED.rds")))

pattern <- "multi_genome.newid.gtf$"

genome <- list.files(path = wd, pattern = pattern, full.names = T)

.gr <- rtracklayer::import(genome)

# .gr %>% as_tibble() %>% dplyr::count(type)

# 5 start_codon      55645
# 6 stop_codon       55603
#  transcript       55609 (protein_coding)

gr <- .gr[which(.gr$type == "transcript")]

gr <- gr[which(gr$transcript_biotype == "protein_coding")]

# Sanity check
# MUST MATCH THE SIZE 
gr %>% as_tibble() %>% dplyr::count(transcript_biotype) # 55609

transcript_features <- gr %>% as_tibble() %>% 
  mutate(gene_coords = paste(start, end, strand, sep = ":")) %>%
  select(seqnames, gene_coords, gene_id, transcript_id, type, transcript_biotype)

nrow(gene_features) # 31171 genes 
nrow(transcript_features) # 55609 isoforms (contained in the cds.all.fa)

out <- list(transcript_features, gene_features)

write_rds(out, file = paste0(wd, "genome_features.rds"))
