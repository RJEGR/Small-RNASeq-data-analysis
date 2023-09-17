
# INSTEAD OF IT, USE CDS.ALL.FA 
# SUBSEQ SEQUENCE BASED ON COORDINATES CLOSE TO UTR-TARGET SITES
# PREPARE KNOW ANNOTATION SEQUENCES TO RERUN BLAST AND RETRIVE GENE ONTOLOGIES

library(tidyverse)

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/FUNCTIONAL_MIR_ANNOT"

pattern <- "/mir_vs_utr_rmdup_RNAhybrid.out.psig.gene_features.rds"

df_drop <- read_rds(paste0(wd, pattern))

df_drop %>%
  drop_na(gene_id) %>%
  mutate(query = strsplit(query, ";")) %>%
  unnest(query) %>% distinct(query)

gene_ids <- df_drop %>% drop_na(gene_id) %>% distinct(gene_id) %>% pull(gene_id)

# fileName <- paste0(gsub(".rds", "", pattern), "_known_annot.ids")

# write_lines(gene_ids, file = paste0(wd, fileName))

# cat Haliotis_rufescens_gca023055435v1rs.xgHalRufe1.0.p.cds.all.fa | seqkit grep -f $file

f <- "Haliotis_rufescens_gca023055435v1rs.xgHalRufe1.0.p.cds.all.fa"

dna_wd <- "~/Documents/MIRNA_HALIOTIS/ENSEMBLE"

f <- list.files(path = dna_wd, pattern = f, full.names = T)

dna <- Biostrings::readDNAStringSet(f)

str(transcript_id <- sapply(strsplit(names(dna), " "), `[`, 1))

str(unique(transcript_id)) # check if unique identifier

head(gene_id <- sapply(strsplit(names(dna), " "), `[`, 4))

str(gene_id <- gsub("gene:","", gene_id))

str(unique(gene_id))

query_ids <- data.frame(gene_id, transcript_id) %>% 
  as_tibble() %>%
  filter(gene_id %in% gene_ids) %>%
  distinct(transcript_id) %>% pull()

str(query_ids) # 15718

names(dna) <- transcript_id

sum(keep <- names(dna) %in% query_ids) # 15718

fileName <- paste0(gsub(".rds", "", pattern), "_known_annot.fa")

Biostrings::writeXStringSet(dna[keep], file = paste0(wd, fileName))

# re-run mir_vs_utr_rmdup_RNAhybrid.out.psig.gene_features_known_annot.fa w/ blast diamont
# uniprot_sprot.ncbi.blastx.tsv

# LOAD uniprot_sprot.ncbi.blastx.tsv
# MAP TO http://current.geneontology.org/ontology/go-basic.obo

# ~/Trinotate-Trinotate-v4.0.0/util/admin/util/obo_to_tab.pl go-basic.obo > obo_to_tab.tsv

url <- "http://current.geneontology.org/ontology/go-basic.obo"

OBO <- readLine(url)

OBO <- read_tsv(url,col_names = F)

OBO %>% head(100) %>% view()

# OBO %>% filter(grepl("^[Term]", X1))

# 
# my ($field, $annot) = split(/\s+/, $_, 2);