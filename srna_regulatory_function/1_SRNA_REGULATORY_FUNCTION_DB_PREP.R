# RICARDO GOMEZ-REYES
# READ GFF FROM GENOME 
# SELECT PROTEIN_CODING COORDS (UTR IS EMBED HERE)
# CREATE A DF OBJECT W/ COLS:
# [1] "seqnames"    "gene_coords" "gene_id"     "description" "type"        "biotype"


rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/FUNCTIONAL_MIR_ANNOT"

# GET gene features from locus ====

pattern <- "multi_genome.newid.gff3$"

genome <- list.files(path = wd, pattern = pattern, full.names = T)

.gr <- rtracklayer::import(genome)

# .gr %>% as_tibble() %>% dplyr::count(type)
# .gr %>% as_tibble() %>% drop_na(description) %>% dplyr::count(biotype)

# ...
# 4 protein_coding 31171 # <---
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

# save gene/protein description to:

write_rds(gene_features, file = paste0(wd, "gene_features.rds"))
