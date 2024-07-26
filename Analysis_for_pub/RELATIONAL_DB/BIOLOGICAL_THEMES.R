
# PROCESS EGGNOGMAPPER + BLAST-DIAMONT RESULTS

# Define biological-phentypic themes
# Calcification <---
# Development
# Growth
# RM (Respiratory Metabolism)
# Group transcriptome to these themes
# Idealy, using only miRNA:mRNA transcripts (~ 170 transcripts)
# Additionally, blast miRNA:mRNA transcripts to biomineralization-genes database

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)


dir <- "~/Documents/MIRNA_HALIOTIS/ENSEMBLE/calcif/"

f <- list.files(dir, pattern = ".outfmt6", full.names = T)

# selecting only biomineralization
f <- f[grepl("Combined_Gastropoda|shell_matrix", f)]

library(tidyverse)

read_outfmt6 <- function(f) {
  
  # seqid = transcript_id
  outfmt6.names <- c("transcript_id", "subject", "identity", "coverage", "mismatches", "gaps", "seq_start", "seq_end", "sub_start", "sub_end", "e", "score")
  
  
  df <- read_tsv(f, col_names = F) 
  
  colnames(df) <- outfmt6.names
  
  df <- df %>% mutate(db = basename(f))
  
  return(df)
  
  
}

SMP_DB <- do.call(rbind, lapply(f, read_outfmt6))

# df %>% ggplot(aes(identity, color = db)) + stat_ecdf()

# merge to functional db
# gene2tr -----
wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/ENSEMBLE/"

gene2tr <- read_rds(paste0(wd, "genome_features.rds"))[[1]] %>% 
  distinct(gene_id, transcript_id)

SMP_DB <- SMP_DB %>% distinct(transcript_id, subject) %>% left_join(gene2tr)

wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/FUNCTIONAL_MIR_ANNOT/"
# 
print(FUNCTIONAL_DB <- read_rds(paste0(wd,"SRNA_FUNCTION_PREDICTED_LONG_EXPRESSED.rds")))
# 
# df <- df %>% left_join(FUNCTIONAL_DB, by = "gene_id")
# 
# df %>% dplyr::select(dplyr::starts_with("SRR")) %>% as("matrix") %>% heatmap(Colv = NA)

# MajorRNA2Name -----

wd <- "~/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out/"

# bind to Loci mirna db

print(.LOCATION_DB <- read_rds(paste0(wd, "RNA_LOCATION_MIR_DB.rds")))

MajorRNA2Name <- .LOCATION_DB %>% distinct(MajorRNA, Name) %>%   dplyr::rename("query" = "Name")

# EGGMAPPER -----
# THIS OUTPUT IS DERIVED FROM ALL miRNA-mRNA TARGETS (~170 GENES) passed through eggnogmapper
# 
# f <- "/Users/cigom/Documents/GitHub/Small-RNASeq-data-analysis/NOG.annotations.tsv"
# NOG.annotations <- read.table(f, sep="\t", stringsAsFactors=FALSE, quote="")
# names(NOG.annotations) <- c("db", "nog", "proteins", "species", "code", "name")
# NOG.annotations <- NOG.annotations %>% distinct(code, name) %>% as_tibble()


ref_path <- "/Users/cigom/Documents/MIRNA_HALIOTIS/ENSEMBLE/ANNOTATIONS/"

NOG.col <- read_rds(paste0(ref_path, '/cogs.rds'))


# NOG.col <- NOG.annotations %>% mutate(clrs = "#CDC") %>%
#   anti_join(NOG.col, by = "code") %>%
#   rbind(NOG.col)
  
NOG.col

dir <- "~/Documents/MIRNA_HALIOTIS/ENSEMBLE/calcif/"

f <- list.files(dir, pattern = "eggnog_mapper.emapper.annotations", full.names = T)

# the eggnog_mapper.emapper.annotations file is the good one
eggnog_db <- read_tsv(f, comment = "##")

# Prior to EDA, bind data bases
# Bind keyids data bases ----
# 
 
FUNCTIONAL_DB
  # mutate(query = strsplit(query, ";")) %>%
  # unnest(query) %>% 
  # distinct(gene_id, query)

# 1) bind transcript_id to gene_id and 
# collapse the transcript_id
DB <- gene2tr %>%
  right_join(eggnog_db, by = c("transcript_id" = "#query" )) %>%
  select(- transcript_id) %>% distinct()

# 1.2) add kegg theme

DB <- DB %>% left_join(NOG.col, by = c("COG_category" = "code"))

DB <- DB %>% rename("COG_name" = "name")

view(DB)

DB %>% distinct(Preferred_name) %>% view() 

# 2) bind miRNA info (ONLY TO SPLIT DIFFEXP THEMES)
# Note: duplicate genes means different target location
# duplicates=  many-to-many relationship 

FUNCTIONAL_DB %>% distinct(gene_id, query) %>%
  mutate(query = strsplit(query, ";")) %>%
  unnest(query) %>% 
  mutate(query = gsub(".mature$", "", query)) %>%
  distinct() %>%
  left_join(MajorRNA2Name) %>%
  select(-query) 
  # left_join(DB, by = "gene_id") %>% view()




# EDA eggnog_mapper ----
# n <- match(eggnog_db$COG_category, NOG.annotations$class)
# y <- table(unlist(strsplit(NOG.annotations$class[n], "")))
# NOG.annotations$description[n]
# y <- data.frame(y)
# names(y) <- c('code', 'Freq')
# sum(y$Freq)


col_palette <- NOG.col %>% distinct(code, clrs) %>% pull(clrs, name = code)

eggnog_db %>%  distinct(`#query`)

plotdf <- eggnog_db %>% count(COG_category) %>% 
  left_join(NOG.col, by = c("COG_category" = "code")) %>%
  drop_na(name) %>% 
  arrange(desc(n))


plotdf %>%
  mutate(name = factor(name, levels = unique(name))) %>%
  ggplot(aes(x = name, y = n, fill = COG_category)) + # 
  geom_col() +
  coord_flip() +
  geom_text(aes(label = n), size = 3, hjust = -0.05, family = "GillSans") +
  scale_fill_manual(values = col_palette) 


# merge to gene_id

eggnog2gene <- eggnog_db %>% distinct(`#query`, COG_category) %>% 
  left_join(NOG.col, by = c("COG_category" = "code")) %>%
  left_join(gene2tr, by = c("#query" = "transcript_id")) %>% 
  distinct(gene_id, name, COG_category) 

# eggnog2gene %>% 
#   count(gene_id, name, COG_category, sort = T) %>%  drop_na(name) %>% 
#   mutate(name = factor(name, levels = unique(name))) %>%
#   ggplot(aes(x = n, y = name)) +
#   ggridges::geom_density_ridges()

eggnog2gene %>% 
  count(name, COG_category, sort = T) %>%  drop_na(name) %>% 
  mutate(name = factor(name, levels = unique(name))) %>%
  ggplot(aes(x = name, y = n, fill = COG_category)) + # 
  geom_col() +
  coord_flip() +
  geom_text(aes(label = n), size = 3, hjust = -0.05, family = "GillSans") +
  scale_fill_manual(values = col_palette)

distinct(FUNCTIONAL_DB, gene_id, query) %>%
  right_join(eggnog2gene, by = "gene_id") %>% 
  filter(name != "Function unknown") %>%
  mutate(query = strsplit(query, ";")) %>%
  unnest(query) %>% 
  distinct(query, name)


eggnog_db %>% distinct(`#query`, GOs) %>%
  left_join(gene2tr, by = c("#query" = "transcript_id")) %>% 
  distinct(gene_id, GOs)


# Enrichment analysis

# clusterprof

# BiocManager::install("clusterProfiler")
library(clusterProfiler)

hrf <- search_kegg_organism("Haliotis rufescens", by='scientific_name')

hrf

kegg_organism = "hrf" # cel hsa

geneList <- eggnog_db %>% 
  mutate(KEGG_ko = strsplit(KEGG_ko, ",")) %>%
  distinct(KEGG_ko, score) %>%
  unnest(KEGG_ko) %>% mutate(KEGG_ko = gsub("ko:K", "", KEGG_ko)) %>%
  arrange(score) %>% filter(KEGG_ko != "-") %>%
  pull(score, name = KEGG_ko)

kk2 <- gseKEGG(geneList = geneList,
  organism     = kegg_organism,
  nPerm        = 10000,
  minGSSize    = 3,
  maxGSSize    = 800,
  pvalueCutoff = 0.05,
  pAdjustMethod = "none",
  keyType       = "kegg")

