
# Q: ARE HOST-GENE FROM INTRAGENIC MIRS RELEVANT FOR DEVELOPMENT?

# FOR INTRAGENIC-LOCI MICRORNA

# 0) LOAD gene2tr DATABASE
# 2) LOAD EGGMAPPER RESULTS FROM miRNA:LOCI transcripts
# 3) LOAD STRING AND PARSE
# 4) JOIN DB
# 5 RUN TOPGO (SPLIT BY DEGS-TIME)


library(tidyverse)

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

read_outfmt6 <- function(f) {
  
  # seqid = transcript_id
  outfmt6.names <- c("transcript_id", "subject", "identity", "coverage", "mismatches", "gaps", "seq_start", "seq_end", "sub_start", "sub_end", "e", "score")
  
  
  df <- read_tsv(f, col_names = F) 
  
  colnames(df) <- outfmt6.names
  
  df <- df %>% mutate(db = basename(f))
  
  return(df)
  
  
}

out_dir <- "~/Documents/MIRNA_HALIOTIS/MIRNA_PUB_2024/"

# 0

dir <- "/Users/cigom/Documents/MIRNA_HALIOTIS/ENSEMBLE/"

genedescr <- read_rds(paste0(dir, "genome_features.rds"))[[2]]

gene2tr <- read_rds(paste0(dir, "genome_features.rds"))[[1]] %>% 
  distinct(gene_id, transcript_id)


dir <- "~/Documents/MIRNA_HALIOTIS/ENSEMBLE/ONLY_miRNA_loci_dir/"

mapper_f <- list.files(dir, pattern = "eggnog_mapper.emapper.annotations", full.names = T)

outfmt_f <- list.files(dir, pattern = "outfmt6", full.names = T)


# MAPPER ----


MAPPER_DB <- read_tsv(mapper_f, comment = "##")

eggNOG_cols <- c("gene_id", "Preferred_name","Description", "GO","PFAMs","BRITE", "CAZy", "COG_category")

MAPPER_DB <- gene2tr %>%
  right_join(MAPPER_DB, by = c("transcript_id" = "#query" )) %>%
  select(- transcript_id, -evalue, -score) %>%
  select_at(vars(contains(eggNOG_cols), starts_with("KEGG"))) %>%
  distinct()

MAPPER_DB %>% distinct(gene_id)

# Bind to the COG_category

dir <- "/Users/cigom/Documents/MIRNA_HALIOTIS/ENSEMBLE/ANNOTATIONS/"

NOG.col <- read_rds(paste0(dir, '/cogs.rds')) %>%
  dplyr::rename("COG_name" = "name")

MAPPER_DB <- MAPPER_DB %>% 
  left_join(NOG.col, by = c("COG_category" = "code"))

GOID <- MAPPER_DB %>% distinct(gene_id, GOs) %>% filter(GOs != "-") # %>% 
  # mutate(GOs = strsplit(GOs, ",")) %>%
  # unnest(GOs) %>%
  # group_by(gene_id) %>%
  # summarise(across(GOs, .fns = paste_col))

# TOPGO ----

print(DB <- read_tsv(file.path(dir, "LONGER_RELATIONAL_DB.tsv")))

DB %>%
  distinct(MajorRNA, gene_id) %>%
  right_join(GOID)

dir <- "~/Documents/MIRNA_HALIOTIS/MIRNA_PUB_2024/"

RES <- read_rds(list.files(path = dir, 
  pattern = "SEQUENCES_MERGED_DESEQ_RES_WGCNA.rds", full.names = T)) %>%
  # left_join(WGCNA) %>%
  mutate(Contrast = ifelse(sign(log2FoldChange) == -1, sampleB, sampleA)) %>%
  mutate(Contrast = paste0(CONTRAST,"_",Contrast))

RES.P <- RES %>% filter( padj < 0.05  & abs(log2FoldChange) > 1)


QUERYDB <- RES.P %>% 
  left_join(LOCIDB, by = "MajorRNA") %>% 
  mutate(dir = ifelse(sign(log2FoldChange) == -1, "-logFC","+logFC"), 
    CONTRAST = paste0(CONTRAST, ":", dir))

# TOPGO


runTopGO <- function(x, y, Nodes = 20, onto = "BP") {
  
  allMIRS <- y %>% filter(CONTRAST == x) %>% pull(padj, name = MajorRNA)
  
  # gene2GO <- gene2GO %>% filter(CONTRAST == x)
  gene2GO <- split(strsplit(gene2GO$GOs, ";") , gene2GO$MajorRNA)
  
  gene2GO <- lapply(gene2GO, unlist)
  
  keep <- names(gene2GO) %in% names(allMIRS)
  
  gene2GO <- gene2GO[keep]
  
  keep <- names(allMIRS) %in% names(gene2GO)
  
  allMIRS <- allMIRS[keep]
  
  description <- "Run topGO enrichment"
  
  library(topGO)
  
  topGOdata <- new("topGOdata", 
    ontology = onto, 
    description = "description",
    allGenes = allMIRS,
    geneSel = function(x) x,
    annot = annFUN.gene2GO,
    gene2GO = gene2GO)
  
  # Boostraping
  if(is.infinite(Nodes)) {
    cat("Using boostraping")
    
    allRes <- list()
    
    for (i in c(20, 50, 75, 100, length(allMIRS))) {
      
      DF <- runtopGO(topGOdata, topNodes = i, conservative = T)
      
      allRes[[i]] <- data.frame(DF, Top = i)
    }
    
    allRes <- do.call(rbind, allRes)
    
  } else {
    
    cat("Using topNodes lim")
    
    allRes <- runtopGO(topGOdata, topNodes = Nodes, conservative =T)
    
    allRes <- data.frame(allRes, Top = Nodes)
  }
  
  out <- data.frame(allRes, GROUP = x)
  
  return(out)
}

runTopGO(, DB, Nodes = 50)

OUT1 <- lapply(CONTRAST2GO1, function(x) runTopGO(x, QUERYDB1, Nodes = 50))


# BLAST

BLAST <- do.call(rbind, lapply(outfmt_f, read_outfmt6))

