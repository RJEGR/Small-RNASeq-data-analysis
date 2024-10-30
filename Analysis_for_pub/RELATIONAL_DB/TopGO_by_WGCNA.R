
# RUN BY WGCNA 

library(tidyverse)

rm(list = ls());

if(!is.null(dev.list())) dev.off()


dir <- "~/Documents/MIRNA_HALIOTIS/MIRNA_PUB_2024/"

print(DB <- read_tsv(file.path(dir, "LONGER_RELATIONAL_DB.tsv")))


RES <- read_rds(list.files(path = dir, 
  pattern = "SEQUENCES_MERGED_DESEQ_RES_WGCNA.rds", full.names = T)) %>%
  # left_join(WGCNA) %>%
  mutate(Contrast = ifelse(sign(log2FoldChange) == -1, sampleB, sampleA)) %>%
  mutate(Contrast = paste0(CONTRAST,"_",Contrast))

RES.P <- RES %>% filter( padj < 0.05  & abs(log2FoldChange) > 1)

# Previous version use Ontologies from trinotate
# wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/FUNCTIONAL_MIR_ANNOT/"
# print(.SRNA2GO <- read_tsv(paste0(wd, "SRNA_REGULATORY_FUNCTION_DB.tsv")))

# upgraded version use from eggnogg mapper, in addition to the subset of expressed genes during larval dev.

GODB <- read_tsv(file.path(dir, "Gene_ontologies_DB.tsv"))

GODB <- GODB %>% drop_na(GOs) %>% filter(GOs != "-") %>% 
  distinct(gene_id, GOs)

# QUERY DEGS PRIOR TO SPLIT
# QUERY ONLY MIRS CONSERVED IN Ctrl

QUERYDB <- RES.P %>% filter(CONTRAST == "CONTRAST_C") %>%
  mutate(HPF = ifelse(sign(log2FoldChange) == -1, "110 hpf", "24 hpf"))

# CAUTION: GENES ARE MULTI-TARGETS

DB %>%
  distinct(gene_id, MajorRNA) %>%
  count(MajorRNA, sort = T)

paste_col <- function(x) { 
  x <- x[!is.na(x)] 
  x <- unique(sort(x))
  x <- paste(x, sep = ';', collapse = ';') 
  
  return(x)
  
}


gene2GO <- DB %>% 
  right_join(QUERYDB, by = "MajorRNA") %>%
  distinct(gene_id, MajorRNA) 


gene2GO <- GODB %>% 
  mutate(GOs = strsplit(GOs, ",")) %>%
  unnest(GOs) %>% 
  left_join(gene2GO) %>%
  distinct(MajorRNA, GOs) %>%
  group_by(MajorRNA) %>%
  summarise(across(GOs, .fns = paste_col))

gene2GO <- QUERYDB %>% distinct(MajorRNA, HPF) %>% 
  right_join(gene2GO)

# split(strsplit(GODB$GOs, ";") , GODB$GROUP)

runTopGO <- function(x, y, Nodes = 20, onto = "BP") {
  
  allMIRS <- y %>% filter(CONTRAST == x) %>% pull(pvalue, name = MajorRNA)
  
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

QUERYDB <- RES.P %>% 
  mutate(dir = ifelse(sign(log2FoldChange) == -1, "-logFC","+logFC"), 
    CONTRAST = paste0(CONTRAST, ":", dir))

CONTRAST2GO <- QUERYDB %>% distinct(CONTRAST) %>% pull()

OUT <- lapply(CONTRAST2GO, function(x) runTopGO(x, QUERYDB, Nodes = Inf))

print(data <- dplyr::bind_rows(OUT) %>% as_tibble())
# 
# write_tsv(data, file = file.path(dir, "Boostrap_topGO.tsv"))

str(GO.IDS <- data %>% distinct(GO.ID) %>% pull() %>% sort())

wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/FUNCTIONAL_MIR_ANNOT/"

orgdb <- "org.Hs.eg.db" # "org.Ce.eg.db"

# semdata <- GOSemSim::godata(orgdb, ont="BP")
# write_rds(semdata, file = paste0(wd, orgdb, ".rds"))

semdata <- read_rds(paste0(wd, orgdb, ".rds"))

SEM <- SEMANTIC_SEARCH(GO.IDS, orgdb = orgdb, semdata = semdata)

data <- data %>% left_join(SEM, by = c("GO.ID" = "go"))

write_tsv(data, file = file.path(dir, "Boostrap_topGO.tsv"))


recode_c <- c("CONTRAST_A:-logFC" = "24 hpf:pH 7.6",
  "CONTRAST_B:+logFC" = "110 hpf:pH 8.0",
  "CONTRAST_B:-logFC" = "110 hpf:pH 7.6",
  "CONTRAST_C:+logFC" = "24 hpf:pH 8.0",
  "CONTRAST_C:-logFC"= "110 hpf:pH 8.0",
  "CONTRAST_D:+logFC"= "24 hpf:pH 7.6",
  "CONTRAST_D:-logFC" = "110 hpf:pH 7.6")


data <- data %>%
  mutate(CONTRAST = dplyr::recode_factor(GROUP, !!!recode_c)) %>%
  mutate(Top = ifelse(!Top %in% c(20, 50, 75, 100), "All", Top)) %>%
  mutate(Top = factor(as.character(Top), levels = c(20, 50, 75, 100, "All"))) %>%
  group_by(CONTRAST, Top) %>%
  mutate(Ratio = size / max(size)) %>%
  filter(Ratio > 0) %>%
  separate(CONTRAST, into = c("f1", "f2"), sep = ":")

write_tsv(data, file = file.path(dir, "Boostrap_topGO.tsv"))


