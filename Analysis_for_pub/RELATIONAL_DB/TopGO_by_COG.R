
# PURPOSE ----
# QUERY THEE DEGS/WGCNA MODULES
# 1) GENERAL EXPRESSION PATTERNS OF LARVAL DEV.
# 2) DEGS FOR Ctrl 24 AND 110 HPF
# 2.1) JOIN Ctrl DEGS to LONGER_RELATIONAL_DB.tsv
# 3) Join Ctrl DEGS TO pH 7.6 DEGS
# PLOT Enrichment Ratio of NOGS
# 4) TOP GO ENRICHMENT
# 4.1) JOIN LONGER_RELATIONAL_DB TO DEGS
# 4.2) RUN TOP GO ENRICHMENT

library(tidyverse)

rm(list = ls());

if(!is.null(dev.list())) dev.off()

dir <- "~/Documents/MIRNA_HALIOTIS/MIRNA_PUB_2024/"

print(DB <- read_tsv(file.path(dir, "LONGER_RELATIONAL_DB.tsv")))

# TOP GO ENRICHMENT by NOG ----
dir <- "~/Documents/MIRNA_HALIOTIS/MIRNA_PUB_2024/"

GODB <- read_tsv(file.path(dir, "Gene_ontologies_DB.tsv"))

GODB <- GODB %>% drop_na(GOs) %>% filter(GOs != "-") %>% 
  distinct(gene_id, GOs)

paste_col <- function(x) { 
  x <- x[!is.na(x)] 
  x <- unique(sort(x))
  x <- paste(x, sep = ';', collapse = ';') 
  
  return(x)
  
}


gene2GO <- DB %>% 
  distinct(gene_id, COG_name) 

gene2GO <- GODB %>% 
  mutate(GOs = strsplit(GOs, ",")) %>%
  unnest(GOs) %>%
  left_join(gene2GO) %>%
  distinct() %>%
  group_by(COG_name, gene_id) %>%
  summarise(across(GOs, .fns = paste_col)) %>%
  drop_na(COG_name, GOs) %>% 
  ungroup()

# split(strsplit(GODB$GOs, ";") , GODB$GROUP)

runTopGO <- function(x, y, Nodes = 20, onto = "BP") {
  
  # allGenes <- y %>% filter(CONTRAST == x) %>% pull(pvalue, name = MajorRNA)
  
  cat("Running ", x)
  
  allGenes <- y %>% filter(COG_name == x) %>% pull(gene_id)
  
  allGenes <- structure(rep(0.05, length(allGenes)), names = allGenes)
  
  # gene2GO <- gene2GO %>% filter(CONTRAST == x)
  
  gene2GO <- split(strsplit(gene2GO$GOs, ";") , gene2GO$gene_id)
  
  gene2GO <- lapply(gene2GO, unlist)
  
  keep <- names(gene2GO) %in% names(allGenes)
  
  gene2GO <- gene2GO[keep]
  
  keep <- names(allGenes) %in% names(gene2GO)
  
  allGenes <- allGenes[keep]
  
  description <- "Run topGO enrichment"
  
  library(topGO)
  
  topGOdata <- new("topGOdata", 
    ontology = onto, 
    description = "description",
    allGenes = allGenes,
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

CONTRAST2GO <- gene2GO %>% distinct(COG_name) %>% pull()

# Omit single gene_id
patt <- "Cell wall/membrane/envelope biogenesis"
CONTRAST2GO <- CONTRAST2GO[!grepl(patt, CONTRAST2GO)]

OUT <- lapply(CONTRAST2GO, function(x) runTopGO(x, gene2GO, Nodes = 25))

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

write_tsv(data, file = file.path(dir, "Boostrap_topGO_by_COG.tsv"))

data %>%
  mutate(classicKS = as.double(classicKS), 
    elimKS = as.double(elimKS)) %>%
  ggplot(aes(y = GROUP, x = parentTerm, fill = -log10(classicKS))) +
  geom_tile(color = 'white', linewidth = 0.2) +
  # facet_grid(~ f1+f2, switch = "y", scales = "free_y") +
  theme_bw(base_family = "GillSans", base_size = 10) +
  theme(axis.text.x = element_text(angle = 90, size = 10, hjust = 1))
