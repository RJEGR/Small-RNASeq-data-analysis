
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
  summarise(across(GOs, .fns = paste_col), n = n()) %>%
  drop_na(COG_name, GOs) %>% 
  ungroup()

runTopGO <- function(x, y, Nodes = 20, onto = "BP") {
  
  # allGenes <- y %>% filter(CONTRAST == x) %>% pull(pvalue, name = MajorRNA)
  
  cat("Running ", x)
  
  allGenes <- y %>% filter(COG_name == x) %>% pull(gene_id)
  
  allGenes <- structure(rep(0.05, length(allGenes)), names = allGenes)
  
  gene2GO <- y %>% filter(COG_name == x)
  
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

omit_single_genes <- gene2GO %>% 
  group_by(COG_name) %>% 
  summarise(ngo = sum(n), ngenes=n()) %>%
  filter(ngenes == 1) %>% pull(COG_name) %>%
  paste(collapse = "|")

CONTRAST2GO <- gene2GO %>% 
  distinct(COG_name) %>% pull()

CONTRAST2GO <- CONTRAST2GO[!grepl(omit_single_genes, CONTRAST2GO)]

OUT <- lapply(CONTRAST2GO, function(x) runTopGO(x, gene2GO, Nodes = 25))

print(data <- dplyr::bind_rows(OUT) %>% as_tibble())

data %>% count(GROUP)
# 
write_tsv(data, file = file.path(dir, "COG2topGO.tsv"))

str(GO.IDS <- data %>% distinct(GO.ID) %>% pull() %>% sort())

wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/FUNCTIONAL_MIR_ANNOT/"

orgdb <- "org.Hs.eg.db" # "org.Ce.eg.db"

# semdata <- GOSemSim::godata(orgdb, ont="BP")
# write_rds(semdata, file = paste0(wd, orgdb, ".rds"))

semdata <- read_rds(paste0(wd, orgdb, ".rds"))

SEM <- SEMANTIC_SEARCH(GO.IDS, orgdb = orgdb, semdata = semdata)

data <- data %>% left_join(SEM, by = c("GO.ID" = "go"))

write_tsv(data, file = file.path(dir, "COG2topGO.tsv"))

data %>%
  group_by(GROUP, parentTerm) %>%
  summarise(score = sum(Annotated)) %>%
  mutate(score = score/sum(score)) %>%
  # mutate(classicKS = as.double(classicKS), 
  #   elimKS = as.double(elimKS)) %>%
  ggplot(aes(y = GROUP, x = parentTerm, fill = score)) +
  geom_tile(color = 'white', linewidth = 0.2) +
  # facet_grid(~ f1+f2, switch = "y", scales = "free_y") +
  theme_bw(base_family = "GillSans", base_size = 10) +
  theme(axis.text.x = element_text(angle = 90, size = 10, hjust = 1))

# due to omit_single_genes
# run semantic grouping

cogs2GO <- gene2GO %>% 
  mutate(GOs = strsplit(GOs, ";")) %>%
  unnest(GOs) %>%
  group_by(COG_name) %>%
  summarise(across(GOs, .fns = paste_col)) 


cogs2GO <- split(strsplit(cogs2GO$GOs, ";") , cogs2GO$COG_name)
cogs2GO <- lapply(cogs2GO, unlist)

OUT <- lapply(cogs2GO, function(x) SEMANTIC_SEARCH(x, semdata = semdata))

print(DATA <- dplyr::bind_rows(OUT, .id = "COG_name") %>% as_tibble())

write_tsv(DATA, file = file.path(dir, "COG2revigo.tsv"))

DATA %>% distinct(parentTerm)

DATA %>% 
  count(COG_name, parentTerm) %>%
  group_by(COG_name) %>%
  mutate(size = n/sum(n)) %>% 
  # filter(size > 0.2) %>%
  ggplot(aes(y = parentTerm, x = size)) +
  geom_col() +
  facet_grid(COG_name ~., scales = "free", space = "free") +
  geom_tile(color = 'white', linewidth = 0.2) +
  # facet_grid(~ f1+f2, switch = "y", scales = "free_y") +
  theme_bw(base_family = "GillSans", base_size = 10) +
  theme(
    strip.text.y = element_text(angle = 0, hjust = 1),
    strip.background = element_rect(fill = 'grey95', color = 'white'),
    # strip.placement = "outside",
    axis.text.y = element_text(size = 5, hjust = 1),
    axis.text.x = element_text(angle = 0, size = 5, hjust = 1)) +
  ggsave(filename = 'COGs2revigo.png', 
    path = dir, width = 10, height = 10, device = png, dpi = 300)

DATA %>%
  group_by(parentTerm) %>%
  # summarise(size = sum(size)) %>%
  mutate(size = size/sum(size)) %>%
  ggplot(aes(y = parentTerm, x = size)) +
  geom_col() +
  facet_grid(COG_name ~., scales = "free", space = "free") +
  geom_tile(color = 'white', linewidth = 0.2) +
  # facet_grid(~ f1+f2, switch = "y", scales = "free_y") +
  theme_bw(base_family = "GillSans", base_size = 10) +
  theme(axis.text.x = element_text(angle = 90, size = 10, hjust = 1))

