
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

# 1) GENERAL EXPRESSION PATTERNS -----
# 
# plotdf <- DB %>% count(COG_category, COG_name, sort = T) %>% 
#   drop_na() 
# 
# 
# plotdf %>%
#   mutate(COG_name = factor(COG_name, levels = unique(COG_name))) %>%
#   ggplot(aes(x = COG_name, y = n, fill = COG_category)) + # 
#   geom_col() +
#   coord_flip() +
#   geom_text(aes(label = n), size = 3, hjust = -0.05, family = "GillSans")
#   # scale_fill_manual(values = col_palette) 
# 

# 2) PREP. DEGS FOR NORMAL 24 AND 110 HPF ----

# 89 DEGS IN NORMAL

dir <-  "~/Documents/MIRNA_HALIOTIS/MIRNA_PUB_2024/"

RES <- read_tsv(list.files(path = dir, 
  pattern = "SEQUENCES_MERGED_DESEQ_RES.tsv", full.names = T)) 

RES.P <- RES %>% 
  filter( padj < 0.05  & abs(log2FoldChange) > 1)

RES.P %>% 
  count(CONTRAST) %>% arrange(CONTRAST)

DEGS <- RES.P %>% 
  filter(CONTRAST == "CONTRAST_C") %>%
  mutate(HPF = ifelse(sign(log2FoldChange) == -1, "24 hpf", "110 hpf")) %>%
  distinct(MajorRNA, HPF, CONTRAST)

DEGS %>% count(HPF, CONTRAST)

# JOIN LOW PH DEGS -----
# 76 MIRS DEGS IN OA

DEGS_D <- RES.P %>% 
  filter(CONTRAST == "CONTRAST_D") %>%
  mutate(HPF = ifelse(sign(log2FoldChange) == -1, "24 hpf", "110 hpf")) %>%
  distinct(MajorRNA, HPF, CONTRAST)

# DEGS %>% ggplot(aes(log2FoldChange, fill = HPF)) + geom_histogram()


# JOIN DEGS TO LONGER_DB ----

plotdf <- DEGS %>% 
  left_join(DB) %>% 
  drop_na(COG_name) %>%
  # distinct(MajorRNA, HPF, CONTRAST) %>% count(HPF, CONTRAST)
  distinct(MajorRNA, HPF, CONTRAST, COG_name) %>% 
  count(COG_name, HPF, CONTRAST, sort = T)


order_n <- plotdf %>% group_by(COG_name) %>% tally(n, sort = T) %>% pull(COG_name)

sum(DEGS_D$MajorRNA %in% DEGS$MajorRNA)

recode_to <- structure(c("pH 8.0","pH 7.6"), names = c("CONTRAST_C", "CONTRAST_D"))


# DEGS_D %>% 
#   # anti_join(DEGS, by = "MajorRNA") %>%
#   left_join(DB) %>% 
#   drop_na(COG_name) %>%
#   distinct(MajorRNA, HPF, CONTRAST, COG_name) %>% 
#   count(COG_name, HPF, CONTRAST, sort = T) %>%
#   rbind(plotdf) %>%
#   # group_by(COG_name, CONTRAST) %>% mutate(frac = n / sum(n)) %>%
#   mutate(COG_name = factor(COG_name, levels = rev(order_n))) %>%
#   mutate(HPF = factor(HPF, levels = c("24 hpf","110 hpf"))) %>%
#   mutate(CONTRAST = dplyr::recode_factor(CONTRAST, !!!recode_to)) %>%
#   ggplot(aes(y = COG_name, x = n, fill = CONTRAST)) + 
#   labs(y = "COG", x = "Number of microRNAs") +
#   facet_grid(~ HPF) +
#   geom_col(position = position_dodge()) +
#   scale_fill_grey() +
#   guides(fill=guide_legend(title = "", nrow = 1)) +
#   theme_bw(base_size = 12, base_family = "GillSans") +
#   theme(legend.position = "top", 
#     strip.background = element_rect(fill = 'grey89', color = 'white'),
#     axis.line.x = element_blank(),
#     axis.line.y = element_blank()) 

plotdf <- DEGS_D %>% 
  # anti_join(DEGS, by = "MajorRNA") %>%
  left_join(DB) %>% 
  drop_na(COG_name) %>%
  distinct(MajorRNA, HPF, CONTRAST, COG_name) %>% 
  count(COG_name, HPF, CONTRAST, sort = T) %>%
  rbind(plotdf) %>%
  group_by(COG_name, CONTRAST) %>% mutate(frac = n / sum(n)) %>%
  mutate(COG_name = factor(COG_name, levels = order_n)) %>%
  mutate(HPF = factor(HPF, levels = c("24 hpf","110 hpf"))) %>%
  mutate(CONTRAST = dplyr::recode_factor(CONTRAST, !!!recode_to))

data_text <- plotdf %>% group_by(COG_name, CONTRAST) %>% summarise(n = sum(n)) %>%
  mutate(n = paste0("(", n, ")"))

plotdf %>%
  ggplot() + 
  labs(y = "NOGS", x = "Enrichment ratio (miRNA target set)") +
  facet_grid(COG_name~ ., scales = "free_y", space = "free_y", switch = "y") + 
  geom_col(aes(y = CONTRAST, x = frac, fill = HPF))  +
  # scale_y_discrete(position = "right") +
  scale_fill_grey() +
  guides(fill=guide_legend(title = "", nrow = 1)) +
  theme_classic(base_size = 12, base_family = "GillSans") +
  theme(legend.position = "top", 
    strip.background = element_rect(fill = 'grey95', color = 'white'),
    axis.line.x = element_blank(),
    # axis.line.y = element_blank(),
    axis.text.y = element_text(size = 7, hjust = 1),
    strip.text.y.left = element_text(angle = 0, hjust = 1),
    strip.placement = "outside") +
  geom_text(data = data_text, 
    aes(y = CONTRAST, x = 1, group = CONTRAST, label = n), size = 2.5,
    hjust = -0.1, vjust = 0, 
    family = "GillSans", position = position_dodge(width = 1)) -> p

# p

ggsave(p, filename = 'NOGS.png', path = dir, width = 8, height = 8, device = png, dpi = 300)


# TOP GO ENRICHMENT ----
dir <- "~/Documents/MIRNA_HALIOTIS/MIRNA_PUB_2024/"

GODB <- read_tsv(file.path(dir, "Gene_ontologies_DB.tsv"))

GODB <- GODB %>% drop_na(GOs) %>% filter(GOs != "-") %>% 
  distinct(gene_id, GOs)

# QUERY DEGS PRIOR TO SPLIT
# QUERY ONLY MIRS CONSERVED IN Ctrl

QUERYDB <- RES.P %>% filter(CONTRAST == "CONTRAST_C") %>%
  mutate(HPF = ifelse(sign(log2FoldChange) == -1, "24 hpf", "110 hpf"))

# CAUTION: GENES ARE MULTI-TARGETS

# DB %>% 
#   distinct(gene_id, MajorRNA) %>%
#   count(MajorRNA, sort = T)

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
  
  
  allRes <- runtopGO(topGOdata, topNodes = Nodes, conservative = T)
  
  allRes <- data.frame(allRes, GROUP = x)
}

QUERYDB <- RES.P %>% 
  mutate(dir = ifelse(sign(log2FoldChange) == -1, "-logFC","logFC"), 
    CONTRAST = paste0(CONTRAST, dir))

CONTRAST2GO <- QUERYDB %>% distinct(CONTRAST) %>% pull()

OUT <- lapply(CONTRAST2GO, function(x) runTopGO(x, QUERYDB, Nodes = 50))

print(data <- dplyr::bind_rows(OUT, .id = "DE") %>% as_tibble())

hist(as.numeric(data$classicFisher))

# Next we want to analyse the distribution of the genes annotated to a GO term of interest. In an enrichment analysis one expects that the genes annotated to a signicantly enriched GO term have higher scores than the average gene' score of the gene universe

goID <- allRes[20, "GO.ID"]

showGroupDensity(topGOdata, goID, ranks = TRUE)

# df <- GroupDensityDF(topGOdata, goID, ranks = TRUE)

geneScore(topGOdata, use.names = TRUE)

# genesInTerm(topGOdata, goID)[[1]] 
#
data %>%
  mutate(classicKS = as.double(classicKS), elimKS = as.double(elimKS)) %>%
  ggplot(aes(x = as.factor(`Rank.in.classicFisher`), y = -log10(elimKS))) + 
  facet_grid(GROUP ~.) +
  geom_path(group = 1)
