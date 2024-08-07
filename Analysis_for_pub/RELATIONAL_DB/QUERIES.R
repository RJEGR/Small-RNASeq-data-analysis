
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

# ----
# Define biological-phentypic themes
# Calcification <---
# Development
# Growth
# RM (Respiratory Metabolism)
# Group transcriptome to these themes
# Idealy, using only miRNA:mRNA transcripts (~ 170 transcripts)
# Additionally, blast miRNA:mRNA transcripts to biomineralization-genes database

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

WGCNA <- read_rds(paste0(dir, "WGCNA_MIRS.rds"))

# or load SEQUENCES_MERGED_DESEQ_RES_WGCNA.rds

RES <- read_tsv(list.files(path = dir, 
  pattern = "SEQUENCES_MERGED_DESEQ_RES.tsv", full.names = T)) %>%
  left_join(WGCNA) %>%
  mutate(Contrast = ifelse(sign(log2FoldChange) == -1, sampleB, sampleA)) %>%
  mutate(Contrast = paste0(CONTRAST,"_",Contrast))

RES.P <- RES %>% filter( padj < 0.05  & abs(log2FoldChange) > 1)

# RES.P %>% count(Contrast) %>% arrange(Contrast)

DEGS <- RES.P %>% 
  filter(CONTRAST == "CONTRAST_C") %>%
  mutate(HPF = ifelse(sign(log2FoldChange) == -1, "110 hpf", "24 hpf")) %>%
  distinct(MajorRNA, HPF, CONTRAST)

DEGS %>% count(HPF, CONTRAST)

# JOIN LOW PH DEGS -----
# 76 MIRS DEGS IN OA

DEGS_D <- RES.P %>% 
  filter(CONTRAST == "CONTRAST_D") %>% 
  mutate(HPF = ifelse(sign(log2FoldChange) == -1, "110 hpf", "24 hpf")) %>%
  distinct(MajorRNA, HPF, CONTRAST)

# DEGS %>% ggplot(aes(log2FoldChange, fill = HPF)) + geom_histogram()


# JOIN DEGS TO LONGER_DB ----

plotdf <- DEGS %>% 
  left_join(DB) %>% 
  drop_na(COG_name) %>%
  # distinct(MajorRNA, HPF, CONTRAST) %>% count(HPF, CONTRAST)
  distinct(MajorRNA, HPF, CONTRAST, COG_name) %>% 
  count(COG_name, HPF, CONTRAST, sort = T)

sum(plotdf$n)

order_n <- plotdf %>% group_by(COG_name) %>% tally(n, sort = T) %>% pull(COG_name)

# sum(DEGS_D$MajorRNA %in% DEGS$MajorRNA)

recode_to <- structure(c("pH 8.0","pH 7.6"), names = c("CONTRAST_C", "CONTRAST_D"))


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
  mutate(CONTRAST = dplyr::recode_factor(CONTRAST, !!!rev(recode_to)))

data_text <- plotdf %>% 
  group_by(COG_name, CONTRAST) %>% 
  summarise(n = sum(n)) %>%
  mutate(n = paste0("(", n, ")")) %>%
  mutate(CONTRAST = factor(CONTRAST, levels = rev(recode_to))) %>%
  mutate(log2FoldChange = -10)


# Find brach of interactions -----

rbind(DEGS, DEGS_D) %>%
  left_join(DB) %>% 
  drop_na(COG_name) %>% 
  distinct(MajorRNA, COG_name) %>%
  mutate(n = 1) %>%
  pivot_wider(names_from = COG_name, values_from = n, values_fill = 0) %>%
  data.frame() -> m


rownames(m) <- m$MajorRNA

m$MajorRNA <- NULL

yhc <- stats::hclust(dist(t(m), method = "binary"), method="ward.D2")

xhc <- stats::hclust(dist(m, method = "binary"), method="ward.D2")


library(ggh4x)

rbind(DEGS, DEGS_D) %>%
  left_join(DB) %>% 
  drop_na(COG_name) %>% 
  distinct(MajorRNA, HPF, CONTRAST, COG_name) %>%
  # mutate(COG_name = gsub("[[:blank:]]|,", ".",COG_name)) %>%
  mutate(COG_name = factor(COG_name, levels = order_n)) %>%
  ggplot(aes(fill = 1, y = COG_name, x = MajorRNA)) +
  geom_tile(color = 'white', size = 0.7, width = 1) +
  # ggh4x::scale_y_dendrogram(hclust = yhc) +
  ggh4x::scale_x_dendrogram(hclust = xhc, label = NULL)

  
# scale_col <- c("#cd201f", "#FFFC00","#00b489","#31759b")


plotdf %>%
  ggplot() + 
  labs(y = "NOGS", x = "Enrichment ratio (miRNA target set)") +
  facet_grid(COG_name~ ., scales = "free_y", space = "free_y", switch = "y") +
  # facet_grid(COG_name~ HPF, scales = "free_y", space = "free_y", switch = "y") +
  geom_col(position = position_stack(reverse = T),
    aes(y = CONTRAST, x = frac, fill = HPF))  +
  geom_text(aes(y = CONTRAST, x = frac, label=CONTRAST), x = 0.01, hjust=0, size = 2.5, family = "GillSans") +
  # scale_y_discrete(position = "right") +
  # scale_fill_grey() +
  scale_fill_manual(values = c("grey89", "gray20")) +
  guides(fill=guide_legend(title = "", nrow = 1)) +
  theme_classic(base_size = 12, base_family = "GillSans") +
  theme(legend.position = "top", 
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    strip.background = element_rect(fill = 'grey95', color = 'white'),
    axis.line.x = element_blank(), axis.text.y = element_blank(),
    # axis.line.y = element_blank(),
    # axis.text.y = element_text(size = 7, hjust = 1),
    strip.text.y.left = element_text(angle = 0, hjust = 1),
    strip.placement = "outside") +
  geom_text(data = data_text, 
    aes(y = CONTRAST, x = 1, group = CONTRAST, label = n), size = 2.5,
    hjust = -0.1, vjust = 0, 
    family = "GillSans", position = position_dodge(width = 1)) -> p

# p

# ggsave(p, filename = 'NOGS.png', path = dir, width = 8, height = 8, device = png, dpi = 300)


# OR CONSIDERING LOG2FC -----

# (optional) 1st How to diffentiated the mirs?
# lets to upset first

query <- RES.P %>% 
  filter(CONTRAST %in% c("CONTRAST_C", "CONTRAST_D")) %>%
  mutate(CONTRAST = dplyr::recode_factor(CONTRAST, !!!rev(recode_to))) %>%
  group_by(MajorRNA) %>%
  summarise(across(CONTRAST, .fns = list), n = n()) %>% arrange(desc(n)) %>%
  filter(n == 1) %>% distinct(MajorRNA) %>% pull()


pdens <- RES.P %>% 
  filter(CONTRAST %in% c("CONTRAST_C", "CONTRAST_D")) %>%
  mutate(CONTRAST = dplyr::recode_factor(CONTRAST, !!!rev(recode_to))) %>%
  left_join(DB, by = "MajorRNA") %>% 
  drop_na(COG_name) %>%
  # if filter the intersected
  filter(MajorRNA %in% query) %>%
  # mutate(COG_name = ifelse(is.na(COG_name), "Unknown", COG_name)) %>%
  mutate(COG_name = factor(COG_name, levels = order_n)) %>%
  # !Note: here invert log2FoldChange for sort  24 and 110 hpf in the plot
  mutate(log2FoldChange = log2FoldChange*-1) %>%
  # drop_na(SMPID) %>%
  ggplot(aes(y = COG_name, x = log2FoldChange, fill = CONTRAST, color = CONTRAST)) + 
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.5) +
  # geom_point() +
  ggridges::geom_density_ridges(
    alpha = 0.5, 
    jittered_points = T,
    position = ggridges::position_points_jitter(width = 0.5, height = 0),
    point_shape = '|', point_size = 1, point_alpha = 1, alpha = 0.2) +
  scale_fill_manual(values = c("black", "gray")) +
  scale_color_manual(values = c("black", "gray")) +
  guides(fill=guide_legend(title = "", nrow = 1), color = "none") +
  annotate("text", x = -9, y = 20, label = "24 hpf", family = "GillSans") +
  annotate("text", x = 9, y = 20, label = "110 hpf", family = "GillSans") +
  annotate("text", x = 0, y = 21, label = "") 

pdens <- pdens +
  labs(y = "Nested Orthologous Gene Group (NOGS)") +
  theme_bw(base_size = 12, base_family = "GillSans") +
  theme(legend.position = "top",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.key.width = unit(0.2, "cm"),
    legend.key.height = unit(0.12, "cm"))

pdens<- pdens + 
  facet_grid(~ CONTRAST) +
  geom_text(data = data_text,
    aes(y = COG_name, x = -14, group = CONTRAST, label = n), 
    size = 2.5, hjust = -0.1, vjust = 0, 
    family = "GillSans", position = position_dodge(width = 1))

ggsave(pdens, filename = 'NOGS_dens_facet.png', path = dir, width = 5.6, height = 8, device = png, dpi = 300)


RES.P %>% 
  filter(CONTRAST %in% c("CONTRAST_A", "CONTRAST_B")) %>%
  mutate(CONTRAST = dplyr::recode_factor(CONTRAST, !!!c("CONTRAST_B" = "110 hpf"))) %>%
  left_join(DB, by = "MajorRNA") %>% 
  drop_na(COG_name) %>%
  # mutate(COG_name = ifelse(is.na(COG_name), "Unknown", COG_name)) %>%
  mutate(COG_name = factor(COG_name, levels = order_n)) %>%
  # drop_na(SMPID) %>%
  ggplot(aes(y = COG_name, x = log2FoldChange, fill = CONTRAST, color = CONTRAST)) + 
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.5) +
  geom_point() +
  facet_grid(~ CONTRAST) +
  scale_fill_manual(values = c("black", "gray")) +
  scale_color_manual(values = c("black", "gray")) +
  ggrepel::geom_text_repel(aes(label = MirGeneDB_ID), 
    size = 2.5, hjust = -0.1, vjust = 0, 
    family = "GillSans", position = position_dodge(width = 1), max.overlaps = 50)

# by wgcna module ----

RES.P %>% 
  filter(CONTRAST %in% c("CONTRAST_C", "CONTRAST_D")) %>%
  mutate(CONTRAST = dplyr::recode_factor(CONTRAST, !!!rev(recode_to))) %>%
  left_join(DB, by = "MajorRNA") %>% 
  drop_na(COG_name) %>%
  distinct(WGCNA, COG_name, MajorRNA) %>%
  count(WGCNA, COG_name) %>%
  mutate(COG_name = factor(COG_name, levels = order_n)) %>%
  ggplot() +
  geom_tile(aes(fill = n, y = COG_name, x = WGCNA),
    color = 'white', size = 0.7, width = 1) +
  theme(axis.text.x = element_text(angle = 45, size = 10, hjust = 1))

# TOP GO ENRICHMENT ----
dir <- "~/Documents/MIRNA_HALIOTIS/MIRNA_PUB_2024/"

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


data %>% ungroup() %>% count(f1, f2)

p <- data %>%
  mutate(classicKS = as.double(classicKS), 
    elimKS = as.double(elimKS)) %>%
  ggplot(aes(x = Top, y = parentTerm, fill = -log10(classicKS))) +
  geom_tile(color = 'white', linewidth = 0.2) +
  facet_grid(~ f1+f2, switch = "y", scales = "free_y") +
  theme_bw(base_family = "GillSans", base_size = 10) +
  # guides(fill = "none") +
  labs(y = "Biological process", x = "Top Enrichment") +
  scale_fill_viridis_c("Size", option = "inferno", direction = -1)

p

p <- data %>%  
  # filter(as.numeric(classicKS) < 0.05) %>%
  group_by(parentTerm, f1, f2) %>% 
  # count()
  summarise(Ratio = size/sum(size)) %>% 
  ggplot(aes(x = f2, y = parentTerm, fill = Ratio)) +
  geom_tile(color = 'white', linewidth = 0.2) +
  facet_grid(~ f1, switch = "y", scales = "free_y") +
  theme_bw(base_family = "GillSans", base_size = 10) +
  # guides(fill = "none") +
  labs(y = "Biological process", x = "") +
  scale_fill_viridis_c("Enrichment ratio", option = "inferno", direction = -1)

p +
  theme(legend.position = "top",
    # axis.text.y.left = element_blank(),
    axis.ticks.y.left = element_blank(),
    strip.background.x = element_rect(fill = 'grey89', color = 'white'),
    strip.text.y.left = element_text(angle = 0, size = 10, hjust = 1),
    strip.background.y = element_rect(fill = 'white', color = 'white'),
    panel.border = element_blank(),
    plot.title = element_text(hjust = 0),
    plot.caption = element_text(hjust = 0),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.y = element_text(angle = 0, size = 10),
    axis.text.x = element_text(angle = 90, size = 10)) 


# Next we want to analyse the distribution of the genes annotated to a GO term of interest. In an enrichment analysis one expects that the genes annotated to a signicantly enriched GO term have higher scores than the average gene' score of the gene universe

# goID <- allRes[20, "GO.ID"]

# showGroupDensity(topGOdata, goID, ranks = TRUE)

# df <- GroupDensityDF(topGOdata, goID, ranks = TRUE)

# geneScore(topGOdata, use.names = TRUE)

# genesInTerm(topGOdata, goID)[[1]] 
#
data %>%
  mutate(classicKS = as.double(classicKS), 
         elimKS = as.double(elimKS)) %>%
  ggplot(aes(x = as.factor(`Rank.in.classicFisher`), y = -log(elimKS))) + 
  geom_boxplot()
  # facet_grid(GROUP ~.) +
  # geom_path(group = 1)
