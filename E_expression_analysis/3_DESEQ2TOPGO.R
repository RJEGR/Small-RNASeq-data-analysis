
# RUN GO ANALYSIS USING:
# 0) PREPARE UPSET 
# 1) SPLIT DOWN FROM UP SIGNIFICANT (PADJ < 0.05) VALUES
# 2) RUN INDEPENDENTLY THE TOP GO ANALYSIS

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

wd <- "~/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out/"

RES.P <- read_tsv(paste0(wd, "DESEQ_RES.tsv")) %>% filter( padj < 0.05)

.UPSETDF <- read_rds(paste0(wd, "UPSETDF.rds"))

EXCLUSIVE_MIRS <- .UPSETDF[[1]]

INTERSECTED_MIRS <- .UPSETDF[[2]]

wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/FUNCTIONAL_MIR_ANNOT/"

print(SRNA2GO <- read_tsv(paste0(wd, "SRNA_REGULATORY_FUNCTION_DB.tsv")))

url <- "https://raw.githubusercontent.com/RJEGR/Cancer_sete_T_assembly/main/functions.R"

source(url)

# Mature are you functional strand used in DE

paste_go <- function(x) { 
  x <- x[!is.na(x)] 
  x <- unique(sort(x))
  x <- paste(x, sep = ';', collapse = ';')
}


SRNA2GO <- SRNA2GO %>%
  filter(predicted == "BOTH") %>%
  mutate(query = strsplit(query, ";")) %>%
  unnest(query) %>% 
  distinct(query, GO.ID) %>%
  group_by(query) %>%
  summarise(
    across(GO.ID, .fns = paste_go), 
    .groups = "drop_last")

SRNA2GO <- SRNA2GO %>% 
  separate(query, into = c("query", "arm"), sep = "[.]") %>%
  filter(arm == "mature")


SRNA2GO <- split(strsplit(SRNA2GO$GO.ID, ";") , SRNA2GO$query)

SRNA2GO <- lapply(SRNA2GO, unlist)

# EXCLUSIVE ====

str(CONTRAST <- EXCLUSIVE_MIRS %>% ungroup() %>% filter(SIGN == "Up") %>% distinct(CONTRAST) %>% pull())

# SANITY CHECK IF SINGLE MIR NOT GO ENRICHMENT IS POSSIBLE  

EXCLUSIVE_MIRS %>%
  group_by(CONTRAST, SIGN) %>%
  summarise(across(Name, .fns = list), n = n())
  

# UP

DF <- list()

for(j in 1:length(CONTRAST)) {
  
  i <- j
  
  cat("\nRunning ",CONTRAST[i], "\n")
  
  # DF2GO <- RES %>% filter(CONTRAST %in% CONTRAST[i])
  
  # DF2GO <- DF2GO %>% filter(padj < 0.05) # filter(log2FoldChange < 0 )
  
  # query.p <- DF2GO %>% pull(padj)
  
  # query.names <- DF2GO %>% pull(Name)
  
  query.names <- EXCLUSIVE_MIRS %>% ungroup() %>% filter(SIGN == "Up") %>% filter(CONTRAST %in% CONTRAST[i]) %>% pull(Name)
  
  query.p <- c(rep(0.05, length(query.names)))
  
  allRes <- GOenrichment(query.p, query.names, SRNA2GO, Nodes = 20)
  
  allRes <- allRes %>% mutate(CONTRAST = CONTRAST[i]) %>% as_tibble()
  
  DF[[i]] <- allRes
  
}

EXCLUSIVE_UP_MIRS <- do.call(rbind, DF) %>% mutate(SIGN = "Up")

# DOWN

str(CONTRAST <- EXCLUSIVE_MIRS %>% ungroup() %>% filter(SIGN == "Down") %>% distinct(CONTRAST) %>% pull())

DF <- list()

for(j in 1:length(CONTRAST)) {
  
  i <- j
  
  cat("\nRunning ",CONTRAST[i], "\n")

  query.names <- EXCLUSIVE_MIRS %>% ungroup() %>% filter(SIGN == "Down") %>% filter(CONTRAST %in% CONTRAST[i]) %>% pull(Name)
  
  query.p <- c(rep(0.05, length(query.names)))
  
  allRes <- GOenrichment(query.p, query.names, SRNA2GO, Nodes = 20)
  
  allRes <- allRes %>% mutate(CONTRAST = CONTRAST[i]) %>% as_tibble()
  
  DF[[i]] <- allRes
  
}

EXCLUSIVE_DOWN_MIRS <- do.call(rbind, DF) %>% mutate(SIGN = "Down")


# Similarity grouping ====

library(GOSemSim)

# hsGO <- godata('org.Hs.eg.db', ont="BP")


GO.ID <- rbind(EXCLUSIVE_DOWN_MIRS, EXCLUSIVE_UP_MIRS) %>% distinct(GO.ID) %>% pull()

hsGO <- read_rds(paste0(wd, '/hsGO_BP.rds'))

termSim <- GOSemSim::termSim(GO.ID, GO.ID,  semData = hsGO, method = "Wang")

cmds <- termSim %>% 
  # use distance metric
  dist(method = "euclidean") %>%
  # compute cMDS
  cmdscale() %>%
  data.frame() %>%
  rownames_to_column(var= 'GO.ID')

cutree <- termSim %>% 
  # use distance metric
  dist(method = "euclidean") %>%
  # compute cMDS
  hclust() %>%
  cutree(., 5)

head(cmds <- cbind(cmds, cutree) %>% mutate(cutree = as.factor(cutree)))

rbind(EXCLUSIVE_DOWN_MIRS, EXCLUSIVE_UP_MIRS) %>%
  left_join(cmds) -> TOPGO_EXCLUSIVE


p <- TOPGO_EXCLUSIVE %>%
  arrange(desc(Term)) %>%
  mutate(Term = factor(Term, levels= unique(Term))) %>%
  separate(CONTRAST, into = c("WRAP","CONTRAST"), sep = "[|]") %>%
  ggplot(aes(y = Term, x = WRAP)) + # , fill = p.adj.ks fill = -log10(p.adj.ks))
  # geom_col() +
  # geom_point(aes(size = Annotated)) +
  geom_tile(aes(fill = Annotated), color = 'white', size = 0.7, width = 1) +
  geom_text(aes(label = Annotated), color = 'black', family = "GillSans", size = 1.5) +
  theme_bw(base_family = "GillSans", base_size = 12) +
  ggh4x::facet_nested( ~ CONTRAST+SIGN, nest_line = F, scales = "free", space = "free") +
  theme_bw(base_family = "GillSans", base_size = 12) +
  ggsci::scale_fill_material() +
  labs(x = "", y = "")

p <- p + theme(legend.position = "top",
  strip.background = element_rect(fill = 'grey89', color = 'white'),
  panel.border = element_blank(),
  plot.title = element_text(hjust = 0),
  plot.caption = element_text(hjust = 0),
  panel.grid.minor.y = element_blank(),
  panel.grid.major.y = element_blank(),
  panel.grid.minor.x = element_blank(),
  panel.grid.major.x = element_blank(),
  # strip.background.y = element_blank(),
  axis.text.y = element_text(angle = 0, size = 5),
  axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10))

ggsave(p, filename = 'DESEQ2TOPGO.png', path = wd, width = 4, height = 10, device = png, dpi = 300)


# INTERSECTED ====

# DOWN INTERSECTED:

str(CONTRAST <- INTERSECTED_MIRS %>% ungroup() %>% filter(SIGN == "Down") %>% distinct(CONTRAST) %>% pull())

DF <- list()

for(j in 1:length(CONTRAST)) {
  
  i <- j
  
  cat("\nRunning ",CONTRAST[i], "\n")
  
  query.names <- INTERSECTED_MIRS %>% ungroup() %>% filter(SIGN == "Down") %>% filter(CONTRAST %in% CONTRAST[i]) %>% pull(Name)
  
  query.p <- c(rep(0.05, length(query.names)))
  
  allRes <- GOenrichment(query.p, query.names, SRNA2GO, Nodes = 20)
  
  allRes <- allRes %>% mutate(CONTRAST = CONTRAST[i]) %>% as_tibble()
  
  DF[[i]] <- allRes
  
}

INTERSECTED_DOWN_MIRS <- do.call(rbind, DF) %>% mutate(SIGN = "Down")


# UP INTERSECTED:

str(CONTRAST <- INTERSECTED_MIRS %>% ungroup() %>% filter(SIGN == "Up") %>% distinct(CONTRAST) %>% pull())

DF <- list()

for(j in 1:length(CONTRAST)) {
  
  i <- j
  
  cat("\nRunning ",CONTRAST[i], "\n")
  
  query.names <- INTERSECTED_MIRS %>% ungroup() %>% filter(SIGN == "Up") %>% filter(CONTRAST %in% CONTRAST[i]) %>% pull(Name)
  
  query.p <- c(rep(0.05, length(query.names)))
  
  allRes <- GOenrichment(query.p, query.names, SRNA2GO, Nodes = 20)
  
  allRes <- allRes %>% mutate(CONTRAST = CONTRAST[i]) %>% as_tibble()
  
  DF[[i]] <- allRes
  
}

INTERSECTED_UP_MIRS <- do.call(rbind, DF) %>% mutate(SIGN = "Up")

rbind(INTERSECTED_DOWN_MIRS, INTERSECTED_UP_MIRS) %>%
  left_join(cmds) -> TOPGO_INTERSECTED


p <- TOPGO_INTERSECTED %>%
  arrange(desc(Term)) %>%
  mutate(Term = factor(Term, levels= unique(Term))) %>%
  separate(CONTRAST, into = c("WRAP","CONTRAST"), sep = "[|]") %>%
  ggplot(aes(y = Term, x = WRAP)) + # , fill = p.adj.ks fill = -log10(p.adj.ks))
  # geom_col() +
  # geom_point(aes(size = Annotated)) +
  geom_tile(aes(fill = Annotated), color = 'white', size = 0.7, width = 1) +
  geom_text(aes(label = Annotated), color = 'black', family = "GillSans", size = 1.5) +
  theme_bw(base_family = "GillSans", base_size = 12) +
  ggh4x::facet_nested( ~ CONTRAST+SIGN, nest_line = F, scales = "free", space = "free") +
  theme_bw(base_family = "GillSans", base_size = 12) +
  ggsci::scale_fill_material() +
  labs(x = "", y = "")

p <- p + theme(legend.position = "top",
  strip.background = element_rect(fill = 'grey89', color = 'white'),
  panel.border = element_blank(),
  plot.title = element_text(hjust = 0),
  plot.caption = element_text(hjust = 0),
  panel.grid.minor.y = element_blank(),
  panel.grid.major.y = element_blank(),
  panel.grid.minor.x = element_blank(),
  panel.grid.major.x = element_blank(),
  # strip.background.y = element_blank(),
  axis.text.y = element_text(angle = 0, size = 5),
  axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10))

ggsave(p, filename = 'DESEQ2TOPGO_EXCLUSIVE.png', path = wd, width = 4, height = 10, device = png, dpi = 300)


# AS DIMENSIONAL SCATTERPLOT


str(GO.ID <- rbind(INTERSECTED_DOWN_MIRS, INTERSECTED_UP_MIRS) %>% distinct(GO.ID) %>% pull())

hsGO <- read_rds(paste0(wd, '/hsGO_BP.rds'))

termSim <- GOSemSim::termSim(GO.ID, GO.ID,  semData = hsGO, method = "Wang")

cmds <- termSim %>% 
  # use distance metric
  dist(method = "euclidean") %>%
  # compute cMDS
  cmdscale() %>%
  data.frame() %>%
  rownames_to_column(var= 'GO.ID')

cutree <- termSim %>% 
  # use distance metric
  dist(method = "euclidean") %>%
  # compute cMDS
  hclust() %>%
  cutree(., 5)

head(cmds <- cbind(cmds, cutree) %>% mutate(cutree = as.factor(cutree)))

TOPGO_INTERSECTED %>%
  ggplot(aes(X1,X2, color = cutree)) +
  geom_point(aes(alpha = p.adj.ks), size = 4) +
  ggrepel::geom_text_repel(aes(label = Term),
    max.overlaps = 20, color = 'black',
    family = "GillSans") +
  theme_classic(base_size = 16, base_family = "GillSans") +
  labs(x = 'Dimension 1', y = 'Dimension 2') 


# ENRICHMENT ANALYSIS: =====

.bioc_packages <- c('biomaRt','GOSemSim',"org.Hs.eg.db") # DESeq2'

sapply(c(.bioc_packages), require, character.only = TRUE)

# PLOT ====
# POSITIVE LFC == UP EXPRESSED IN EXPERIMENTAL
# NEGATIVE LFC == UP EXPRESSED IN CONTROL

# recode_to <- paste(LETTERS[1:4], ")", sep = "")

recode_to <- c("24 HPF: Ctrl|pH", "110 HPF: Ctrl|pH", "Ctrl: 24|110 HPF", "pH: 24|110")
recode_to <- structure(recode_to, names = CONTRAST)

# res.p <- prep_DE_data(allRes, alpha = 0.05, lfcThreshold = 2) %>% drop_na(cc)

sigfc <- "Sign and FC";pv <- "Sign";fc <- "FC"

colors_fc <- c("red2",  "#4169E1", "forestgreen", "grey30")

colors_fc <- structure(colors_fc, names = c(sigfc, pv, fc, "NS"))

allRes$Color <- 'NS'
allRes[which(abs(allRes$log2FoldChange) > 2), 'Color'] <- fc
allRes[which(abs(allRes$padj) <= 0.05), 'Color'] <- pv
allRes[which(allRes$padj <= 0.05 & abs(allRes$log2FoldChange) > 2), 'Color'] <- sigfc

colors_fc <- structure(colors_fc, names = c(sigfc, pv, fc, "NS"))

p <- allRes  %>%
  dplyr::mutate(CONTRAST = dplyr::recode_factor(CONTRAST, !!!recode_to))
ggplot(aes(y = -log10(pvalue), x = log2FoldChange, color = Color)) +
  geom_point()  

p <- p +
  scale_color_manual(name = "", values = colors_fc) +
  labs(x= expression(Log[2] ~ "Fold Change"), y = expression(-Log[10] ~ "P")) +
  theme_bw(base_family = "GillSans", base_size = 18) +
  theme(legend.position = "top")  +
  geom_abline(slope = 0, intercept = -log10(0.05), linetype="dashed", alpha=0.5) +
  geom_vline(xintercept = 0, linetype="dashed", alpha=0.5) 

p + facet_grid(~ CONTRAST)

allRes %>% 
  mutate(g = sign(log2FoldChange)) %>% 
  dplyr::count(Color, g, CONTRAST) %>%
  # mutate(g = ifelse(g == "1", "CON_CANCER", "SIN_CANCER")) %>%
  filter(Color != 'NS') %>%
  group_by(g, CONTRAST) %>%
  mutate(pct = n / sum(n)) %>%
  ggplot() +
  facet_grid(g ~ .) +
  geom_col(aes(x = CONTRAST, y = pct, fill = Color), width = 0.5) +
  scale_fill_manual(name = "", values = colors_fc[-4]) +
  theme_bw(base_family = "GillSans", base_size = 18) +
  labs(y = '% Transcripts', x = '') +
  theme(legend.position = 'top') + coord_flip()

