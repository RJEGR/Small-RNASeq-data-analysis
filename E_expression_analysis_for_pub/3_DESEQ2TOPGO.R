
# RUN GO ANALYSIS USING:
# 0) PREPARE UPSET 
# 1) SPLIT DOWN FROM UP SIGNIFICANT (PADJ < 0.05) VALUES
# 2) RUN INDEPENDENTLY THE TOP GO ANALYSIS
# 3) PLOT SOME STUFFS


# NOTE:

# POSITIVE LFC == UP EXPRESSED IN EXPERIMENTAL (OR sampleB)
# NEGATIVE LFC == UP EXPRESSED IN CONTROL (OR sampleA)
# view baseMeanA and baseMeanB to contrast the expression values


rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

wd <- "~/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out/"

RES.P <- read_tsv(paste0(wd, "DESEQ_RES.tsv")) %>% filter( padj < 0.05 & abs(log2FoldChange) > 1)

RES.P %>% dplyr::count(CONTRAST)

.UPSETDF <- read_rds(paste0(wd, "UPSETDF.rds"))

EXCLUSIVE_MIRS <- .UPSETDF[[1]] %>% ungroup()

INTERSECTED_MIRS <- .UPSETDF[[2]]  %>% ungroup()

EXCLUSIVE_MIRS %>% dplyr::count(CONTRAST, SIGN, sort = T) 

INTERSECTED_MIRS %>% dplyr::count(CONTRAST, SIGN, sort = T)

wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/FUNCTIONAL_MIR_ANNOT/"

print(.SRNA2GO <- read_tsv(paste0(wd, "SRNA_REGULATORY_FUNCTION_DB.tsv")))

url <- "https://raw.githubusercontent.com/RJEGR/Cancer_sete_T_assembly/main/functions.R"

source(url)

# Mature are you functional strand used in DE

paste_go <- function(x) { 
  x <- x[!is.na(x)] 
  x <- unique(sort(x))
  x <- paste(x, sep = ';', collapse = ';')
}


SRNA2GO <- .SRNA2GO %>%
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

RES.P <- SRNA2GO %>% 
  dplyr::rename("Name" = "query") %>%
  dplyr::select(Name, GO.ID) %>%
  right_join(RES.P)

RES.P <- .SRNA2GO %>%
  filter(predicted == "BOTH") %>%
  mutate(query = strsplit(query, ";")) %>%
  unnest(query) %>% 
  distinct(query, gene_id) %>%
  group_by(query) %>%
  summarise(
    across(gene_id, .fns = paste_go), 
    .groups = "drop_last") %>%
  separate(query, into = c("query", "arm"), sep = "[.]") %>%
  filter(arm == "mature") %>%
  dplyr::rename("Name" = "query") %>%
  dplyr::select(Name, gene_id) %>%
  right_join(RES.P)
  
write_tsv(RES.P, paste0(wd, "DESEQ_RES_P.tsv"))

# RES.P %>% 

SRNA2GO <- split(strsplit(SRNA2GO$GO.ID, ";") , SRNA2GO$query)

SRNA2GO <- lapply(SRNA2GO, unlist)

# INTERSECTED ====

# 38+46=84 mirs shared under low and control pH at 24 HPF and 110 HPF, respectively

INTERSECTED_MIRS %>% dplyr::count(CONTRAST, SIGN, sort = T)

INTERSECTED_MIRS %>% distinct(Name)
  
str(WHICH_SIGN <- INTERSECTED_MIRS %>% ungroup() %>% distinct(SIGN) %>% pull())

DF <- list()

for(j in 1:length(WHICH_SIGN)) {
  
  i <- j
  
  cat("\nRunning ",WHICH_SIGN[i], "\n")
  
  query.names <- INTERSECTED_MIRS %>% ungroup() %>% filter(SIGN %in% WHICH_SIGN[i]) %>% distinct(Name) %>% pull(Name)
  
  query.p <- c(rep(0.05, length(query.names)))
  
  allRes <- GOenrichment(query.p, query.names, SRNA2GO, Nodes = 20)
  
  allRes <- allRes %>% mutate(SIGN = WHICH_SIGN[i]) %>% as_tibble()
  
  DF[[i]] <- allRes
  
}

INTERSECTED_MIRS_ALLRes <- do.call(rbind, DF) %>% mutate(CONTRAST = "Both pH")

view(INTERSECTED_MIRS_ALLRes)


INTERSECTED_MIRS %>% distinct(Name) %>%
  left_join(RES.P) %>% group_by(Name) %>% view()


INTERSECTED_MIRS_ALLRes %>%
  group_by(SIGN) %>% arrange(Term) %>%
  # mutate(Term = factor(Term, levels= unique(Term))) %>%
  mutate(Term = fct_reorder(Term, Annotated)) %>%
  ggplot(aes(y = Term, x = Annotated, fill = p.adj.ks)) + # , fill = p.adj.ks fill = -log10(p.adj.ks))
  geom_col() +
  # geom_text(aes(label = Annotated), color = 'black', family = "GillSans", size = 3, angle = 0, hjust = -0.5) +
  theme_bw(base_family = "GillSans", base_size = 12) +
  ggh4x::facet_nested( SIGN~., nest_line = F, scales = "free", space = "free") +
  theme_bw(base_family = "GillSans", base_size = 12) +
  ggsci::scale_fill_material("blue-grey") +
  labs(x = "miRs anotados", y = "") +
  theme(strip.background = element_rect(fill = 'grey89', color = 'white'),
    # legend.position = "top",
    panel.border = element_blank(),
    plot.title = element_text(hjust = 0),
    plot.caption = element_text(hjust = 0),
    panel.grid = element_blank(),
    # axis.text.y.right = element_text(angle = 0, hjust = 1, vjust = 0, size = 2.5),
    axis.text.y = element_text(angle = 0, size = 7),
    axis.text.x = element_text(angle = 0, size = 10)) -> p

ggsave(p, filename = 'DESEQ2TOPGO_INTERSECTED.png', path = wd, width = 7, height = 4, device = png, dpi = 300)

# INTERSECTED_MIRS_ALLRes
# SEMANTIC ====

library(GOSemSim)

# hsGO <- godata('org.Hs.eg.db', ont="BP")

# rbind(TOPGO_INTERSECTED, TOPGO_EXCLUSIVE)


str(GO.ID <- INTERSECTED_MIRS_ALLRes %>% distinct(GO.ID) %>% pull())

# hsGO <- GOSemSim::godata('org.Hs.eg.db', ont="BP")

hsGO <- read_rds(paste0(wd, '/hsGO_BP.rds'))

termSim <- GOSemSim::termSim(GO.ID, GO.ID,  semData = hsGO, method = "Wang")

# heatmap(termSim)

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
  cutree(., 7)

# Sanity check

identical(cmds$GO.ID, names(cutree))

head(cmds <- cbind(cmds, cutree) %>% mutate(cutree = as.factor(cutree)))

INTERSECTED_MIRS_ALLRes <- INTERSECTED_MIRS_ALLRes %>% left_join(cmds)


# EXCLUSIVE ====

# CONTRAST SIGN           n
# 1 pH 7.6   B) 110 HPF    18
# 2 pH 8.0   A) 24 HPF      9
# 3 pH 8.0   B) 110 HPF     1

WHICH_SIGN <- "A) 24 HPF"

str(CONTRAST <- EXCLUSIVE_MIRS %>% ungroup() %>% filter(SIGN == WHICH_SIGN) %>% distinct(CONTRAST) %>% pull())

CONTRAST <- as.character(CONTRAST)

EXCLUSIVE_MIRS %>% unnest(CONTRAST_DE) %>% filter(SIGN == WHICH_SIGN) %>% dplyr::count(SIGN, CONTRAST)

str(query.names <- EXCLUSIVE_MIRS %>% unnest(CONTRAST_DE) %>% 
    filter(SIGN == WHICH_SIGN) %>% distinct(Name) %>% pull())

query.p <- c(rep(0.05, length(query.names)))

EXCLUSIVE_MIRS_ALLRes_24HPF <- GOenrichment(query.p, query.names, SRNA2GO, Nodes = 3)

EXCLUSIVE_MIRS_ALLRes_24HPF <- EXCLUSIVE_MIRS_ALLRes_24HPF %>% mutate(SIGN = WHICH_SIGN, CONTRAST = CONTRAST) %>% as_tibble()

view(EXCLUSIVE_MIRS_ALLRes_24HPF)

# RES.P %>%
#   filter(Name %in% query.names) %>%
#   left_join(distinct(.SRNA2GO, gene_id, query)) %>%
#   view()

.SRNA2GO %>% 
  mutate(query = strsplit(query, ";")) %>%
  unnest(query) %>%
  separate(query, into = c("query", "arm"), sep = "[.]") %>% 
  filter(arm == "mature") %>%
  filter(query %in% query.names & predicted == "BOTH") %>% 
  group_by(query)

# str(CONTRAST <- EXCLUSIVE_MIRS %>% ungroup() %>% distinct(CONTRAST) %>% pull())

# SANITY CHECK: IF SINGLE MIR NOT GO ENRICHMENT IS POSSIBLE  

drop.names <- EXCLUSIVE_MIRS %>%
  group_by(CONTRAST, SIGN) %>%
  summarise(across(Name, .fns = list), n = n()) %>%
  filter(n == 1) %>%
  unnest(Name) %>% 
  pull(Name)


WHICH_SIGN <- "B) 110 HPF"

str(CONTRAST <- EXCLUSIVE_MIRS %>% ungroup() %>% filter(SIGN == WHICH_SIGN) %>% distinct(CONTRAST) %>% pull())

CONTRAST <- as.character(CONTRAST)

EXCLUSIVE_MIRS %>% unnest(CONTRAST_DE) %>% filter(SIGN == WHICH_SIGN) %>% dplyr::count(SIGN, CONTRAST)

str(query.names <- EXCLUSIVE_MIRS %>% unnest(CONTRAST_DE) %>% 
    filter(SIGN == WHICH_SIGN & Name != drop.names) %>% distinct(Name) %>% pull())

query.p <- c(rep(0.05, length(query.names)))

EXCLUSIVE_MIRS_ALLRes_110HPF <- GOenrichment(query.p, query.names, SRNA2GO, Nodes = 10)

EXCLUSIVE_MIRS_ALLRes_110HPF <- EXCLUSIVE_MIRS_ALLRes_110HPF %>% mutate(SIGN = WHICH_SIGN, CONTRAST = CONTRAST[1]) %>% as_tibble()

view(rbind(INTERSECTED_MIRS_ALLRes, EXCLUSIVE_MIRS_ALLRes_24HPF, EXCLUSIVE_MIRS_ALLRes_110HPF))

INTERSECTED_MIRS_ALLRes


print(EXCLUSIVE_MIRS_ALLRes_110HPF <- do.call(rbind, DF) %>% mutate(SIGN = WHICH_SIGN))

# DOWN

WHICH_SIGN <- "B) 110 HPF"

str(CONTRAST <- EXCLUSIVE_MIRS %>% ungroup() %>% filter(SIGN == WHICH_SIGN) %>% distinct(CONTRAST) %>% pull())

DF <- list()

for(j in 1:length(CONTRAST)) {
  
  i <- j
  
  cat("\nRunning ",CONTRAST[i], "\n")
  
  query.names <- EXCLUSIVE_MIRS %>% ungroup() %>% filter(SIGN == WHICH_SIGN) %>% filter(CONTRAST %in% CONTRAST[i]) %>% pull(Name)
  
  query.p <- c(rep(0.05, length(query.names)))
  
  allRes <- GOenrichment(query.p, query.names, SRNA2GO, Nodes = 20)
  
  allRes <- allRes %>% mutate(CONTRAST = CONTRAST[i]) %>% as_tibble()
  
  DF[[i]] <- allRes
  
}

EXCLUSIVE_DOWN_MIRS <- do.call(rbind, DF) %>% mutate(SIGN = WHICH_SIGN)

rbind(EXCLUSIVE_DOWN_MIRS, EXCLUSIVE_UP_MIRS) -> TOPGO_EXCLUSIVE

p <- TOPGO_EXCLUSIVE %>%
  # filter(p.adj.ks < 0.05) %>%
  arrange(Term) %>%
  mutate(Term = factor(Term, levels= unique(Term))) %>%
  separate(CONTRAST, into = c("WRAP","CONTRAST"), sep = "[|]") %>%
  ggplot(aes(y = Term, x = WRAP)) + # , fill = p.adj.ks fill = -log10(p.adj.ks))
  # geom_col() +
  # geom_point(aes(size = Annotated)) +
  geom_tile(aes(fill = Annotated), color = 'white', size = 0.7, width = 1) +
  geom_text(aes(label = Annotated), color = 'black', family = "GillSans", size = 1.5, angle = 90) +
  theme_bw(base_family = "GillSans", base_size = 12) +
  ggh4x::facet_nested( ~ CONTRAST+SIGN, nest_line = F, scales = "free", space = "free") +
  theme_bw(base_family = "GillSans", base_size = 12) +
  ggsci::scale_fill_material("blue-grey") +
  labs(x = "", y = "")

p <- p + theme(legend.position = "none",
  strip.background = element_rect(fill = 'grey89', color = 'white'),
  panel.border = element_blank(),
  plot.title = element_text(hjust = 0),
  plot.caption = element_text(hjust = 0),
  panel.grid.minor.y = element_blank(),
  panel.grid.major.y = element_blank(),
  panel.grid.minor.x = element_blank(),
  panel.grid.major.x = element_blank(),
  # strip.background.y = element_blank(),
  axis.text.y = element_text(angle = 15, size = 5),
  axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10))

ggsave(p, filename = 'DESEQ2TOPGO_EXCLUSIVE.png', path = wd, width = 3.5, height = 10, device = png, dpi = 300)

# or plot by bars:
# split enrichment annotated by: if 

TOPGO_EXCLUSIVE %>% distinct(Term)

p <- TOPGO_EXCLUSIVE %>%
  # separate(CONTRAST, into = c("WRAP","CONTRAST"), sep = "[|]") %>%
  mutate(Annotated = ifelse(SIGN == "Down", -Annotated, Annotated)) %>%
  mutate(SIGN = factor(SIGN, levels = c("Up","Down"))) %>%
  mutate(ordering = -as.numeric(SIGN) + Annotated, Term = fct_reorder(Term, ordering, .desc = T)) %>%
  ggplot(aes(y = Term, x = Annotated, fill = SIGN, color = SIGN)) + 
  # geom_col() +
  geom_segment(aes(x = Annotated, xend = 0, yend = Term)) +
  ggh4x::facet_nested(CONTRAST+WRAP+SIGN ~ ., nest_line = F, scales = "free") +
  ggsci::scale_fill_aaas() +
  ggsci::scale_color_aaas() +
  theme_bw(base_family = "GillSans", base_size = 12) +
  labs(x = "", y = "") +
  theme(legend.position = "none",
    strip.background = element_rect(fill = 'grey89', color = 'white'),
    panel.border = element_blank(),
    plot.title = element_text(hjust = 0),
    plot.caption = element_text(hjust = 0),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_line(linetype = "dashed", linewidth = 0.5),
    axis.text.y = element_text(angle = 0, size = 5),
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10))

ggsave(p, filename = 'DESEQ2TOPGO_EXCLUSIVE_2.png', path = wd, width = 3.5, height = 10, device = png, dpi = 300)

TOPGO_EXCLUSIVE %>%
  mutate(Annotated = ifelse(SIGN == "A) 24 HPF B)", -Annotated, Annotated)) %>%
  mutate(SIGN = factor(SIGN, levels = c("A) 24 HPF","B) 110 HPF"))) %>%
  mutate(ordering = -as.numeric(SIGN) + Annotated, Term = fct_reorder(Term, ordering, .desc = T)) %>%
  ggplot(aes(y = Term, x = Annotated, fill = SIGN, color = SIGN)) + 
  # geom_col() +
  geom_segment(aes(x = Annotated, xend = 0, yend = Term)) +
  ggh4x::facet_nested(CONTRAST+SIGN ~ ., nest_line = F, scales = "free") +
  ggsci::scale_fill_aaas() +
  ggsci::scale_color_aaas() +
  theme_bw(base_family = "GillSans", base_size = 12) +
  labs(x = "", y = "") +
  theme(legend.position = "none",
    strip.background = element_rect(fill = 'grey89', color = 'white'),
    panel.border = element_blank(),
    plot.title = element_text(hjust = 0),
    plot.caption = element_text(hjust = 0),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_line(linetype = "dashed", linewidth = 0.5),
    axis.text.y = element_text(angle = 0, size = 5),
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10))


out <- list(TOPGO_EXCLUSIVE, TOPGO_INTERSECTED) 

names(out) <- c("EXCLUSIVE","INTERSECTED")

write_rds(out, file = paste0(wd, "DESEQ2TOPGO.rds"))


# QUIT ===
