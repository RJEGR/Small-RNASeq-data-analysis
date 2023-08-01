
# RUN GO ANALYSIS USING:
# 0) PREPARE UPSET 
# 1) SPLIT DOWN FROM UP SIGNIFICANT (PADJ < 0.05) VALUES
# 2) RUN INDEPENDENTLY THE TOP GO ANALYSIS
# 3) PLOT SOME STUFFS

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

wd <- "~/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out/"

RES.P <- read_tsv(paste0(wd, "DESEQ_RES_P.tsv")) %>% filter( padj < 0.05)

.UPSETDF <- read_rds(paste0(wd, "UPSETDF.rds"))

EXCLUSIVE_MIRS <- .UPSETDF[[1]] %>% ungroup()

INTERSECTED_MIRS <- .UPSETDF[[2]]  %>% ungroup()


EXCLUSIVE_MIRS %>% distinct(Family)

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

SRNA2GO <- split(strsplit(SRNA2GO$GO.ID, ";") , SRNA2GO$query)

SRNA2GO <- lapply(SRNA2GO, unlist)

# EXCLUSIVE ====

str(CONTRAST <- EXCLUSIVE_MIRS %>% ungroup() %>% filter(SIGN == "Up") %>% distinct(CONTRAST) %>% pull())

# SANITY CHECK: IF SINGLE MIR NOT GO ENRICHMENT IS POSSIBLE  

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

print(EXCLUSIVE_UP_MIRS <- do.call(rbind, DF) %>% mutate(SIGN = "Up"))

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

rbind(INTERSECTED_DOWN_MIRS, INTERSECTED_UP_MIRS) -> TOPGO_INTERSECTED



p <- TOPGO_INTERSECTED %>%
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

ggsave(p, filename = 'DESEQ2TOPGO_INTERSECTED.png', path = wd, width = 3.5, height = 10, device = png, dpi = 300)


out <- list(TOPGO_EXCLUSIVE, TOPGO_INTERSECTED) 

names(out) <- c("EXCLUSIVE","INTERSECTED")

write_rds(out, file = paste0(wd, "DESEQ2TOPGO.rds"))


# QUIT ===
