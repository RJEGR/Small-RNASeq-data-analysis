# LOAD
# RNA_LOCATION_DB
# DESEQ2MIRGENE2SRNA_REGULATORY
# CREATE 2 DB W/
# 1) Name/Family (list), gene_id, description (protein), group/contrast
# 2) Name/Family (list), go.id, term, group/contrast
# SAVE DF W/ NCOLS AND 147 ROWS MIRS


rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)


wd <- "~/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out/"

print(DB <- read_tsv(paste0(wd, "/RNA_LOCATION_DB.tsv")) %>% filter(SRNAtype == "miR"))

print(RES.P <- read_tsv(paste0(wd, "DESEQ_RES_P.tsv")) %>% filter( padj < 0.05 & abs(log2FoldChange) > 1))

.UPSETDF <- read_rds(paste0(wd, "UPSETDF.rds"))

EXCLUSIVE_MIRS <- .UPSETDF[[1]] %>% ungroup()

INTERSECTED_MIRS <- .UPSETDF[[2]]  %>% ungroup()

EXCLUSIVE_MIRS %>% count(CONTRAST, SIGN, sort = T) # ???

INTERSECTED_MIRS %>% count(CONTRAST, SIGN, sort = T)

wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/FUNCTIONAL_MIR_ANNOT/"

print(.SRNA2GO <- read_tsv(paste0(wd, "SRNA_REGULATORY_FUNCTION_DB.tsv")))

# 1) ====

print(DB1 <- .SRNA2GO %>% 
  mutate(query = strsplit(query, ";")) %>%
  unnest(query) %>%
  separate(query, into = c("query", "arm"), sep = "[.]") %>% 
  filter(arm == "mature") %>%
  filter(predicted == "BOTH") %>%
  distinct(query, gene_id) %>%
  dplyr::rename("Name" = "query") %>%
  right_join(distinct(RES.P, Name, Family)) %>% 
  group_by(gene_id) %>%
  summarise(n_mirs = n(), across(Family, .fns = paste_go), .groups = "drop_last") %>%
  left_join(distinct(.SRNA2GO, gene_id, description)))


# DB1 %>% view()

N_TARGETS <- .SRNA2GO %>% 
  mutate(query = strsplit(query, ";")) %>%
  unnest(query) %>%
  separate(query, into = c("query", "arm"), sep = "[.]") %>% 
  filter(arm == "mature") %>%
  filter(predicted == "BOTH") %>% 
  group_by(query) %>%
  summarise(n_targets = n()) %>%
  dplyr::rename("Name" = "query") %>%
  right_join(distinct(RES.P, Name, Family))

# GENERATE GO BY MIR ===
# read_tsv(paste0(wd, "SRNA2GO.tsv")) %>%

SRNA2GO <- .SRNA2GO %>% 
  mutate(query = strsplit(query, ";")) %>%
  unnest(query) %>%
  separate(query, into = c("query", "arm"), sep = "[.]") %>% 
  filter(arm == "mature") %>%
  filter(predicted == "BOTH") %>%
  dplyr::select(GO.ID, query)
  
SRNA2GO <- split(strsplit(SRNA2GO$GO.ID, ";") , SRNA2GO$query)

SRNA2GO <- lapply(SRNA2GO, unlist)

str(query.names <- RES.P %>% ungroup() %>% distinct(Name) %>% pull(Name))

# str(lapply(SRNA2GO, length))

str(query.names <- query.names[query.names %in% names(SRNA2GO)])

cat("\nUsing ",length(query.names), " Names...\n")

query.p <- RES.P %>% group_by(Name) %>% sample_n(1) %>% pull(pvalue, name = Name)

# allResg <- GOenrichment(query.p, query.names, SRNA2GO, Nodes = 7)

allRes <- list()

for (i in 1:length(query.names)) {
  
  q <- query.names[i]
  
  cat("\nUsing ",q, " Names...\n")
  
  
  # p <- 0.05
  p <- query.p[i]
  
  df <- GOenrichment(p, q, SRNA2GO, Nodes = 3)
  
  allRes[[i]] <- data.frame(df, Name = q)
}


DF <- do.call(rbind, allRes) %>% as_tibble() %>% left_join(distinct(RES.P, Name, Family))

# DF <- DF %>% mutate(SIGN = WHICH_SIGN, CONTRAST = CONTRAST) 


DB2 <- DF %>% 
  group_by(GO.ID, Term) %>%
  summarise(n = n(),
    across(Family, .fns = paste_go), 
    .groups = "drop_last")


DB2 %>% view()
DB1 %>% view()


# RUN TOPGO USING CONTRAST

WHICH_SIGN <- "A) 24 HPF"

EXCLUSIVE_MIRS %>% unnest(CONTRAST_DE) %>% filter(SIGN == WHICH_SIGN) %>% count(SIGN, CONTRAST)

str(CONTRAST <- EXCLUSIVE_MIRS %>% ungroup() %>% filter(SIGN == WHICH_SIGN) %>% distinct(CONTRAST) %>% pull())

str(query.names <- EXCLUSIVE_MIRS %>% ungroup() %>% filter(SIGN == WHICH_SIGN) %>% distinct(Name) %>% pull(Name))

str(query.names <- query.names[query.names %in% names(SRNA2GO)])

cat("\nUsing ",length(query.names), " Names...\n")

# query.p <- c(rep(0.05, length(query.names)))

query.p <- RES.P %>% group_by(Name) %>% sample_n(1) %>% pull(pvalue, name = Name)

str(query.p <- query.p[names(query.p) %in% query.names])

allRes1 <- GOenrichment(query.p, query.names, SRNA2GO, Nodes = 12)

allRes1 <- allRes1 %>% mutate(SIGN = WHICH_SIGN, CONTRAST = CONTRAST)

# B)

# drop.names <- EXCLUSIVE_MIRS %>%
#   group_by(CONTRAST, SIGN) %>%
#   summarise(across(Name, .fns = list), n = n()) %>%
#   filter(n == 1) %>%
#   unnest(Name) %>% 
#   pull(Name)

WHICH_SIGN <- "B) 110 HPF"

EXCLUSIVE_MIRS %>% unnest(CONTRAST_DE) %>% filter(SIGN == WHICH_SIGN) %>% count(SIGN, CONTRAST)

str(CONTRAST <- EXCLUSIVE_MIRS %>% ungroup() %>% filter(SIGN == WHICH_SIGN) %>% distinct(CONTRAST) %>% pull())

str(query.names <- EXCLUSIVE_MIRS %>% ungroup() %>% filter(SIGN == WHICH_SIGN) %>% distinct(Name) %>% pull(Name))

str(query.names <- query.names[query.names %in% names(SRNA2GO)])

cat("\nUsing ",length(query.names), " Names...\n")

# query.p <- c(rep(0.05, length(query.names)))

query.p <- RES.P %>% group_by(Name) %>% sample_n(1) %>% pull(pvalue, name = Name)

str(query.p <- query.p[names(query.p) %in% query.names])

# allRes1 <- GOenrichment(query.p, query.names, SRNA2GO, Nodes = 12)

# allRes1 <- allRes1 %>% mutate(SIGN = WHICH_SIGN, CONTRAST = CONTRAST)

DF < list()

for(j in 1:length(CONTRAST)) {
  
  i <- j
  
  cat("\nRunning ",CONTRAST[i], "\n")
  
  query.names <- EXCLUSIVE_MIRS %>% ungroup() %>% 
    filter(SIGN == WHICH_SIGN) %>% 
    filter(CONTRAST %in% CONTRAST[i]) %>% distinct() %>% pull(Name)
  
  # RES.P %>% filter(Name %in% query.names) %>% filter(CONTRAST %in% CONTRAST[i]) %>% view()
  
  str(query.names <- query.names[query.names %in% names(SRNA2GO)])
  
  query.p <- RES.P %>% group_by(Name) %>% sample_n(1) %>% pull(pvalue, name = Name)
  
  str(query.p <- query.p[names(query.p) %in% query.names])
  
  
  cat("\nUsing ",length(query.names), " Names...\n")

  allRes <- GOenrichment(query.p, query.names, SRNA2GO, Nodes = 10)
  
  allRes <- allRes %>% mutate(CONTRAST = CONTRAST[i]) %>% as_tibble()
  
  DF[[i]] <- allRes
  
}

