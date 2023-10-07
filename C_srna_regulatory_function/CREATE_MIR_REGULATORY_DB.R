# LOAD
# RNA_LOCATION_DB
# DESEQ2MIRGENE2SRNA_REGULATORY
# CREATE 2 DB W/
# 1) Name/Family (list), gene_id, description (protein), group/contrast
# 2) Name/Family (list), go.id, term, group/contrast
# SAVE DF W/ NCOLS AND 147 ROWS MIRS

# The GO terms that should be used for the action of the miRNA on gene expression
# should be Biological Process terms OR Molecular Function terms

# https://wiki.geneontology.org/index.php/MicroRNA_GO_annotation_manual

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)


wd <- "~/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out/"

print(DB <- read_tsv(paste0(wd, "/RNA_LOCATION_DB.tsv")) %>% filter(SRNAtype == "miR"))

print(RES.P <- read_tsv(paste0(wd, "DESEQ_RES.tsv")) %>% filter( padj < 0.05 & abs(log2FoldChange) > 1))

.UPSETDF <- read_rds(paste0(wd, "UPSETDF.rds"))

EXCLUSIVE_MIRS <- .UPSETDF[[1]] %>% ungroup()

INTERSECTED_MIRS <- .UPSETDF[[2]]  %>% ungroup()

EXCLUSIVE_MIRS %>% count(CONTRAST, SIGN, sort = T) # ???

INTERSECTED_MIRS %>% count(CONTRAST, SIGN, sort = T)

wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/FUNCTIONAL_MIR_ANNOT/"

print(.SRNA2GO <- read_tsv(paste0(wd, "SRNA_REGULATORY_FUNCTION_DB.tsv")))

# BIND SWISSPROT TO SRNA2GO
ref_path <- "/Users/cigom/Documents/MIRNA_HALIOTIS/ENSEMBLE/ANNOTATIONS/"

annot_f <- list.files(path = ref_path, pattern = "Trinotate.xls", full.names = T)

annot <- read_tsv(annot_f, na = ".")

names(annot)[c(1,2)] <- c("gene_id", "transcript_id")

annot <- annot %>% drop_na(gene_ontology_BLASTX)

blastx_df <- split_blast(annot, hit = "BLASTX") # HARD TO RUN THIS TIME

DB <- blastx_df %>% dplyr::filter(domain == "Eukaryota")

# blastx_df %>% dplyr::count(genus, sort = T)

DB <- DB %>% 
  group_by(gene) %>%
  arrange(desc(identity)) %>%
  sample_n(1) %>%
  ungroup()

# DB %>% dplyr::count(genus, sort = T)

genome_feat <- read_rds("/Users/cigom/Documents/MIRNA_HALIOTIS/ENSEMBLE/genome_features.rds")[[1]]

genome_feat <- genome_feat %>% dplyr::select(gene_id, transcript_id)

DB <- DB %>% 
  dplyr::select(-align, -identity, -evalue, -gene_id, -protein,-lineage, -domain, -genus) %>%
  right_join(genome_feat, by = "transcript_id") 


# UNNEST MIRS ====

SRNA2GO_DE <- .SRNA2GO %>%
  filter(predicted == "BOTH") %>%
  mutate(query = strsplit(query, ";")) %>%
  unnest(query) %>%
  separate(query, into = c("query", "arm"), sep = "[.]") %>%
  filter(arm == "mature") %>%
  dplyr::rename("Name" = "query") %>%
  dplyr::select(-n, -n_rnas, -arm) %>%
  distinct() %>%
  right_join(distinct(RES.P, Name, Family))



# HOW POSSIBLE TARGETS REGUALTED BY DEGS? 167

SRNA2GO_DE <- DB %>% distinct(gene_id, uniprot, name) %>%
  right_join(SRNA2GO_DE, by = "gene_id") %>%
  distinct(Family, description, GO.ID, uniprot, name)

distinct(SRNA2GO_DE, description) # 167

# SRNA2GO_DE <- distinct(SRNA2GO_DE, Family, description, GO.ID )

# view(SRNA2GO_DE)

SRNA2GO_DE %>% dplyr::count(Family,sort = T)


# REMOVE uncharacterized VARIANT ?
# SRNA2GO_DE %>% filter(grepl("uncharacterized", description))

wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/FUNCTIONAL_MIR_ANNOT/"

write_tsv(SRNA2GO_DE, paste0(wd, "DESEQ2SRNATARGET.tsv"))

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

.SRNA2GO %>% 
  mutate(query = strsplit(query, ";")) %>%
  unnest(query) %>%
  separate(query, into = c("query", "arm"), sep = "[.]") %>% 
  filter(arm == "mature") %>%
  filter(predicted == "BOTH") %>% 
  group_by(query) %>%
  summarise(n_targets = n(), gene_id = paste_go(gene_id)) %>% view()

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

# DEFINE WHICH MIRS?

RES.OUT <- RES.P # %>% filter(CONTRAST %in% c("CONTRAST_A","CONTRAST_B"))

RES.OUT <- RES.P %>% filter(!grepl("Cluster", Family))

str(query.names <- RES.OUT %>% ungroup() %>% distinct(Name) %>% pull(Name))

# str(lapply(SRNA2GO, length))

str(query.names <- query.names[query.names %in% names(SRNA2GO)])

cat("\nUsing ",length(query.names), " QUERY Names...\n")

query.p <- RES.OUT %>% 
  group_by(Name) %>% sample_n(1) %>% pull(pvalue, name = Name)

query.p <- query.p[match(query.names, names(query.p))]

identical(names(query.p), query.names)

# allResg <- GOenrichment(query.p, query.names, SRNA2GO, Nodes = 7)

# BIOLOGICAL PROCES.

allRes <- list()

for (i in 1:length(query.names)) {
  
  q <- query.names[i]
  
  cat("\nUsing ",q, " Names...\n")
  
  n <- lapply(SRNA2GO[names(SRNA2GO) %in% q], length)[[1]]
  
  cat("\nUsing ",n, " GO terms\n")
  
  
  # p <- 0.05
  p <- query.p[i]
  
  df <- GOenrichment(p, q, SRNA2GO, Nodes = 10, onto = "BP")
  
  allRes[[i]] <- data.frame(df, Name = q)
}


BPDF <- do.call(rbind, allRes) %>% as_tibble() %>% 
  left_join(distinct(RES.P, Name, Family)) %>% mutate(Onto = "Biological Process")

# CELLULAR COMPOMENTS (OMIT)
# 
# allRes <- list()
# 
# for (i in 1:length(query.names)) {
#   
#   q <- query.names[i]
#   
#   cat("\nUsing ",q, " Names...\n")
#   
#   n <- lapply(SRNA2GO[names(SRNA2GO) %in% q], length)[[1]]
#   
#   cat("\nUsing ",n, " GO terms\n")
#   
#   
#   # p <- 0.05
#   p <- query.p[i]
#   
#   df <- GOenrichment(p, q, SRNA2GO, Nodes = 3, onto = "CC")
#   
#   allRes[[i]] <- data.frame(df, Name = q)
# }
# 
# CCDF <- do.call(rbind, allRes) %>% as_tibble() %>% left_join(distinct(RES.P, Name, Family)) %>% mutate(Onto = "Cellular Component")

# MOLLECULAR FUNCTION

allRes <- list()

for (i in 1:length(query.names)) {
  
  q <- query.names[i]
  
  cat("\nUsing ",q, " Names...\n")
  
  n <- lapply(SRNA2GO[names(SRNA2GO) %in% q], length)[[1]]
  
  cat("\nUsing ",n, " GO terms\n")
  
  
  # p <- 0.05
  p <- query.p[i]
  
  df <- GOenrichment(p, q, SRNA2GO, Nodes = 10, onto = "MF")
  
  allRes[[i]] <- data.frame(df, Name = q)
}

MFDF <- do.call(rbind, allRes) %>% as_tibble() %>% left_join(distinct(RES.P, Name, Family)) %>% mutate(Onto = "Molecular Function")

DF <- rbind(BPDF, MFDF)

DF %>% group_by(Family, Onto) %>% dplyr::count(Term, sort = T)

# 

BPDF

str(GO.ID <- BPDF %>% distinct(GO.ID) %>% pull())

# hsGO <- GOSemSim::godata('org.Hs.eg.db', ont="BP")

wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/FUNCTIONAL_MIR_ANNOT/"


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

hclust <- termSim %>% 
  # use distance metric
  dist(method = "euclidean") %>%
  # compute cMDS
  hclust()

cutree <- hclust %>% cutree(., 7)

go_order <- hclust$labels[hclust$order]

recode_to <- BPDF %>% filter(GO.ID %in% go_order) %>% distinct(GO.ID, Term)

recode_to <- structure(recode_to$Term, names = recode_to$GO.ID)

identical(sort(names(recode_to)),sort(go_order))

go_order <- recode_to[match(go_order, names(recode_to))]

identical(names(go_order),  hclust$labels[hclust$order])


# Sanity check

identical(cmds$GO.ID, names(cutree))

head(cmds <- cbind(cmds, cutree) %>% mutate(cutree = as.factor(cutree)))

BPDF <- BPDF %>% left_join(cmds)

BPDF %>% 
  # filter(Family %in% "MIR-190") %>%
  # arrange(desc(Term)) %>%
  # mutate(Term = factor(Term, levels = unique(Term))) %>%
  # mutate(Term = fct_reorder(Term, Annotated))
  ggplot(aes(x = Family, y = GO.ID, fill = Annotated) ) +
  geom_tile(color = 'white', linewidth = 0.2) +
  scale_x_discrete(position = 'top') +
  ggh4x::scale_y_dendrogram(hclust = hclust, position = "left", labels = NULL) +
  guides(y.sec = guide_axis_manual(labels = go_order, label_size = 5)) +
  theme_bw(base_size = 10, base_family = "GillSans") +
  labs(x = '', y = '') +
  facet_grid(Onto ~ ., scales = "free") +
  theme(
    legend.position = "none",
    # axis.text.y = element_blank(),
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(angle = 75, hjust = -0.05, vjust = 1),
    strip.background = element_rect(fill = 'grey89', color = 'white'),
    panel.border = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank()
    ) 

# GENERATE HEATMAP OF PROCESS ====

wd <- "~/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out/"

dds <- read_rds(paste0(wd, "/DDS_DESEQ2.rds"))

barplot(COUNT <- DESeq2::counts(dds, normalized = F, replaced = F)[query.names,])

barplot(COUNT <- DESeq2::varianceStabilizingTransformation(round(COUNT)))


sample_dist = dist(t(COUNT), method='euclidean')

hc_samples = hclust(sample_dist, method='complete')

hc_order <- hc_samples$labels[hc_samples$order]

srna_dist = dist(COUNT, method='euclidean')

hc_srna = hclust(srna_dist, method='complete')

srna_order <- hc_srna$labels[hc_srna$order]

recode_to <- RES.OUT %>% filter(Name %in% srna_order) %>% distinct(Name, Family)

recode_to <- structure(recode_to$Family, names = recode_to$Name)

identical(sort(names(recode_to)),sort(srna_order))

srna_order <- recode_to[match(srna_order, names(recode_to))]

identical(names(srna_order),  hc_srna$labels[hc_srna$order])

dplyr::count(BPDF, Term,sort = T)

.colData <- as_tibble(SummarizedExperiment::colData(dds)) %>% dplyr::select(LIBRARY_ID, colName,hpf, pH)

COUNT %>% 
  as_tibble(rownames = 'Name') %>%
  pivot_longer(-Name, names_to = "LIBRARY_ID") %>%
  left_join(.colData) %>%
  mutate(LIBRARY_ID = factor(LIBRARY_ID, levels = hc_order)) -> COUNT_LONG

library(ggh4x)

COUNT_LONG %>%
  ggplot(aes(x = LIBRARY_ID, y = Name, fill = value)) + 
  geom_tile(color = 'white', linewidth = 0.2) +
  scale_fill_viridis_c(option = "B", name = "Log2", direction = -1, na.value = "white") +
  ggh4x::facet_nested( ~ hpf+pH, nest_line = F, scales = "free", space = "free") +
  # scale_x_discrete(position = 'top') +
  # ggh4x::scale_x_dendrogram(hclust = hc_samples, position = 'top') +
  ggh4x::scale_y_dendrogram(hclust = hc_srna, position = "left", labels = NULL) +
  guides(y.sec = guide_axis_manual(labels = srna_order, label_size = 12)) +
  theme_bw(base_size = 12, base_family = "GillSans") +
  labs(x = '', y = '') +
  theme(
    axis.text.x = element_text(angle = 90, hjust = -0.15, vjust = 1))
  


# DF <- DF %>% mutate(SIGN = WHICH_SIGN, CONTRAST = CONTRAST) 

# QUIT ====

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

DF <- list()

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


allRes2 <- do.call(rbind, DF) %>% as_tibble() %>% mutate(SIGN = WHICH_SIGN)

view(allRes1)
view(allRes2)

identical(names(allRes1), names(allRes2))
