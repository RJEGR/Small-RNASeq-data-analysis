
# UPGRADE VERSION OF 3_DESEQ2TOPGO.R

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

wd <- "~/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out/"

RES.P <- read_tsv(paste0(wd, "DESEQ_RES.tsv")) %>% filter( padj < 0.05 & abs(log2FoldChange) > 1)

RES.P %>% dplyr::count(CONTRAST)

wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/FUNCTIONAL_MIR_ANNOT/"

print(.SRNA2GO <- read_tsv(paste0(wd, "SRNA_REGULATORY_FUNCTION_DB.tsv")))

URL <- "https://raw.githubusercontent.com/RJEGR/Small-RNASeq-data-analysis/master/FUNCTIONS.R"

source(URL)

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

.SRNA2GO <- SRNA2GO

SRNA2GO <- split(strsplit(SRNA2GO$GO.ID, ";") , SRNA2GO$query)

SRNA2GO <- lapply(SRNA2GO, unlist)


# IF MAPPING 

# org.Dr.eg.db # ZEBRAFISH
# org.Ce.eg.db # WORM
# org.Xl.eg.db # XENOPUS

WHICH_DB <- "org.Ce.eg.db"

# BiocManager::install(WHICH_DB)

# MAPPING <- GOSemSim::godata(WHICH_DB, ont="BP") # TAKE A LONG (3- 5 min)

# write_rds(MAPPING, paste0(wd, WHICH_DB, ".rds"))

MAPPING <- read_rds(paste0(wd, WHICH_DB, ".rds"))

# 24 HPF UP EXPRESSED BY OA ====
WCONTRAST <- "CONTRAST_A"

query.names <- RES.P %>% filter(CONTRAST %in% WCONTRAST) %>% 
  filter(log2FoldChange < 0) %>%
  distinct(Name) %>% pull()

names_under_OA_24 <- query.names

str(query.names <- query.names[query.names %in% names(SRNA2GO)])

cat("\nUsing ",length(query.names), " QUERY Names...\n")

query.p <- RES.P %>% 
  group_by(Name) %>% sample_n(1) %>% 
  pull(pvalue, name = Name)

query.p <- query.p[match(query.names, names(query.p))]

identical(names(query.p), query.names)

allRes <- list()

for (i in 1:length(query.names)) {
  
  q <- query.names[i]
  
  cat("\nUsing ",q, " Names...\n")
  
  n <- lapply(SRNA2GO[names(SRNA2GO) %in% q], length)[[1]]
  
  cat("\nUsing ",n, " GO terms\n")
  
  p <- query.p[i]
  
  df <- GOenrichment(p, q, SRNA2GO, Nodes = 10, onto = "BP") # mapping = NULL
  
  allRes[[i]] <- data.frame(df, Name = q)
}

DF1 <- do.call(rbind, allRes) %>% as_tibble() %>% mutate(facet = WCONTRAST)

# df <- GOenrichment(query.p, query.names, SRNA2GO, Nodes = 10, onto = "BP", mapping = NULL)

# 110 HPF UP EXPRESSED BY OA ====
WCONTRAST <- "CONTRAST_B"

str(query.names <- RES.P %>% filter(CONTRAST %in% WCONTRAST) %>% 
  filter(log2FoldChange < 0) %>%
  distinct(Name) %>% pull())

names_under_OA_110 <- query.names

str(query.names <- query.names[query.names %in% names(SRNA2GO)])

cat("\nUsing ",length(query.names), " QUERY Names...\n")

query.p <- RES.P %>% 
  group_by(Name) %>% sample_n(1) %>% 
  pull(pvalue, name = Name)

query.p <- query.p[match(query.names, names(query.p))]

identical(names(query.p), query.names)

allRes <- list()

for (i in 1:length(query.names)) {
  
  q <- query.names[i]
  
  cat("\nUsing ",q, " Names...\n")
  
  n <- lapply(SRNA2GO[names(SRNA2GO) %in% q], length)[[1]]
  
  cat("\nUsing ",n, " GO terms\n")
  
  p <- query.p[i]
  
  df <- GOenrichment(p, q, SRNA2GO, Nodes = 10, onto = "BP", mapping = NULL)
  
  allRes[[i]] <- data.frame(df, Name = q)
}

DF2 <- do.call(rbind, allRes) %>% as_tibble() %>% mutate(facet = WCONTRAST)


# 110 HPF DOWN EXPRESSED BY OA ====
WCONTRAST <- "CONTRAST_B"

str(query.names <- RES.P %>% 
    filter(CONTRAST %in% WCONTRAST) %>% 
    filter(log2FoldChange > 0) %>%
    # filter(!grepl("Cluster", Family)) %>%
    distinct(Name) %>% pull())

names_down_under_OA_110 <- query.names

str(query.names <- query.names[query.names %in% names(SRNA2GO)])

cat("\nUsing ",length(query.names), " QUERY Names...\n")

query.p <- RES.P %>% 
  group_by(Name) %>% sample_n(1) %>% 
  pull(pvalue, name = Name)

query.p <- query.p[match(query.names, names(query.p))]

identical(names(query.p), query.names)

allRes <- list()

for (i in 1:length(query.names)) {
  
  q <- query.names[i]
  
  cat("\nUsing ",q, " Names...\n")
  
  n <- lapply(SRNA2GO[names(SRNA2GO) %in% q], length)[[1]]
  
  cat("\nUsing ",n, " GO terms\n")
  
  p <- query.p[i]
  
  df <- GOenrichment(p, q, SRNA2GO, Nodes = 10, onto = "BP", mapping = NULL)
  
  allRes[[i]] <- data.frame(df, Name = q)
}

DF3 <- do.call(rbind, allRes) %>% as_tibble() %>% mutate(facet = WCONTRAST)

recode_to_contrat <- structure(c("24 HPF", "110 HPF"), names = c("CONTRAST_A", "CONTRAST_B"))

recode_to_de <- structure(c("pH 7.6", "pH 8.0"), names = c("Up", "Down"))

scale_col <- structure(c("#3B4992FF","#EE0000FF"), names = c("pH 8.0", "pH 7.6"))


rbind(DF1, DF2, DF3) %>%
  dplyr::mutate(facet = dplyr::recode_factor(facet, !!!recode_to_contrat, .ordered = T)) %>%
  group_by(facet, col, Term) %>%
  count(Term, sort = T) %>%
  group_by(facet, col) %>%
  arrange(n) %>% 
  ungroup() %>% 
  mutate(order = row_number()) %>%
  mutate(Annotated = ifelse(col == "Up", -n, n)) %>%
  mutate(col = factor(col, levels = c("Up","Down"))) %>%
  dplyr::mutate(col = dplyr::recode_factor(col, !!!recode_to_de)) %>%
  # mutate(ordering = -as.numeric(col) + Annotated, Term = fct_reorder(Term, ordering, .desc = T)) %>%
  mutate(Term = fct_reorder(Term, order, .desc = F)) %>%
  ggplot(aes(y = Term, x = Annotated, fill = col, color = col)) + # 
  # geom_col() +
  geom_segment(aes(x = Annotated, xend = 0, yend = Term)) +
  ggh4x::facet_nested(facet+col~ ., nest_line = F, scales = "free_y", space = "free_y") +
  scale_color_manual("", values = scale_col) +
  scale_fill_manual("", values =  scale_col) +
  theme_bw(base_family = "GillSans", base_size = 12) +
  labs(x = "", y = "") +
  theme(legend.position = "top",
    strip.background = element_rect(fill = 'grey89', color = 'white'),
    panel.border = element_blank(),
    plot.title = element_text(hjust = 0),
    plot.caption = element_text(hjust = 0),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_line(linetype = "dashed", linewidth = 0.5),
    axis.text.y = element_text(angle = 0, size = 5),
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10)) -> p

ggsave(p, filename = 'DESEQ2TOPGO_OA_RESPONSE.png', path = wd, width = 3.5, height = 10, device = png, dpi = 300)

# TEST REVIGO INSTEAD ====

SRNA2GO <- .SRNA2GO %>% mutate(DE = NA)

SRNA2GO <- SRNA2GO %>% mutate(DE = ifelse(query %in% names_under_OA_24, "24 HPF", DE))
SRNA2GO <- SRNA2GO %>% mutate(DE = ifelse(query %in% names_under_OA_110, "110 HPF", DE))
SRNA2GO <- SRNA2GO %>% mutate(DE = ifelse(query %in% names_down_under_OA_110, "-110 HPF", DE))

DEDF <- SRNA2GO %>% drop_na(DE) %>% distinct(query, DE) %>% dplyr::rename("Name" = "query")

GODF <- SRNA2GO %>% 
  drop_na(DE) %>%
  mutate(GO.ID = strsplit(GO.ID, ";")) %>%
  unnest(GO.ID) %>%
  distinct(DE, GO.ID) %>% # <-- OPTIONAL
  group_by(DE) %>%
  summarise(
    across(GO.ID, .fns = paste_go), 
    .groups = "drop_last")

GODF <- split(strsplit(GODF$GO.ID, ";") , GODF$DE)

GODF <- lapply(GODF, unlist)

orgdb <- "org.Ce.eg.db"

semdata <- read_rds(paste0(wd, orgdb, ".rds"))

OUT <- lapply(GODF, function(x) SEMANTIC_SEARCH(x, semdata = semdata))

print(data <- dplyr::bind_rows(OUT, .id = "DE") %>% as_tibble())

write_tsv(data, file = paste0(wd, "DESEQ2REVIGO_BY_GROUP.tsv"))

which_dup <- data %>% distinct(DE ,parentTerm) %>% dplyr::count(parentTerm) %>% filter(n == 2) %>% pull(parentTerm)

recode_to <- structure(c("24 hpf", "110 hpf", "110 hpf"),names = c("24 HPF", "110 HPF", "-110 HPF"))

data %>%
  # filter(!parentTerm %in% which_dup) %>%
  group_by(DE ,parentTerm) %>%
  summarise(size = sum(size)) %>%
  mutate(size = size / max(size)) %>%
  filter(size > 0) %>%
  ungroup() %>%
  mutate(size = ifelse(DE %in% "-110 HPF", -size, size)) %>%
  dplyr::mutate(DE = dplyr::recode_factor(DE, !!!recode_to)) %>%
  # mutate(DE = factor(DE, levels = c("24 HPF", "110 HPF"))) %>%
  mutate(parentTerm = fct_reorder2(parentTerm, DE, size, .desc = F)) %>%
  ggplot(aes(y = parentTerm, x = size, fill = DE, color = DE)) + # 
  geom_vline(xintercept = 0, linetype="dashed", alpha=0.5) +
  geom_segment(aes(x = size, xend = 0, yend = parentTerm), size = 3) +
  ggh4x::facet_nested(DE~ ., nest_line = F, scales = "free_y", space = "free_y") +
  theme_bw(base_family = "GillSans", base_size = 12) +
  labs(x = "", y = "") +
  scale_color_manual("", values = c("#DADADA", "#D4DBC2")) +
  scale_fill_manual("", values =  c("#DADADA", "#D4DBC2")) +
  theme(legend.position = "top",
    strip.background = element_blank(),
    strip.text = element_blank(),
    panel.border = element_blank(),
    plot.title = element_text(hjust = 0),
    plot.caption = element_text(hjust = 0),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_line(linetype = "dashed", linewidth = 0.5),
    axis.text.y = element_text(angle = 0, size = 10),
    axis.text.x = element_text(size = 10)) -> p


ggsave(p, filename = 'DESEQ2REVIGO_OA_RESPONSE.png', path = wd, width = 6, height = 6, device = png, dpi = 300)

# OR BY MIR ====

str(SORT_MIRS <- c("MIR-278","LET-7","MIR-133",
  "Cluster_55760","Cluster_55776","MIR-2","MIR-315", "MIR-153", "BANTAM", "MIR-190", "MIR-2722", "MIR-1988", 
  "MIR-92", "miR-277b", "MIR-216"))

GODF <- SRNA2GO %>% drop_na(DE) 

GODF <- split(strsplit(GODF$GO.ID, ";") , GODF$query)

GODF <- lapply(GODF, unlist)

OUT <- lapply(GODF, function(x) SEMANTIC_SEARCH(x, semdata = semdata))

print(data2 <- dplyr::bind_rows(OUT, .id = "Name") %>% as_tibble())

data2 <- data2 %>% left_join(distinct(RES.P, Name, Family))

write_tsv(data2, file = paste0(wd, "DESEQ2REVIGO_BY_MIR.tsv"))

# CONTINUE here =====

data2 <- read_tsv(paste0(wd, "DESEQ2REVIGO_BY_MIR.tsv"))


which_dup <- data2 %>% distinct(Name ,parentTerm) %>% dplyr::count(parentTerm) %>% filter(n == 2) %>% pull(parentTerm)

recode_to <- structure(c("A) 24 hpf", "B) 110 hpf", "B) 110 hpf"),names = c("24 HPF", "110 HPF", "-110 HPF"))


sum_data2 <- data2 %>%
  # filter(!parentTerm %in% which_dup) %>%
  left_join(DEDF) %>%
  left_join(distinct(RES.P, Name, Family)) %>%
  group_by(Family ,parentTerm, DE) %>%
  summarise(size = sum(size)) %>%
  group_by(Family) %>%
  mutate(size = size / max(size)) %>%
  mutate(size = ifelse(DE %in% "110 HPF", -size, size)) %>%
  mutate(size = ifelse(DE %in% "24 HPF", -size, size)) %>%
  dplyr::mutate(DE = dplyr::recode_factor(DE, !!!recode_to)) 

# DUPLICATE MIR-133 IN THE 24 HPF
sum_data2 <- sum_data2 %>% filter(Family == "MIR-133") %>% mutate(DE = "A) 24 hpf", size = -size) %>% rbind(sum_data2)

sum_data2 %>%
  # filter(abs(size) > 0) %>%
  ungroup() %>%
  mutate(parentTerm = fct_reorder2(parentTerm, Family, size, .desc = F)) %>%
  dplyr::mutate(DE = factor(DE, levels = c("A) 24 hpf", "B) 110 hpf"))) %>%
  ggplot(aes(y = parentTerm, x = Family, fill = size, color = size)) +
  labs(x = "", y = "") +
  ggh4x::facet_nested(DE~ ., nest_line = F, scales = "free_y", space = "free_y") +
  geom_tile(color = 'white', size = 0.7, width = 1) +
  ggsci::scale_fill_gsea(name = "",reverse = T) +
  ggh4x::scale_x_dendrogram(hclust = hc_srna, position = "bottom", labels = NULL) +
  guides(x.sec = guide_axis_manual(labels = srna_order, label_size = 10)) +
  guides(fill = guide_colorbar(barwidth = unit(2.5, "in"),
    barheight = unit(0.1, "in"), label.position = "top",
    alignd = 0.5,
    ticks.colour = "black", ticks.linewidth = 0.5,
    frame.colour = "black", frame.linewidth = 0.5,
    label.theme = element_text(size = 10))) +
  theme_bw(base_family = "GillSans", base_size = 12) +
    theme(legend.position = "top",
      strip.background = element_rect(fill = 'grey89', color = 'white'),
      panel.border = element_blank(),
      plot.title = element_text(hjust = 0),
      plot.caption = element_text(hjust = 0),
      panel.grid.minor.y = element_blank(),
      # axis.ticks = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.text.y = element_text(angle = 0, size = 10),
      axis.text.x.top = element_text(angle = 90, size = 10, hjust = 0, vjust = 0)) -> p2

# p2

ggsave(p2, filename = 'DESEQ2REVIGO_OA_RESPONSE_HEATMAP.png', path = wd, width = 7, height = 10, device = png, dpi = 300)

# create supp. table ====

data2 %>%
  group_by(Family, parentTerm) %>%
  # distinct(Name, parentTerm) %>%
  summarise(
    across(term, .fns = paste_go), 
    .groups = "drop_last") %>% view()

# COMPLEMENT W/ USING GENE_ID AS Y AXIS

print(.SRNA2GO <- read_tsv(paste0(wd, "SRNA_REGULATORY_FUNCTION_DB.tsv")))


DF <- .SRNA2GO %>%
  filter(predicted == "BOTH") %>%
  mutate(query = strsplit(query, ";")) %>%
  unnest(query) %>% 
  distinct(query, gene_id) %>%
  separate(query, into = c("query", "arm"), sep = "[.]") %>%
  filter(arm == "mature") %>%
  dplyr::rename("Name" = "query" ) %>%
  right_join(distinct(RES.P, Name, Family))
  

DF <- DF %>% mutate(DE = NA)

DF <- DF %>% mutate(DE = ifelse(Name %in% names_under_OA_24, "24 HPF", DE))
DF <- DF %>% mutate(DE = ifelse(Name %in% names_under_OA_110, "110 HPF", DE))
DF <- DF %>% mutate(DE = ifelse(Name %in% names_down_under_OA_110, "-110 HPF", DE))

DF <- DF %>% distinct(Family, gene_id, DE)

DF <- DF %>% drop_na(DE) %>% filter(Family %in% SORT_MIRS) # <- IF USING ONLY DEGS (198 targets)

MAT <- DF %>%  
  distinct(Family,gene_id) %>%
  mutate(size = 1) %>%
  pivot_wider(names_from = Family, values_from = size, values_fill = 0) 

head(MAT <- as(MAT, "matrix"))

rownames(MAT) <- MAT[,1]

dim(MAT <- MAT[,2:ncol(MAT)] )

gene_dist = dist(MAT, method='euclidean')
hc_genes = hclust(gene_dist, method='complete')
hc_order <- hc_genes$labels[hc_genes$order]

srna_dist = dist(t(MAT), method='euclidean')
hc_srna = hclust(srna_dist, method='complete')
srna_order <- hc_srna$labels[hc_srna$order]

plot(hc_srna)

library(ggh4x)

DF %>% distinct(gene_id)

DF %>%
  mutate(size = 1, facet = "C) Gene target predicted") %>%
  # mutate(gene_id = factor(gene_id, levels = hc_order)) %>%
  ggplot(aes(x = gene_id, y = Family, fill = size)) +
  labs(x = "", y = "") +
  ggh4x::facet_nested(~ facet, nest_line = F) + 
  geom_tile(color = 'white', fill = "grey50", size = 0.7, width = 1) + # #2196F3
  ggh4x::scale_x_dendrogram(hclust = hc_genes, position = "top", labels = NULL) +
  guides(x.sec = guide_axis_manual(labels = hc_order, label_size = 3.5)) +
  ggh4x::scale_y_dendrogram(hclust = hc_srna, position = "right", labels = NULL) +
  # guides(y.sec = guide_axis_manual(labels = srna_order, label_size = 3.5)) +
  theme_bw(base_family = "GillSans", base_size = 12) +
  theme(legend.position = "none",
    strip.background = element_rect(fill = 'grey89', color = 'white'),
    panel.border = element_blank(),
    plot.title = element_text(hjust = 0),
    plot.caption = element_text(hjust = 0), axis.ticks = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.y = element_text(angle = 0, size = 12),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 10)) -> p2
  
# ggsave(p, filename = 'TARGETS_HEATMAP_DEGS.png', path = wd, width = 2, height = 6, device = png, dpi = 300)

# SEE 2_DESEQ2BARPLOT TO BIND HEATMAP OF MULTIPLE TARGET PROCESS TO DEGS.

psave <- p+p2 + plot_layout(widths = c(2, 8))

ggsave(psave, filename = 'TARGETS_HEATMAP2BARPLOT_DEGS_EN.png', path = wd, width = 10, height = 3.5, device = png, dpi = 300)

