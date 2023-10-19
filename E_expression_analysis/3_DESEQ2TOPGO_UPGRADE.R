
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

DF1 <- do.call(rbind, allRes) %>% as_tibble() %>% mutate(facet = CONTRAST)

# df <- GOenrichment(query.p, query.names, SRNA2GO, Nodes = 10, onto = "BP", mapping = NULL)

# 110 HPF UP EXPRESSED BY OA ====
WCONTRAST <- "CONTRAST_B"

str(query.names <- RES.P %>% filter(CONTRAST %in% WCONTRAST) %>% 
  filter(log2FoldChange < 0) %>%
  distinct(Name) %>% pull())

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

DF3 <- do.call(rbind, allRes) %>% as_tibble() %>% mutate(facet = WCONTRAST, col = "Down")

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

