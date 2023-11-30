
# RICARDO GOMEZ-REYES
# DESEQ2REVIGO
rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

wd <- "~/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out/"

RES.P <- read_tsv(paste0(wd, "DESEQ_RES.tsv")) %>% filter( padj < 0.05 & abs(log2FoldChange) > 1)

RES.P %>% dplyr::count(CONTRAST)

semdata <- read_rds(paste0(wd, orgdb, ".rds"))

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

# SRNA2GO <- split(strsplit(SRNA2GO$GO.ID, ";") , SRNA2GO$query)

# SRNA2GO <- lapply(SRNA2GO, unlist)


# 2)
# FILTER OUT BY OA RESPONSE
# log2FoldChange < 0 MEANS EITHER, 24 OR 110 HPF UNDER LOW PH

str(filter.out <- RES.P %>% 
  filter(CONTRAST %in% c("CONTRAST_A", "CONTRAST_B")) %>%
  filter(log2FoldChange < 0) %>%
  distinct(Name) %>% pull())

# BY DEVELOPMENT
# log2FoldChange > 0 MEANS EITHER, CTRL OR LOW PH DUIRNG 24 HPF
# USING INTERSECTED MIRS TO SUBSET

str(query.names <- RES.P %>% 
    filter(CONTRAST %in% c("CONTRAST_C", "CONTRAST_D")) %>% 
    filter(log2FoldChange > 0) %>%
    dplyr::count(Name, sort = T) %>% # keeping if are presented in both, low and control pH
    filter(n == 2) %>% 
    filter(!Name %in% filter.out) %>% # MUST MATCH 46 
    pull(Name))

SRNA2GO <- SRNA2GO %>% mutate(DE = NA)

SRNA2GO <- SRNA2GO %>% mutate(DE = ifelse(query %in% query.names, "24 HPF", DE))
  
# log2FoldChange < 0 MEANS EITHER, CTRL OR LOW PH DUIRNG 110 HPF

str(query.names <- RES.P %>% 
    filter(CONTRAST %in% c("CONTRAST_C", "CONTRAST_D")) %>% 
    filter(log2FoldChange < 0) %>%
    dplyr::count(Name, sort = T) %>% # keeping if are presented in both, low and control pH
    filter(n == 2) %>% 
    filter(!Name %in% filter.out) %>% # MUST MATCH 38 INTERSECTED
    pull(Name))

SRNA2GO <- SRNA2GO %>% mutate(DE = ifelse(query %in% query.names, "110 HPF", DE))


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

SEMANTIC_SEARCH <- function(x, orgdb = "org.Ce.eg.db", semdata = semdata) {
  
  
  require(rrvgo)
  
  # semdata <- read_rds(paste0(wd, orgdb, ".rds"))
  
  x <- sort(x)
  
  SimMatrix <- calculateSimMatrix(x, 
    orgdb = orgdb,
    ont="BP", 
    semdata = semdata,
    method = 'Wang')
  
  data <- reduceSimMatrix(SimMatrix, threshold = 0.9, orgdb = orgdb) 
  
  y <- cmdscale(as.matrix(as.dist(1 - SimMatrix)), eig = TRUE, k = 2)
  
  data <- cbind(as.data.frame(y$points), data[match(rownames(y$points), data$go),])
  
  return(data)
}

str(GODF)

OUT <- lapply(GODF, SEMANTIC_SEARCH)

OUT

print(data <- dplyr::bind_rows(OUT, .id = "DE") %>% as_tibble())


p <- ggplot2::ggplot(data, aes(x = V1, y = V2,
  color = as.factor(cluster))) +
  ggplot2::geom_point(ggplot2::aes(size = size), alpha = 0.5) +
  ggplot2::scale_size_continuous(range = c(0, 10)) +
  # ggplot2::scale_x_continuous(name = "") +
  # ggplot2::scale_y_continuous(name = "") +
  ggplot2::theme_bw(base_family='GillSans', base_size = 14) +
  ggplot2::theme(legend.position = 'none',
    axis.line.x = ggplot2::element_blank(),
    axis.line.y = ggplot2::element_blank(),
    strip.background = element_rect(fill = 'grey86', color = 'white'),
    panel.border = element_blank()
  )

p <- p + facet_wrap(~DE)

data_subset <- data %>% distinct(parentTerm, .keep_all = T)

p + ggrepel::geom_text_repel(aes(label = parentTerm),
  data = data_subset,
  family = 'GillSans',
  max.overlaps = 10,
  box.padding = grid::unit(1, "lines"), size = 5) +
  labs(x = 'Dimension 1', y = 'Dimension 2') -> p


data

data %>%
  group_by(DE) %>%
  mutate(size = size / max(size)) %>%
  mutate(parentTerm = fct_reorder(parentTerm, size, .desc = T)) %>%
  ungroup() %>%
  mutate(DE = factor(DE, levels = c("24 HPF", "110 HPF"))) %>%
  ggplot(aes(y = parentTerm, x = DE, size = size, color = DE)) +
  geom_point(shape = 1) +
  # scale_color_manual('', values = structure(scale_col, names = scale_col) ) +
  # scale_fill_manual('', values = structure(scale_col, names = scale_col) ) +
  theme_bw(base_family = "GillSans", base_size = 12) +
  guides(color = "none") +
  theme(
    legend.position = "top",
    # axis.text.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 8),
    strip.background = element_rect(fill = 'grey89', color = 'white'),
    panel.border = element_blank(),
    plot.title = element_text(hjust = 0),
    plot.caption = element_text(hjust = 0),
    panel.grid.minor.x = element_blank(),
    panel.grid.major = element_blank(),
    # panel.grid.major = element_line(linewidth = 0.5, color = "grey89", linetype = "dashed"),
    panel.grid.major.x = element_blank())

# 3)

c("#DADADA", "#D4DBC2")

# WHICH BIOL. PROCESS ARE REGULATED BY DE-MIRS, BUT ARE SIMILAR PARENT TERM.

which_proc <- data %>% distinct(DE ,parentTerm) %>% dplyr::count(parentTerm) %>% filter(n == 2) %>% pull(parentTerm)

recode_to <- structure(c("A) 24 HPF", "B) 110 HPF"),names = c("24 HPF", "110 HPF"))

data %>%
  filter(!parentTerm %in% which_proc) %>%
  group_by(DE ,parentTerm) %>%
  summarise(size = sum(size)) %>%
  # group_by(DE) %>%
  mutate(size = size / max(size)) %>%
  ungroup() %>%
  # mutate(DE = factor(DE, levels = c("24 HPF", "110 HPF"))) %>%
  dplyr::mutate(DE = dplyr::recode_factor(DE, !!!recode_to)) %>%
  mutate(parentTerm = fct_reorder2(parentTerm, DE, size, .desc = F)) %>%
  ggplot(aes(y = parentTerm, x = size, fill = DE, color = DE)) + # 
  geom_segment(aes(x = size, xend = 0, yend = parentTerm), size = 4) +
  ggh4x::facet_nested(DE~ ., nest_line = F, scales = "free_y", space = "free_y") +
  theme_bw(base_family = "GillSans", base_size = 12) +
  labs(x = "", y = "") +
  scale_color_manual("", values = c("#DADADA", "#D4DBC2")) +
  scale_fill_manual("", values =  c("#DADADA", "#D4DBC2")) +
  theme(legend.position = "none",
    strip.background = element_rect(fill = 'grey89', color = 'white'),
    panel.border = element_blank(),
    plot.title = element_text(hjust = 0),
    plot.caption = element_text(hjust = 0),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_line(linetype = "dashed", linewidth = 0.5),
    axis.text.y = element_text(angle = 0, size = 10),
    axis.text.x = element_text(size = 10)) -> p

ggsave(p, filename = 'DESEQ2REVIGO_UP_BY_DEV.png', path = wd, width = 5, height = 4.5, device = png, dpi = 300)


data %>%
  filter(parentTerm %in% which_proc) %>%
  group_by(DE ,parentTerm) %>%
  summarise(size = sum(size)) %>%
  # group_by(DE) %>%
  mutate(size = size / max(size)) %>%
  ungroup() %>%
  mutate(DE = factor(DE, levels = c("24 HPF", "110 HPF"))) %>%
  mutate(parentTerm = fct_reorder2(parentTerm, DE, size, .desc = F)) %>%
  ggplot(aes(y = parentTerm, x = size, fill = DE, color = DE)) + # 
  geom_segment(aes(x = size, xend = 0, yend = parentTerm), size = 3) +
  ggh4x::facet_nested(~ DE, nest_line = F, scales = "free_y", space = "free_y") +
  theme_bw(base_family = "GillSans", base_size = 12) +
  labs(x = "", y = "") +
  scale_color_manual("", values = c("#DADADA", "#D4DBC2")) +
  scale_fill_manual("", values =  c("#DADADA", "#D4DBC2")) +
  theme(legend.position = "none",
    strip.background = element_rect(fill = 'grey89', color = 'white'),
    panel.border = element_blank(),
    plot.title = element_text(hjust = 0),
    plot.caption = element_text(hjust = 0),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_line(linetype = "dashed", linewidth = 0.5),
    axis.text.y = element_text(angle = 0, size = 10),
    axis.text.x = element_text(size = 10)) -> p

ggsave(p, filename = 'DESEQ2REVIGO_UP_BY_DEV_2.png', path = wd, width = 7, height = 3, device = png, dpi = 300)


#


