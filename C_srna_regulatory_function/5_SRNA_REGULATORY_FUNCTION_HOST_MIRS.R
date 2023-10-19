
# RICARDO GOMEZ-REYES
# SEARCH FOR HOST REGULATING MIR-MRNA
#


path <- "/Users/cigom/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out/"

.DB <- read_rds(paste0(path, "RNA_LOCATION_MIR_DB.rds"))

wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/FUNCTIONAL_MIR_ANNOT/"

# print(.SRNA2GO <- read_tsv(paste0(wd, "DESEQ2SRNATARGET.tsv"))) # <- HERE UNIPROT ID

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
    .groups = "drop_last") %>%
  separate(query, into = c("query", "arm"), sep = "[.]") %>%
  filter(arm == "mature") %>%
  rename("query" = "Name")
  


# SRNA2GO <- split(strsplit(SRNA2GO$GO.ID, ";") , SRNA2GO$query)

# SRNA2GO <- lapply(SRNA2GO, unlist)


# HOST GENE REGULATING ??

LOCI <- .SRNA2GO %>%
  filter(predicted == "BOTH") %>%
  mutate(query = strsplit(query, ";")) %>%
  unnest(query) %>% 
  distinct(query, gene_id) %>%
  separate(query, into = c("query", "arm"), sep = "[.]") %>%
  filter(arm == "mature") %>%
  rename("query" = "Name", "gene_id" = "target") %>% 
  dplyr::select(-arm)
  
# FROM THOSE INTRAGENIC (66 MIRS), ARE ANY REGULATING HOST?. R: ZERO

LOCI2TARGETDB <- .DB %>% distinct(Name, gene_id) %>% drop_na(gene_id) %>% right_join(LOCI) 

LOCI2TARGETDB %>% filter(gene_id  ==  target)

LOCI2TARGETDB %>% distinct(gene_id, target) %>% filter(gene_id  ==  target)

write_tsv(LOCI2TARGETDB, file = paste0(wd, "LOCI2TARGETDB.tsv"))

# GO ENRICHMENT BY WGCNA MODULE ====
# USING REVIGO
# https://bioconductor.org/packages/release/bioc/html/rrvgo.html

GODF <- .DB %>% distinct(Name, WGCNA) %>%
  right_join(SRNA2GO) %>%
  mutate(GO.ID = strsplit(GO.ID, ";")) %>%
  unnest(GO.ID) %>%
  distinct(WGCNA, GO.ID) %>% # <-- OPTIONAL
  group_by(WGCNA) %>%
  summarise(
    across(GO.ID, .fns = paste_go), 
    .groups = "drop_last")

GODF <- split(strsplit(GODF$GO.ID, ";") , GODF$WGCNA)

GODF <- lapply(GODF, unlist)


SEMANTIC_SEARCH <- function(x, orgdb = "org.Ce.eg.db") {

  
  require(rrvgo)
  
  semdata <- read_rds(paste0(wd, orgdb, ".rds"))
  
  x <- sort(x)
  
  SimMatrix <- calculateSimMatrix(x, 
    orgdb = orgdb,
    ont="BP", 
    semdata = semdata,
    method = 'Wang')
  
  data <- reduceSimMatrix(SimMatrix, threshold = 0.9, orgdb = orgdb) 
  
  y <- cmdscale(as.matrix(as.dist(1 - SimMatrix)), eig = TRUE, k = 2)
  
  data <- cbind(as.data.frame(y$points), 
    data[match(rownames(y$points), data$go),])
  
  return(data)
}

str(GODF)

OUT <- lapply(GODF, SEMANTIC_SEARCH)

OUT

head(data <- dplyr::bind_rows(OUT, .id = "WGCNA") %>% as_tibble())

write_tsv(data, file = paste0(wd, "WGCNA_MIRS2SEMANTIC.tsv"))

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

p <- p + facet_wrap(~WGCNA)

data_subset <- data %>% distinct(parentTerm, .keep_all = T)

p + ggrepel::geom_text_repel(aes(label = parentTerm),
  data = data_subset,
  family = 'GillSans',
  max.overlaps = 10,
  box.padding = grid::unit(1, "lines"), size = 5) +
  labs(x = 'Dimension 1', y = 'Dimension 2') -> p

mod_membr_1 <- c("grey", "black", "yellow", "pink", "turquoise")

mod_membr_2 <- c("brown", "red", "green", "blue")

scale_col <- data %>% distinct(WGCNA) %>% pull()

data %>%
  group_by(WGCNA) %>%
  mutate(size = size / max(size)) %>%
  mutate(parentTerm = fct_reorder(parentTerm, size, .desc = T)) %>%
  # arrange(size) %>% 
  ungroup() %>%
  # mutate(order = row_number()) %>%
  # mutate(parentTerm = fct_reorder(parentTerm, order, .desc = F)) %>%
  mutate(WGCNA = factor(WGCNA, levels = c(mod_membr_1, mod_membr_2))) %>%
  ggplot(aes(y = parentTerm, x = WGCNA, size = size, color = WGCNA)) +
  geom_point(shape = 1) +
  scale_color_manual('', values = structure(scale_col, names = scale_col) ) +
  scale_fill_manual('', values = structure(scale_col, names = scale_col) ) +
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
    panel.grid.major.x = element_blank()) -> p


ggsave(p, filename = 'WGCNA_MIRS2SEMANTICH.png', path = path, width = 7.3, height = 10, device = png, dpi = 300)
