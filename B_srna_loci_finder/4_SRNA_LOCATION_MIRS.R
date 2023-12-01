# RICARDO GOMEZ-REYES, 2023
# SRNA_LOCATION_MIRS.R
# REPLACE DATA FROM WGCNA FOR ONLY MIRS (6_WGCNA_MIRS.R)

wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out/"

.DB <- read_tsv(paste0(wd, "/RNA_LOCATION_DB.tsv")) %>% dplyr::filter(SRNAtype == "miR") %>%
  dplyr::select(-WGCNA)

WGCNA <- read_rds(paste0(wd, "/WGCNA_MIRS.rds"))

# read_tsv(paste0(wd, "/SRNA2MIRGENEDB.tsv")) %>%
  
# path <- "/Users/cigom/Documents/MIRNA_HALIOTIS/"

FAMILYDB <- read_tsv(paste0(wd, "/DESEQ_RES.tsv")) %>% distinct(Name, Family)

module_members_1 <- c("grey", "black", "yellow", "pink", "turquoise")

module_members_2 <- c("brown", "red", "green", "blue")

DB <- WGCNA %>% right_join(FAMILYDB) %>% left_join(.DB, by = "Name")
  
other_nc <- c("snRNA", "rRNA", "tRNA", "lncRNA", "srpRNA")

other_intra <- c("exon", "UTR")


DB <- DB %>% mutate(Chrom = gsub("^JALGQA010000", "", Chrom)) %>% 
  mutate(Chrom = gsub(".1$", "", Chrom)) %>% 
  mutate(biotype_best_rank = ifelse(biotype_best_rank %in% other_nc, "Other ncRNA", biotype_best_rank)) %>%
  mutate(biotype_best_rank = ifelse(biotype_best_rank %in% other_intra, "Intragenic", biotype_best_rank))

write_rds(DB, file = paste0(path, "RNA_LOCATION_MIR_DB.rds"))

str(query.chr <- DB %>% distinct(Chrom) %>% pull())

# VIZ

pattern <- "multi_genome.newid.gff3$" # <- to get sequence-region per chromosome

wd_genome <- "/Users/cigom/Documents/MIRNA_HALIOTIS/ENSEMBLE/"

f <- list.files(path = wd_genome, pattern = pattern, full.names = T)

.G <- rtracklayer::import(f)

length(G <- .G[.G$type == "region"])

G <- G %>% as_tibble() %>% dplyr::select(seqnames, start,  end)

G <- G %>% mutate(Chrom = gsub("^JALGQA010000", "", as.character(seqnames))) %>% 
  mutate(Chrom = gsub(".1$", "", Chrom)) %>%
  filter(Chrom %in% query.chr ) %>%
  dplyr::rename("Start" = "start", "End" = "end")

col_recode <- structure(c("#f44336", "#FFC107", "#2196F3"), names = c("Other ncRNA", "Intergenic", "Intragenic"))

scale_col <- DB %>% pull(WGCNA)



library(gggenes)


ggplot() +
  geom_segment(data = G, aes(x = Start, xend = End, yend = Chrom, y = Chrom), linewidth = 0.5, color="#E7DFD5", linetype="dashed") +
  geom_gene_arrow(data = DB, aes(xmin = Start, xmax = End, y = Chrom, fill = WGCNA, color = WGCNA)) +
  facet_grid(~ WGCNA) +
  theme_genes() +
  scale_x_continuous("Tamaño (Nucleótidos)", 
    labels = scales::number_format(scale = 1/1000000, suffix = " Gb")) +
  labs(y = "Andamios del genoma") +
  scale_color_manual('', values = structure(scale_col, names = scale_col) ) +
  scale_fill_manual('', values = structure(scale_col, names = scale_col) ) +
  theme_bw(base_family = "GillSans", base_size = 12) +
  theme(
    legend.position = "top",
    # axis.text.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 7),
    strip.background = element_rect(fill = 'grey89', color = 'white'),
    panel.border = element_blank(),
    plot.title = element_text(hjust = 0),
    plot.caption = element_text(hjust = 0),
    panel.grid.minor.x = element_blank(),
    panel.grid.major = element_blank(),
    # panel.grid.major = element_line(linewidth = 0.5, color = "grey89", linetype = "dashed"),
    panel.grid.major.x = element_blank())

