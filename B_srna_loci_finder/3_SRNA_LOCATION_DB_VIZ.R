
# RICARDO GOMEZ-REYES  
# AFTER RUN 2_SRNA_LOCATION_DB_PREP.R
# EVALUATE THE SOURCE OF SRNAS
# TRY CIRCOS

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

library(GenomicRanges)

wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out"

print(DB <- read_tsv(paste0(wd, "/RNA_LOCATION_DB.tsv")))

print(DE <- read_tsv(paste0(wd, "/DESEQ_RES.tsv")) )


pattern <- "multi_genome.newid.gff3$" # <- to get sequence-region per chromosome

wd_genome <- "/Users/cigom/Documents/MIRNA_HALIOTIS/ENSEMBLE/"

f <- list.files(path = wd_genome, pattern = pattern, full.names = T)

G <- rtracklayer::import(f)

length(G <- G[G$type == "region"])

G <- G %>% as_tibble() %>% dplyr::select(seqnames, start,  end)

# view(DB)

# 

DB %>% count(biotype_best_rank, sort = T)

DB %>% count(biotype, sort = T) 

DB %>% count(type, sort = T)

DB %>% count(SRNAtype, sort = T)

other_nc <- c("snRNA", "rRNA", "tRNA", "lncRNA", "srpRNA")

# CIRCOS ====

# view(head(DB))

library(gggenes)

other_intra <- c("exon", "UTR")

DB <- DB %>% mutate(Chrom = gsub("^JALGQA010000", "", Chrom)) %>% 
  mutate(Chrom = gsub(".1$", "", Chrom)) %>% 
  mutate(biotype_best_rank = ifelse(biotype_best_rank %in% other_nc, "Other ncRNA", biotype_best_rank)) %>%
  mutate(biotype_best_rank = ifelse(biotype_best_rank %in% other_intra, "Intragenic", biotype_best_rank))


query.chr <- DB %>% filter(SRNAtype == "miR") %>% distinct(Chrom) %>% pull()


# DB %>%
#   # filter(Chrom %in% query.chr ) %>%
#   filter(SRNAtype != "siR") %>%
#   ggplot(aes(xmin = Start, xmax = End, y = Chrom)) +
#   geom_gene_arrow(aes(fill = SRNAtype, color = SRNAtype)) +
#   geom_feature_label(aes(x = Start+(End-Start)/2, y = Chrom, label = biotype_best_rank, forward = Strand)) +
#   theme_genes() 

#  
# DB %>%
#   filter(SRNAtype != "siR") %>%
#   mutate(SRNAtype = ifelse(SRNAtype == "miR", "a", "b")) %>%
#   ggplot(aes(xmin = Start, xmax = End, y = Chrom)) +
#   geom_gene_arrow(aes(fill = biotype_best_rank, color = biotype_best_rank)) +
#   geom_feature_label(aes(x = Start+(End-Start)/2, y = Chrom, label = SRNAtype, forward = Strand), family = "GillSans") + 
#   facet_grid(~ biotype_best_rank) +
#   theme_genes() +
#   scale_x_continuous("Chromosome width (Nucleotides)", 
#     labels = scales::number_format(scale = 1/1000000, suffix = " Gb")) +
#   see::scale_fill_material() +
#   see::scale_color_material() +
#   theme_minimal()

# see::material_colors()

# scales::show_col(see::material_colors())



col_recode <- structure(c("#E7DFD5", "#FFC107", "#2196F3"), names = c("siR", "piR", "miR"))

DB_SIRS <- DB %>% 
  filter(SRNAtype == "siR" & SampleFreq > 3 & UniqueReads > 100) %>%
  filter(between(x = as.numeric(DicerCall), 24, 30))

DB_SIRS %>% count(DicerCall)

G <- G %>% mutate(Chrom = gsub("^JALGQA010000", "", as.character(seqnames))) %>% 
  mutate(Chrom = gsub(".1$", "", Chrom)) %>%
  filter(Chrom %in% query.chr ) %>%
  dplyr::rename("Start" = "start", "End" = "end")


DB <- DB %>% 
  filter(SRNAtype != "siR") %>%
  rbind(DB_SIRS) %>% # Using only siRNA from filtered criteria
  filter(Chrom %in% query.chr ) %>%
  mutate(SRNAtype = factor(SRNAtype, levels = rev(names(col_recode))))

p <- ggplot() +
  geom_segment(data = G, aes(x = Start, xend = End, yend = Chrom, y = Chrom), linewidth = 0.5, color="#E7DFD5", linetype="dashed") +
  geom_gene_arrow(data = DB, aes(xmin = Start, xmax = End, y = Chrom, fill = SRNAtype, color = SRNAtype)) +
  facet_grid(~ biotype_best_rank) +
  theme_genes() +
  scale_x_continuous("Tamaño (Nucleótidos)", 
    labels = scales::number_format(scale = 1/1000000, suffix = " Gb")) +
  labs(y = "Andamios del genoma") +
  scale_color_manual("",values = col_recode) +
  scale_fill_manual("",values = col_recode) +
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


ggplot2::ggsave(p, filename = "SRNA_LOCATION2.png", path = wd, width = 7, height = 7, device = png, dpi = 300)

# ONLY MIRS BY CATEGORY ====


DBMIR <- DB %>% 
  filter(SRNAtype == "miR") %>%
  filter(Chrom %in% query.chr ) 

DB_FAM <- DE %>% distinct(Name, Family)


DE <- DE %>% filter( padj < 0.05 & abs(log2FoldChange) > 1) 

WHICH_DE <- DE %>% filter(CONTRAST %in% c("CONTRAST_A", "CONTRAST_B")) %>% distinct(Family) %>% pull()

DBMIR <- DBMIR %>%
  left_join(DB_FAM, by = "Name") %>% 
  mutate(facet = ifelse(grepl("^Cluster", Family), "A)", "B)")) 

DBMIR %>% dplyr::count(facet)

key_recode <- structure(c("Intergénico", "Intragénico", "Otros ncRNA"), 
  names = c("Intergenic", "Intragenic", "Other ncRNA"))

DBMIR <- DBMIR %>%
  dplyr::mutate(biotype_best_rank = dplyr::recode_factor(biotype_best_rank, !!!key_recode)) 

# col_recode <- structure(c("#f44336", "#FFC107", "#2196F3"), names = c("Other ncRNA", "Intergenic", "Intragenic"))

col_recode <- structure(c("#f44336", "#FFC107", "#2196F3"), names = c("Otros ncRNA", "Intergénico", "Intragénico"))

p <- ggplot() +
  geom_segment(data = G, aes(x = Start, xend = End, yend = Chrom, y = Chrom), linewidth = 0.5, color="#E7DFD5", linetype="dashed") +
  geom_gene_arrow(data = DBMIR, aes(xmin = Start, xmax = End, y = Chrom, fill = biotype_best_rank, color = biotype_best_rank)) +
  # facet_grid(~ facet) +
  theme_genes() +
  scale_x_continuous("Tamaño (Nucleótidos)", 
    labels = scales::number_format(scale = 1/1000000, suffix = " Gb")) +
  labs(y = "Andamios del genoma") +
  scale_color_manual("",values = col_recode) +
  scale_fill_manual("",values = col_recode) +
  theme_bw(base_family = "GillSans", base_size = 10) +
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

ggplot2::ggsave(p, filename = "SRNA_LOCATION_MIRS_ES.png", path = wd, width = 7, height = 5, device = png, dpi = 300)

write_tsv(DBMIR, file = paste0(wd, "MIRLOCI_BIOTYPES_best_rank.tsv"))

# and by barplot

DBMIR %>% 
  group_by(facet, biotype_best_rank) %>%
  # dplyr::count()
  summarise(Reads = sum(Reads), n = n()) %>% 
  ggplot(aes(x = biotype_best_rank, y = n, fill = facet)) +
  geom_col(position = position_identity()) +
  scale_color_manual("",values = col_recode) +
  # scale_fill_manual("",values = col_recode) +
  scale_y_continuous("Freq. of reads",labels = scales::comma) +
  theme_bw(base_family = "GillSans", base_size = 11) +
  theme(strip.background = element_rect(fill = 'grey89', color = 'white'),
    legend.position = 'top',
    panel.border = element_blank(),
    panel.grid.minor = element_blank())

# IN RESPONSE TO ACIDIFICATION

WHICH_DE <- DE %>% 
  filter(CONTRAST %in% c("CONTRAST_A", "CONTRAST_B")) %>% 
  dplyr::distinct(CONTRAST, Family) %>%
  pull(Family)


DBMIR %>% 
  filter(Family %in% WHICH_DE) %>%
  group_by(biotype_best_rank) %>%
  # dplyr::count()
  summarise(Reads = sum(Reads), n = n()) 

# IN RESPONSE TO DEV
# 110 HPF
WHICH_DE <- DE %>% 
  filter(log2FoldChange > 0) %>%
  filter(CONTRAST %in% c("CONTRAST_C", "CONTRAST_D")) %>% 
  dplyr::distinct(CONTRAST, Family) %>%
  dplyr::count(Family, sort = T) %>% 
  filter(n > 1) %>% pull(Family)

DBMIR %>% 
  filter(Family %in% WHICH_DE) %>%
  group_by(biotype_best_rank) %>%
  # dplyr::count()
  summarise(Reads = sum(Reads), n = n()) 

# p

df <- DB %>% 
  mutate(width = nchar(MajorRNA)) %>%
  mutate(biotype_best_rank = ifelse(biotype_best_rank %in% other_nc, "Other ncRNA", biotype_best_rank)) %>%
  group_by(SRNAtype, biotype_best_rank, width) %>% 
  summarise(Reads = sum(Reads), n = n()) %>%
  group_by(SRNAtype) %>%
  mutate(reads_pct = Reads/sum(Reads), n_freq = n/sum(n))

df %>% summarise(sum(Reads), sum(n))

df %>% summarise(sum(reads_pct), sum(n_freq))

df %>% group_by(SRNAtype, biotype_best_rank) %>% 
  summarise(values_from = sum(n)) %>% 
  pivot_wider(names_from = SRNAtype, values_from = values_from) %>% view()

df %>% group_by(SRNAtype, biotype_best_rank) %>% 
  summarise(values_from = sum(Reads)) %>%  
  pivot_wider(names_from = SRNAtype, values_from = values_from) %>% view()

# (doi.org/10.1186/1471-2164-11-533): Roughly half of known miRNA genes are located within previously annotated protein-coding regions ("intragenic miRNAs"). A high-confidence set of predicted mRNA targets of intragenic miRNAs also shared many of these features with the host genes. Approximately 20% of intragenic miRNAs were predicted to target their host mRNA transcript.

df %>% filter(SRNAtype == "miR") %>% group_by(biotype_best_rank) %>% summarise(sum(n_freq))

p <- df %>%
  ggplot(aes(x = biotype_best_rank, y = Reads, fill = SRNAtype)) +
  geom_col(position = position_dodge()) +
  scale_color_manual("",values = col_recode) +
  scale_fill_manual("",values = col_recode) +
  scale_y_continuous("Freq. of reads",labels = scales::comma) +
  theme_bw(base_family = "GillSans", base_size = 11) +
  theme(strip.background = element_rect(fill = 'grey89', color = 'white'),
    legend.position = 'top',
    panel.border = element_blank(),
    panel.grid.minor = element_blank())

# ggplot2::ggsave(p, filename = "SRNA_LOCATION.png", path = wd, width = 5, height = 5, device = png, dpi = 300)

df %>%
  ggplot(aes(x = width, y = Reads, fill = SRNAtype)) +
  geom_col() +
  labs(x = "Read length (nt)") +
  see::scale_fill_material() +
  scale_y_continuous("Number of reads",labels = scales::comma) +
  theme_bw(base_family = "GillSans", base_size = 11) +
  theme(strip.background = element_rect(fill = 'grey89', color = 'white'),
    legend.position = 'top',
    panel.border = element_blank(),
    panel.grid.minor = element_blank()) +
  facet_wrap(~ biotype_best_rank, scales = "free_y", nrow = 1)
# facet_grid(~ SRNAtype, scales = "free", space = "free")


DB %>% filter(SRNAtype == "miR") %>% dplyr::count(MIRNA)  

DB %>% filter(SRNAtype == "miR" & MIRNA == "N") %>% group_by(biotype_best_rank) %>% dplyr::count(MIRNA)

DB %>% drop_na(KnownRNAs) %>% filter(SRNAtype == "miR" & MIRNA == "Y") %>% group_by(biotype_best_rank) %>% dplyr::count(MIRNA)

DB %>% filter(is.na(KnownRNAs)) %>% filter(SRNAtype == "miR" & MIRNA == "Y") %>% group_by(biotype_best_rank) %>% dplyr::count(MIRNA)


DB %>% 
  filter(is.na(KnownRNAs)) %>% filter(SRNAtype == "miR" & MIRNA == "Y") %>%
  group_by(SRNAtype, biotype_best_rank) %>% 
  summarise(values_from = sum(Reads)) %>%  
  pivot_wider(names_from = SRNAtype, values_from = values_from) 


library(ggbio)

# DB <- # genomic range

ggplot() +
  ggbio::layout_circle(DB, aes(fill = SRNAtype, y = Reads), geom = "rect") +
  theme(legend.position = "top")


# PLOT BY TRANPOSON TYPE? ====
# PLOT of piRNA cluster expression across developmental stages
# SEE FIG 7 ADN 8 FROM 10.1080/15476286.2017.1349048. .
# Praher D, Zimmermann B, Genikhovich G, Columbus-Shenkar Y, Modepalli V, Aharoni R, Moran Y, Technau U. Characterization of the piRNA pathway during development of the sea anemone Nematostella vectensis. RNA Biol. 2017 Dec 2;14(12):1727-1741. doi: 10.1080/15476286.2017.1349048. Epub 2017 Sep 13. PMID: 28783426; PMCID: PMC5731801.

#
DB %>%
  filter(!grepl("exon|Intergenic|Intragenic|UTR", biotype)) %>%
  mutate(width = nchar(MajorRNA)) %>%
  # mutate(biotype_best_rank = ifelse(biotype_best_rank %in% other_nc, "Other ncRNA", biotype_best_rank)) %>%
  group_by(SRNAtype,biotype, type) %>%
  summarise(Reads = sum(Reads), n = n())  %>% view()


# ANALYSIS OF CANONICAL OR ISOFORM (NEILSEN ET AL 2012) ?

DB %>% 
  drop_na(KnownRNAs) %>%
  ggplot(aes(MajorRNAReads/Reads)) +
  geom_histogram()

#

DB %>% 
  # drop_na(KnownRNAs) %>%
  ggplot(aes(FracTop, MajorRNAReads/Reads, color = Strand)) +
  facet_grid(~ SRNAtype) +
  geom_point()


# ANY CUTOFF?

DB %>% count(SRNAtype)

DB %>% 
  # filter()
  group_by(SRNAtype) %>%
  count(type, sort = T) %>%
  ggplot()
  
  
# 

.out %>% 
  # filter(Name %in% which_mirs) %>%
  count(Name, sort = T) 

# gggenes::example_dummies %>%
#   ggplot(aes(xmin = start, xmax = end, y = molecule,fill = gene)) +
#   geom_gene_arrow()

which_cluster <- "Cluster_2825"

cluster_bind <- assembly[assembly$ID == which_cluster] %>% as_tibble() %>% 
  rename("Name"="ID") %>%
  mutate(biotype = type, gene_id = DicerCall, transcript_id = MIRNA) %>% 
  select(any_of(names(.out)))


.out %>% 
  filter(Name == which_cluster) %>% 
  # filter(!type %in% c("mRNA", "gene", "ncRNA")) %>%
  rbind(cluster_bind) %>%
  ggplot(aes(xmin = start, xmax = end, y = 1, forward = T, fill = biotype)) +
  geom_gene_arrow() +
  facet_grid(type ~.)


  
  

