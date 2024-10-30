# Calculate distance-distance microRNA loci
# Plot proportion of intragenic/intergenic
# # By 1) Total, 24/110 (development), and acidification
# Merge dataviz to 

# intergenic regions (i.e. independent or canonical transcription)
# Intragenic regions (i.e including introns and exons)


rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

library(GenomicRanges)

dir <-  "~/Documents/MIRNA_HALIOTIS/MIRNA_PUB_2024"

RES.P <- read_rds(paste0(dir, "/SEQUENCES_MERGED_DESEQ_RES_WGCNA.rds")) %>%
  mutate(Contrast = ifelse(sign(log2FoldChange) == -1, sampleB, sampleA)) %>%
  mutate(Contrast = paste0(CONTRAST,"_",Contrast)) %>%
  filter( padj < 0.05  & abs(log2FoldChange) > 1) 

RES.P %>% dplyr::count(CONTRAST)

MIRDB <- read_rds(paste0(dir, "/SEQUENCES_MERGED_DESEQ_RES_WGCNA.rds")) %>% distinct(MirGeneDB_ID, MajorRNA, WGCNA)

MIRDB <- MIRDB %>% mutate(Mir = ifelse(grepl("Cluster", MirGeneDB_ID),"Novel", "Known"))




# Main coords ----

pattern <- "multi_genome.newid.gff3$" # <- to get sequence-region per chromosome

wd_genome <- "/Users/cigom/Documents/MIRNA_HALIOTIS/ENSEMBLE/"

f <- list.files(path = wd_genome, pattern = pattern, full.names = T)

G <- rtracklayer::import(f)

length(G <- G[G$type == "region"])

G <- G %>% as_tibble() %>% dplyr::select(seqnames, start,  end)

G <- G %>% mutate(Chrom = gsub("^JALGQA010000", "", seqnames)) %>% 
  mutate(Chrom = gsub(".1$", "", Chrom)) 

# Mir loci ----
wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out/"

print(.DB <- read_tsv(paste0(wd, "/RNA_LOCATION_DB.tsv")) %>%   filter(SRNAtype == "miR") )


query.chr <- .DB %>% distinct(Chrom) %>% pull()


G  <- G %>% 
  filter(seqnames %in% query.chr) %>%
  mutate(Chrom = gsub("^JALGQA010000", "", seqnames)) %>% 
  mutate(Chrom = gsub(".1$", "", Chrom)) %>% 
  dplyr::rename("Start" = "start", "End" = "end")


other_nc <- c("snRNA", "rRNA", "tRNA", "lncRNA", "srpRNA")
other_intra <- c("exon", "UTR")

.DB <- .DB %>% mutate(Chrom = gsub("^JALGQA010000", "", Chrom)) %>% 
  mutate(Chrom = gsub(".1$", "", Chrom)) %>% 
  # mutate(biotype_best_rank = ifelse(biotype_best_rank %in% other_nc, "Other ncRNA", biotype_best_rank)) %>%
  mutate(biotype_best_rank = ifelse(biotype_best_rank %in% other_nc, "Intergenic", biotype_best_rank)) %>%
  mutate(biotype_best_rank = ifelse(biotype_best_rank %in% other_intra, "Intragenic", biotype_best_rank))


# There is a MajorRNA cluster hosted in two separated regions:

.DB %>% filter(MajorRNA %in% "UCGGUGGGACUUUCGUUCGUCU") 

DB <- .DB %>% 
  select(biotype_best_rank, MajorRNA) %>%
  distinct()


DF1 <- DB %>% dplyr::count(biotype_best_rank) %>% mutate(Frac = n/118, col = "Total")

# .DB %>% distinct(MajorRNA, biotype_best_rank, type, biotype) %>% dplyr::count(biotype_best_rank, type, biotype) %>% view()


"#E7DFD5"
col_recode <- structure(c("#FFC107", "#2196F3","red"), 
  names = c("Intergenic", "Intragenic", "Other ncRNA"))


base_size_theme <- 9

# Merge 

DF2 <- RES.P %>% 
  filter(CONTRAST != "CONTRAST_D") %>%
  select(MajorRNA, Contrast) %>% 
  mutate(Contrast = ifelse(grepl("_Low",Contrast), "Acidification", Contrast)) %>%
  mutate(Contrast = ifelse(grepl("_Low",Contrast), "Acidification", Contrast)) %>%
  mutate(Contrast = ifelse(grepl("CONTRAST_B_Control",Contrast), "Acidification", Contrast)) %>%
  mutate(Contrast = ifelse(grepl("CONTRAST_C",Contrast), "Development", Contrast)) %>%
  distinct() %>% left_join(DB) %>%
  group_by(Contrast) %>%
  dplyr::count(biotype_best_rank) %>% 
  mutate(Frac = n/sum(n), col = Contrast) %>% 
  ungroup() %>% select(-Contrast)

rbind(DF1, DF2) %>%
  ggplot(aes(x = Frac, y = col, fill = biotype_best_rank)) +
  geom_col(position = position_stack(reverse = T)) +
  scale_color_manual("",values = col_recode) +
  scale_fill_manual("",values = col_recode) +
  ylab("") + xlab("Fraction of microRNAs") +
  theme_bw(base_size = base_size_theme, base_family = "GillSans") +
  theme( legend.position = 'top',
    # axis.text.y = element_text(size = 7),
    axis.text.y = element_blank(), 
    axis.ticks.y = element_blank(),
    # legend.key.width = unit(2.5, "mm"),
    # legend.key.height = unit(0.25, "mm"),
    panel.grid = element_blank()) -> p1

data_text <- data.frame(biotype_best_rank = "biotype_best_rank", n = 1, Frac = 0.1, col = c("Total", "Development","Acidification"), 
  label=c("Total (118)","Development (89)","Acidification (16)"))

p1 <- p1 + geom_text(data = data_text, aes(label = label), x = 0.01, hjust=0, size = 3, family = "GillSans")

ggplot2::ggsave(p1, filename = "LOCI_PROPORTION.png", path = dir, width = 2.5, height = 1.5, device = png, dpi = 300)

# p1 <- p1 + ggplot2::scale_x_continuous(position = "top")


# 
# Calculate Distance from the nearest CDS for each miRNA that is located outside of a CDS
# Or  the distance to the closest upstream exon and transcriptional orienta- tion.

.DB %>%
  select(Chrom, MajorRNA, End, biotype_best_rank) %>%
  # left_join(G) %>% mutate(End = End/end)
  arrange(Chrom, End) %>%
  group_by(Chrom) %>%
  mutate(Distance = 1 - min(End)/End) %>%
  mutate(Distance2 = (End-min(End))/1000000) %>%
  filter(Distance2 > 0) %>% count(biotype_best_rank, Distance2) %>% view()
  ggplot(aes(Distance2, fill = biotype_best_rank)) + 
  geom_histogram( position="identity", alpha = 1) +
  ylab("Number of\nmicroRNAs") +
  scale_x_continuous("Distance", 
    labels = scales::number_format(scale = 1, suffix = " Gb")) +
  scale_color_manual("",values = col_recode) +
  scale_fill_manual("",values = col_recode) +
  theme_bw(base_size = base_size_theme, base_family = "GillSans") +
  theme( legend.position = 'none',
    axis.text.y = element_text(size = 7),
    legend.key.width = unit(2.5, "mm"),
    legend.key.height = unit(0.25, "mm"),
    panel.grid = element_blank()) -> p2
p2
# ggsave(p2, filename = 'Nearest_distance.png', path = dir, width = 3.5, height = 2.5, device = png, dpi = 500)

# p+ facet_grid(biotype_best_rank ~.)

# Loci plot -----
MIRDB <- MIRDB %>% right_join(.DB)

p3 <- ggplot() +
  geom_segment(data = G, aes(x = Start, xend = End, yend = Chrom, y = Chrom), linewidth = 0.5, color="#E7DFD5", linetype="dashed") +
  geom_point(data = .DB, 
    aes(y = Chrom, x = Start, color = biotype_best_rank), shape = 108, size = 3) +
  labs(y = "Genomic scaffold") +
  scale_color_manual("",values = col_recode) +
  theme_bw(base_family = "GillSans", base_size = base_size_theme) +
  theme(
    legend.position = "top",
    legend.key.width = unit(2.5, "mm"),
    legend.key.height = unit(0.05, "mm"),
    # axis.text.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 6),
    strip.background = element_rect(fill = 'grey89', color = 'white'),
    # panel.border = element_blank(),
    plot.title = element_text(hjust = 0),
    plot.caption = element_text(hjust = 0),
    panel.grid.minor.x = element_blank(),
    panel.grid.major = element_blank(),
    # panel.grid.major = element_line(linewidth = 0.5, color = "grey89", linetype = "dashed"),
    panel.grid.major.x = element_blank())

p3 <- p3 +   scale_x_continuous("Position", 
  labels = scales::number_format(scale = 1/1000000, suffix = " Gb"), position = "top") 

library(patchwork)

# heights <- c("p1" = 0.05, "spacer" = -0.025, "p2" = 0.1, "spacer" = -0.05, "p3" = 0.5)

# p <- p1/ plot_spacer() / p2 / plot_spacer() / p3 +  plot_layout(heights = heights)

# p <- p1/ p2 / p3 + plot_layout(heights = c(0.07, 0.10, 0.50))

p <- p3/ p2 / p1 + plot_layout(heights = c(0.50, 0.10, 0.07))

ggsave(p, filename = 'LOCI_PLOT2.png', path = dir, width = 2.7, height = 7, device = png, dpi = 500)

# ggplot2::ggsave(p, filename = "SRNA_LOCATION_MIRS_ES.png", path = dir, width = 7, height = 5, device = png, dpi = 300)


# log2fc plot

RES.P %>% 
  left_join(DB) %>%
  filter(CONTRAST %in% c("CONTRAST_C", "CONTRAST_D")) %>%
  # ggplot() +
  ggplot(aes(x = log2FoldChange, y = Contrast, fill = biotype_best_rank)) + 
  facet_grid(CONTRAST ~., scales = "free_y", space = "free_y") +
  ggridges::geom_density_ridges(
    jittered_points = T,
    position = ggridges::position_points_jitter(width = 0.05, height = 0),
    point_shape = '|', point_size = 1, point_alpha = 0.5, alpha = 0.5) +
  scale_color_manual("",values = col_recode) +
  scale_fill_manual("",values = col_recode) +
  theme_bw()


# 2


RES.P %>% 
  left_join(DB) %>%
  filter(CONTRAST %in% c("CONTRAST_C", "CONTRAST_D")) %>%
  # ggplot() +
  ggplot(aes(x = log2FoldChange, y = CONTRAST, fill = biotype_best_rank)) + 
  facet_grid(CONTRAST ~., scales = "free_y", space = "free_y") +
  ggridges::geom_density_ridges(
    jittered_points = T,
    position = ggridges::position_points_jitter(width = 0.05, height = 0),
    point_shape = '|', point_size = 1, point_alpha = 0.5, alpha = 0.5) +
  scale_color_manual("",values = col_recode) +
  scale_fill_manual("",values = col_recode) +
  theme_bw()

  
RES.P %>% 
  left_join(DB) %>%
  filter(CONTRAST %in% "CONTRAST_C") %>% 
  ggplot() +
  geom_point(aes(y = log2FoldChange, x = log2(baseMean), color = biotype_best_rank)) +
  theme_classic() +
  scale_color_manual("",values = col_recode) +
  scale_fill_manual("",values = col_recode)
