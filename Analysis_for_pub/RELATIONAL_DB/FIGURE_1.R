
# CREATE A PANEL OF PLOTS FOR FIG. 1
#
rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

# PANEL A: READ DISTRIUTION ----
MTD <- read_tsv('~/Documents/MIRNA_HALIOTIS/METADATA.tsv')

scale_col <- c("#cd201f", "#FFFC00","#00b489","#31759b")

recode_to <- c(`Control` = "pH 8.0", `Low` = "pH 7.6")

MTD <- MTD %>%
  dplyr::rename("sample_id" = "LIBRARY_ID") %>%
  dplyr::mutate(hpf = paste0(hpf, " HPF")) %>%
  dplyr::mutate(hpf = factor(hpf, levels = c("24 HPF", "110 HPF"))) %>%
  dplyr::mutate(pH = dplyr::recode_factor(pH, !!!recode_to))

path <- '~/Documents/MIRNA_HALIOTIS/PROFILING_BY_READ_LENGTH/'

path_out <- "~/Documents/MIRNA_HALIOTIS/MIRNA_PUB_2024/"

recode_to <- c(`miRNAs` = "miRNAs", `Unknown`= "Unknown",`rRNA` = "rRNA", `tRNA` = "tRNA", `Artifacts` = "Artifacts")

out <- read_rds(paste0(path, "prof_by_read_length_summary.rds")) %>% 
  left_join(MTD) %>%
  mutate(first_nuc = recode(first_nuc, `T` = "U")) %>%
  dplyr::mutate(rnatype = dplyr::recode_factor(rnatype, !!!recode_to))


# SANITY CHECK

out %>% group_by(sample_id) %>% tally(n)


# 1) Bottom -----

xlab <-  "Read Length (Nucleotides)" # "Read Length (nt)"
ylab <- "Reads (Millions)"

recode_to <- c(`24 HPF`= "B) 24 HPF", `110 HPF` = "C) 110 HPF")

width_col <- 0.8
base_t_text_size <- 5

out %>% 
  dplyr::mutate(hpf = dplyr::recode_factor(hpf, !!!recode_to)) %>%
  group_by(hpf, Length, rnatype) %>%
  summarise(n = sum(n)) %>%
  group_by(Length) %>% mutate(pct = n / sum(n)) %>% 
  # summarise(n = sum(pct)) %>% # scheck
  ggplot(aes(x = Length, y = n, fill = rnatype)) + 
  geom_col(width = width_col, color = "black", linewidth = 0.2) +
  # geom_segment(aes(xend = Length, yend = 0,  color = rnatype), linewidth = 1) +
  facet_grid(hpf ~ ., scales = "free_y", switch = "y") +
  scale_y_continuous(ylab, labels = scales::number_format(scale = 1/1000000, suffix = "")) +
  scale_x_continuous(xlab, breaks = seq(min(out$Length),max(out$Length), by = 4)) +
  see::scale_fill_pizza(reverse = T) +
  see::scale_color_pizza(reverse = T) -> bottom

bottom <- bottom + theme_classic(base_family = "GillSans", base_size = base_t_text_size) +
  theme(strip.background = element_rect(fill = 'grey89', color = 'white'),
    legend.position = 'none',
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.placement = "outside",
    axis.text.x = element_text(size = base_t_text_size))

# 2) top -----

ylab <- "Frac. of reads"
xlab <- ""


top <- out %>% 
  group_by(Length, rnatype) %>% summarise(n = sum(n)) %>%
  group_by(Length) %>% mutate(pct = n / sum(n)) %>%
  mutate(sample_id = "A) Global") %>%
  ggplot(aes(x = Length, y = pct, fill = rnatype)) + 
  geom_col(width = width_col, color = "black", linewidth = 0.2) +
  # geom_segment(aes(xend = Length, yend = 0, color = rnatype), linewidth = 1) +
  facet_grid(sample_id ~ ., scales = "free_y", switch = "y") +
  scale_y_continuous(ylab, labels = scales::percent) +
  scale_x_continuous(xlab, breaks = seq(min(out$Length),max(out$Length), by = 2)) +
  see::scale_fill_pizza(reverse = T) + 
  see::scale_color_pizza(reverse = T) +
  guides(fill = guide_legend(title = "", label_size = 3)) +
  theme_classic(base_family = "GillSans", base_size = base_t_text_size) +
  theme(strip.background = element_rect(fill = 'grey86', color = 'white'),
    legend.position = 'top',
    legend.key.width = unit(0.2, "cm"),
    legend.key.height = unit(0.2, "cm"),
    panel.border = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    strip.placement = "outside",
    panel.grid.major = element_blank())


library(patchwork)


ps <- top / plot_spacer() / bottom + plot_layout(heights = c(0.5, -0.35, 1))

# ggsave(ps, filename = 'FIGURE_2.png', path = path_out, width = 2, height = 1.7, device = png, dpi = 300)

# PANEL B): Loci plot -----

wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/MIRGENEDB_20230314/"

print(MIRGENEDB <- read_tsv(paste0(wd, "/MIRGENEDB_2.1.tsv")))

MIRGENEDB <- MIRGENEDB %>% 
  select(MirGeneDB_ID, `Node_of_origin_(family)` ) %>%
  distinct() %>% mutate(DBsource = "MIRGENEDB_2.1") %>%
  rename("KnownRNAs" = "MirGeneDB_ID", "Node_of_origin" = "Node_of_origin_(family)")

dir <- "~/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out/"

print(.KnownRNAsDB <- read_rds(paste0(dir, "RNA_LOCATION_MIR_DB.rds")))

# KnownRNAsDB %>% distinct(Name, biotype_best_rank) %>% count(biotype_best_rank)

KnownRNAsDB <- .KnownRNAsDB %>%
  mutate(KnownRNAs = strsplit(KnownRNAs, ";")) %>%
  dplyr::select(Name, KnownRNAs, MajorRNA) %>%
  unnest(KnownRNAs) %>%
  filter(!grepl("^piR-",KnownRNAs)) %>% # <- ALSO EXCLUDE OVERLAPED MIRS/PIRS
  mutate(KnownRNAs = gsub("_[3-5]p$","", KnownRNAs)) %>%
  distinct(Name, MajorRNA, KnownRNAs) %>%
  # separate(KnownRNAs, into = c("KnownMirID", "arm"), sep = "_") %>% 
  # mutate(sp = sapply(strsplit(KnownMirID, "-"), `[`, 1)) %>%
  left_join(MIRGENEDB, by = "KnownRNAs")


KnownRNAsDB <- KnownRNAsDB %>% distinct(Name, MajorRNA, Node_of_origin, DBsource)

# RECODE 7 MIRS FAMILY FROM MIRBASE/MOLLUSC DB (FROM 4_MIR_PAIRWISEALIGNMENT.R)

KnownMirID <- c(
  `Cluster_39774` = "MIR-242", 
  `Cluster_17642` = "MIR-1989b-5p", 
  `Cluster_38139` = "MIR-745", 
  `Cluster_45860` = "MIR-10492a-5p", 
  `Cluster_47716` = "MIR-277b-3p", 
  `Cluster_37147` = "MIR-184-3p", 
  `Cluster_27861` = "MIR-252-5p")

recode_majorRNA <- KnownRNAsDB %>% filter(Name %in% names(recode_mir)) %>% pull(MajorRNA, name = Name)

MIRBASEDB <- data.frame(KnownMirID, 
  # arm = KnownMirID,
  MajorRNA = recode_majorRNA, 
  # sp = "Mollusca",
  Node_of_origin = "Lophotrochozoa",
  DBsource = "Mirbase") %>% 
  as_tibble(rownames = "Name") %>% select(-KnownMirID)


KnownRNAsDB <- rbind(KnownRNAsDB, MIRBASEDB)

KnownRNAsDB <- .KnownRNAsDB %>% 
  distinct(Name, biotype_best_rank) %>% 
  right_join(KnownRNAsDB)
  
KnownRNAsDB %>% distinct(MajorRNA)

# FRAC. BY KNOWN MIR SPECIES GROUP
DF <- KnownRNAsDB %>% 
  # distinct(MajorRNA,Node_of_origin, DBsource) %>%
  # drop_na(Node_of_origin) %>%
  mutate(Node_of_origin = ifelse(is.na(Node_of_origin), "Novel", Node_of_origin)) %>%
  group_by(DBsource) %>%
  count(Node_of_origin, DBsource, biotype_best_rank, sort = T) %>%
  mutate(pct = n/147)


DF %>% ungroup() %>%
  distinct(Node_of_origin, pct)

DF %>% 
  ungroup() %>%
  mutate(factor = biotype_best_rank) %>%
  # pivot_longer(cols = c('DBsource','Node_of_origin', 'biotype_best_rank')) %>%
  ggplot(aes(x = "name", y = pct, fill = factor)) + 
  # facet_grid(~ value) +
  geom_col() +
  geom_text(aes(label = factor), stat = 'identity', angle = 0, 
    position = 'identity', check_overlap = TRUE)

# KnownRNAsDB <- KnownRNAsDB %>%
#   mutate(MirGeneDB_ID = ifelse(is.na(MirGeneDB_ID), MajorRNA, MirGeneDB_ID)) %>%
#   dplyr::mutate(MirGeneDB_ID = dplyr::recode(MirGeneDB_ID, !!!recode_mir)) %>% 
#   drop_na()
