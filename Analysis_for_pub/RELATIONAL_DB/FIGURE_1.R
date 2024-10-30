
# CREATE A PANEL OF PLOTS FOR FIG. 1
#
rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

# PANEL A: READ DISTRIUTION ----
MTD <- read_tsv('~/Documents/MIRNA_HALIOTIS/METADATA.tsv')

# scale_col <- c("#cd201f", "#FFFC00","#00b489","#31759b")

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

recode_to <- c(`24 HPF`= "24 hpf", `110 HPF` = "110 hpf")

width_col <- 0.8
base_t_text_size <- 6

out %>% 
  dplyr::mutate(hpf = dplyr::recode_factor(hpf, !!!recode_to)) %>%
  group_by(hpf, Length, rnatype) %>%
  summarise(n = sum(n)) %>%
  group_by(Length) %>% mutate(pct = n / sum(n)) %>% 
  # summarise(n = sum(pct)) %>% # scheck
  ggplot(aes(x = Length, y = n, fill = rnatype)) + 
  # geom_vline(xintercept = 25, linetype = "dotted", size = 0.2) +
  geom_col(width = width_col, linewidth = 0.2) +
  # geom_segment(aes(xend = Length, yend = 0,  color = rnatype), linewidth = 1) +
  # facet_grid(hpf ~ ., scales = "free_y", switch = "y") +
  facet_wrap(~ hpf, nrow = 2, scales = "free_y") +
  scale_y_continuous(ylab, labels = scales::number_format(scale = 1/1000000, suffix = "")) +
  scale_x_continuous(xlab, breaks = seq(min(out$Length),max(out$Length), by = 4)) +
  see::scale_fill_pizza(reverse = T) +
  see::scale_color_pizza(reverse = T) -> bottom

bottom <- bottom + 
  guides(fill = guide_legend(title = "", label_size = 1, nrow = 2)) +
  theme_bw(base_family = "GillSans", base_size = base_t_text_size) +
  theme(
    strip.background = element_rect(fill = 'white', color = 'white'),
    strip.text = element_text(color = "black",hjust = 1),
    legend.position = 'top',
    legend.key.width = unit(0.12, "cm"),
    legend.key.height = unit(0.12, "cm"),
    # panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.placement = "inside"
    # axis.text.x = element_text(size = 3),
    )


# bottom <- bottom +
#   annotate("text", y = 16E6, x = 26, 
#     size = 0.75,
#     label = "24 hpf", color = "black", family = "GillSans")

ggsave(bottom, filename = 'FIGURE_1_PANEL_A.png', path = path_out, 
  width = 0.7, height = 1, device = png, dpi = 500)


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

KnownRNAsDB <- KnownRNAsDB %>% mutate(Name = recode(Name, !!!KnownMirID)) %>% 
  mutate(Node_of_origin = ifelse(Name %in% KnownMirID, "Mollusca", Node_of_origin)) %>%
  mutate(DBsource = ifelse(Name %in% KnownMirID, "Mirbase", DBsource)) 

KnownRNAsDB %>%
  filter(Name %in% KnownMirID)

# recode_majorRNA <- KnownRNAsDB %>% filter(Name %in% names(KnownMirID)) %>% pull(MajorRNA, name = Name)
# 
# MIRBASEDB <- data.frame(KnownMirID, 
#   # arm = KnownMirID,
#   MajorRNA = recode_majorRNA, 
#   # sp = "Mollusca",
#   Node_of_origin = "Lophotrochozoa",
#   DBsource = "Mirbase") %>% 
#   as_tibble(rownames = "Name") %>% select(-KnownMirID)

# names(MIRBASEDB) %in% names(KnownRNAsDB)

# KnownRNAsDB <- rbind(KnownRNAsDB, MIRBASEDB)

KnownRNAsDB <- KnownRNAsDB %>% distinct(MajorRNA, Node_of_origin, DBsource)

KnownRNAsDB %>% count(MajorRNA, sort = T)

KnownRNAsDB <- .KnownRNAsDB %>% 
  distinct(MajorRNA, biotype_best_rank) %>% 
  right_join(KnownRNAsDB)
  
KnownRNAsDB %>% count(MajorRNA, sort = T) # less 1 cause Cluster_45662 Cluster_45657 are close 

KnownRNAsDB %>% 
  mutate(Node_of_origin = ifelse(grepl("", Node_of_origin)))

# FRAC. BY KNOWN MIR SPECIES GROUP
DF <- KnownRNAsDB %>% 
  # distinct(MajorRNA,Node_of_origin, DBsource) %>%
  # drop_na(Node_of_origin) %>%
  mutate(Node_of_origin = ifelse(is.na(Node_of_origin), "Novel", Node_of_origin)) %>%
  group_by(DBsource) %>%
  count(Node_of_origin, DBsource, biotype_best_rank, sort = T) %>%
  mutate(pct = n/147)

DF %>% tally(n)


# 1) Plot Node_of_origin by DBsource
DVIZ <- DF %>% drop_na(DBsource) %>%
  group_by(Node_of_origin) %>%
  summarise(pct = sum(pct), n = sum(n)) %>%
  ungroup() %>% arrange(desc(pct)) %>% 
  mutate(Node_of_origin = factor(Node_of_origin, levels=unique(Node_of_origin))) %>%
  mutate(label = paste0(Node_of_origin, " (", n, ")")) 

scale_y_disc <- DVIZ %>% pull(label, name = Node_of_origin)

pr <- DVIZ %>% 
  mutate(facet = "% of known microRNAs") %>%
  ggplot(aes(y = Node_of_origin, x = pct)) +
  facet_wrap(~ facet, scales = "free_x",  switch = "x") +
  geom_col(fill = "grey89") +
  labs(x = "", y = "") +
  scale_y_discrete(labels = scale_y_disc) +
  geom_text(aes(label=label), x = 0.01, hjust=0, size = 2, family = "GillSans") +
  scale_x_continuous(labels = scales::percent) +
  theme_classic(base_family = "GillSans", base_size = 5) +
  theme(
    strip.text = element_text(color = "black",hjust = 1),
    strip.placement = "outside",
    strip.background = element_blank(),
    legend.position = 'none',
    legend.key.width = unit(0.2, "cm"),
    legend.key.height = unit(0.2, "cm"),
    panel.border = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.text.x = element_text(size = 5),
    # axis.line.x = element_blank(),
    panel.grid.major = element_blank())

pr

# 2) Plot Loci and number of known mirbase

DF %>% drop_na(DBsource) 

cols <- c('DBsource', 'biotype_best_rank')

DFlong <- DF %>% 
  mutate(DBsource = ifelse(is.na(DBsource), "Novel", DBsource)) %>%
  group_by(DBsource, biotype_best_rank) %>%
  summarise(pct = sum(pct), n = sum(n)) %>%
  ungroup() %>% arrange(desc(pct)) %>% 
  pivot_longer(cols = all_of(cols), names_to = "x",values_to = "fill")
  
dat_text <- DFlong %>% group_by(x, fill) %>% sample_n(1)

# fill = 

pl <- DF %>% 
  mutate(DBsource = ifelse(is.na(DBsource), "Novel", "Known")) %>%
  ungroup() %>% arrange(desc(pct)) %>% 
  ggplot() + 
  # facet_grid(DBsource ~ ., scales = "free_y", space = "free", switch = "y") +
  geom_col(aes(x = "x", y = n, fill = DBsource), 
    width = width_col, linewidth = 0.2) +
  # facet_grid( DBsource ~., scales = "free", space = "free") +
  labs(x = "", y = "", fill = "") +
  scale_fill_manual(values = c("grey89", "black")) +
  theme_classic(base_family = "GillSans", base_size = 5) +
  theme(
    strip.background = element_blank(),
    strip.text = element_blank(),
    legend.position = 'none',
    legend.key.width = unit(0.2, "cm"),
    legend.key.height = unit(0.2, "cm"),
    panel.border = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    # axis.text.y = element_text(size = 5),
    axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
    panel.grid.major = element_blank())

pl <- pl + annotate("text", x = 1, y = 20, size = 2,
  label = "Novel (47)", angle = 90, color = "white", family = "GillSans", fontface="bold") +
  annotate("text", x = 1, y = 70, size = 2,
    label = "Known (70)", angle = 90, color = "black",
    family = "GillSans", fontface="bold")

# Precisison values 3) -----

# CALL FOR VARIANT AND PRECISION IDENTIFICATION:

pd <- .KnownRNAsDB %>%
  mutate(x = ifelse(is.na(KnownRNAs), "Novel (47)", "Known (70)")) %>%
  mutate(PRECISION = (UniqueReads+MajorRNAReads)/Reads) %>%
  # mutate(PRECISION = MajorRNAReads/Reads) %>%
  ggplot(aes(y = x, x = PRECISION, fill = stat(x))) +
  # facet_grid(x ~ ., scales = "free_y", space = "free", switch = "y") +
  facet_wrap(~ x, nrow = 2, scales = "free_y") +
  labs(y = "", x = "microRNA Precision") +
  scale_fill_viridis_c(option = "C") +
  # xlim(0,1) +
  ggridges::geom_density_ridges_gradient(
    jittered_points = T,
    position = ggridges::position_points_jitter(width = 0.05, height = 0),
    point_shape = '|', point_size = 1, point_alpha = 0.5, alpha = 0.2) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
  theme_classic(base_family = "GillSans", base_size = 7) +
  theme(
    strip.text = element_text(color = "black",hjust = 1),
    strip.background = element_blank(),
    legend.position = 'none',
    legend.key.width = unit(0.2, "cm"),
    legend.key.height = unit(0.2, "cm"),
    panel.border = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.text.x = element_text(size = 5),
    # axis.line.x = element_blank(),
    panel.grid.major = element_blank())

# Bind w/ escoresdf from STRUCVIZ.R
pd <- .KnownRNAsDB %>%
  left_join(escoresdf) %>%
  mutate(x = ifelse(is.na(KnownRNAs), "Novel (47)", "Known (70)")) %>%
  mutate(PRECISION = (UniqueReads+MajorRNAReads)/Reads) %>%
  select(PRECISION, escore, x) %>%
  pivot_longer(-x) %>% 
  mutate(name = ifelse(name == "escore", "Minimum free-energy (MFE)", "microRNA Precision")) %>%
  ggplot(aes(y = x, x = value)) + 
  # facet_grid(x ~ ., scales = "free_y", space = "free", switch = "y") +
  facet_wrap(~ name, nrow = 2, scales = "free_x",  switch = "x") +
  labs(y = "", x = "") +
  # xlim(0,1) +
  ggridges::geom_density_ridges(
    jittered_points = T,
    position = ggridges::position_points_jitter(width = 0.05, height = 0),
    point_shape = '|', point_size = 0.5, point_alpha = 0.5, alpha = 0.2) +
  theme_classic(base_family = "GillSans", base_size = 5) +
  theme(
    strip.text = element_text(color = "black",hjust = 1),
    strip.placement = "outside",
    strip.background = element_blank(),
    legend.position = 'none',
    legend.key.width = unit(0.2, "cm"),
    legend.key.height = unit(0.2, "cm"),
    panel.border = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.text.x = element_text(size = 5),
    # axis.line.x = element_blank(),
    panel.grid.major = element_blank())


pd

p1 <- pl +  plot_spacer() + pr + plot_layout(widths = c(0.3, -0.1, 1))

pss <- pd +  plot_spacer() + p1 + plot_layout(widths = c(0.2, -0.025, 0.9))

ggsave(pss, filename = 'FIGURE_1_PANEL_B.png', path = path_out, 
  width = 4.5, height = 2.5, device = png, dpi = 300)

pss <- pd +  plot_spacer() + pr + plot_layout(nrow = 1, widths = c(2.2, -0.5, 2.5), guides = "collect")

f1 <- bottom +  plot_spacer() + pss + plot_layout(widths = c(1.5, -0.15, 3))


ggsave(f1, filename = 'FIGURE_1.png', path = path_out, 
  width = 3.5, height = 2.5, device = png, dpi = 300)
