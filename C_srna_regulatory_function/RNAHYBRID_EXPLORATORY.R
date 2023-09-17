

# READ RNAHYBRYD CSV ====


# pattern = "mir_vs_utr_rmdup_RNAhybrid.out.tsv$"

pattern <- "mir_vs_utr_rmdup_RNAhybrid.out.psig.tsv"


f <- list.files(path = wd, pattern = pattern, full.names = T)

# df <- read_delim(f, col_names = F, delim = "\t")
# colNames <- c("target", "query", "mfe", "pval", "pos", "lenT", "lenQ")
# colnames(df) <- colNames

# df <- read_tsv(f)

# str(targetsids <- df %>% distinct(target) %>% pull(target))

# write_lines(targetsids, file = paste0(gsub(".tsv", "", f), ".ids"))

# df <- df %>% filter(pval < 0.05)

# write_tsv(df, file = paste0(gsub(".tsv", "", f), ".psig.tsv"))

# AFTER CLEAN pval < 0.05, JUST LOAD .out.psig.tsv file:

df <- read_tsv(f)


# Q: how many unique binding sites? =====

nrow(df %>% distinct(target)) #  37,663

nrow(df %>% distinct(query)) # 261 from 393 clusters from mir.fasta (I.E. ONLY MATURE AND STAR SEQS)

# mir.fasta: 
# This is a FASTA formatted file containing hairpin, mature miRNA, and miRNA* sequences derived from ShortStack's identification of MIRNA loci.

# Q: how many were 5' and 3' UTR? ====

str(df %>% distinct(target) %>% pull()) #  37 663

str(three_prime_utr <- sapply(strsplit(x, " "), `[`, 1)) # three_prime_utr ids 

howtargets <- df %>% 
  distinct(target) %>% # comment this if wish to know multi-binding sites 
  mutate(UTR = ifelse(target %in% three_prime_utr, "three_prime_utr", "five_prime_utr"))

howtargets %>% count(UTR)


# Explanation: 
# 8,969 three_prime_utr ids with predicted target site for miR binding 
# 28, 694 five_prime_utr ids with predicted target site for miR binding


howmirs <- df %>% 
  # distinct(target, query) %>%
  mutate(MIRNA = "hairpin") %>%
  mutate(MIRNA = ifelse(grepl("mature", query), "mature", MIRNA)) %>%
  mutate(MIRNA = ifelse(grepl("star", query), "star", MIRNA)) 


howmirs %>%
  mutate(UTR = ifelse(target %in% three_prime_utr, 
    "three_prime_utr", "five_prime_utr")) %>%
  distinct(target ,.keep_all = T) %>% # comment this if wish to know multi-binding sites
  group_by(UTR) %>%
  count(MIRNA)

# By energy?

p1 <- howmirs %>%
  mutate(UTR = ifelse(target %in% three_prime_utr, 
    "three_prime_utr", "five_prime_utr")) %>%
  ggplot(aes(y = -log10(pval), x = mfe, color = UTR)) + 
  # facet_grid(~ MIRNA) +
  geom_point(shape = 1, alpha = 0.5) +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color="grey") +
  geom_vline(xintercept = -35, linetype = "dashed", color="grey") +
  theme_bw(base_family = "GillSans", base_size = 12) +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  theme(strip.background = element_rect(fill = 'white', color = 'white'),
    legend.position = 'none',
    panel.border = element_blank(),
    panel.grid.minor = element_blank()) +
  labs(x = "mfe (kcal/mol)")

p2 <- howmirs %>% 
  mutate(UTR = ifelse(target %in% three_prime_utr, 
  "three_prime_utr", "five_prime_utr")) %>%
  ggplot(aes(mfe, color = UTR)) + 
  stat_ecdf(linewidth = 2) + #  linetype = "dashed"
  geom_rug() +
  scale_y_continuous("f(mfe)",labels = scales::comma) +
  theme_bw(base_family = "GillSans", base_size = 12) +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  theme(strip.background = element_rect(fill = 'white', color = 'white'),
    legend.position = 'right',
    panel.border = element_blank(),
    panel.grid.minor = element_blank()) 
 
library(patchwork)
 

ps <- p1 + p2 + plot_layout(widths = c(5, 6))

ggsave(ps, filename = 'RNAhybrid.png', path = wd, width = 10, height = 4, 
  device = png, dpi = 300)


# Q: How many mirs per target:? =====

mir_sites <- howmirs %>%
  mutate(UTR = ifelse(target %in% three_prime_utr, 
    "three_prime_utr", "five_prime_utr")) %>%
  group_by(target, UTR) %>%
  count(MIRNA, sort = T) %>%
  group_by(target)

# Santity check

identical(sum(mir_sites$n), nrow(df)) # 90918

mir_sites %>% 
  # group_by(UTR, MIRNA, n) %>%
  # tally(n) %>% view() 
  group_by(UTR, MIRNA) %>%
  rstatix::get_summary_stats(n, type = "common")

p1 <- mir_sites %>% 
  group_by(n) %>%
  dplyr::count(UTR, MIRNA) %>%
  ggplot(aes(y = nn, x = n, color = UTR, fill = UTR)) +
  geom_col() +
  coord_flip() +
  facet_grid(~ MIRNA) +
  theme_bw(base_size = 16, base_family = "GillSans") +
  labs(x = "binding miRs per target") +
  scale_y_continuous("Target gene", labels = scales::comma) +
  theme(legend.position = "top") +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  scale_fill_manual(values = c("#00AFBB", "#E7B800")) 

p2 <- mir_sites %>% 
  group_by(n) %>%
  dplyr::count(UTR, MIRNA) %>% 
  # pivot_wider(names_from = MIRNA, values_from = nn)
  ggplot(aes(nn, color = UTR)) +
  facet_grid(~ MIRNA) +
  stat_ecdf(linewidth = 2, alpha = 0.5) +
  theme_bw(base_size = 16, base_family = "GillSans") +
  labs(x = "Target gene") +
  scale_y_continuous("Freq.", labels = scales::comma) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) 

library(patchwork)

p1/p2

# Q: how many chromosomes were binding sites for mirs? ====
# 
# chromosome_df <- df %>% 
#   separate(col = target, into = c("seqnames", "ranges"), sep = "_") %>%
#   count(seqnames)

chromosome_df <- howmirs %>%
  mutate(UTR = ifelse(target %in% three_prime_utr, 
    "three_prime_utr", "five_prime_utr")) %>%
  separate(col = target, into = c("seqnames", "ranges"), sep = "_") %>%
  group_by(MIRNA, UTR) %>%
  count(seqnames) %>%
  ungroup()

length(unique(chromosome_df$seqnames)) # 243

chromosome_df %>% pull(n) %>% sum() # 90918

chromosome_df <- .gr[.gr$type == "region"] %>% 
  as_tibble() %>% select(seqnames, width) %>%
  left_join(chromosome_df) %>% drop_na()

# How many 

# mith? zero <--

chromosome_df %>% filter(seqnames %in% "JALGQA010000616.1")

p3 <- chromosome_df %>%
  ggplot(aes(width, n, color = UTR)) +
  facet_grid(~ MIRNA) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", linetype="dashed", size = 0.5, alpha=0.5,
    se = F, na.rm = TRUE) +
  ggpubr::stat_cor(method = "spearman", cor.coef.name = "R", p.accuracy = 0.001) +
  theme_bw(base_size = 16, base_family = "GillSans") +
  theme(legend.position = "top") +
  scale_y_continuous("Binding sites (UTR targets)", labels = scales::comma) +
  scale_x_continuous("Chromosome width (Nucleotides)", 
    labels = scales::number_format(scale = 1/1000000, suffix = " Gb")) +
  theme(legend.position = "top") +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  scale_fill_manual(values = c("#00AFBB", "#E7B800")) 

library(patchwork)

# p1/p2/p3

# Join features to miRNA:mRNA target prediction ====

# Q: How many target annotation from features? ====

str(targets <- df %>% distinct(target) %>% pull(target)) # 37 663

nrow(utr_source %>% filter(target %in% targets)) # 8 970

llist <- function(x) {
  x <- paste(x, sep = ';', collapse = ';')
}

# x <- "Cluster_286.mature::JALGQA010000001.1:13045384-13045405(+)"
# sapply(strsplit(x, "[.]"), `[`, 1)

# sapply(strsplit(x, "::"), `[`, 1) 


df_drop <- df %>%
  mutate(query = sapply(strsplit(query, "::"), `[`, 1) ) %>%
  group_by(target) %>%
  summarise(across(query, .fns = llist), .groups = "drop_last") %>% 
  ungroup()

df_drop <- utr_source %>% right_join(df_drop)
  

saveRDS(df_drop, file = paste0(wd, "/mir_vs_utr_rmdup_RNAhybrid.out.psig.gene_features.rds"))


# KARYOPLOTER (OMIT) ====
# seqnames      ranges strand

df <- df %>% 
  separate(col = target, into = c("seqnames", "ranges"), sep = "_") %>%
  separate(col = ranges, into = c("ranges", "strand"), sep = ":") %>%
  separate(col = ranges, into = c("start", "end"), sep = "-")

df <- df %>% mutate_at(c("start", "end", "mfe", "pos"), as.numeric)

ranges_ <- IRanges(start = df$start, end = df$end) # Those UTR ranges

seqnames_ <- df$seqnames

grr <- GRanges(Rle(seqnames_), 
  ranges = ranges_, strand = df$strand, 
  mfe = df$mfe, pval = df$pval, query = df$query)


# Q: how many targets were binding sites for mirs? ====
# Ex.

howmirs %>% filter(target %in% "JALGQA010000017.1_22665876-22665941:-") 
howmirs %>% filter(query %in% "Cluster_47717.mature::JALGQA010000024.1:8020797-8020820(+)") %>% view()
# Contrast w/ mir_sites

target_df <- howmirs %>%
  mutate(UTR = ifelse(target %in% three_prime_utr, 
    "three_prime_utr", "five_prime_utr")) %>%
  group_by(MIRNA, UTR) %>%
  count(query, sort = T) %>%
  group_by(query)

p1 <- target_df %>%
  ggplot(aes(n, color = UTR, fill = UTR)) +
  facet_grid(~ MIRNA) +
  geom_histogram() +
  theme_bw(base_size = 12, base_family = "GillSans") +
  scale_x_continuous("Number of target", labels = scales::comma) +
  scale_y_continuous("Number of miRs", labels = scales::comma) +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  scale_fill_manual(values = c("#00AFBB", "#E7B800")) +
  theme(strip.background = element_rect(fill = 'white', color = 'white'),
    legend.position = 'top',
    panel.border = element_blank(),
    panel.grid.minor = element_blank()) 

p2 <- mir_sites %>% 
  ggplot(aes(n, color = UTR, fill = UTR)) +
  facet_grid(~ MIRNA) +
  geom_histogram() +
  theme_bw(base_size = 12, base_family = "GillSans") +
  labs(x = "miRs") +
  scale_x_continuous("Number of miRs", labels = scales::comma) +
  scale_y_continuous("Number of targets", labels = scales::comma) +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  scale_fill_manual(values = c("#00AFBB", "#E7B800")) +
  theme(strip.background = element_rect(fill = 'white', color = 'white'),
    legend.position = 'none',
    panel.border = element_blank(),
    panel.grid.minor = element_blank()) 

ps <- p1/p2

ggsave(ps, filename = 'RNAhybrid_2.png', path = wd, width = 5, height = 5, 
  device = png, dpi = 300)

# DATAVIZ W/
# SE PUEDE CONVERTIR A RANGOS CON IRANGE::RANGES <----
# ALGO COMO LO REALIZADO EN 1/2_SRNA_LOCATION_DB.R


.g <- howmirs %>%
  mutate(UTR = ifelse(target %in% three_prime_utr, 
    "three_prime_utr", "five_prime_utr")) %>%
  separate(col = target, into = c("seqnames", "ranges"), sep = "_") %>%
  separate(col = ranges, into = c("ranges", "strand"), sep = ":") %>%
  separate(col = ranges, into = c("start", "end"), sep = "-") %>%
  mutate_at(c("start", "end"), as.integer)
# group_by(MIRNA, UTR, query) %>%
# count(target) %>%
# ungroup()


.g %>% group_by(start, end)


target_gr <- GRanges(Rle(.g$seqnames), 
  ranges =  IRanges(start = .g$start, end = .g$end), 
  strand = .g$strand, 
  source = "RNAHybrid",
  UTR = .g$UTR,
  MIRNA = .g$MIRNA,
  mfe = .g$mfe, 
  pval = .g$pval)

# Remember:

length(target_gr) # 90 918

length(IRanges::reduce(target_gr)) # 33 663

length(range(target_gr)) # 404

IRanges::reduce(target_gr)



IRanges::resize(target_gr, width(target_gr), fix = "center") %>% 
  as_tibble() %>% group_by(start, end)

# IF WHICH TRACK ALL type OF SOURCE:
# length(x <- .gr[!.gr$type == "region"])
# target_df <- subsetByOverlaps(x = x, ranges = range(target_gr), type="any")


# TEST UTR FROM BIOMART (OMIT) ====== 
library(biomaRt)
library(tidyverse)

## list the available Ensembl marts and use Ensembl Genes

listMarts()
# biomart                version
# 1 ENSEMBL_MART_ENSEMBL      Ensembl Genes 109
# 2   ENSEMBL_MART_MOUSE      Mouse strains 109
# 3     ENSEMBL_MART_SNP  Ensembl Variation 109
# 4 ENSEMBL_MART_FUNCGEN Ensembl Regulation 109

listMarts(host="uswest.ensembl.org")

mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL")

datasets <- listDatasets(mart)
datasets %>% view()
## list the available datasets in this Mart
searchDatasets(mart = ensembl, pattern = "rufescens")

datasets <- listDatasets(mart = ensembl)

datasets %>% 
  filter(str_detect(description, 'gca02305543'))
