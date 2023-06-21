

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)


library(tidyverse)

library(GenomicRanges)

wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/FUNCTIONAL_MIR_ANNOT"

# GET gene features from locus ====

pattern <- "multi_genome.newid.gff3$"

genome <- list.files(path = wd, pattern = pattern, full.names = T)

.gr <- rtracklayer::import(genome)

# .gr %>% as_tibble() %>% dplyr::count(type)
# .gr %>% as_tibble() %>% drop_na(description) %>% dplyr::count(biotype)

# ...
# 4 protein_coding 31171 # <---
# ...
#

gr <- .gr[which(.gr$type == "gene")]

gr <- gr[which(gr$biotype == "protein_coding")]

# Sanity check
# 31 171 of protein_coding with gene/protein description:

gr %>% as_tibble() %>% drop_na(description) %>% dplyr::count(biotype) 


which_cols <- c("seqnames", "gene_coords", "gene_id","description", "type", "biotype")

gene_features <- gr %>% as_tibble() %>% 
  drop_na(description) %>%
  mutate(gene_coords = paste(start, end, strand, sep = ":")) %>%
  select_at(all_of((which_cols)))

# gene_features <- gene_features %>% distinct(gene_id, description)

# from gene feature to UTR (flat information)

utr_f <- list.files(path = wd, full.names = T, pattern = "three_prime_utr.ids")

str(x <- read_lines(utr_f)) # 54 432 <-- i.e the N sequences from three_prime_utr.fa

str(gene_id <- sapply(strsplit(x, " "), `[`, 2)) # LOC*

str(target <- sapply(strsplit(x, " "), `[`, 1)) # three_prime_utr ids 

nrow(utr_source <- data.frame(target, gene_id) %>% as_tibble())

nrow(utr_source <- utr_source %>% distinct(target, gene_id)) # 31,654

# Las columnas start:end del objeto gene_features corresponden 
# a las coordenadas de todo el gene/protein_coding de cada gene_id

nrow(utr_source <- utr_source %>% left_join(gene_features)) # 31, 654

# READ RNAHYBRYD CSV ====


# pattern = "mir_vs_utr_rmdup_RNAhybrid.out.tsv$"

pattern <- "mir_vs_utr_rmdup_RNAhybrid.out.psig.tsv"


f <- list.files(path = wd, pattern = pattern, full.names = T)

# df <- read_delim(f, col_names = F, delim = "\t")
# colNames <- c("target", "query", "mfe", "pval", "pos", "lenT", "lenQ")
# colnames(df) <- colNames

df <- read_tsv(f)

# str(targetsids <- df %>% distinct(target) %>% pull(target))

# write_lines(targetsids, file = paste0(gsub(".tsv", "", f), ".ids"))

# df <- df %>% filter(pval < 0.05)

# write_tsv(df, file = paste0(gsub(".tsv", "", f), ".psig.tsv"))


# Q: how many unique binding sites? =====

nrow(df %>% distinct(target)) #  37,663
nrow(df %>% distinct(query)) # 261 from 393 clusters from mir.fasta


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

# mir.fasta: 
# This is a FASTA formatted file containing hairpin, mature miRNA, and miRNA* sequences derived from ShortStack's identification of MIRNA loci.

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
  facet_grid(~ MIRNA) +
  geom_point(shape = 1, alpha = 0.5) +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color="grey") +
  geom_vline(xintercept = -35, linetype = "dashed", color="grey") +
  theme_bw(base_size = 16, base_family = "GillSans") +
  labs(x = "mfe (kcal/mol)") +
  theme(legend.position = "top")


# Q: How many mirs per target:? =====

mir_sites <- howmirs %>%
  mutate(UTR = ifelse(target %in% three_prime_utr, 
    "three_prime_utr", "five_prime_utr")) %>%
  group_by(target, UTR) %>%
  count(MIRNA)

sum(mir_sites$n)

mir_sites %>% 
  # group_by(UTR, MIRNA, n) %>%
  # tally(n) %>% view() 
  group_by(UTR, MIRNA) %>%
  rstatix::get_summary_stats(n, type = "common")

p2 <- mir_sites %>% 
  ggplot(aes(n, color = UTR, fill = UTR)) +
  facet_grid(~ MIRNA) +
  geom_histogram() +
  # stat_ecdf(linewidth = 2, alpha = 0.5) +
  theme_bw(base_size = 16, base_family = "GillSans") +
  labs(x = "binding miRs per target") +
  scale_y_continuous("Target gene", labels = scales::comma) +
  theme(legend.position = "none")

library(patchwork)

# p1/p2

# Q: how many chromosomes were binding sites for mirs? ====

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
  theme(legend.position = "none")

library(patchwork)

p1/p2/p3

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



# Find overlaping features (OMIT) ====
# 


ov <- findOverlaps(grr, gr, minoverlap = 2)

to_index <- unique(subjectHits(ov)) # position vector where coordinates from subject is found

from_index <- unique(queryHits(ov)) # position vector where coordinates from query is found

gr[sort(to_index)] %>%  as_tibble() %>% drop_na(description) %>% select(seqnames, start, end, description)

grr[sort(from_index)]

# Rsubread::flattenGTF(genome, GTF.featureType = "three_prime_UTR", GTF.attrType = "gene_id")

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
