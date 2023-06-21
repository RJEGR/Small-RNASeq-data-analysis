# TARGETSCAN

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/FUNCTIONAL_MIR_ANNOT"

f <- "mature_star_mir_vs_mir_vs_utr_rmdup_RNAhybrid.out.psig_targetscan.out"

f <- list.files(path = wd, pattern = f, full.names = T)

df <- read_tsv(f)

df %>% distinct(a_Gene_ID) # 861

df %>% distinct(miRNA_family_ID) # 262

df %>%
  mutate(MIRNA = "hairpin") %>%
  mutate(MIRNA = ifelse(grepl("mature", miRNA_family_ID), "mature", MIRNA)) %>%
  mutate(MIRNA = ifelse(grepl("star", miRNA_family_ID), "star", MIRNA)) %>%
  count(MIRNA)

# OVERLAP WITH RNAHYBRID
