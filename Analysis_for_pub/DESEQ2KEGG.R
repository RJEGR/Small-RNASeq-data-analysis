
# LOAD gene TRINOTATE 
# BIND INFO TO SRNA_REGULATORY_FUNCTION_DB.tsv
# LOAD DESEQ DATA
# USE TARGETED TRANSCRIPT_IDS AND DIFF-EXP MIRS INFO TO GET KEGG THEMES

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

# LOAD EGGS

# download.file("https://raw.githubusercontent.com/RJEGR/infovis/master/NOG.annotations.tsv", "NOG.annotations.tsv")

f <- "/Users/cigom/Documents/GitHub/Small-RNASeq-data-analysis/NOG.annotations.tsv"

egg <- read.table(f, sep="\t", stringsAsFactors=FALSE, quote="")

names(egg) <- c("db", "nog", "proteins", "species", "class", "description")

library(tidyverse)

ref_path <- "/Users/cigom/Documents/MIRNA_HALIOTIS/ENSEMBLE/ANNOTATIONS/"

cogs <- read_rds(paste0(ref_path, '/cogs.rds'))

annot_f <- list.files(path = ref_path, pattern = "Trinotate.xls", full.names = T)

annot <- read_tsv(annot_f, na = ".")

names(annot)[1] <- "gene_id"

annot <- annot %>% drop_na(gene_ontology_BLASTX)

wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/ENSEMBLE/"

genome_feat <- read_rds(paste0(wd, "genome_features.rds"))[[1]]

annot <- dplyr::select(annot, transcript_id, prot_id, eggnog, Kegg) %>% 
  left_join(genome_feat, by = "transcript_id")

# OR

blastx_df <- split_blast(annot, hit = "BLASTX") # HARD TO RUN THIS TIME

DB <- blastx_df %>% dplyr::filter(domain == "Eukaryota")

blastx_df %>% dplyr::count(genus, sort = T)


DB <- DB %>% 
    group_by(gene) %>%
    arrange(desc(identity)) %>%
    sample_n(1) %>%
    ungroup()

DB %>% dplyr::count(genus, sort = T)

genome_feat <- read_rds("/Users/cigom/Documents/MIRNA_HALIOTIS/ENSEMBLE/genome_features.rds")[[1]]

genome_feat <- genome_feat %>% dplyr::select(gene_id, transcript_id)

names(DB)[c(1,2)] <- c("gene_id", "transcript_id")

DB %>% dplyr::select(-align, -identity, -evalue, -gene_id) %>%
  left_join(genome_feat, by = "transcript_id")

# EXAMPLE:

ids <- annot %>% drop_na(eggnog) %>% sample_n(100) %>% pull(transcript_id)

kegg_df <- get_eggnog(x = annot, ids = ids, by = "transcript_id")


path_out <- "~/Documents/MIRNA_HALIOTIS/MIRNA_PUB_2024/"

ids <- read_tsv(paste0(path_out, "REFBASED_MODE_COUNT.tsv")) %>% pull(gene_id)

kegg_df <- get_eggnog(x = annot, ids = ids, by = "gene_id")

annot %>% filter(gene_id %in% ids) %>%
  split_kegg("Kegg") %>%
  dplyr::select(gene, Kegg) %>%
  # filter(grepl('KEGG', Kegg)) %>%
  separate(col = Kegg, sep = ":", into = c("db", "sp", "kegg")) %>%
  distinct(gene, kegg, .keep_all = T) -> Kegg


col_palette <- cogs %>% distinct(code, clrs)

col_palette <- structure(col_palette$clrs, names = col_palette$code)


cogs %>% left_join(kegg_df) %>% 
  arrange(desc(Freq)) %>%
  mutate(name = factor(name, levels = unique(name))) %>%
  ggplot(aes(x = name, y = Freq, fill = code)) +
  geom_col() +
  coord_flip() +
  geom_text(aes(label = Freq), size = 3, hjust = -0.05, family = "GillSans") +
  scale_fill_manual(values = col_palette)

# SPLIT BY MIRS TARGETS:
# BECAUSE NON EGGNOG FOR MIRS TARGET, LETS USE KEGG:

wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/FUNCTIONAL_MIR_ANNOT/"

# print(.SRNA2GO <- read_tsv(paste0(wd, "SRNA_REGULATORY_FUNCTION_DB.tsv")))

print(RES.P <- read_tsv(paste0(wd, "DESEQ_RES_P.tsv")))


DB <- RES.P %>%
  mutate(gene_id = strsplit(gene_id, ";")) %>%
  unnest(gene_id) %>%
  distinct(Family, gene_id) %>%
  # pull(gene_id)
  left_join(annot, by = "gene_id") %>% 
  distinct()

DB %>%
  split_kegg("Kegg") %>%
  dplyr::select(gene, Kegg) %>%
  # filter(grepl('KEGG', Kegg)) %>%
  separate(col = Kegg, sep = ":", into = c("db", "sp", "kegg")) %>%
  distinct(gene, kegg, .keep_all = T) -> Kegg

Kegg %>% dplyr::count(db)

# SEARCH FOR KEGG Entry/	NCBI-GeneID
# https://www.genome.jp/kegg-bin/get_htext?ko00001.keg
# try to download https://www.genome.jp/kegg/
# Ex.dre:100003958
# https://www.genome.jp/entry/dre:100003958
