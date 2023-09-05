
# LOAD gene TRINOTATE 
# BIND INFO TO SRNA_REGULATORY_FUNCTION_DB.tsv
# LOAD DESEQ DATA
# USE TARGETED TRANSCRIPT_IDS AND DIFF-EXP MIRS INFO TO GET KEGG THEMES

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

# LOAD EGGS

download.file("https://raw.githubusercontent.com/RJEGR/infovis/master/NOG.annotations.tsv", "NOG.annotations.tsv")

egg <- read.table("NOG.annotations.tsv", sep="\t", stringsAsFactors=FALSE, quote="")

names(egg) <- c("db", "nog", "proteins", "species", "class", "description")

library(tidyverse)

ref_path <- "/Users/cigom/Documents/MIRNA_HALIOTIS/ENSEMBLE/ANNOTATIONS/"

cogs <- read_rds(paste0(ref_path, '/cogs.rds'))

annot_f <- list.files(path = ref_path, pattern = "Trinotate.xls", full.names = T)

annot <- read_tsv(annot_f, na = ".")

# EXAMPLE:

ids <- annot %>% drop_na(eggnog) %>% sample_n(100) %>% pull(transcript_id)

kegg_df <- get_eggnog(x = annot, ids = ids)


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

