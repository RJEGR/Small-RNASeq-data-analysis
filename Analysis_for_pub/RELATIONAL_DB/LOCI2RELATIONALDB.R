# INTRAGENIC INSPECTION
# collectively. Such intragenic (Polycistrons as cooperative functional units.) co- regulated miRNAs expressed from the same clus- ter have a tendency to target the same gene or target different genes in the same pathways. this reinforces the network-regulating roles of these miRNAs, as exem- plified by the co-targeting of the actin cytoskeleton pathway by members of the miR-200 family. (doi:10.1038/nrg.2016.134).

# LOAD LOCI DATABASE
# BIND WITH RELATIONALDB

# COUNT NUMBER OF INTRAGENIC AGAINST INTERGENIC TARGETDE GREE


library(tidyverse)

dir <- "~/Documents/MIRNA_HALIOTIS/MIRNA_PUB_2024/"

print(DB <- read_tsv(file.path(dir, "LONGER_RELATIONAL_DB.tsv")))

DB <- DB %>% 
  mutate(Loci_geneid = strsplit(Loci_geneid, ";")) %>%
  unnest(Loci_geneid) %>%
  mutate(STRINGID = strsplit(STRINGID, ";")) %>%
  unnest(STRINGID)

col_recode <- structure(c("#FFC107", "#2196F3","red"), 
  names = c("Intergenic", "Intragenic", "Other ncRNA"))
  
DB %>% 
  mutate(GENESET = ifelse(is.na(STRINGID), gene_id, STRINGID)) %>%
  drop_na(GENESET) %>% 
  # distinct(MajorRNA, STRINGID, biotype_best_rank) %>% 
  count(GENESET, biotype_best_rank, sort = T) %>%
  arrange(desc(n)) %>%
  mutate(GENESET = factor(GENESET, levels = unique(GENESET))) %>%
  ggplot(aes(x = GENESET, y = n, color = biotype_best_rank)) + 
  geom_point(shape = 21) +
  geom_step(aes(group = biotype_best_rank)) +
  theme_classic(base_family = "GillSans", base_size = 12) +
  labs(y = "microRNA degree", x = "Gene target") +
  theme(legend.position = 'top', 
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 7)) +
  scale_color_manual("",values = col_recode) +
  scale_fill_manual("",values = col_recode) 

dir <- "/Users/cigom/Documents/MIRNA_HALIOTIS/ENSEMBLE/"

genedescr <- read_rds(paste0(dir, "genome_features.rds"))[[2]]

genedescr <- genedescr %>% 
  mutate(gene_id = strsplit(gene_id, ";")) %>%
  unnest(gene_id)

DB %>% distinct(MajorRNAID, Loci_geneid) %>% drop_na() %>% 
  left_join(genedescr, by = c("Loci_geneid" = "gene_id")) %>%
  select(MajorRNAID, Loci_geneid, description) %>%
  view()

