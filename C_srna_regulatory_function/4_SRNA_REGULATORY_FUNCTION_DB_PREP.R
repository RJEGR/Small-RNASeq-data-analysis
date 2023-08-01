# JOIN TRANSCRIPTOMIC DATA RELATED TO MIR-TARGETS:
#
# 1) READ SRNA_FUNCTION_PREDICTED.rds AND PULL GENE_ID
# 2) RETRIVE BACK TRANSCRIPT IDS (XM_*) OR gene_id (MSTRG*) USING COLUMN ref_gene_id 
# FROM THE TRANSCRIPTOME ASSEMBLY (transcripts.gtf) 
# 3) FILTER IDS FROM EXPRESSION MATRIX (gene_count_matrix.csv)
# 4) CORRELATE MIRS EXPRESSION W/ TARGET GENE FOUND IN TRANSCRIPTOME

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/FUNCTIONAL_MIR_ANNOT/"

# AFTER RUN 3_DESEQ2TOPGO.R
# LOAD DESEQ_RES_P WHICH INCLUDE GENE_ID TARGETED BY BOTH, RNAHYBRID AND TARGETSCAN:

print(RES.P <- read_tsv(paste0(wd, "DESEQ_RES_P.tsv")))

# nrow(out <- read_rds(paste0(wd, "SRNA_FUNCTION_PREDICTED.rds")))

str(query.ids <- RES.P %>% distinct(gene_id) %>% pull()) # 8465

# str(query.ids <- paste(query.ids, collapse = "|"))

  # mutate(query = strsplit(query, ";")) %>%
  # unnest(query) %>% 
  # distinct(query, GO.ID) %>%
  # group_by(query) %>%
  # summarise(
  #   across(GO.ID, .fns = paste_go), 
  #   .groups = "drop_last")

GTF <- "transcripts.gtf"

wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/RNA_SEQ_ANALYSIS/ASSEMBLY/STRINGTIE/REFBASED_MODE/"

f <- list.files(path = wd, pattern = "METADATA_RNASEQ", full.names = T)

.colData <- read_csv(f)

f <- list.files(path = wd, pattern = GTF, full.names = T)

GTF <- rtracklayer::import(f)

GTF2DF <- GTF %>% as_tibble() %>% distinct(gene_id, transcript_id, ref_gene_id)

GTF2DF <- GTF2DF %>% filter(!is.na(ref_gene_id))

nrow(GTF2DF <- GTF2DF %>% filter(ref_gene_id %in% query.ids)) # 17532


GTF2DF %>% count(ref_gene_id, sort = T)

# TRANSCRIPT ARE UNIQUE IDENTIFIER:
# gene_id ID AND ref_gene_id CLUSTERS TRANSCRIPT IDENTIFERS (SAME INFO)
# gene_id == ref_gene_id, EX. MSTRG.44204 == LOC124140464 W/ 41 DIFF. transcript_id (ISOFORMS)
# 

str(query2.ids <- GTF2DF %>% distinct(ref_gene_id) %>% pull() %>% sort()) # 8465

# SANITY CHECK

any(sort(query.ids) ==  query2.ids) # TRUE

# 3) ====
# WHICH ID RETRIVE BACK BETTER DETAILS?

# grep -c "XM_" transcript_count_matrix.csv 55 596 from 193 407
# rep -c "XM_" gene_count_matrix.csv Zero

# grep -c "LOC" transcript_count_matrix.csv 360 from 193 407
# grep -c "LOC" gene_count_matrix.csv 6077 from 58 593

# grep -c "MSTRG" transcript_count_matrix.csv 123 403
# grep -c "MSTRG" gene_count_matrix.csv 46 472

wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/RNA_SEQ_ANALYSIS/ASSEMBLY/STRINGTIE/REFBASED_MODE/QUANTIFICATION/"

f <- list.files(path = wd, pattern = "gene_count_matrix.csv", full.names = T)

.COUNT <- read_csv(f)

COUNT <- GTF2DF %>% 
  distinct(gene_id, ref_gene_id) %>% 
  left_join(.COUNT) %>%
  rename("assembled_id" = "gene_id", "gene_id" = "ref_gene_id")

View(COUNT)

DB <- RES.P %>% 
  mutate(gene_id = strsplit(gene_id, ";")) %>%
  unnest(gene_id) %>% 
  distinct(Name, gene_id) %>%
  group_by(gene_id) %>%
  summarise(
    across(Name, .fns = list), n_srnas = n(), .groups = "drop_last") %>%
  arrange(desc(n_srnas))

print(DB <- DB %>% left_join(COUNT))


# DB %>% unnest(Name) %>% view()

# View(DB)

head(COUNT <- DB %>% dplyr::select(starts_with("SRR")) %>% as(., "matrix"))

COUNT[is.na(COUNT)] <- 0

rownames(COUNT) <- DB$gene_id

COUNT <- log2(COUNT+1)

heatmap(COUNT)

# 
sample_dist = dist(t(COUNT), method='euclidean')

hc_samples = hclust(sample_dist, method='complete')

hc_order <- hc_samples$labels[hc_samples$order]

gene_dist = dist(COUNT, method='euclidean')

hc_gene = hclust(gene_dist, method='complete')

gene_order <- hc_gene$labels[hc_gene$order]

which_cols <- DB %>% dplyr::select(starts_with("SRR")) %>% names()

recode_to <- .colData %>% arrange(dpf) %>% distinct(Isolated) %>% pull()

pattern <- "Day Post-Fertilization Sample|Days Post-Fertilization Sample|Days Post-Fertilization"

recode_to <- structure(gsub(pattern, "DPF", recode_to ), names = recode_to)


DB_LONG <- DB %>% 
  pivot_longer(all_of(which_cols),names_to = "LIBRARY_ID") %>%
  # filter(value > 0) %>%
  mutate(value = ifelse(value == 0, NA, value)) %>%
  left_join(.colData) %>%
  dplyr::mutate(Isolated = dplyr::recode_factor(Isolated, !!!recode_to))
  # group_by(Isolated, gene_id) %>%
  # summarise(value = sum(value))
  # mutate(Isolated = factor(Isolated, levels = facet_level))

library(ggh4x)

DB_LONG %>%
  ggplot(aes(x = LIBRARY_ID, y = gene_id, fill = log2(value+1))) + 
  geom_tile(color = 'white', linewidth = 0.2) +
  scale_fill_viridis_c(option = "B", name = "Log2", direction = -1, na.value = "white") +
  ggh4x::facet_nested( ~ Isolated, nest_line = F, scales = "free", space = "free") +
  # scale_x_discrete(position = 'top') +
  # ggh4x::scale_x_dendrogram(hclust = hc_samples, position = 'top') +
  ggh4x::scale_y_dendrogram(hclust = hc_gene, position = "left", labels = NULL) +
  guides(y.sec = ggh4x::guide_axis_manual(labels = gene_order, breaks = gene_order, label_size = 3.5)) +
  theme_bw(base_size = 10, base_family = "GillSans") +
  labs(x = '', y = '') +
  theme(
    legend.position = "top",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    # axis.text.x = element_text(angle = 90, hjust = -0.15, vjust = 1),
    strip.background = element_rect(fill = 'grey89', color = 'white'),
    panel.border = element_blank(),
    plot.title = element_text(hjust = 0),
    plot.caption = element_text(hjust = 0),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank()) -> pright

pleft <- DB %>% dplyr::select(gene_id, n_srnas) %>% 
  # arrange(desc(match(gene_id, gene_order))) %>%
  mutate(gene_id = factor(gene_id, levels = gene_order)) %>%
  ggplot(aes(y = gene_id, x = n_srnas)) +
  geom_col() +
  theme_bw(base_size = 10, base_family = "GillSans") +
  labs(x = 'Number of miRs', y = '') +
  theme(
    legend.position = "top",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.background = element_rect(fill = 'grey89', color = 'white'),
    panel.border = element_blank(),
    plot.title = element_text(hjust = 0),
    plot.caption = element_text(hjust = 0),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank())

library(patchwork)

p <- pright + plot_spacer() + pleft + plot_layout(widths = c(7, -1.5 , 3.5)) 

ggsave(p, filename = "DESEQ2HEATMAP_TRANSCRIPTOME.png", path = wd, width = 5, height = 10, device = png, dpi = 300)


# COMO COLAPSAR INFORMACION DE NAME PARA IDENTIFICAR EL TARGET EN UN HEATMAP U OTRO DATAVIZ???

DB %>% unnest(Name) %>% 
  # dplyr::select(all_of(c("Name",which_cols))) %>% 
  right_join(RES.P %>% dplyr::select(Name, Family) %>% filter(Family == "MIR-278")) %>%
  view()
