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

# print(RES.P <- read_tsv(paste0(wd, "DESEQ_RES_P.tsv")))

# wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/FUNCTIONAL_MIR_ANNOT/"

print(out <- read_rds(paste0(wd, "SRNA_FUNCTION_PREDICTED.rds")))

out <- out %>% filter(predicted == "BOTH")

str(query.ids <- out %>% distinct(gene_id) %>% pull())
  
# str(query.ids <- RES.P %>% 
#     mutate(gene_id = strsplit(gene_id, ";")) %>%
#     unnest(gene_id) %>%
#     distinct(gene_id) %>% pull()) 

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

print(GTF2DF <- GTF %>% 
    as_tibble() %>% 
    distinct(gene_id, transcript_id, ref_gene_id) %>% 
    filter(!is.na(ref_gene_id))) # A tibble: 70,001 × 3

# nrow(GTF2DF <- GTF2DF %>% filter(ref_gene_id %in% query.ids)) # 17532

GTF2DF %>% dplyr::count(ref_gene_id, sort = T)

# TRANSCRIPT ARE UNIQUE IDENTIFIER:
# gene_id ID AND ref_gene_id CLUSTERS TRANSCRIPT IDENTIFERS (SAME INFO)
# gene_id == ref_gene_id, EX. MSTRG.44204 == LOC124140464 W/ 41 DIFF. transcript_id (ISOFORMS)
# 

str(query2.ids <- GTF2DF %>% distinct(ref_gene_id) %>% pull() %>% sort()) # 8465

# SANITY CHECK
length(query.ids) # 168 predicted targets

any(sort(query.ids) %in%  query2.ids) # TRUE

sum(sort(query.ids) %in%  query2.ids) # 155 expressed  of  168 predicted targets

write_rds(out %>% filter(gene_id %in% query2.ids), 
  file = paste0(wd, "SRNA_FUNCTION_PREDICTED_EXPRESSED.rds"))


# 3) USING TRANSCRITPOME DATA ====
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

.COUNT %>% distinct(gene_id) %>% nrow() # 58592 genes assembled

any(GTF2DF$gene_id %in% .COUNT$gene_id)

GTF2DF <- GTF2DF[GTF2DF$gene_id %in% .COUNT$gene_id,]

COUNT <- GTF2DF %>% 
  distinct(gene_id, ref_gene_id) %>%
  right_join(.COUNT) %>%
  rename("assembled_id" = "gene_id", "gene_id" = "ref_gene_id")

which_cols <- COUNT %>% dplyr::select(dplyr::starts_with("SRR")) %>% names()

# optional, change transcript(ie. isoform ) to gene mode

COUNT <- COUNT %>%
  group_by(gene_id) %>%
  summarise_at(vars(all_of(which_cols)), sum) 

# filter out Mantle samples?
keep <- .colData %>% filter(!Tissue %in% "Mantle") %>% pull(LIBRARY_ID)

keep <- c("gene_id",which_cols[which_cols %in% keep])

COUNT <- COUNT %>% dplyr::select(dplyr::starts_with(keep))

# View(COUNT)

DB <- out %>%
  mutate(query = strsplit(query, ";")) %>%
  unnest(query) %>%
  group_by(gene_id) %>%
  summarise(
    across(query, .fns = list), n_srnas = n(), .groups = "drop_last") %>%
  arrange(desc(n_srnas))

# DB <- RES.P %>%
#   mutate(gene_id = strsplit(gene_id, ";")) %>%
#   unnest(gene_id) %>%
#   distinct(Name, gene_id) %>%
#   group_by(gene_id) %>%
#   summarise(
#     across(Name, .fns = list), n_srnas = n(), .groups = "drop_last") %>%
#   arrange(desc(n_srnas))

print(DB <- DB %>% left_join(COUNT))

DB <- DB %>% drop_na()
# DB %>% unnest(Name) %>% view()

# View(DB)

head(COUNT <- DB %>% dplyr::select(starts_with("SRR")) %>% as(., "matrix"))

sum(keep <- rowSums(COUNT) > 1)

dim(COUNT <- COUNT[keep,])

COUNT[is.na(COUNT)] <- 0 # add pseudo-count

rownames(COUNT) <- DB$gene_id[keep]

colSums(COUNT)

COUNT <- DESeq2::varianceStabilizingTransformation(round(COUNT))
# 
# rowSums(COUNT)[1]
# 
# COUNT[1,]
# 
# z_scores(COUNT[1,])

z_scores <- function(x) {(x-mean(x))/sd(x)}

COUNT <- apply(COUNT, 1, z_scores)

COUNT <- t(COUNT)

# COUNT <- log2(COUNT+1)
COUNT[is.na(COUNT)] <- 0

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


DB_LONG <- COUNT %>% as_tibble(rownames = "gene_id") %>%
  pivot_longer(all_of(which_cols),names_to = "LIBRARY_ID") %>%
  # filter(value > 0) %>%
  # mutate(value = ifelse(value == 0, NA, value)) %>%
  left_join(.colData) %>%
  dplyr::mutate(Isolated = dplyr::recode_factor(Isolated, !!!recode_to))
  # group_by(Isolated, gene_id) %>%
  # summarise(value = sum(value))
  # mutate(Isolated = factor(Isolated, levels = facet_level))


lo = floor(min(DB_LONG$value))
up = ceiling(max(DB_LONG$value))
mid = (lo + up)/2


library(ggh4x)

DB_LONG %>%
  # mutate(value = ifelse(value == 0, NA, value)) %>%
  ggplot(aes(x = LIBRARY_ID, y = gene_id, fill = value)) + 
  geom_tile(color = 'white', linewidth = 0.2) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
    na.value = "white", midpoint = mid, limit = c(lo, up),
    name = NULL) +
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
    # axis.text.x = element_blank(),
    # axis.ticks.x = element_blank(),
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

p
ggsave(p, filename = "DESEQ2HEATMAP_TRANSCRIPTOME.png", path = wd, width = 5, height = 10, device = png, dpi = 300)


# COMO COLAPSAR INFORMACION DE NAME (I.E MIRS) PARA IDENTIFICAR EL TARGET EN UN HEATMAP U OTRO DATAVIZ???

# GAD, 
# CORRELACIONAR CONTEOS TOTAL DE MIRS POR TARGET
# GENERAR UNA TABLA QUE CONTENGA TRES COLUMNAS:
# gene_id, EXPRESSION, MIR-TARGET-EXPRESSION
# 1) LEER CONTEOS DE MICRORNAS Y SUMAR POR TARGET
# 2) LEER CONTEOS DE GENE-EXPRESSION RNA-SEQ Y SUMAR INDEPENDIENTE DE LA ETAPA DE DESARROLLO
# 3) VERIFICAR SI ALGUNO DE LOS TARGET SE ENCUENTRA EXPRESADO DURANTE EL DESARROLLO (VISUAL)
# 4) HACER JOIN DE LOS CONTEOS DE LOS MIR-TARGET-EXPRESSION Y GENE-EXPRESSION
# 5) GRAFICAR UN SCATTERPLOT CON COEFICIENTES PEARSON 
# (OPTIONAL) CALCULAR MATRIZ DE CORRELACION POR ETAPA DE DESARROLLO?

wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out/"

dim(MIR_COUNT <- read_rds(paste0(wd, "COUNT.rds"))) # must match 147 rows x 12 cols

wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/FUNCTIONAL_MIR_ANNOT/"

TARGET <- read_tsv(paste0(wd, "LOCI2TARGETDB.tsv"))

TARGET %>% distinct(Name) %>% nrow() # must match 147 mirs

bind_srnas <- function(x) {
  x <- paste(x, sep = ';', collapse = ';')
}

MIR_EXPRESSION <- TARGET %>% 
  group_by(Name) %>%
  summarise(across(target, .fns = bind_srnas), .groups = "drop_last") %>%
  right_join(rowSums(MIR_COUNT) %>% as_tibble(rownames = 'Name')) %>% view()
  mutate(target = strsplit(target, ";")) %>%
  unnest(target) %>%
  group_by(target) %>%
  summarise(across(Name, .fns = bind_srnas), n_mirs= n(), mir_exp = sum(value),
    .groups = "drop_last") %>%
  dplyr::rename("gene_id" = "target")

# 2) LEER CONTEOS DE GENE-EXPRESSION RNA-SEQ Y SUMAR INDEPENDIENTE DE LA ETAPA DE DESARROLLO

GTF <- "transcripts.gtf"

wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/RNA_SEQ_ANALYSIS/ASSEMBLY/STRINGTIE/REFBASED_MODE/"

f <- list.files(path = wd, pattern = "METADATA_RNASEQ", full.names = T)

.colData <- read_csv(f)

f <- list.files(path = wd, pattern = GTF, full.names = T)

GTF <- rtracklayer::import(f)

print(GTF2DF <- GTF %>% as_tibble() %>% 
    distinct(gene_id, transcript_id, ref_gene_id) %>% 
    filter(!is.na(ref_gene_id))) # A tibble: 70,001 × 3

wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/RNA_SEQ_ANALYSIS/ASSEMBLY/STRINGTIE/REFBASED_MODE/QUANTIFICATION/"

f <- list.files(path = wd, pattern = "gene_count_matrix.csv", full.names = T)

.COUNT <- read_csv(f)

.COUNT %>% distinct(gene_id) %>% nrow() # 58592 genes assembled

any(GTF2DF$gene_id %in% .COUNT$gene_id)

GTF2DF <- GTF2DF[GTF2DF$gene_id %in% .COUNT$gene_id,]

COUNT <- GTF2DF %>% 
  distinct(gene_id, ref_gene_id) %>%
  right_join(.COUNT) %>%
  rename("assembled_id" = "gene_id", "gene_id" = "ref_gene_id")

# 3) VERIFICAR SI ALGUNO DE LOS TARGET SE ENCUENTRA EXPRESADO DURANTE EL DESARROLLO (VISUAL) ====

# how many targets?

TARGET %>% distinct(target) %>% nrow() # 168

q.targets <- TARGET %>% distinct(target) %>% pull()

COUNT %>% filter(gene_id %in% q.targets) %>% nrow() # 155 targets are usually expressed during larval dev.

COUNT <- COUNT %>% filter(gene_id %in% q.targets)

which_cols <- COUNT %>% select_if(is.double) %>% names()

z_scores <- function(x) {(x-mean(x))/sd(x)}

M <- COUNT %>% dplyr::select(starts_with("SRR")) %>% as(., "matrix")

M <- apply(M, 1, z_scores)

M[is.na(M)] <- 0

heatmap(M)

dpfLevel <- c("0.5 dpf","1 dpf", "6 dpf","7 dpf","21 dpf", "30 dpf","Female" )

LONGERCOUNT <- COUNT %>% 
  filter(gene_id %in% q.targets) %>%
  pivot_longer(all_of(which_cols),names_to = "LIBRARY_ID") %>%
  filter(value > 0) %>%
  # mutate(value = ifelse(value == 0, NA, value)) %>%
  left_join(.colData) %>%
  group_by(gene_id, CONTRAST_B) %>%
  summarise(value = sum(value)) %>%
  mutate(CONTRAST_B = factor(CONTRAST_B, levels = dpfLevel))

LONGERCOUNT %>%
  mutate(value = z_scores(value)) %>%
  ggplot(aes(x = CONTRAST_B, y = gene_id, fill = value)) + 
  geom_tile(color = 'white', linewidth = 0.2) +
  scale_fill_viridis_c(option = "B", name = "z", direction = -1, na.value = "white")

# 4) HACER JOIN DE LOS CONTEOS DE LOS MIR-TARGET-EXPRESSION Y GENE-EXPRESSION ====
# 5) GRAFICAR UN SCATTERPLOT CON COEFICIENTES PEARSON 


GENE_EXPRESSION <- LONGERCOUNT %>% summarise(gene_exp = sum(value))

GENE_EXPRESSION %>% left_join(MIR_EXPRESSION) %>%
  arrange(desc(n_mirs)) %>%
  mutate(gene_exp = log10(gene_exp), mir_exp = log10(mir_exp)) %>%
  ggplot(aes(gene_exp, n_mirs)) + geom_point() +
  geom_smooth(method = "lm", linetype="dashed", size = 0.5, alpha=0.5,
    se = F, na.rm = TRUE) +
  ggpubr::stat_cor(method = "spearman", cor.coef.name = "R", p.accuracy = 0.001) +
  theme_bw(base_size = 16, base_family = "GillSans")

# (OPTIONAL) CALCULAR MATRIZ DE CORRELACION POR ETAPA DE DESARROLLO?
library(rstatix)

# M <- COUNT %>% dplyr::select(starts_with("SRR")) # %>% as(., "matrix")
# 
# M <- apply(M, 1, z_scores)
# 
# M[is.na(M)] <- 0
# 

cor.mat <- COUNT %>% 
  # mutate_at(vars(contains(c("SRR"))), z_scores) %>%
  left_join(MIR_EXPRESSION) %>%
  dplyr::select(contains(c("SRR", "n_mirs", "mir_exp"))) %>% 
  rstatix::cor_mat()

cor.mat %>%
  cor_reorder() %>%
  pull_lower_triangle() %>%
  cor_plot(label = TRUE)
