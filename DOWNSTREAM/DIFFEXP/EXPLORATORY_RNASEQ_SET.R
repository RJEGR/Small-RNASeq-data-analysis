
# HACE FALTA CORRER TRINOTATE: 
# REVISAR. PORQUE HAY SOLO 149,366 transcritos (transcript.fa) y la matriz de conteos tiene 193,406
# i.e REHACER EL ARCHIVO DE SECUENCIAS USANDO EL GTF ADECUADO
# /home/rvazquez/RNA_SEQ_ANALYSIS/ASSEMBLY/STRINGTIE/REFBASED_MODE/ANNOTATION

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)


# 

dir <- "~/Documents/MIRNA_HALIOTIS/"

wd <- paste0(dir, "RNA_SEQ_ANALYSIS/ASSEMBLY/STRINGTIE/REFBASED_MODE/")

count_f <- list.files(path = paste0(wd, "QUANTIFICATION"), pattern = "transcript_count_matrix", full.names = T)

# blastp_f <- list.files(path = paste0(wd, "ANNOTATION/OUTPUTS/"), pattern = ".rds", full.names = T)

# Cargamos archivos en R  ====

load(paste0(wd, "ANNOTATION/OUTPUTS/annot.Rdata"))

str(query.ids <- blastp_df %>% distinct(transcript) %>% pull(transcript)) # 1 348

dim(count <- read_csv(count_f)) # caution: na = c("", "NA")

count <- count %>% filter(transcript_id %in% query.ids) # sample_n(1000) # SUBSAMPLE

count[is.na(count)] <- 0

query.ids <- count$transcript_id

count <- count %>% select_if(is.double) %>% as(., "matrix")

rownames(count) <- query.ids

# 1) Filter data by removing low-abundance genes ----
# 

by_count <- 1; by_freq <- 2

keep <- rowSums(count > by_count) >= by_freq

sum(keep) # N transcripts

nrow(count <- count[keep,])

count <- round(count) + 1

# 2) Format metadata and count matrix ----

# Usamos objeto colData para elaborar el diseño experimental

colData <- read_csv(list.files(path = wd, pattern = "METADATA_RNASEQ.csv", full.names = T))

colData[is.na(colData)] <- "Low CO2"

colData <- colData %>% arrange(match(LIBRARY_ID, colnames(count)))

any(colnames(count) == colData$LIBRARY_ID) # sanity check

colData <- mutate_if(colData, is.character, as.factor)

# Usamos objeto colData para elaborar el diseño experimental

f_col <- "CONTRAST_A"

names(colData)[names(colData) %in% f_col] <- "Design"

# Si es posible, especificar las el factor "control" como nivel de referencia:

colData <- colData %>% mutate(Design = relevel(Design, ref = "Low CO2"))

# DESEQ2 =====

library(DESeq2)

ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = count,
  colData = colData,
  design = ~ Design )

dds <- estimateSizeFactors(ddsFullCountTable) 

# sizeFactors(dds)

dds <- estimateDispersions(dds)

# boxplot(log2(counts(ddsFullCountTable)+0.5))
# boxplot(log2(counts(dds, normalized=TRUE)+0.5))

dds <- nbinomWaldTest(dds) 

res <- results(dds)

summary(res)

# out of 1023 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 10, 0.98%
# LFC < 0 (down)     : 29, 2.8%
# outliers [1]       : 668, 65%
# low counts [2]     : 0, 0%
# (mean count < 4)

res <- res %>% as_tibble(rownames = "transcript") %>% arrange(padj)

write_rds(res, file = paste0(wd, "Low-CO2-High-CO2-DEGS-orfs.rds"))

# plotCounts(dds, gene= "MSTRG.1734.1", intgroup="Design")

# BY RULE:

# POSITIVE LFC == UP EXPRESSED in High CO2 (but down-regulated in Ref Group)
# NEGATIVE LFC == UP EXPRESSED in Low CO2 (Ref Group)

sigfc <- "Sign and FC";pv <- "Sign";fc <- "FC"

res$cc <- 'NS'
res[which(abs(res$log2FoldChange) > 2), 'cc'] <- fc
res[which(abs(res$padj) <= 0.05), 'cc'] <- pv
res[which(res$padj <= 0.05 & abs(res$log2FoldChange) > 2), 'cc'] <- sigfc

# Join annotation ====

res <- read_rds(paste0(wd, "Low-CO2-High-CO2-DEGS-orfs.rds"))

blastp_df <- blastp_df %>% group_by(transcript) %>% arrange(identity) %>% sample_n(1) %>% ungroup()

res <- blastp_df %>% right_join(res)

genes <- res$transcript

gene1 <- genes[which.min(res$log2FoldChange)]
gene2 <- genes[which.max(res$log2FoldChange)]

d1 <- plotCounts(dds, gene=gene1, intgroup="Design", returnData=TRUE) %>% mutate(transcript = gene1)
d2 <- plotCounts(dds, gene=gene2, intgroup="Design", returnData=TRUE) %>% mutate(transcript = gene2)

# 
go_df2 <- go_df %>% filter(ontology == "biological_process") %>% group_by(transcript) %>% sample_n(1)

df <- rbind(d1,d2) %>% left_join(go_df2)

ggplot(df, aes(x=Design, y=count, group = Design)) +
  facet_grid(~ name) +
  geom_jitter(position=position_jitter(w=0.1,h=0), size = 5, alpha = 0.5) +
  scale_y_log10() +
  stat_summary(fun.data = "mean_cl_boot", colour = "red", linewidth = 0.7, size = 1) +
  labs(y = "Gene count (Log10)", x = "Treatment") +
  theme_bw(base_family = "GillSans", base_size = 20)


# Subset/subseq sequence =====
library(Biostrings)

# f <- list.files(path = paste0(wd, "ANNOTATION/OUTPUTS/"), pattern = "transcripts.fa.transdecoder.renamed.cds$", full.names = T)

f <- list.files(path = wd, pattern = "transcripts.fa", full.names = T)

seqs <- Biostrings::readDNAStringSet(f, format="fasta")

# query.ids <- res %>% distinct(protein) %>% pull()

query.ids <- res %>% distinct(transcript) %>% pull()

sum(keep <- names(seqs) %in% query.ids)

seqs <- seqs[keep]

gff_f <- list.files(path = paste0(wd, "ANNOTATION/OUTPUTS/"), pattern = "prime_UTR.gff3$", full.names = T)

gff <- rtracklayer::import(gff_f) %>% as_tibble()

# No hay consistencia de las coordenadas y los ids de las secuencias (seqs)

coords <- gff %>% filter(seqnames %in% query.ids)

coords <- coords %>% mutate(Parent = unlist(Parent)) %>% select(seqnames, Parent, start, end, strand, type)

coords %>% dplyr::count(strand)

coords %>% arrange(desc(strand)) %>% distinct(Parent)

coords %>% filter(seqnames %in% "MSTRG.15170.2")

seqs[names(seqs) %in% "MSTRG.15170.2"]

# subseq

n <- nrow(coords)

pb <- txtProgressBar(min = 0, max = n, style = 3, width = 50, char = "=")

out <- list()

for(j in 1:n) {
  
  i <- j
  
  p1 <- coords[i,] %>% mutate(p = as.numeric(paste0(strand, start))) %>% pull(p)
  
  p2 <- coords[i,] %>% mutate(p = as.numeric(paste0(strand, end))) %>% pull(p)
  
  seqname <- coords[i,] %>% pull(seqnames)
  
  x <- seqs[names(seqs) %in% seqname]
  
  x <- as.character(unlist(x))
  
  x <- unlist(strsplit(x, ""))
  
  subsq <- paste(x[p1:p2], collapse = "")
  
  out[[i]] <- fasta <- c(rbind(paste0(">", seqname), subsq))
  
  setTxtProgressBar(pb, i)
  
}

unlist(out)

write(fasta, file= paste0(path, "MajorRNA.fasta"))






# namef <- paste(levels(colData$Design), collapse = "-", sep = "-")

# "Low-CO2-High-CO2"

# ShortRead::writeFasta(seqs, file = paste0(wd, "ANNOTATION/OUTPUTS/", "Low-CO2-High-CO2-DEGS.cds"))



# (OMIT) =====
# Generar un ciclo para hacer contrate multiple entre los estadios y el resto de las etapas

# for (i in unique(colData$dpf)) {}

DESeqDataSetFromInputs <- function(colData, count, f_col = ...) {
  
  
  if(!is.matrix(count)) stop()

  colData <- colData %>% drop_na(Design)
  
  sum(keep <- colnames(count) %in% colData$LIBRARY_ID)
  
  head(count <- count[,keep])
  
  # 1) Filter data by removing low-abundance genes ----
  
  by_count <- 1; by_freq <- 2
  
  keep <- rowSums(count > by_count) >= by_freq
  
  sum(keep) # N transcripts
  
  nrow(count <- count[keep,])
  
  head(count <- round(count))
  
  # 2) Format metadata and count matrix ----
  
  colData <- colData %>% arrange(match(LIBRARY_ID, colnames(count)))
  
  any(colnames(count) == colData$LIBRARY_ID) # sanity check
  
  # Si es posible, especificar las el factor "control" como nivel de referencia:
  
  colData <- mutate_if(colData, is.character, as.factor)
  
  reflev <- "Control" # "Low CO2"
  
  colData <- colData %>% mutate(Design = relevel(Design, ref = reflev))
  
  library(DESeq2)
  library(tidyverse)
  
  ddsFullCountTable <- DESeqDataSetFromMatrix(
    countData = count,
    colData = colData,
    design = ~ Design )
  
  # 4) estimateSizeF for DESeq  ----
  
  # dds <- estimateSizeFactors(ddsFullCountTable) 
  
  # out <- DESeq(dds)
  
  dds <- estimateSizeFactors(ddsFullCountTable) 
  
  dds <- estimateDispersions(dds)
  
  dds <- nbinomWaldTest(dds) 
  
  # res <- results(dds)
  
  
  write_rds(dds, file = paste0(path, f_col, '_DDS.rds'))
  
  return(dds)
}

for (i in unique(colData$dpf)) {
  
  CONTRAST <- colData$dpf %in% unique(colData$dpf)[i]
  
  colData$CONTRAST <- "Control"
  
  colData$CONTRAST[CONTRAST] <- "Contrast"
  
  f_col <- "CONTRAST"
  
  names(colData)[names(colData) %in% f_col] <- "Design"
  
  DESeqDataSetFromInputs(colData, count, f_col = "CONTRAST_A")
}

CONTRAST <- colData$dpf %in% unique(colData$dpf)[1]

colData$Design <- "Control"

colData$Design[CONTRAST] <- "Contrast"

.colData <- colData

colData <- colData %>% drop_na(Design)

sum(keep <- colnames(count) %in% colData$LIBRARY_ID)

head(count <- count[,keep])

# 1) Filter data by removing low-abundance genes ----
 
by_count <- 1; by_freq <- 2

keep <- rowSums(count > by_count) >= by_freq

sum(keep) # N transcripts

nrow(count <- count[keep,])

head(count <- round(count))

# 2) Format metadata and count matrix ----

colData <- colData %>% arrange(match(LIBRARY_ID, colnames(count)))

any(colnames(count) == colData$LIBRARY_ID) # sanity check

# Si es posible, especificar las el factor "control" como nivel de referencia:

colData <- mutate_if(colData, is.character, as.factor)

reflev <- "Control" # "Low CO2"

colData <- colData %>% mutate(Design = relevel(Design, ref = reflev))

# DESEQ2 ====

library(DESeq2)

ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = count,
  colData = colData,
  design = ~ Design )

# dds <- estimateSizeFactors(ddsFullCountTable)

# out <- DESeq(dds) # to low also !!!

dds <- estimateSizeFactors(ddsFullCountTable) 

dds <- estimateDispersions(dds)

dds <- nbinomWaldTest(dds) # <-- to low!!! 

res <- results(dds)

summary(res)
contrast <- levels(colData(dds)$Design)


get_res(dds, contrast)

# DATAVIZ

res <- res %>% as_tibble(rownames = "transcript")


sigfc <- "Sign and FC";pv <- "Sign";fc <- "FC"

colors_fc <- c("red2",  "#4169E1", "forestgreen", "grey30")

colors_fc <- structure(colors_fc, names = c(sigfc, pv, fc, "NS"))

res$cc <- 'NS'
res[which(abs(res$log2FoldChange) > 2), 'cc'] <- fc
res[which(abs(res$padj) <= 0.05), 'cc'] <- pv
res[which(res$padj <= 0.05 & abs(res$log2FoldChange) > 2), 'cc'] <- sigfc

p <- res  %>%
  ggplot(aes(y = -log10(pvalue), x = log2FoldChange, color = cc)) +
  geom_point()  

p <- p +
  scale_color_manual(name = "", values = colors_fc) +
  labs(x= expression(Log[2] ~ "Fold Change"), y = expression(-Log[10] ~ "P")) +
  theme_bw(base_family = "GillSans", base_size = 18) +
  theme(legend.position = "top")  +
  geom_abline(slope = 0, intercept = -log10(0.05), linetype="dashed", alpha=0.5) +
  geom_vline(xintercept = 0, linetype="dashed", alpha=0.5) 


# POSITIVE LFC == UP EXPRESSED IN CANCER
# NEGATIVE LFC == UP EXPRESSED IN CONTROL

res %>% mutate(g = sign(log2FoldChange)) %>% 
  dplyr::count(cc, g) %>%
  mutate(g = ifelse(g == "1", "High CO2", "Low CO2")) %>%
  filter(cc != 'NS') %>%
  group_by(g) %>%
  mutate(pct = n / sum(n)) %>%
  ggplot() +
  geom_col(aes(x = g, y = pct, fill = cc), width = 0.5) +
  scale_fill_manual(name = "", values = colors_fc[-4]) +
  theme_bw(base_family = "GillSans", base_size = 18) +
  labs(y = '% Transcripts', x = '') +
  theme(legend.position = 'top') + coord_flip()

# PREVALENCE ----

apply(.count, 1, function(x) sum(x > 0)) %>% table()

prevelancedf = apply(.count, 1, function(x) sum(x > 0))

data.frame(Prevalence = prevelancedf, 
  TotalAbundance = rowSums(.count)) %>% # mean_se
  as_tibble(rownames = "Name") %>%
  arrange(desc(TotalAbundance)) -> prevelancedf

prevelancedf %>% 
  dplyr::count(Prevalence) %>% 
  ggplot(aes(Prevalence, n)) + geom_col() +
  theme_minimal(base_family = "GillSans", base_size = 18) + 
  scale_y_continuous("Number of transcript", labels = scales::comma) 

# BOXPLOT ====

qprobs <- function(x) { 
  x <- x[x > 1]
  quantile(x, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
}

apply(log2(.count+1), 2, qprobs) %>% 
  t() %>%
  as_tibble(rownames = 'LIBRARY_ID') %>%
  left_join(.colData) -> probs_df

probs_df %>%
  ggplot(., 
    aes(x = LIBRARY_ID, ymin = `5%`, lower = `25%`,
      middle = `50%`, upper = `75%`, ymax = `95%`)) +
  geom_errorbar(width = 0.3, position = position_dodge(0.6)) +
  geom_boxplot(aes(fill = Design), width = 0.5, stat = 'identity', 
    position = position_dodge(0.6)) +
  labs(y = expression(log[2]~ 'Abundance'), x = '') +
  theme_classic(base_family = "GillSans", base_size = 18) +
  theme(
    legend.position = 'top',
    panel.border = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank()) -> ptop

ptop <- ptop + facet_grid(~ Design, scales = 'free') +
  theme(
    strip.background = element_rect(fill = 'grey', color = 'white'))

# BOTTOM ====

apply(.count, 2, function(x) sum(x > 0)) -> Total_genes

# How singletones are per sample? ----

apply(.count, 2, function(x) sum(x > 0)) -> filtered_genes

cbind(as_tibble(Total_genes, rownames = 'name'), as_tibble(filtered_genes)) -> n_genes

names(n_genes) <- c('LIBRARY_ID','Raw', 'Filt')

n_genes %>% mutate(pct_genes_retained = Filt/Raw) -> n_genes

n_genes %>%
  left_join(.colData) %>%
  mutate(LIBRARY_ID = factor(LIBRARY_ID, levels = unique(LIBRARY_ID))) %>%
  ggplot() + 
  geom_col(aes(x = LIBRARY_ID, y = Raw, fill = Design)) + 
  labs(x = 'Samples', y = 'Total Transcript') +
  scale_y_continuous(labels = scales::comma) +
  theme_minimal(base_family = "GillSans", base_size = 18) +
  theme(legend.position = 'none',
    axis.text.x = element_text(angle = 45,
      hjust = 1, vjust = 1, size = 10)) -> pbottom

pbottom + facet_grid(~ dpf, scales = 'free', space = "free") 

pbottom <- pbottom + facet_grid(~ Design, scales = 'free') +
  theme(
    strip.background = element_blank(), 
    strip.text = element_blank())

library(patchwork)

ptop / pbottom + patchwork::plot_layout(heights = c(1,1.2))

# 6) PCA =====

# ncol(data <- log2(count+1))

vst <- DESeq2::vst(ddsFullCountTable) # vst if cols > 10

ncol(data <- assay(vst))

PCA = prcomp(t(data), center = FALSE, scale. = FALSE)

percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)

sd_ratio <- sqrt(percentVar[2] / percentVar[1])

PCAdf <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2])

PCAdf %>%
  mutate(LIBRARY_ID = rownames(.)) %>%
  left_join(colData) %>%
  arrange(dpf) %>% mutate(Isolated = factor(Isolated, levels = unique(Isolated))) %>%
  ggplot(., aes(PC1, PC2)) +
  # coord_fixed(ratio = sd_ratio) +
  geom_abline(slope = 0, intercept = 0, linetype="dashed", alpha=0.5) +
  geom_vline(xintercept = -74, linetype="dashed", alpha=0.5) +
  # ggforce::geom_mark_ellipse(aes(group = as.factor(dpf)),
    # fill = 'grey', con.colour = 'grey') +
  geom_point(size = 12, alpha = 0.7, aes(color = Isolated)) +
  geom_text( family = "GillSans",
    mapping = aes(label = paste0(dpf, " dpf")), size = 3.5) +
  labs(caption = '', color = "") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme_minimal(base_family = "GillSans", base_size = 18) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = 'top')

# 5) (OPTIONAL) Correlation heatmap ----

sample_cor = cor(assay(vst), method='pearson', use='pairwise.complete.obs')
sample_dist = dist(t(assay(vst)), method='euclidean')
hc_samples = hclust(sample_dist, method='complete')


hc_order <- hc_samples$labels[hc_samples$order]

# heatmap(sample_cor, col = cm.colors(12))

sample_cor %>% 
  as_tibble(rownames = 'LIBRARY_ID') %>%
  pivot_longer(cols = colnames(sample_cor), values_to = 'cor') %>%
  left_join(colData) %>% 
  mutate(LIBRARY_ID = factor(LIBRARY_ID, levels = hc_order)) -> sample_cor_long

library(ggh4x)

sample_cor_long %>%
  ggplot(aes(x = LIBRARY_ID, y = name, fill = cor)) + 
  geom_tile(color = 'white', linewidth = 0.2) +
  # geom_raster() + 
  # geom_text(aes(label = cor), color = 'white') +
  scale_fill_viridis_c(option = "B", name = "Pearson", direction = -1) +
  scale_x_discrete(position = 'top') +
  ggh4x::scale_x_dendrogram(hclust = hc_samples, position = 'top') +
  ggh4x::scale_y_dendrogram(hclust = hc_samples) +
  theme_classic(base_size = 12, base_family = "GillSans") +
  labs(x = '', y = '') +
  theme(axis.text.x = element_text(angle = 50,
    hjust = -0.15, vjust = 1)) # -> pheat

