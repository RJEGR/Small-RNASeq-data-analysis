

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

path <- "~/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out/"

list.files(path = path, pattern = "txt")

count_f <- list.files(path = path, pattern = "Counts.txt", full.names = T)

res_f <- list.files(path = path, pattern = "Results.txt", full.names = T)

MTD_f <- list.files(path = path, pattern = "METADATA.tsv", full.names = T)

MTD <- read_tsv(MTD_f)

RESULTS <- read_tsv(res_f)

COUNTS <- read_tsv(count_f)

colNames <- COUNTS %>% select_if(is.double) %>% names()

names(COUNTS)[names(COUNTS) %in% colNames] <- gsub(".clean.newid.subset", "", colNames)

colNames <- COUNTS %>% select_if(is.double) %>% names()

# #d73027

pHpalette <- c(`Low`="#ad1f1f", `Control`= "#4575b4")


# Sanity check

identical(RESULTS$Locus, COUNTS$Coords)
identical(RESULTS$Name, COUNTS$Name)

# cbind()

# UNIQ AND MULTIMAP READS ====

# Reads: Number of aligned sRNA-seq reads that overlap this locus.
# UniqueReads: Number of uniquely aligned (e.g. not multi-mapping) reads that overlap this locus.

RESULTS %>%
  mutate(MIRNA = ifelse(!is.na(KnownRNAs), "Y", MIRNA)) %>%
  count(MIRNA)

RESULTS <- RESULTS %>% 
  mutate(MIRNA = ifelse(!is.na(KnownRNAs), "Y", MIRNA)) %>%
  mutate(MIRNA = factor(MIRNA, levels = c("Y", "N")))

RESULTS %>% 
  ggplot(aes(Reads, UniqueReads)) + 
  facet_grid(~ Strand) +
  scale_y_log10() + scale_x_log10() +
  geom_point(data = subset(RESULTS, MIRNA == "N"), color = "red",alpha = 0.1) + 
  geom_point(data = subset(RESULTS, MIRNA == "Y"), color = "blue")
  #stat_smooth(method = "lm")

RESULTS %>% distinct(MajorRNA)

RESULTS %>% 
  mutate(KnownRNAs = ifelse(is.na(KnownRNAs) & MIRNA == "Y", Name, KnownRNAs)) %>%
  mutate(MIRNA = ifelse(!is.na(KnownRNAs), "Y", MIRNA)) %>%
  select(Name, Chrom, KnownRNAs, MIRNA) %>%
  drop_na(KnownRNAs) %>%
  mutate(KnownRNAs = ifelse(grepl("Cluster_", KnownRNAs), "Novel", "known")) %>%
  count(KnownRNAs)


RESULTS %>% 
  filter(MIRNA == "Y") %>% distinct(Chrom) 

RESULTS %>% 
  # filter(MIRNA == "Y") %>%
  # filter(Chrom == "JALGQA010000001.1") %>%
  pivot_longer(cols = c("Start", "End"), names_to = "names", values_to = "Coords") %>%
  ggplot(aes(x = Coords, y = log10(Reads))) + 
  theme_minimal() + facet_grid( MIRNA ~ .) +
  geom_line(colour="grey50", lineend = "round") +
  geom_area(fill="grey", alpha=0.6) +
  geom_hline(yintercept = 2, colour = "grey60", size = 0.5)
  


# PREVALENCE ----

nrow(raw_count <- COUNTS[colNames]) # 66,226

head(raw_count <- as.data.frame(raw_count, row.names = COUNTS$Name))

rownames(raw_count) <- COUNTS$Name

apply(raw_count, 1, function(x) sum(x > 0)) %>% table()

prevelancedf = apply(raw_count, 1, function(x) sum(x > 0))

# mean_se = apply(raw_count, 1, function(x) mean_se(x)) %>% do.call(rbind, .)

data.frame(Prevalence = prevelancedf, 
  TotalAbundance = rowSums(raw_count)) %>% # mean_se
  as_tibble(rownames = "Name") %>%
  arrange(desc(TotalAbundance)) -> prevelancedf

prevelancedf <- prevelancedf %>% left_join(RESULTS %>% select(Name, MIRNA, KnownRNAs, Strand, DicerCall))

prevelancedf %>% 
  count(Prevalence) %>% 
  ggplot(aes(Prevalence, n)) + geom_col() +
  theme_classic(base_family = "GillSans") + 
  scale_y_continuous("Number of sRNAs", labels = scales::comma) +
  scale_x_continuous(breaks = 1:12) -> ps

ggsave(ps, filename = 'prevalence_hist.png', 
  path = path, width = 3, height = 2, device = png)

prevelancedf %>% 
  arrange(Prevalence) %>%
  mutate(Prevalence = paste0(Prevalence, ' Samples')) %>%
  mutate(Prevalence = factor(Prevalence, levels = unique(Prevalence))) %>%
  mutate(TotalAbundance = log2(TotalAbundance+1)) %>% # edgeR::cpm(count) 
  ggplot(aes(TotalAbundance)) + geom_histogram() + 
  facet_wrap(~ Prevalence, scales = 'free_y') -> p1

dat_text <- prevelancedf %>% group_by(Prevalence) %>% tally() %>% 
  mutate(cumsum = cumsum(n)) %>% arrange(Prevalence) %>%
  mutate(Prevalence = paste0(Prevalence, ' Samples')) %>%
  mutate(Prevalence = factor(Prevalence, levels = unique(Prevalence)))

p1 <- p1 + geom_text(
  data    = dat_text, family = "GillSans",
  mapping = aes(x = -Inf, y = -Inf, label = paste0(n, " sRNAs")),
  hjust   = -1, vjust   = -2.5, size = 2.5) + 
  theme_classic(base_size = 7, base_family = "GillSans") +
  labs(x = expression(~Log[2]~('TotalAbundance'~+1)), y = "") +
  scale_y_continuous("Number of sRNAs", labels = scales::comma)


ggsave(p1, filename = 'prevalence_hist_facet.png', 
  path = path, width = 7, height = 4, device = png)


prevelancedf %>% 
  mutate(MIRNA = factor(MIRNA, levels = c("Y","N") )) %>%
  # filter(MIRNA == "Y") %>%
  ggplot(aes(TotalAbundance, Prevalence/12)) +
  # geom_point(aes(color = MIRNA), alpha = 0.2) +
  stat_density2d(aes(alpha = ..density..), geom = "tile", contour = FALSE) +
  scale_x_log10("Total Abundance (log10 scale)", labels = scales::comma) +  
  scale_y_continuous("Prevalence [Frac. Samples]", 
    labels = scales::percent_format(scale = 100)) +
  theme_classic(base_family = "GillSans") +
  theme(legend.position="top") -> ps
  # facet_grid(~ MIRNA) 
  # scale_color_manual(values = c("red", "grey"))


ggsave(ps, filename = 'prevalence.png', 
  path = path, width = 4, height = 4, device = png)


prevelancedf %>% 
  mutate(Prevalence = factor(Prevalence)) %>%
  # filter(MIRNA == "Y") %>%
  ggplot(aes(log10(TotalAbundance), color  = Prevalence, fill = Prevalence)) +
  # scale_x_continuous() +
  stat_ecdf() 
  # geom_density(alpha = 0.3)


# BOXPLOT ====
# 
# prevelancedf = apply(raw_count, 1, function(x) sum(x > 0))
# 
# mean_se = apply(raw_count, 1, function(x) mean_se(x)) %>% do.call(rbind, .)

qprobs <- function(x) { 
  x <- x[x > 1]
  quantile(x, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
}

apply(log2(raw_count+1), 2, qprobs) %>% 
  t() %>%
  as_tibble(rownames = 'LIBRARY_ID') %>%
  left_join(MTD) -> probs_df

probs_df %>%
  ggplot(., 
    aes(x = LIBRARY_ID, ymin = `5%`, lower = `25%`,
      middle = `50%`, upper = `75%`, ymax = `95%`)) +
  geom_errorbar(width = 0.3, position = position_dodge(0.6)) +
  geom_boxplot(aes(fill = pH), width = 0.5, stat = 'identity', 
    position = position_dodge(0.6)) +
  labs(y = expression(log[2]~ 'Abundance'), x = '') +
  theme_classic(base_family = "GillSans") +
  scale_color_manual("", values = rev(pHpalette)) +
  scale_fill_manual("", values = rev(pHpalette)) +
  theme(
    legend.position = 'top',
    panel.border = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank()) -> ptop

ptop <- ptop + facet_grid(~ hpf, scales = 'free') +
  theme(
    strip.background = element_rect(fill = 'grey', color = 'white'))

# BOTTOM ====

apply(raw_count, 2, function(x) sum(x > 0)) -> Total_genes

# Filter data by removing low-abundance genes

keep <- rowSums(edgeR::cpm(raw_count) > 1) >= 2


nrow(count <- raw_count[keep,])

# How singletones are per sample? ----

apply(count, 2, function(x) sum(x > 0)) -> filtered_genes

cbind(as_tibble(Total_genes, rownames = 'name'), as_tibble(filtered_genes)) -> n_genes

names(n_genes) <- c('LIBRARY_ID','Raw', 'Filt')

n_genes %>% mutate(pct_genes_retained = Filt/Raw) -> n_genes

n_genes %>%
  left_join(MTD) %>%
  mutate(LIBRARY_ID = factor(LIBRARY_ID, levels = unique(LIBRARY_ID))) %>%
  ggplot() + 
  geom_col(aes(x = LIBRARY_ID, y = Raw, fill = pH)) + 
  # geom_errorbar(aes(x = LIBRARY_ID, y = Filt, group = pH, ymin = Filt, ymax = Filt), color = 'black')
  labs(x = 'Samples', y = 'Total sRNAs') +
  scale_y_continuous(labels = scales::comma) +
  scale_color_manual("", values = rev(pHpalette)) +
  scale_fill_manual("", values = rev(pHpalette)) +
  theme_classic(base_family = "GillSans") +
  theme(legend.position = 'none',
    axis.text.x = element_text(angle = 45,
      hjust = 1, vjust = 1, size = 10)) -> pbottom
  
pbottom <- pbottom + facet_grid(~ hpf, scales = 'free') +
  theme(
    strip.background = element_blank(), 
    strip.text = element_blank())

library(patchwork)

ps <- ptop / pbottom + patchwork::plot_layout(heights = c(1,1.2))


ggsave(ps, filename = "transcripts_and_reads_plots.png", 
  path = path, width = 5, height = 4.5, device = png)


# PCA =====

# ncol(data <- log2(count+1))

ncol(data <- edgeR::cpm(raw_count+1))

PCA = prcomp(t(data), center = FALSE, scale. = FALSE)

percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)

sd_ratio <- sqrt(percentVar[2] / percentVar[1])

PCAdf <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2])

# library(mclust)

# d_clust <- mclust::Mclust(as.matrix(PCAdf), G=1:4, modelNames = mclust.options("emModelNames"))
# plot(d_clust)
# d_clust$BIC
# k <- d_clust$G
# names(k) <- d_clust$modelName

PCAdf %>% 
  dist(method = "euclidean") %>% 
  hclust() %>% 
  cutree(., k) %>% 
  as_tibble(rownames = 'LIBRARY_ID') %>% 
  mutate(cluster = paste0('C', value)) %>% 
  dplyr::select(-value) -> hclust_res


PCAdf %>%
  mutate(LIBRARY_ID = rownames(.)) %>%
  # mutate(g = substr(sample_id, 1,1)) %>%
  # left_join(hclust_res) %>%
  left_join(MTD) %>%
  ggplot(., aes(PC1, PC2)) +
  # coord_fixed(ratio = sd_ratio)
  geom_abline(slope = 0, intercept = 0, linetype="dashed", alpha=0.5) +
  geom_vline(xintercept = 0, linetype="dashed", alpha=0.5) +
  # ggforce::geom_mark_ellipse(aes(group = as.factor(cluster)),
  #   fill = 'grey', con.colour = 'grey') +
  geom_point(size = 7, alpha = 0.7, aes(color = pH)) +
  geom_text( family = "GillSans",
    mapping = aes(label = paste0(hpf, " hpf")), size = 2.5) +
  labs(caption = '') +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  scale_color_manual("", values = rev(pHpalette)) +
  scale_fill_manual("", values = rev(pHpalette)) +
  theme_classic(base_family = "GillSans") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = 'top') -> pcaplot

pcaplot

ggsave(pcaplot, 
  filename = "PCA.png", path = path, 
  width = 10, height = 7)

# Correlation heatmap ----

sample_cor = cor(data, method='pearson', use='pairwise.complete.obs')
sample_dist = dist(t(data), method='euclidean')
hc_samples = hclust(sample_dist, method='complete')


hc_order <- hc_samples$labels[hc_samples$order]

heatmap(sample_cor, col = cm.colors(12))

sample_cor %>% 
  as_tibble(rownames = 'LIBRARY_ID') %>%
  pivot_longer(cols = colnames(sample_cor), values_to = 'cor') %>%
  left_join(MTD) %>% 
  mutate(LIBRARY_ID = factor(LIBRARY_ID, levels = hc_order)) -> sample_cor_long

# sample_cor_long %>% distinct(sample_id, Diagnosis) %>% mutate(col = ifelse(g %in% 'C', 'red', 'blue')) -> coldf 

n <- length(unique(sample_cor_long$pH))

library(ggh4x)

sample_cor_long %>%
  ggplot(aes(x = LIBRARY_ID, y = name, fill = cor)) + 
  geom_tile(color = 'white', size = 0.2) +
  # geom_raster() + 
  # geom_text(aes(label = cor), color = 'white') +
  scale_fill_viridis_c(option = "B", name = "Pearson", direction = 1) +
  scale_x_discrete(position = 'top') +
  # ggh4x::scale_x_dendrogram(hclust = hc_samples, position = 'top') +
  # ggh4x::scale_y_dendrogram(hclust = hc_samples) +
  theme_classic(base_size = 12, base_family = "GillSans") +
  labs(x = '', y = '') +
  theme(axis.text.x = element_text(angle = 45,
    hjust = -0.15, vjust = 1)) -> pheat

# pheat +  ggh4x::scale_x_dendrogram(hclust = hc_samples, position = 'top') -> pheat

ggsave(pheat, 
  filename = "cor_data_matrix.png", path = path, 
  width = 6, height = 5, device = png)



# TEST DESEQ2 ----

# mtd <- read_tsv(list.files(path = path, pattern = 'METADATA', full.names = T))

# nrow(raw_count <- COUNTS[colNames])

# x) Choose contrast  ----

colData <- MTD

f_col <- "CONTRAST_A"

ref_group <- "Control"

names(colData)[names(colData) %in% f_col] <- "Design"

# 0) Remove single samples ----

sam <- table(colData$Design)

sam <- names(sam[sam > 1])

colData <- colData %>% drop_na(Design) %>% filter(Design %in% sam)

x <- colData %>% pull(Design, name = LIBRARY_ID)

# 1) According to metadata, filter samples ----

keep_cols <- names(raw_count) %in% colData$LIBRARY_ID

sum(keep_cols)

nrow(count <- raw_count[,keep_cols])

# dim(count <- raw_count[names(raw_count) %in% names(x)])

count <- as(count, "matrix")

rownames(count) <- COUNTS$Coords

# 2) Filter data by removing low-abundance genes ----

keep_rows <- rowSums(edgeR::cpm(count) > 1) >= 2

sum(keep_rows) # N transcripts

nrow(count <- count[keep_rows,])

dim(count <- round(count))


apply(count, 1, function(x) sum(x > 0)) %>% table()


# 3) Format metadata and count matrix ----

colData <- colData %>% arrange(match(LIBRARY_ID, names(count)))

sum(colnames(count) == colData$LIBRARY_ID) # sanity check


colData <- mutate_if(colData, is.character, as.factor)

# using relevel, just specifying the reference level:

colData <- colData %>% mutate(Design = relevel(Design, ref = ref_group))

library(DESeq2)
library(tidyverse)


table(colData$Design)

ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = count,
  colData = colData,
  design = ~ Design )

# 4) Run DESeq in the following functions order ----

dds <- estimateSizeFactors(ddsFullCountTable) 

# dds <- estimateDispersions(dds)
# dds <- nbinomWaldTest(dds) # use Wald test when contrasting a unique pair (ex. control vs cancer)

# Note: also run as a single batch 
# dds <- estimateSizeFactors(ddsFullCountTable)

dds <- DESeq(dds)

# 5) Test transformation for downstream analyses ----

# In order to test for differential expression, 
# we operate on raw counts and use discrete distributions ... 
# However for other downstream analyses – e.g. for visualization or clustering – it might be useful to work with transformed versions of the count data. Ex:

vst <- vst(dds) # wrapped version for varianceStabilizingTransformation

ntr <- DESeq2::normTransform(dds)

# 
# DESeq2::plotPCA(ntr, intgroup = "Design")
#
# DESeq2::plotPCA(vst, intgroup = "Design")

# BiocManager::install("vsn")

# VERIFICAR LA SD DE LOS DATOS TRANSFORMADOS, normTransform y vsd. Encontrar que transformacion disminuye la sd (linea roja) de los datos.

# ntr mantiene un sd debajo de 1 para la mayoria de los transcritos (eje x)

# vsn::meanSdPlot(assay(vst))
# vsn::meanSdPlot(assay(dds))
# vsn::meanSdPlot(assay(ntr)) 

# En este caso, ntr funciona mejor para estos datos


# 6 Generate results ----

get_res <- function(dds, contrast, alpha_cutoff = 0.1, lfc_cutoff = 0) {
  
  sA <- contrast[1]
  sB <- contrast[2]
  
  contrast <- as.character(DESeq2::design(dds))[2]
  
  keepA <- as.data.frame(colData(dds))[,contrast] == sA
  keepB <- as.data.frame(colData(dds))[,contrast] == sB
  
  contrast <- c(contrast, sA, sB)
  
  res = results(dds, contrast = contrast, alpha = alpha_cutoff, lfcThreshold = lfc_cutoff)
  
  
  baseMeanA <- rowMeans(DESeq2::counts(dds,normalized=TRUE)[,keepA])
  baseMeanB <- rowMeans(DESeq2::counts(dds,normalized=TRUE)[,keepB])
  
  res %>%
    as.data.frame(.) %>%
    cbind(baseMeanA, baseMeanB, .) %>%
    cbind(sampleA = sA, sampleB = sB, .) %>%
    as_tibble(rownames = "ids") %>%
    mutate(padj = ifelse(is.na(padj), 1, padj)) %>%
    mutate_at(vars(!matches("ids|sample|pvalue|padj")),
      round ,digits = 2)
}

contrast <- levels(colData(dds)$Design)

res <- get_res(dds, contrast) # # take a long

# write_rds(dds, file = paste0(path, '/CONTRAST_A_DDS.rds'))

# 7) multiple contrat from designs against control -----

DESeqDataSetFromInputs <- function(colData, count, f_col = ...) {
  
  
  if(!is.matrix(count)) stop()
  
  names(colData)[names(colData) %in% f_col] <- "Design"
  
  # 0) Remove single samples:
  
  sam <- table(colData$Design)
  
  sam <- names(sam[sam > 1])
  
  colData <- colData %>% drop_na(Design) %>% filter(Design %in% sam)
  
  x <- colData %>% pull(Design, name = LIBRARY_ID)
  
  # 1) According to colData, filter samples ----
  keep_cols <- colnames(count) %in% names(x)
  
  ncol(count <- count[,keep_cols])
  
  # 2) Filter data by removing low-abundance genes ----
  
  keep_rows <- rowSums(edgeR::cpm(count) > 1) >= 2
  
  nrow(count <- count[keep_rows,])
  
  count <- round(count)
  
  cat("Data to process: ", dim(count), "\n")
  
  # 3) Format metadata and count matrix ----
  
  colData <- colData %>% arrange(match(LIBRARY_ID, names(count)))
  
  sum(names(count) == colData$LIBRARY_ID) # sanity check
  
  colData <- mutate_if(colData, is.character, as.factor)
  
  # using relevel, just specifying the reference level:
  
  colData <- colData %>% mutate(Design = relevel(Design, ref = "Control"))
  
  library(DESeq2)
  library(tidyverse)
  
  ddsFullCountTable <- DESeqDataSetFromMatrix(
    countData = count,
    colData = colData,
    design = ~ Design )
  
  # 4) estimateSizeF for DESeq  ----
  
  dds <- estimateSizeFactors(ddsFullCountTable) 
  
  out <- DESeq(dds)
  
  write_rds(out, file = paste0(path, f_col, '_DDS.rds'))
  
  return(out)
}

nrow(raw_count <- COUNTS[colNames]) # 66,226

raw_count <- as(raw_count, "matrix")

rownames(raw_count) <- COUNTS$Name

# Warning, count data must be a matrix class object !!


dds_A <- DESeqDataSetFromInputs(mtd, raw_count, f_col = "CONTRAST_A") # 62,817 transcripts and 56 samples
dds_B <- DESeqDataSetFromInputs(mtd, raw_count, f_col = "CONTRAST_B") # 21,006 transcripts and 56 samples
dds_C <- DESeqDataSetFromInputs(mtd, raw_count, f_col = "CONTRAST_D") # 44,167 transcripts and 56 samples
dds_D <- DESeqDataSetFromInputs(mtd, raw_count, f_col = "CONTRAST_C") # 110533 transcripts and 56 samples

head(assay(dds_C))

# Generate results 
get_res <- function(dds, contrast, alpha_cutoff = 0.1, weighted_p = FALSE) {
  
  require(DESeq2)
  require(tidyverse)
  require(IHW)
  
  sA <- contrast[1]
  sB <- contrast[2]
  
  contrast <- as.character(DESeq2::design(dds))[2]
  
  keepA <- as.data.frame(colData(dds))[,contrast] == sA
  keepB <- as.data.frame(colData(dds))[,contrast] == sB
  
  contrast <- c(contrast, sA, sB)
  
  res <- results(dds, contrast, alpha = 0.1, lfcThreshold = 0)
  
  
  baseMeanA <- rowMeans(DESeq2::counts(dds,normalized=TRUE)[,keepA])
  baseMeanB <- rowMeans(DESeq2::counts(dds,normalized=TRUE)[,keepB])
  
  
  res <- res %>%
    as.data.frame(.) %>%
    cbind(baseMeanA, baseMeanB, .) %>%
    cbind(sampleA = sA, sampleB = sB, .) %>%
    mutate(padj = ifelse(is.na(padj), 1, padj)) %>%
    mutate_at(vars(!matches("transcript_id|sample|pvalue|padj")),
      round ,digits = 2)
  
  if(weighted_p) {
    
    #  take the set of p-values and then calculate the adjusted p-values using a independent Hypothesis Weighting
    
    ihwRes <- ihw(pvalue ~ baseMean,  data = res, alpha = 0.1) # toma menos de 1 minuto.
    
    weighted_pvalue <- as.data.frame(ihwRes)$weighted_pvalue
    
    
    res <- res %>% cbind(., weighted_pvalue = weighted_pvalue) %>% 
      as_tibble(rownames = "transcript_id")
    
  } else
  
  return(res)
}

dds2res <- function(dds_path_f, ...) {
  
  dds <- read_rds(dds_path_f)
  
  contrast <- levels(colData(dds)$Design)
  
  out <- list()
  
  for(i in 2:length(contrast)) {
    j <- i
    cat('\nContrast', contrast[1],' and ', contrast[j],'\n')
    res <- get_res(dds, contrast[c(1,j)], alpha_cutoff = 0.1)
    out[[j]] <- res
    
  }
  
  do.call(rbind, out) -> res
  
  return(res)
  
}

# data viz ----
rds_f <- list.files(path, pattern = 'CONTRAST_',  full.names = TRUE)

out <- lapply(rds_f, dds2res)

names(out) <- gsub("_DDS.rds", "", basename(rds_f))

res <- bind_rows(out, .id = 'Contrast') %>% as_tibble(rownames = "Locus")

write_rds(res, file = paste0(path, "MULTIPLE_CONTRAST_RES.rds"))

