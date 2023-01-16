
rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

path <- '~/Documents/MIRNA_HALIOTIS/CLEAN_INPUT/'


mtd <- read_tsv(list.files(path = path, pattern = 'METADATA', full.names = T))

fileName <- list.files(path = path, pattern = 'miRNAs_expressed_all_samples', full.names = T)

count <- read_tsv(file = fileName)

whichCols <- !grepl('norm', names(count))

count <- count[whichCols]

# count %>% distinct(`#miRNA`)

count <- count %>% filter(total > 0)

names(count)

summary <- count[,1:4]

summary %>% distinct(`#miRNA`)

count <- as.data.frame(count[,5:16])

row_id <- paste0(1:nrow(summary),"-" ,summary$`#miRNA`)

rownames(count) <- row_id


dim(count <- count[rowSums(edgeR::cpm(count) > 5 ) >= 2, ])

# Multiple contrast ----

colors_fc <- c("red2",  "#4169E1", "forestgreen", "grey30")
  
colData <- mtd %>% arrange(match(colName, names(count)))

colData <- mutate_if(colData, is.character, as.factor)
colData <- mutate_if(colData, is.double, as.factor)

# using relevel, just specifying the reference level:

colData$pH <- relevel(colData$pH, ref = "Control")

levels(colData$pH)

library(DESeq2)

count <- round(count)

table(colData$pH)

ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = count,
  colData = colData,
  design = ~ pH) # if not rep use design = ~ 1

dds <- DESeq(ddsFullCountTable)

# write_rds(dds, file = paste0(path, '/Cancer_vs_Control_dds.rds'))

# dds <- read_rds("~/Documents/DOCTORADO/human_cancer_dataset/DiffExp/Cancer_vs_Control_dds.rds")

contrast <- levels(colData(dds)$pH)

get_res <- function(dds, contrast, alpha_cutoff = 0.1) {
  
  sA <- contrast[1]
  sB <- contrast[2]
  
  contrast <- as.character(DESeq2::design(dds))[2]
  
  keepA <- as.data.frame(colData(dds))[,contrast] == sA
  keepB <- as.data.frame(colData(dds))[,contrast] == sB
  
  contrast <- c(contrast, sA, sB)
  
  res = results(dds, contrast, alpha = alpha_cutoff)
  
  
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

prep_DE_data <- function(res, alpha, lfcThreshold) {
  
  FC_cols <- c('logFC','log2FoldChange')
  pv_cols <- c('P.Value','pvalue', 'PValue')
  pvad_cols <- c('padj', 'FDR', 'adj.P.Val')
  
  sam <- c('sampleA',  'sampleB')
  samv <- c('baseMeanA', 'baseMeanB')
  
  rename_to <- c('logFC', 'pvalue', 'padj')
  # 
  # sigfc <- "<b>P</b>-value & Log<sub>2</sub> FC"
  # pv <- "<b>P</b>-value"
  # fc <- "Log<sub>2</sub> FC"
  # 
  # sigfc <- expression(p - value ~ and ~ log[2] ~ FC)
  # pv <- "p-value"
  # fc <- expression(Log[2] ~ FC)
  # 
  sigfc <- "p - value ~ and ~ log[2] ~ FC"
  pv <- "p-value"
  fc <- "Log[2] ~ FC"
  
  
  sample_NS <- function(res, n) {
    
    x <- filter(res, !(cc %in% c('NS', fc)))
    
    y <- sample_n(filter(res, cc == fc), n)
    
    z <- sample_n(filter(res, cc == 'NS'), n) 
    
    dat_sampled <- rbind(x,y,z)
    
    return(dat_sampled)
  }
  
  res %>%
    # drop_na() %>%
    select_at(vars(ids, contains(sam), contains(samv), contains(FC_cols), 
      contains(pv_cols),
      contains(pvad_cols))) %>%
    rename_at(vars(contains(FC_cols),
      contains(pv_cols),
      contains(pvad_cols)), ~ rename_to) %>%
    arrange(pvalue) -> res
  
  
  res$cc <- 'NS'
  res[which(abs(res$logFC) > lfcThreshold), 'cc'] <- fc
  res[which(abs(res$padj) <= alpha), 'cc'] <- pv
  res[which(res$padj <= alpha & abs(res$logFC) > lfcThreshold), 'cc'] <- sigfc
  
  up <- 'Up-regulated'
  down <- 'Down-regulated'
  
  res$lfcT <- 'Basal'
  res[which(res$padj < alpha & res$logFC > lfcThreshold), 'lfcT'] <- up
  res[which(res$padj < alpha & res$logFC < -lfcThreshold), 'lfcT'] <- down
  
  res %>%
    # sample_NS(.,1000) %>%
    mutate(cc = factor(cc, levels = c(sigfc, pv, fc, "NS"))) %>%
    mutate(lfcT = factor(lfcT, levels = c(up, down, 'Basal')))
  
}

res <- get_res(dds, contrast)

res.p <- prep_DE_data(res, alpha = 0.05, lfcThreshold = 2) %>% drop_na(cc)


logfc_in <- 2

padj_in <- 0.05 

# volcano ----

res.p %>%
  mutate(pvalue = -log10(pvalue)) %>%
  ggplot(aes(x = logFC, y = pvalue)) +
  geom_point(aes(color = cc), alpha = 3/5) +
  scale_color_manual(name = "", values = colors_fc) + 
  labs(x= expression(Log[2] ~ "Fold Change"), 
    y = expression(-Log[10] ~ "P")) +
  theme_bw(base_family = "GillSans", base_size = 12) +
  theme(legend.position = "top") -> pvol

pvol + geom_abline(slope = 0, intercept = -log10(padj_in), linetype="dashed", alpha=0.5) +
  geom_vline(xintercept = 0, linetype="dashed", alpha=0.5) -> pvol

ggsave(pvol,  filename = "volcano.png", path = path, width = 8, height = 5)


res.p %>%
  mutate(g = ifelse(logFC > 0, 'Control', 'Low pH')) %>%
  group_by(cc,g ) %>% tally() %>%
  filter(cc != 'NS') %>%
  ggplot() +
  geom_bar(aes(x = g, y = n, fill = cc), 
    stat = 'identity', position = position_dodge2()) +
  scale_fill_manual(name = "", values = colors_fc[-4]) +
  theme_bw(base_family = "GillSans", base_size = 12) +
  labs(y = '# microRNAs', x = '') +
  theme(legend.position = 'top') -> pbar


res.p %>%
  filter(cc != 'NS') %>%
  filter(cc ==  'p - value ~ and ~ log[2] ~ FC') %>%
  ggplot() +
  geom_histogram(aes(padj, fill = cc)) +
  # facet_wrap( cc ~ ., scales = 'free_x') +
  scale_fill_manual(name = "", values = colors_fc[-4]) +
  theme_bw(base_family = "GillSans", base_size = 12) +
  theme(legend.position = 'none') +
  labs(y = '# microRNAs') -> padjp

library(patchwork)

pbar / padjp -> savep

ggsave(savep, 
  filename = "diffExp_bar.png", path = path, 
  width = 5, height = 5)

# Heatmap ----

res %>% distinct(sampleA, sampleB)

res %>% arrange(log2FoldChange)

# res %>% sample_n(1000) %>%
#   ggplot() +
#   geom_point(aes(y = log2FoldChange, x = baseMean, color = -log10(padj)))


alpha = 0.05 

lfcThreshold = 2

# Queremos genes que fueron superiores al umbral lfcThreshold & < alpha == FC+sig

res %>% 
  mutate(cc = NA) %>% 
  mutate(cc = ifelse(padj < alpha & abs(log2FoldChange) > lfcThreshold, 'sigfc', cc)) %>%
  drop_na(cc) -> res.p

up <- 'Up-regulated'
down <- 'Down-regulated'

res.p %>%
  mutate(lfcT = 'Basal') %>%
  mutate(lfcT = ifelse(log2FoldChange > lfcThreshold, up, lfcT)) %>%
  mutate(lfcT = ifelse(log2FoldChange < -lfcThreshold, down, lfcT)) -> res.p


res.p %>%
  # filter(cc ==  'p - value ~ and ~ log[2] ~ FC') %>%
  pull(ids) -> diffexpList

keep <- rownames(count) %in% diffexpList

dim(count[keep,] -> DEcount)

sample_cor = cor(DEcount, method='pearson', use='pairwise.complete.obs')

# superheat::superheat()

dist.method <- 'euclidean'
linkage.method <- 'complete'

m <- log2(DEcount+1)



hc_samples <- hclust(dist(t(m), method = dist.method), 
  method = linkage.method)

# sample_dist = dist(t(DEcount), method='euclidean')
# sample_dist = dist(sample_cor, method='euclidean')

hc_sam_order <- hc_samples$labels[hc_samples$order]

hc_genes <- hclust(dist(m, method = dist.method), method = linkage.method)

hc_genes_ord <- hc_genes$labels[hc_genes$order]

DEcount %>% 
  as_tibble(rownames = 'id') %>%
  pivot_longer(cols = colnames(DEcount), values_to = 'count', names_to = 'colName') %>%
  # mutate(sample_id = factor(sample_id, levels = hc_sam_order)) %>%
  # mutate(id = factor(id, levels = hc_genes_ord)) %>%
  left_join(mtd) %>%
  filter(count > 0) -> countLong

countLong %>% 
  # distinct(colName, Diagnosis) %>% 
  # mutate(col = ifelse(g1 %in% 'C', 'red', 'blue')) %>%
  arrange(match(colName, hc_sam_order)) -> coldf 

coldf %>% mutate()

n <- length(unique(coldf$pH))

# getPalette <- RColorBrewer::brewer.pal(n, 'Paired')

getPalette <- c( "#4169E1", "red2")

axis_col <- structure(getPalette, names = unique(coldf$pH))
axis_col <- getPalette[match(mtd$pH, names(axis_col))]
axis_col <- structure(axis_col, names = unique(coldf$colName))



library(ggh4x)

countLong %>%
  # sample_n(1000) %>%
  ggplot(aes(x = colName, y = id, fill = log2(count+1))) + 
  geom_raster() + 
  theme_classic(base_size = 12, base_family = "GillSans") +
  scale_fill_viridis_c(name = expression(Log[2]~'(count+1)')) +
  labs(x = '', y = '') +
  ggh4x::scale_y_dendrogram(hclust = hc_genes) +
  ggh4x::scale_x_dendrogram(hclust = hc_samples, position = 'top') +
  theme(axis.text.x = element_text(angle = 90, 
    hjust = 1, vjust = 1, size = 7, color = axis_col), # 
    axis.ticks.length = unit(10, "pt"),
    # legend.position = 'top',
    panel.border = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank()) -> pheat


pheat + guides(fill = guide_colorbar(barheight = unit(4, "in"),
  ticks.colour = "black",
  frame.colour = "black",
  label.theme = element_text(size = 12))) -> pheat
