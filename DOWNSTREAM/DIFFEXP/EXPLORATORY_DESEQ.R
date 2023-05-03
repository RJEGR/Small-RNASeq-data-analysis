
rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

library(DESeq2)

source("~/Documents/GitHub/EDAT/DE_methods.R")

path <- '~/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out/'

mtd <- read_tsv(list.files(path = path, pattern = 'METADATA', full.names = T))

# fileName <- list.files(path = path, pattern = 'miRNAs_expressed_all_samples', full.names = T)
fileName <- list.files(path = path, pattern = "Counts.txt", full.names = T)

Countdf <- read_tsv(file = fileName)

# whichCols <- !grepl('norm', names(count))

# count <- count[whichCols]

Countdf %>% dplyr::count(MIRNA)

# count %>% distinct(`#miRNA`)

count <- Countdf %>% dplyr::filter(total > 0)

names(count)

summary <- count[,1:4]

count %>% distinct(`#miRNA`)

summary %>% dplyr::count(`#miRNA`, sort = T)

summary %>% filter(`#miRNA`== 'dre-miR-430a-3p')

count <- as.data.frame(count[,5:16])

row_id <- paste0(1:nrow(summary),"-" ,summary$`#miRNA`)

rownames(count) <- row_id

dim(count)

dim(count0 <- count[rowSums(edgeR::cpm(count) > 1 ) >= 2, ])

# DESEQ2 -----
# Select samples and run deseq2 

which_sam <- function(mtd, f_col = "hpf", f = 24, ref = "Control") {
  
  samples_id <- mtd$colName
  
  # mtd %>% filter_at(vars(starts_with(f_col)), any_vars((. %% f) == 0))
  
  subset_mtd <- filter(mtd, if_any(starts_with(f_col), ~ (.x == f) == 0))
  
  x <- subset_mtd %>% pull(pH, name = colName)
  
  y <- samples_id %in% names(x)
  
  y[y == FALSE] <- "X"
  
  y <- structure(y, names = samples_id)
  
  y[names(y) %in% names(x)] <- as.double(relevel(factor(x), ref = ref))
  
  y <- paste(y, collapse = '')
  
  return(y)
  
}


# Ex: g <- "000XX1XXXXX1XXXXXXX2XXX2XXXXXXXXXXXX"
  
g1 <- which_sam(mtd, f_col = "hpf", f = 24)
g2 <- which_sam(mtd, f_col = "hpf", f = 110)

# which_sam(mtd, f_col = "pH", f = "low")
# which_sam(mtd, f_col = "pH", f = "Control", ref = 24)  # not work!
# count, g

res1 <- run_DESEQ2(count0, g1)
res2 <- run_DESEQ2(count0, g2)

# DESEQ DATAVIZ ----
alpha <- 0.05

res1.p <- prep_DE_data(res1, alpha = 0.05, lfcThreshold = 2) %>% mutate(hpf = "24 hpf")
res2.p <- prep_DE_data(res2, alpha = 0.05, lfcThreshold = 2)  %>% mutate(hpf = "110 hpf")

res.p <- rbind(res1.p , res2.p)

# Volcano ----

res.p %>%
  mutate(hpf = factor(hpf, levels = c("24 hpf", "110 hpf"))) %>%
  mutate(pvalue = -log10(pvalue)) %>%
  ggplot(aes(x = logFC, y = pvalue)) +
  geom_point(aes(color = cc), alpha = 3/5) +
  scale_color_manual(name = "", values = colors_fc) + 
  labs(x= expression(Log[2] ~ "Fold Change"), 
    y = expression(-Log[10] ~ "P")) +
  theme_bw(base_family = "GillSans", base_size = 12) +
  theme(legend.position = "top")  +
  geom_abline(slope = 0, intercept = -log10(alpha), linetype="dashed", alpha=0.5) +
  geom_vline(xintercept = 0, linetype="dashed", alpha=0.5) -> pvol

pvol <- pvol + facet_wrap(~ hpf, scales = 'free_y')

ggsave(pvol,  filename = "volcano.png", path = path, width = 6, height = 3)

res.p %>% arrange(logFC)

# count0[rownames(count0) %in% '844-ami-miR-193b-3p',]

# SI ES NEGATIVO ES EN LOW PH

res.p %>%
  mutate(g = ifelse(logFC > 0, 'Control', 'Low pH')) %>%
  group_by(cc,g, hpf) %>%
  tally() %>%
  mutate(hpf = factor(hpf, levels = c("24 hpf", "110 hpf"))) %>%
  filter(cc != 'NS') %>%
  ggplot() +
  geom_bar(aes(x = g, y = n, fill = cc), 
    stat = 'identity', position = position_dodge2()) +
  scale_fill_manual(name = "", values = colors_fc[-4]) +
  theme_bw(base_family = "GillSans", base_size = 12) +
  labs(y = '# microRNAs', x = '') +
  theme(legend.position = 'top') -> pbar

pbar <- pbar + facet_wrap(~ hpf, scales = 'free_y')
  
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

#how up and down genes ----

title <- 'Differentially expressed features'

subtitle <- expression(p[adj] < 0.05 ~ and ~ abs ~ (log[2] ~ FC) > 2)

table(res.p$cc)

res.p %>%
  # filter(cc ==  'p - value ~ and ~ log[2] ~ FC') %>%
  group_by(hpf, sampleA, lfcT) %>%
  tally() %>% mutate(pct = n / nrow(count0)) %>%
  mutate(hpf = factor(hpf, levels = c("24 hpf", "110 hpf"))) %>%
  ggplot(aes(x = hpf, y = pct, fill = lfcT)) +
  geom_col(color = 'black') +
  scale_fill_manual(name = "", values = c('#1f78b4', '#a6cee3', "grey80")) +
  theme_bw(base_family = "GillSans", base_size = 12) +
  labs(y = '# microRNAs', x = '', subtitle = subtitle, title = title) +
  theme(legend.position = 'top') -> psave

psave <- psave + coord_flip()

ggsave(psave, filename = "up_down_genes.png", path = path,
  width = 4, height = 2.5)


# Heatmap ----

res2.p %>%
  filter(cc ==  'Log[2] ~ FC') %>%
  pull(ids) -> diffexpList2

res1.p %>%
  filter(cc == 'p - value ~ and ~ log[2] ~ FC') %>%
  pull(ids) -> diffexpList1

str(diffexpList <- c(diffexpList1, diffexpList2))
str(diffexpList <- unique(diffexpList))

# 
# res.p %>%
#   filter(cc ==  'p - value ~ and ~ log[2] ~ FC') %>%
#   group_by(hpf) %>%
#   arrange(padj) %>%
#   # slice_head(n = 20) %>%
#   pull(ids) -> diffexpList

keep <- rownames(count0) %in% diffexpList

dim(count0[keep,] -> DEcount)

sample_cor = cor(DEcount, method='pearson', use='pairwise.complete.obs')

# superheat::superheat()

dist.method <- 'euclidean'

linkage.method <- 'complete'

# m <- log2(DEcount+1)

m <- edgeR::cpm(DEcount)



hc_samples <- hclust(dist(t(m), method = dist.method), 
  method = linkage.method)

# sample_dist = dist(t(DEcount), method='euclidean')
# sample_dist = dist(sample_cor, method='euclidean')

hc_sam_order <- hc_samples$labels[hc_samples$order]

hc_genes <- hclust(dist(m, method = dist.method), method = linkage.method)

hc_genes_ord <- hc_genes$labels[hc_genes$order]

DEcount %>% 
  # edgeR::cpm(DEcount)
  as_tibble(rownames = 'id') %>%
  pivot_longer(cols = colnames(DEcount), values_to = 'count', names_to = 'colName') %>%
  mutate(colName = factor(colName, levels = hc_sam_order)) %>%
  mutate(id = factor(id, levels = hc_genes_ord)) %>%
  left_join(mtd) %>%
  filter(count > 0) -> countLong

?rstatix::get_summary_stats()

countLong %>%
  group_by(id, hpf, pH) %>% 
  # rstatix::get_summary_stats(type = "mean_se") # TO LATE
  summarise(count = mean(count), n = n()) -> countLong

countLong %>% 
  # distinct(colName, Diagnosis) %>% 
  # mutate(col = ifelse(g1 %in% 'C', 'red', 'blue')) %>%
  arrange(match(colName, hc_sam_order)) -> coldf 

n <- length(unique(coldf$pH))

# getPalette <- RColorBrewer::brewer.pal(n, 'Paired')

getPalette <- c( "#4169E1", "red2")

axis_col <- structure(getPalette, names = unique(coldf$pH))
axis_col <- getPalette[match(mtd$pH, names(axis_col))]
axis_col <- structure(axis_col, names = unique(coldf$colName))



library(ggh4x)

countLong %>%
  ggplot(aes(x = pH, y = id, fill = log2(count+1))) + 
  facet_grid(~ hpf) +
  geom_raster() + 
  theme_classic(base_size = 12, base_family = "GillSans") +
  scale_fill_viridis_c(name = expression(Log[2]~'(count+1)')) +
  labs(x = '', y = '') +
  # ggh4x::facet_nested(~ hpf+pH,nest_line = T, scales = 'free_x') +
  ggh4x::scale_y_dendrogram(hclust = hc_genes) +
  # ggh4x::scale_x_dendrogram(hclust = hc_samples, position = 'top') +
  theme(axis.text.x = element_text(angle = 45, 
    hjust = 1, vjust = 1, size = 14), # color = axis_col
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

# pheat + facet_grid(~ hpf+pH, scales = 'free_x')

ggsave(pheat, 
  filename = "24hpf_contrast_control_vs_lowpH_hetmap.png", path = path, 
  width = 8, height = 8)


# (OMIT) Multiple contrast ----

# Filter data for multiple contrast

# mtd <- mtd[mtd$hpf == 110,]

# mtd <- mtd[mtd$hpf == 24,]

# head(count0 <- count0[names(count0) %in% mtd$colName])

colors_fc <- c("red2",  "#4169E1", "forestgreen", "grey30")
  
colData <- mtd %>% arrange(match(colName, names(count0)))

# colData <- colData %>% mutate(g = paste(pH,hpf, sep = ''))

colData <- mutate_if(colData, is.character, as.factor)
colData <- mutate_if(colData, is.double, as.factor)

# using relevel, just specifying the reference level:

colData$pH <- relevel(colData$pH, ref = "Control")

levels(colData$pH)

count0 <- round(count0)

table(colData$pH)

ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = count0,
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
