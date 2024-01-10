# EXPLORE DATA TRANSFORMATION AND NORMALIZATION
# GENERATE HEATMAP OF UP/DOWN NORMALIZED DATA:

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

path <- "~/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out/"

count_f <- list.files(path = path, pattern = "Counts.txt", full.names = T)

DB_f <- list.files(path = path, pattern = "RNA_LOCATION_DB.tsv", full.names = T)

MTD_f <- list.files(path = path, pattern = "METADATA.tsv", full.names = T)

# wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/FUNCTIONAL_MIR_ANNOT/"

# SRNA2GO <- read_tsv(paste0(wd, "SRNA2GO.tsv")) # "SRNA_REGULATORY_FUNCTION_DB.tsv"

.colData <- read_tsv(MTD_f)

DB <- read_tsv(DB_f)

str(query.ids <- DB %>% filter(SRNAtype == "miR") %>% distinct(Name) %>% pull())

COUNTS <- read_tsv(count_f)

colNames <- gsub(".clean.newid.subset", "", names(COUNTS))

colnames(COUNTS) <- colNames

# OPTIONAL:

COUNTS <- COUNTS %>% filter(Name %in% query.ids)

# THEN:

rowNames <- COUNTS$Name

COUNTS <- COUNTS %>% select_if(is.double) %>% as(., "matrix")

rownames(COUNTS) <- rowNames

.COUNTS <- COUNTS

# rowSums(COUNTS) %>% as_tibble(rownames = "Name")

# DB %>% select(Name, Locus_type, MajorRNA, SRNAtype, Reads) %>% filter(SRNAtype == "miR")

# 1) Filter data by removing low-abundance genes ----

by_count <- 1; by_freq <- 2

keep <- rowSums(COUNTS > by_count) >= by_freq

sum(keep) # N transcripts

nrow(COUNTS <- COUNTS[keep,])

COUNTS <- round(COUNTS)

# DESeq2::varianceStabilizingTransformation(COUNTS)

library(DESeq2)

ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = COUNTS,
  colData = .colData,
  design = ~ 1 )

dds <- estimateSizeFactors(ddsFullCountTable) 

dds <- estimateDispersions(dds)

# However for other downstream analyses – e.g. for visualization or clustering – it might be useful to work with transformed versions of the count data. Ex:

vst <- DESeq2::varianceStabilizingTransformation(dds) # vst if cols > 10

# vst <- DESeq2::vst(dds) # vst if cols > 10

ntr <- DESeq2::normTransform(dds)

DESeq2::plotPCA(ntr, intgroup = "pH")
DESeq2::plotPCA(vst, intgroup = "pH")

raw_df <- vsn::meanSdPlot(assay(dds), plot = F)

vst_df <- vsn::meanSdPlot(assay(vst), plot = F)
ntr_df <- vsn::meanSdPlot(assay(ntr), plot = F)

rbind(data.frame(py = vst_df$sd, px = vst_df$rank, col = "vst"),
  data.frame(py = ntr_df$sd, px = ntr_df$rank, col = "ntr"),
  data.frame(py = raw_df$sd, px = raw_df$rank, col = "raw")) %>%
  filter(col != "raw") %>%
  ggplot(aes(px, py, color = col)) +
  labs(x = "Ranks", y = "sd", color = "") +
  geom_line(orientation = NA, position = position_identity(), size = 2) +
  theme_bw(base_family = "GillSans", base_size = 20) +
  theme(legend.position = "top")


# USING VST FOR HEATMAP

head(DATA <- assay(vst))

# heatmap(log2(DATA))

# 


sample_cor = cor(DATA, method='pearson', use='pairwise.complete.obs')

# sample_dist = dist(sample_cor, method='euclidean')
# hc_samples = hclust(sample_dist, method='complete')

# plot(hc_samples)

# 
# hc_order <- hc_samples$labels[hc_samples$order]

# hc_samples <- as.hclust(reorder(as.dendrogram(hc_samples), 12:1))

h <- heatmap(sample_cor, col = cm.colors(12), keep.dendro = T)

hc_samples <- as.hclust(h$Rowv)

hc_order <- hc_samples$labels[h$rowInd]

sample_cor %>% 
  as_tibble(rownames = 'LIBRARY_ID') %>%
  pivot_longer(cols = colnames(sample_cor), values_to = 'cor') %>%
  left_join(.colData) %>% 
  mutate(LIBRARY_ID = factor(LIBRARY_ID, levels = rev(hc_order))) -> sample_cor_long

library(ggh4x)

P <- sample_cor_long %>%
  ggplot(aes(x = LIBRARY_ID, y = name, fill = cor)) +  
  geom_tile(linewidth = 0.2) +
  ggsci::scale_fill_material(name = "", "blue-grey") +
  ggsci::scale_color_material(name = "", "blue-grey") +
  scale_x_discrete(position = 'bottom') +
  ggh4x::scale_y_dendrogram(hclust = hc_samples, position = "left", labels = NULL) +
  guides(y.sec = guide_axis_manual(labels = hc_order, label_size = 5, label_family = "GillSans")) +
  # ggh4x::scale_x_dendrogram(hclust = hc_samples, position = "top", labels = NULL) +
  theme_bw(base_size = 7, base_family = "GillSans") +
  labs(x = '', y = '') +
  theme(
    legend.position = "bottom",
    # axis.text.x = element_blank(),
    # axis.ticks.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 5),
    strip.background = element_rect(fill = 'grey89', color = 'white'),
    panel.border = element_blank(),
    plot.title = element_text(hjust = 0),
    plot.caption = element_text(hjust = 0),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank()) 

P <- P + guides(
  fill = guide_colorbar(barwidth = unit(1.5, "in"),
    barheight = unit(0.05, "in"), label.position = "bottom",
    alignd = 0.5,
    ticks.colour = "black", ticks.linewidth = 0.5,
    frame.colour = "black", frame.linewidth = 0.5,
    label.theme = element_text(family = "GillSans", size = 7)))

P

# TOP

recode_to <- c(`110` = "110 hpf", `24` = "24 hpf")

TOPDF <- sample_cor_long %>%
  distinct(LIBRARY_ID, hpf, pH) %>%
  dplyr::mutate(hpf = dplyr::recode_factor(hpf, !!!recode_to)) %>%
  mutate(label = ifelse(pH %in% "Low", "*", "")) %>%
  mutate(y = 1)

topplot <- TOPDF %>%
  ggplot(aes(y = y, x = LIBRARY_ID, color = as.factor(hpf))) +
  geom_point(shape = 15, size = 2) +
  geom_text(aes(label = label),  vjust = -0.7, hjust = 0.5, size = 1.5, family =  "GillSans", color = "#d73027") +
  ggh4x::scale_x_dendrogram(hclust = hc_samples, position = 'top', labels = NULL) +
  # ggh4x::guide_dendro()
  # guides(x.sec = guide_axis_manual(labels = hc_order, label_size = 3.5)) +
  theme_bw(base_family = "GillSans", base_size = 10) +
  # see::scale_color_pizza(name = "", reverse = T) +
  scale_color_manual("", values = c("#DADADA", "#D4DBC2")) + # "#4575b4", "#d73027"
  theme(legend.position = 'top',
    panel.border = element_blank(), 
    plot.background = element_rect(fill='transparent', color = 'transparent'),
    plot.margin = unit(c(0,0,0,0), "pt"),
    panel.grid.minor = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    panel.grid.major = element_blank()) 


library(patchwork)

psave <- topplot/ plot_spacer() /P + plot_layout(heights = c(0.6, -0.5, 5))

psave

ggsave(psave, filename = 'SAMPLE_HEATMAP_MIR.png', path = path, width = 2.5, height = 3.5, device = png, dpi = 300)

# PCA
vst <- DESeq2::varianceStabilizingTransformation(dds) # vst if cols > 10

ncol(data <- assay(vst))

PCA = prcomp(t(data), center = T, scale. = FALSE)

percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)

sd_ratio <- sqrt(percentVar[2] / percentVar[1])

PCAdf <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2])

PCAdf %>%
  mutate(LIBRARY_ID = rownames(.)) %>%
  left_join(.colData) %>%
  ggplot(., aes(PC1, PC2)) +
  # coord_fixed(ratio = sd_ratio) +
  geom_abline(slope = 0, intercept = 0, linetype="dashed", alpha=0.5) +
  geom_vline(xintercept = 0, linetype="dashed", alpha=0.5) +
  ggforce::geom_mark_ellipse(aes(group = as.factor(hpf), label = as.factor(hpf)),
    fill = 'grey', colour = NA) +
  geom_point(size = 7, alpha = 0.7, aes(color = pH)) +
  geom_text( family = "GillSans",
    mapping = aes(label = pH), size = 2.5) +
  labs(caption = '', color = "") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme_minimal(base_family = "GillSans", base_size = 18) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = 'none')


# heatmap of z-score ====

# including only degs

wd <- "~/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out/"

RES.P <- read_tsv(paste0(wd, "DESEQ_RES.tsv")) %>% filter( padj < 0.05 & abs(log2FoldChange) > 1) # abs(log2FoldChange) > 2 &

RES.P <- RES.P %>% filter(CONTRAST %in% c("CONTRAST_A","CONTRAST_B")) # ,  "CONTRAST_D"

str(query.names <- RES.P %>% distinct(Name) %>% pull())

dim(M <- assay(vst))

datExpr <- M %>% as_tibble(rownames = 'Name') %>%
  right_join(distinct(RES.P, Name, Family), by = "Name") %>%
  group_by(Family) %>%
  summarise_at(vars(colnames(M)), sum)

# BY SAMPLE

# datExpr <- datExpr %>% 
#   pivot_longer(-Family, names_to = "LIBRARY_ID") %>%
#   left_join(.colData, by = "LIBRARY_ID") %>%
#   group_by(Family, colName) %>% summarise(value = mean(value)) %>%
#   pivot_wider(names_from = colName, values_from = value) %>%
#   ungroup()

which_cols <- datExpr %>% dplyr::select_if(is.double) %>% colnames()

M <- datExpr %>% dplyr::select(all_of(which_cols))

M <- as(M, "matrix")

rownames(M) <- datExpr$Family

head(M)

heatmap(t(apply(M, 1, z_scores)))

# 
# plot(density(M[1,]), type = "b", pch = 19)
# 
# abline(v = mean(M[1,]), col="red", lwd=3, lty=2)
# 
# abline(col="blue", lwd=3, lty=2, v = mean(M[1,])+sd(M[1,]))
# abline(col="blue", lwd=3, lty=2, v = mean(M[1,])-sd(M[1,]))
# 
# hist(apply(t(M), 1, z_scores)[1,])

# keep_mirs <- rownames(M) %in% query.names

# dim(M <- M[keep_mirs,])

# head(M <- t(scale(t(M))))

z_scores <- function(x) {(x-mean(x))/sd(x)}

heatmap(M <- apply(M, 1, z_scores))


# head(M1 <- apply(M[!grepl("24", rownames(M)),], 2, z_scores))
# head(M2 <- apply(M[grepl("24", rownames(M)),], 2, z_scores))
# 
# 
# dim(M <- rbind(M1, M2))

# head(M1 <- t(scale(t(M[,!grepl("24", colnames(M))]))))
# head(M2 <- t(scale(t(M[,grepl("24", colnames(M))]))))
# 
# dim(M <- cbind(M1, M2))

# scale to zscore by groups

# head(M1 <- t(scale(t(M[,!grepl("24", colnames(M))]))))
# head(M2 <- t(scale(t(M[,grepl("24", colnames(M))]))))
# dim(M <- cbind(M1, M2))
# M[is.na(M)] <- 0

h <- heatmap(t(M), col = cm.colors(12), keep.dendro = T)

hc_samples <- as.hclust(h$Colv)
plot(hc_samples)
hc_order <- hc_samples$labels[h$colInd]

hc_mirs <- as.hclust(h$Rowv)
plot(hc_mirs)
order_mirs <- hc_mirs$labels[h$rowInd]

summary(M)


DATA <- M %>% 
  as_tibble(rownames = 'LIBRARY_ID') %>%
  pivot_longer(cols = colnames(M), values_to = 'fill', names_to = "Name") %>%
  left_join(.colData, by = "LIBRARY_ID") %>%
  mutate(LIBRARY_ID = factor(LIBRARY_ID, levels = rev(hc_order)))


DATA %>%
  dplyr::mutate(pH = dplyr::recode_factor(pH, !!!c(`Low` = "pH 7.6",`Control` = "pH 8.0"))) %>%
  dplyr::mutate(hpf = dplyr::recode_factor(hpf, !!!c(`24` = "24 hpf", `110` = "110 hpf"))) %>%
  ggplot(aes(x = LIBRARY_ID, y = Name, fill = fill)) +  
  geom_tile(linewidth = 0.2) +
  # facet_grid(~ hpf+pH, scales = "free_x") +
  ggsci::scale_fill_gsea(name = "", reverse = T) +
  ggsci::scale_color_gsea(name = "", reverse = T) +
  # scale_x_discrete(position = 'bottom') +
  # ggh4x::scale_x_dendrogram(hclust = hc_samples, position = "top", labels = NULL) 
  # guides(x.sec = guide_axis_manual(labels = hc_order, label_size = 5, label_family = "GillSans")) 
  ggh4x::scale_y_dendrogram(hclust = hc_mirs, position = "left", labels = NULL) +
  guides(y.sec = guide_axis_manual(labels = order_mirs, label_size = 5, label_family = "GillSans")) +
  theme_bw(base_size = 7, base_family = "GillSans") +
  labs(x = '', y = '') +
  theme(
    legend.position = "bottom",
    # axis.text.x = element_blank(),
    # axis.ticks.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 5),
    strip.background = element_rect(fill = 'grey89', color = 'white'),
    panel.border = element_blank(),
    plot.title = element_text(hjust = 0),
    plot.caption = element_text(hjust = 0),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank()) -> P

P <- P + ggh4x::facet_nested(~ hpf+pH, nest_line = T, scales = "free_x", space = "free_x") +
  theme(strip.background = element_rect(fill = 'white', color = 'white')) 

P <- P + guides(
  fill = guide_colorbar(barwidth = unit(1.5, "in"),
    barheight = unit(0.05, "in"), label.position = "bottom",
    alignd = 0.5,
    ticks.colour = "black", ticks.linewidth = 0.5,
    frame.colour = "black", frame.linewidth = 0.5,
    label.theme = element_text(family = "GillSans", size = 7)))

P

recode_to <- c(`24` = "24 hpf", `110` = "110 hpf")

TOPDF <- DATA %>%
  distinct(LIBRARY_ID, hpf, pH) %>%
  dplyr::mutate(hpf = dplyr::recode_factor(hpf, !!!recode_to)) %>%
  mutate(label = ifelse(pH %in% "Low", "*", "")) %>%
  mutate(y = 1)

topplot <- TOPDF %>%
  ggplot(aes(y = y, x = LIBRARY_ID)) + # color = as.factor(hpf)
  # geom_point(shape = 15, size = 2) +
  # geom_text(aes(label = label),  vjust = -0.7, hjust = 0.5, size = 1.5, family =  "GillSans", color = "#d73027") +
  ggh4x::scale_x_dendrogram(hclust = hc_samples, position = 'top', labels = NULL) +
  # ggh4x::guide_dendro()
  # guides(x.sec = guide_axis_manual(labels = hc_order, label_size = 3.5)) +
  theme_bw(base_family = "GillSans", base_size = 10) +
  # see::scale_color_pizza(name = "", reverse = T) +
  scale_color_manual("", values = c("#DADADA", "#D4DBC2")) + # "#4575b4", "#d73027"
  theme(legend.position = 'top',
    panel.border = element_blank(), 
    plot.background = element_rect(fill='transparent', color = 'transparent'),
    plot.margin = unit(c(0,0,0,0), "pt"),
    panel.grid.minor = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    panel.grid.major = element_blank()) 


library(patchwork)

psave <- topplot/ plot_spacer() /P + plot_layout(heights = c(0.5, -0.5, 5))

# psave

ggsave(psave, filename = 'DESEQ2_OADEGS_HEATMAP_MIR.png', path = path, width = 3, height = 3.5, device = png, dpi = 300)

# by log2fc
.RES.P <- read_tsv(paste0(wd, "DESEQ_RES.tsv")) %>% 
  filter(CONTRAST %in% c("CONTRAST_A","CONTRAST_B"))
  # filter( padj < 0.05 & abs(log2FoldChange) > 1)

subset <- .RES.P %>% filter(Name %in% query.names) %>%
  dplyr::mutate(CONTRAST = dplyr::recode_factor(CONTRAST, !!!c(`CONTRAST_A` = "24 hpf",
    `CONTRAST_B` = "110 hpf"))) %>%
  mutate(cc = ifelse(padj < 0.05, "p-value","")) 

.RES.P %>%
  dplyr::mutate(CONTRAST = dplyr::recode_factor(CONTRAST, !!!c(`CONTRAST_A` = "24 hpf",
    `CONTRAST_B` = "110 hpf"))) %>%
  mutate(pHcc = ifelse(log2FoldChange < 0, "pH 7.6","pH 8.0")) %>%
  ggplot(aes(x = log10(baseMean), y = log2FoldChange)) + 
  facet_grid(~ CONTRAST) +
  geom_abline(slope = 0, intercept = 0, linetype="dashed", alpha=0.5) +
  # geom_abline(slope = 0, intercept = -1, linetype="dashed", alpha=0.5) +
  # geom_abline(slope = 0, intercept = 1, linetype="dashed", alpha=0.5) +
  # geom_rect( aes(xmin=-Inf, xmax = Inf, ymin = 1, ymax = Inf), fill = "grey89") +
  # geom_rect( aes(xmin=-Inf, xmax = Inf, ymin = -1, ymax = -Inf), fill = "grey89") +
  # annotate("text", x = 5, y = 20, label = "pH 8.0", family = "GillSans", color = "black") +
  # annotate("text", x = 5, y = -10, label = "pH 7.6", family = "GillSans",  color = "black") +
  geom_point(shape = 1, alpha = 1, aes(color = pHcc)) +
  ggrepel::geom_text_repel(data = subset, aes(label = Family), size = 2.5, family = "GillSans", max.overlaps = 50) +
  scale_color_manual("", values = c(`pH 8.0` = "#4169E1", `pH 7.6` = "#d73027")) +
  theme_bw(base_family = "GillSans", base_size = 12) +
  theme(legend.position = "top",
    panel.border = element_blank(),
    strip.background = element_rect(fill = 'grey89', color = 'white'),
    plot.title = element_text(hjust = 0),
    plot.caption = element_text(hjust = 0),
    panel.grid.minor.y = element_blank(),
    # panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    strip.background.y = element_blank(),
    # axis.text.y.right = element_text(angle = 0, hjust = 1, vjust = 0, size = 2.5),
    # axis.text.y = element_text(angle = 0, size = 5),
    axis.text.x = element_text(angle = 0)) -> p

ggsave(p, filename = 'DESEQ2_MEANvs_log2fc_MIR.png', path = path, width = 6, height = 4.5, device = png, dpi = 300)
