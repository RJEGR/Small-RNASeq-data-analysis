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
sample_dist = dist(sample_cor, method='euclidean')
hc_samples = hclust(sample_dist, method='complete')

hc_order <- hc_samples$labels[hc_samples$order]

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

recode_to <- c(`24` = "24 hpf", `110` = "110 hpf")

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


