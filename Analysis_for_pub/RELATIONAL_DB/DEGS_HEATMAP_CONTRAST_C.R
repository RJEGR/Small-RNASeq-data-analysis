# Plot z-score heatmap for vital process for 24 and 110 hpf
# group 1st COGs/GOs by modules
# split modules  24 from 110 hpf
# plot heatmap
# This code is similar from DEGS_HEATMAP_ZSCORE



rm(list = ls());

if(!is.null(dev.list())) dev.off()


library(DESeq2)

library(tidyverse)


dir <-  "~/Documents/MIRNA_HALIOTIS/MIRNA_PUB_2024"

print(DB <- read_tsv(file.path(dir, "LONGER_RELATIONAL_DB.tsv")))

dds <- read_rds(list.files(path = dir, pattern = "DDS_DESEQ2.rds", full.names = T))

colData <- data.frame(colData(dds))

WHICH_CONTRAST <- c("CONTRAST_C") # "CONTRAST_D"

RES <- read_rds(paste0(dir, "/SEQUENCES_MERGED_DESEQ_RES_WGCNA.rds")) %>%
  mutate(Contrast = ifelse(sign(log2FoldChange) == -1, sampleB, sampleA)) %>%
  mutate(Contrast = paste0(CONTRAST,"_",Contrast))

RES.P <- RES %>% 
  filter( padj < 0.05  & abs(log2FoldChange) > 1) %>%
  filter(CONTRAST %in% WHICH_CONTRAST)

# Evals if both, or upexpr in contrast C and D

which_contrast <- function(x) { 
  x <- x[!is.na(x)] 
  n <- length(unique(x))
  x <- unique(x)
  
  if(n > 1) {
    x <- "BOTH"
  } else
    x <- paste(x, sep = '|', collapse = '|') }

RES.P %>%
  mutate(HPF = ifelse(sign(log2FoldChange) == -1, "110 hpf", "24 hpf")) %>%
  group_by(MajorRNA, HPF) %>%
  summarise(across(CONTRAST, .fns = which_contrast)) %>%
  ungroup() %>% count(CONTRAST)

which_mirs <- RES.P %>% 
  # filter(CONTRAST_DE %in% WHICH_CONTRAST) %>%
  distinct(MajorRNA) %>% pull()

str(which_mirs)  


# divide the counts by the size factors or normalization factors before
barplot(cnts <- DESeq2::counts(dds, normalized = F, replaced = F)[which_mirs,])
barplot(cnts <- DESeq2::varianceStabilizingTransformation(round(cnts)))

na_fill <- function(x) { ifelse(is.na(x), 0, x)} 

z_scores <- function(x) {(x-mean(x))/sd(x)}

# c1 <- cnts[,grepl("HR24", colnames(cnts))]
# c2 <- cnts[,!grepl("HR24", colnames(cnts))]

colnames(c1 <- cnts[,grepl("8[1-3]", colnames(cnts))])
colnames(c2 <- cnts[,!grepl("8[1-3]", colnames(cnts))])

M1 <- apply(c1, 1, z_scores)
M2 <- apply(c2, 1, z_scores)

# dim(M <- rbind(M1, M2))

dim(M <- M1)

# M <- apply(M, 1, na_fill)

# heatmap(M <- apply(cnts, 2, z_scores))

h <- heatmap(t(M), keep.dendro = T)


hc_samples <- as.hclust(h$Colv)
plot(hc_samples)
hc_order <- hc_samples$labels[h$colInd]

hc_mirs <- as.hclust(h$Rowv)
# plot(hc_mirs)
order_mirs <- hc_mirs$labels[h$rowInd]



# DATA <- M %>% 
#   as_tibble(rownames = 'MajorRNA') %>%
#   pivot_longer(cols = colnames(M), values_to = 'fill', names_to = "LIBRARY_ID") %>%
#   left_join(colData, by = "LIBRARY_ID") %>%
#   mutate(LIBRARY_ID = factor(LIBRARY_ID, levels = rev(hc_order)))

DATA <- M %>% 
  as_tibble(rownames = 'LIBRARY_ID') %>%
  pivot_longer(cols = colnames(M), values_to = 'fill', names_to = "MajorRNA") %>%
  left_join(colData, by = "LIBRARY_ID") %>%
  mutate(LIBRARY_ID = factor(LIBRARY_ID, levels = rev(hc_order))) %>%
  mutate(MajorRNA = factor(MajorRNA, levels = hc_mirs$labels[h$rowInd]))
  

labels <- RES.P %>% distinct(MajorRNA, MirGeneDB_ID)

labels <- structure(labels$MirGeneDB_ID, names = labels$MajorRNA)

labels <- labels[match(order_mirs, names(labels))]

identical(order_mirs, names(labels))

order_mirs <- labels

# Add cogs

COG_names <- DB %>%
  filter(MajorRNA %in% names(order_mirs)) %>%
  drop_na(COG_name) %>% 
  dplyr::count(MajorRNA, COG_name, sort = T) %>%
  group_by(MajorRNA) %>% slice_head(n = 1) %>%
  pull(COG_name, name = MajorRNA)

COG_names <- COG_names[match(names(order_mirs), names(COG_names))]

identical(names(COG_names), names(order_mirs)) # TRUE

order_mirs <- structure(paste0(order_mirs, " (", COG_names, ")"), names = names(order_mirs))

# Then plot

lo = floor(min(DATA$fill))
up = ceiling(max(DATA$fill))
mid = (lo + up)/2

library(ggh4x)

DATA %>%
  dplyr::mutate(pH = dplyr::recode_factor(pH, !!!c(`Control` = "pH 8.0", `Low` = "pH 7.6"))) %>%
  dplyr::mutate(hpf = dplyr::recode_factor(hpf, !!!c(`24` = "24 hpf", `110` = "110 hpf"))) %>%
  mutate(LIBRARY_ID = gsub("HR", "", LIBRARY_ID)) %>%
  ggplot(aes(x = LIBRARY_ID, y = MajorRNA, fill = fill)) +  
  geom_tile(color = 'white', linewidth = 0.5, width = 1) +
  # facet_grid(~ hpf+pH, scales = "free_x") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
    na.value = "white", midpoint = mid, limit = c(lo, up),
    name = NULL) +
  scale_y_discrete(labels = order_mirs, position = "left") +
  # ggh4x::scale_y_dendrogram(hclust = hc_mirs, position = NULL, labels = NULL) +
  # guides(y.sec = guide_axis_manual(labels = order_mirs, label_size = 8, label_family = "GillSans")) +
  theme_bw(base_size = 10, base_family = "GillSans") +
  labs(x = '', y = '') +
  theme(
    legend.position = "top",
    axis.text.y = element_text(size = 2.5),
    # axis.ticks.y = element_blank(),
    axis.ticks.y = element_line(size=0.05),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5),
    axis.ticks.x = element_line(size=0.05),
    # strip.background = element_rect(fill = 'grey89', color = 'white'),
    panel.border = element_blank(),
    plot.title = element_text(hjust = 0),
    plot.caption = element_text(hjust = 0),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank()) -> P

P <- P + 
  ggh4x::facet_nested(~ pH+hpf, nest_line = T, scales = "free", space = "free") +
  theme(
    strip.text.x = element_text(angle = 0, hjust = 1),
    strip.placement = "outside",
    strip.background = element_rect(fill = 'white', color = 'white')) 

P <- P + guides(
  fill = guide_colorbar(barwidth = unit(1.5, "in"),
    barheight = unit(0.05, "in"), 
    label.position = "top",
    title = "Z-score",
    title.position  = "top",
    title.theme = element_text(size = 7, family = "GillSans", hjust = 1),
    alignd = 0.8,
    ticks.colour = "black", ticks.linewidth = 0.25,
    frame.colour = "black", frame.linewidth = 0.25,
    label.theme = element_text(size = 7, family = "GillSans")))
# P

TOPDF <- RES.P %>% 
  mutate(HPF = ifelse(sign(log2FoldChange) == -1, "110 hpf", "24 hpf")) %>%
  group_by(MajorRNA, HPF, WGCNA) %>%
  summarise(across(CONTRAST, .fns = which_contrast)) %>%
  # distinct(MajorRNA, WGCNA, HPF) %>%
  mutate(shape = ifelse(CONTRAST == "CONTRAST_D", "pH 7.6", "")) %>%
  mutate(y = 0.5) 

# TOPDF %>% group_by(CONTRAST, HPF) %>% count()

col <- RES.P %>% distinct(WGCNA) %>% pull(WGCNA, name = WGCNA)

LEFTP <- TOPDF %>%
  mutate(HPF = factor(HPF)) %>%
  ggplot(aes(x = y, y = MajorRNA, color = WGCNA)) + # color = as.factor(hpf)
  # facet_grid(~ HPF, scales = "free", space = "free") +
  ggh4x::scale_y_dendrogram(hclust = hc_mirs, position = 'left', labels = NULL) +
  # guides(y.sec = guide_axis_manual(labels = order_mirs, 
  #   label_size = 2.5, label_family = "GillSans")) +
  geom_point(aes(shape = shape), size =0.7) +
  scale_shape_manual(values = c(15,19, 8)) +
  theme_bw(base_family = "GillSans", base_size = 10) +
  # scale_color_manual(values = c("grey79", "black")) +
  scale_color_manual(values = col) +
  guides(color=guide_legend(title = "", nrow = 1), 
    shape = guide_legend(title = "", nrow = 1)) +
  theme(legend.position = 'none',
    panel.border = element_blank(), 
    plot.background = element_rect(fill='transparent', color = 'transparent'),
    panel.grid.minor = element_blank(),
    axis.ticks.x = element_blank(),
    # axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.title = element_blank(),
    panel.grid.major = element_blank()) 

library(patchwork)

ps <- LEFTP + plot_spacer() + P + 
  plot_layout(nrow = 1, widths = c(1, -0.9999, 7)) +
  theme(plot.margin = unit(c(.2,.2,.2,0), "cm"))

ggsave(ps, filename = 'DEGS_HEATMAP_CONTRAST_C.png', path = dir, width = 3.5, height = 5, device = png, dpi = 500)
