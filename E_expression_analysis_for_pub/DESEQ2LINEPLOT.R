# GENERATE AVERAGE-LINE-PLOT USING DIFF EXP. MIRS AS INPUTS
# FOCUS ON MIRS UP-EXPRESSED UNDER LOW PH (i.e CONTRAST A AND B)
# OR ONLY IN EXCLUSIVE MIRS
# COMPARISON FROM CONTRAST C AND D SHOW HOW REGULATORY NETWORK DURING DEVELOPMENT (AT 110 HPF) CHANGES IN REPONSE TO LOW PH

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

# wd <- "~/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out/"
# 

wd <-  "~/Documents/MIRNA_HALIOTIS/MIRNA_PUB_2024/res"


dds <- read_rds(paste0(wd, "/SEQUENCES_MERGED_DDS_DESEQ2.rds"))

# RES.P <- read_tsv(paste0(wd, "DESEQ_RES_P.tsv")) %>% filter( padj < 0.05) 


RES.P <- read_tsv(list.files(path = wd, pattern = "DESEQ_RES.tsv", full.names = T)) %>% 
  filter( padj < 0.05  & abs(log2FoldChange) > 1) %>%
  mutate(star = ifelse(padj <.001, "***", 
    ifelse(padj <.01, "**",
      ifelse(padj <.05, "*", "")))) 

RES.P %>% dplyr::count(CONTRAST)

.colData <- as_tibble(SummarizedExperiment::colData(dds)) # %>% dplyr::select(contains(c("LIBRARY_ID","hpf", "pH")))


# WHICH MIRS REGARDING ACIDIFICATION
WHICH_CONT <- c("CONTRAST_A", "CONTRAST_B")

query.genes <- RES.P %>% filter(CONTRAST %in% WHICH_CONT) %>% pull(MajorRNA)

sort.genes <- RES.P %>% filter(CONTRAST %in% WHICH_CONT) %>% pull(Name)
  # group_by(CONTRAST) %>% arrange(log2FoldChange) %>% pull(Name)

cnts <- DESeq2::counts(dds, normalized = T, replaced = F)[query.genes,]

# heatmap(cnts)

z_scores <- function(x) {(x-mean(x))/sd(x)}

heatmap(cnts <- apply(cnts, 1, z_scores))


DF <- t(cnts) %>% as_tibble(rownames = "MajorRNA") %>% 
  pivot_longer(-MajorRNA, names_to = "LIBRARY_ID", values_to = "count")

DF <- RES.P %>% distinct(Name, MajorRNA) %>% 
  right_join(DF, by = "MajorRNA") %>% 
  left_join(.colData, by = "LIBRARY_ID") %>%
  mutate(Name = factor(Name, levels = unique(sort.genes)))
  # drop_na(any_of(intgroup)) %>%
  # mutate(CONTRAST = intgroup)

# names(DF)[names(DF) %in% intgroup] <- "Design"

fun.data.trend <- "mean_se" # "mean_cl_boot", "mean_sdl"


DF %>%
  # arrange(count) %>% 
  # group_by(Family, LIBRARY_ID, Design, CONTRAST) %>% sample_n(1) %>%
  # dplyr::mutate(CONTRAST = dplyr::recode(CONTRAST, !!!recode_to)) %>%
  # dplyr::mutate(Design = dplyr::recode(Design, !!!recode_to_AB)) %>%
  # mutate(CONTRAST = factor(CONTRAST, levels = recode_to[1:2])) %>%
  # mutate(Design = factor(Design, levels = recode_to_AB)) %>%
  # mutate(Name = factor(Name, levels = unique(names(LOGFC_LABELLER)))) %>%
  ggplot(aes(x = as.factor(pH), y=count, group = as.factor(hpf), color = as.factor(hpf))) +
  facet_wrap(~ Name) +
  # facet_wrap(CONTRAST ~ Family, labeller = labeller(Family = LOGFC_LABELLER), scales = "free_y") +
  # ggh4x::facet_nested_wrap(~ CONTRAST+Family, labeller = labeller(Family = LOGFC_LABELLER), scales = "free_y", nest_line = F) +
  geom_jitter(position=position_jitter(w=0.1,h=0), size = 1, alpha = 0.5) +
  stat_summary(fun.data = fun.data.trend, linewidth = 0.7, size = 0.7, alpha = 0.7) +
  stat_summary(fun = mean, geom = "line") +
  # labs(y = ylab, x = "") +
  theme_bw(base_family = "GillSans", base_size = 11) +
  theme(strip.background = element_rect(fill = 'grey89', color = 'white'),
    panel.border = element_blank(),
    plot.title = element_text(hjust = 0),
    plot.caption = element_text(hjust = 0),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    # axis.text.y.right = element_text(angle = 0, hjust = 1, vjust = 0, size = 2.5),
    axis.text.y = element_text(angle = 0, size = 7),
    axis.text.x = element_text(angle = 0, size = 10)) -> p


p
# 
ggsave(p, filename = 'DESEQ2UPEXPRESSED_MIRS_IN_REPONSE_TO_LOW_PH.png', path = wd, width = 7, height = 5, device = png, dpi = 300)

# exit



h <- heatmap(t(cnts), col = cm.colors(12), keep.dendro = T)

hc_samples <- as.hclust(h$Colv)
plot(hc_samples)
hc_order <- hc_samples$labels[h$colInd]

hc_mirs <- as.hclust(h$Rowv)
plot(hc_mirs)
order_mirs <- hc_mirs$labels[h$rowInd]


DATA <- t(cnts) %>% 
  as_tibble(rownames = 'MajorRNA') %>%
  pivot_longer(-MajorRNA, values_to = 'fill', names_to = "LIBRARY_ID") %>%
  left_join(.colData, by = "LIBRARY_ID") %>%
  mutate(LIBRARY_ID = factor(LIBRARY_ID, levels = rev(hc_order))) %>%
  mutate(MajorRNA = factor(LIBRARY_ID, levels = unique(order_mirs))) 

library(ggh4x)

DATA %>%
  dplyr::mutate(pH = dplyr::recode_factor(pH, !!!c(`Low` = "pH 7.6",`Control` = "pH 8.0"))) %>%
  dplyr::mutate(hpf = dplyr::recode_factor(hpf, !!!c(`24` = "24 hpf", `110` = "110 hpf"))) %>%
  ggplot(aes(x = LIBRARY_ID, y = MajorRNA, fill = fill)) +  
  geom_tile(linewidth = 0.2) +
  # facet_grid(~ hpf+pH, scales = "free_x") +
  ggsci::scale_fill_gsea(name = "", reverse = T) +
  ggsci::scale_color_gsea(name = "", reverse = T) +
  # scale_x_discrete(position = 'bottom') +
  # ggh4x::scale_x_dendrogram(hclust = hc_samples, position = "top", labels = NULL) 
  # guides(x.sec = guide_axis_manual(labels = hc_order, label_size = 5, label_family = "GillSans")) 
  # ggh4x::scale_y_dendrogram(hclust = hc_mirs, position = "left", labels = NULL) +
  # guides(y.sec = guide_axis_manual(labels = order_mirs, label_size = 5, label_family = "GillSans")) +
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

