
# LOAD DATA RESULTS FROM DESEQ ANALYSIS
# 1) FILTER DATA BY PADJ < 0.05 (NOT USE LOG2FC)
# 2) GENERATE BARPLOT OF LOG2FC W/ SIGNIFICANCE FOR KNOWN AND NOVEL MIRS
# 3) SPLIT EXCLUSIVE FROM INTERSECTED MIRS
# 4) SAVE 

# NOTE:

# POSITIVE LFC == UP EXPRESSED IN EXPERIMENTAL (OR sampleB)
# NEGATIVE LFC == UP EXPRESSED IN CONTROL (OR sampleA)
# view baseMeanA and baseMeanB to contrast the expression values

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

wd <-  "~/Documents/MIRNA_HALIOTIS/MIRNA_PUB_2024"

RES <- read_tsv(list.files(path = wd, pattern = "DESEQ_RES.tsv", full.names = T)) # %>% mutate(log2FoldChange = log2FoldChange * -1)

url <- "https://raw.githubusercontent.com/RJEGR/Cancer_sete_T_assembly/main/functions.R"

source(url)

CONTRAST <- RES %>% distinct(CONTRAST) %>% pull()

recode_to <- c("24 HPF | Ctrl : Low pH", "110 HPF | Ctrl : Low pH", "Ctrl pH | 24 : 110 HPF", "Low pH | 24 : 110 HPF")

recode_to <- structure(recode_to, names = CONTRAST)

recode_fc <- structure(c("Up","Down"), names = c(1,-1))

RES.P <- RES %>% filter( padj < 0.05 & abs(log2FoldChange) > 1) # abs(log2FoldChange) > 2 &

RES.P %>% distinct(Name)

RES.P <- RES.P %>%
  mutate(star = ifelse(padj <.001, "***", 
    ifelse(padj <.01, "**",
      ifelse(padj <.05, "*", "")))) %>%
  mutate(SIGN = sign(log2FoldChange)) %>%
  mutate(CONTRAST_DE = CONTRAST) %>%
  dplyr::mutate(CONTRAST = dplyr::recode_factor(CONTRAST, !!!recode_to)) %>%
  dplyr::mutate(SIGN = dplyr::recode_factor(SIGN, !!!recode_fc)) %>%
  separate(CONTRAST, into = c("WRAP","CONTRAST"), sep = "[|]")

# ggsave(p, filename = 'DESEQ2BARPLOT.png', path = wd, width = 5, height = 10, device = png, dpi = 300)

# 2) BY DEVELOPMENT ====

# INCLUDE VOLCANO ====


WHICH_CONTRAST <- c("CONTRAST_C", "CONTRAST_D")

recode_to <-  structure(c("A) pH 8.0", "B) pH 7.6"), names = WHICH_CONTRAST)

RES_CC <- prep_DE_data(RES, alpha = 0.05, lfcThreshold = 1)

colors_fc <- c("forestgreen",  "#4169E1", "red2", "grey")

c("#DADADA", "#D4DBC2")

upplot <- RES_CC %>% 
  filter(CONTRAST %in% WHICH_CONTRAST) %>%
  dplyr::mutate(SIGN = sign(logFC)) %>%
  dplyr::mutate(CONTRAST = dplyr::recode_factor(CONTRAST, !!!recode_to)) %>%
  # filter(cc != "NS") %>%
  mutate(logFC = logFC*-1) %>% # reverse sort of HPF
  ggplot(aes(x = logFC, y = -log10(padj))) +
  facet_grid(~ CONTRAST) +
  geom_rect(
    aes(xmin=-15, xmax = -1, ymin = 1, ymax = Inf), fill = '#DADADA') +
  geom_rect(
    aes(xmin=1, xmax = 15, ymin = 1, ymax = Inf), fill = '#D4DBC2') +
  geom_point(aes(color = cc), alpha = 3/5) +
  scale_color_manual(name = "", values = colors_fc, labels = scales::parse_format()) +
  labs(x= expression(Log[2] ~ "Fold Change"), 
    y = expression(-Log[10] ~ "padj")) +
  theme_bw(base_family = "GillSans", base_size = 10) +
  geom_abline(slope = 0, intercept = -log10(0.05), linetype="dashed", alpha=0.5) +
  geom_vline(xintercept = 1, linetype="dashed", alpha=0.5) +
  geom_vline(xintercept = -1, linetype="dashed", alpha=0.5) +
  annotate("text", x = -10, y = 90, label = "24 hpf", family = "GillSans") +
  annotate("text", x = 10, y = 90, label = "110 hpf", family = "GillSans") 
  # xlim(-15, 30)

upplot <- upplot + theme(legend.position = "top",
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
  axis.text.x = element_text(angle = 0))

# ggsave(upplot, filename = 'DESEQ2VOLCANO_CONTRAST_C_D.png', path = wd, width = 5, height = 3, device = png, dpi = 300)

recode_fc <- structure(c("24 hpf","110 hpf"), names = c(1,-1))

recode_to <-  structure(c("pH 8.0", "pH 7.6"), names = WHICH_CONTRAST)

# RES.P %>% filter(CONTRAST_DE %in% WHICH_CONTRAST) %>% distinct(Name)

UPSETDF <- RES.P %>% 
  filter(CONTRAST_DE %in% WHICH_CONTRAST) %>%
  mutate(SIGN = sign(log2FoldChange)) %>%
  dplyr::mutate(SIGN = dplyr::recode_factor(SIGN, !!!recode_fc)) %>%
  dplyr::mutate(CONTRAST_DE = dplyr::recode_factor(CONTRAST_DE, !!!recode_to)) %>%
  group_by(Name, SIGN) %>%
  summarise(across(CONTRAST_DE, .fns = list), n = n())

# UPSETDF <- UPSETDF %>%
#   ungroup() %>%
#   left_join(distinct(RES.P, Name, Family)) %>%
#   distinct(Family, .keep_all = T)

library(ggupset)

UPSETDF %>%
  mutate(facet = "C) Intersected") %>%
  ggplot(aes(x = CONTRAST_DE, fill = SIGN)) +
  geom_bar(position = position_dodge(width = 1), color = "black", linewidth = 0.2) +
  # geom_segment(aes(xend = CONTRAST_DE, yend = n, y = Inf, color = SIGN), linewidth = 1) +
  geom_text(stat='count', aes(label = after_stat(count)), 
    position = position_dodge(width = 1), vjust = -0.2, family = "GillSans", size = 2.5) +
  ggupset::scale_x_upset(order_by = "degree", reverse = F) +
  theme_bw(base_family = "GillSans", base_size = 10) +
  ggupset::theme_combmatrix(combmatrix.panel.point.color.fill = "black",
    combmatrix.panel.line.size = 0, base_family = "GillSans") +
  axis_combmatrix(levels = c("pH 8.0", "pH 7.6"), clip = "off") +
  labs(x = '', y = 'miRNAs number') +
  scale_color_manual("", values = c("#DADADA", "#D4DBC2")) +
  scale_fill_manual("", values =  c("#DADADA", "#D4DBC2")) +
  guides(fill = guide_legend(title = "", nrow = 1)) +
  ylim(c(0,40)) -> p2

p2 <- p2 + theme(legend.position = "top",
  legend.key.width = unit(0.15, "cm"),
  legend.key.height = unit(0.05, "cm"),
  panel.border = element_blank(),
  plot.title = element_text(hjust = 0),
  plot.caption = element_text(hjust = 0),
  panel.grid.minor.y = element_blank(),
  panel.grid.major.y = element_blank(),
  panel.grid.major.x = element_blank(),
  panel.grid.minor.x = element_blank(),
  strip.background.y = element_blank())

p2 <- p2 + facet_grid(cols = vars(facet), scales = "free") +
  theme(strip.background = element_rect(fill = 'grey89', color = 'white'))


# ggsave(p2, filename = 'DESEQ2UPSET_CONTRAST_C_D.png', path = wd, width = 5, height = 3.5, device = png, dpi = 300)

psave <- upplot + p2 + plot_layout(width = c(6, 3)) 

path_out <- "~/Documents/MIRNA_HALIOTIS/MIRNA_PUB_2024"

ggsave(psave, filename = 'DESEQ2PLOT_CONTRAST_C_D.png', path = path_out, width = 7, height = 3.2, device = png, dpi = 300)


# using density func


pdens <- RES_CC %>%
  filter(CONTRAST %in% WHICH_CONTRAST) %>%
  dplyr::mutate(SIGN = sign(logFC)) %>%
  dplyr::mutate(CONTRAST = dplyr::recode_factor(CONTRAST, !!!recode_to)) %>%
  ggplot(aes(x = log2(baseMean), y = logFC)) +
  facet_grid(~ CONTRAST, scales = "free") +
  geom_point(aes(alpha = -log10(padj)), color = 'grey') +
  # stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  # scale_fill_viridis_c() +
  geom_density_2d(aes(color = ..level..), linewidth = 0.5) +
  scale_color_viridis_c() +
  # stat_density_2d(
  #   geom = "raster",
  #   aes(fill = after_stat(density)),
  #   contour = FALSE
  # ) +
  geom_abline(slope = 0, intercept = 0, linetype="dashed", alpha=0.5) +
  theme_bw(base_family = "GillSans", base_size = 12) +
  theme(legend.position = "top",
    panel.border = element_blank(),
    strip.background = element_rect(fill = 'grey89', color = 'white'),
    plot.title = element_text(hjust = 0),
    plot.caption = element_text(hjust = 0),
    # panel.grid.major = element_blank(),
    panel.grid = element_blank(),
    strip.background.y = element_blank(),
    axis.text.x = element_text(angle = 0))



SUBSET_RES <- RES_CC %>%
  mutate(cc = ifelse(padj < 0.05, "p-value","")) %>%
  # left_join(ANNOT, by = "transcript_id") %>%
  filter(abs(logFC) > 5) %>%
  dplyr::mutate(CONTRAST = dplyr::recode_factor(CONTRAST, !!!recode_to))


pdens +
  # ylim(c(-10,10))
  geom_rect(
    aes(ymin=-30, ymax = -1, xmin = 1, xmax = Inf), fill = '#DADADA', alpha = 0.02) +
  geom_rect(
    aes(ymin=1, ymax = 30, xmin = 1, xmax = Inf), fill = '#DADADA', alpha = 0.02) +
  ggrepel::geom_text_repel(data = SUBSET_RES, aes(label = Name),
    size = 5, family = "GillSans", max.overlaps = 100)


library(ggdensity)

# https://jamesotto852.github.io/ggdensity/

RES_CC %>%
  filter(CONTRAST %in% WHICH_CONTRAST) %>%
  dplyr::mutate(SIGN = sign(logFC)) %>%
  dplyr::mutate(CONTRAST = dplyr::recode_factor(CONTRAST, !!!recode_to)) %>%
  ggplot(aes(x = log2(baseMean), y = logFC, color = CONTRAST)) +
  facet_grid(~ CONTRAST, scales = "free") +
  geom_point() +
  geom_hdr_lines() +
  scale_color_manual(values = c("black", "black")) +
  # geom_rug() +
  # scale_fill_viridis_c() +
  # geom_density_2d_filled(aes(colour = CONTRAST), linewidth = 0.5)
  # geom_density_2d(aes(colour = CONTRAST), linewidth = 0.5) +
  geom_abline(slope = 0, intercept = 0, linetype="dashed", alpha=0.5) +
  theme_bw(base_family = "GillSans", base_size = 12) +
  theme(legend.position = "top",
    panel.border = element_blank(),
    strip.background = element_rect(fill = 'grey89', color = 'white'),
    plot.title = element_text(hjust = 0),
    plot.caption = element_text(hjust = 0),
    # panel.grid.major = element_blank(),
    panel.grid = element_blank(),
    strip.background.y = element_blank(),
    axis.text.x = element_text(angle = 0)) +
  ggrepel::geom_text_repel(data = SUBSET_RES, aes(label = Name),
    size = 3, family = "GillSans", max.overlaps = 100)
