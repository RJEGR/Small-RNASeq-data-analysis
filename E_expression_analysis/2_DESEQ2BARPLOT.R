
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

wd <- "~/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out/"

RES <- read_tsv(paste0(wd, "DESEQ_RES.tsv")) # %>% mutate(log2FoldChange = log2FoldChange * -1)

url <- "https://raw.githubusercontent.com/RJEGR/Cancer_sete_T_assembly/main/functions.R"

source(url)

CONTRAST <- RES %>% distinct(CONTRAST) %>% pull()

recode_to <- c("24 HPF | Ctrl : Low pH", "110 HPF | Ctrl : Low pH", "Ctrl pH | 24 : 110 HPF", "Low pH | 24 : 110 HPF")

recode_to <- structure(recode_to, names = CONTRAST)

recode_fc <- structure(c("Up","Down"), names = c(1,-1))

RES.P <- RES %>% filter( padj < 0.05 & abs(log2FoldChange) > 1) # abs(log2FoldChange) > 2 &

RES.P <- RES.P %>%
  mutate(star = ifelse(padj <.001, "***", 
    ifelse(padj <.01, "**",
      ifelse(padj <.05, "*", "")))) %>%
  mutate(SIGN = sign(log2FoldChange)) %>%
  mutate(CONTRAST_DE = CONTRAST) %>%
  dplyr::mutate(CONTRAST = dplyr::recode_factor(CONTRAST, !!!recode_to)) %>%
  dplyr::mutate(SIGN = dplyr::recode_factor(SIGN, !!!recode_fc)) %>%
  separate(CONTRAST, into = c("WRAP","CONTRAST"), sep = "[|]")

write_tsv(RES.P, paste0(wd, "DESEQ_RES_P.tsv"))

# 1) =====


RES.P %>% 
  # drop_na(padj)
  mutate(YWRAP = ifelse(grepl("Cluster", Family), "B", "A")) %>%
  mutate(
    ymin = (abs(log2FoldChange) - lfcSE) * sign(log2FoldChange),
    ymax = (abs(log2FoldChange) + lfcSE) * sign(log2FoldChange),
    y_star = ymax + (0.15+lfcSE)* sign(log2FoldChange)) %>% 
  ggplot(aes(x = Family, y = log2FoldChange, fill = SIGN)) + 
  geom_bar(stat = "identity", width = 0.5, 
    position = position_identity()) +
  coord_flip() +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.25,
    position = position_identity(), color = "black") + 
  geom_text(aes(y = y_star, label=star), 
    vjust=  .7, color="black", position = position_identity(), family = "GillSans", size = 1.5) +
  ggsci::scale_fill_aaas() +
  guides(color = "none") +
  labs(x = NULL, y = "Log fold change") +
  guides(fill = guide_legend(title = "", nrow = 1)) +
  theme_bw(base_family = "GillSans", base_size = 12) +
  # theme_classic(base_size = 16, base_family = "GillSans") + 
  # geom_abline(slope = 0, intercept = 0, linetype = "dashed", alpha = 0.5) +
  theme(legend.position = "top",
    panel.border = element_blank(),
    plot.title = element_text(hjust = 0),
    plot.caption = element_text(hjust = 0),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    strip.background.y = element_blank(),
    # axis.text.y.right = element_text(angle = 0, hjust = 1, vjust = 0, size = 2.5),
    axis.text.y = element_text(angle = 0, size = 5),
    axis.text.x = element_text(angle = 0)) -> p

p <- p +  ggh4x::facet_nested( YWRAP ~ CONTRAST_DE, nest_line = F, scales = "free_y", space = "free_y") +
  theme(strip.background = element_rect(fill = 'grey89', color = 'white')) 
  
ggsave(p, filename = 'DESEQ2BARPLOT.png', path = wd, width = 5, height = 10, device = png, dpi = 300)

# ONLY CONTRAST A AND B ====

str(SORT_MIRS <- c("MIR-278","MIR-133","LET-7",
  "Cluster_55760","Cluster_55776","MIR-2","MIR-315", "MIR-153", "BANTAM", "MIR-190", "MIR-2722", "MIR-1988", 
  "MIR-92", "MIR-277B", "MIR-216"))


recode_c <- structure(c("A) 24 HPF","B) 110 HPF"), names = c("CONTRAST_A", "CONTRAST_B"))
recode_fc <- structure(c("pH 8.0","pH 7.6"), names = c("Up","Down"))


RES.P %>%
  filter(CONTRAST_DE %in% c("CONTRAST_A", "CONTRAST_B")) %>%
  dplyr::mutate(CONTRAST_DE = dplyr::recode_factor(CONTRAST_DE, !!!recode_c)) %>%
  dplyr::mutate(SIGN = dplyr::recode_factor(SIGN, !!!recode_fc)) %>% 
  arrange(desc(log2FoldChange)) %>%
  mutate(Family = factor(Family, levels = unique(Family))) %>%
  mutate(
    ymin = (abs(log2FoldChange) - lfcSE) * sign(log2FoldChange),
    ymax = (abs(log2FoldChange) + lfcSE) * sign(log2FoldChange),
    y_star = ymax + (0.15+lfcSE)* sign(log2FoldChange)) %>% 
  mutate(Family = factor(Family, levels = srna_order)) %>%
  ggplot(aes(x = Family, y = log2FoldChange, fill = SIGN)) + 
  geom_bar(stat = "identity", width = 0.5, 
    position = position_identity()) +
  coord_flip() +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.25,
    position = position_identity(), color = "black") + 
  geom_text(aes(y = y_star, label=star), 
    vjust=  .7, color="black", position = position_identity(), family = "GillSans", size = 1.5) +
  ggsci::scale_fill_aaas() +
  guides(color = "top") +
  labs(x = NULL, y = "Log fold change") +
  guides(fill = guide_legend(title = "", nrow = 1)) +
  theme_bw(base_family = "GillSans", base_size = 12) +
  theme(legend.position = "top",
    strip.background = element_rect(fill = 'grey89', color = 'white'),
    panel.border = element_blank(),
    plot.title = element_text(hjust = 0),
    plot.caption = element_text(hjust = 0),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.y = element_text(angle = 0, size = 6),
    axis.text.x = element_text(angle = 0)) +  
  ggh4x::facet_nested(~ CONTRAST_DE , nest_line = F, scales = "free_y", space = "free_y") -> p


ggsave(p, filename = 'DESEQ2BARPLOT_ONLY_DEGS.png', path = wd, width = 4, height = 5, device = png, dpi = 300)

# SEE DESEQ2TOPGO_UPGRADE TO BIND HEATMAP OF MULTIPLE TARGET PROCESS.

# 2 ====
# not only unique

RES.P %>% 
  mutate(SIGN = sign(log2FoldChange)) %>%
  dplyr::count(CONTRAST, WRAP, SIGN, sampleB) %>%
  dplyr::mutate(CONTRAST = dplyr::recode_factor(CONTRAST, !!!recode_to)) %>%
  dplyr::mutate(SIGN = dplyr::recode_factor(SIGN, !!!recode_fc)) 
  # separate(CONTRAST, into = c("WRAP","CONTRAST"), sep = "[|]")


library(ggupset)

# UPSETDF <- RES %>% 
#   filter(padj < 0.05) %>%
#   mutate(SIGN = sign(log2FoldChange)) %>%
#   dplyr::mutate(CONTRAST = dplyr::recode(CONTRAST, !!!recode_to)) %>%
#   dplyr::mutate(SIGN = dplyr::recode_factor(SIGN, !!!recode_fc)) %>%
#   # separate(CONTRAST, into = c("CONTRAST","WRAP"), sep = "[|]") %>%
#   group_by(Name, SIGN) %>%
#   summarise(across(CONTRAST, .fns = list), n = n()) 

# 1) omit

recode_fc <- structure(c("pH 7.6","pH 8.0"), names = c(-1,1))

recode_to <-  structure(c("24 HPF", "110 HPF"), names = c("CONTRAST_A", "CONTRAST_B"))

UPSETDF <- RES.P %>% 
  filter(CONTRAST_DE %in% c("CONTRAST_A", "CONTRAST_B")) %>%
  mutate(SIGN = sign(log2FoldChange)) %>%
  dplyr::mutate(SIGN = dplyr::recode_factor(SIGN, !!!recode_fc)) %>%
  dplyr::mutate(CONTRAST_DE = dplyr::recode_factor(CONTRAST_DE, !!!recode_to)) %>%
  group_by(Name, SIGN) %>%
  summarise(across(CONTRAST_DE, .fns = list), n = n()) %>%
  mutate(SIGN = factor(SIGN, levels = rev(recode_fc)))


UPSETDF <- UPSETDF %>%
  ungroup() %>%
  left_join(distinct(RES.P, Name, Family)) %>%
  distinct(Family, .keep_all = T)

# UPSETDF <- RES.P %>%
#   # unite("CONTRAST", CONTRAST:WRAP, sep = " | ") %>%
#   group_by(Name, SIGN) %>%
#   summarise(across(CONTRAST_DE, .fns = list), n = n())

# Levels <- RES.P %>% unite("CONTRAST", CONTRAST:WRAP, sep = " | ") %>% distinct(CONTRAST) %>% pull()

col <- c("#3B4992FF","#EE0000FF")  

UPSETDF %>%
  ggplot(aes(x = CONTRAST_DE, fill = SIGN)) +
  geom_bar(position = position_dodge(width = 1)) +
  geom_text(stat='count', aes(label = after_stat(count)), 
    position = position_dodge(width = 1), vjust = -0.5, family = "GillSans", size = 3.5) +
  scale_x_upset(order_by = "degree", reverse = F) +
  theme_bw(base_family = "GillSans") +
  theme_combmatrix(combmatrix.panel.point.color.fill = "black",
    combmatrix.panel.line.size = 0, base_family = "GillSans") +
  axis_combmatrix(levels = c("24 HPF", "110 HPF")) +
  labs(x = '', y = 'Number of miRs') +
  scale_color_manual("", values = col) +
  scale_fill_manual("", values =  col) +
  guides(fill = guide_legend(title = "", nrow = 1)) -> p1

p1 <- p1 + theme(legend.position = "top",
  panel.border = element_blank(),
  plot.title = element_text(hjust = 0),
  plot.caption = element_text(hjust = 0),
  panel.grid.minor.y = element_blank(),
  panel.grid.major.y = element_blank(),
   panel.grid.major.x = element_blank(),
  panel.grid.minor.x = element_blank(),
  strip.background.y = element_blank())

# needs manual eddition if facetwrawp
# p1 + facet_grid(cols = vars(SIGN), scales = "free") + theme(strip.background = element_rect(fill = 'grey89', color = 'white'))

ggsave(p1, filename = 'DESEQ2UPSET_CONTRAST_A_B.png', path = wd, width = 5, height = 5, device = png, dpi = 300)


# 2) BY DEVELOPMENT ====

# INCLUDE VOLCANO ====


WHICH_CONTRAST <- c("CONTRAST_C", "CONTRAST_D")

recode_to <-  structure(c("A) pH 8.0", "B) pH 7.6"), names = WHICH_CONTRAST)

RES_CC <- prep_DE_data(RES, alpha = 0.05, lfcThreshold = 2)

colors_fc <- c("forestgreen",  "#4169E1", "red2")

c("#DADADA", "#D4DBC2")

upplot <- RES_CC %>% 
  filter(CONTRAST %in% WHICH_CONTRAST) %>%
  dplyr::mutate(SIGN = sign(logFC)) %>%
  dplyr::mutate(CONTRAST = dplyr::recode_factor(CONTRAST, !!!recode_to)) %>%
  filter(cc != "NS") %>%
  mutate(logFC = logFC*-1) %>% # reverse sort of HPF
  ggplot(aes(x = logFC, y = -log10(padj))) +
  facet_grid(~ CONTRAST) +
  geom_rect(
    aes(xmin=-15, xmax = -1, ymin = 1, ymax = Inf), fill = '#DADADA') +
  geom_rect(
    aes(xmin=1, xmax = 30, ymin = 1, ymax = Inf), fill = '#D4DBC2') +
  geom_point(aes(color = cc), alpha = 3/5) +
  scale_color_manual(name = "", values = colors_fc) +
  labs(x= expression(Log[2] ~ "Fold Change"), 
    y = expression(-Log[10] ~ "padj")) +
  theme_bw(base_family = "GillSans", base_size = 12) +
  geom_abline(slope = 0, intercept = -log10(0.05), linetype="dashed", alpha=0.5) +
  geom_vline(xintercept = 1, linetype="dashed", alpha=0.5) +
  geom_vline(xintercept = -1, linetype="dashed", alpha=0.5) +
  annotate("text", x = -10, y = 150, label = "24 hpf", family = "GillSans") +
  annotate("text", x = 10, y = 150, label = "110 hpf", family = "GillSans") +
  xlim(-15, 30)

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

UPSETDF <- UPSETDF %>%
  ungroup() %>%
  left_join(distinct(RES.P, Name, Family)) %>%
  distinct(Family, .keep_all = T)

library(ggupset)

UPSETDF %>%
  mutate(facet = "C)") %>%
  ggplot(aes(x = CONTRAST_DE, fill = SIGN)) +
  geom_bar(position = position_dodge(width = 1)) +
  geom_text(stat='count', aes(label = after_stat(count)), 
    position = position_dodge(width = 1), vjust = -0.2, family = "GillSans", size = 3.5) +
  ggupset::scale_x_upset(order_by = "degree", reverse = F) +
  theme_bw(base_family = "GillSans", base_size = 12) +
  ggupset::theme_combmatrix(combmatrix.panel.point.color.fill = "black",
    combmatrix.panel.line.size = 0, base_family = "GillSans") +
  axis_combmatrix(levels = c("pH 8.0", "pH 7.6"), clip = "off") +
  labs(x = '', y = 'Número de miRs') +
  scale_color_manual("", values = c("#DADADA", "#D4DBC2")) +
  scale_fill_manual("", values =  c("#DADADA", "#D4DBC2")) +
  guides(fill = guide_legend(title = "", nrow = 1)) +
  ylim(c(0,50)) -> p2

p2 <- p2 + theme(legend.position = "top",
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


ggsave(psave, filename = 'DESEQ2PLOT_CONTRAST_C_D.png', path = wd, width = 7, height = 3.2, device = png, dpi = 300)


# 3) SPLIT EXCLUSIVE FROM INTERSECTED MIRS ----
# INCLUDE ALL CONTRAST

UPSETDF <- RES.P %>% 
  mutate(SIGN = sign(log2FoldChange)) %>%
  # dplyr::mutate(SIGN = dplyr::recode_factor(SIGN, !!!recode_fc)) %>%
  # dplyr::mutate(CONTRAST_DE = dplyr::recode_factor(CONTRAST_DE, !!!recode_to)) %>%
  group_by(Name, SIGN) %>%
  summarise(across(CONTRAST, .fns = list), n = n())

UPSETDF %>% group_by(n) %>% tally()

UPSETDF %>% group_by(n, SIGN) %>% tally()

# Sanity check if unique:

str(query.ids <- UPSETDF %>% filter(n == 1) %>% pull(Name))

# str(query.ids <- EXCLUSIVE_MIRS %>% count(Name, sort = T) %>% filter(n == 2) %>% pull(Name))

RES.P %>% count(sampleA, sampleB)

RES.P %>%
  unite("CONTRAST", CONTRAST:WRAP, sep = " | ") %>%
  filter(Name %in% query.ids) %>% # view()
  # dplyr::select(Name, baseMeanA, baseMeanB, CONTRAST, WRAP, SIGN) %>%
  pivot_longer(cols = c("baseMeanA", "baseMeanB")) %>%
  ggplot(aes(x = value, y = Name, fill = name)) +
  facet_grid(SIGN ~ CONTRAST) +
  geom_col(position = position_dodge2())


EXCLUSIVE_MIRS <- UPSETDF %>% filter(n == 1) %>%
  mutate(CONTRAST = unlist(CONTRAST)) %>%
  left_join(RES.P %>% dplyr::distinct(Name, Family))
  # dplyr::select(-n)

EXCLUSIVE_MIRS %>% ungroup() %>% count(SIGN, CONTRAST)  %>% view()

INTERSECTED_MIRS <- UPSETDF %>% filter(n != 1) %>%
  unnest(CONTRAST) %>%
  # dplyr::rename("CONTRAST" = "CONTRAST_DE") %>%
  left_join(RES.P %>% dplyr::distinct(Name, Family))

out <- list(EXCLUSIVE_MIRS, INTERSECTED_MIRS)

write_rds(out, file = paste0(wd, "UPSETDF.rds"))


# HEATMAP OF DEGS (z-score) ====

dds <- read_rds(paste0(wd, "/DDS_DESEQ2.rds"))

which_mirs <- RES.P %>% 
  # filter(CONTRAST_DE %in% WHICH_CONTRAST) %>%
  distinct(Name) %>% pull()

str(which_mirs)  
  
# divide the counts by the size factors or normalization factors before
barplot(cnts <- DESeq2::counts(dds, normalized = T, replaced = F)[which_mirs,])

barplot(cnts <- DESeq2::varianceStabilizingTransformation(round(cnts)))

# as z score

cnts <- t(scale(t(cnts)))

heatmap(cnts)

#

RES %>% 
  # filter(CONTRAST %in% WHICH_CONTRAST) %>%
  dplyr::mutate(CONTRAST = dplyr::recode_factor(CONTRAST, !!!recode_to)) %>%
  # filter(cc != "NS") %>%
  mutate(log2FoldChange = log2FoldChange*-1) %>% # reverse sort of HPF
  ggplot(aes(x = log10(baseMean), y = log2FoldChange)) +
  facet_grid(~ CONTRAST) +
  # geom_rect(
  #   aes(ymin=-15, ymax = -1, xmin = 1, xmax = Inf), fill = '#DADADA') +
  # geom_rect(
  #   aes(ymin=1, ymax = 30, xmin = 1, xmax = Inf), fill = '#D4DBC2') +
  geom_point(aes(alpha = -log10(padj)), color = 'black') +
  # scale_color_manual(name = "", values = colors_fc) +
  labs(y= expression(Log[10] ~ "BaseMean")) +
  theme_bw(base_family = "GillSans", base_size = 12) + 
  geom_abline(slope = 0, intercept = -2, linetype="dashed", alpha=0.5) +
  geom_abline(slope = 0, intercept = 2, linetype="dashed", alpha=0.5) +
  geom_density_2d(aes(color = ..level..), linewidth = 0.5) +
  scale_color_viridis_c() +
  theme(legend.position = "top",
    panel.border = element_blank(),
    strip.background = element_rect(fill = 'grey89', color = 'white'),
    plot.title = element_text(hjust = 0),
    plot.caption = element_text(hjust = 0),
    # panel.grid.major = element_blank(),
    panel.grid = element_blank(),
    strip.background.y = element_blank(),
    axis.text.x = element_text(angle = 0)) +
  # geom_vline(xintercept = 1, linetype="dashed", alpha=0.5) +
  # geom_vline(xintercept = -1, linetype="dashed", alpha=0.5) +
  annotate("text", x = -10, y = 150, label = "24 hpf", family = "GillSans") +
  annotate("text", x = 10, y = 150, label = "110 hpf", family = "GillSans") +
  xlim(-15, 30)

RES %>% 
  filter(CONTRAST %in% WHICH_CONTRAST) %>%
  dplyr::mutate(CONTRAST = dplyr::recode_factor(CONTRAST, !!!recode_to)) %>%
  mutate(log2FoldChange = log2FoldChange*-1) %>% # reverse sort of HPF
  ggplot(aes(x = log10(baseMean), y = log2FoldChange)) +
  geom_point(aes(alpha = -log10(padj)), color = 'black') +
  # geom_density_2d(aes(color = CONTRAST), linewidth = 0.5) +
  theme_bw(base_family = "GillSans", base_size = 12) + 
  facet_grid(~ CONTRAST) +
  stat_density_2d(
    geom = "raster",
    aes(fill = after_stat(density)),
    contour = FALSE
  ) +
  geom_abline(slope = 0, intercept = -2, linetype="dashed", alpha=0.5) +
  geom_abline(slope = 0, intercept = 2, linetype="dashed", alpha=0.5) +
  theme(legend.position = "top",
    panel.border = element_blank(),
    strip.background = element_rect(fill = 'grey89', color = 'white'),
    plot.title = element_text(hjust = 0),
    plot.caption = element_text(hjust = 0),
    # panel.grid.major = element_blank(),
    panel.grid = element_blank(),
    strip.background.y = element_blank(),
    axis.text.x = element_text(angle = 0)) 
  
  
  
