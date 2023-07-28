
# LOAD DATA RESULTS FROM DESEQ ANALYSIS
# FILTER DATA BY PADJ < 0.05 (NOT USE LOG2FC)
# GENERATE BARPLOT OF LOG2FC W/ SIGNIFICANCE FOR KNOWN AND NOVEL MIRS
# SUBSET EXCLUSIVE DE PER CONTRAST 
# VIZ UPSET

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

wd <- "~/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out/"

RES <- read_tsv(paste0(wd, "DESEQ_RES.tsv"))

# wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/FUNCTIONAL_MIR_ANNOT/"

# print(SRNA2GO <- read_tsv(paste0(wd, "SRNA_REGULATORY_FUNCTION_DB.tsv")))

url <- "https://raw.githubusercontent.com/RJEGR/Cancer_sete_T_assembly/main/functions.R"

source(url)

CONTRAST <- RES %>% distinct(CONTRAST) %>% pull()

recode_to <- c("24 HPF | Ctrl : Low pH", "110 HPF | Ctrl : Low pH", "Ctrl pH | 24 : 110 HPF", "Low pH | 24 : 110 HPF")

recode_to <- structure(recode_to, names = CONTRAST)

recode_fc <- structure(c("Up","Down"), names = c(1,-1))

RES.P <- RES %>% filter( padj < 0.05) # abs(log2FoldChange) > 2 &


# 1) =====


RES.P <- RES.P %>%
  mutate(star = ifelse(padj <.001, "***", 
  ifelse(padj <.01, "**",
    ifelse(padj <.05, "*", "")))) %>%
  mutate(SIGN = sign(log2FoldChange)) %>%
  dplyr::mutate(CONTRAST = dplyr::recode_factor(CONTRAST, !!!recode_to)) %>%
  dplyr::mutate(SIGN = dplyr::recode_factor(SIGN, !!!recode_fc)) %>%
  separate(CONTRAST, into = c("WRAP","CONTRAST"), sep = "[|]")


RES.P %>% 
  # drop_na(padj)
  mutate(YWRAP = ifelse(grepl("Cluster", Family), "B", "A")) %>%
  mutate(
    ymin = (abs(log2FoldChange) - lfcSE) * sign(log2FoldChange),
    ymax = (abs(log2FoldChange) + lfcSE) * sign(log2FoldChange),
    y_star = ymax + (0.15+lfcSE)* sign(log2FoldChange)) %>% 
  # mutate(Phylum = factor(Phylum,  levels = labels)) %>%
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

p <- p +  ggh4x::facet_nested( YWRAP ~ CONTRAST+WRAP, nest_line = F, scales = "free_y", space = "free_y") +
  theme(strip.background = element_rect(fill = 'grey89', color = 'white')) 
  
ggsave(p, filename = 'DESEQ2BARPLOT.png', path = wd, width = 5, height = 10, device = png, dpi = 300)



# upset ====

RES.P %>% 
  mutate(SIGN = sign(log2FoldChange)) %>%
  dplyr::count(CONTRAST, WRAP, SIGN, sampleB, sampleA) 
  # dplyr::mutate(CONTRAST = dplyr::recode_factor(CONTRAST, !!!recode_to)) %>%
  # dplyr::mutate(SIGN = dplyr::recode_factor(SIGN, !!!recode_fc)) %>%
  # separate(CONTRAST, into = c("WRAP","CONTRAST"), sep = "[|]")

library(ggupset)


RES %>% 
  filter(padj < 0.05) %>%
  mutate(SIGN = sign(log2FoldChange)) %>%
  dplyr::mutate(CONTRAST = dplyr::recode(CONTRAST, !!!recode_to)) %>%
  dplyr::mutate(SIGN = dplyr::recode_factor(SIGN, !!!recode_fc)) %>%
  # separate(CONTRAST, into = c("WRAP","CONTRAST"), sep = "[|]") %>%
  group_by(Name, SIGN) %>%
  summarise(across(CONTRAST, .fns = list), n = n()) %>% 
  ggplot(aes(x = CONTRAST, fill = SIGN)) +
  geom_bar(position = position_dodge(width = 1)) +
  geom_text(stat='count', aes(label = after_stat(count)), 
    position = position_dodge(width = 1), vjust = -0.5, family = "GillSans", size = 2.5) +
  scale_x_upset(order_by = "freq", reverse = F) +
  theme_bw(base_family = "GillSans") +
  theme_combmatrix(combmatrix.panel.point.color.fill = "black",
    combmatrix.panel.line.size = 0, base_family = "GillSans") +
  labs(x = '', y = 'Number of Intersected miRs') +
  ggsci::scale_fill_aaas() +
  ggsci::scale_color_aaas() +
  guides(fill = guide_legend(title = "", nrow = 1)) -> p

p <- p + theme(legend.position = "top",
  panel.border = element_blank(),
  plot.title = element_text(hjust = 0),
  plot.caption = element_text(hjust = 0),
  panel.grid.minor.y = element_blank(),
  panel.grid.major.y = element_blank(),
   panel.grid.major.x = element_blank(),
  panel.grid.minor.x = element_blank(),
  strip.background.y = element_blank())

ggsave(p, filename = 'DESEQ2UPSET.png', path = wd, width = 5, height = 5.2, device = png, dpi = 300)


# Generate single / exclusive transcripts per condition ----
  
RES %>% 
  filter(abs(log2FoldChange) > 2) %>%
  mutate(SIGN = sign(log2FoldChange)) %>%
  group_by(Name, SIGN) %>%
  summarise(n_contrast = length(CONTRAST), 
    across(CONTRAST, .fns = list)) -> distinct_trans

distinct_trans %>% group_by(n_contrast) %>% tally()

