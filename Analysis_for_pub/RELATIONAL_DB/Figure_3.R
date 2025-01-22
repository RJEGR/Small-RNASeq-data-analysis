
library(tidyverse)

rm(list = ls());

if(!is.null(dev.list())) dev.off()

dir <- "~/Documents/MIRNA_HALIOTIS/MIRNA_PUB_2024/"

print(DB <- read_tsv(file.path(dir, "LONGER_RELATIONAL_DB.tsv")))

RES <- read_tsv(list.files(path = dir, 
  pattern = "SEQUENCES_MERGED_DESEQ_RES.tsv", full.names = T)) %>%
  mutate(Contrast = ifelse(sign(log2FoldChange) == -1, sampleB, sampleA)) %>%
  mutate(Contrast = paste0(CONTRAST,"_",Contrast))

# RES.P <- RES %>% filter( padj < 0.05  & abs(log2FoldChange) > 1)

RES.p <- RES %>% 
  filter(CONTRAST %in% c("CONTRAST_A","CONTRAST_B")) %>%
  # mutate(CONTRAST = dplyr::recode_factor(CONTRAST, !!!recode_to)) %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1)

QUERIES <- RES.p %>% 
  distinct(MajorRNA, MirGeneDB_ID) %>%
  pull(MajorRNA, name = MirGeneDB_ID)


# QUERIES

DF1 <- DB %>% 
  # mutate(STRINGID = ifelse(is.na(STRINGID), gene_id, STRINGID)) %>% 
  distinct(STRINGID, gene_id, MajorRNA, MajorRNAID) %>%
  count(STRINGID, gene_id, sort = T) %>%
  dplyr::rename("MIR_degree" = "n")



paste_col <- function(x) { 
  x <- x[!is.na(x)] 
  x <- unique(sort(x))
  x <- paste(x, sep = ';', collapse = ';')
}

DF2 <- DB %>% 
  filter(MajorRNA %in% QUERIES) %>%
  distinct(gene_id, MajorRNAID) %>%
  mutate(Contrast = ifelse(MajorRNAID %in% names(QUERIES[1:3]), "24 hpf", "110 hpf")) %>%
  group_by(gene_id, Contrast) %>%
  summarise(across(MajorRNAID, .fns = paste_col), OA_MIR_degree = n()) %>%
  arrange(desc(OA_MIR_degree)) %>%
  mutate(OA_MIR_degree) %>%
  mutate(OA_MIR_degree = ifelse(is.na(OA_MIR_degree), 0, OA_MIR_degree))
# dplyr::rename("preferred_name" = "STRINGID")


NodeDF <- DF2 %>% right_join(DF1 , by = "gene_id") %>% 
  mutate(STRINGID = ifelse(is.na(STRINGID), gene_id, STRINGID)) %>%
  mutate(OA_MIR_degree = ifelse(is.na(OA_MIR_degree), 0, OA_MIR_degree)) %>% ungroup()

# add COG_cat or biological theme

print(BT <- read_tsv(file.path(dir, "Supplementary_tables - Biological_themes.tsv"), col_names = T))

BT <- BT %>% select(-microRNA_degree, -STRINGID)

DataViz <- DB %>% 
  left_join(BT, by = "gene_id") %>%
  mutate(y_axis = Bilogical_pathway) %>% # or COG_category
  # mutate(y_axis = ifelse(is.na(y_axis), COG_category, y_axis)) %>%
  # drop_na(COG_category) %>%
  distinct(gene_id, y_axis) %>%
  # drop_na(y_axis) %>%
  # count(gene_id, COG_category, sort = T) %>% filter(n >1 )
  right_join(NodeDF, by = "gene_id")

# DataViz <- DataViz %>% mutate(y_axis = ifelse(is.na(y_axis), "Unknown", y_axis)) 

DataViz <- DataViz %>% drop_na(y_axis)

DataViz <- DataViz %>%
  mutate(Frac = OA_MIR_degree) %>%
  filter(Frac > 0) %>%
  group_by(Contrast, y_axis) %>% 
  summarise(OA_MIR_degree = sum(OA_MIR_degree), MIR_degree = sum(MIR_degree)) %>%
  group_by(Contrast) %>%
  arrange(desc(OA_MIR_degree), .by_group = T) %>% 
  mutate(Label = as.character(y_axis), row_number = row_number()) %>% 
  mutate(Label = paste(Label, Contrast, sep = "__")) %>%
  mutate(Label = factor(Label, levels = rev(unique(Label)))) %>%
  # mutate(Label = factor(paste(Label, row_number, sep = "__"), levels = rev(paste(Label, row_number, sep = "__")))) %>%
  mutate(Contrast = factor(Contrast, levels = c("24 hpf", "110 hpf")))


DataViz %>%
  ggplot(aes(y = Label, x = OA_MIR_degree)) +
  # facet_grid( COG_name ~ ., space = "free_y", scales = "free_y") +
  ggforce::facet_col(~ Contrast, space = "free", scales = "free_y") +
  # facet_wrap(~ COG_name, ncol = 1, scales = "free_y") +
  geom_point(shape = 21, size = 1, stroke = 1.2, fill = "white") +
  geom_segment(aes(xend = 0, yend = Label), linewidth = 1.3) +
  scale_y_discrete("",labels = function(x) gsub("__.+$", "", x)) +
  labs(x = "MicroRNA degree") +
  theme_bw(base_family = "GillSans", base_size = 10) +
  theme(
    strip.background = element_rect(fill = 'white', color = 'white'),
    strip.text = element_text(color = "black",hjust = 0),
    legend.position = 'top',
    # panel.border = element_blank(),
    plot.margin = unit(c(0,0,0,0), "pt"),
    panel.grid.minor = element_blank(),
    # axis.ticks.y = element_blank(),
    # axis.text.y = element_blank(),
    # axis.title = element_blank(),
    panel.grid.major = element_blank()) -> P1


# 2 -----

RESviz <- RES.p %>% distinct(MajorRNA, MirGeneDB_ID, log2FoldChange) %>%
  # mutate(Contrast = ifelse(MirGeneDB_ID == "MIR-133-3p" & log2FoldChange > 0, NA, "")) %>% drop_na(Contrast) %>%
  filter(log2FoldChange < 0) %>%
  mutate(Contrast = ifelse(MirGeneDB_ID %in% names(QUERIES[1:3]), "24 hpf", "110 hpf")) %>%
  mutate(Contrast = factor(Contrast, levels = c("24 hpf", "110 hpf")))



DB %>% 
  filter(MajorRNA %in% QUERIES) %>%
  left_join(BT, by = "gene_id") %>%
  mutate(y_axis = Bilogical_pathway) %>% # or COG_category
  # mutate(y_axis = ifelse(is.na(y_axis), COG_category, y_axis)) %>%
  # drop_na(COG_category) %>%
  distinct(MajorRNA, y_axis) %>%
  right_join(RESviz, by = "MajorRNA") %>%
  drop_na(y_axis) %>%
  mutate(Label = paste(y_axis, Contrast, sep = "__")) %>%
  mutate(Label = factor(Label, levels = levels(DataViz$Label))) %>%
  ggplot(aes(y = Label, x = log2FoldChange*-1)) +
  scale_y_discrete("", labels = function(x) gsub("__.+$", "", x)) +
  # geom_jitter(alpha = 0.5) +
  stat_boxplot(geom ='errorbar', width = 0.3, linetype="dashed") +
  geom_point(shape = 21, size = 1, stroke = 1.2, fill = "white", alpha = 0.5) +
  # geom_segment(aes(xend = 0, yend = COG_name), linewidth = 0.5, color="#E7DFD5") +
  # geom_boxplot(width = 0.3, outlier.alpha = 0, linetype="dashed") +
  ggforce::facet_col(~ Contrast, space = "free", scales = "free_y") +
  theme_bw(base_family = "GillSans", base_size = 10) +
  # scale_color_manual(values = c("black", "black")) +
  labs(x = "Log2FoldChange") +
  theme(
    strip.background = element_rect(fill = 'white', color = 'white'),
    strip.text = element_text(color = "black",hjust = 0),
    legend.position = 'top',
    # panel.border = element_blank(),
    plot.margin = unit(c(0,0,0,0), "pt"),
    panel.grid.minor = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    # axis.title = element_blank(),
    panel.grid.major = element_blank()) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.5) -> P2




library(patchwork)

ps <- P1 + plot_spacer() + P2 + plot_layout(widths = c(5, -0.1, 4))

# ps

ggsave(ps, filename = 'Degree_and_log2fc.png', path = dir, width = 6, height = 4, device = png, dpi = 300)

