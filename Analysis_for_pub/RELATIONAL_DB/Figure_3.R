
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


QUERIES_PHCONTROL <- RES %>% 
  filter(CONTRAST %in% c("CONTRAST_C")) %>%
  # mutate(CONTRAST = dplyr::recode_factor(CONTRAST, !!!recode_to)) %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>% 
  distinct(MajorRNA, MirGeneDB_ID) %>%
  pull(MajorRNA, name = MirGeneDB_ID)


RES.p <- RES %>% 
  filter(CONTRAST %in% c("CONTRAST_A","CONTRAST_B")) %>%
  # mutate(CONTRAST = dplyr::recode_factor(CONTRAST, !!!recode_to)) %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1)

QUERIES <- RES.p %>% 
  distinct(MajorRNA, MirGeneDB_ID) %>%
  pull(MajorRNA, name = MirGeneDB_ID)


# QUERIES

DF1 <- DB %>% 
  filter(MajorRNA %in% QUERIES_PHCONTROL) %>%
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

BT <- BT %>% select(-STRINGID)

# rename ?
BT %>% distinct(Bilogical_pathway)

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

# 
# DataViz %>%
#   mutate(OA_MIR_degree = OA_MIR_degree/sum(OA_MIR_degree)) %>%
#   mutate(MIR_degree = MIR_degree/sum(MIR_degree)) %>%
#   ggplot(aes(y = Label)) +
#   # facet_grid( COG_name ~ ., space = "free_y", scales = "free_y") +
#   ggforce::facet_col(~ Contrast, space = "free", scales = "free_y") +
#   # facet_wrap(~ COG_name, ncol = 1, scales = "free_y") +
#   geom_point(aes(x = OA_MIR_degree), shape = 21, size = 1, stroke = 1.2, fill = "white", color = "black") +
#   geom_segment(aes(x = OA_MIR_degree, xend = 0, yend = Label), linewidth = 1.3, color = "black") +
#   scale_y_discrete("",labels = function(x) gsub("__.+$", "", x)) +
#   labs(x = "MicroRNA degree") +
#   theme_bw(base_family = "GillSans", base_size = 10) +
#   theme(
#     strip.background = element_rect(fill = 'white', color = 'white'),
#     strip.text = element_text(color = "black",hjust = 0),
#     legend.position = 'top',
#     # panel.border = element_blank(),
#     plot.margin = unit(c(0,0,0,0), "pt"),
#     panel.grid.minor = element_blank(),
#     # axis.ticks.y = element_blank(),
#     # axis.text.y = element_blank(),
#     # axis.title = element_blank(),
#     panel.grid.major = element_blank()) -> P1

# Including degree of development

DataViz %>%
  mutate(OA_MIR_degree = OA_MIR_degree/sum(OA_MIR_degree) ) %>%
  mutate(MIR_degree = MIR_degree/sum(MIR_degree) * -1) %>%
  pivot_longer(cols = c("MIR_degree","OA_MIR_degree")) %>%
  mutate(name = ifelse(name %in% "MIR_degree", "pH 8.0", "pH 7.6")) %>%
  ggplot(aes(y = Label, color = name, x = value)) +
  ggforce::facet_col(~ Contrast, space = "free", scales = "free_y") +
  # geom_point(shape = 21, size = 1, stroke = 1.2, fill = "white", position = position_dodge(width = 0.5)) +
  # replace segment to linerange to position dodge
  # geom_linerange(aes(xmin = 0, xmax =  value), linewidth = 1.3, position = position_dodge(width = 1)) +
  geom_segment(aes(x = value, xend = 0, yend = Label), linewidth = 1.3) +
  scale_y_discrete("",labels = function(x) gsub("__.+$", "", x)) +
  labs(x = "MicroRNA degree (Fraction)") +
  theme_bw(base_family = "GillSans", base_size = 10) +
  scale_color_manual(values = c("black", "gray86")) +
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

P1 <- P1 + 
  # theme(axis.text.y = element_blank(), axis.ticks.x = element_blank()) +
  guides(color=guide_legend(title = "", nrow = 1), fill = "none") +
  scale_x_continuous("MicroRNA density",labels = scales::percent)

# of to stack 100% bars per bio.theme?

DataViz %>%
  # mutate(OA_MIR_degree = OA_MIR_degree/sum(OA_MIR_degree)) %>%
  # mutate(MIR_degree = MIR_degree/sum(MIR_degree)) %>%
  mutate(OA_MIR_degree = OA_MIR_degree/(OA_MIR_degree+MIR_degree)) %>%
  mutate(MIR_degree = MIR_degree/(MIR_degree-OA_MIR_degree)) %>%
  pivot_longer(cols = c("MIR_degree","OA_MIR_degree")) %>%
  mutate(name = ifelse(name %in% "MIR_degree", "pH 8.0", "pH 7.6")) %>%
  ggplot(aes(y = Label, color = name, x = value)) +
  ggforce::facet_col(~ Contrast, space = "free", scales = "free_y") +
  # geom_point(shape = 21, size = 1, stroke = 1.2, fill = "white", position = position_dodge(width = 0.5)) +
  geom_col()



# 2 -----

RESviz <- RES.p %>% distinct(MajorRNA, MirGeneDB_ID, log2FoldChange) %>%
  filter(abs(log2FoldChange) > 0) %>% # abs(log2FoldChange)
  mutate(Contrast = ifelse(MirGeneDB_ID %in% names(QUERIES[1:3]), "24 hpf", "110 hpf")) %>%
  mutate(Contrast = factor(Contrast, levels = c("24 hpf", "110 hpf")))



DB %>% 
  filter(MajorRNA %in% QUERIES) %>%
  left_join(BT, by = "gene_id") %>%
  mutate(y_axis = Bilogical_pathway) %>% # or COG_category
  # mutate(y_axis = ifelse(is.na(y_axis), COG_category, y_axis)) %>%
  # drop_na(COG_category) %>%
  distinct(MajorRNA, y_axis) %>%
  right_join(RESviz %>% filter(log2FoldChange < 0), by = "MajorRNA") %>%
  drop_na(y_axis) %>%
  mutate(Label = paste(y_axis, Contrast, sep = "__")) %>%
  mutate(Label = factor(Label, levels = levels(DataViz$Label))) %>%
  ggplot(aes(y = Label, x = log2FoldChange*-1)) +
  scale_y_discrete("", labels = function(x) gsub("__.+$", "", x)) +
  stat_boxplot(geom ='errorbar', width = 0.3, linetype="dashed") +
  geom_point(shape = 21, size = 1, stroke = 1.2, fill = "white", alpha = 0.7, color = "black") +
  # scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
  ggforce::facet_col(~ Contrast, space = "free", scales = "free_y") +
  theme_bw(base_family = "GillSans", base_size = 10) +
  # scale_color_manual(values = c("black", "black")) +
  labs(x = "Log2FoldChange") +
  # geom_vline(xintercept = 0, linetype = "dashed", size = 0.5)
  theme(
    strip.background = element_rect(fill = 'white', color = 'white'),
    strip.text = element_text(color = "white",hjust = 0),
    legend.position = 'top',
    # panel.border = element_blank(),
    plot.margin = unit(c(0,0,0,0), "pt"),
    panel.grid.minor = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    # axis.title = element_blank(),
    panel.grid.major = element_blank()) -> P2


# Or by ggridges

RES %>% 
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
  filter(CONTRAST %in% c("CONTRAST_C", "CONTRAST_D")) %>%
  drop_na(MajorRNA) %>%
  mutate(HPF = ifelse(log2FoldChange > 0, "24 hpf", "110 hpf")) %>%
  mutate(HPF = factor(HPF, levels = c("24 hpf", "110 hpf"))) %>%
  mutate(CONTRAST = dplyr::recode_factor(CONTRAST, !!!recode_to)) %>%
  filter(MajorRNA %in% QUERIES) %>%
  select(log2FoldChange, MajorRNA, CONTRAST, HPF) %>%
  left_join(DB, by = "MajorRNA") %>%
  left_join(BT, by = "gene_id") %>%
  drop_na(Bilogical_pathway) %>%
  mutate(Label = paste(Bilogical_pathway, HPF, sep = "__")) %>%
  mutate(Label = factor(Label, levels = levels(DataViz$Label))) %>%
  drop_na(Label) %>%
  ggplot(aes(y = Label, x = abs(log2FoldChange), fill = CONTRAST, color = CONTRAST)) +
  ggridges::geom_density_ridges(
    jittered_points = T,
    position = ggridges::position_points_jitter(width = 0.05, height = 0),
    point_shape = '|', point_size = 0.5, point_alpha = 0.5, alpha = 0.2) +
  ggforce::facet_col(~ HPF, space = "free", scales = "free_y") +
  theme_bw(base_family = "GillSans", base_size = 10) +
  # scale_color_manual(values = c("black", "black")) +
  labs(x = "Log2FoldChange")

library(patchwork)

ps <- P1 + plot_spacer() + P2 + plot_layout(widths = c(4, -0.2, 3))

# ps

ggsave(ps, filename = 'Degree_and_log2fc.png', path = dir, width = 6, height = 4, device = png, dpi = 300)

# left panel Heatmap ----


recode_to <- structure(c("pH 8.0","pH 7.6", "24 hpf", "110 hpf"), names = c("CONTRAST_C", "CONTRAST_D", "CONTRAST_A", "CONTRAST_B"))

MajorRNALevs <- RES %>% 
  filter(CONTRAST %in% c("CONTRAST_C", "CONTRAST_D")) %>%
  mutate(CONTRAST = dplyr::recode_factor(CONTRAST, !!!recode_to)) %>%
  select(log2FoldChange, MajorRNA, CONTRAST) %>% 
  pivot_wider(names_from = CONTRAST, values_from = log2FoldChange) %>% arrange(desc(`pH 8.0`)) %>% pull(MajorRNA)

# Instead of log2fc use basemean per sample group
RES %>% 
  filter(CONTRAST %in% c("CONTRAST_C", "CONTRAST_D")) %>%
  drop_na(MajorRNA) %>%
  mutate(CONTRAST = dplyr::recode_factor(CONTRAST, !!!recode_to)) %>%
  filter(MajorRNA %in% QUERIES) %>%
  # select(log2FoldChange, MajorRNA, CONTRAST) %>%
  mutate(MajorRNA = factor(MajorRNA, levels = MajorRNALevs)) %>%
  # pivot_wider(names_from = CONTRAST, values_from = log2FoldChange) %>% arrange(desc(`pH 8.0`)) %>% 
  ggplot(aes(x = MajorRNA, y = log2FoldChange, color = CONTRAST, group = CONTRAST)) +
  geom_line(orientation = "x") +
  # geom_point(shape = 21, size = 1, stroke = 1.2, fill = "white") +
  scale_color_manual(values = c("gray86","black")) +
  theme_bw(base_size = 14, base_family = "GillSans") +
  theme(legend.position = "top",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.key.width = unit(0.2, "cm"),
    legend.key.height = unit(0.12, "cm")) +
  guides(color=guide_legend(title = "", nrow = 1), fill = "none") +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.5) +
  annotate("text", x = 10, y = 12, label = "24 hpf", family = "GillSans") +
  annotate("text", x = 10, y = -5, label = "110 hpf", family = "GillSans") +
  labs(x = "") -> P3


P3 <- P3 + 
  # scale_x_discrete(breaks = SIGNMIRS, labels = names(SIGNMIRS)) +
  # theme(axis.text.x = element_text(angle = 90, size = 8, hjust = 1, vjust = 0.5, color = "gray")) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

# Instead of log2fc use basemean per sample group

SIGNMIRS <- RES.p %>%  filter(!CONTRAST %in% c("CONTRAST_C", "CONTRAST_D")) %>% 
  distinct(MajorRNA, MirGeneDB_ID) %>% pull(MajorRNA, name = MirGeneDB_ID)


RES %>% 
  # filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
  filter(CONTRAST %in% c("CONTRAST_C", "CONTRAST_D")) %>%
  drop_na(MajorRNA) %>%
  mutate(CONTRAST = dplyr::recode_factor(CONTRAST, !!!recode_to)) %>%
  filter(MajorRNA %in% QUERIES) %>%
  select(baseMeanA, baseMeanB, MajorRNA, CONTRAST) %>%
  mutate(MajorRNA = factor(MajorRNA, levels = MajorRNALevs)) %>%
  pivot_longer(cols = c("baseMeanA","baseMeanB"), names_to = "HPF", values_to = "baseMean") %>%
  mutate(HPF = ifelse(HPF %in% "baseMeanA", "24 hpf", "110 hpf")) %>%
  mutate(HPF = factor(HPF, levels = c("24 hpf", "110 hpf"))) %>%
  mutate(trend = paste(MajorRNA,HPF, sep = "__")) %>%
  ggplot(aes(x = MajorRNA, y = log10(baseMean), color = CONTRAST, group = CONTRAST)) +
  geom_point(position = position_dodge(width = 0.5), shape = 21, size = 1, stroke = 1.2) +
  geom_line(aes(group = trend), position = position_dodge(width = 0.5)) +
  ggforce::facet_col(~ HPF, space = "free", scales = "free_y") +
  scale_color_manual(values = c("gray86","black")) +
  theme_bw(base_size = 14, base_family = "GillSans") +
  theme(legend.position = "top",
    strip.background = element_rect(fill = 'white', color = 'white'),
    strip.text = element_text(color = "black",hjust = 0),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.key.width = unit(0.2, "cm"),
    legend.key.height = unit(0.12, "cm")) +
  guides(color=guide_legend(title = "", nrow = 1), fill = "none") +
  scale_x_discrete(breaks = SIGNMIRS, labels = names(SIGNMIRS)) +
  theme(axis.text.x = element_text(angle = 90, size = 8, hjust = 1, vjust = 0.5, color = "gray"))



# 

DB %>% 
  filter(MajorRNA %in% QUERIES) %>%
  left_join(BT, by = "gene_id") %>%
  mutate(y_axis = Bilogical_pathway) %>% # or COG_category
  # mutate(y_axis = ifelse(is.na(y_axis), COG_category, y_axis)) %>%
  drop_na(y_axis) %>%
  distinct(MajorRNA, y_axis) %>%
  right_join(RESviz, by = "MajorRNA") %>%
  # count(MajorRNA, y_axis, Contrast, sort = T) %>%
  mutate(Label = paste(y_axis, Contrast, sep = "__")) %>%
  mutate(Label = factor(Label, levels = levels(DataViz$Label))) %>%
  mutate(MajorRNACol = ifelse(MajorRNA %in% SIGNMIRS, "black", "gray86")) %>%
  mutate(MajorRNA = factor(MajorRNA, levels = MajorRNALevs)) %>%
  ggplot() +
  scale_y_discrete("", labels = function(x) gsub("__.+$", "", x)) +
  geom_tile(aes(fill = factor(sign(log2FoldChange)), y = Label, x = MajorRNA), 
    color = 'white', size = 0.7, width = 1, height = 0.5) +
  theme_bw(base_size = 14, base_family = "GillSans") +
  ggforce::facet_col(~ Contrast, space = "free", scales = "free_y") +
  # scale_fill_manual(values = c("black", "gray86")) +
  guides(color=guide_legend(title = "", nrow = 1), fill = "none") +
  labs(y = "NOGS Category") +
  theme(
    # legend.position = "bottom",
    strip.background = element_rect(fill = 'white', color = 'white'),
    strip.text = element_text(color = "black",hjust = 0),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.key.width = unit(0.2, "cm"),
    legend.key.height = unit(0.12, "cm")) +
  theme(axis.text.x = element_text(angle = 90, size = 8, hjust = 1, vjust = 0.5, color = "gray"),
    axis.text.y = element_text(size = 8)) -> P4


P4 <- P4 + scale_x_discrete(breaks = SIGNMIRS, labels = names(SIGNMIRS))

library(patchwork)

PR <- P3/ plot_spacer() / P4 + plot_layout(heights = c(2.5, -0.9, 6.5))


PL <- P4 + plot_spacer() + P1 + plot_layout(widths = c(5, -0.1, 2))

PL

# Count numbers

DB %>% 
  # filter(MajorRNA %in% QUERIES[1:3]) %>%
  select(-Contrast) %>% distinct() %>%
  left_join(BT, by = "gene_id") %>%
  mutate(y_axis = Bilogical_pathway) %>% # or COG_category
  drop_na(COG_category) %>%
  mutate(COG_category = paste0(COG_category, ",",COG_name)) %>%
  mutate(y_axis = ifelse(is.na(y_axis), COG_category, y_axis)) %>%
  distinct(MajorRNA, y_axis) %>%
  right_join(RESviz, by = "MajorRNA") %>% 
  # filter(Contrast == "110 hpf") %>% 
  count(y_axis, sort = T) 


DB %>% 
  filter(MajorRNA %in% "UAUCACAGCCAGCUUUGAUGAGCU") %>%
  drop_na(STRINGID) %>%
  select(-Contrast) %>% distinct() %>%
  left_join(BT, by = "gene_id") %>% 
  # view()
  mutate(y_axis = Bilogical_pathway) %>% # or COG_category
  drop_na(COG_category) %>%
  mutate(COG_category = paste0(COG_category, ",",COG_name)) %>%
  mutate(y_axis = ifelse(is.na(y_axis), COG_category, y_axis)) %>%
  count(MajorRNA, y_axis, sort = T) 

  
DB %>% 
    filter(MajorRNA %in% QUERIES) %>%
  drop_na(Contrast) %>%
    select(-Contrast) %>% distinct() %>%
    left_join(BT, by = "gene_id") %>% 
  filter(grepl("ankyrin|glycoprotein|ubiquitin", `ENSEMBL description`)) %>% 
  count(gene_id, `ENSEMBL description`, STRINGID, Bilogical_pathway) %>%
  view()
  

# Intergenic / intragenic

# 
col_recode <- structure(c("#FFC107", "#2196F3","red"),
  names = c("Intergenic", "Intragenic", "Other ncRNA"))

DB %>%
  filter(MajorRNA %in% QUERIES) %>%
  left_join(BT, by = "gene_id") %>%
  view()
  distinct(gene_id, STRINGID, Bilogical_pathway, MajorRNA, MajorRNAID, biotype_best_rank) %>%
  count(Bilogical_pathway, biotype_best_rank, sort = T) %>%
  drop_na(Bilogical_pathway) %>%
  dplyr::rename("MIR_degree" = "n") %>%
  # group_by(Bilogical_pathway) %>% mutate(MIR_degree = MIR_degree/sum(MIR_degree)) %>%
  ggplot(aes(x = MIR_degree, y = Bilogical_pathway, fill = biotype_best_rank)) +
  geom_col(position = position_dodge(width = 1)) +
  scale_color_manual("",values = col_recode) +
  scale_fill_manual("",values = col_recode)