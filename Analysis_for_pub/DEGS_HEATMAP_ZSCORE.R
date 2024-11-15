
# PLOT HEATMAP OF DEGS, USING Z SCORE

rm(list = ls());

if(!is.null(dev.list())) dev.off()


library(DESeq2)

library(tidyverse)

wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out/"

print(.LOCIDB <- read_tsv(paste0(wd, "/RNA_LOCATION_DB.tsv")) %>%   filter(SRNAtype == "miR") )


other_nc <- c("snRNA", "rRNA", "tRNA", "lncRNA", "srpRNA")
other_intra <- c("exon", "UTR")

.LOCIDB <- .LOCIDB %>% mutate(Chrom = gsub("^JALGQA010000", "", Chrom)) %>% 
  mutate(Chrom = gsub(".1$", "", Chrom)) %>% 
  # mutate(biotype_best_rank = ifelse(biotype_best_rank %in% other_nc, "Other ncRNA", biotype_best_rank)) %>%
  mutate(biotype_best_rank = ifelse(biotype_best_rank %in% other_nc, "Intergenic", biotype_best_rank)) %>%
  mutate(biotype_best_rank = ifelse(biotype_best_rank %in% other_intra, "Intragenic", biotype_best_rank))


.LOCIDB <- .LOCIDB %>% distinct(MajorRNA, biotype_best_rank)

dir <-  "~/Documents/MIRNA_HALIOTIS/MIRNA_PUB_2024"


print(DB <- read_tsv(file.path(dir, "LONGER_RELATIONAL_DB.tsv")))

DF <- DB %>% 
  drop_na(STRINGID) %>% 
  mutate(STRINGID = strsplit(STRINGID, ";")) %>%
  unnest(STRINGID) %>% 
  mutate(SMPID = ifelse(!is.na(SMPID), "(*) ", "")) %>% 
  mutate(STRINGID_label = paste0(SMPID, STRINGID)) %>%
  mutate(Contrast = strsplit(Contrast, ";")) %>%
  unnest(Contrast) %>%
  left_join(.LOCIDB) %>%
  filter(grepl("CONTRAST_C|CONTRAST_D", Contrast)) %>%
  distinct(biotype_best_rank, Contrast, MajorRNA, STRINGID_label, STRINGID) %>%
  dplyr::count(biotype_best_rank, Contrast, STRINGID_label, STRINGID) %>% 
  dplyr::rename("preferred_name" = "STRINGID", "degree" = "n") %>% 
  # dplyr::rename("protein_protein" = "degree") %>%
  arrange(desc(degree)) %>%
  mutate(preferred_name = factor(STRINGID_label, levels = unique(STRINGID_label)))

recode_c <- c(
  "CONTRAST_D_Control" = "24 hpf:pH 7.6",
  "CONTRAST_D_Competent" = "110 hpf:pH 7.6",
  "CONTRAST_C_Control" = "24 hpf:pH 8.0",
  "CONTRAST_C_Competent" = "110 hpf:pH 8.0",
  "CONTRAST_B_Low"= "110 hpf:pH 7.6",
  "CONTRAST_A_Low"= "24 hpf:pH 7.6",
  "CONTRAST_B_Control" = "110 hpf:pH 8.0")


DF %>%
  filter(grepl("CONTRAST_C", Contrast)) %>%
  drop_na() %>%
  mutate(Contrast = dplyr::recode_factor(Contrast, !!!recode_c)) %>%
  separate(Contrast, into = c("hpf", "pH"), sep = ":") %>%
  ggplot(aes(x = preferred_name, y = degree, group = pH, color = biotype_best_rank)) + 
  geom_point(shape = 13) + 
  # geom_step() +
  facet_grid(hpf ~.) +
  theme_classic(base_family = "GillSans", base_size = 12) +
  theme(legend.position = 'top', 
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 10))

dds <- read_rds(list.files(path = dir, pattern = "DDS_DESEQ2.rds", full.names = T))

colData <- data.frame(colData(dds))

WHICH_CONTRAST <- c("CONTRAST_A", "CONTRAST_B")

RES <- read_rds(paste0(dir, "/SEQUENCES_MERGED_DESEQ_RES_WGCNA.rds")) %>%
  mutate(Contrast = ifelse(sign(log2FoldChange) == -1, sampleB, sampleA)) %>%
  mutate(Contrast = paste0(CONTRAST,"_",Contrast))

RES.P <- RES %>% 
  filter( padj < 0.05  & abs(log2FoldChange) > 1) %>%
  filter(CONTRAST %in% WHICH_CONTRAST)

which_mirs <- RES.P %>% 
  # filter(CONTRAST_DE %in% WHICH_CONTRAST) %>%
  distinct(MajorRNA) %>% pull()

str(which_mirs)  



# divide the counts by the size factors or normalization factors before
barplot(cnts <- DESeq2::counts(dds, normalized = F, replaced = F)[which_mirs,])
barplot(cnts <- DESeq2::varianceStabilizingTransformation(round(cnts)))


z_scores <- function(x) {(x-mean(x))/sd(x)}

c1 <- cnts[,grepl("HR24", colnames(cnts))]
c2 <- cnts[,!grepl("HR24", colnames(cnts))]

M1 <- apply(c1, 1, z_scores)
M2 <- apply(c2, 1, z_scores)


dim(M <- rbind(M1, M2))

# heatmap(M <- apply(cnts, 1, z_scores))

h <- heatmap(t(M), col = cm.colors(12), keep.dendro = T)

hc_samples <- as.hclust(h$Colv)
plot(hc_samples)
hc_order <- hc_samples$labels[h$colInd]

hc_mirs <- as.hclust(h$Rowv)
plot(hc_mirs)
order_mirs <- hc_mirs$labels[h$rowInd]
plot(hc_mirs)
summary(M)


DATA <- M %>% 
  as_tibble(rownames = 'LIBRARY_ID') %>%
  pivot_longer(cols = colnames(M), values_to = 'fill', names_to = "MajorRNA") %>%
  left_join(colData, by = "LIBRARY_ID") %>%
  mutate(LIBRARY_ID = factor(LIBRARY_ID, levels = rev(hc_order)))


labels <- RES.P %>% select(MajorRNA, MirGeneDB_ID)
  # arrange(desc(Name)) %>%
  # group_by(Name) %>%
  # mutate(Label = Name, row_number = row_number(Label)) %>%
  # mutate(Label = factor(paste(Label, row_number, sep = "__"), 
  #   levels = rev(paste(Label, row_number, sep = "__"))))


labels <- structure(labels$MirGeneDB_ID, names = labels$MajorRNA)

labels <- labels[match(order_mirs, names(labels))]

identical(order_mirs, names(labels))

order_mirs <- labels

COG_names <- DB %>%
  filter(MajorRNA %in% names(order_mirs)) %>%
  drop_na(COG_name) %>% 
  dplyr::count(MajorRNA, COG_name, sort = T) %>%
  group_by(MajorRNA) %>% slice_head(n = 1) %>%
  summarise(across(COG_name, .fns = paste_col)) %>%
  pull(COG_name, name = MajorRNA)
  

# identical(names(COG_names), names(order_mirs)) # FALSE

COG_names <- COG_names[match(names(order_mirs), names(COG_names))]

identical(names(COG_names), names(order_mirs)) # TRUE

order_mirs <- structure(paste0(order_mirs, " (", COG_names, ")"), names = names(order_mirs))

# order_mirs <- paste0(order_mirs, " (", COG_names, ")")


lo = floor(min(DATA$fill))
up = ceiling(max(DATA$fill))
mid = (lo + up)/2

library(ggh4x)

DATA %>%
  dplyr::mutate(pH = dplyr::recode_factor(pH, !!!c(`Control` = "pH 8.0", `Low` = "pH 7.6"))) %>%
  dplyr::mutate(hpf = dplyr::recode_factor(hpf, !!!c(`24` = "24 hpf", `110` = "110 hpf"))) %>%
  mutate(LIBRARY_ID = gsub("HR", "", LIBRARY_ID)) %>%
  ggplot(aes(x = LIBRARY_ID, y = MajorRNA, fill = fill)) +  
  geom_tile(color = 'white', linewidth = 0.7, width = 1) +
  # facet_grid(~ hpf+pH, scales = "free_x") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
    na.value = "white", midpoint = mid, limit = c(lo, up),
    name = NULL) +
  ggh4x::scale_y_dendrogram(hclust = hc_mirs, position = "left", labels = NULL) +
  guides(y.sec = guide_axis_manual(labels = order_mirs, label_size = 8, label_family = "GillSans")) +
  theme_bw(base_size = 10, base_family = "GillSans") +
  labs(x = '', y = '') +
  theme(
    legend.position = "top",
    # axis.text.x = element_blank(),
    # axis.ticks.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7.5),
    # strip.background = element_rect(fill = 'grey89', color = 'white'),
    panel.border = element_blank(),
    plot.title = element_text(hjust = 0),
    plot.caption = element_text(hjust = 0),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank()) -> P

P <- P + 
  ggh4x::facet_nested(~ hpf+pH, nest_line = T, scales = "free_x", space = "free_x") +
  theme(
    strip.text.x = element_text(angle = 0, hjust = 1),
    strip.placement = "outside",
    strip.background = element_rect(fill = 'white', color = 'white')) 

P <- P + guides(
  fill = guide_colorbar(barwidth = unit(2.5, "in"),
    barheight = unit(0.07, "in"), 
    label.position = "top",
    title = "Z-score",
    title.position  = "top",
    title.theme = element_text(size = 7, family = "GillSans", hjust = 1),
    alignd = 0.8,
    ticks.colour = "black", ticks.linewidth = 0.25,
    frame.colour = "black", frame.linewidth = 0.25,
    label.theme = element_text(size = 7, family = "GillSans")))



recode_to <- c(`24` = "24 hpf", `110` = "110 hpf")

TOPDF <- DATA %>%
  distinct(LIBRARY_ID, hpf, pH) %>%
  dplyr::mutate(hpf = dplyr::recode_factor(hpf, !!!recode_to)) %>%
  # mutate(label = ifelse(pH %in% "Low", "*", "")) %>%
  # mutate(shape = ifelse(pH == "Low", "pH 7.6", "")) %>%
  mutate(y = 1)

topplot <- TOPDF %>%
  ggplot(aes(y = y, x = LIBRARY_ID)) + # color = as.factor(hpf)
  # geom_point(shape = 15, size = 2) +
  # geom_point(aes(shape = shape, color = shape), alpha = 0.5) +
  # geom_text(aes(label = label),  vjust = -0.7, hjust = 0.5, size = 1.5, family =  "GillSans", color = "#d73027") +
  ggh4x::scale_x_dendrogram(hclust = hc_samples, position = 'bottom', labels = NULL) +
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

psave <- P/ plot_spacer() /topplot + plot_layout(heights = c(3, -0.6, 0.1))

# psave

ggsave(psave, filename = 'DESEQ2_DEGS_HEATMAP_MIR.png', path = dir, 
  width = 6.5, height = 4.2, device = png, dpi = 300)

# Merge w/ COGs ------
paste_col <- function(x) { 
  x <- x[!is.na(x)] 
  x <- unique(sort(x))
  x <- paste(x, sep = ';', collapse = ';') 
  
  return(x)
  
}

dir <- "~/Documents/MIRNA_HALIOTIS/MIRNA_PUB_2024/"

print(DB <- read_tsv(file.path(dir, "LONGER_RELATIONAL_DB.tsv")))


DB %>%
  filter(MajorRNA %in% which_mirs) %>%
  drop_na(COG_name) %>% 
  group_by(MajorRNA) %>%
  summarise(across(COG_category, .fns = paste_col), n = n()) %>%
  left_join(DATA) %>%
  dplyr::mutate(pH = dplyr::recode_factor(pH, !!!c(`Control` = "pH 8.0", `Low` = "pH 7.6"))) %>%
  dplyr::mutate(hpf = dplyr::recode_factor(hpf, !!!c(`24` = "24 hpf", `110` = "110 hpf"))) %>%
  ggplot(aes(x = LIBRARY_ID, y = MajorRNA, fill = fill)) +  
  geom_tile(color = 'white', linewidth = 0.7, width = 1) +
  # facet_grid(~ hpf+pH, scales = "free_x") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
    na.value = "white", midpoint = mid, limit = c(lo, up),
    name = NULL) +
  ggh4x::scale_y_dendrogram(hclust = hc_mirs, position = "left", labels = NULL) +
  guides(y.sec = guide_axis_manual(labels = order_mirs, label_size = 8, label_family = "GillSans")) +
  theme_bw(base_size = 8, base_family = "GillSans") +
  ggh4x::facet_nested(~ hpf+pH, nest_line = T, scales = "free_x", space = "free_x") +
  theme(strip.background = element_rect(fill = 'white', color = 'white'))
  
  
DATA %>% 
  left_join(DB)
