
# 1) GENERATE SOME PLOTS AS HEATMAP OF RELEVANT TERMS FROM EACH GROUP COMPARED:
# 2) EX. MAKE HEATMAP OF KNOWN-NOVEL MIRS W SEMANTIC SIMILARITY GROUPING (SPLIT BY DOWN AND UP EXPRESSED)
# 2.1) TOP LEGEND INCLUDE EXPERIMENTAL DESIGN
# 2.2) LEFT OR RIGHT LABEL INCLUDE SEMANTIC SIMILARITY GROUPING

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

wd <- "~/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out/"

.UPSETDF <- read_rds(paste0(wd, "UPSETDF.rds"))

EXCLUSIVE_IDS <- .UPSETDF[[1]] %>% distinct(Name) %>% pull()


wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/FUNCTIONAL_MIR_ANNOT/"

# 
# .OUT <- read_rds(paste0(wd, "DESEQ2TOPGO.rds"))
# 
# TOPGO_EXCLUSIVE <- .OUT$EXCLUSIVE
# 
# TOPGO_INTERSECTED <- .OUT$INTERSECTED

RES.P <- read_tsv(paste0(wd, "DESEQ_RES_P.tsv")) 

wd <- "~/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out/"

.COUNT <- read_rds(paste0(wd, "COUNT.rds"))

recode_to <- structure(c("Ctrl pH", "Low pH"), names = c("Control", "Low"))

.colData <- read_tsv(paste0(wd, "METADATA.tsv")) %>%
  dplyr::mutate(hpf = paste0(hpf, " HPF")) %>%
  dplyr::mutate(hpf = factor(hpf, levels = c("24 HPF", "110 HPF"))) %>%
  dplyr::mutate(pH = dplyr::recode_factor(pH, !!!recode_to))

# 3) HEATMAP OF DATA:
# ALL OR EXCLUSIVE?

# str(query.ids <- UPSETDF %>% filter( n == 1) %>% distinct(Name) %>% pull())

RES.P %>% filter(abs(log2FoldChange) > 1) %>% count(WRAP, CONTRAST, SIGN)

RES.P <- RES.P %>% filter(abs(log2FoldChange) > 1)


# UPEXPRESSED 

str(query.ids <- RES.P %>% filter(SIGN == "Up") %>% distinct(Name) %>% pull()) # 46

sum(KEEP <- rownames(.COUNT) %in% EXCLUSIVE_IDS) # OR query.ids

nrow(COUNT <- .COUNT[KEEP,])

COUNT <- DESeq2::varianceStabilizingTransformation(COUNT)

# COUNT <- edgeR::cpm(COUNT)

colSums(COUNT)

# COUNT <- COUNT %>% 
#   as_tibble(rownames = "Name") %>%
#   left_join(RES.P %>% distinct(Name, Family)) 


sample_dist = dist(t(COUNT), method='euclidean')

hc_samples = hclust(sample_dist, method='complete')

hc_order <- hc_samples$labels[hc_samples$order]

srna_dist = dist(COUNT, method='euclidean')

hc_srna = hclust(srna_dist, method='complete')

srna_order <- hc_srna$labels[hc_srna$order]

recode_to <- RES.P %>% filter(Name %in% srna_order) %>% distinct(Name, Family)

recode_to <- structure(recode_to$Family, names = recode_to$Name)

identical(sort(names(recode_to)),sort(srna_order))

srna_order <- recode_to[match(srna_order, names(recode_to))]

identical(names(srna_order),  hc_srna$labels[hc_srna$order])


sum(tag <- srna_order %in% EXCLUSIVE_IDS)

srna_order[tag] <- paste0(srna_order[tag], "*") 

# heatmap(COUNT, col = cm.colors(12))

COUNT %>% 
  as_tibble(rownames = 'Name') %>%
  pivot_longer(-Name, names_to = "LIBRARY_ID") %>%
  left_join(.colData) %>%
  mutate(LIBRARY_ID = factor(LIBRARY_ID, levels = hc_order)) -> COUNT_LONG


COUNT_LONG %>%
  right_join(distinct(RES.P, Name, Family), by = "Name") %>% 
  mutate(Family = ifelse(!grepl("Cluster", Family), toupper(Family), Family)) %>%
  group_by(LIBRARY_ID, Family) %>% 
  summarise(value = sum(value)) 

library(ggh4x)

p1 <- COUNT_LONG %>%
  ggplot(aes(x = LIBRARY_ID, y = Name, fill = value)) + 
  geom_tile(color = 'white', linewidth = 0.2) +
  scale_fill_viridis_c(option = "B", name = "Log2", direction = -1, na.value = "white") +
  # ggsci::scale_fill_gsea(name = "", reverse = T, na.value = "white") +
  ggh4x::facet_nested( ~ hpf+pH, nest_line = F, scales = "free", space = "free") +
  # scale_x_discrete(position = 'top') +
  # ggh4x::scale_x_dendrogram(hclust = hc_samples, position = 'top') +
  ggh4x::scale_y_dendrogram(hclust = hc_srna, position = "left", labels = NULL) +
  guides(y.sec = guide_axis_manual(labels = srna_order, label_size = 3.5)) +
  theme_bw(base_size = 10, base_family = "GillSans") +
  labs(x = '', y = '') +
  theme(
    legend.position = "top",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    # axis.text.x = element_text(angle = 90, hjust = -0.15, vjust = 1),
    strip.background = element_rect(fill = 'grey89', color = 'white'),
    panel.border = element_blank(),
    plot.title = element_text(hjust = 0),
    plot.caption = element_text(hjust = 0),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank()) 

# DOWN EXP

str(query.ids <- RES.P %>% filter(SIGN == "Down") %>% distinct(Name) %>% pull()) # 59

sum(KEEP <- rownames(.COUNT) %in% query.ids)

nrow(COUNT <- .COUNT[KEEP,])

COUNT <- DESeq2::varianceStabilizingTransformation(COUNT)

# COUNT <- edgeR::cpm(COUNT)

colSums(COUNT)

# COUNT <- COUNT %>% 
#   as_tibble(rownames = "Name") %>%
#   left_join(RES.P %>% distinct(Name, Family)) 


sample_dist = dist(t(COUNT), method='euclidean')

hc_samples = hclust(sample_dist, method='complete')

hc_order <- hc_samples$labels[hc_samples$order]

srna_dist = dist(COUNT, method='euclidean')

hc_srna = hclust(srna_dist, method='complete')

srna_order <- hc_srna$labels[hc_srna$order]

recode_to <- RES.P %>% filter(Name %in% srna_order) %>% distinct(Name, Family)

recode_to <- structure(recode_to$Family, names = recode_to$Name)

identical(sort(names(recode_to)),sort(srna_order))

srna_order <- recode_to[match(srna_order, names(recode_to))]

identical(names(srna_order),  hc_srna$labels[hc_srna$order])

sum(tag <- srna_order %in% EXCLUSIVE_IDS)

srna_order[tag] <- paste0(srna_order[tag], "*") 

# heatmap(COUNT, col = cm.colors(12))

COUNT %>% 
  as_tibble(rownames = 'Name') %>%
  pivot_longer(-Name, names_to = "LIBRARY_ID") %>%
  left_join(.colData) %>%
  mutate(LIBRARY_ID = factor(LIBRARY_ID, levels = hc_order)) -> COUNT_LONG

library(ggh4x)

# scale_fill_gradient2(low = "red", mid = "white", high = "blue", 
#   midpoint = 0, na.value = "white",
#   breaks = seq(-5, 0, 5), 
#   limits = c(-5, 5)) +
  
p2 <- COUNT_LONG %>%
  ggplot(aes(x = LIBRARY_ID, y = Name, fill = log2(value))) + 
  geom_tile(color = 'white', linewidth = 0.2) +
  scale_fill_viridis_c(option = "B", name = "Log2", direction = -1, na.value = "white") +
  ggh4x::facet_nested( ~ hpf+pH, nest_line = F, scales = "free", space = "free") +
  # scale_x_discrete(position = 'top') +
  # ggh4x::scale_x_dendrogram(hclust = hc_samples, position = 'top') +
  ggh4x::scale_y_dendrogram(hclust = hc_srna, position = "left", labels = NULL) +
  guides(y.sec = guide_axis_manual(labels = srna_order, label_size = 3.5)) +
  theme_bw(base_size = 10, base_family = "GillSans") +
  labs(x = '', y = '') +
  theme(
    legend.position = "none",
    # axis.text.y = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = -0.15, vjust = 1),
    # strip.background = element_rect(fill = 'grey89', color = 'white'),
    strip.background = element_blank(), 
    strip.text = element_blank(),
    panel.border = element_blank(),
    plot.title = element_text(hjust = 0),
    plot.caption = element_text(hjust = 0),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank()) 


library(patchwork)


p <- p1 / plot_spacer() / p2 + plot_layout(heights = c(3, -0.3 , 3))

# p <- p1 / p2 + patchwork::plot_layout(heights = c(1.2,1.2))


ggsave(p, filename = "DESEQ2HEATMAP.png", 
  path = wd, width = 3.5, height = 10, device = png, dpi = 300)


# SINGLE HEATMAP

str(query.ids <- RES.P %>%
  # filter(!grepl("Cluster", Family)) %>%
  distinct(Name) %>% pull()) # 65

sum(KEEP <- rownames(.COUNT) %in% query.ids)

nrow(COUNT <- .COUNT[KEEP,])

COUNT <- DESeq2::varianceStabilizingTransformation(COUNT)

# COUNT <- edgeR::cpm(COUNT)

colSums(COUNT)

# COUNT <- COUNT %>% 
#   as_tibble(rownames = "Name") %>%
#   left_join(RES.P %>% distinct(Name, Family)) 


sample_dist = dist(t(COUNT), method='euclidean')

hc_samples = hclust(sample_dist, method='complete')

hc_order <- hc_samples$labels[hc_samples$order]

srna_dist = dist(COUNT, method='euclidean')

hc_srna = hclust(srna_dist, method='complete')

srna_order <- hc_srna$labels[hc_srna$order]

recode_to <- RES.P %>% filter(Name %in% srna_order) %>% distinct(Name, Family)

recode_to <- structure(recode_to$Family, names = recode_to$Name)

identical(sort(names(recode_to)),sort(srna_order))

srna_order <- recode_to[match(srna_order, names(recode_to))]

identical(names(srna_order),  hc_srna$labels[hc_srna$order])

sum(tag <- srna_order %in% EXCLUSIVE_IDS)

srna_order[tag] <- paste0(srna_order[tag], "*") 

# heatmap(COUNT, col = cm.colors(12))

COUNT %>% 
  as_tibble(rownames = 'Name') %>%
  pivot_longer(-Name, names_to = "LIBRARY_ID") %>%
  left_join(.colData) %>%
  mutate(LIBRARY_ID = factor(LIBRARY_ID, levels = hc_order)) -> COUNT_LONG

library(ggh4x)

p3 <- COUNT_LONG %>%
  ggplot(aes(x = LIBRARY_ID, y = Name, fill = log2(value))) + 
  geom_tile(color = 'white', linewidth = 0.2) +
  scale_fill_viridis_c(option = "B", name = "", direction = -1, na.value = "white") +
  # ggh4x::facet_nested( ~ hpf+pH, nest_line = F, scales = "free", space = "free") +
  # scale_x_discrete(position = 'top') +
  ggh4x::scale_x_dendrogram(hclust = hc_samples, position = 'top', labels = NULL) +
  ggh4x::scale_y_dendrogram(hclust = hc_srna, position = "left", labels = NULL) +
  guides(y.sec = guide_axis_manual(labels = srna_order, label_size = 3.5),
    x.sec = guide_axis_manual(labels = hc_order, label_size = 10)) +
  theme_bw(base_size = 10, base_family = "GillSans") +
  labs(x = '', y = '') +
  theme(
    legend.position = "top",
    # axis.text.y = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1.1),
    # strip.background = element_rect(fill = 'grey89', color = 'white'),
    strip.background = element_blank(), 
    strip.text = element_blank(),
    panel.border = element_blank(),
    plot.title = element_text(hjust = 0),
    plot.caption = element_text(hjust = 0),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank()) 

p3

# 2) LEFT OR RIGHT LABEL INCLUDE SEMANTIC SIMILARITY GROUPING

RES.P %>% distinct(GO.ID)


# Similarity grouping ====

library(GOSemSim)

# hsGO <- godata('org.Hs.eg.db', ont="BP")

# rbind(TOPGO_INTERSECTED, TOPGO_EXCLUSIVE)


str(GO.ID <- rbind(TOPGO_INTERSECTED, TOPGO_EXCLUSIVE) %>% distinct(GO.ID) %>% pull())

# org.Dr.eg.db # ZEBRAFISH
# org.Ce.eg.db # WORM
# org.Xl.eg.db # XENOPUS

WHICH_DB <- "org.Ce.eg.db"

BiocManager::install(WHICH_DB)

GOSEM <- GOSemSim::godata(WHICH_DB, ont="BP")

# hsGO <- GOSemSim::godata('org.Hs.eg.db', ont="BP")

# hsGO <- read_rds(paste0(wd, '/hsGO_BP.rds'))

termSim <- GOSemSim::termSim(GO.ID, GO.ID,  semData = GOSEM, method = "Wang")

# heatmap(termSim)

cmds <- termSim %>% 
  # use distance metric
  dist(method = "euclidean") %>%
  # compute cMDS
  cmdscale() %>%
  data.frame() %>%
  rownames_to_column(var= 'GO.ID')

cutree <- termSim %>% 
  # use distance metric
  dist(method = "euclidean") %>%
  # compute cMDS
  hclust() %>%
  cutree(., 6)

# Sanity check

identical(cmds$GO.ID, names(cutree))

head(cmds <- cbind(cmds, cutree) %>% mutate(cutree = as.factor(cutree)))

# AS DIMENSIONAL SCATTERPLOT
TOPGO <- rbind(TOPGO_INTERSECTED, TOPGO_EXCLUSIVE) %>% left_join(cmds)

# separate(CONTRAST, into = c("WRAP","CONTRAST"), sep = "[|]")

TEXTDF <- TOPGO %>% group_by(CONTRAST, SIGN) %>% arrange(p.adj.ks) %>% sample_n(1)

ggplot(data = TOPGO, aes(X1, X2)) +
  ggh4x::facet_nested(SIGN ~ CONTRAST, nest_line = F, scales = "free", space = "free") +
  ggforce::geom_mark_circle(aes(group = SIGN, label = SIGN), fill = 'grey89', color = NA) +
  geom_point(data = TOPGO, mapping = aes(X1, X2, color = cutree, alpha = p.adj.ks, size = Annotated)) +
  ggrepel::geom_text_repel(data = TEXTDF, mapping = aes(X1, X2, label = Term), 
    max.overlaps = 20, color = 'black', family = "GillSans") +
  theme_classic(base_size = 16, base_family = "GillSans") +
  labs(x = 'Dimension 1', y = 'Dimension 2') 


# ENRICHMENT ANALYSIS: =====

.bioc_packages <- c('biomaRt','GOSemSim',"org.Hs.eg.db") # DESeq2'

sapply(c(.bioc_packages), require, character.only = TRUE)

# PLOT ====
# POSITIVE LFC == UP EXPRESSED IN EXPERIMENTAL
# NEGATIVE LFC == UP EXPRESSED IN CONTROL

# recode_to <- paste(LETTERS[1:4], ")", sep = "")

recode_to <- c("24 HPF: Ctrl|pH", "110 HPF: Ctrl|pH", "Ctrl: 24|110 HPF", "pH: 24|110")
recode_to <- structure(recode_to, names = CONTRAST)

# res.p <- prep_DE_data(allRes, alpha = 0.05, lfcThreshold = 2) %>% drop_na(cc)

sigfc <- "Sign and FC";pv <- "Sign";fc <- "FC"

colors_fc <- c("red2",  "#4169E1", "forestgreen", "grey30")

colors_fc <- structure(colors_fc, names = c(sigfc, pv, fc, "NS"))

allRes$Color <- 'NS'
allRes[which(abs(allRes$log2FoldChange) > 2), 'Color'] <- fc
allRes[which(abs(allRes$padj) <= 0.05), 'Color'] <- pv
allRes[which(allRes$padj <= 0.05 & abs(allRes$log2FoldChange) > 2), 'Color'] <- sigfc

colors_fc <- structure(colors_fc, names = c(sigfc, pv, fc, "NS"))

p <- allRes  %>%
  dplyr::mutate(CONTRAST = dplyr::recode_factor(CONTRAST, !!!recode_to))
ggplot(aes(y = -log10(pvalue), x = log2FoldChange, color = Color)) +
  geom_point()  

p <- p +
  scale_color_manual(name = "", values = colors_fc) +
  labs(x= expression(Log[2] ~ "Fold Change"), y = expression(-Log[10] ~ "P")) +
  theme_bw(base_family = "GillSans", base_size = 18) +
  theme(legend.position = "top")  +
  geom_abline(slope = 0, intercept = -log10(0.05), linetype="dashed", alpha=0.5) +
  geom_vline(xintercept = 0, linetype="dashed", alpha=0.5) 

p + facet_grid(~ CONTRAST)

allRes %>% 
  mutate(g = sign(log2FoldChange)) %>% 
  dplyr::count(Color, g, CONTRAST) %>%
  # mutate(g = ifelse(g == "1", "CON_CANCER", "SIN_CANCER")) %>%
  filter(Color != 'NS') %>%
  group_by(g, CONTRAST) %>%
  mutate(pct = n / sum(n)) %>%
  ggplot() +
  facet_grid(g ~ .) +
  geom_col(aes(x = CONTRAST, y = pct, fill = Color), width = 0.5) +
  scale_fill_manual(name = "", values = colors_fc[-4]) +
  theme_bw(base_family = "GillSans", base_size = 18) +
  labs(y = '% Transcripts', x = '') +
  theme(legend.position = 'top') + coord_flip()






