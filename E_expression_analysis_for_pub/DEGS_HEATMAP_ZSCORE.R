
# PLOT HEATMAP OF DEGS, USING Z SCORE

library(DESeq2)


wd <-  "~/Documents/MIRNA_HALIOTIS/MIRNA_PUB_2024"

dds <- read_rds(list.files(path = wd, pattern = "DDS_DESEQ2.rds", full.names = T))

which_mirs <- RES.P %>% 
  # filter(CONTRAST_DE %in% WHICH_CONTRAST) %>%
  distinct(MajorRNA) %>% pull()

str(which_mirs)  

# divide the counts by the size factors or normalization factors before
barplot(cnts <- DESeq2::counts(dds, normalized = T, replaced = F)[which_mirs,])

barplot(cnts <- DESeq2::varianceStabilizingTransformation(round(cnts)))


z_scores <- function(x) {(x-mean(x))/sd(x)}

heatmap(M <- apply(M, 1, z_scores))



h <- heatmap(t(M), col = cm.colors(12), keep.dendro = T)

hc_samples <- as.hclust(h$Colv)
plot(hc_samples)
hc_order <- hc_samples$labels[h$colInd]

hc_mirs <- as.hclust(h$Rowv)
plot(hc_mirs)
order_mirs <- hc_mirs$labels[h$rowInd]

summary(M)


DATA <- M %>% 
  as_tibble(rownames = 'LIBRARY_ID') %>%
  pivot_longer(cols = colnames(M), values_to = 'fill', names_to = "Name") %>%
  left_join(.colData, by = "LIBRARY_ID") %>%
  mutate(LIBRARY_ID = factor(LIBRARY_ID, levels = rev(hc_order)))


DATA %>%
  dplyr::mutate(pH = dplyr::recode_factor(pH, !!!c(`Low` = "pH 7.6",`Control` = "pH 8.0"))) %>%
  dplyr::mutate(hpf = dplyr::recode_factor(hpf, !!!c(`24` = "24 hpf", `110` = "110 hpf"))) %>%
  ggplot(aes(x = LIBRARY_ID, y = Name, fill = fill)) +  
  geom_tile(linewidth = 0.2) +
  # facet_grid(~ hpf+pH, scales = "free_x") +
  ggsci::scale_fill_gsea(name = "", reverse = T) +
  ggsci::scale_color_gsea(name = "", reverse = T) +
  # scale_x_discrete(position = 'bottom') +
  # ggh4x::scale_x_dendrogram(hclust = hc_samples, position = "top", labels = NULL) 
  # guides(x.sec = guide_axis_manual(labels = hc_order, label_size = 5, label_family = "GillSans")) 
  ggh4x::scale_y_dendrogram(hclust = hc_mirs, position = "left", labels = NULL) +
  guides(y.sec = guide_axis_manual(labels = order_mirs, label_size = 5, label_family = "GillSans")) +
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

recode_to <- c(`24` = "24 hpf", `110` = "110 hpf")

TOPDF <- DATA %>%
  distinct(LIBRARY_ID, hpf, pH) %>%
  dplyr::mutate(hpf = dplyr::recode_factor(hpf, !!!recode_to)) %>%
  mutate(label = ifelse(pH %in% "Low", "*", "")) %>%
  mutate(y = 1)

topplot <- TOPDF %>%
  ggplot(aes(y = y, x = LIBRARY_ID)) + # color = as.factor(hpf)
  # geom_point(shape = 15, size = 2) +
  # geom_text(aes(label = label),  vjust = -0.7, hjust = 0.5, size = 1.5, family =  "GillSans", color = "#d73027") +
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

psave <- topplot/ plot_spacer() /P + plot_layout(heights = c(0.5, -0.5, 5))

# psave

ggsave(psave, filename = 'DESEQ2_OADEGS_HEATMAP_MIR.png', path = path, width = 3, height = 3.5, device = png, dpi = 300)
