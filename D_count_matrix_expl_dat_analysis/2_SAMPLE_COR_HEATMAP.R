
# OLD VERSION, 
# 

library(tidyverse)

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

path <- "~/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out/"

# 
pHpalette <- c(`Experimental`="#ad1f1f", `Control`= "#4575b4")
#

out <- read_rds(paste0(path, "Counts_and_mtd.rds"))

dim(COUNTS <- out[[1]])

dim(MTD <- out[[2]])

# SAMPLE-SAMPLE CORRELATION HEATMAP ====

ncol(data <- edgeR::cpm(COUNTS+1))

# 1.1) Prior, calculate stats: (zero output)

library(rstatix)

# data %>%
#   as_tibble(rownames = "id") %>%
#   pivot_longer(-id) %>%
#   mutate(name = substr(name, 1,nchar(name)-1)) %>%
#   group_by(id, name) %>%
#   summarise(value = sum(value)) %>%
#   pivot_wider(names_from = name, values_from = value) %>%
#   # cor_test()
#   ungroup() %>%
#   select(-id)
#   

# cor_stats <- data %>%
#   cor_mat(method='pearson') %>%
#   # cor_get_pval()
#   cor_gather()

# summary(cor_stats$p)

sample_cor = cor(data, method='pearson', use='pairwise.complete.obs')

sample_dist = dist(sample_cor, method='euclidean')

hc_samples = hclust(sample_dist, method='complete')

hc_order <- hc_samples$labels[hc_samples$order]

# heatmap(sample_cor, col = cm.colors(12))

sample_cor %>% 
  as_tibble(rownames = 'LIBRARY_ID') %>%
  pivot_longer(cols = colnames(sample_cor), values_to = 'cor') %>%
  left_join(MTD) %>% 
  mutate(LIBRARY_ID = factor(LIBRARY_ID, levels = hc_order)) -> sample_cor_long

library(ggh4x)

sample_cor_long %>%
  ggplot(aes(x = LIBRARY_ID, y = name, fill = cor)) + 
  # geom_tile(color = 'white', size = 0.2) +
  geom_raster() +
  # geom_text(aes(label = cor), color = 'white') +
  scale_fill_viridis_c(option = "B", name = "", direction = -1) +
  # scale_x_discrete(position = 'top') +
  # scale_y_discrete(position = "right") +
  ggh4x::scale_x_dendrogram(hclust = hc_samples, position = 'top') +
  ggh4x::scale_y_dendrogram(hclust = hc_samples, position = "left") +
  guides(fill = guide_legend(title = "", nrow = 1)) +
  labs(x = '', y = '') +
  theme_bw(base_family = "GillSans", base_size = 10) +
  theme(legend.position = 'top',
    axis.text.x = element_text(angle = 90,
      hjust = -0.15, vjust = 1)) -> p

p

p <- p + theme(strip.background = element_rect(fill = 'white', color = 'white'),
  panel.border = element_blank(),
  panel.grid.minor = element_blank(),
  axis.ticks.x = element_blank(),
  panel.grid.major = element_blank())

ggsave(p, filename = 'HEATMAP_COR_MAT.png', path = path, width = 4, 
  height = 4.2, device = png, dpi = 300)

