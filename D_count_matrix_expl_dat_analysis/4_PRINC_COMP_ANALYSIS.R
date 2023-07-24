
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



# PCA =====

# ncol(data <- log2(count+1))

ncol(data <- log2(COUNTS+1))

PCA = prcomp(t(data), center = T, scale. = FALSE)

percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)

sd_ratio <- sqrt(percentVar[2] / percentVar[1])

PCAdf <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2])

# Visualize eigenvalues 
# Show the percentage of variances explained by each principal component.

barplot(PCA$sdev)

# library(mclust)

# d_clust <- mclust::Mclust(as.matrix(PCAdf), G=1:4, modelNames = mclust.options("emModelNames"))
# plot(d_clust)
# d_clust$BIC
# k <- d_clust$G
# names(k) <- d_clust$modelName

k <- 4

PCAdf %>% 
  dist(method = "euclidean") %>% 
  hclust() %>% 
  cutree(., k) %>% 
  as_tibble(rownames = 'LIBRARY_ID') %>% 
  mutate(cluster = paste0('C', value)) %>% 
  dplyr::select(-value) -> hclust_res


PCAdf %>%
  mutate(LIBRARY_ID = rownames(.)) %>%
  # mutate(g = substr(sample_id, 1,1)) %>%
  # left_join(hclust_res) %>%
  left_join(MTD) %>%
  mutate(sample_group = substr(LIBRARY_ID, 1,nchar(LIBRARY_ID)-1)) %>%
  ggplot(., aes(PC1, PC2)) +
  # coord_fixed(ratio = sd_ratio) +
  geom_abline(slope = 0, intercept = 0, linetype="dashed", alpha=0.5) +
  geom_vline(xintercept = 0, linetype="dashed", alpha=0.5) +
  ggforce::geom_mark_ellipse(aes(group = as.factor(sample_group)),
    fill = 'grey89', color = NA) +
  geom_point(size = 7, alpha = 0.7, aes(color = pH)) +
  geom_text( family = "GillSans",
    mapping = aes(label = paste0(hpf, " hpf")), size = 2.5) +
  labs(caption = '') +
  ylim(-250, 250) +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  scale_color_manual("", values = rev(pHpalette)) +
  scale_fill_manual("", values = rev(pHpalette)) +
  theme_classic(base_family = "GillSans", base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = 'top') -> pcaplot

# pcaplot

ggsave(pcaplot, 
  filename = "PCA.png", path = path, 
  width = 5, height = 5, device = png, dpi = 300)

# 

