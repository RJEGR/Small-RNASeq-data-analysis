
rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

# MIRDEEP2 -----

path <- '~/Documents/MIRNA_HALIOTIS/CLEAN_INPUT/'

mtd <- read_tsv(list.files(path = path, pattern = 'METADATA', full.names = T))


fileName <- list.files(path = path, pattern = 'miRNAs_expressed_all_samples', full.names = T)

count <- read_tsv(file = fileName)

whichCols <- !grepl('norm', names(count))

count <- count[whichCols]

count %>% distinct(`#miRNA`)

count <- count %>% filter(total > 0)

names(count)

summary <- count[,1:4]

count <- count[,5:16]

# PCA

# PCA <- prcomp(t(log2(count+1)), scale. = FALSE) 

PCA <- prcomp(t(edgeR::cpm(count)), scale. = FALSE) 

percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)

sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dtvis <- data.frame(PC1 = PCA$x[,1],  PC2 = PCA$x[,2]) %>% as_tibble(rownames = 'colName')
# mtd %>% distinct(group, .keep_all = T))

dtvis <- mtd %>% left_join(dtvis) %>% mutate(label = paste0(hpf, "-", pH))

# g <- dtvis %>% filter(grepl('Developmetal', Sample.type )) %>% pull(stage)
# n <- length(unique(g))
# grid.col <- c(ggsci::pal_d3()(10), ggsci::pal_aaas()(n-10))
# grid.col <- structure(grid.col, names = stagesLev)

unique(dtvis$label)

sortLabel <- c("24-Control","24-low","110-Control", "110-low")


dtvis %>%
  mutate(label = factor(label,levels = sortLabel)) %>%
  ggplot(., aes(PC1, PC2)) +
  geom_point(aes(color = pH), size = 10, alpha = 0.9) +
  geom_abline(slope = 0, intercept = 0, linetype="dashed", alpha=0.5) +
  geom_vline(xintercept = 0, linetype="dashed", alpha=0.5) +
  geom_text(aes(label = as.factor(hpf)), family = "GillSans") +
  ggforce::geom_mark_hull(aes(group = as.factor(label)), fill = "grey", color = NA) +
  # color = "black", linetype = 2
  labs(caption = 'cpm (x)') +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme_bw(base_family = "GillSans", base_size = 16) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = 'top') +
  # scale_color_manual(values = grid.col) +
  # coord_fixed(ratio = sd_ratio) +
  guides(fill = "none") 
# scale_color_manual('', values = getPalette) -> psave
# guides(color = FALSE, fill = FALSE)

ggsave(psave, filename = 'digitalExp_sRNA_biogenesis_dev_PCA.png', path = path, width = 12, height = 6.5)

# SHORTSTACKS -----


path2 <- "/Users/cigom/Documents/MIRNA_HALIOTIS/SHORTSTACKS"
fileName <- list.files(path = path2, pattern = 'Counts.txt', full.names = T)

nrow(count <- read_tsv(file = fileName)) # 918,958

nrow(x <- count %>% filter(main > 1)) # 543,291

keep <- grepl('clean.newid.subset', names(x))

x <- x[keep]

names(x) <- gsub(".clean.newid.subset", "", names(x))

PCA <- prcomp(t(edgeR::cpm(x)), scale. = FALSE) 

percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)

sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dtvis <- data.frame(PC1 = PCA$x[,1],  PC2 = PCA$x[,2]) %>% as_tibble(rownames = 'File')

# mtd %>% distinct(group, .keep_all = T))

dtvis <- mtd %>% left_join(dtvis) %>% mutate(label = paste0(hpf, "+", pH))

# sortLabel <- c("24+Control","24+low","110-Control", "110-low")


psave <- dtvis %>%
  # mutate(label = factor(label,levels = sortLabel)) %>%
  ggplot(., aes(PC1, PC2)) +
  geom_point(aes(color = pH), size = 10, alpha = 0.9) +
  geom_abline(slope = 0, intercept = 0, linetype="dashed", alpha=0.5) +
  geom_vline(xintercept = 0, linetype="dashed", alpha=0.5) +
  geom_text(aes(label = as.factor(hpf)), family = "GillSans") +
  ggforce::geom_mark_ellipse(aes(group = as.factor(label)), fill = "grey", color = NA) +
  # color = "black", linetype = 2
  labs(caption = 'cpm (x)') +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme_bw(base_family = "GillSans", base_size = 16) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = 'top') +
  # scale_color_manual(values = grid.col) +
  # coord_fixed(ratio = sd_ratio) +
  guides(fill = "none") 

ggsave(psave, filename = 'PCA.png', path = path2, device = png)




