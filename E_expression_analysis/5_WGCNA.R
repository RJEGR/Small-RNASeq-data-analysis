
# LETS CLUSTER MIRS TO OVERVIEW PANEL OF RESPONSE TO ACIDIFICATION AND DEVELOPMENTAL STAGES:
# EX. https://doi.org/10.1038/s41598-020-75945-2 (LOOK GENES HUBS SEACH M&M)


library(WGCNA)
library(flashClust)
library(tidyverse)


rm(list = ls());

if(!is.null(dev.list())) dev.off()

path <- "~/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out/"

count_f <- list.files(path = path, pattern = "Counts.txt", full.names = T)

count <- read_tsv(count_f)


which_cols <- count %>% select_if(is.double) %>% colnames()


datExpr <- count %>% select(all_of(which_cols))

datExpr <- as(datExpr, "matrix")

rownames(datExpr) <- count$Name

.datExpr <- datExpr

datExpr <- log2(datExpr+1)

datExpr <- t(datExpr) # log2(count+1) # 

str(datExpr)

cat("\n:::::\n")

gsg = goodSamplesGenes(datExpr, verbose = 3)

gsg$allOK

if (!gsg$allOK) {
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse= ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse=", ")))
  datExpr= datExpr[gsg$goodSamples, gsg$goodGenes]
}


# Simulated-05-NetworkConstruction.pdf (Horvath and Langfelder, 2011)

# When you have a lot of genes use the following code
k=softConnectivity(datE=datExpr,power=6)
# Plot a histogram of k and a scale free topology plot
sizeGrWindow(10,5)
par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k, main="Check scale free topology\n")

# Debemos señalar que no es necesario que una red satisfaga una topología libre de escala; La topología libre de escala puede no satisfacerse si los datos están compuestos por grupos de muestras globalmente muy distintos (por ejemplo, diferentes tipos de tejidos).

# DEBIDO A LA CANTIDAD DE DATOS MASIVOS. NO PODEMOS GENERAR: Average linkage hierachical clustering with adjacency-based dissimilarity


# 1) Detect max power ----

max_power <- 30

powers = c(c(1:10), seq(from = 10, to = max_power, by=1)) 

#powers = unique(powers)

allowWGCNAThreads()

cor_method =  "cor" # by default WGCNA::cor(method =  'pearson') is used, "bicor"

corOptionsList = list(use ='p') # maxPOutliers = 0.05, blocksize = 20000

sft = pickSoftThreshold(datExpr, 
  powerVector = powers, 
  corFnc = cor_method,
  corOptions = corOptionsList,
  verbose = 5, 
  networkType = "signed")

saveRDS(sft, file = paste0(path, 'SoftThreshold_',cor_method, '.rds'))

sft <- read_rds(paste0(path, 'SoftThreshold_',cor_method, '.rds'))

soft_values <- abs(sign(sft$fitIndices[,3])*sft$fitIndices[,2])

soft_values <- round(soft_values, digits = 2)

hist(soft_values)

power_pct <- quantile(soft_values, probs = 0.95)

softPower <- sft$fitIndices[,1][which(soft_values >= power_pct)]

meanK <- sft$fitIndices[softPower,5]

hist(sft$fitIndices[,5])

softPower <- min(softPower)

cat("\nsoftPower value", softPower, '\n')


title1 = 'Scale Free Topology Model Fit,signed R^2'
title2 = 'Mean Connectivity'

caption = paste0("Lowest power for which the scale free topology index reaches the ", power_pct*100, " %")

sft$fitIndices %>% 
  mutate(scale = -sign(slope)*SFT.R.sq) %>%
  select(Power, mean.k., scale) %>% pivot_longer(-Power) %>%
  mutate(name = ifelse(name %in% 'scale', title1, title2)) %>%
  ggplot(aes(y = Power, x = value)) +
  geom_text(aes(label = Power), size = 2) +
  geom_abline(slope = 0, intercept = softPower, linetype="dashed", alpha=0.5) +
  # geom_vline(xintercept = min(meanK), linetype="dashed", alpha=0.5) +
  labs(y = 'Soft Threshold (power)', x = '', 
    caption = caption) +
  facet_grid(~name, scales = 'free_x', switch = "x") +
  # scale_x_continuous(position = "top") +
  theme_light(base_family = "GillSans",base_size = 16) -> psave

ggsave(psave, path = path, filename = 'SoftThreshold.png', width = 7, height = 4)

# The soft thresholding, 
# is a value used to power the correlation of the genes to that threshold. 
# The assumption on that by raising the correlation to a power will reduce the noise of 
# the correlations in the adjacency matrix


# 2) Construct a gene co-expression matrix and generate modules ----
# Using Blockwise construction

# Call the network topology analysis function

cor_method = 'cor'

filename <- paste0(path, 'SoftThreshold_',cor_method, '.rds')

sft <- readRDS(filename)

soft_values <- abs(sign(sft$fitIndices[,3])*sft$fitIndices[,2])

soft_values <- round(soft_values, digits = 2)

power_pct <- quantile(soft_values, probs = 0.95)

softPower <- sft$fitIndices[,1][which(soft_values >= power_pct)]

softPower <- min(softPower)

allowWGCNAThreads()

wd <- paste0(path, Sys.Date())

system(paste0('mkdir ', wd))

setwd(wd)

# RUN: 

bwnet <- blockwiseModules(datExpr, 
  maxBlockSize = 5000,
  power = 5, # softPower, 
  TOMType = "signed", 
  networkType = "signed",
  minModuleSize = 30,
  corType = "bicor",
  reassignThreshold = 0, 
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  saveTOMs = TRUE,
  saveTOMFileBase = "TOM-blockwise",
  verbose = 3)

saveRDS(bwnet, "bwnet.rds")

# 3) Output exploratory ----


bwnet <- readRDS(paste0(path, "2023-06-26/bwnet.rds"))

bwmodules = labels2colors(bwnet$colors)

names(bwmodules) <- names(bwnet$colors)

table(bwmodules)


# (TEST) Plot the dendrogram and the module colors underneath for block 1

# i <- 8
# 
# plotDendroAndColors(bwnet$dendrograms[[i]], bwmodules[bwnet$blockGenes[[i]]],
#   "Module colors", main = "Gene dendrogram and module colors in block i",
#   dendroLabels = FALSE, hang = 0.03,
#   addGuide = TRUE, guideHang = 0.05)



reads <- rowSums(.datExpr)

Total <- sum(reads)

# Sanity check:

identical(names(colSums(datExpr)), names(bwmodules))

data.frame(reads, bwmodules) %>% 
  as_tibble() %>% 
  group_by(bwmodules) %>% 
  summarise(n = n(), reads = sum(reads)) %>%
  dplyr::rename('module' = 'bwmodules') -> stats

# Prep binary datTraits

datTraits <- gsub(".clean.newid.subset", "", rownames(datExpr))

HR11076 <- grepl('HR11076', datTraits)

HR1108 <- grepl('HR1108', datTraits)

HR2476 <- grepl('HR2476', datTraits)

HR248 <- grepl('HR248', datTraits)

datTraits <- data.frame(HR11076, HR1108, HR2476, HR248)

datTraits <- 1*datTraits

rownames(datTraits) <- rownames(datExpr)

# Recalculate MEs with color labels

MEs0 = moduleEigengenes(datExpr, bwmodules)$eigengenes

MEs = orderMEs(MEs0)

names(MEs) <- str_replace_all(names(MEs), '^ME', '')

moduleTraitCor = cor(MEs, datTraits, use= "p")

moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(datTraits))


moduleTraitCor %>% as_tibble(rownames = 'module') %>% 
  pivot_longer(-module, values_to = 'moduleTraitCor') -> df1

moduleTraitPvalue %>% as_tibble(rownames = 'module') %>% 
  pivot_longer(-module, values_to = 'corPvalueStudent') %>%
  right_join(df1) -> df1

hclust <- hclust(dist(moduleTraitCor), "complete")

# # BASED ON THE dendogram from the heatmap, WE REALIZE THAN "tan" (2n clade), "magenta" (2n clade), "greenyellow" (1st clade) are within same cluster than grey module in the dendo.


# plot(hclust)

# SPLIT CLUSTERS BY MIRS/PIRS/SIRS

SRNAS <- read_rds(paste0(path, "KNOWN_CLUSTERS_MIRS_PIRS.rds"))

str(which_pirs <- SRNAS %>% filter(grepl("piR", KnownRNAs)) %>% distinct(Name) %>% pull())
str(which_mirs <- SRNAS %>% filter(!grepl("piR", KnownRNAs)) %>% distinct(Name) %>% pull())


bwmodules %>% 
  as_tibble(., rownames = 'Name') %>%
  mutate(biotype = ifelse(Name %in% which_mirs, 'miR', 
    ifelse(Name %in% which_pirs, 'piR', 'siRs'))) %>%
  dplyr::rename('module' = 'value') -> bwModuleCol

# ADD COUNT INFO:
# reads <- rowSums(datExpr)
# Total <- sum(reads)

bwModuleCol <- data.frame(reads) %>% 
  as_tibble(rownames = "Name") %>% 
  right_join(bwModuleCol) 
  

# bwModuleCol %>%  group_by(module) %>% count(sort = T) %>% left_join(stats)

bwModuleCol %>% 
  group_by(module, biotype) %>% 
  summarise(reads = sum(reads), n = n()) %>% 
  arrange(desc(n)) -> bwModuleDF

bwModuleDF %>% 
  mutate(module = factor(module, levels = hclust$labels[hclust$order])) -> bwModuleDF

bwModuleDF %>% 
  group_by(module) %>%
  # ungroup() %>%
  mutate(cluster_frac = n / sum(n), 
    reads_frac = reads / sum(reads)) -> bwModuleDF

# bwModuleDF %>% tally(reads_frac)


df1 %>%
  mutate(star = ifelse(corPvalueStudent <.001, "***", 
    ifelse(corPvalueStudent <.01, "**",
      ifelse(corPvalueStudent <.05, "*", "")))) -> df1

df1 %>%
  mutate(facet = ifelse(name %in% c('HR11076', 'HR2476'), 'Low pH', 'Control')) %>%
  mutate(moduleTraitCor = round(moduleTraitCor, 2)) %>%
  mutate(star = ifelse(star != '', paste0(moduleTraitCor, '(', star,')'), moduleTraitCor)) %>%
  # mutate(star = ifelse(star != '', paste0(moduleTraitCor, '(', star,')'), '')) %>%
  ggplot(aes(y = module, x = name, fill = moduleTraitCor)) +
  geom_tile(color = 'white', size = 0.7, width = 1) +
  # geom_raster() +
  geom_text(aes(label = star),  vjust = 0.5, hjust = 0.5, size= 4, family =  "GillSans") +
  ggsci::scale_fill_gsea(name = "", reverse = T, na.value = "white") +
  # scale_fill_viridis_c(name = "Membership", na.value = "white") +
  ggh4x::scale_y_dendrogram(hclust = hclust) +
  labs(x = '', y = 'Module') +
  guides(fill = guide_colorbar(barwidth = unit(3.5, "in"),
    barheight = unit(0.1, "in"), label.position = "top",
    alignd = 0.5,
    ticks.colour = "black", ticks.linewidth = 0.5,
    frame.colour = "black", frame.linewidth = 0.5,
    label.theme = element_text(size = 10))) +
  theme_classic(base_size = 12, base_family = "GillSans") +
  theme(legend.position = "top",
    strip.background = element_rect(fill = 'white', color = 'white')) +
  facet_wrap(~ facet, scales = 'free_x') -> p1 

p1 <- p1 + theme(
  axis.line.x = element_blank(),
  axis.line.y = element_blank(),
  axis.text.y = element_text(hjust = 1.2),
  axis.ticks.length = unit(5, "pt"))

p1 <- p1 + theme(panel.spacing.x = unit(-0.5, "mm"))


bwModuleDF %>% 
  ggplot(aes(y = module, fill = biotype)) +
  # facet_wrap(biotype ~ ., nrow = 1, scales = "free_x") +
  scale_x_continuous("Frac. of reads", labels = scales::percent) +
  geom_col(aes(x = reads_frac), width = 0.95, position = position_stack(reverse = TRUE)) +
  # geom_col(aes(x = reads_frac), width = 0.95, fill = "grey")
  scale_fill_manual(name = '', values = c("#303960", "#647687", "#E7DFD5")) + # grey90
  theme_classic(base_size = 12, base_family = "GillSans") +
  theme(legend.position = "top",
    strip.background = element_rect(fill = 'white', color = 'white'),
    axis.title.y = element_blank(), 
    axis.text.y= element_blank(),
    axis.ticks.y =element_blank(), 
    axis.line.y = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.length = unit(5, "pt")) -> p2


library(patchwork)

# p1 + plot_spacer() + p2 + plot_layout(widths = c(5,-0.5, 10))  #& theme(plot.margin = 0)

ps <- p1 + p2 + plot_layout(widths = c(6, 5)) + labs(caption = '* corPvalueStudent < 0.05 ') 

# ps

ggsave(ps, filename = 'WGCNA.png', 
  path = path, width = 8, height = 4.5, device = png, dpi = 300)


bwModuleDF %>%
  ungroup() %>%
  mutate(x = n / sum(n), y = reads / sum(reads)) %>%
  ggplot(aes(x, y, color = biotype)) +
  # facet_wrap(~module) +
  theme_classic(base_size = 12, base_family = "GillSans") +
  scale_color_manual(name = '', values = c("#303960", "#647687", "grey90")) +
  geom_point()

bwModuleDF %>% view()

df1 %>% filter(star != "") %>% View()

