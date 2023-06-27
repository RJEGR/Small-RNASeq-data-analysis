
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


# detect max power ----

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
  networkType = "unsigned")

saveRDS(sft, file = paste0(path, 'SoftThreshold_',cor_method, '.rds'))

# Construct a gene co-expression matrix and generate modules ----
# 

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

# Blockwise construction ----

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

# START: 

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
