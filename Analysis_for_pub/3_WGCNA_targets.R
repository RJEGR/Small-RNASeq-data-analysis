# WGCNA FOR Target
# this version use the majorRNA sequences summarization
# 


library(WGCNA)
library(flashClust)
library(tidyverse)


library(tidyverse)

rm(list = ls());

if(!is.null(dev.list())) dev.off()

wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/FUNCTIONAL_MIR_ANNOT/"


print(FUNCTIONAL_DB <- read_rds(paste0(wd,"SRNA_FUNCTION_PREDICTED_LONG_EXPRESSED.rds")))

is_na <- function(x) ifelse(is.na(x), 0, x)

keep_expressed <- FUNCTIONAL_DB %>% 
  dplyr::select(starts_with("SRR")) %>% 
  mutate(across(where(is.double), ~is_na(.))) %>% 
  rowSums()

# FUNCTIONAL_DB <- FUNCTIONAL_DB[keep_expressed > 1,]

datExpr <- FUNCTIONAL_DB[keep_expressed > 1,] %>% 
  dplyr::select(starts_with("SRR")) %>%
  as("matrix")

rownames(datExpr) <- FUNCTIONAL_DB$gene_id
# Nota: usar transformacion ablanda los datos
# datExpr <- log2(datExpr+1)
# z_scores <- function(x) {(x-mean(x))/sd(x)}
# datExpr <- apply(datExpr, 1, z_scores)

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

# When you have a lot of genes use the following code
k=softConnectivity(datE=datExpr,power=6)
# Plot a histogram of k and a scale free topology plot

scaleFreePlot(k, main="Check scale free topology\n")

max_power <- 20

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

# sft <- read_rds(paste0(path, 'SoftThreshold_',cor_method, '.rds'))

soft_values <- abs(sign(sft$fitIndices[,3])*sft$fitIndices[,2])

soft_values <- round(soft_values, digits = 2)

hist(soft_values)

# power_pct <- quantile(soft_values, probs = 0.95)

power_pct <- soft_values[which.max(soft_values)]

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
  facet_grid(~name, scales = 'free_x', switch = "x") +
  geom_text(aes(label = Power), size = 2) +
  geom_abline(slope = 0, intercept = softPower, linetype="dashed", alpha=0.5) +
  # geom_vline(xintercept = min(meanK), linetype="dashed", alpha=0.5) +
  labs(y = 'Soft Threshold (power)', x = '', 
    caption = caption) +
  # scale_x_continuous(position = "top") +
  theme_light(base_family = "GillSans",base_size = 16) -> psave

psave

# Build a adjacency "correlation" matrix

enableWGCNAThreads()

# specify network type

adjacency <- adjacency(datExpr, 
  power = softPower, 
  type = "unsigned")


TOM <- TOMsimilarity(adjacency, TOMType = "unsigned") # specify network type


dissTOM = 1 - TOM

heatmap(TOM)

# Generate Modules ----
# Generate a clustered gene tree

geneTree = flashClust(as.dist(dissTOM), method="average")

plot(geneTree, xlab="", sub="", 
  main= "Gene Clustering on TOM-based dissimilarity", labels= FALSE, hang=0.04)

#This sets the minimum number of genes to cluster into a module

#   Modules of co-expressed miRNAs were then determined by using the dynamic tree cut procedure with a minimum module size cut-off of 5. This cut-off was chosen considering the miRNA transcriptome’s small size and the fact that a single miRNA can target multiple RNA transcripts. With this procedure, 18 different modules of EigenmiRNA were identified (https://doi.org/10.1038/s41598-020-)

minClusterSize <- softPower

dynamicMods <- cutreeDynamic(dendro= geneTree, 
  distM = dissTOM,
  method = "hybrid",
  deepSplit = 2, 
  cutHeight = 0.8,
  pamRespectsDendro = FALSE,
  minClusterSize = minClusterSize)

dynamicColors = labels2colors(dynamicMods)
names(dynamicColors) <- colnames(datExpr)
MEList = moduleEigengenes(datExpr, 
  colors= dynamicColors,
  softPower = softPower)

MEs = MEList$eigengenes

MEDiss= 1 - cor(MEs)

METree = flashClust(as.dist(MEDiss), method= "average")

# Set a threhold for merging modules. In this example we are not merging so MEDissThres=0.0

MEDissThres = 0

title <- paste0("Dynamic Tree Cut", "Merged dynamic\n(cutHeight=", MEDissThres,")")

merge = mergeCloseModules(datExpr, dynamicColors, cutHeight= MEDissThres, verbose =3)

mergedColors = merge$colors

mergedMEs = merge$newMEs

#plot dendrogram with module colors below it
plotDendroAndColors(geneTree, dynamicColors, c("Modules"), dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang=0.05)

plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), title, dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang=0.05)

moduleColors <- mergedColors

colorOrder <- c("grey", standardColors(50))

moduleLabels <- match(moduleColors, colorOrder) - 1


MEs <- mergedMEs


# how modules where obtained:
nm <- table(moduleColors)


cat('Number of mudules obtained\n :', length(nm))

print(nm)

WGCNA <- data.frame(gene_id = names(moduleColors), WGCNA_target = moduleColors)

WGCNA <- as_tibble(WGCNA)

FUNCTIONAL_DB <- WGCNA %>% right_join(FUNCTIONAL_DB) %>% distinct()

write_rds(FUNCTIONAL_DB, file = paste0(wd, "SRNA_FUNCTION_PREDICTED_LONG_EXPRESSED.rds"))
