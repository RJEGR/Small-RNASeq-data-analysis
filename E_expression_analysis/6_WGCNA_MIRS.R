# WGCNA FOR ONLY MIRS
# 

# wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/FUNCTIONAL_MIR_ANNOT/"
# read_tsv(paste0(wd, "DESEQ2SRNATARGET.tsv"))

library(WGCNA)
library(flashClust)
library(tidyverse)


rm(list = ls());

if(!is.null(dev.list())) dev.off()

path <- "~/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out/"

head(COUNT <- read_rds(paste0(path, "COUNT.rds"))) # ONLY MIRS

DB <- read_tsv(paste0(path, "DESEQ_RES.tsv")) %>% distinct(Name, Family)

# MIRGENEDB <- read_tsv(paste0(path, "SRNA2MIRGENEDB.tsv")) %>% distinct(Name, Family, arm)
# MIRGENEDB <- MIRGENEDB %>% mutate(Fam = Family) %>% tidyr::unite("MIR", c("Fam", "arm"), sep = "-") 
# DB <- DB %>% left_join(MIRGENEDB)


# COLLAPSE BY FAM

DB %>% dplyr::count(Name, sort = T)

COUNT <- COUNT %>% as_tibble(rownames = "Name") %>%
  right_join(DB %>% distinct(Name, Family), by = "Name") %>%
  group_by(Family) %>%
  summarise_at(vars(colnames(COUNTS)), sum) 

which_cols <- COUNT %>% select_if(is.double) %>% colnames()

datExpr <- COUNT %>% select(all_of(which_cols))

datExpr <- as(datExpr, "matrix")

rownames(datExpr) <- COUNT$Family

# THEN,

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

# When you have a lot of genes use the following code
k=softConnectivity(datE=datExpr,power=6)
# Plot a histogram of k and a scale free topology plot
sizeGrWindow(10,5)
par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k, main="Check scale free topology\n")

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

# sft <- read_rds(paste0(path, 'SoftThreshold_',cor_method, '.rds'))

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

psave

