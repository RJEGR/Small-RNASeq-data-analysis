# THIS IS THE UPGRADE VERSION FOR D.E ANALYSIS
# 1) KEEP ONLY MIRS USING DB AS GUIDE
# 2) RUN DESEQ TEST USING MULTIPLE CONTRAST COMPARISON:

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

path <- "~/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out/"

count_f <- list.files(path = path, pattern = "Counts.txt", full.names = T)

DB_f <- list.files(path = path, pattern = "RNA_LOCATION_DB.tsv", full.names = T)

MTD_f <- list.files(path = path, pattern = "METADATA.tsv", full.names = T)

# wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/FUNCTIONAL_MIR_ANNOT/"

# SRNA2GO <- read_tsv(paste0(wd, "SRNA2GO.tsv")) # "SRNA_REGULATORY_FUNCTION_DB.tsv"

.colData <- read_tsv(MTD_f)

DB <- read_tsv(DB_f)

str(query.ids <- DB %>% filter(SRNAtype == "miR") %>% distinct(Name) %>% pull())

COUNTS <- read_tsv(count_f)

colNames <- gsub(".clean.newid.subset", "", names(COUNTS))

colnames(COUNTS) <- colNames

# OPTIONAL:

COUNTS <- COUNTS %>% filter(Name %in% query.ids)

# THEN:

rowNames <- COUNTS$Name

COUNTS <- COUNTS %>% select_if(is.double) %>% as(., "matrix")

rownames(COUNTS) <- rowNames

.COUNTS <- COUNTS

rowSums(COUNTS) %>% as_tibble(rownames = "Name")

DB %>% select(Name, Locus_type, MajorRNA, SRNAtype, Reads) %>% filter(SRNAtype == "miR")

# 1) Filter data by removing low-abundance genes ----

by_count <- 1; by_freq <- 2

keep <- rowSums(COUNTS > by_count) >= by_freq

sum(keep) # N transcripts

nrow(COUNTS <- COUNTS[keep,])

COUNTS <- round(COUNTS)

# 2) Format metadata and count matrix ----

colData <- .colData %>% arrange(match(LIBRARY_ID, colnames(COUNTS)))

any(colnames(COUNTS) == colData$LIBRARY_ID) # sanity check

colData <- mutate_if(colData, is.character, as.factor)

# Using colData to create the experimental design

f_col <- "CONTRAST_A"

names(colData)[names(colData) %in% f_col] <- "Design"

# DROP NA COLS

colData <- colData %>% drop_na(Design)

keep <- colnames(COUNTS) %in% colData$LIBRARY_ID 

head(COUNTS <- COUNTS[,keep])

# Si es posible, especificar  el factor "control" como nivel de referencia:

colData <- colData %>% mutate(Design = relevel(Design, ref = "Control"))

# DESEQ2 =====
library(DESeq2)

ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = COUNTS,
  colData = colData,
  design = ~ Design )

# In order to test for differential expression, we operate on raw counts and use discrete distributions ...

# 3) Run DESeq in the following functions order ----
# Negative Binomial GLM

# # For the analysis we need to estimate the effective library size to normalize for.
# Imagine, a gene has the same number of counts in two samples. But the library size 
# (total number of reads) was twice as high in the first sample. Then we would conclude that the gene was higher expressed in the second sample. 
# You can estimate the size factors from the count data using the function estimateSizeFactors():

dds <- estimateSizeFactors(ddsFullCountTable) 

# sizeFactors(dds)

# Next we need to estimate for each condition a function that allows to predict the 
# dispersion (= variance across samples). The core assumption of this method is that the  # mean is a good predictor of the variance, i.e., that genes with a similar expression level also have similar variance across the samples:

dds <- estimateDispersions(dds)

# boxplot(log2(counts(ddsFullCountTable)+0.5))
# boxplot(log2(counts(dds, normalized=TRUE)+0.5))

# Now we can test for differential expression of the genes between the two groups 
# by calling the function nbinomWaldTest(). We provide the dds and the names of 
# the groups of our samples to this function. This might take a few minutes:

dds <- nbinomWaldTest(dds)

# (Above) use Wald test when contrasting a unique pair (ex. group1 vs group2)

res <- results(dds)

summary(res)

# OUTPUT DF ====

contrast <- levels(colData(dds)$Design)

res <- get_res(dds, contrast)

res.p <- prep_DE_data(res, alpha = 0.05, lfcThreshold = 2) %>% drop_na(cc)

res.p

# Note: also run as a single batch 

# dds <- DESeq(ddsFullCountTable);res <- results(dds)

# CONCACT

CONTRAST <- .colData %>% select(starts_with("CONTRAST")) %>% names()

run_contrast_DE <- function(COUNTS, colData, CONTRAST = NULL, ref = NULL) {
  
  # CONTRAST: Column in colData with character vector of Design
  
  names(colData)[1] <- "LIBRARY_ID"
  
  names(colData)[names(colData) %in% CONTRAST] <- "Design"
  
  colData <- colData %>% drop_na(Design)
  
  # any(colnames(COUNTS) == colData$LIBRARY_ID) # sanity check
  
  colData <- mutate_if(colData, is.character, as.factor)
  
  keep <- colnames(COUNTS) %in% colData$LIBRARY_ID 
  
  COUNTS <- COUNTS[,keep]
  
  colData <- colData %>% mutate(Design = relevel(Design, ref = "Control"))
  
  require(DESeq2)
  
  ddsFullCountTable <- DESeqDataSetFromMatrix(
    countData = COUNTS,
    colData = colData,
    design = ~ Design )
  
  dds <- estimateSizeFactors(ddsFullCountTable) 
  
  dds <- estimateDispersions(dds)
  
  dds <- nbinomWaldTest(dds)
  
  # return(dds)
  
  contrast <- levels(colData(dds)$Design)
  
  res <- get_res(dds, contrast)
  
  return(res)
}

# run_contrast_DE(.COUNTS, .colData, CONTRAST = CONTRAST[1])

any(colnames(.COUNTS) == .colData$LIBRARY_ID)

out <- list()

for (j in 1:length(CONTRAST)) {
  
  i <- j
  
  cat("\nRunning ",CONTRAST[i], "\n")
  

  out [[i]] <- run_contrast_DE(.COUNTS, .colData, CONTRAST = CONTRAST[i]) %>% mutate(CONTRAST = CONTRAST[i])
}

do.call(rbind, out) -> allRes

# POSITIVE LFC == UP EXPRESSED IN EXPERIMENTAL
# NEGATIVE LFC == UP EXPRESSED IN CONTROL

# recode_to <- paste(LETTERS[1:4], ")", sep = "")

recode_to <- c("24 HPF: Ctrl|pH", "110 HPF: Ctrl|pH", "Ctrl: 24|110 HPF", "pH: 24|110")
recode_to <- structure(recode_to, names = CONTRAST)

allRes <- allRes %>% 
  dplyr::mutate(CONTRAST = dplyr::recode_factor(CONTRAST, !!!recode_to))

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
  
allRes %>%
