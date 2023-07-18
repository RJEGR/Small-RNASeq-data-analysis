# Differential Expression Analysis

 ([Example](https://www.nature.com/articles/nature25027)) The resulting annotated small-RNA loci (Supplementary Data 1) were analysed for differential expression (interface versus parasite stem) using DESeq236, with a log2 fold threshold of 1, alternative hypothesis of ‘greaterAbs’, and alpha of 0.05. P values were adjusted for multiple testing using the Benjamini–Hochberg procedure, and loci with an adjusted P value of ≤0.05 (equivalent to an FDR of ≤0.05) were denoted upregulated in interfaces relative to parasite stem. Among the upregulated loci, those annotated by ShortStack as miRNAs deriving from the C. campestris genome which produced either a 21- or 22-nucleotide mature miRNA (Supplementary Data 2) were retained and further analysed. The predicted secondary structures and observed small-RNA-seq read coverage was visualized (Supplementary Data 3, 4) using strucVis (version 0.3; https://github.com/MikeAxtell/strucVis).
 
```R
library(tidyverse)

path <- '~/Documents/MIRNA_HALIOTIS/CLEAN_INPUT/'

mtd <- read_tsv(list.files(path = path, pattern = 'METADATA', full.names = T))

fileName <- list.files(path = path, pattern = 'miRNAs_expressed_all_samples', full.names = T)

count <- read_tsv(file = fileName)

```
