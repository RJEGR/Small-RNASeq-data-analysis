

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

path <- "~/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out//"

list.files(path = path, pattern = "txt")

count_f <- list.files(path = path, pattern = "Counts.txt", full.names = T)

res_f <- list.files(path = path, pattern = "Results.txt", full.names = T)


RESULTS <- read_tsv(res_f)

COUNTS <- read_tsv(count_f)

colNames <- COUNTS %>% select_if(is.double) %>% names()

names(COUNTS)[names(COUNTS) %in% colNames] <- gsub(".clean.newid.subset", "", colNames)

colNames <- COUNTS %>% select_if(is.double) %>% names()

# Sanity check

identical(RESULTS$Locus, COUNTS$Coords)
identical(RESULTS$Name, COUNTS$Name)

# cbind()

# PREVALENCE ----
nrow(raw_count <- COUNTS[colNames]) # 66,226

head(raw_count <- as.data.frame(raw_count, row.names = COUNTS$Name))

apply(raw_count, 1, function(x) sum(x > 0)) %>% table()

prevelancedf = apply(raw_count, 1, function(x) sum(x > 0))

mean_se = apply(raw_count, 1, function(x) mean_se(x)) %>% do.call(rbind, .)

data.frame(Prevalence = prevelancedf, 
  TotalAbundance = rowSums(raw_count),
  mean_se) %>% 
  as_tibble(rownames = "Name") %>%
  arrange(desc(TotalAbundance)) -> prevelancedf

prevelancedf <- prevelancedf %>% left_join(RESULTS %>% select(Name, MIRNA))

prevelancedf %>% 
  mutate(MIRNA = factor(MIRNA, levels = c("Y","N") )) %>%
  # filter(MIRNA == "Y") %>%
  ggplot(aes(TotalAbundance, Prevalence/12)) +
  geom_point(aes(color = MIRNA), alpha = 0.2) +
  # stat_density2d(aes(alpha = ..density..), geom = "tile", contour = FALSE) +
  scale_x_log10("Total Abundance (log10 scale)", labels = scales::comma) +  
  scale_y_continuous("Prevalence [Frac. Samples]", 
    labels = scales::percent_format(scale = 100)) +
  theme_classic(base_family = "GillSans") +
  theme(legend.position="top") + 
  facet_grid(~ MIRNA) 
  # scale_color_manual(values = c("red", "grey"))

prevelancedf %>% 
  mutate(MIRNA = factor(MIRNA, levels = c("Y","N") )) %>%
  # filter(MIRNA == "Y") %>%
  ggplot(aes(log10(TotalAbundance), color  = MIRNA, fill = MIRNA)) +
  # scale_x_continuous() +
  stat_ecdf() 
  # geom_density(alpha = 0.3)


prevelancedf %>% 
  arrange(Prevalence) %>%
  mutate(Prevalence = paste0(Prevalence, ' Samples')) %>%
  mutate(Prevalence = factor(Prevalence, levels = unique(Prevalence))) %>%
  mutate(TotalAbundance = log2(TotalAbundance+1)) %>% # edgeR::cpm(count) 
  ggplot(aes(TotalAbundance)) + geom_histogram() + 
  facet_wrap(~ Prevalence, scales = 'free_y') -> p1

dat_text <- prevelancedf %>% group_by(Prevalence) %>% tally() %>% 
  mutate(cumsum = cumsum(n)) %>% arrange(Prevalence) %>%
  mutate(Prevalence = paste0(Prevalence, ' Samples')) %>%
  mutate(Prevalence = factor(Prevalence, levels = unique(Prevalence)))

p1 <- p1 + geom_text(
  data    = dat_text, family = "GillSans",
  mapping = aes(x = -Inf, y = -Inf, label = paste0(n, " genes")),
  hjust   = -1, vjust   = -2, label.size = 0.2) + 
  theme_classic(base_size = 7, base_family = "GillSans") +
  labs(x = expression(~Log[2]~('TotalAbundance'~+1)), y = "")


# TEST DESEQ2 ----

mtd <- read_tsv(list.files(path = path, pattern = 'METADATA', full.names = T))

nrow(raw_count <- COUNTS[colNames])

# x) Choose contrast  ----

colData <- mtd

f_col <- "CONTRAST_A"

ref_group <- "Control"

names(colData)[names(colData) %in% f_col] <- "Design"

# 0) Remove single samples ----

sam <- table(colData$Design)

sam <- names(sam[sam > 1])

colData <- colData %>% drop_na(Design) %>% filter(Design %in% sam)

x <- colData %>% pull(Design, name = LIBRARY_ID)

# 1) According to metadata, filter samples ----

keep_cols <- names(raw_count) %in% colData$LIBRARY_ID

sum(keep_cols)

nrow(count <- raw_count[,keep_cols])

# dim(count <- raw_count[names(raw_count) %in% names(x)])

count <- as(count, "matrix")

rownames(count) <- COUNTS$Coords

# 2) Filter data by removing low-abundance genes ----

keep_rows <- rowSums(edgeR::cpm(count) > 1) >= 2

sum(keep_rows) # N transcripts

nrow(count <- count[keep_rows,])

dim(count <- round(count))


apply(count, 1, function(x) sum(x > 0)) %>% table()


# 3) Format metadata and count matrix ----

colData <- colData %>% arrange(match(LIBRARY_ID, names(count)))

sum(colnames(count) == colData$LIBRARY_ID) # sanity check


colData <- mutate_if(colData, is.character, as.factor)

# using relevel, just specifying the reference level:

colData <- colData %>% mutate(Design = relevel(Design, ref = ref_group))

library(DESeq2)
library(tidyverse)


table(colData$Design)

ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = count,
  colData = colData,
  design = ~ Design )

# 4) Run DESeq in the following functions order ----

dds <- estimateSizeFactors(ddsFullCountTable) 

# dds <- estimateDispersions(dds)
# dds <- nbinomWaldTest(dds) # use Wald test when contrasting a unique pair (ex. control vs cancer)

# Note: also run as a single batch 
# dds <- estimateSizeFactors(ddsFullCountTable)

dds <- DESeq(dds)

# 5) Test transformation for downstream analyses ----

# In order to test for differential expression, 
# we operate on raw counts and use discrete distributions ... 
# However for other downstream analyses – e.g. for visualization or clustering – it might be useful to work with transformed versions of the count data. Ex:

vst <- vst(dds) # wrapped version for varianceStabilizingTransformation

ntr <- DESeq2::normTransform(dds)

# 
# DESeq2::plotPCA(ntr, intgroup = "Design")
#
# DESeq2::plotPCA(vst, intgroup = "Design")

# BiocManager::install("vsn")

# VERIFICAR LA SD DE LOS DATOS TRANSFORMADOS, normTransform y vsd. Encontrar que transformacion disminuye la sd (linea roja) de los datos.

# ntr mantiene un sd debajo de 1 para la mayoria de los transcritos (eje x)

# vsn::meanSdPlot(assay(vst))
# vsn::meanSdPlot(assay(dds))
# vsn::meanSdPlot(assay(ntr)) 

# En este caso, ntr funciona mejor para estos datos


# 6 Generate results ----

get_res <- function(dds, contrast, alpha_cutoff = 0.1, lfc_cutoff = 0) {
  
  sA <- contrast[1]
  sB <- contrast[2]
  
  contrast <- as.character(DESeq2::design(dds))[2]
  
  keepA <- as.data.frame(colData(dds))[,contrast] == sA
  keepB <- as.data.frame(colData(dds))[,contrast] == sB
  
  contrast <- c(contrast, sA, sB)
  
  res = results(dds, contrast = contrast, alpha = alpha_cutoff, lfcThreshold = lfc_cutoff)
  
  
  baseMeanA <- rowMeans(DESeq2::counts(dds,normalized=TRUE)[,keepA])
  baseMeanB <- rowMeans(DESeq2::counts(dds,normalized=TRUE)[,keepB])
  
  res %>%
    as.data.frame(.) %>%
    cbind(baseMeanA, baseMeanB, .) %>%
    cbind(sampleA = sA, sampleB = sB, .) %>%
    as_tibble(rownames = "ids") %>%
    mutate(padj = ifelse(is.na(padj), 1, padj)) %>%
    mutate_at(vars(!matches("ids|sample|pvalue|padj")),
      round ,digits = 2)
}

contrast <- levels(colData(dds)$Design)

res <- get_res(dds, contrast) # # take a long

# write_rds(dds, file = paste0(path, '/CONTRAST_A_DDS.rds'))

# 7) multiple contrat from designs against control -----

DESeqDataSetFromInputs <- function(colData, count, f_col = ...) {
  
  names(colData)[names(colData) %in% f_col] <- "Design"
  
  # 0) Remove single samples:
  
  sam <- table(colData$Design)
  
  sam <- names(sam[sam > 1])
  
  colData <- colData %>% drop_na(Design) %>% filter(Design %in% sam)
  
  x <- colData %>% pull(Design, name = LIBRARY_ID)
  
  # 1) According to colData, filter samples ----
  
  ncol(count <- count[names(count) %in% names(x)])
  
  # 2) Filter data by removing low-abundance genes ----
  
  keep <- rowSums(edgeR::cpm(count) > 1) >= 2
  
  nrow(count <- count[keep,])
  
  count <- round(count)
  
  cat("Data to process: ", dim(count), "\n")
  
  # 3) Format metadata and count matrix ----
  
  colData <- colData %>% arrange(match(LIBRARY_ID, names(count)))
  
  sum(names(count) == colData$LIBRARY_ID) # sanity check
  
  colData <- mutate_if(colData, is.character, as.factor)
  
  # using relevel, just specifying the reference level:
  
  colData <- colData %>% mutate(Design = relevel(Design, ref = "Control"))
  
  library(DESeq2)
  library(tidyverse)
  
  ddsFullCountTable <- DESeqDataSetFromMatrix(
    countData = count,
    colData = colData,
    design = ~ Design )
  
  # 4) estimateSizeF for DESeq  ----
  
  dds <- estimateSizeFactors(ddsFullCountTable) 
  
  out <- DESeq(dds)
  
  write_rds(out, file = paste0(path, f_col, '_DDS.rds'))
  
  return(out)
}

count <- as(raw_count, "matrix")

rownames(count) <- COUNTS$Coords

dds_A <- DESeqDataSetFromInputs(mtd, raw_count, f_col = "CONTRAST_A") # 62,817 transcripts and 56 samples
dds_B <- DESeqDataSetFromInputs(mtd, raw_count, f_col = "CONTRAST_B") # 21,006 transcripts and 56 samples
dds_C <- DESeqDataSetFromInputs(mtd, raw_count, f_col = "CONTRAST_D") # 44,167 transcripts and 56 samples
dds_D <- DESeqDataSetFromInputs(mtd, raw_count, f_col = "CONTRAST_C") # 110533 transcripts and 56 samples


# Generate results 
get_res <- function(dds, contrast, alpha_cutoff = 0.1, weighted_p = FALSE) {
  
  require(DESeq2)
  require(tidyverse)
  require(IHW)
  
  sA <- contrast[1]
  sB <- contrast[2]
  
  contrast <- as.character(DESeq2::design(dds))[2]
  
  keepA <- as.data.frame(colData(dds))[,contrast] == sA
  keepB <- as.data.frame(colData(dds))[,contrast] == sB
  
  contrast <- c(contrast, sA, sB)
  
  res = results(dds, contrast, alpha = 0.1, lfcThreshold = 0)
  
  
  baseMeanA <- rowMeans(DESeq2::counts(dds,normalized=TRUE)[,keepA])
  baseMeanB <- rowMeans(DESeq2::counts(dds,normalized=TRUE)[,keepB])
  
  
  res <- res %>%
    as.data.frame(.) %>%
    cbind(baseMeanA, baseMeanB, .) %>%
    cbind(sampleA = sA, sampleB = sB, .) %>%
    mutate(padj = ifelse(is.na(padj), 1, padj)) %>%
    mutate_at(vars(!matches("transcript_id|sample|pvalue|padj")),
      round ,digits = 2)
  
  if(weighted_p) {
    
    #  take the set of p-values and then calculate the adjusted p-values using a independent Hypothesis Weighting
    
    ihwRes <- ihw(pvalue ~ baseMean,  data = res, alpha = 0.1) # toma menos de 1 minuto.
    
    weighted_pvalue <- as.data.frame(ihwRes)$weighted_pvalue
    
    
    res <- res %>% cbind(., weighted_pvalue = weighted_pvalue) %>% 
      as_tibble(rownames = "transcript_id")
    
  } else
  
  return(res)
}

dds2res <- function(dds_path_f, ...) {
  
  dds <- read_rds(dds_path_f)
  
  contrast <- levels(colData(dds)$Design)
  
  out <- list()
  
  for(i in 2:length(contrast)) {
    j <- i
    cat('\nContrast', contrast[1],' and ', contrast[j],'\n')
    res <- get_res(dds, contrast[c(1,j)], alpha_cutoff = 0.1)
    out[[j]] <- res
    
  }
  
  do.call(rbind, out) -> res
  
  return(res)
  
}


rds_f <- list.files(path, pattern = 'CONTRAST_',  full.names = TRUE)

out <- lapply(rds_f, dds2res)

names(out) <- gsub("_DDS.rds", "", basename(rds_f))

res <- bind_rows(out, .id = 'Contrast')

do.call(rbind, out) -> res

