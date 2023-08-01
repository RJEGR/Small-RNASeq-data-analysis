# THIS IS THE UPGRADE VERSION FOR D.E ANALYSIS
# 1) KEEP ONLY MIRS USING DB AS GUIDE (OPTIONAL)
# 2) RUN DESEQ TEST USING MULTIPLE CONTRAST COMPARISON:
# 3) BIND DE OUTPUT TO MIRGENDB AND SRNA_REGULATORY DB
# 4) SAVE AS RDATA FILE 
# 
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

write_rds(.COUNTS, paste0(path, "COUNT.rds"))

# rowSums(COUNTS) %>% as_tibble(rownames = "Name")

# DB %>% select(Name, Locus_type, MajorRNA, SRNAtype, Reads) %>% filter(SRNAtype == "miR")

# 1) Filter data by removing low-abundance genes ----

by_count <- 1; by_freq <- 2

keep <- rowSums(COUNTS > by_count) >= by_freq

sum(keep) # N transcripts

nrow(COUNTS <- COUNTS[keep,])

COUNTS <- round(COUNTS)


# In order to test for differential expression, we operate on raw counts and use discrete distributions ...

# 2) Run Multiple Contrast Comparison =====

library(DESeq2)

CONTRAST <- .colData %>% dplyr::select(starts_with("CONTRAST")) %>% names()

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

do.call(rbind, out) -> RES

#

# wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/MIRGENEDB_20230314/"

RES <- read_tsv(paste0(path, "SRNA2MIRGENEDB.tsv")) %>% 
  distinct(Name, Family) %>%
  right_join(RES)

RES <- RES %>% mutate(Family = ifelse(is.na(Family), Name, Family))

write_tsv(RES, file = paste0(path, "DESEQ_RES.tsv"))


RES %>% 
  filter(abs(log2FoldChange) > 0 & padj < 0.05) %>%
  mutate(REGULATED = sign(log2FoldChange)) %>%
  group_by(REGULATED, Name) %>%
  summarise(
    across(CONTRAST, .fns = paste_go), 
    n = n(),
    .groups = "drop_last")

