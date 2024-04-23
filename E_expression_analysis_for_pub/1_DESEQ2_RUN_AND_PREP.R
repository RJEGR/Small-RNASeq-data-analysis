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

COUNTS <- COUNTS %>% filter(Name %in% query.ids)

# RE-GROUP MICRORNAS IF MAJORRNA SEQUENCE ARE IDENTICAL
# THE OUTPUT WILL BE A UPGRADED DATASET FOR D.E AND CO-EXPRESSION ANALYSIS FOR PUBLICATION

# DR.RICARDO GOMEZ-REYES

which_cols <- COUNTS %>% select_if(is.double) %>% names()

COUNTS <- DB %>% filter(SRNAtype == "miR") %>% 
  select(Name, MajorRNA) %>% 
  full_join(COUNTS, by = "Name") %>%
  group_by(MajorRNA) %>%
  summarise_at(vars(all_of(which_cols)), sum) 

# THEN:

rowNames <- COUNTS$MajorRNA

COUNTS <- COUNTS %>% select_if(is.double) %>% as(., "matrix")

rownames(COUNTS) <- rowNames

.COUNTS <- COUNTS

path_out <- "~/Documents/MIRNA_HALIOTIS/MIRNA_PUB_2024"

write_rds(.COUNTS, file.path(path_out, "IDENTICAL_SEQUENCES_MERGED_COUNT.rds"))

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

RES <- RES %>% dplyr::rename("MajorRNA" = "Name")

# deal w/ duplication ids due to random seed sequences and 5p/3p label

dup_majorRNA <- read_tsv(paste0(path, "SRNA2MIRGENEDB.tsv")) %>% 
  mutate(MirGeneDB_ID = paste0(Family,"-", arm)) %>%
  distinct(MajorRNA, Family, MirGeneDB_ID) %>% 
  dplyr::count(MajorRNA, sort = T) %>% 
  filter(n > 1) %>% pull(MajorRNA)


RES <- read_tsv(paste0(path, "SRNA2MIRGENEDB.tsv")) %>% 
  mutate(MirGeneDB_ID = paste0(Family,"-", arm)) %>%
  mutate(MirGeneDB_ID = ifelse(MajorRNA %in% dup_majorRNA, Family, MirGeneDB_ID)) %>%
  distinct(MajorRNA, MirGeneDB_ID) %>%
  right_join(RES)

RES <- RES %>%
  group_by(CONTRAST) %>%
  arrange(desc(MajorRNA), .by_group = T) %>%
  mutate(Name = paste0("Cluster_", row_number(MajorRNA)))



view(RES)

# RECODE 7 MIRS FAMILY FROM MIRBASE/MOLLUSC DB (FROM 4_MIR_PAIRWISEALIGNMENT.R)

recode_mir <- c(
  `Cluster_39774` = "MIR-242", 
  `Cluster_17642` = "MIR-1989b-5p", 
  `Cluster_38139` = "MIR-745", 
  `Cluster_45860` = "MIR-10492a-5p", 
  `Cluster_47716` = "MIR-277b-3p", 
  `Cluster_37147` = "MIR-184-3p", 
  `Cluster_27861` = "MIR-252-5p")


recode_majorRNA <- DB %>% filter(Name %in% names(recode_mir)) %>% pull(MajorRNA, name = Name)

recode_mir <- structure(recode_mir, names = recode_majorRNA)

RES <- RES %>% 
  mutate(MirGeneDB_ID = ifelse(is.na(MirGeneDB_ID), MajorRNA, MirGeneDB_ID)) %>%
  dplyr::mutate(MirGeneDB_ID = dplyr::recode(MirGeneDB_ID, !!!recode_mir))

RES <- RES %>% mutate(Name = ifelse(grepl("^MIR", MirGeneDB_ID), MirGeneDB_ID, Name))

write_tsv(RES, file = file.path(path_out, "SEQUENCES_MERGED_DESEQ_RES.tsv"))

# CREATE A DDS OBJ.

ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = .COUNTS,
  colData = .colData,
  design = ~ colName )

dds <- estimateSizeFactors(ddsFullCountTable) 

write_rds(dds,file = file.path(path_out, "SEQUENCES_MERGED_DDS_DESEQ2.rds"))

