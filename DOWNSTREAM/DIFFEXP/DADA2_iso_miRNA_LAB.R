# OMIT
# due to isoMIRNAs lest try dada2 to find single sequence variants
# revisar codigo filterAndTriming_learnError_test.R AND coi_dada2_mock.R (OMEGA_ALFA function)

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

# library(tidyverse)
library(dada2)
# problemas con R CMD INSTALL Matrix_1.5-3.tar.gz
# packageVersion("Matrix")
# library(Matrix)

packageVersion("dada2")

# matrix needs https://github.com/fxcoudert/gfortran-for-macOS/releases


data_path <- "~/Documents/MIRNA_HALIOTIS/RAW_DATA/usftp21.novogene.com/01.RawData/"


cat("\nInput files available in dir:", data_path, "\n")

subdir <- list.files(path = data_path, pattern = 'HR*', full.names = T)

list.files(subdir, pattern = "fq.gz")

date <- format(Sys.time(), "%Y%m%d%H%M%S")

run <- paste0("run","_", date)

# ================
# Outputs in `pwd`:
# ================

out_dir = "~/Documents/MIRNA_HALIOTIS/"
out_path <- file.path(out_dir, "isomiRNA_lab", run)
system(command = paste0("mkdir -p ", out_path), intern = F)


fnFs <- sort(list.files(path = subdir, pattern = 'fq.gz', full.names = T))

# subset to test 
fnFs <- fnFs[1]

cat("\nForward read files to process:\n")

fnFs


# seqs1 <- Biostrings::readDNAStringSet(fnFs[1], format="fastq")

# plotQualityProfile(fnFs[1])

source('~/Documents/GitHub/metagenomics/plotQualityProfile.R')

# plotQP(fnFs)

## Filtrar y recortar ----

# Crear directorio y nombres de archivo
sam_name <- sapply(strsplit(basename(fnFs), "[.]"), '[', 1)

filtFs <- file.path(out_path, paste0(sam_name, "_F_filt.fastq.gz"))

# Ahora si Filtrar y recortar ----

# Parameters FilterAndTrim

# adapter_3p <- "AGATCGGAAGAGCACACGTCT"

out <- filterAndTrim(fwd = fnFs, filt = filtFs,
  truncLen = 0, 
  # truncQ = 0, # (default - no truncation) 
  maxEE= Inf, 
  minLen = 15,
  trimLeft = 0, trimRight = 20,
  maxN = 5, # 10 % from 50SE sequencing
  rm.phix = FALSE, compress = T, multithread = F)


pct_trim <- 1 - out[,2]/out[,1]

pct_trim

# Dereplication ----
derepFs <- derepFastq(filtFs, verbose = T)

# Encountered 1,116,758 unique sequences from 19,164,919 total sequences read.

# Learn Error ----
MAX_CONSIST <- 20

threads <- NULL

errF <- learn_Err(filtFs)

source('~/Documents/GitHub/BM/filteringAndTriming_learnError_test.R')

plotErrors(errF, nominalQ=TRUE)

# dada algorithm
OMEGA_A <- 1e-40

dadaFs <- dada2(derepFs, OMEGA_A, errF)

propagateCol <- names(dadaFs)

nrow(mergers <- mergePairs(dadaFs, derepFs, 
  minOverlap = 12, maxMismatch = 1,
  verbose=TRUE, returnRejects = TRUE, 
  propagateCol=propagateCol))



dim(seqtab <- makeSequenceTable(mergers))

table(nchar(getSequences(seqtab)))

