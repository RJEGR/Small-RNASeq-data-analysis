
# Using a biochemically-based miRNA target prediction


rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

# BiocManager::install("scanMiR")

library(scanMiR)
library(tidyverse)

# 0) Input queries ----
# Using a miRNA sequences predicted by shortstacks

wd <- "~/Documents/MIRNA_HALIOTIS/"

query_path <- paste0(wd, "/SHORTSTACKS/ShortStack_20230315_out/")

ref_path <- paste0(wd, "/RNA_SEQ_ANALYSIS/ASSEMBLY/STRINGTIE/REFBASED_MODE/ANNOTATION/")

res_f <- list.files(path = query_path, pattern = "Results.txt", full.names = T)

RESULTS <- read_tsv(res_f)

RESULTS %>% dplyr::count(MIRNA)

#  pull(MajorRNA, name = Name) 

MajorRNA <- RESULTS %>% filter(MIRNA == "Y") %>% pull(MajorRNA) 

MajorRNA <- sample(MajorRNA, 10)


# 1) Input targets ----

# 1.1) larvae rna-seq assembly

# f <- " /Users/cigom/Documents/MIRNA_HALIOTIS/RNA_SEQ_ANALYSIS/ASSEMBLY/STRINGTIE/REFBASED_MODE/transcript.fa"

# 1.2) Genomic five_prime_utr and three_prime_utr flanked fasta format

# f <- '/Users/cigom/Documents/MIRNA_HALIOTIS/ENSEMBLE/utr_rmdup.fa'

# `seqkit rmdup -s -P -D duplicated_detail.txt -d duplicates.fasta -o utr_rmdup.fa`

# cat transcripts.fa.transdecoder.cds | seqkit rmdup -s -P -D duplicated_detail.txt -d transcripts.fa.transdecoder.dups.cds -o transcripts.fa.transdecoder.rmdup.cds


# 108812 CDS
# 108812 exon
# 87680 five_prime_UTR
# 108812 gene
# 108812 mRNA
# 100219 three_prime_UTR


gff_f <- list.files(path = ref_path, pattern = "transcripts.fa.transdecoder.gff3", full.names = T)

# nrow(ape::read.gff(gff_f))

gff <- read_tsv(gff_f, comment = "#", col_names = F)

gffNames <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")

names(gff) <- gffNames

seqname_id <- gff %>% 
  filter(grepl("prime_UTR",feature)) %>% 
  arrange(seqname) %>% distinct(seqname) %>% pull()

f <- list.files(path = ref_path, pattern = "transcripts.fa.transdecoder.rmdup.cds", full.names = T)


seqs <- Biostrings::readDNAStringSet(f)

str(keep <- names(seqs))

keep <- sapply(strsplit(keep, " "), `[`, 1)

str(keep <- unique(keep))

sum(keep <- keep %in% seqname_id)

subseqs <- sample(seqs, 10)

# 2) run scan ----
# matches2 <- findSeedMatches(seqs, SampleKdModel, verbose = FALSE)

matchesdf <- findSeedMatches(subseqs, MajorRNA, verbose = F, keepMatchSeq = T)


MajorRNA <- gsub("U", "T", MajorRNA)

hist(matchesdf$p3.score)

table(matchesdf$type)

n <- 2

sum(keep <- as.character(matchesdf$miRNA) %in% MajorRNA[n])

m <- matchesdf[keep,]

m 

seqnames <- subseqs[names(subseqs) %in% as.character(m@seqnames)[5]]

viewTargetAlignment(m[5,], MajorRNA[n], seqs = seqnames)

?scanMiR::plotKdModel()

matchesdf

matchesdf %>% view()

matchesdf <- as.data.frame(matchesdf)

mrna_id <- as.character(matchesdf$seqnames)

# mrna_id <- sapply(strsplit(mrna_id, " "), `[`, 1)

mrna_seq <- seqs[names(seqs) %in% mrna_id]
mrna_seq <- as.character(mrna_seq)
names(mrna_seq) <- NULL

mirna_seq <- as.character(matchesdf$miRNA)

label <- matchesdf$p3.score
# mirna_id <- ...

head(data.frame(mirna_seq, mrna_id, mrna_seq, label))
# Output # mirna_id        mirna_seq       mrna_id mrna_seq        label   split




# run scan
# load a sample transcript
data("SampleTranscript")

matches <- findSeedMatches(SampleTranscript, miRNA, verbose = FALSE)

matches

viewTargetAlignment(matches[1], miRNA, SampleTranscript)

# Using a KdModel
# KdModel collections corresponding to all human, mouse and rat mirbase miRNAs can be obtained through the scanMiRData package.

# load sample KdModel
data("SampleKdModel")

# run scan

matches <- findSeedMatches(SampleTranscript, SampleKdModel, verbose = FALSE)
matches


# by ORF length

library(Biostrings)

# generate set of random sequences
seqs <- DNAStringSet(getRandomSeq(length = 1000, n = 10))

# add vector of ORF lengths
mcols(seqs)$ORF.length <- sample(500:800, length(seqs))

# run scan
# matches2 <- findSeedMatches(seqs, SampleKdModel, verbose = FALSE)
matches2 <- findSeedMatches(seqs, miRNA, verbose = FALSE)
head(matches2)

viewTargetAlignment(matches[1], SampleKdModel, SampleTranscript)

# Creating a KdModel object ---- 
# not working
# The scanMiRData package contains KdModel collections corresponding to all human, mouse and rat mirbase miRNAs.

scanMiR::getkd

?getKdModel()

dim(kd <- dummyKdData())
str(kd)


mod3 <- getKdModel(kd=kd, mirseq="TTAATGCTAATCGTGATAGGGGTT", name = "my-miRNA")


summary(mod3)


# devtools::install_github("ETHZ-INS/scanMiRData")

# library(scanMiRData)

# getKdModels(species = c("hsa", "mmu", "rno"), categories = NULL)

w_categories <- c("Low-confidence", "Poorly conserved", "Conserved across mammals", "Conserved across vertebrates")

mods <- getKdModels("mmu", categories = w_categories)

# summary(mods)

dim(kd <- dummyKdData(KdModelList(mods, makeUnique = T)[[1]]))
str(kd)
# dim(kd <- dummyKdData())

mirseq <- "TTAATGCTAATCGTGATAGGGGTT" #"CUAUACAAUCUACUGUCUUUCC"
mirseq <- gsub("U", "T", mirseq)

getKdModel2(kd=kd, mirseq=mirseq, name = "my-miRNA")


getKdModel2 <- function(kd, mirseq = NULL, name = NULL, conservation = NA_integer_, 
  ...) {
    if (is.character(kd) && length(kd) == 1) {
      if (is.null(name)) 
        name <- gsub("\\.txt$|\\.csv$", "", gsub("_kds", 
          "", basename(kd)))
      kd <- read.delim(kd, header = TRUE, stringsAsFactors = FALSE)[, 
        c(1, 2, 4)]
    }
    if (is.null(mirseq) && !is.null(kd$mirseq)) 
      mirseq <- as.character(kd$mirseq[1])
    if (is.null(name) && !is.null(kd$mir)) 
      name <- as.character(kd$mir[1])
    stopifnot(!is.null(name) && !is.null(mirseq))
    if (!("X12mer" %in% colnames(kd)) && "12mer" %in% colnames(kd)) 
      colnames(kd) <- gsub("^12mer$", "X12mer", colnames(kd))
    kd <- kd[, c("X12mer", "log_kd")]
    seed <- reverseComplement(DNAString(substr(mirseq, 2, 8)))
    seed <- paste0(as.character(seed), "A")
    w <- grep("X|N", kd$X12mer, invert = TRUE)
    pwm <- Biostrings::consensusMatrix(as.character(rep(kd$X12mer[w], 
      floor((exp(-kd$log_kd[w]))/3))), as.prob = TRUE, width = 12L)
    fields <- c("mer8", "fl.score")
    if (!all(fields %in% colnames(kd))) 
      kd <- .prep12mers(kd, seed = seed)
    fields <- c(fields, "log_kd")
    if (!all(fields %in% colnames(kd))) 
      stop("Malformed `kd` data.frame.")
    co <- t(vapply(split(kd[, c("log_kd", "fl.score")], kd$mer8), 
      FUN.VALUE = numeric(2), FUN = function(x) {
        .lm.fit(cbind(1, x$fl.score), x$log_kd)$coefficients
      }))
    
    cat("\n co dim", nrow(co)[1]," ", ncol(co)[1], "\n")
    
    fitted <- co[kd$mer8, 1] + co[kd$mer8, 2] * kd$fl.score
    if (!is.na(conservation)) {
      co <- .conservation_levels()
      if (!is.numeric(conservation)) {
        if (!(conservation %in% co)) {
          warning("Unknown conservation level - will be set to NA")
          conservation <- NA_integer_
        }
        else {
          conservation <- as.integer(names(co)[which(co == 
              conservation)])
        }
      }
      else {
        if (!(as.character(conservation) %in% names(co))) 
          warning("Unknown conservation level.")
      }
    }
    new("KdModel", list(mer8 = as.integer(round(co[, 1] * 1000)), 
      fl = as.integer(round(co[, 2] * 1000)), name = name, 
      mirseq = mirseq, canonical.seed = seed, pwm = pwm, conservation = conservation, 
      cor = cor(fitted, kd$log_kd), mae = median(abs(kd$log_kd - 
          fitted)), ...))
}


