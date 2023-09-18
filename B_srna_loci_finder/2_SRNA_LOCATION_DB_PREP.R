
# RICARDO GOMEZ-REYES  
# AFTER RUN 1_SRNA_LOCATION_DB_PREP.R
# 

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

library(GenomicRanges)

wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out"

out <- read_rds(paste0(wd, "/SRNA_LOCATION_OUT.rds"))

print(.out <- do.call(rbind, out) %>% as_tibble())

# 1) HIERARCHICAL CATEGORIES FOR THE GENOMIC ANNOTATIONS ==== 

# SIMILAR TO ISAAC MARTINEZ UGALDE (CEI ABREU STUDENT)
# GENOMIC ANNOTATION MUST BE CURATED BASED ON 
# OVERLAPPING FEATURES ARE HIERARCHICALY CATEGORIZED FOR EACH GENOMIC REGION.
# THEN, ARE USED TO IDENTIFY MIRNA/SRNAS SOURCE
# 1) USING SETDIFF (FROM GENOMICRANGES; Lawrence et al., 2013) TO OBTAIN THE NON-INTERSECTED REGION AND THE INTERSECTED
# 1.2) GENOMIC CATEGORIES ARE SUMMARIZED IN biotypeL OBJECT
# 1.3) THE LEVEL OF SRNA PRODUCTION IN EACH GENOMIC REGION WAS OBTAINED USING subsetByOverlaps


# 1.1) CURATE NA biotype ====
# CURATE NA biotype based on type 

.out %>% filter(is.na(biotype)) %>% count(type, sort = T)

.out <- .out %>% 
  mutate(biotype = ifelse(type == "exon", type, biotype)) %>% 
  mutate(biotype = ifelse(type == "CDS", type, biotype)) %>% 
  mutate(biotype = ifelse(grepl("prime_UTR", type), "UTR", biotype))

# Sanity check
.out %>% filter(is.na(biotype)) %>% count(type, sort = T) # MUST BE EMPTY


# 1.2) MASK INTRAGENIC ====
# (DOI: 10.1111/mec.14973) miRNA intronic (intragenic) locations were those that intersected with a gene but were not annotated as either an exon or an UTR 
  
.out %>% count(biotype, sort = T) 

.out %>% filter(biotype == "protein_coding") %>% count(type, sort = T)
.out %>% filter(biotype == "misc_RNA") %>% count(type, sort = T)
.out %>% filter(biotype == "pseudogene") %>% count(type, sort = T)

recode_to <- c(`protein_coding` = "Intragenic", `misc_RNA`= "Intragenic",`pseudogene` = "Intragenic", `mRNA` = "Intragenic", `gene` = "Intragenic")

.out <- .out %>% 
  # mutate(SRNA_CLASS = biotype) %>%
  dplyr::mutate(biotype = dplyr::recode_factor(biotype, !!!recode_to)) %>%
  dplyr::mutate(type = dplyr::recode_factor(type, !!!recode_to))

known_repeat_ranks <- c(
  "Simple_repeat","LTR/ERV1", "LTR/Gypsy", "DNA/TcMar-Tc2", "Low_complexity",
  "LINE/CR1", "DNA/TcMar-Mariner", "DNA/hAT-Charlie", "Satellite", 
  "LTR/ERVL-MaLR", "DNA/TcMar-Tigger", "DNA/MULE-MuDR", "Unknown", 
  "Satellite/centr", "LINE/RTE-X", "SINE/MIR", "LINE/L2", "LINE/L1") 

recode_to <- structure(rep("Repeat element", 
  length(known_repeat_ranks)), names =known_repeat_ranks)

.out <- .out %>% 
  # mutate(SRNA_CLASS = biotype) %>%
  dplyr::mutate(biotype = dplyr::recode_factor(biotype, !!!recode_to))


.out %>% count(biotype, sort = T) 

# WHAT ABOUT NON CODING BIOTYPES??
# ARE ALSO INTRAGENIC? MAY BE TRUE!
# OR BETTER, REMOVE FROM DOWNSTREAM ANALYSIS!

.out %>% filter(biotype == "rRNA") %>% count(type, sort = T)
.out %>% filter(biotype == "lncRNA") %>% count(type, sort = T)
.out %>% filter(biotype == "tRNA") %>% count(type, sort = T)

.out %>% filter(biotype == "Simple_repeat") %>% count(type, sort = T)
.out %>% filter(biotype == "srpRNA") %>% count(type, sort = T)
.out %>% filter(biotype == "snoRNA") %>% count(type, sort = T)
.out %>% filter(biotype == "snRNA") %>% count(type, sort = T)
.out %>% filter(biotype == "Unknown") %>% count(type, sort = T)

# 1.3) PREP RANK BIOTYPE HIERARCHIC CATEGORIZE ====

# (doi.org/10.1101/2022.12.30.522274) Non redundant (nr) genomic annotation mask may be build to reduce the Intra-range features (i.e overlapped annotations) and assign mapped sequences to unique/single annotations. Using preference (or hierarchical) annotation selection as follow: rRNA > tRNA, Transposable elements (TE) > protein-coding exon, other ncRNAs, introns, pseudogenes ... 


known_non_coding_ranks <- c("tRNA", "rRNA","lncRNA", "snoRNA", "srpRNA", "snRNA")

primary_ranks <- c("Repeat element", "UTR","exon", "CDS")

last_ranks <- c("Intergenic", "Intragenic")

biotypeL <- c(primary_ranks, known_non_coding_ranks, last_ranks)

.out <- .out %>% mutate(biotype = factor(biotype, levels = biotypeL))

# sanity check

identical(levels(.out$biotype), biotypeL) # MUST BE TRUE

pick_rank <- function(x) {
  
  x <- unique(x[!is.na(x)])
  x <- levels(x)[sort(as.integer(x))][1]
}

# USE THEN IN STEP 4)

paste_headers <- function(x) { 
  x <- x[!is.na(x)] 
  n <- length(x)
  x <- unique(x)
  
  x <- paste(x, sep = '|', collapse = '|') 
  x <- paste(n, x, sep = '|', collapse = '|') }


# out %>% drop_na(gene_id)

# MERGE W/ SHORTSTACKS OUTPUT =====

# SRNA-SEQ TRANSCRIPTOME

wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out"

pattern <- "Results.gff3"

f <- list.files(path = wd, pattern = pattern, full.names = T)

assembly <- rtracklayer::import(f)

# SPLIT TRUE MIRS/PIRS

SRNAS <- read_rds(paste0(wd, "/SPLIT_SHORTSTACKS_MIRS_PIRS_SIRS.rds"))

SRNAS %>% count(SRNAtype)

str(which_pirs <- SRNAS %>% filter(SRNAtype == "piR") %>% distinct(Name) %>% pull())
str(which_mirs <- SRNAS %>% filter(SRNAtype == "miR") %>% distinct(Name) %>% pull())

# PREPARE DB FOR VIZ ====
# 1)
# TO KEEP TRACK LOCUS TYPE PER CLUSTER
# Note: be aware score (Reads) cols for *.mature/*.star locus are documented in assembly (262 rows) but not in res
# 
# assembly %>% as_tibble() %>% filter(grepl(".mature|.star", ID)) %>% select()
  
.DB <- assembly %>% as_tibble() %>% select(ID, type) %>% dplyr::rename("Locus_type"="type", "Name"="ID")

# 2)  BIND SEVERAL OUTPUT  ====
# TO KEEP TRACK READS COUNT PER CLUSTER

f <- list.files(path = wd, pattern = "Results.txt", full.names = T)

DB <- read_tsv(f) %>% 
  select(Locus, Name, Chrom, Start, End, Strand, DicerCall, Reads, UniqueReads, FracTop, 
    MajorRNA, MajorRNAReads, MIRNA, KnownRNAs) %>%
  left_join(.DB, by = "Name")


# 3) BIND KNOWN RNA LABEL =====

# SRNAS %>% mutate(RNAType = ifelse(grepl("piR", KnownRNAs)))

.DB <- SRNAS %>% select(-MajorRNA, -KnownRNAs) %>% dplyr::rename("NKnownRNAs"="n")

DB <- DB %>% left_join(.DB, by = "Name")

# 4) BIND LOCATION DB ====


bind_types <- function(x) { x <- paste(unique(x[!is.na(x)]), sep = '|', collapse = "|") }

bind_ids <- function(x) { 
  
  x <- x[!is.na(x)]
  
  x <- unique(x)
  x <- paste(x, sep = ';', collapse = ';') }

.DB <- .out %>%
  group_by(Name) %>%
  mutate(biotype_best_rank = biotype) %>%
  summarise(
    across(biotype_best_rank, .fns = pick_rank),
    across(type, .fns = bind_types), 
    across(biotype, .fns = bind_types),
    across(gene_id, .fns = bind_ids),
    across(transcript_id, .fns = bind_ids),
    .groups = "drop_last") %>% 
  ungroup()

DB <- DB %>% left_join(.DB, by = "Name")

# 5) ADD COUNT-MATRIX-BASED CLUSTERING  ====

# LOAD MODULES

bwnet <- readRDS(paste0(wd, "/2023-06-26/bwnet.rds"))

bwmodules = WGCNA::labels2colors(bwnet$colors)

names(bwmodules) <- names(bwnet$colors)

DB$WGCNA <- NA

# Sanity check

any(names(bwmodules) %in% DB$Name)

DB$WGCNA  <- bwmodules[match(DB$Name, names(bwmodules))]

# DB %>% head() %>% view()

write_tsv(DB, file = paste0(wd, "/RNA_LOCATION_DB.tsv"))

# FURTHER: ADD TARGET: ====

# ADD INNER JOIN ANTIJOIN FROM TARGET-SCAN AND RNAHYBRID:

# head(read_tsv(paste0(wd, "/RNA_LOCATION_DB.tsv")))


