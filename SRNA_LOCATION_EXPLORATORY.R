# RICARDO GOMEZ-REYES  

# THE PURPOSE IS DEFINE IF THOSE SRNA LOCI ARE INTERGENIC OR NOT
# I.E FOR EACH CLUSTER (66226) ANNOTATE GENOMIC LOCATION (INTRONIC, EXONIC, ETC.)
# AIM: FIND GENOMIC COODS OVERLAPS BETWEEN SRNA LOCI AND ANNOT. AS INTERGENIC, TRANSCRIPT UNIT (TU)/GENE HOST

# OUTPUTS WILL GENERATE OBJECT AS DATA.FRAME FORMAT W/ FOLLOW COLS:
# seqnames start end strand type  biotype Name gene_id transcript_id

# GENE_ID IS CONCORDANCE W/ mir_vs_utr_rmdup_RNAhybrid.out.psig.gene_features.rds

# REFERENCE: 
# doi: 10.1101/gr.2722704: These findings strongly suggest that miRNAs are transcribed in parallel with their host transcripts, and that the two different transcription classes of miRNAs (`exonic' and `intronic') identified here may require slightly different mechanisms of biogenesis.

# (DOI: 10.1111/mec.14973) miRNA locations that were characterized as intergenic were those that did not intersect with any annotated gene in the genome, while intronic locations were those that intersected with a gene but were not annotated as either an exon or an UTR 


# (doi.org/10.1186/1471-2164-11-533): Roughly half of known miRNA genes are located within previously annotated protein-coding regions ("intragenic miRNAs"). A high-confidence set of predicted mRNA targets of intragenic miRNAs also shared many of these features with the host genes. Approximately 20% of intragenic miRNAs were predicted to target their host mRNA transcript.

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

library(GenomicRanges)

# 1) SRNA-SEQ TRANSCRIPTOME ====

wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out"

pattern <- "Results.gff3"

f <- list.files(path = wd, pattern = pattern, full.names = T)

assembly <- rtracklayer::import(f)


# Remember:

length(assembly) # 66488

length(IRanges::reduce(assembly)) # 66226

length(range(assembly)) # 528

# IRanges::ranges(assembly, use.mcols = T) %>% as_tibble() %>% distinct(start, end)

# IRanges::resize(assembly, width(assembly), fix = "center") %>% as_tibble() %>% distinct(start, end)

# length(findOverlaps(range(assembly)))

# 2) LOAD GENOME REFERENCE ====
#
pattern <- "multi_genome.newid.gff3$" # <- to get sequence-region per chromosome

wd_genome <- "/Users/cigom/Documents/MIRNA_HALIOTIS/ENSEMBLE/"

f <- list.files(path = wd_genome, pattern = pattern, full.names = T)

genome <- rtracklayer::import(f)

length(g <- genome[!genome$type == "region"])

# length(g)

table(g$biotype)

table(g$type)

# 3) LOAD REPEAT MASK GENOME ANNOTATION ====

wd_rmask <- "/Users/cigom/Documents/MIRNA_HALIOTIS/REPEAT_MASKER_OUT/"

pattern <- "multi_genome.newid_repeat_masked.gff3" # "multi_genome.newid.fa.out.gff"

mask_f <- list.files(path = wd_rmask, pattern = pattern, full.names = T)

maskg <- rtracklayer::import(mask_f)

table(maskg$type) 

maskdf <- maskg %>% as_tibble()

f <- list.files(path = wd_rmask, pattern = "multi_genome.newid.fa.out$", full.names = T)

mask_features <- read_table(f, col_names = F, skip = 3) %>% as_tibble()

colNames <- c("bit_score", "perc_div", "perc_del", "perc_ins", "seqnames", "start", "end", "left", "strand",
  "matching_repeat", "repeat_family", "start_r", "end_r", "left", "ID")

colnames(mask_features) <- colNames

nrow(maskdf)

nrow(mask_features)

identical(as.character(maskdf$ID), as.character(mask_features$ID))

mask_features <- mask_features %>% select(matching_repeat, repeat_family )

maskdf <- cbind(maskdf, mask_features)

maskdf <- data.frame(seqnames = maskdf$seqnames, 
  start = maskdf$start, end = maskdf$end, 
  strand = maskdf$strand, 
  source = maskdf$source,
  type = maskdf$repeat_family,
  biotype = maskdf$matching_repeat,
  score = maskdf$score, 
  phase = maskdf$phase,
  # description = maskdf$type, 
  gene_id = NA, transcript_id = NA)

# 3.1) BIND REPEAT MASK TO GENOME: ====

which_names <- maskdf %>% names()

.g <- genome %>% as_tibble() %>% select(all_of(which_names))

.g <- as_tibble(rbind(.g, maskdf))

.g <- GRanges(Rle(.g$seqnames), 
  ranges =  IRanges(start = .g$start, end = .g$end), 
  strand = .g$strand, 
  source = .g$source,
  type = .g$type,
  biotype = .g$biotype,
  score = .g$score, 
  phase = .g$phase,
  # description = .g$type, 
  gene_id = .g$gene_id, transcript_id = .g$transcript_id)

fileout <- paste0(wd_rmask, "multi_genome.newid_repeat_masked.gff3")

write(rtracklayer::export(.g, format = "gff3"), file = fileout)

.g <- rtracklayer::import(fileout)

length(g <- .g[!.g$type == "region"])

# 
# 4) =====
# Use cluster 6 as example:

subsetByOverlaps(x = g, ranges = assembly[6,], type="any") %>% 
  as_tibble() %>% 
  mutate(Name = "Cluster_6") %>%
  select(seqnames, start, end, strand, type, biotype, Name, gene_id, transcript_id) 

# group_by(biotype) %>%
# summarise(across(type, .fns = paste_headers), .groups = "drop_last") 

g

out <- list()

max <- length(assembly)

pb <- txtProgressBar(min = 0, max = max, style = 3, width = 50, char = "=")

for(i in 1:max) { 
  
  j <- i
  
  Name <- assembly[j,]$ID
  
  ranges <- assembly[j,]
  
  L <- subsetByOverlaps(x = g, ranges = ranges, type="any")
  
  if(length(L) == 0) {
    
    t <- "Intergenic"
    
    s <- start(ranges)
    e <- end(ranges)
    
    seqnames <- unfactor(seqnames(ranges)@values)
    
    strand_ <- unfactor(strand(ranges)@values)
    
    gene_id <- NA
    transcript_id <- NA
    
    out[[j]] <- data.frame(seqnames, start = s, end = e, strand = strand_, 
      type = t, biotype = t, Name, gene_id, transcript_id)
  } else
    
    out[[j]] <- L %>% as_tibble() %>% mutate(Name = Name) %>%
    select(seqnames, start, end, strand, type, biotype, Name, gene_id, transcript_id) 
  
  setTxtProgressBar(pb, i)
  
}

# df <- lapply(file, read.table) # leer todas las tablas

paste_headers <- function(x) { 
  
  x <- unique(x)
  x <- paste(x, sep = '|', collapse = '|') }

length(out) == max

head(out)


# df <- do.call(rbind, lapply(out, function(x) { setNames(x, names(out[[6]])) }))
# df <- df %>% as_tibble()

# df <- do.call(rbind, out) %>% as_tibble()

write_rds(out, file = paste0(wd, "/SRNA_LOCATION_OUT.rds"))

head(read_rds(paste0(wd, "/SRNA_LOCATION_OUT.rds")))

df %>% count(type)

df %>% count(Name)

df %>%
  # filter(type != "Intergenic")
  group_by(Name, biotype) %>%
  summarise(across(type, .fns = paste_headers), .groups = "drop_last")

wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out"

write_rds(df, file = paste0(wd, "/SRNA_LOCATION.rds"))

# ALSO INCLUDE INTERGENIC TERM:
# https://support.bioconductor.org/p/54789/
### Removing the strand to perform the reduce operation
genes <- range(TU[strand(TU) != "*"])

strand(genes) <- '*'

### Now, let's find the gaps
genic <- IRanges::reduce(genes)

### Because of teh reduce, there should be no overlaping genic regions
# MUST BE ZERO
length(findOverlaps(genic))

### Finding the intergenic regions
intergenic <- gaps(genic)

length(findOverlaps(intergenic))


# OR

length(findOverlaps(IRanges::setdiff(g,genic)))

# 

# PROCESS FEATURES:
# 1) SET INTERGENIC REGIONS ====
# gaps(x) produces a set of ranges covering the positions in [start(x), end(x)] 
# that are not covered by any range in x. 
# Given coding sequence addresses and exon intervals, this can be used to enumerate introns.
?"intra-range-methods"
?"inter-range-methods"



# EX:
# genome %>% as_tibble() %>% head() %>% view()
genome %>% as_tibble() %>% 
  filter(type != "region") %>% 
  group_by(seqnames, start) %>% 
  count(type, sort = T)  

genome %>% 
  as_tibble() %>% 
  filter(seqnames == "JALGQA010000004.1" & start == 9291699) %>% view()

subsetByOverlaps(x = genome, ranges = CDS[1,1])

length(CDS <- genome[genome$type == "CDS"])

length(CDS <- CDS[seqnames(CDS) %in% query.locus])

seqlevels(CDS) <- unique(unfactor(seqnames(CDS)@values))

ranges <- IRanges::ranges(CDS)

# NUMBER OF CDS PER CHROMOSOME?

CDS %>% as_tibble() %>% count(seqnames, sort = T) 

# gaps(ranges)


gaps(CDS)[1]

# CDS

subsetByOverlaps(x = genome, ranges = mirs_true) %>% 
  as_tibble() %>% 
  # group_by(seqnames, start, end) %>%
  # summarise(across(type, .fns = paste_headers), .groups = "drop_last") %>%
  count(type, sort = T) 

# 46718-46755

# g3 <- genome[!genome$type == "region"]

# length(g3 <- g3[seqnames(g3) %in% query.locus])

# seqlevels(g3) <- unique(unfactor(seqnames(g3)@values))

# g3 %>% as_tibble() %>% count(type, sort = T)

## split into a GRangesList
## where each element has all ranges for one chromosome

# g3 <- split(g3, seqnames(g3))
## apply a function to the ranges of each chromosome
# lapply(g3, IRanges::ranges)

ranges <- IRanges::ranges(g3)

gaps(g3, start = start(ranges), end = end(ranges))

# g3 <- subsetByOverlaps(x = g3, ranges = mirs_true) # <- Overlaps se hace despues de reducir o de conocer el tamano del gen source (no necesario)

g3 %>% as_tibble() %>% count(type, sort = T)

IRanges::reduce(g3) # <- NECESITAMOS REDUCIR, PERO CONSERVAR LOS FEATURES:

IRanges::ranges(g3, use.mcols = T) #, use.mcols = T) 

IRanges::resize(g3, width(g3), fix = "center") # %>% as_tibble() %>% count(type, sort = T)

x <- IRanges::resize(g3, width(g3), fix = "center") %>%
  as_tibble() %>%
  group_by(seqnames, start, end) %>%
  summarise(across(type, .fns = paste_headers), .groups = "drop_last") 

source_annot_df <- GRanges(Rle(x$seqnames), 
  ranges =  IRanges(start = x$start, end = x$end), type = x$type)




# GenomicRanges::setdiff(g3, g3)

# %>% as_tibble() %>% count(type, sort = T)


