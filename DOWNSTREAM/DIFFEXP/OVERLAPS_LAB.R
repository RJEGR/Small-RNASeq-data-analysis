# TEST OVERLAPS IRANGE

# GENOMIC SOURCE ===
# /Users/cigom/Documents/MIRNA_HALIOTIS/RNA_SEQ_ANALYSIS/REPEAT_MASKER_OUT
# multi_genome.newid.fa.out.gff

# /Users/cigom/Documents/MIRNA_HALIOTIS/ENSEMBLE
# Haliotis_rufescens_gca023055435v1rs.xgHalRufe1.0.p.56.gtf
# multi_genome.newid.gtf


# TRANSCRIPTOMIC SOURCE ===

# /Users/cigom/Documents/MIRNA_HALIOTIS/RNA_SEQ_ANALYSIS/ASSEMBLY/STRINGTIE/REFBASED_MODE/ANNOTATION
# transcripts.fa.transdecoder.gff3

# /Users/cigom/Documents/MIRNA_HALIOTIS/RNA_SEQ_ANALYSIS/ASSEMBLY/STRINGTIE/REFBASED_MODE
# transcripts.gtf

# SRNA SOURCE ===

# path <- '/Users/cigom/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out'
# Results.gff3	knownRNAs.gff3

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(GenomicRanges)
library(tidyverse)

read_gene_features <- function(f) {
  
  require(GenomicRanges)
  require(tidyverse)
  
  col_names <- c("seqname", 
    "source", 
    "feature", 
    "start", 
    "end", 
    "score", 
    "strand", 
    "frame", 
    "attribute")
  
  df <- read_tsv(f, col_names = F, comment = "#", na = ".")
  
  names(df) <- col_names
  
  # start_ <- df$start
  # 
  # end_ <- df$end
  # 
  ranges_ <- IRanges(start = df$start, end = df$end)
  
  seqnames_ <- df$seqname
  
  # score_ <- df$score
  # 
  # strand_ <- df$strand
  
  # strand_ <- gsub("[.]", "*", strand_)
  
  gr <- GRanges(Rle(seqnames_), 
    ranges = ranges_, strand = df$strand, score = df$score, 
    attribute = df$attribute, feature = df$feature)
  
  # genome(gr) <- basename(f)
  
  return(gr)
}

# RANGE 0 ====
# GENOMIC COORDS

path <- '/Users/cigom/Documents/MIRNA_HALIOTIS/ENSEMBLE/'
# 
#pattern <- "multi_genome.newid.gff3$" # <- to get sequence-region per chromosome

pattern <- "multi_genome.newid.gtf"

f <- list.files(path = path, pattern = pattern, full.names = T)

# gr <- read_gene_features(f)

# gr %>% as_tibble() %>% group_by(type) %>% count(gene_biotype) %>% view()

gr <- rtracklayer::import(f)

# gr <- gr[gr$type == "region"]

# Let's transform to up and downstream from the gene to define region to find miRNA source
# This step is to find sRNA targets sites in further steps (so omit now)

srna_promoters <- promoters(gr[gr$type == 'exon'], upstream = 1e3, downstream = 1e3)

# filter by something (Ex. feature) ?:


# gr[seqnames(gr) %in% "JALGQA010000001.1"]

# seqlengths(gr)
# seqinfo(gr)

#===== (Omit)


which_att <- c("gene_id", "gene_biotype")

.attr <- df$attribute

split_attr <- function(s, which_att = c("gene_id", "gene_biotype")) {
  
  # s <- df$attribute
  
  s <- sapply(strsplit(s, ";"), `[`, )
  
  keep <- grepl(paste0(which_att, collapse = "|"), s)
  
  s <- s[keep]
  
  s <- gsub(paste0(which_att, collapse = "|"), "", s)
  
  bind <- data.frame(t(s))
  
  names(bind) <- which_att
  
  return(bind)
}

df_sub <- df %>% filter(seqname %in% "JALGQA010000001.1")

.attr <- df_sub$attribute

att_df <- lapply(.attr, split_attr) # <- take long time if rows > 1000
att_df <- do.call(rbind, att_df)

df_sub <- cbind(df_sub, att_df) %>% as_tibble()

df_sub %>% count(gene_biotype)

which_feat <- "transcript"

df_sub %>% filter(grepl(which_feat, feature))

df_sub %>% group_by(feature) %>% distinct(start, end) %>% view()


df_sub %>% filter(start >= 44946 & end <= 72500)

# https://github.com/wilkox/gggenes
# https://bernatgel.github.io/karyoploter_tutorial/

library(gggenes)

set.seed(202355)

df_sub %>% 
  filter(start >= 44946 & end <= 72500) %>% view()
  # group_by(gene_biotype) %>% distinct(start, end) %>%
  ggplot(aes(xmin = start, xmax = end, y = feature, forward = F, fill = gene_biotype)) +
  geom_gene_arrow() 
# facet_grid(~ seqname)


# df_sub %>% filter(start >= 44946 & end <= 72500) %>%
df %>% group_by(seqname, feature) %>% sample_n(1) %>%
  # distinct(feature, .keep_all = T) %>% filter(feature != "gene") %>%
  # group_by(gene_biotype) %>% distinct(start, end) %>%
  ggplot(aes(xmin = start, xmax = end, y = seqname, forward = F, fill = feature)) +
  geom_gene_arrow() 
  # facet_grid(~ seqname)


gggenes::example_dummies %>%
  ggplot(aes(xmin = start, xmax = end, y = molecule,fill = gene)) +
  geom_gene_arrow()

# RANGE 1 ====

# After running shortstacks, lets to used the GenomicRanges to get the coordinates of the sRNA source within the genome. After obtaining the coordinates we used the GenomicRanges (Lawrence et al., 2013) resize function with the parameters fix = "center", width=100, this allowed us to get a window of interaction between sRNAs and genome

path <- '/Users/cigom/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out'

srna_f <- list.files(path = path, pattern = "Results.gff3", full.names = T)
res_f <- list.files(path = path, pattern = "Results.txt", full.names = T)

# MTD_f <- list.files(path = path, pattern = "METADATA.tsv", full.names = T)
# MTD <- read_tsv(MTD_f)

# The 'score' column in the Results.gff3 format stores the number of sRNA-seq aligned reads at that locus.

# gr2 <- read_gene_features(srna_f) 

gr2 <- rtracklayer::import(srna_f)

# FILTER BY TRUE MIRNS (NOVEL AND KNOWN)

query.ids <- read_tsv(res_f) %>%
  mutate(KnownRNAs = ifelse(is.na(KnownRNAs) & MIRNA == "Y", Name, KnownRNAs)) %>%
  mutate(MIRNA = ifelse(!is.na(KnownRNAs), "Y", MIRNA)) %>%
  drop_na(KnownRNAs) %>% pull(Name)
  # mutate(MIRNA = factor(MIRNA, levels = c("Y", "N")))

# gr2 <- gr2[gr2$ID %in% query.ids]


# Q: EN CUANTOS SCAFFOLDS DEL GENOMA DEL ABULON () SE ANOTARON FUENTES DE SRNAS FUNCIONALES ====

sum(seqlevels(gr2) %in% seqlevels(gr)) # 221 / 616 cromosomas

length(gr) # 1,613,804 (GTF) y 616 (GFF3)

length(gr <- gr[seqnames(gr) %in% seqlevels(gr2)]) # 221 GFF3

length(gr2) # 66,488

not_founded <- gr2[!seqnames(gr2) %in% seqlevels(gr)]

length(gr2 <- gr2[seqnames(gr2) %in% seqlevels(gr)]) # 66,449

sum(seqlevels(gr2) %in% seqlevels(gr)) # 205 / 221 scaffolds

srna_to_genome_ov <- findOverlaps(gr2, gr, minoverlap = 1)

# FIND OVERLAPS BY CHROMOSOME

## optional, if you want a genomic order of the chromosomes
gr0 = sortSeqlevels(gr2)

REPEAT_MASKER_OUT## split into a GRangesList
## where each element has all ranges for one chromosome
grl = split(gr0, seqnames(gr0))

x <- grl[[2]]



find_chr_overlaps <- function(query_gr, subject_gr) {
  
  gr <- subject_gr
  
  gr <- sortSeqlevels(gr)
  
  gr <- gr[seqnames(gr) == as.character(seqnames(query_gr)@values)]
  
  fo <- findOverlaps(query_gr, gr, minoverlap = 1)
  
  return(fo)
  
}


## apply a function to the ranges of each chromosome
res = lapply(grl, find_chr_overlaps)


res[[1]]

res %>% as_tibble()

findOverlaps(gr2, gr, minoverlap = 1)

# overlaps functions are strand-specific, although findOverlaps has an ignore.strand option.

hist(countOverlaps(gr2, gr, minoverlap = 1))

table(countOverlaps(gr2, gr))

# Hits object with 55030 hits and 0 metadata columns:

# Q: DONDE SE ORIGINAN LOS SRNAS? ====

to_index <- unique(subjectHits(srna_to_genome_ov)) # position vector where coordinates from subject is found

df <- gr[sort(to_index)] %>% 
  as_tibble() %>%
  group_by(type) %>%
  count(gene_biotype) %>% 
  mutate(type = as.character(type)) %>%
  # mutate(gene_biotype = ifelse(gene_biotype %in% "protein_coding", type, gene_biotype)) %>%
  drop_na() %>%
  arrange(desc(n)) %>%
  mutate(gene_biotype = factor(gene_biotype, levels = unique(gene_biotype)))



df %>%
  ggplot() +
  geom_col(aes(y = gene_biotype, x = n), fill = "grey20") +
  geom_col(data = subset(df, gene_biotype == "protein_coding"), 
    aes(y = gene_biotype, x = n,fill = type), 
    position = position_stack(reverse = T)) +
  # geom_text(aes(label = n), size = 2, hjust = -0.05, family = "GillSans") +
  scale_x_continuous(labels = scales::comma) +
  theme_bw(base_family = "GillSans", base_size = 14) +
  ggsci::scale_fill_aaas() +
  theme(
    # axis.text.x = element_text(angle = 45, 
    # hjust = 1, vjust = 1, size = 14), # color = axis_col
    axis.ticks.length = unit(10, "pt"),
    legend.position = 'top',
    panel.border = element_blank(),
    # axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank()) 

# 2)

miRNAs_loc <- c("MIRNA_hairpin", "mature_miRNA", "miRNA-star")

# Q: Que tipo de srnas fueron sobrelapados en regiones diferentes de "protein_coding"? ======

# .gr <- gr[sort(to_index)]

keep <- which(!gr$gene_biotype == "protein_coding" )

.gr <- gr[keep]

table(.gr$gene_biotype)

ov <- findOverlaps(gr2, .gr, minoverlap = 1)

from_index <- unique(queryHits(ov)) # position vector where coordinates from query is found

gr2[sort(from_index)] %>% 
  as_tibble() %>% 
  group_by(type, strand) %>% 
  summarise(n = n(), TotalReads = sum(score)) %>% # view()
  mutate(col = type) %>% 
  ungroup() %>% mutate_if(is.factor, as.character) %>%
  mutate(col = ifelse(col %in% miRNAs_loc, "miRNA_locus", col)) %>%
  mutate(col = ifelse(grepl("siRNA", col), "siRNA_locus", col)) %>%
  mutate(pct = n / sum(n)) %>%
  ggplot() +
  facet_grid(~ col, scales = "free", space = "free") +
  geom_col(aes(y = n, x = type, fill = strand), 
    position = position_dodge(width = 0.5), width = 0.4) +
  theme_bw(base_family = "GillSans", base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, 
    hjust = 1, vjust = 1, size = 14))

# OMIT =====
# For a GRangesList, overlap detection reports a hit at the element level, i.e., when any range within an element overlaps a query range. This semantic is convenient, for example, when counting the total number of RNA-seq read pairs overlapping the exonic regions of each transcript. In that case, both the reads and the transcripts are GRangesList objects.

GRangesList(gr2, gr)

# overlaps <- countOverlaps(gr2, gr)

# gr2 %>% head(10)

# head(fo <- findOverlaps(gr2, gr))
# head(countOverlaps(gr2, gr))


queryHits(fo) # position vector where coordinates from query is found

subjectHits(fo) # position vector where coordinates from subject is found

# ranges from the query for which we found a hit in the subject

index = queryHits(fo)

ovlps <- gr[index,]

srna_df[index,]



# GRangesList(gr, gr2)

coverage(gr)




# strand_ <- gsub("[.]", "*", strand_)

miRNAs_loc <- c("MIRNA_hairpin", "mature_miRNA", "miRNA-star")

table(srna_df$X3)

gr2 %>% as_tibble() %>% 
  group_by(feature, strand) %>% 
  # tally(score, sort = T) %>%
  # group_by(X3, X7) %>% 
  summarise(n = n(), TotalReads = sum(score)) %>%
  mutate(col = feature) %>% 
  mutate(col = ifelse(col %in% miRNAs_loc, "miRNA_locus", col)) %>%
  mutate(col = ifelse(grepl("siRNA", col), "siRNA_locus", col)) %>%
  ungroup() %>%
  mutate(pct = n / sum(n)) %>%
  ggplot() +
  facet_wrap(~ col, scales = "free") +
  geom_col(aes(y = n, x = feature, fill = strand), 
    position = position_dodge(width = 0.5), width = 0.4) +
  theme_bw(base_family = "GillSans", base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, 
    hjust = 1, vjust = 1, size = 14))


# THIS STEP IS TO FILTER ACCORDING TO SRNA UP-XPRESSED BY CONDITION

srna_df %>% select(X9) %>% 
  separate(col =X9, sep = ";", into = c("Name", "DicerCall", "MIRNA")) %>%
  count(MIRNA)


# RANGE 2 =====

path <- "/Users/cigom/Documents/MIRNA_HALIOTIS/RNA_SEQ_ANALYSIS/REPEAT_MASKER_OUT"

mask_f <- list.files(path = path, pattern = "multi_genome.newid.fa.out.gff", full.names = T)

# gr3 <- read_gene_features(mask_f) 
gr3 <- rtracklayer::import(mask_f)

# mask_df <- read_tsv(mask_f, col_names = F, comment = "#")

# SEARCH OVERLAPS

# gr %over% gr2
# gr1[gr1 %over% gr2]

srna_to_genome_ov

# te_to_genome_ov <- findOverlaps(gr3, gr, minoverlap = 1)

srna_to_te_ov <- findOverlaps(gr2, gr3, minoverlap = 1)

# Hits object with 2656 hits and 0 metadata columns:

unique(queryHits(fo)) # position vector where coordinates from query is found

unique(subjectHits(fo)) # position vector where coordinates from subject is found

# ranges from the query for which we found a hit in the subject

index = queryHits(fo)

ovlps <- gr[index,]

srna_df[index,] %>% 
  group_by(X3, X7) %>% 
  summarise(n = n(), TotalReads = sum(X6)) %>%
  mutate(col = X3) %>% 
  mutate(col = ifelse(col %in% miRNAs_loc, "miRNA_locus", col)) %>%
  mutate(col = ifelse(grepl("siRNA", col), "siRNA_locus", col)) %>%
  ungroup() %>%
  mutate(pct = n / sum(n)) %>%
  ggplot() +
  geom_col(aes(y = n, x = col, fill = X7), 
    position = position_dodge(width = 0.5), width = 0.4) +
  theme_bw(base_family = "GillSans", base_size = 12)

# 

coverage(gr)
# maxPos <- which.max(ctcfCoverage10)
# 
# > roi <- resize(IRanges(maxPos, width = 1), 5000, “center”)
# 
# > roiCoverage <- ctcfCoverage$chr10[roi]

# RANGE 3 ====

path <- "/Users/cigom/Documents/MIRNA_HALIOTIS/RNA_SEQ_ANALYSIS/ASSEMBLY/STRINGTIE/REFBASED_MODE/"


orf_f <- list.files(path = paste0(path, "ANNOTATION"), 
  pattern = "transcripts.fa.transdecoder.gff3", full.names = T)

orf_df <- read_tsv(df_f, col_names = F, comment = "#")

gtf_f <- list.files(path = path, pattern = "^transcripts.gtf", full.names = T)

gtf_df <- read_tsv(gtf_f, col_names = F, comment = "#")

srna_df %>% select(X9) %>% 
  separate(col =X9, sep = ";", into = c("Name", "DicerCall", "MIRNA")) %>%
  count(MIRNA)


orf_df %>% distinct(X1) %>% filter(!grepl("LOC", X1))

orf_df %>% head() %>% view()
gtf_df %>% filter(grepl("LOC", X9)) %>% head() %>% view()

start_ <- df$X4

end_ <- df$X5

score_ <- mask_df$X6

ranges_ <- IRanges(start = start_, end = end_)

seqnames_ <- mask_df$X1

strand_ <- mask_df$X7

strand_ <- gsub("[.]", "*", strand_)

gr3 <- GRanges(Rle(seqnames_), ranges = ranges_, strand = strand_)


# WHICH MITOCONDRIAL
# JALGQA010000616.1 <- MITOCHONDR

srna_df %>% filter(X1 %in% "JALGQA010000616.1") %>%
  group_by(X3) %>%
  summarise(n = n(), TotalReads = sum(X6))

mask_df %>%
  filter(X1 %in% "JALGQA010000616.1") %>%
  group_by(X3) %>%
  summarise(n = n(), TotalReads = sum(X6))


# karyoploteR ====

library(karyoploteR)
# https://bernatgel.github.io/karyoploter_tutorial//Tutorial/PlotDensity/PlotDensity.html
# PLOT DENSITY AND REGIONS
# cd /usr/local/lib
# sudo ln -s /usr/local/gfortran/lib/libquadmath.0.dylib .
# sudo ln -s /usr/local/gfortran/lib/libgfortran.5.dylib

# srna_to_genome_ov <- findOverlaps(gr2, gr, minoverlap = 1)

# to_index <- unique(subjectHits(srna_to_genome_ov)) # WHICH IN THE GENOME 
# from_index <- unique(queryHits(srna_to_genome_ov)) # WHICH IN THE SRNAS

# genome <- gr[sort(to_index)]

# query.locus <- seqnames(genome)@values[seqnames(genome)@lengths > 1]

# gr2 <- gr2[gr2$ID %in% query.ids]

# query.locus <- seqnames(gr2)@values[seqnames(gr2)@lengths > 1]

viz <- gr2 %>% as_tibble() %>% 
  group_by(seqnames) %>%
  summarise(n = n(), TotalReads = sum(score)) %>% 
  arrange(desc(n)) 

viz %>%
  ggplot(aes(log10(n), log10(TotalReads))) + geom_point()

query.locus <- viz %>% pull(seqnames) %>% head(30)

# GRange object with seqnames, ranges and strand cols at least

genome <- gr[seqnames(gr) %in% unique(unfactor(query.locus))]

seqlevels(genome) <- unique(unfactor(seqnames(genome)@values))

chr.names <- as.character(GenomeInfoDb::seqnames(genome))
any(duplicated(chr.names))

kp <- plotKaryotype(genome = genome, chromosomes = query.locus) 

# kp <- plotKaryotype(genome, plot.type=2, chromosomes = "JALGQA010000015.1")

data <- gr2 # gr2[gr2$ID %in% query.ids]

kpPlotDensity(kp, data=data)

kpPlotRegions(kp, data=data, data.panel=1, avoid.overlapping = T)

kp <- plotKaryotype(genome, plot.type=2, chromosomes = "JALGQA010000015.1")

kpPlotCoverage(kp, data=data)

kpPlotDensity(kp, data=data, data.panel=2)
