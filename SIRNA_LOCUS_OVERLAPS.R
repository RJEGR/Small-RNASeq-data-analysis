# SEARCH OVERLAPINGS DISTANCES BASED ON MIRS (NOVEL AND KNOWN) AND PIRS ANNOTATION:
# SEE (Marco et al 2014, 10.1093/nar/gkt534)

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

library(GenomicRanges)

library(karyoploteR)

# SPLITTED DATA:

wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out"

SRNAS <- read_rds(paste0(wd, "/KNOWN_CLUSTERS_MIRS_PIRS.rds"))

str(which_pirs <- SRNAS %>% filter(grepl("piR", KnownRNAs)) %>% distinct(Name) %>% pull())
str(which_mirs <- SRNAS %>% filter(!grepl("piR", KnownRNAs)) %>% distinct(Name) %>% pull())

# SRNA-SEQ TRANSCRIPTOME ====

pattern <- "Results.gff3"

f <- list.files(path = wd, pattern = pattern, full.names = T)

assembly <- rtracklayer::import(f)

mirs_true <- assembly[which(assembly$MIRNA == "Y")]

# GENOME REFERENCE ====
#
pattern <- "multi_genome.newid.gff3$" # <- to get sequence-region per chromosome

wd_genome <- "/Users/cigom/Documents/MIRNA_HALIOTIS/ENSEMBLE/"

f <- list.files(path = wd_genome, pattern = pattern, full.names = T)

genome <- rtracklayer::import(f)

region <- genome[genome$type == "region"]

chr.names <- as.character(GenomeInfoDb::seqnames(region))

any(duplicated(chr.names)) # MUST BE FALSE

# EXAMPLE 1): ====
# SET 131 COORDINATE FROM HAIRPINS AND OVERLAPS INTO MATURE AND STAR COORDINATES 

mirs_true <- subsetByOverlaps(x = assembly, ranges = mirs_true)

str(chr.names <- unique(unfactor(seqnames(mirs_true)@values)))

seqlevels(mirs_true) <- chr.names

# PROCESS GENOME ACCORDING TO MIR-LOCI ====

query.locus <- unique(seqlevels(mirs_true))

length(g <- region[seqnames(region) %in% query.locus])

seqlevels(g) <- unique(unfactor(seqnames(g)@values))

chr.names <- as.character(GenomeInfoDb::seqnames(g))

any(duplicated(chr.names))

mirs_true <- subsetByOverlaps(x = assembly, ranges = mirs_true)


# SEE INTER-microRNA DISTANCE =====

str(query_names <- c(which_pirs, which_mirs))

assembly %>% as_tibble() %>% 
  # filter(seqnames %in% query.locus) %>%
  filter(!type %in% c("mature_miRNA", "miRNA-star", "MIRNA_hairpin")) %>%
  filter(!ID %in% query_names) %>%
  count(type, sort = T)

str(siRNA_locus <- assembly %>% as_tibble() %>%
    filter(!type %in% c("mature_miRNA", "miRNA-star", "MIRNA_hairpin")) %>%
    filter(!ID %in% query_names) %>% distinct(ID) %>% pull())

# SELECT CHROMOSOME-RELATED MIR-LOCI

length(g2 <- assembly[seqnames(assembly) %in% query.locus])

# SELECT ONLY siRNA_locus
length(g2 <- g2[which(g2$ID %in% siRNA_locus)])

seqlevels(g2) <- unique(unfactor(seqnames(g2)@values))

# NUMBER OF siRNA_locus (NOT MIRS NOR PIRS) RELATED TO MIR-LOCI

g2 %>% as_tibble() %>% count(type, sort = T) 

# # NUMBER OF siRNA_locus (NOT MIRS NOR PIRS) RELATED TO MIR-LOCI

# IRanges::reduce(g2) # <- those are already reduced (i.e not overlapped) !
# IRanges::ranges(g2)

g2 <- GenomicRanges::resize(g2, width(g2), fix = "center")

mirs_true <- GenomicRanges::resize(mirs_true, width(mirs_true), fix = "center")

# PREPARE mycytobands.txt ====
# Cytoband file format is used to define the chromosome ideograms for a reference genome,
# Including follow colNames:
# chr	start	end	name	gieStain
# The gieStain levels are the ones used at UCSC: gneg, gpos25, gpos75, gpos100, gvar, acen, stalk, etc.

# https://software.broadinstitute.org/software/igv/Cytoband

cytobands <- g3 %>% as_tibble() %>% select(seqnames, start, end, strand, type, biotype)

cytobands <- GRanges(Rle(cytobands$seqnames), 
  ranges =  IRanges(start = cytobands$start, end = cytobands$end), 
  strand = cytobands$strand, type = cytobands$type, cytobands = cytobands$biotype)

# VIEW INTER-DISTANCE BETWEEN MIR-LOCI ====

# plotDefaultPlotParams(plot.type=1)

pp <- getDefaultPlotParams(1)

# pp$data2outmargin <- 100
# pp$leftmargin <- 0.2
# kp <- plotKaryotype(g, plot.type = 6, plot.params = pp)

# png(filename = paste0(wd, "/INTER-DISTANCE-MIR-LOCI.png"), width = 15, height = 15)

kp <- plotKaryotype(genome = g, plot.type=1, cex = 0.5,  
  # cytobands = cytobands, # Include cytobands relaled to genomic features
  plot.params = pp, main = "INTER-DISTANCE MIR-LOCI")

kpDataBackground(kp, color = "#FFFFFFAA", r0 = -0.5, r1 = 0)

kpAddBaseNumbers(kp, add.units = T, units = "Mb")

kpPlotRegions(kp,g2, col = "#647FA4", r0 = 0,r1 = 0.5)

kpPlotRegions(kp,mirs_true, col="red", r0 = 0,r1 = 0.5)


dev.off()

# kpPlotRegions(kp, g3, border="red", r0 = 0,r1 = 0.5)

# CONCLUSION:
# 1)  We (may) define a cluster of microRNAs as a group of microRNA precursors with an inter-microRNA distance of <10 kb on the same genomic strand (Marco et al 2014, 10.1093/nar/gkt534)

# LABEL BY MIR (Ex.) ====

text_df <- SRNAS %>% filter(!grepl("piR", KnownRNAs) & !KnownRNAs == "Novel") %>%
  mutate(KnownRNAs = ifelse(n > 1, paste0("Conserved|", n), KnownRNAs)) %>%
  left_join(assembly %>% as_tibble(), by = c("Name"="ID"))

text_df <- GRanges(Rle(text_df$seqnames), 
    ranges =  IRanges(start = text_df$start, end = text_df$end), 
  strand = text_df$strand, score = text_df$score, KnownRNAs = text_df$KnownRNAs)

kpPlotMarkers(kp, data = mirs_true, labels = text_df$KnownRNAs)

# LABEL BY MIRS-SOURCE: ====


# PROCESS FEATURES:
# 1) SET INTERGENIC REGIONS
# gaps(x) produces a set of ranges covering the positions in [start(x), end(x)] 
# that are not covered by any range in x. 
# Given coding sequence addresses and exon intervals, this can be used to enumerate introns.

g3 <- genome[!genome$type == "region"]

length(g3 <- g3[seqnames(g3) %in% query.locus])

seqlevels(g3) <- unique(unfactor(seqnames(g3)@values))

# g3 %>% as_tibble() %>% count(type, sort = T)

## split into a GRangesList
## where each element has all ranges for one chromosome

# g3 <- split(g3, seqnames(g3))
## apply a function to the ranges of each chromosome
# lapply(g3, IRanges::ranges)

ranges <- IRanges::ranges(g3)

gaps(g3, start = start(ranges), end = end(ranges))

# g3 <- subsetByOverlaps(x = g3, ranges = mirs_true) # <- Overlaps se hace despues de reducir o de conocer el tamano del gen source

g3 %>% as_tibble() %>% count(type, sort = T)

IRanges::reduce(g3) # <- NECESITAMOS REDUCIR, PERO CONSERVAR LOS FEATURES:

IRanges::ranges(g3, use.mcols = T) #, use.mcols = T) 

IRanges::resize(g3, width(g3), fix = "center") # %>% as_tibble() %>% count(type, sort = T)

paste_headers <- function(x) { 
  
  x <- unique(x)
  x <- paste(x, sep = '|', collapse = '|') }

x <- IRanges::resize(g3, width(g3), fix = "center") %>%
  as_tibble() %>%
  group_by(seqnames, start, end) %>%
  summarise(across(type, .fns = paste_headers), .groups = "drop_last") 

source_annot_df <- GRanges(Rle(x$seqnames), 
  ranges =  IRanges(start = x$start, end = x$end), type = x$type)


# GenomicRanges::setdiff(g3, g3)

# %>% as_tibble() %>% count(type, sort = T)

# path <- "/Users/cigom/Documents/MIRNA_HALIOTIS/RNA_SEQ_ANALYSIS/REPEAT_MASKER_OUT"

# mask_f <- list.files(path = path, pattern = "multi_genome.newid.fa.out.gff", full.names = T)

# gr3 <- rtracklayer::import(mask_f)



kpPlotMarkers(kp, data = mirs_true, labels = source_annot_df$type)
# ggbio::autoplot(mirs_true, aes(color = type, fill = type))

# ggplot(mirs_true) + 
#   ggbio::layout_circle(aes(fill = type, y = score), geom = "rect") +
#   theme(legend.position = "top")


# omit ====

kp <- plotKaryotype(genome = g, plot.type=4,  labels.plotter = NULL, chromosomes = chr.names) 

kpPlotDensity(kp, data=mirs_true, r0=0, r1=0.5)

nreps <- length(mirs_true)

pb <- txtProgressBar(min = 0, max = nreps, style = 3, width = 50, char = "=")

for(i in 1:nreps) {
  
  kpPlotRegions(kp, mirs_true[i], r0 = (i-1)*(0.8/nreps), r1 = (i)*(0.8/nreps), 
    col="black", avoid.overlapping = T)
  
  setTxtProgressBar(pb, i)
  
}

# IF ZOOM

# kpPlotDensity(kp, mirs_true, window.size = 0.5e6, data.panel="ideogram", col="#3388FF", border="#3388FF", r0=0.5, r1=0)

# kpAddBaseNumbers(kp, tick.dist = 50000, add.units = TRUE)

# kpPlotDensity(kp, mirs_true, window.size = 0.5e6, data.panel="ideogram", col="#3388FF", border="#3388FF")

kpPlotRegions(kp, data=mirs_true, 
  data.panel= 2, layer.margin = 0.05, 
  border="red", 
  r0=0.9, r1=1,
  avoid.overlapping = T)
# 
# kpPlotCoverage(kp, mirs_true, ymax = 20, r0=0.8,  r1=1, col="#CCCCFF")
# kpAxis(kp, ymin = 0, ymax= 20, numticks = 2, r0 = 0.8, r1=1)
# 

