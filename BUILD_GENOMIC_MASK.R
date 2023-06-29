

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)


# TEST https://github.com/MikeAxtell/sRNA_Viewer

# Non redundant (nr) genomic annotation mask may be build to reduce the Intra-range features (i.e overlapped annotations) and assign mapped sequences to unique/single annotations. Using proference (or hierarchical) annotation selection as follow: rRNA > tRNA, Transposable elements (TE) > protein-coding exon, other ncRNAs, introns, pseudogenes ... (Cei Abreu,2023)



# LA VERSION DE ENSEMBLE CONTIENE INFORMACION SOBRE LAS REGIONES UTR, LO QUE LA VERSION DE NCBI NO TIENE:
# (from ensemble) xgHalRufe1.0.p region 1 94,228,061
# (from ncbi) xgHalRufe1.0.p region 1 94,228,061

library(GenomicRanges)

path <- "~/Documents/MIRNA_HALIOTIS/ENSEMBLE/"

gff_f <- list.files(path, pattern = ".gff3", full.names = T)

# PRIOR TO READ GFF3 FILES
# Lines beginning with '##' are directives (sometimes called pragmas or meta-data) and provide meta-information about the document as a whole

ensemble_gff <- read_tsv(gff_f, comment = "#", col_names = F)

names(ensemble_gff) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")


ensemble_gff %>% count(feature) 
ensemble_gff %>% count(source)
# ensemble_gff <- ape::read.gff(gff_f)

# Note -----
#Inter-range features: 
# input
#         (--------)
#              (----------------)
#                      (-----)
#                                     (---)
#                                     
# range   (-------------------------------)
# 
# reduce  (----------------------)    (---)
# 
# disjoin (---)(---)(-)(-----)(--)    (---)
# 
# GenomicRanges::setdiff(range(input),input)

# Intra-range features:
# input
#                                 (--)
#                         (----)
#                         .    .
# resize                  (--------)
# resize(fix="end")   (--------)

# Continue -------

# For this purpose, single sequence-annotation positions was obtained using the `findOverlaps` function
# from `GenomicRanges` R package (Lawrence et al., 2013). In order to avoid
# overlaps, we count according to the overlap of the central nucleotide of each read,
# using the resize function, with parameter fix=”center” from the GenomicRanges
# R package (Tesis de Isaac Martínez Ugalde)

gr <- head(ensemble_gff, 500)

gr <- GRanges(gr); #rm(ensemble_gff)

genome(gr) <- "xgHalRufe1.0.p"

# Intra range transformations: shift(), narrow(), resize(), flank() 
# Range-based set operations: union(), intersect(), setdiff(), punion(), pintersect()
# Finding/counting overlapping ranges: findOverlaps(), countOverlaps()
# ?GenomicRanges::resize(x, width, fix="center")
# GenomicRanges::setdiff

resize(gr,width = 10)

findOverlaps(gr)
hist(countOverlaps(gr))
# GenomicRanges::resize(gr, width, fix="center")

# https://r-crash-course.github.io/16-bioconductor/
plotRanges_gg <- function(ir) {
  bins <- disjointBins(IRanges(start(ir), end(ir) + 1))
  dat <- cbind(as.data.frame(ir), bin = bins)
  ggplot2::ggplot(dat) +
    ggplot2::geom_rect(ggplot2::aes(xmin = start - 0.5, xmax = end + 0.5,
      ymin = bin, ymax = bin + 0.9), fill = "#FFFFFF", colour = "black") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.title.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_line(colour = "black"))
}


plotRanges_gg(ranges(gr))


# BiocManager::install("ggbio")
# https://github.com/lawremi/ggbio

library(ggbio)