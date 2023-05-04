# EXPLORING IRANGE

# =====

library(IRanges)

ir1 <- IRanges(start = c(1,4,8), end = c(3,7,10))

ir2 <- IRanges(start = c(3,4), width = 3)

ov <- findOverlaps(ir1, ir2)

ov

disjoin(ir1)

library(gggenomes)

# https://github.com/thackl/gggenomes/
# http://genomicsclass.github.io/book/

devtools::install_github("thackl/thacklr")
devtools::install_github("thackl/gggenomes")




#

# Bioconductor and GitHub packages

# BiocManager::install(c("Homo.sapiens",
#   "GenomicFeatures",
#   "genomicsclass/ERBS",
#   "genomicsclass/ph525x"))
# 
# # load GM12878 and HepG2 objects from ERBS package
# library(ERBS)
# data(GM12878)
# data(HepG2)
# 
# x = HepG2[order(HepG2),]
# seqnames(x)     # demonstrate usefulness of Rle type
# as.character(seqnames(x))

# findOverlaps and %over% ========
# We will demonstrate two commonly used methods for comparing GRanges objects. First we build two sets of ranges:

library(GenomicRanges)

# GRangesList()

gr <- GRanges("chrZ", IRanges(start=c(5,10),end=c(35,45)),
  strand="+", seqlengths=c(chrZ=100L))

gr2 <- GRanges("chrZ",IRanges(11:13,51:53))
gr1 <- GRanges("chrZ",IRanges(c(1,11,21,31,41),width=5),strand="*")

GRangesList(gr, gr2)

fo <- findOverlaps(gr1, gr2)

queryHits(fo) # position vector where coordinates from query is found
subjectHits(fo) # position vector where coordinates from subject is found


gr1 %over% gr2
gr1[gr1 %over% gr2]

# ranges from the query for which we found a hit in the subject

index = queryHits(fo)
erbs = gr1[index,]

# Extract only ranges

granges(erbs)


