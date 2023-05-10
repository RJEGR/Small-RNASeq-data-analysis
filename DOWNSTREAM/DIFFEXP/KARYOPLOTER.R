


# Using Karyotyper to viz srna data from shortstacks



rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(GenomicRanges)
library(tidyverse)

path <- '/Users/cigom/Documents/MIRNA_HALIOTIS/ENSEMBLE/'
# 
pattern <- "multi_genome.newid.gff3$" # <- to get sequence-region per chromosome

fgff <- list.files(path = path, pattern = pattern, full.names = T)

pattern <- "multi_genome.newid.gtf"

fgtf <- list.files(path = path, pattern = pattern, full.names = T)

path <- '/Users/cigom/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out'

srna_f <- list.files(path = path, pattern = "Results.gff3", full.names = T)

path <- "/Users/cigom/Documents/MIRNA_HALIOTIS/RNA_SEQ_ANALYSIS/REPEAT_MASKER_OUT"

mask_f <- list.files(path = path, pattern = "multi_genome.newid.fa.out.gff", full.names = T)

# LOAD DATA ==== 

gr <- rtracklayer::import(fgff)

gr <- gr[gr$type == "region"]

gr2 <- rtracklayer::import(srna_f)

gr3 <- rtracklayer::import(mask_f)


# (pending) upgrade cytoband information is available from gtf
# https://bernatgel.github.io/karyoploter_tutorial//Tutorial/CustomGenomes/CustomGenomes.html

# cytoband <- rtracklayer::import(fgtf)

# PREPARE LOCUS TO VIZ

# gr2 %>% as_tibble() %>% head() %>% view()

viz <- gr2 %>% as_tibble() %>% 
  # filter(seqnames %in% query.locus[1:105]) %>%
  # filter(MIRNA == "Y") %>% # SWITCH OR NOT
  group_by(seqnames) %>%
  summarise(n = n(), TotalReads = sum(score)) %>% 
  arrange(desc(n))

# viz %>% view()

viz %>% ggplot(aes(n, TotalReads, color = MIRNA)) + geom_point() + scale_y_log10() + scale_x_log10() 

query.locus <- viz %>% distinct(seqnames) %>% pull()

str(query.locus <- unfactor(query.locus))

length(genome <- gr[seqnames(gr) %in% query.locus])

seqlevels(genome) <- unique(unfactor(seqnames(genome)@values))

chr.names <- as.character(GenomeInfoDb::seqnames(genome))

any(duplicated(chr.names))

library(karyoploteR)

# https://bernatgel.github.io/karyoploter_tutorial//Tutorial/PlotTypes/PlotTypes.html



# data <- gr2[which(gr2$MIRNA == "Y")]

# Example 1 (global viz)

data <- gr2

# data %>% as_tibble() %>% filter(score > 10)

query.locus <- unique(seqlevels(data))

length(genome <- gr[seqnames(gr) %in% query.locus])

seqlevels(genome) <- unique(unfactor(seqnames(genome)@values))

chr.names <- as.character(GenomeInfoDb::seqnames(genome))

any(duplicated(chr.names))

kp <- plotKaryotype(genome = genome, plot.type=3,  labels.plotter = NULL, chromosomes = query.locus) 

kpPlotDensity(kp, data=data, r0=0, r1=0.5)

# kpPlotRegions(kp, data=data, r0=0, r1=0.5, data.panel = 2)

y <- log10(data$score)
ymax <- ceiling(max(abs(range(y))))
ymin <- ceiling(min(abs(range(y))))

kpPoints(kp, data, y = log10(data$score), data.panel = 2, ymin = ymin, ymax = ymax)

kpAxis(kp, data.panel=2, ymin = ymin, ymax = ymax)

kpAddLabels(kp, labels = "Log10(ab)", srt=90,  data.panel = 2, srt=90, pos=1, label.margin = 0.05)

# Example 2: zooming

# range(data)

ldat <- split(data, seqnames(data))

zoom <- range(ldat[[1]])

kp <- plotKaryotype(genome = genome, plot.type=3, zoom = zoom[1]) 

kpPlotDensity(kp, data=data)

y <- log10(ldat[[1]]$score)

ymax <- ceiling(max(abs(range(y))))
ymin <- ceiling(min(abs(range(y))))

kpPoints(kp, data, y = log10(data$score), data.panel = 2, ymin = ymin, ymax = ymax)

kpPlotRegions(kp, data=data, data.panel=1, avoid.overlapping = T,
  r0=0, r1=0.5, num.layers = 10, clipping = F)


# Example 3 (make stacked version of clusters)

nreps <- 100 

# ldat[[1]]

chr <- "JALGQA010000616.1" # "JALGQA010000014.1"

# zoom <- range(data[1:nreps])
dat <- ldat[[which(seqlevels(ldat) == chr)]]

zoom <- range(dat)

nreps <- 100 # length(dat)

kp <- plotKaryotype(genome = genome, plot.type=1, zoom = zoom[1]) 

pb <- txtProgressBar(min = 0, max = nreps, style = 3, width = 50, char = "=")

for(i in 1:nreps) {
  
  kpPlotRegions(kp, dat[i], r0 = (i-1)*(0.8/nreps), r1 = (i)*(0.8/nreps), 
    col="black", avoid.overlapping = T)
  
  setTxtProgressBar(pb, i)
  
}

close(pb)

kpPlotRegions(kp, data=dat, 
  data.panel= 2, layer.margin = 0.05, 
  border="red", 
  r0=0.9, r1=1,
  avoid.overlapping = T)

#Example 2: Do the same with a single bigger set of possibly overlapping regions

# kpPlotRegions(kp, dat, r0 = 0, r1 = 0.8, col="#AAAAAA")

kp <- plotKaryotype(genome = genome, plot.type=1, zoom = zoom[1]) 

kpPlotCoverage(kp, dat, ymax = 20, r0=0.8,  r1=1, col="#CCCCFF")
kpAxis(kp, ymin = 0, ymax= 20, numticks = 2, r0 = 0.8, r1=1)


# Example 4 (using bam as input)
# https://bernatgel.github.io/karyoploter_tutorial//Tutorial/PlotBAMCoverage/PlotBAMCoverage.html

path <- '/Users/cigom/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out/'

bam_f <- list.files(pattern = ".subset.bam$", path = path, full.names = T)

kp <- plotKaryotype(genome = genome,  plot.type=4,  labels.plotter = NULL) 

# kpAddBaseNumbers(kp, tick.dist = 25000, add.units = TRUE)

n <- length(bam_f)

pb <- txtProgressBar(min = 0, max = n, style = 3, width = 50, char = "=")

col <- c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3")
cols <- rep(col, each = 3)

for(i in 1:n) {
 
  kpPlotBAMDensity(kp, data=bam_f[i], col=cols[i], border=NA, r0 = (i-1)*(0.8/n), r1 = (i)*(0.8/n), normalize = F)
  # kpPlotBAMCoverage
  kpAxis(kp, ymax=kp$latest.plot$computed.values$max.coverage, r0 = (i-1)*(0.8/n), r1 = (i)*(0.8/n))
  
  kpAddLabels(kp, gsub(".clean.newid.subset.bam", "", basename(bam_f[i])), 
    r0 = (i-1)*(0.8/n), r1 = (i)*(0.8/n), label.margin = 0.03)
  
  setTxtProgressBar(pb, i)
  
}

add_legend("top", legend=c("110-Low", "110-Control", "24-Low", "24-Control"), 
  pch=12, 
  col=col,
  horiz=TRUE, bty='n', cex=0.7)

