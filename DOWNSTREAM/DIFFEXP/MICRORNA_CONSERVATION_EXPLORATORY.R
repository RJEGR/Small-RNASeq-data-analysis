# RICARDO GOMEZ-REYES

# MICRORNA CONSERVATION
# Identificación de microRNAs conocidos en la etapa trocófora (pre-competente) y competente (veliger tardía) en el abulón rojo Haliotis rufescens a través de técnicas de secuenciación masiva.


rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

# 2) ==== 
# TREE-BASED FAMILY ANALYSIS OF DE NOVO/KNOWN MIRS:
# /Users/cigom/Documents/MIRNA_HALIOTIS/MIRGENEDB_20230314
# /Users/cigom/Documents/MIRNA_HALIOTIS/PIRNA_DB

# Generate multiple-sequence-alignment:

# 1) as fig 2. from James Tarver paper: 

# 2) AS Sempere, Wheeler and collaborators (35,36). 
# USING Clustal /MAFFT 6.85 /MUSCLE. AND manually refined the alignments with RALEE/Gblocks/ and reconstructed evolutionary trees with standard phylogenetic methods: neighbor-joining (40) and maximum likelihood (41) using DECIPHER. (EX. 10.1093/nar/gkt534)

# 3) using msa inside R

# BiocManager::install("msa")

WD <- "/Users/cigom/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out/"

f <- "KNOWN_AND_NOVEL_MIRS_MajorRNA.fasta"

# f <- "KNOWN_AND_NOVEL_MIRS_MajorRNA.aln$"

f <- list.files(path = WD, pattern = f, full.names = T)

.bioc_packages <- c("Biostrings", "DECIPHER", "phangorn", "tidyverse")

sapply(c(.bioc_packages), require, character.only = TRUE)

print(seqs <- Biostrings::readRNAStringSet(f))

names(seqs) <- sapply(strsplit(names(seqs), " "), `[`, 1)

# seqs <- DECIPHER::RemoveGaps(seqs)

# hist(width(seqs))

# Construct phylogenetic tree
# https://bioconductor.org/packages/release/bioc/vignettes/DECIPHER/inst/doc/FindingNonCodingRNAs.pdf

# use RNAStringSet or DNAStringSet to turn RNA/DNA and viseversa:
#
# fasta=KNOWN_AND_NOVEL_MIRS_MajorRNA.fasta
# mafft --thread 2 $fasta > ${fasta%.*}.aln
# mafft --globalpair --thread 2 $fasta > ${fasta%.*}.aln2
  
# alignment <- DECIPHER::AlignSeqs(RNAStringSet(seqs))

alignment <- seqs

# P <- PredictDBN(alignment, type="structures")
# BrowseSeqs(alignment, patterns=P )

out <- ggmsa::tidy_msa(alignment) %>% as_tibble() %>% rename("Name"="name") 
  # separate(name, into = c("Name", "KnownRNA"), sep = " ")

# see::social_colors()
# scales::show_col(see::social_colors())

ps <- out %>%
  ggplot(aes(x = position, y = Name)) +
  geom_tile(aes(fill = character)) + # size = 0.1, width = 0.95, color = "black"
  # geom_text(aes(label = character), vjust = 0.5, hjust = 0.5, size= 2.5, family =  "GillSans") +
  theme_classic(base_size = 7, base_family = "GillSans") +
  scale_fill_manual("", values = c("white", "#cd201f", "#FFFC00","#00b489","#31759b")) +
  theme(
    legend.position = "none", 
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.text.y = element_text(hjust = 1),
    axis.ticks.length = unit(5, "pt"))
  

# ggsave(ps, filename = 'ALIGMENT_TREE.png', 
  # path = wd, width = 6, height = 12, device = png, dpi = 300)

# ADD LEFT-TREE AND TOP MOTIF ANALYSIS:

# EX:
# (doi:10.1093/gbe/evy096 ) This alignment was run under a GTR+G model in PhyloBayes.

# The phangorn R package is then used to construct a phylogenetic tree. 
# follow lines construct a neighbor-joining tree, 
# and then fit a GTR+G+I (Generalized time-reversible with Gamma rate variation) maximum likelihood tree 
# using the neighbor-joining tree as a starting point

phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")

dm <- phangorn::dist.ml(phangAlign)

treeNJ <- phangorn::NJ(dm) # Note, tip order != sequence order

# plot(treeNJ)

fit <- phangorn::pml(treeNJ, data=phangAlign)

fitGTR <- stats::update(fit, k = 4, inv = 0.2)

fitGTR <- optim.pml(fitGTR)

# fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE, 
#   control = pml.control(trace = 0))

fitGTR$tree
class(fitGTR$tree)

is.ultrametric(fitGTR$tree)

as.hclust(fitGTR$tree)
?as.hclust.phylo

stats::as.dendrogram(fitGTR$tree)
ps + ggh4x::scale_y_dendrogram(hclust = fitGTR)

# 

topdf <- out %>% 
  group_by(position) %>% 
  dplyr::count(character) %>%
  mutate(Freq = n/sum(n))

topdf %>% tally(Freq)

top <- topdf %>%
  ggplot(aes(y = Freq, x = position, fill = character)) +
  geom_col() +
  labs(x = "") +
  theme_classic(base_size = 7, base_family = "GillSans") +
  scale_fill_manual("", values = c("white", "#cd201f", "#FFFC00","#00b489","#31759b")) +
  theme(
    legend.position = "top", 
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.text.y = element_text(hjust = 1))


library(patchwork)

top/ps + plot_layout(widths = c(2, 2), heights = c(2,10))


ggsave(p, filename = "XXX", 
  path = wd, width = 3.5, height = 10, device = png, dpi = 300)


# OMIT =====
detach("package:phangorn", unload=TRUE)

tree <- fitGTR$tree

x <- ggtree::as_data_frame(tree)

# JOIN TO ...AND THEN:
x %>% tidytree::as.treedata()

# OMIT ====

# CONTRAST AGAINST A PRECURSOR: PROFILE (NOT) :

f <- "mir.fa"

f <- list.files(path = WD, pattern = f, full.names = T)

seqs2 <- Biostrings::readRNAStringSet(f)

keep <- !grepl("mature|star", names(seqs2))

seqs2 <- seqs2[keep]

# perform the alignment
# aligned <- AlignProfiles(seqs, seqs2)


# 
# dots <- matrix(0, width(alignment)[1], width(alignment)[1])
# 
# evidence <- PredictDBN(alignment, type="evidence", threshold=0, verbose=FALSE)
# 
# dots[evidence[, 1:2]] <- evidence[, 3]
# 
# dots[P[, 2:1]] <- 1
# 
# heatmap(dots)



# clade-specific mir family (mirtrace):

# https://raw.githubusercontent.com/friedlanderlab/mirtrace/87a50b8813f2e663e2e3ce5cdf6ca73550317083/src/lib/curated/clade-specific_miRNA_families_of_animal_clades.txt

# (https://doi.org/10.1186/s13059-018-1588-9) According to the previously curated clade-specific miRNA family numbers (Additional file 1), the reference catalog of clade-specific miRNA sequences is obtained from miRBase v21 by selecting the mature miRNA sequences with an ID that contains the family numbers. For example, the primate-specific miRNA family 580 yields five sequences with miRBase v21 IDs hsa-miR-580-5p, hsa-miR-580-3p, mml-miR-580, ptr-miR-580, and ppy-miR-580. A read is identified as clade-specific miRNA if its first 20 nt have an exact match to a reference sequence.

# The mature sequences of the qualified miRNA precursors were used as reference sequences. We mapped the reads of each sample to the collected species-specific miRNA sequences.

# SEE MIRGENE 
# TREEMirGeneDB version 2.0 (http://mirgenedb.org), which now contains high-quality annotations of 10 899 bona fide and consistently named miRNAs constituting 1275 miRNA families from 45 species, representing every major metazoan group, including many well-established and emerging invertebrate and vertebrate model organisms 
f <-  "https://raw.githubusercontent.com/sinanugur/MirMachine/19fab398b12e3cae04c8c51949852609ad39bbf8/mirmachine/meta/tree.newick"

plot(ape::read.tree(f))


# OMIT=====
# SEARCH START NUCLEOTIDE PREFERENCY IN PIRS AND MIRS (FROM SHORTSTACKS)

# SEARCH MOTIF CONSERVATION USING MEME OR OTHER R PACKAGE FOR MIRS
#  https://omarwagih.github.io/ggseqlogo/#installation
devtools::install_github("omarwagih/ggseqlogo")

# or http://yulab-smu.top/ggmsa/articles/ggmsa.html (maybe or not)
devtools::install_github("YuLab-SMU/ggmsa")
miRNA_sequences <- system.file("extdata", "seedSample.fa", package = "ggmsa")

miRNA_sequences <- Biostrings::readRNAStringSet(miRNA_sequences)

library(ggmsa)

ggmsa(miRNA_sequences,  seq_name = TRUE, color = "Clustal", font = "DroidSansMono") + 
  geom_seqlogo(color = "Chemistry_AA")

# tidy_msa(msa, start = start, end = end)

