
# FORMAT MajorRNA miRNA sequences (Deduplicate sequences and headers)
# FORMAT mirgeneDB (Deduplicate sequences and headers)
# LOAD MajorRNA sequences (118)
# LOAD mirgeneDB sequences (17599)
# 
# RUN alignment
# Estimate specificity of clade-specific miRNa detection as Kang et al., 2018 (mirtrace)
# To estimate how many reads are likely to be ientified as clade-specific miRNAs by chance, a binomial model is adopted. 

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)


.bioc_packages <- c("Biostrings", "DECIPHER", "phangorn", "tidyverse")

sapply(c(.bioc_packages), require, character.only = TRUE)

# 0. DATABASE
wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/MIRGENEDB_20230314/"

dim(MIRGENEDB <- read_tsv(paste0(wd, "/MIRGENEDB_2.1.tsv")))

MIRGENEDB <- MIRGENEDB %>% select_at(vars(contains(c("MirGeneDB_ID", "Family","Node_of_origin_(family)", "Seed"))))

# 1. LOAD SEQUENCES -----
# MAJORRNAS

wd <-  "~/Documents/MIRNA_HALIOTIS/MIRNA_PUB_2024/"

f <- list.files(path = wd, pattern = "DESEQ_RES.tsv", full.names = T)

DESEQDB <- read_tsv(f) %>% distinct(MajorRNA, MirGeneDB_ID) %>% 
  dplyr::rename("seqs" = "MajorRNA") 


DESEQDB <- DESEQDB %>%
  # separate(MirGeneDB_ID, into = c("MirGeneDB_ID", "arm"), sep = "_") %>%
  mutate(Family = "Abalone", 
        `Node_of_origin_(family)` = "Abalone",
       Seed = substr(seqs, 2, 8)) %>%
  filter(grepl("Cluster_", MirGeneDB_ID))

# are 47 novel majorRNA sequences

# Using tip labels from cluster/clade analysis made from all the mirgenedb+majorRNAs seqs.

tip_labels <- read_rds(file.path(wd, "tip_labels.rds"))

# fasta_prep <- read_tsv(f) %>% distinct(MajorRNA, MirGeneDB_ID) 
  
# seqs <- fasta_prep %>% pull(MajorRNA)

# headers <- fasta_prep %>% mutate(MirGeneDB_ID = paste0(">", MirGeneDB_ID)) %>% pull(MirGeneDB_ID) 

# fasta <- c(rbind(headers, seqs))

# write(fasta, file= paste0(wd, "MIRS_MajorRNA.fasta"))

f <- "^MIRS_MajorRNA.fasta$"

f <- list.files(path = wd, pattern = f, full.names = T)

print(seqs <- Biostrings::readRNAStringSet(f))

# names(seqs) <- sapply(strsplit(names(seqs), " "), `[`, 1)

# LOAD sequences from MirgeneDB ----

wd <-  "~/Documents/MIRNA_HALIOTIS/MIRGENEDB_20230314/"

f <- list.files(path = wd, pattern = "ALL-mat.fa", full.names = T)

print(mirgenedb <- Biostrings::readRNAStringSet(f))

# a priori we known abalone mirs identified in this research is clustered together 
# the clade 

mirgenedb <- mirgenedb[mirgenedb %in% tip_labels]


table(sapply(strsplit(names(mirgenedb), "-"), `[`, 1))

mirgenedb <- unique(mirgenedb)

MIRGENEDB <- data.frame(seqs = mirgenedb) %>% 
  as_tibble(rownames = "Name") %>%
  separate(Name, into = c("MirGeneDB_ID", "arm"), sep = "_") %>%
  left_join(MIRGENEDB, by = "MirGeneDB_ID") %>%
  unite("MirGeneDB_ID",MirGeneDB_ID:arm)


MIRGENEDB <- DESEQDB %>%
  select(names(MIRGENEDB)) %>%
  rbind(MIRGENEDB)

# 
# mirgenedb

seqs <- c(seqs, mirgenedb)

print(dup <- table(names(seqs))[table(names(seqs))>1])

seqs[names(seqs) %in% names(dup)]

# Deduplicate ----

seqs <- unique(seqs)

# deduplicat names by using seqs by names

names(seqs) <- seqs

# Test methods and models

# subseqs <- seqs[sample(nrow(seqs), 100)]

# seqs must match w/ updated mirgenedb

# data.frame(seqs = seqs) %>% 
#   as_tibble() %>%
#   left_join(MIRGENEDB, by = "seqs")

# Generate multiple-sequence-alignment: ----
# using ClustalW or muscle, (Time-computer consuming)
# As Sempere, Wheeler , et al. 
# using Clustal /MAFFT 6.85 /MUSCLE. AND 
# manually refined the alignments with RALEE/Gblocks/

# library(msa)
msa_align <- msa::msa(RNAStringSet(seqs), method = "ClustalW")

RNAStringSet(msa_align)

# OR using refined alignment
# The phangorn R package is then used to construct a phylogenetic tree. 
# follow lines construct a neighbor-joining tree, 
# and then fit a GTR+G+I (Generalized time-reversible with Gamma rate variation) maximum likelihood tree 
# using the neighbor-joining tree as a starting point


alignment <- DECIPHER::AlignSeqs(RNAStringSet(seqs))

# This method is time-computer consuming
phangAlign <- phyDat(as(alignment, "matrix"), type="DNA") 
dm <- phangorn::dist.ml(phangAlign)
treeNJ <- phangorn::NJ(dm) # Note, tip order != sequence order
# plot(treeNJ)
fit <- phangorn::pml(treeNJ, data=phangAlign) # this step is fast
fitGTR <- stats::update(fit, k = 4, inv = 0.2) # this step is fast
fitGTR <- optim.pml(fitGTR) # this step is fast

# saveRDS(fitGTR, file.path(wd, "fitGTR.rds")) # this is the full output of mirgenedb

# Instead, use this method
# Get/calculate dendogram group -----
# so fast!!
# alignment_bin <- as.DNAbin(DNAStringSet(alignment))
# distmat <- dist.dna(alignment_bin, as.matrix = TRUE, pairwise.deletion = TRUE)
# 
# clust <- hclust(distmat, method = "single") # <-- run

# clust <- DECIPHER::IdClusters(distmat,
#   method = "UPGMA",
#   cutoff = 0.01,
#   showPlot = TRUE,
#   type = "both",
#   myXStringSet = NULL,
#   # model = MODELS,
#   collapse = 0,
#   processors = NULL,
#   verbose = TRUE)


# Join clades to MIRGENEDB ----

hist(fitGTR$tree$edge.length)

dist_mat <- ape::cophenetic.phylo(fitGTR$tree)

hc_pw_seq_align <- hclust(as.dist(dist_mat), method="single")


# Find special clusters:
library(dynamicTreeCut)

clusters <- cutreeDynamic(hc_pw_seq_align, distM = as.matrix(dist_mat), method = "tree")
# we need to sort them to the order of the dendrogram:
clusters <- clusters[order.dendrogram(as.dendrogram(hc_pw_seq_align))]

clusters_numbers <- unique(clusters) - (0 %in% clusters)
n_clusters <- length(clusters_numbers)
print(n_clusters)
# Screen which clusters/clades are for abalone 

data.frame(seqs = seqs) %>% 
  as_tibble() %>% mutate(Cluster = clusters) %>%
  left_join(MIRGENEDB, by = "seqs") %>%
  mutate(`Node_of_origin_(family)` = ifelse(is.na(`Node_of_origin_(family)`), "New", `Node_of_origin_(family)`)) %>%
  # count(Cluster, `Node_of_origin_(family)`) %>% # cluster 0 and 12 are clade for new
  filter(Cluster %in% c(0, 12)) %>%
  pull(seqs) -> tip_labels


# write_rds(tip_labels, file.path(wd, "tip_labels.rds"))



tree <- fitGTR$tree

# tree_f <- ape::drop.tip(tree, tree$tip.label[!tree$tip.label %in% tip_labels])



hc <- hclust(as.dist(cophenetic.phylo(tree)), method="ward.D2")

library(ggtree)


tree <- as_data_frame(tree) %>%
  as_tibble() %>%
  left_join(MIRGENEDB, by = c("label"="seqs")) %>%
  mutate(isNovel = ifelse(grepl("Cluster_", MirGeneDB_ID), "True", "False")) %>%
  tidytree::as.treedata()



ggtree(tree, layout="circular") + 
  geom_nodepoint(color="blue", alpha=2/4, size=2.5) +
  geom_tiplab(aes(label = `Node_of_origin_(family)`, color = isNovel), 
    offset = 1, family= "GillSans", align = F) 
  # geom_label(aes(x=branch, label=`Node_of_origin_(family)`, fill = isNovel),
  #   color = "white", family = "GillSans", hjust = -.7)


# plot as sequences 

subalign <- RNAStringSet(msa_align) # or alignment

# subalign <- subalign[names(subalign) %in% tip_labels]

out <- ggmsa::tidy_msa(subalign) %>% as_tibble() %>% dplyr::rename("seqs"="name") 

ps <- out %>%
  left_join(MIRGENEDB, by = "seqs") %>%
  ggplot(aes(x = position, y = seqs)) +
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

  #

ps 
  # facet_grid(`Node_of_origin_(family)` ~., scales = "free") 
  # ggh4x::scale_y_dendrogram(hclust = hc, position = 'right', labels = NULL)
  # guides(y.sec = guide_axis_manual(labels = hc_order, label_size = 5, label_family = "GillSans")) 

# out <- ggmsa::tidy_msa(alignment) %>% as_tibble() %>% dplyr::rename("Name"="name") 


# Regulatory tree ----
# Include regulatory function tree and compare tree to tree relationship

# https://www.datanovia.com/en/lessons/comparing-cluster-dendrograms-in-r/
  
dir <- "~/Documents/MIRNA_HALIOTIS/MIRNA_PUB_2024/"
TARGETDB <- read_tsv(file.path(dir, "SEQUENCES_MERGED_MIRNA_TARGET_DB.tsv"))

m <- TARGETDB %>% distinct(MajorRNA, gene_id) %>% mutate(n = 1) %>%
  pivot_wider(names_from = gene_id, values_from = n, values_fill = 0) %>%
  data.frame()

rownames(m) <- m$MajorRNA

m$MajorRNA <- NULL

d <- dist(m, method = "binary")
hc_regulatory_mir <- hclust(d, method="ward.D2")

# plot(hc_regulatory_mir)

tree <- fitGTR$tree

tree_f <- ape::drop.tip(tree, tree$tip.label[!tree$tip.label %in% hc_regulatory_mir$labels]) 




dist(cophenetic.phylo(tree_f), method = "binary")

hc_pw_seq_align <- hclust(as.dist(cophenetic.phylo(tree_f)), method="ward.D2")

library(dendextend)

# Create two dendrograms
dend1 <- as.dendrogram (hc_regulatory_mir)
dend2 <- as.dendrogram (hc_pw_seq_align)

# Create a list to hold dendrograms
dend_list <- dendlist(dend1, dend2)
  


#  alignment quality vals
untangle(method = "ladderize") # 0.6520217
untangle(method = "step1side") #  0.3670336 <- most linear
untangle(method = "random", R = 10) # 0.4658022


# Align and plot two dendrograms side by side
dendlist(dend1, dend2) %>%
  untangle(method = "step1side") %>% # Find the best alignment layout
  tanglegram(
    highlight_distinct_edges = FALSE,
    common_subtrees_color_lines = TRUE,
    common_subtrees_color_branches = TRUE) 


# Compute alignment quality. Lower value = good alignment quality
dendlist(dend1, dend2) %>%
  untangle(method = "step1side") %>%
  entanglement() # Alignment quality

# Baker’s Gamma Index (see baker’s paper from 1974) is a measure of association (similarity) between two trees of Hierarchical clustering (dendrograms). It is defined as the rank correlation between the stages at which pairs of objects combine in each of the two tree

cor_bakers_gamma(dend1, dend2)

# dend_diff(dend1, dend2)

# Entanglement is a measure between 1 (full entanglement) and 0 (no entanglement). T

# 0.3670336

cors <- cor.dendlist(dend_list, method = "common_nodes")
# Print correlation matrix
round(cors, 2)

# Visualize the correlation matrix using corrplot package
library(corrplot)
# corrplot(cors, "pie", "lower")

# PHYLOGENETIC TREE BASED ON KNOWN (MIRBASE) MICRORNAS AND NOVEL MICRORNA IDENTIFIED IN ABALONE
#  using the ‘phytools’ package91 (v1.9.16) and used to perform an Ancestral State Reconstruction (ASR). All the species were coded as LAS or NSA. The ASR was then conducted using the ‘ape’92 package (v5.7-1) and plotted using the ‘phytools’ package (v1.9.16) in R (v4.2.2). Both equal rate (ER) and all rate different (ARD) models were tested, and the all-rate different model was retained based on the log-likelihood value.

# https://peerj.com/articles/16505/#aff-1

# RICARDO GOMEZ-REYES

# MICRORNA CONSERVATION
# Identificación de microRNAs conocidos en la etapa trocófora (pre-competente) y competente (veliger tardía) en el abulón rojo Haliotis rufescens a través de técnicas de secuenciación masiva.


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

# f <- "KNOWN_AND_NOVEL_MIRS_MajorRNA.fasta"

f <- "^MajorRNA.fasta$"

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

# alignment <- seqs

# P <- PredictDBN(alignment, type="structures")
# BrowseSeqs(alignment, patterns=P )

out <- ggmsa::tidy_msa(alignment) %>% as_tibble() %>% dplyr::rename("Name"="name") 
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
#   control = pml.control(trace = F))

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

library(ggtree)

x <- as_data_frame(tree)

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


# (10.1093/bioinformatics/btab708) approach to capture the essential sequence motifs and hairpin loops representing a non-coding RNA famil
# align the sequences

