
# MICRORNA CONSERVATION

# /Users/cigom/Documents/MIRNA_HALIOTIS/MIRGENEDB_20230314
# /Users/cigom/Documents/MIRNA_HALIOTIS/PIRNA_DB


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

# or generate its own miR multiple-sequence-alignment:

# 1) as fig 2. from James Tarver paper: 

# 2) AS Sempere, Wheeler and collaborators (35,36). 
# USING Clustal /MAFFT 6.85 /MUSCLE. AND manually refined the alignments with RALEE/Gblocks/ and reconstructed evolutionary trees with standard phylogenetic methods: neighbor-joining (40) and maximum likelihood (41) using DECIPHER. (EX. 10.1093/nar/gkt534)

# 3) using msa inside R

# BiocManager::install("msa")

WD <- "/Users/cigom/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out/"

f <- "KNOWN_AND_NOVEL_MIRS_MajorRNA.fasta"

# f <- "mir.fasta"


f <- list.files(path = WD, pattern = f, full.names = T)

.bioc_packages <- c("Biostrings", "DECIPHER", "phangorn", "tidyverse")

sapply(c(.bioc_packages), require, character.only = TRUE)

seqs <- Biostrings::readRNAStringSet(f)

# seqs <- Biostrings::readDNAStringSet(f)

keep <- !grepl("mature|star", names(seqs))

seqs <- seqs[keep]

# Biostrings::replaceAt(seqs, "U")

seqs <- DECIPHER::RemoveGaps(seqs)

hist(width(seqs))

# Construct phylogenetic tree
# https://bioconductor.org/packages/release/bioc/vignettes/DECIPHER/inst/doc/FindingNonCodingRNAs.pdf

alignment <- DECIPHER::AlignSeqs(DNAStringSet(seqs))

P <- PredictDBN(alignment, type="structures")

BrowseSeqs(alignment, patterns=P )



ggmsa::tidy_msa(alignment) %>%
  as_tibble() %>%
  mutate(facet = ifelse(grepl("Mir|Highly_conserved", name), "Known", "Novel")) %>%
  filter(facet == "Novel") %>%
  ggplot(aes(x = position, y = name)) +
  geom_tile(aes(fill = character)) + # size = 0.1, width = 0.95, color = "black"
  facet_grid(facet ~., scales = "free_y") +
  geom_text(aes(label = character), vjust = 0.5, hjust = 0.5, size= 2.5, family =  "GillSans") +
  theme_classic(base_size = 12, base_family = "GillSans") +
  see::scale_fill_pizza() 

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

# The phangorn R package is then used to construct a phylogenetic tree. 
# follow lines construct a neighbor-joining tree, 
# and then fit a GTR+G+I (Generalized time-reversible with Gamma rate variation) maximum likelihood tree 
# using the neighbor-joining tree as a starting point

phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")

dm <- dist.ml(phangAlign)

treeNJ <- NJ(dm) # Note, tip order != sequence order

fit = pml(treeNJ, data=phangAlign)

fitGTR <- update(fit, k=4, inv=0.2)

fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE, control = pml.control(trace = 0))

# save(fitGTR,fit, treeNJ, dm,phangAlign, file = paste0(dir, "dna-sequences-fitGTR.RData"))

plot(fitGTR$tree)

detach("package:phangorn", unload=TRUE)

library(tidyverse)

tree <- fitGTR$tree

x <- ggtree::as_data_frame(fitGTR$tree)

# JOIN TO ...AND THEN:
x %>% tidytree::as.treedata()
