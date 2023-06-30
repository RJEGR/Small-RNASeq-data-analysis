# RICARDO GOMEZ-REYES

# MICRORNA CONSERVATION


rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)


wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/MIRGENEDB_20230314/"

dim(MIRGENEDB <- read_tsv(paste0(wd, "/MIRGENEDB_2.1.tsv")))

MIRGENEDB <- MIRGENEDB %>% select_at(vars(contains(c("MirGeneDB_ID", "Family","Node_of_origin_", "Seed"))))

view(MIRGENEDB)


# 1) ====
# COUNT THE NUMBER OF HOMOLOGY-BASED:SHORTSTACKS PREDICTED MIRS:

wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out"

head(DB <- read_tsv(paste0(wd, "/RNA_LOCATION_DB.tsv")))

nrow(DB <- DB %>% filter(SRNAtype == "miR"))


DB <- DB %>% 
  drop_na(KnownRNAs) %>%
  mutate(KnownRNAs = strsplit(KnownRNAs, ";")) %>%
  unnest(KnownRNAs) %>%
  select(Name, KnownRNAs, DicerCall, Reads) %>%
  mutate(sep = KnownRNAs) %>%
  separate(sep, into = c("MirGeneDB_ID", "arm"), sep = "_") %>% 
  mutate(sp = sapply(strsplit(KnownRNAs, "-"), `[`, 1)) %>%
  filter(sp != "piR") %>% # <- ALSO EXCLUDE OVERLAPED MIRS/PIRS
  left_join(MIRGENEDB, by = "MirGeneDB_ID")

DB %>% distinct(sp, Family)

DB %>% filter(`Node_of_origin_(family)` == "Mollusca") %>% view()

DB %>% count(arm, `Node_of_origin_(family)`, sort = T)
# DB %>% count(`Node_of_origin_(locus)`) %>% View()

# DB %>% filter(is.na(`Node_of_origin_(family)`))

DB %>%
  # filter(`Node_of_origin_(family)` == "Mollusca") %>%
  group_by(sp, arm, `Node_of_origin_(family)`) %>%
  summarise(n = n()) %>% # Reads = sum(unique(Reads)), 
  # group_by(`Node_of_origin_(family)`) %>%
  mutate(sp = fct_reorder(sp, n)) %>%
  ggplot(aes(x = sp, y = n, fill = arm)) +
  # facet_grid(~ `Node_of_origin_(family)`, scales = "free", space = "free") +
  geom_col(position = "dodge2")
  

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

