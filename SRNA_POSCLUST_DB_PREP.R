
# RICARDO GOMEZ-REYES
# PART OF STEP: SRNA POS-CLUSTERING AND NOVEL CRITERIA TO IDENTIFY FUNCTIONAL SIRNAS/PIRNAS:

# PREPARE INFO REGARGIND:
# 1) SAMPLE FREQUENCY (I.E. PREVALENCE, SPREAD OR CO-OCURRANCE OF SRNAS IN THE SAMPLES)
# 2) SRNA ABUNDANCES (ROWSUM OF SRNAS. IT WOULD OR WOULDNT BE CORRELATED TO SAMPLE FREQ.)
# 3) SRNA-SEQUENCE GROUPING (I.E. BASED ON SEQ. SIMILARITY GROUP IF SEQS > 70 % Sequence Identity)

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

path <- "~/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out/"

out <- read_rds(paste0(path, "Counts_and_mtd.rds"))

dim(COUNTS <- out[[1]])

dim(MTD <- out[[2]])

head(.DB <- read_tsv(paste0(path, "/RNA_LOCATION_DB.tsv")))

MIRGENEDB <- read_tsv(paste0(path, "SRNA2MIRGENEDB.tsv")) %>%
  tidyr::unite("Family", c("Family", "arm"), sep = "-")

# head(sort_m <- sort(rowSums(COUNTS), decreasing = T))

# any(names(sort_m) %in% rownames(COUNTS))

# head(COUNTS[match(rownames(COUNTS), names(sort_m) ),])

# ORDER SRNA MATRIX BY DECREASING COUNT (Freq. and total read count):

# head(SORTED_MAT <- COUNTS[order(rowSums(COUNTS),decreasing=T),])

.GCOUNTS <- COUNTS %>%
  as_tibble(rownames = "Name") %>%
  right_join(.DB %>% select(Name, MajorRNA)) %>%
  group_by(MajorRNA) %>%
  summarise_at(vars(colnames(COUNTS)), sum) 


GCOUNTS <- .GCOUNTS %>% select_if(is.double) %>% as(., "matrix")

rownames(GCOUNTS) <- .GCOUNTS$MajorRNA

head(GCOUNTS <- GCOUNTS[order(rowSums(GCOUNTS),decreasing=T),])

READS_AND_FREQ <- GCOUNTS %>%
  as_tibble(rownames = "Name") %>%
  pivot_longer(-Name, values_to = "Reads") %>%
  mutate(REP = substr(name, 1,nchar(name)-1)) %>%
  filter(Reads > 0) %>%
  group_by(Name, REP) %>%
  summarise(Reads = sum(Reads), Freq = n()) %>%
  dplyr::rename("MajorRNA" = "Name")

# READS_AND_FREQ %>% group_by(REP) %>% rstatix::cor_test(Reads, Freq)

# READS_AND_FREQ %>% select(-Freq) %>% pivot_wider(names_from = REP, values_from = Reads, values_fill = 0)

# READS_AND_FREQ %>% select(-Reads) %>% pivot_wider(names_from = REP, values_from = Freq, values_fill = 0)

# ADD SAMPLES FREQ. ====

DB <- READS_AND_FREQ %>% 
  summarise(Freq = sum(Freq)) %>% # Reads = sum(Reads), 
  left_join(.DB) %>% 
  dplyr::rename("SampleFreq" = "Freq") %>%
  arrange(match(MajorRNA, rownames(GCOUNTS)))
  # arrange(match(Name, rownames(SORTED_MAT)))

# UPGRADE DB

which_cols <- names(DB)[!names(DB) %in% names(.DB)]

DB <- DB %>% dplyr::select(any_of(c(names(.DB), which_cols)))

write_tsv(DB, paste0(path, "/RNA_LOCATION_DB.tsv"))

DB  %>% 
  filter(SampleFreq == 1) %>%
  mutate(SampleFreq = SampleFreq/12) %>%
  ggplot(aes(x = Reads, y = SampleFreq, color = SRNAtype)) + 
  # geom_point()
  stat_ecdf()

# Sanity check:

sum(.DB$Reads)

sum(DB$Reads) # 83 808 529 <- consistente

DB  %>% 
  group_by(SampleFreq, SRNAtype) %>%
  summarise(Frac = sum(Reads)/sum(DB$Reads)) %>%
  ggplot(aes(x = as.factor(SampleFreq), y = Frac, fill = SRNAtype)) + 
  geom_col() +
  facet_grid(~ SRNAtype, scales = "free", space = "free") +
  scale_y_continuous(labels = scales::percent) 

DB  %>% 
  # filter(SampleFreq >=2) %>%
  group_by(SampleFreq, SRNAtype, DicerCall) %>%
  summarise(Frac = sum(Reads)/sum(DB$Reads)) %>%
  ggplot(aes(x = as.factor(DicerCall), y = Frac, fill = SRNAtype)) + 
  geom_col() +
  facet_grid(SampleFreq ~ SRNAtype, scales = "free", space = "free") +
  scale_y_continuous(labels = scales::percent) 

# KEEP W/ SAMPLE READS > 100 & SAMPLE FREQ. 

# VALIDATE MIRS PREVALENCE == 12:

str(query.rna <- DB %>% filter(SRNAtype == "miR") %>% pull(MajorRNA))

GCOUNTS %>% 
  as_tibble(rownames = "MajorRNA") %>%
  filter(MajorRNA %in% query.rna)

#  ggplot(aes(x = LIBRARY_ID, y = name, fill = cor))

READS_AND_FREQ %>% 
  ungroup() %>%
  filter(MajorRNA %in% query.rna) %>%
  # filter(Freq == 2) %>%
  ggplot(aes(x = REP, y = MajorRNA, fill = as.factor(Freq))) +
  geom_raster() +
  theme_bw() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.background = element_rect(fill = 'white', color = 'white'),
    panel.border = element_blank(),
    panel.grid.major = element_blank())

str(query.ids <- DB %>% filter(SRNAtype == "miR") %>% pull(Name))

any(query.ids %in% rownames(COUNTS))

dim(M <- COUNTS[match(query.ids, rownames(COUNTS)),])

heatmap(edgeR::cpm(M))


# # Generate multiple-sequence-alignment:
# SUBSEQ BY FOLLOW CRITERIO:

# USING MSA (in R)


# 1) as fig 2. from James Tarver paper: 

# 2) AS Sempere, Wheeler and collaborators (35,36). 
# USING Clustal /MAFFT 6.85 /MUSCLE. AND manually refined the alignments with RALEE/Gblocks/ and reconstructed evolutionary trees with standard phylogenetic methods: neighbor-joining (40) and maximum likelihood (41) using DECIPHER. (EX. 10.1093/nar/gkt534)

.DB %>% dplyr::count(SRNAtype)

DB <- .DB %>% filter(SRNAtype == "miR")

str(seqs <- DB %>% arrange(MajorRNA) %>% pull(MajorRNA)) # %>% head(100)

str(seqs.name <- DB %>% arrange(MajorRNA) %>% select(Name) %>% pull(Name))

seqs <- bioseq::rna(seqs)

# bioseq::seq_consensus(bioseq::rna(seqs))

library(msa)
library(bioseq)

# method=c("ClustalW", "ClustalOmega", "Muscle")

align <- msa::msa(RNAStringSet(seqs), method = "Muscle")

.align <- msaConvert(align)$seq

print(bioseq::rna(.align))

# PREP SRNA-SEQUENCE POS-CLUSTERING: 

table(SEQ_CLUSTERS <- bioseq::seq_cluster(bioseq::rna(.align), threshold = 0.01, method = "single"))


print(bioseq::rna(.align)[SEQ_CLUSTERS])

msa <- RNAStringSet(.align)

names(msa) <- seqs.name

out <- ggmsa::tidy_msa(msa) %>% as_tibble() %>% dplyr::rename("Name"="name") %>%
  left_join(MIRGENEDB %>% distinct(Name, Family)) %>% 
  mutate(Family = ifelse(is.na(Family), Name, Family))


x_bin <- ape::as.DNAbin(as(DNAStringSet(msa), "matrix"))
x_bin <- as.matrix(x_bin)
x_dist <- ape::dist.dna(x_bin, model = "raw")

hc_seqs <- stats::hclust(x_dist, method = "complete")

seqs_order <- hc_seqs$labels[hc_seqs$order]

recode_to <- out %>% distinct(Name, Family)

recode_to <- structure(recode_to$Family, names = recode_to$Name)

identical(sort(names(recode_to)),sort(seqs_order))

seqs_order <- recode_to[match(seqs_order, names(recode_to))]

identical(names(seqs_order),  hc_seqs$labels[hc_seqs$order])


ps <- out %>%
  # filter(grepl("^MIR-10", Family)) %>%
  mutate(character = ifelse(character == "-", NA, character)) %>%
  drop_na(character) %>%
  ggplot(aes(x = position, y = Name)) +
  geom_tile(aes(fill = character), color = 'white', linewidth = 0.2) +
  # geom_text(aes(label = character), vjust = 0.5, hjust = 0.5, size= 2.5, family =  "GillSans") +
  theme_classic(base_size = 7, base_family = "GillSans") +
  scale_x_continuous(breaks = seq(1, 41, by = 2)) +
  ggh4x::scale_y_dendrogram(hclust = hc_seqs, position = "left", labels = NULL) +
  guides(y.sec = guide_axis_manual(labels = seqs_order, label_size = 3.5)) +
  scale_fill_manual("Nuc.", values = c("#cd201f", "#FFFC00","#00b489","#31759b"), 
    na.value = "white") +
  theme(
    legend.position = "none", 
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.text.y = element_text(hjust = 1),
    axis.ticks.length = unit(5, "pt"))

topdf <- out %>% 
  mutate(character = ifelse(character == "-", NA, character)) %>%
  drop_na(character) %>%
  group_by(position) %>% 
  dplyr::count(character) %>%
  mutate(Freq = n/sum(n))

topdf %>% tally(Freq)

top <- topdf %>%
  ggplot(aes(y = Freq, x = position, fill = character)) +
  geom_col() +
  labs(x = "") +
  theme_classic(base_size = 7, base_family = "GillSans") +
  scale_fill_manual("Nuc.", values = c("#cd201f", "#FFFC00","#00b489","#31759b"), 
    na.value = "white") +
  theme(
    legend.position = "top", 
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.text.y = element_text(hjust = 1))


library(patchwork)

p <- top / plot_spacer() / ps + plot_layout(heights = c(1, -0.5, 7), widths = c(2, 2, 2))

# p <- top/ps + plot_layout(widths = c(2, 2), heights = c(2,10))

ggsave(p, filename = "ALIGMENT_TREE.png", 
  path = path, width = 3.5, height = 8, device = png, dpi = 300)


# EX.

library(Biostrings)

s1 <- DNAString("AGTATAGATGATAGAT")
s2 <- DNAString("AGTAGATAGATGGATGATAGATA")

palign1 <- pairwiseAlignment(s1, s2)
palign1
pid(palign1)

# palign2 <-
#   pairwiseAlignment(s1, s2,
#     substitutionMatrix =
#       nucleotideSubstitutionMatrix(match = 2, mismatch = 10, baseOnly = TRUE))
# palign2
# pid(palign2, type = "PID4")


# W/ DECIPHER NOT WELL ALIGMENT:
# print(seqs <- Biostrings::RNAStringSet(seqs))
# Biostrings::pairwiseAlignment(seqs, seqs)
# alignedDNA <- DECIPHER::AlignSeqs(RNAStringSet(seqs))
