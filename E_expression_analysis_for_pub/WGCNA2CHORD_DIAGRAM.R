# RICARDO GOMEZ-REYES
# CHORD DIAGRAM
# USING WGCNA RESULTS AND RNA_LOCATION

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

wd <- "~/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out/"

# RES.P <- read_tsv(paste0(wd, "DESEQ_RES.tsv")) %>% filter( padj < 0.05  & abs(log2FoldChange) > 1)

print(DB <- read_tsv(paste0(wd, "/RNA_LOCATION_DB.tsv")))

DATA <- DB %>% filter(SRNAtype == "miR") %>% count(biotype_best_rank, WGCNA)

other_intra <- c("exon", "UTR")

DATA <- DATA %>%
  mutate(biotype_best_rank = ifelse(biotype_best_rank %in% other_intra, "Intragenic", biotype_best_rank))

MODULE_LOW_PH <- c("green", "black", "yellow", "red")

MODULE_IN_DEV <- c("brown", "pink", "turquoise", "blue", "grey")

grid.col <- structure(c(MODULE_LOW_PH, MODULE_IN_DEV), names = c(MODULE_LOW_PH, MODULE_IN_DEV))

library(circlize)

dev.off()  

circos.clear()

DATA %>%
  chordDiagramFromDataFrame(directional = TRUE, 
    # col = colmat, 
    grid.col = grid.col,
    link.sort = T,
    # annotationTrack = "grid", 
    big.gap = 10, small.gap = 1,
    preAllocateTracks = list(track.height = 0.1),
    link.target.prop = FALSE)

# OMIT: ====

DATA <- RES.P %>%  mutate(SIGN = sign(log2FoldChange)) %>% select(Name, CONTRAST,sampleB, SIGN)


DATA %>% 
  mutate(SIGN = paste0(SIGN, "_", CONTRAST)) %>%
  count(SIGN) %>% pull(SIGN)

DATA <- DATA %>%  
  mutate(sampleB = ifelse(SIGN == 1 & CONTRAST %in% c("CONTRAST_C"), "110 HPF (Ctr)", sampleB)) %>%
  mutate(sampleB = ifelse(SIGN == -1 & CONTRAST %in% c("CONTRAST_D"), "110 HPF (Low)", sampleB)) %>%
  mutate(sampleB = ifelse(SIGN == -1 & CONTRAST %in% c("CONTRAST_A"), "pH 7.6", sampleB)) %>%
  mutate(sampleB = ifelse(SIGN ==  1 & CONTRAST %in% c("CONTRAST_A", "CONTRAST_B"), "pH 8.0", sampleB)) %>%
  distinct(Name, sampleB)


DATA <- DATA %>% left_join(select(DB, Name, biotype_best_rank, WGCNA), by = "Name") 



DATA %>%
  dplyr::count(biotype_best_rank,sampleB, WGCNA,sort = T) %>%
  chordDiagramFromDataFrame(directional = TRUE, 
    # annotationTrack = "grid", 
    big.gap = 10, small.gap = 1,
    preAllocateTracks = list(track.height = 0.1),
    link.target.prop = FALSE)

# PREP MATRIX

# DATA <- RES.P %>% 
#   mutate(SIGN = sign(log2FoldChange)) %>% 
#   select(Name, Family, CONTRAST, SIGN) %>% 
#   pivot_wider(names_from = CONTRAST, values_from = SIGN, values_fill = 0)

DATA %>%
  with(., table(Name, CONTRAST_A)) 
  # chordDiagramFromMatrix()

recode_fc <- structure(c("pH 7.6","pH 8.0"), names = c(1,-1))

recode_to <-  structure(c("24 HPF", "110 HPF"), names = c("CONTRAST_A", "CONTRAST_B"))
# 
# RES.P %>% 
#   filter(CONTRAST %in% c("CONTRAST_C", "CONTRAST_D")) %>%
#   mutate(SIGN = sign(log2FoldChange)) %>% 
#   # select(Name, Family, CONTRAST, SIGN) %>%
#   dplyr::mutate(CONTRAST = dplyr::recode_factor(CONTRAST, !!!recode_to)) %>%
#   dplyr::mutate(SIGN = dplyr::recode_factor(SIGN, !!!recode_fc)) %>%
#   with(., table(SIGN, CONTRAST)) %>%
#   chordDiagramFromMatrix()


RES.P %>%  
  mutate(SIGN = sign(log2FoldChange)) %>% 
  dplyr::count(CONTRAST,sampleB, SIGN) %>%
  chordDiagramFromDataFrame()
  # select(Name, Family, CONTRAST, SIGN) %>%
  # mutate(SIGN = ifelse(CONTRAST %in% c("CONTRAST_C", "CONTRAST_D") & SIGN == -1, ""))
  with(., table(Name, sampleB))

EXCLUSIVE_MIRS %>%  
  with(., table(SIGN, CONTRAST)) %>%
  chordDiagramFromMatrix()


library(circlize)

dev.off()  

circos.clear()

chordDiagram(mat, col = colmat, grid.col = state_col2,
  directional = TRUE, annotationTrack = "grid", 
  big.gap = 10, small.gap = 1,
  preAllocateTracks = list(track.height = 0.1),
  link.target.prop = FALSE)