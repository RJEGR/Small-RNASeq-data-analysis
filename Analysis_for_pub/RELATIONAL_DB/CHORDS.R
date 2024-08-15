# Ch


library(tidyverse)

rm(list = ls());

if(!is.null(dev.list())) dev.off()

dir <- "~/Documents/MIRNA_HALIOTIS/MIRNA_PUB_2024/"

print(DB <- read_tsv(file.path(dir, "LONGER_RELATIONAL_DB.tsv")))

dir <-  "~/Documents/MIRNA_HALIOTIS/MIRNA_PUB_2024/"

WGCNA <- read_rds(paste0(dir, "WGCNA_MIRS.rds"))


DB <- DB %>% left_join(WGCNA)

# or load SEQUENCES_MERGED_DESEQ_RES_WGCNA.rds

RES <- read_tsv(list.files(path = dir, 
  pattern = "SEQUENCES_MERGED_DESEQ_RES.tsv", full.names = T)) %>%
  left_join(WGCNA) %>%
  mutate(Contrast = ifelse(sign(log2FoldChange) == -1, sampleB, sampleA)) %>%
  mutate(Contrast = paste0(CONTRAST,"_",Contrast))

RES.P <- RES %>% filter( padj < 0.05  & abs(log2FoldChange) > 1)

# 
DEGS <- RES.P %>% 
  filter(CONTRAST == "CONTRAST_C") %>%
  mutate(HPF = ifelse(sign(log2FoldChange) == -1, "110 hpf", "24 hpf")) %>%
  distinct(MajorRNA, HPF, CONTRAST)

DEGS %>% count(HPF, CONTRAST)

# JOIN LOW PH DEGS -----
# 76 MIRS DEGS IN OA

DEGS_D <- RES.P %>% 
  filter(CONTRAST == "CONTRAST_D") %>% 
  mutate(HPF = ifelse(sign(log2FoldChange) == -1, "110 hpf", "24 hpf")) %>%
  distinct(MajorRNA, HPF, CONTRAST)


mat <- rbind(DEGS, DEGS_D) %>%
  group_by(MajorRNA, HPF) %>%
  summarise(across(CONTRAST, .fns = which_contrast), n = n()) %>%
  left_join(DB) %>% drop_na(COG_name) %>%
  left_join(WGCNA) %>%
  filter(CONTRAST == "CONTRAST_D") %>% # CONTRAST == BOTH ""
  filter(HPF == "110 hpf") %>%
  with(., table(COG_category, WGCNA))

# mat <- DB %>%
#   drop_na(COG_name) %>%
#   with(., table(COG_category, WGCNA))


library(circlize)

co_ocurrance <- rowSums(mat >= 1)
co_ocurrance <- names(co_ocurrance)[co_ocurrance > 1]

co_module <- rowSums(mat)

co_module <- co_module[names(co_module) %in% co_ocurrance]

print(keep <- co_module[co_module/sum(co_module) > 0.1])

mat <- data.frame(mat)

str(col_mat <- data.frame(mat, stringsAsFactors = F))

# table(col_mat$Freq)

# print(threshold <- sum(col_mat$Freq)*0.05)

# sum(keep <- col_mat$Freq >= threshold)


col_mat$Freq <- NULL

col_mat$cols <- as.character(col_mat$WGCNA)

keep <- col_mat$COG_category %in% names(keep)

col_mat$cols[!keep] <- "#00000000"

# col_mat$cols[!keep] <- "#00000000"

grid.col <- unique(as.character(col_mat$WGCNA))
grid.col <- structure(grid.col, names = grid.col)

COG_col <- unique(as.character(col_mat$COG_category))
COG_col <- structure(rep("gray", length(unique(COG_col))), names= COG_col)


dim(col_mat) == dim(mat)  # to make sure it is a matrix

filename <-paste0(dir, "/Chord_110_Contrast_D.png")
# filename <-paste0(dir, "/Chord_110.png")

width= 1500;height=1500;res = 300
# width= 2000;height=2000;res = 300
# 
png(filename, width = width, height = height, res = res)

circos.clear()

circos.par(start.degree = 0, 
  gap.degree = 4, 
  track.margin = c(-0.01, 0.01), 
  points.overflow.warning = FALSE
  # gap.after = c(rep(5, nrow(mat)-1), 15, rep(5, ncol(mat)-1), 15)
)



mat %>%
  chordDiagram(
    order = c(levels(mat$WGCNA), levels(mat$COG_category)),
    grid.col = c(grid.col, COG_col),
    col = col_mat$cols,
    directional = 1,  # to draw direction from module to COGs set 1
    link.sort = T, link.decreasing = T, transparency = 0,
    diffHeight = mm_h(5), target.prop.height = mm_h(1.5),
    annotationTrack = c("grid"), annotationTrackHeight = c(0.03, 0.01),
    # direction.type = c("arrows", "diffHeight"),
    preAllocateTracks = 1,
    small.gap = 10, big.gap = 15)

circos.track(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  xplot = get.cell.meta.data("xplot")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  
  if(abs(xplot[2] - xplot[1]) < 10) {
    circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise", 
      font = par("font"),
      niceFacing = TRUE, adj = c(0, 0.5), col = "black")
  } else {
    circos.text(mean(xlim), ylim[1], sector.name, facing = "inside", 
      font = par("font"),
      niceFacing = TRUE, adj = c(0.5, 0), col= "black")
  }
}, bg.border = NA)

dev.off()

