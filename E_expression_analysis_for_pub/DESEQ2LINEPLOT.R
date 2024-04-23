# GENERATE AVERAGE-LINE-PLOT USING DIFF EXP. MIRS AS INPUTS
# FOCUS ON MIRS UP-EXPRESSED UNDER LOW PH (i.e CONTRAST A AND B)
# OR ONLY IN EXCLUSIVE MIRS
# COMPARISON FROM CONTRAST C AND D SHOW HOW REGULATORY NETWORK DURING DEVELOPMENT (AT 110 HPF) CHANGES IN REPONSE TO LOW PH

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

wd <- "~/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out/"

dds <- read_rds(paste0(wd, "/DDS_DESEQ2.rds"))

# RES.P <- read_tsv(paste0(wd, "DESEQ_RES_P.tsv")) %>% filter( padj < 0.05) 


RES.P <- read_tsv(paste0(wd, "DESEQ_RES.tsv")) %>% filter( padj < 0.05  & abs(log2FoldChange) > 1) %>%
  mutate(star = ifelse(padj <.001, "***", 
    ifelse(padj <.01, "**",
      ifelse(padj <.05, "*", "")))) 

RES.P %>% dplyr::count(CONTRAST)

RES.P %>% select(Family, log2FoldChange, CONTRAST) %>% 
  group_by(Family, CONTRAST) %>% 
  filter(CONTRAST %in% c("CONTRAST_A", "CONTRAST_B")) %>%
  summarise(log2FoldChange = sum(log2FoldChange)) %>%
  view()


# OMIT
# Only mirs upexpressed under low pH

RES.P <- RES.P %>% filter(log2FoldChange < 0)

RES.P %>% dplyr::count(CONTRAST)

# read_tsv(paste0(wd, "SRNA2MIRGENEDB.tsv"))

CONTRAST <- as_tibble(SummarizedExperiment::colData(dds)) %>% dplyr::select(starts_with("CONTRAST")) %>% names()

recode_to <- structure(c("24 HPF", "110 HPF", "Ctrl pH", "Low pH"), names = CONTRAST)

out <- list()

for (j in 1:length(CONTRAST)) {
  
  i <- j
  
  cat("\nRunning ",CONTRAST[i], "\n")
  
  intgroup <- CONTRAST[i]
  
  .colData <- as_tibble(SummarizedExperiment::colData(dds)) %>% 
    dplyr::select(contains(c("LIBRARY_ID",intgroup))) %>% drop_na(any_of(intgroup))
  
  which_group <- RES.P %>% filter(CONTRAST == intgroup)
  
  query.genes <- which_group %>% distinct(Name) %>% pull()
  
  cat("\nUsing ",length(query.genes), " DE genes\n")
  
  # DESeq2 performs an internal normalization where geometric mean is calculated for each gene across all samples.
  # i.e. median of ratios (https://doi.org/10.1186/gb-2010-11-10-r106)
  # Details: https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html
  
  cnts <- DESeq2::counts(dds, normalized = T, replaced = F)[query.genes,]
  
  # cnts <- DESeq2::varianceStabilizingTransformation(round(cnts))
  
  
  DF <- cnts %>% as_tibble(rownames = "Name") %>% 
    pivot_longer(-Name, names_to = "LIBRARY_ID", values_to = "count")
  
  DF <- RES.P %>% distinct(Name, Family) %>% 
    right_join(DF, by = "Name") %>% 
    left_join(.colData, by = "LIBRARY_ID") %>%
    drop_na(any_of(intgroup)) %>%
    mutate(CONTRAST = intgroup)
  
  names(DF)[names(DF) %in% intgroup] <- "Design"
  
  out[[i]] <-  DF

}

do.call(rbind, out) -> DF

DF %>% distinct(Name, CONTRAST) %>% dplyr::count(CONTRAST)

# CONTRAST A AND B) ====

intgroup <- CONTRAST[1:2]

which_group <- RES.P %>% filter(CONTRAST %in% intgroup)

which_group <- which_group %>%
  arrange(log2FoldChange) %>%
  group_by(Family) %>% 
  sample_n(1)

LOGFC_LABELLER <- which_group  %>% 
  arrange(CONTRAST) %>%
  mutate(log2FoldChange = paste0(Family, " (",log2FoldChange,") ", star)) %>%
  dplyr::select(Family, log2FoldChange) %>%
  pull(log2FoldChange, Family)

# downg <- which_group %>% filter(log2FoldChange < 0) %>% pull(Name)

str(LOGFC_LABELLER)

fun.data.trend <- "mean_se" # "mean_cl_boot", "mean_sdl"

recode_to_AB <- structure(c("pH 8.0", "pH 7.6"), names = c("Control", "Low"))

ylab <- expression('Expresión'~ '(' *log[2]* ')')  

DF %>%
  filter(CONTRAST %in% intgroup) %>% 
  arrange(count) %>% group_by(Family, LIBRARY_ID, Design, CONTRAST) %>% sample_n(1) %>%
  dplyr::mutate(CONTRAST = dplyr::recode(CONTRAST, !!!recode_to)) %>%
  dplyr::mutate(Design = dplyr::recode(Design, !!!recode_to_AB)) %>%
  mutate(CONTRAST = factor(CONTRAST, levels = recode_to[1:2])) %>%
  mutate(Design = factor(Design, levels = recode_to_AB)) %>%
  # mutate(Name = factor(Name, levels = unique(names(LOGFC_LABELLER)))) %>%
  ggplot(aes(x = Design, y=log2(round(count)+1), group = Family)) +
  # facet_wrap(CONTRAST ~ Family, labeller = labeller(Family = LOGFC_LABELLER), scales = "free_y") +
  ggh4x::facet_nested_wrap(~ CONTRAST+Family, labeller = labeller(Family = LOGFC_LABELLER), scales = "free_y", nest_line = F) +
  geom_jitter(position=position_jitter(w=0.1,h=0), size = 1, alpha = 0.5) +
  stat_summary(fun.data = fun.data.trend, colour = "red", linewidth = 0.7, size = 0.7, alpha = 0.7) +
  stat_summary(fun = mean, geom = "line", colour = "red") +
  labs(y = ylab, x = "") +
  theme_bw(base_family = "GillSans", base_size = 11) +
  theme(strip.background = element_rect(fill = 'grey89', color = 'white'),
    panel.border = element_blank(),
    plot.title = element_text(hjust = 0),
    plot.caption = element_text(hjust = 0),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    # axis.text.y.right = element_text(angle = 0, hjust = 1, vjust = 0, size = 2.5),
    axis.text.y = element_text(angle = 0, size = 7),
    axis.text.x = element_text(angle = 0, size = 10)) -> p


# 
ggsave(p, filename = 'DESEQ2UPEXPRESSED_MIRS_IN_REPONSE_TO_LOW_PH.png', path = wd, width = 7, height = 5, device = png, dpi = 300)

# 2) CONTRAST C AND D) ====

intgroup <- CONTRAST[3:4]

RES.P <- RES.P %>% filter(log2FoldChange <= -1)

which_group <- RES.P %>% filter(CONTRAST %in% intgroup)

which_unique <- which_group %>%
  group_by(Name) %>%
  summarise(
    across(CONTRAST, .fns = list), 
    n = n(),
    .groups = "drop_last") %>% filter(n == 1) %>% pull(Name)

which_group <- which_group %>% filter(Name %in% which_unique)

intgroup <- CONTRAST[4]

LOGFC_LABELLER <- which_group  %>% 
  filter(CONTRAST %in% intgroup) %>%
  arrange(CONTRAST) %>%
  mutate(log2FoldChange = paste0(Family, " (",log2FoldChange,")")) %>%
  dplyr::select(Family, log2FoldChange) %>%
  pull(log2FoldChange, Family)

str(LOGFC_LABELLER)

fun.data.trend <- "mean_se" # "mean_cl_boot", "mean_sdl"

recode_to_CD <- structure(c("110 HPF", "24 HPF"), names = c("Competent", "Control"))

DF %>%
  filter(CONTRAST %in% intgroup) %>% 
  filter(Name %in% which_unique) %>%
  dplyr::mutate(CONTRAST = dplyr::recode_factor(CONTRAST, !!!recode_to)) %>%
  dplyr::mutate(Design = dplyr::recode_factor(Design, !!!recode_to_CD)) %>%
  mutate(Name = factor(Name, levels = unique(names(LOGFC_LABELLER)))) %>%
  ggplot(aes(x = Design, y=log2(round(count)+1), group = Family)) +
  # facet_wrap(CONTRAST ~ Family, labeller = labeller(Family = LOGFC_LABELLER), scales = "free_y") +
  ggh4x::facet_nested_wrap(~ CONTRAST+Family, labeller = labeller(Family = LOGFC_LABELLER), scales = "free_y", nest_line = F) +
  geom_jitter(position=position_jitter(w=0.1,h=0), size = 1, alpha = 0.5) +
  stat_summary(fun.data = fun.data.trend, colour = "red", linewidth = 0.7, size = 0.7, alpha = 0.7) +
  stat_summary(fun = mean, geom = "line", colour = "red") +
  labs(y = "Expression (Log2)", x = "") +
  theme_bw(base_family = "GillSans", base_size = 11) +
  theme(strip.background = element_rect(fill = 'grey89', color = 'white'),
    panel.border = element_blank(),
    plot.title = element_text(hjust = 0),
    plot.caption = element_text(hjust = 0),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    # axis.text.y.right = element_text(angle = 0, hjust = 1, vjust = 0, size = 2.5),
    axis.text.y = element_text(angle = 0, size = 7),
    axis.text.x = element_text(angle = 90, size = 10)) -> p


# 
ggsave(p, filename = 'DESEQ2UPEXPRESSED_MIRS_IN_DEVELOPMENT_IN_REPONSE_TO_LOW_PH.png', path = wd, width = 5, height = 4.5, device = png, dpi = 300)

# UPGRADE ====
# EXCLUSIVE FROM THE 2_DESEQ2BARPLOT.R

.UPSETDF <- read_rds(paste0(wd, "UPSETDF.rds"))

EXCLUSIVE_MIRS <- .UPSETDF[[1]] %>% ungroup() %>%  
  mutate(Family = ifelse(!grepl("Cluster", Family), toupper(Family), Family))

INTERSECTED_MIRS <- .UPSETDF[[2]]  %>% ungroup() %>%  
  mutate(Family = ifelse(!grepl("Cluster", Family), toupper(Family), Family))

EXCLUSIVE_MIRS %>% 
  # distinct(Family, SIGN, CONTRAST) %>% 
  dplyr::count(CONTRAST, SIGN, sort = T) 

INTERSECTED_MIRS %>% 
  # distinct(Family, SIGN, CONTRAST) %>% 
  dplyr::count(CONTRAST, SIGN, sort = T) 


RES.P %>% 
  # filter(Name %in% which_unique) %>% 
  mutate(SIGN = sign(log2FoldChange)) %>%
  distinct(Name, CONTRAST, SIGN) %>% 
  dplyr::count(CONTRAST, SIGN, sort = T) %>% view()

str(which_unique <- EXCLUSIVE_MIRS %>% 
    filter(CONTRAST %in% c("CONTRAST_A", "CONTRAST_B")) %>%
    distinct(Name) %>% pull())


str(which_unique <- RES.P %>% 
    filter(CONTRAST %in% c("CONTRAST_A", "CONTRAST_B")) %>%
    distinct(Name) %>% pull())

RES.P %>%
  mutate(SIGN = sign(log2FoldChange)) %>%
  distinct(Name, Family, CONTRAST, SIGN) %>%
  mutate(Family = ifelse(Name %in% which_unique, paste0(Family, " *"), Family)) %>%
  pivot_wider(names_from = CONTRAST, values_from = SIGN, values_fill = 0) %>% view()
  # dplyr::count(CONTRAST, SIGN, sort = T)
#   mutate(CONTRAST = ifelse(SIGN == -1 & CONTRAST %in% c("CONTRAST_C"), "110 HPF (Ctr)", CONTRAST)) %>%
#   mutate(sampleB = ifelse(SIGN == -1 & CONTRAST %in% c("CONTRAST_D"), "110 HPF (Low)", sampleB)) %>%
#   mutate(sampleB = ifelse(SIGN == -1 & CONTRAST %in% c("CONTRAST_A"), "pH 7.6", sampleB)) %>%
#   mutate(sampleB = ifelse(SIGN ==  1 & CONTRAST %in% c("CONTRAST_A", "CONTRAST_B"), "pH 8.0", sampleB))
# dplyr::count(CONTRAST, SIGN, sort = T)

  
# divide the counts by the size factors or normalization factors before
barplot(cnts <- DESeq2::counts(dds, normalized = T, replaced = F)[which_unique,])

barplot(cnts <- DESeq2::varianceStabilizingTransformation(round(cnts)))

# as z score

cnts <- t(scale(t(cnts)))

DF <- cnts %>% as_tibble(rownames = "Name") %>% 
  pivot_longer(-Name, names_to = "LIBRARY_ID", values_to = "count")


MTD <- as_tibble(SummarizedExperiment::colData(dds)) %>% dplyr::select(LIBRARY_ID, colName,hpf, pH)

DF <- RES.P %>% distinct(Name, Family) %>% 
  right_join(DF, by = "Name") %>% 
  mutate(Family = ifelse(!grepl("Cluster", Family), toupper(Family), Family)) %>%
  group_by(LIBRARY_ID, Family) %>% 
  summarise(count = sum(count)) %>%
  left_join(MTD, by = "LIBRARY_ID")

# SPLIT MIRS IN RESPONSE TO LOW-PH DURING 24 AND 110 HPF:

# SORT_AB <- RES.P %>% 
#   mutate(SIGN = sign(log2FoldChange)) %>%
#   distinct(Name, Family, CONTRAST, SIGN) %>% 
#   # mutate(Family = ifelse(Name %in% which_unique, paste0(Family, " *"), Family)) %>%
#   mutate(Family = ifelse(!grepl("Cluster", Family), toupper(Family), Family)) %>%
#   filter(CONTRAST %in% c("CONTRAST_A", "CONTRAST_B")) %>%
#   arrange(Family) %>%
#   distinct(Family) %>% pull()


unique(DF$Family)

recode_to_hpf <- structure(c("24 HPF", "110 HPF"), names = c("24", "110"))

recode_to_pH <- structure(c("pH 8.0", "pH 7.6"), names = c("Control", "Low"))


scale_col_pH <- structure(c("#3B4992FF","#EE0000FF"), names = c("pH 8.0", "pH 7.6"))


fun.data.trend <- "mean_se" # "mean_cl_boot", "mean_sdl"


str(SORT_MIRS <- c("MIR-278","LET-7","MIR-133",
  "Cluster_55760","Cluster_55776","MIR-2","MIR-315", "MIR-153", "BANTAM", "MIR-190", "MIR-2722", "MIR-1988", 
  "MIR-92", "MIR-277B", "MIR-216"))



DF %>%
  ungroup() %>%
  dplyr::mutate(hpf = dplyr::recode_factor(hpf, !!!recode_to_hpf, .ordered = T)) %>%
  dplyr::mutate(pH = dplyr::recode_factor(pH, !!!recode_to_pH)) %>%
  mutate(Family = factor(Family, levels = SORT_MIRS)) %>%
  # mutate(y = log2(round(count)+1)) %>%
  mutate(y = count) %>%
  ggplot(aes(x = hpf, y = y, color = pH, fill = pH, group = pH)) +
  facet_wrap(~ Family, scales = "free_y") +
  geom_jitter(position=position_jitter(w=0.1,h=0), size = 1, alpha = 0.5) +
  stat_summary(fun.data = fun.data.trend, linewidth = 0.7, size = 0.7, alpha = 0.7) +
  stat_summary(fun = mean, geom = "line") +
  labs(y = "Read Count", x = "") +
  scale_color_manual("", values = scale_col_pH) +
  scale_fill_manual("", values =  scale_col_pH) +
  theme_bw(base_family = "GillSans", base_size = 11) +
  theme(strip.background = element_rect(fill = 'grey89', color = 'white'),
    legend.position = "top",
    panel.border = element_blank(),
    plot.title = element_text(hjust = 0),
    plot.caption = element_text(hjust = 0),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.y = element_text(angle = 0, size = 7),
    axis.text.x = element_text(angle = 90, size = 10)) -> p

ggsave(p, filename = 'DESEQ2LINEPLOT_ONLY_DEGS.png', path = wd, width = 5, height = 5, device = png, dpi = 300)


DF %>%
  dplyr::mutate(hpf = dplyr::recode_factor(hpf, !!!recode_to_hpf)) %>%
  dplyr::mutate(pH = dplyr::recode_factor(pH, !!!recode_to_pH)) %>%
  group_by(hpf, pH, Family) %>% summarise(mean = mean(count), sd = sd(count)) %>% view()

# out <- list()

# 
# for (i in query.genes) {
#   
#   Name <- i
#   
#   # which(res$transcript_id == gene)
#   
#   df <- plotCounts(dds, gene=Name, intgroup = intgroup, returnData=TRUE, normalized = T, transform = F)
#   
#   df <- data.frame(df, Name)  # %>% as_tibble(rownames = "LIBRARY_ID") 
#   
#   out[[i]] <-  df
# }

# head(DF <- do.call(rbind, out))
