
# RUN GO ANALYSIS USING:
# 0) PREPARE UPSET 
# 1) SPLIT DOWN FROM UP SIGNIFICANT (PADJ < 0.05) VALUES
# 2) RUN INDEPENDENTLY THE TOP GO ANALYSIS

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

wd <- "~/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out/"

RES <- read_tsv(paste0(wd, "DESEQ_RES.tsv"))

wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/FUNCTIONAL_MIR_ANNOT/"

print(SRNA2GO <- read_tsv(paste0(wd, "SRNA_REGULATORY_FUNCTION_DB.tsv")))

url <- "https://raw.githubusercontent.com/RJEGR/Cancer_sete_T_assembly/main/functions.R"

source(url)

# Mature are you functional strand used in DE

SRNA2GO <- SRNA2GO %>%
  filter(predicted == "BOTH") %>%
  mutate(query = strsplit(query, ";")) %>%
  unnest(query) %>% 
  distinct(query, GO.ID) %>%
  group_by(query) %>%
  summarise(
    across(GO.ID, .fns = paste_go), 
    .groups = "drop_last")


SRNA2GO <- SRNA2GO %>% 
  separate(query, into = c("query", "arm"), sep = "[.]") %>%
  filter(arm == "mature")


SRNA2GO <- split(strsplit(SRNA2GO$GO.ID, ";") , SRNA2GO$query)

SRNA2GO <- lapply(SRNA2GO, unlist)

DF <- list()

for(j in 1:length(CONTRAST)) {
  
  i <- j
  
  cat("\nRunning ",CONTRAST[i], "\n")
  
  DF2GO <- RES %>% filter(CONTRAST %in% CONTRAST[i])
  
  DF2GO <- DF2GO %>% filter(log2FoldChange > 0 & padj < 0.05) # filter(log2FoldChange < 0 )
  
  query.p <- DF2GO %>% pull(padj)
  
  query.names <- DF2GO %>% pull(Name)
  
  allRes <- GOenrichment(query.p, query.names, SRNA2GO, Nodes = 20)
  
  allRes <- allRes %>% mutate(CONTRAST = CONTRAST[i]) %>% as_tibble()
  
  DF[[i]] <- allRes
  
}

do.call(rbind, DF) -> allRes

allRes %>% as_tibble() %>%
  arrange(Annotated) %>%
  dplyr::mutate(CONTRAST = dplyr::recode_factor(CONTRAST, !!!recode_to)) %>%
  mutate(Term = factor(Term, levels= unique(Term))) %>%
  ggplot(aes(x = Term, y = Annotated, fill = -log10(p.adj.ks))) + # , fill = p.adj.ks
  coord_flip() +
  geom_col() +
  theme_bw(base_family = "GillSans", base_size = 12) +
  facet_grid( CONTRAST ~ ., scales = "free") +
  theme_bw(base_family = "GillSans", base_size = 11) +
  theme(strip.background = element_rect(fill = 'grey86', color = 'white'),
    legend.position = 'top',
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  ggsci::scale_fill_material()


# Format to bind the DB

allRes %>% filter(is.na(log2FoldChange))

allRes %>% 
  
  res.p <- prep_DE_data(allRes, alpha = 0.05, lfcThreshold = 2)

allRes %>% pivot_wider()

DB %>% left_join(allRes %>% select(log2FoldChange, padj, CONTRAST))

# ENRICHMENT ANALYSIS: =====

.bioc_packages <- c('biomaRt','GOSemSim',"org.Hs.eg.db") # DESeq2'

sapply(c(.bioc_packages), require, character.only = TRUE)

# PLOT ====
# POSITIVE LFC == UP EXPRESSED IN EXPERIMENTAL
# NEGATIVE LFC == UP EXPRESSED IN CONTROL

# recode_to <- paste(LETTERS[1:4], ")", sep = "")

recode_to <- c("24 HPF: Ctrl|pH", "110 HPF: Ctrl|pH", "Ctrl: 24|110 HPF", "pH: 24|110")
recode_to <- structure(recode_to, names = CONTRAST)

# res.p <- prep_DE_data(allRes, alpha = 0.05, lfcThreshold = 2) %>% drop_na(cc)

sigfc <- "Sign and FC";pv <- "Sign";fc <- "FC"

colors_fc <- c("red2",  "#4169E1", "forestgreen", "grey30")

colors_fc <- structure(colors_fc, names = c(sigfc, pv, fc, "NS"))

allRes$Color <- 'NS'
allRes[which(abs(allRes$log2FoldChange) > 2), 'Color'] <- fc
allRes[which(abs(allRes$padj) <= 0.05), 'Color'] <- pv
allRes[which(allRes$padj <= 0.05 & abs(allRes$log2FoldChange) > 2), 'Color'] <- sigfc

colors_fc <- structure(colors_fc, names = c(sigfc, pv, fc, "NS"))

p <- allRes  %>%
  dplyr::mutate(CONTRAST = dplyr::recode_factor(CONTRAST, !!!recode_to))
ggplot(aes(y = -log10(pvalue), x = log2FoldChange, color = Color)) +
  geom_point()  

p <- p +
  scale_color_manual(name = "", values = colors_fc) +
  labs(x= expression(Log[2] ~ "Fold Change"), y = expression(-Log[10] ~ "P")) +
  theme_bw(base_family = "GillSans", base_size = 18) +
  theme(legend.position = "top")  +
  geom_abline(slope = 0, intercept = -log10(0.05), linetype="dashed", alpha=0.5) +
  geom_vline(xintercept = 0, linetype="dashed", alpha=0.5) 

p + facet_grid(~ CONTRAST)

allRes %>% 
  mutate(g = sign(log2FoldChange)) %>% 
  dplyr::count(Color, g, CONTRAST) %>%
  # mutate(g = ifelse(g == "1", "CON_CANCER", "SIN_CANCER")) %>%
  filter(Color != 'NS') %>%
  group_by(g, CONTRAST) %>%
  mutate(pct = n / sum(n)) %>%
  ggplot() +
  facet_grid(g ~ .) +
  geom_col(aes(x = CONTRAST, y = pct, fill = Color), width = 0.5) +
  scale_fill_manual(name = "", values = colors_fc[-4]) +
  theme_bw(base_family = "GillSans", base_size = 18) +
  labs(y = '% Transcripts', x = '') +
  theme(legend.position = 'top') + coord_flip()

