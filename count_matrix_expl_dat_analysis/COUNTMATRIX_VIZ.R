
# THIS IS THE UPGRADE VERSION PREV. TO DE ANALYSIS


rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

path <- "~/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out/"

# 
pHpalette <- c(`Experimental`="#ad1f1f", `Control`= "#4575b4")
#

out <- read_rds(paste0(path, "Counts_and_mtd.rds"))

dim(COUNTS <- out[[1]])

dim(MTD <- out[[2]])

# head(sort_m <- sort(rowSums(COUNTS), decreasing = T))

# any(names(sort_m) %in% rownames(COUNTS))

# head(COUNTS[match(rownames(COUNTS), names(sort_m) ),])

# ORDER SRNA MATRIX BY DECREASING COUNT (Freq. and total read count):

head(SORTED_MAT <- COUNTS[order(rowSums(COUNTS),decreasing=T),])

READS_AND_FREQ <- COUNTS %>%
  as_tibble(rownames = "Name") %>%
  pivot_longer(-Name, values_to = "Reads") %>%
  mutate(REP = substr(name, 1,nchar(name)-1)) %>%
  filter(Reads > 0) %>%
  group_by(Name, REP) %>%
  summarise(Reads = sum(Reads), Freq = n())

READS_AND_FREQ %>% group_by(REP) %>% cor_test(Reads, Freq)
READS_AND_FREQ %>% select(-Freq) %>% pivot_wider(names_from = REP, values_from = Reads, values_fill = 0)
READS_AND_FREQ %>% select(-Reads) %>% pivot_wider(names_from = REP, values_from = Freq, values_fill = 0)


head(DB <- read_tsv(paste0(path, "/RNA_LOCATION_DB.tsv")))

DB <- READS_AND_FREQ %>% summarise(Reads = sum(Reads), Freq = sum(Freq)) %>%
  left_join(DB) %>% rename("SampleFreq" = "Freq") %>%
  arrange(match(Name, rownames(SORTED_MAT)))

DB  %>% ggplot(aes(x = Reads, y = SampleFreq, color = SRNAtype)) + 
  stat_ecdf()

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

# TRY SEQUENCE SIMILARITY ANALYSIS AND 

seqs <- DB %>% distinct(MajorRNA) %>% pull() %>% head(100)

bioseq::rna(seqs)

# print(seqs <- Biostrings::RNAStringSet(seqs))

# Biostrings::pairwiseAlignment(seqs, seqs)
  
# alignedDNA <- DECIPHER::AlignSeqs(RNAStringSet(seqs))

bioseq::seq_consensus(bioseq::rna(seqs))



library(msa)

# method=c("ClustalW", "ClustalOmega", "Muscle")

align <- msa::msa(RNAStringSet(seqs), method = "ClustalW")

align <- msaConvert(align)$seq

print(bioseq::rna(align))

bioseq::seq_cluster(bioseq::rna(align))


# READS_AND_FREQ %>% ungroup() %>% 

# PREVALENCE BY REP ----

nrow(raw_count <- COUNTS[colNames]) # 66,226

head(raw_count <- as.data.frame(raw_count, row.names = COUNTS$Name))

rownames(raw_count) <- COUNTS$Name

apply(raw_count, 1, function(x) sum(x > 0)) %>% table()

prevelancedf = apply(raw_count, 1, function(x) sum(x > 0))

# mean_se = apply(raw_count, 1, function(x) mean_se(x)) %>% do.call(rbind, .)

data.frame(Prevalence = prevelancedf, 
  TotalAbundance = rowSums(raw_count)) %>% # mean_se
  as_tibble(rownames = "Name") %>%
  arrange(desc(TotalAbundance)) -> prevelancedf

prevelancedf <- prevelancedf %>% left_join(RESULTS %>% select(Name, MIRNA, KnownRNAs, Strand, DicerCall))

prevelancedf %>% 
  count(Prevalence) %>% 
  ggplot(aes(Prevalence, n)) + geom_col() +
  theme_classic(base_family = "GillSans") + 
  scale_y_continuous("Number of sRNAs", labels = scales::comma) +
  scale_x_continuous(breaks = 1:12) -> ps

ggsave(ps, filename = 'prevalence_hist.png', 
  path = path, width = 3, height = 2, device = png)

prevelancedf %>% 
  arrange(Prevalence) %>%
  mutate(Prevalence = paste0(Prevalence, ' Samples')) %>%
  mutate(Prevalence = factor(Prevalence, levels = unique(Prevalence))) %>%
  mutate(TotalAbundance = log2(TotalAbundance+1)) %>% # edgeR::cpm(count) 
  ggplot(aes(TotalAbundance)) + geom_histogram() + 
  facet_wrap(~ Prevalence, scales = 'free_y') -> p1

dat_text <- prevelancedf %>% group_by(Prevalence) %>% tally() %>% 
  mutate(cumsum = cumsum(n)) %>% arrange(Prevalence) %>%
  mutate(Prevalence = paste0(Prevalence, ' Samples')) %>%
  mutate(Prevalence = factor(Prevalence, levels = unique(Prevalence)))

p1 <- p1 + geom_text(
  data    = dat_text, family = "GillSans",
  mapping = aes(x = -Inf, y = -Inf, label = paste0(n, " sRNAs")),
  hjust   = -1, vjust   = -2.5, size = 2.5) + 
  theme_classic(base_size = 7, base_family = "GillSans") +
  labs(x = expression(~Log[2]~('TotalAbundance'~+1)), y = "") +
  scale_y_continuous("Number of sRNAs", labels = scales::comma)


ggsave(p1, filename = 'prevalence_hist_facet.png', 
  path = path, width = 7, height = 4, device = png)


prevelancedf %>% 
  mutate(MIRNA = factor(MIRNA, levels = c("Y","N") )) %>%
  # filter(MIRNA == "Y") %>%
  ggplot(aes(TotalAbundance, Prevalence/12)) +
  # geom_point(aes(color = MIRNA), alpha = 0.2) +
  stat_density2d(aes(alpha = ..density..), geom = "tile", contour = FALSE) +
  scale_x_log10("Total Abundance (log10 scale)", labels = scales::comma) +  
  scale_y_continuous("Prevalence [Frac. Samples]", 
    labels = scales::percent_format(scale = 100)) +
  theme_classic(base_family = "GillSans") +
  theme(legend.position="top") -> ps
# facet_grid(~ MIRNA) 
# scale_color_manual(values = c("red", "grey"))


ggsave(ps, filename = 'prevalence.png', 
  path = path, width = 4, height = 4, device = png)


prevelancedf %>% 
  mutate(Prevalence = factor(Prevalence)) %>%
  # filter(MIRNA == "Y") %>%
  ggplot(aes(log10(TotalAbundance), color  = Prevalence, fill = Prevalence)) +
  # scale_x_continuous() +
  stat_ecdf() 
# geom_density(alpha = 0.3)



