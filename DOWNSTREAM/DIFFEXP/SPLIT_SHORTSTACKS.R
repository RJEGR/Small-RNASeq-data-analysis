
rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

path <- "~/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out/"

list.files(path = path, pattern = "txt")

res_f <- list.files(path = path, pattern = "Results.txt", full.names = T)

RESULTS <- read_tsv(res_f)


RESULTS %>% count(MIRNA) # 131 True miRs

# Number of novel and know miRs annotated

RES <- RESULTS %>% 
  select(Name, Chrom, KnownRNAs, MIRNA, Reads, UniqueReads, MajorRNA) %>%
  mutate(KnownRNAs = ifelse(is.na(KnownRNAs) & MIRNA == "Y", Name, KnownRNAs)) %>%
  mutate(MIRNA = ifelse(!is.na(KnownRNAs), "Y", MIRNA)) %>%
  drop_na(KnownRNAs) %>%
  mutate(Type = ifelse(grepl("Cluster_", KnownRNAs), "Novel", "known"))

RES %>% count(Type)


RES %>% head() %>% view()

# PREPARE FASTA  =====

paste_headers <- function(x) {
  x <- paste(x, sep = '|', collapse = '|')
  # x <- list(x)
  # x <- unlist(x)
}

fasta_prep <- RES %>%
  group_by(MajorRNA) %>%
  summarise(across(Name, .fns = paste_headers), .groups = "drop_last") %>% 
  ungroup() %>% filter( nchar(MajorRNA) > 17) %>% arrange(MajorRNA)

seqs <- fasta_prep %>% pull(MajorRNA)

headers <- fasta_prep %>% mutate(Name = paste0(">", Name)) %>% pull(Name) 

fasta <- c(rbind(headers, seqs))

write(fasta, file= paste0(path, "MajorRNA.fasta"))

#
# scp -r MajorRNA.fasta rvazquez@200.23.162.234://home/rvazquez/MIRS_FUNCTIONAL_ANNOT/
# 

# FURTHER STATS ====

# Number of sRNAs?

# Algunas secuencias de piRs son anotadas como miRs, ejemplo: 
# piR-bgl-168596;Esc-Mir-2-o20-v2_3p > Cluster_44687

RES <- RES %>%
  mutate(sRNA = "Novel") %>%
  mutate(sRNA = ifelse(grepl("piR-",KnownRNAs), "piR", sRNA)) %>%
  mutate(sRNA = ifelse(grepl("Mir-",KnownRNAs), "miR", sRNA))

RES %>% count(sRNA)
 
RES %>% 
  mutate(x = log10(Reads), y = log10(UniqueReads)) %>%
  ggplot(aes(x, y, color = sRNA, shape = Type)) +
  geom_point(size = 5, alpha = 0.7) +
  labs(caption = "Diff piR pattern expression from novel and miRs")

# separate know MIRS and Novel from piRs

str(query.seqs <- RES %>% 
    # filter(sRNA != "piR") %>% 
    distinct(Name) %>% 
    pull(Name))


# RES <- RESULTS %>% filter(Name %in% query.seqs)

which_cols <- c("Short", 21:30, "Long")

out <- RESULTS %>% 
  pivot_longer(cols = all_of(which_cols), names_to = "length") %>%
  select(Name, Chrom, length) %>% 
  filter(Name %in% query.seqs) %>%
  left_join(RES) 

out %>%
  group_by(length, sRNA) %>% tally() %>%
  ggplot(aes(x = length, y = n, fill = sRNA)) +
  geom_col(width = 0.85,position = "dodge2")


RESULTS %>% 
  filter(Name %in% query.seqs) %>%
  distinct(MajorRNA) %>% 
  arrange(MajorRNA) %>% 
  mutate(n = nchar(MajorRNA)) %>% 
  ggplot(aes(n)) + geom_histogram()

# revisar la anotacion de los 978 sRNAs (i.e 281 secuencias unicos) 
# y cotejar su anotacion

df <- RES %>% 
  drop_na(KnownRNAs) %>%
  mutate(KnownRNAs = strsplit(as.character(KnownRNAs), ";")) %>% 
  unnest(KnownRNAs) %>% distinct(KnownRNAs ,.keep_all = T)
  # filter(Name %in% query.seqs) 

df %>% count(sRNA)


# ENCONTRAR TABLA DE ESPECIE Y FUNCION CONOCIDA DE KNOWN MIRNAS: https://mirgenedb.org/

# sRNA2annot <- split(df$KnownRNAs, df$Name)



# MULTI-MAP DETAILS

# Reads: Number of aligned sRNA-seq reads that overlap this locus.
# UniqueReads: Number of uniquely aligned (e.g. not multi-mapping) reads that overlap this locus.

# MajorRNA: Sequence of the single most abundant RNA sequence at the locus.
# MajorRNAReads: Number of reads for the MajorRNA aligned to this locus.

RESULTS %>% filter(Length > 200) %>% view()

RESULTS %>%
  select(Name, Length, Reads, UniqueReads, MajorRNAReads) %>%
  pivot_longer(cols = c(Reads, UniqueReads, MajorRNAReads)) %>%
  ggplot(aes(Length, log10(value), color=name)) + geom_point()
  # ggplot(aes(log10(value), fill = name)) +  geom_histogram()
