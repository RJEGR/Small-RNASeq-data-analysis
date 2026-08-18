# =============================================================================
# consolidate_mirna_tibble.R
# Consolidates ShortStack MIRNA loci into a single tibble carrying
#   hairpin / mature / miRNA* sequences (mir.fasta)
#   + MFE of folding (strucVis/*.txt)
#   + KnownRNAsDB annotation (RNA_LOCATION_MIR_DB.rds)
# Output: MIRNA_CONSOLIDATED.tsv / .rds  (147 loci x 39 columns)
# =============================================================================

library(tidyverse)
library(Biostrings)

dir <- "~/Documents/Transcriptomics/MicroRNA_target_validation/"

# ---- 0. inputs --------------------------------------------------------------
KnownRNAsDB <- read_rds(paste0(dir, "RNA_LOCATION_MIR_DB.rds"))
RESULTS     <- read_tsv(paste0(dir, "Shortstacks_outputs/Results.txt"),
                        show_col_types = FALSE)

# ---- 1. MIRS (unchanged logic from Supplementary_data.R) --------------------
SEQS <- RESULTS %>% select(Name, MajorRNA)

paste_headers <- function(x) paste(x, sep = "|", collapse = "|")

SPLIT <- RESULTS %>%
  select(Name, KnownRNAs, MIRNA) %>%
  drop_na(KnownRNAs) %>%
  mutate(KnownRNAs = strsplit(KnownRNAs, ";")) %>%
  unnest(KnownRNAs) %>%
  mutate(biotype = ifelse(grepl("^piR-", KnownRNAs), "piRNA", "Mir"))

MIRS <- SPLIT %>%
  filter(biotype == "Mir") %>%
  group_by(Name) %>%
  summarise(across(KnownRNAs, .fns = paste_headers), n = n(), .groups = "drop")

MIRS <- SPLIT %>%
  filter(MIRNA == "Y" & biotype == "piRNA") %>%
  group_by(Name) %>%
  summarise(across(KnownRNAs, .fns = paste_headers), n = n(), .groups = "drop") %>%
  filter(!Name %in% MIRS$Name) %>%
  rbind(MIRS) %>%
  left_join(SEQS, by = "Name")

NOVEL_MIRS <- RESULTS %>%
  filter(is.na(KnownRNAs), MIRNA == "Y") %>%
  select(Name, MajorRNA)

MIRS <- NOVEL_MIRS %>%
  mutate(KnownRNAs = "Novel_miR", n = 1) %>%
  select(names(MIRS)) %>%
  rbind(MIRS)

# ---- 2. escoresdf : MFE of folding, line 7 of each strucVis txt -------------
f <- list.files(paste0(dir, "Shortstacks_outputs/strucVis"),
                pattern = "txt$", full.names = TRUE)

read_escore <- function(f) {
  x         <- read_lines(f, n_max = 1, skip = 6)
  split_str <- sapply(strsplit(x, " "), `[`, 2)
  tibble(Name      = gsub(".txt", "", basename(f)),
         mfe_score = as.double(gsub("\\(|\\)", "", split_str)))
}

escoresdf <- map_dfr(f, read_escore)

# ---- 3. mir.fasta -> hairpin / mature / star, wide by locus -----------------
dna <- readDNAStringSet(paste0(dir, "Shortstacks_outputs/mir.fasta"))

FA <- tibble(header = names(dna), seq = as.character(dna)) %>%
  separate(header, into = c("id", "coord"), sep = "::", extra = "merge") %>%
  mutate(part = case_when(str_detect(id, "\\.mature$") ~ "mature",
                          str_detect(id, "\\.star$")   ~ "star",
                          TRUE                         ~ "hairpin"),
         Name = str_remove(id, "\\.(mature|star)$")) %>%
  select(Name, part, seq, coord) %>%
  pivot_wider(names_from = part, values_from = c(seq, coord)) %>%
  dplyr::rename(hairpin_seq   = seq_hairpin,   mature_seq   = seq_mature,   star_seq   = seq_star,
         hairpin_coord = coord_hairpin, mature_coord = coord_mature, star_coord = coord_star) %>%
  mutate(hairpin_len = nchar(hairpin_seq),
         mature_len  = nchar(mature_seq),
         star_len    = nchar(star_seq))

# ---- 4. which arm does MajorRNA come from? ---------------------------------
# MajorRNA is RNA (U), mir.fasta is DNA (T). Exact match -> mature / star;
# otherwise place MajorRNA inside the hairpin and assign by overlap (isomiR).
rc <- function(s) as.character(reverseComplement(DNAString(s)))

classify_arm <- function(MajorRNA, hairpin, mature, star) {
  if (is.na(hairpin)) return(c(NA_character_, NA_character_, NA_real_))
  q <- gsub("U", "T", MajorRNA)
  if (identical(q, mature)) return(c("mature", "exact", str_locate(hairpin, fixed(mature))[1]))
  if (identical(q, star))   return(c("star",   "exact", str_locate(hairpin, fixed(star))[1]))
  p <- str_locate(hairpin, fixed(q))[1]; strand <- "sense"
  if (is.na(p)) { p <- str_locate(hairpin, fixed(rc(q)))[1]; strand <- "antisense" }
  if (is.na(p)) return(c("unplaced", "no_match", NA_real_))
  a <- p; b <- p + nchar(q)
  ov <- function(sub) { i <- str_locate(hairpin, fixed(sub))[1]
  if (is.na(i)) 0 else max(0, min(b, i + nchar(sub)) - max(a, i)) }
  om <- ov(mature); os <- ov(star)
  arm <- if (om >= os && om > 0) "mature" else if (os > 0) "star" else "unplaced"
  c(arm, paste0("variant_", strand), p)
}

# ---- 5. consolidate --------------------------------------------------------
CONSOLIDATED <- MIRS %>%
  distinct(Name, MajorRNA) %>%
  left_join(escoresdf,  by = "Name") %>%
  left_join(KnownRNAsDB, by = c("Name", "MajorRNA")) %>%
  left_join(FA,          by = "Name") %>%
  mutate(.arm = pmap(list(MajorRNA, hairpin_seq, mature_seq, star_seq), classify_arm),
         MajorRNA_arm         = map_chr(.arm, 1),
         MajorRNA_match       = map_chr(.arm, 2),
         MajorRNA_hairpin_pos = as.numeric(map_chr(.arm, 3))) %>%
  select(-.arm) %>%
  select(Name, Family, KnownRNAs, MIRNA, Locus, Chrom, Start, End, Strand, DicerCall,
         Locus_type, SRNAtype, biotype, biotype_best_rank, type, gene_id, transcript_id, WGCNA,
         Reads, UniqueReads, MajorRNAReads, FracTop, SampleFreq, Freq, NKnownRNAs,
         MajorRNA, MajorRNA_arm, MajorRNA_match, MajorRNA_hairpin_pos,
         hairpin_seq, hairpin_len, hairpin_coord, mfe_score,
         mature_seq, mature_len, mature_coord,
         star_seq, star_len, star_coord)

CONSOLIDATED

write_tsv(CONSOLIDATED, paste0(dir, "MIRNA_CONSOLIDATED.tsv"))
# write_rds(CONSOLIDATED, paste0(dir, "MIRNA_CONSOLIDATED.rds"))