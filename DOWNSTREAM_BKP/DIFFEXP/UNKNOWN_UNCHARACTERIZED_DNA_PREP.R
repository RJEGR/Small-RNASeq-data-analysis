# OBSOLETO: ESTO YA NO ES NECESARIO (AGOSTO 2023)
# SUBSEQ SEQUENCE BASED ON COORDINATES CLOSE TO UTR-TARGET SITES
# PREPARE UNKNOWN + UNCHARACTERIZED SEQUENCES TO RERUN BLAST AGAINST UNIPRONT/NR

# UNNEST DATA
df_drop <- read_rds(paste0(wd, "/mir_vs_utr_rmdup_RNAhybrid.out.psig.gene_features.rds"))

df_drop %>%
  drop_na(gene_id) %>%
  mutate(query = strsplit(query, ";")) %>%
  unnest(query) %>% distinct(query)

df_drop %>% drop_na(gene_id) 

# GET BACK GENE ONTOLOGY FROM TRANSCRIPTOMIC DATA USING gene_id

# OR RERUN BLASTX (AGAINST UNIPROT OR TRANSCRIPTOME) USING THOSE GENE_ID LOCUS 

transcr_wd <- "~/Documents/MIRNA_HALIOTIS/RNA_SEQ_ANALYSIS/ASSEMBLY/STRINGTIE/REFBASED_MODE/ANNOTATION/OUTPUTS/"

f <- list.files(path = transcr_wd, pattern = "annot.Rdata", full.names = T)

load(f)

# filter(blastx_df, grepl("MSTRG.40173", gene))

filter(blastx_df, grepl("LOC", gene)) %>%
  mutate(gene = gsub("_df_pt", "", gene)) %>%
  filter(gene %in% unique(df_drop$gene_id))

filter(blastp_df, grepl("LOC", gene)) %>% 
  mutate(gene = gsub("_df_pt", "", gene)) %>% 
  filter(gene %in% unique(df_drop$gene_id))


head(unique(df_drop$gene_id))

# NO HAY IDENTIDAD DE LOS IDS DEL GENOMA EN EL TRANSCRIPTOMA, HABRA QUE CORRER BLAST VS UNIPROT OR TRANSCRIPTOME
# O USAR COORDENADAS PARA DEVOLVER LA IDENTIDAD


# PREPARE COORDS FOR BLASTX (AGAINST UNIPROT OR TRANSCRIPTOME) ====
# 1) USING NON ANNOTTED SEQUENCES 

unch_df <- df_drop %>% filter(grepl("uncharacterized", description))
unknown_df <- df_drop %>% filter(is.na(gene_id))

coords_df <- rbind(unch_df, unknown_df) %>%
  separate(col = target, into = c("seqnames", "ranges"), sep = "_") %>%
  separate(col = ranges, into = c("ranges", "strand"), sep = ":") %>%
  separate(col = ranges, into = c("start", "end"), sep = "-") 

coords_df <- coords_df %>% mutate_at(c("start", "end"), as.numeric)

ranges_ <- IRanges(start = coords_df$start, end = coords_df$end) # Those UTR ranges

seqnames_ <- coords_df$seqnames

grr <- GRanges(Rle(seqnames_), 
  ranges = ranges_, strand = coords_df$strand, query = coords_df$query)

nrow(coords_df) # 31 545 unknown + uncharacterized sites

# Transform to single UTR regions 

# grr <- IRanges::reduce(grr) # 27 925 reduced unknown + uncharacterized sites

which_chr <- levels(seqnames(grr))

grr <- split(grr, seqnames(grr))

# Range the exon region (CDS+three_prime_utr) using subsetByOverlaps
# 2)


gr <- .gr[is.na(.gr$gene_id)]

length(gr <- gr[seqnames(gr) %in% which_chr]) # 1456696

seqlevels(gr) <- unique(unfactor(seqnames(gr)@values))

# gr %>% as_tibble() %>% dplyr::count(type) # 664 202

# IRanges::reduce(gr)

# promoters(grr, upstream = 1e3)
## split into a GRangesList
## where each element has all ranges for one chromosome

x <- grr[[1]]

fo <- subsetByOverlaps(x = gr, ranges = x)  

library(gggenomes)
library(gggenes)

fo %>% as_tibble() %>% dplyr::count(type)

fo %>% as_tibble() %>%
  filter(grepl("mRNA|prime_UTR", type)) %>%
  ggplot(aes(xmin = start, xmax = end, y = seq_id, fill = type)) +
  geom_gene_arrow() +
  facet_wrap(~ type, scales = "free", ncol = 1) +
  scale_fill_brewer(palette = "Set3")


# keep mRNA

fo <- fo[fo$type == "mRNA"]

IRanges::reduce(fo) %>% as_tibble() %>%
  ggplot(aes(xmin = start, xmax = end, y = seq_id)) +
  geom_gene_arrow() +
  # facet_wrap(~ type, scales = "free", ncol = 1) +
  scale_fill_brewer(palette = "Set3")

fo <- IRanges::reduce(fo)


# load fasta (toplevel)
f <- "Haliotis_rufescens_gca023055435v1rs.xgHalRufe1.0.p.dna.toplevel.fa"
dna_wd <- "~/Documents/MIRNA_HALIOTIS/ENSEMBLE"

f <- list.files(path = dna_wd, pattern = f, full.names = T)

dna <- Biostrings::readDNAStringSet(f)

head(query_dna <- sapply(strsplit(names(dna), " "), `[`, 1))

# 3)

keep <- query_dna %in% as.character(seqnames(fo)@values)

seq_id <- fo %>% as_tibble() %>%
  mutate(gene_coords = paste(start, end, strand, sep = ":")) %>%
  mutate(seq_id = paste(seq_id,gene_coords, sep = "_")) %>%
  pull(seq_id)

s <- start(fo)
e <- end(fo)
w <- width(fo)

pb <- txtProgressBar(min = 0, max = length(fo), style = 3, width = 50, char = "=")

out <- list()

for(i in 1:length(fo)) {
  
  j <- i
  
  out[[j]] <- Biostrings::subseq(dna[keep], start = s[j], end = e[j])
  
  # names(out[[j]]) <- seq_id[[j]]
  
  setTxtProgressBar(pb, j)
  
}

# TEST AS SINGLE FUNCTION

out <- find_genes_overlaps(query_gr = x, subject_gr =  gr)


## apply a function to the ranges of each chromosome
out <- vapply(grr, find_chr_overlaps, subject_gr = gr)


# out <- do.call(c, out)

seqs <- as.character(out)

headers <- paste0(">", names(out))

fasta <- c(rbind(headers, seqs))

write(fasta, file= paste0(wd, "/test.fasta"))



??Biostrings::subseq()

# Get only genes types

find_genes_overlaps <- function(query_gr, subject_gr) {
  
  gr <- subject_gr
  
  gr <- sortSeqlevels(gr)
  
  gr <- gr[seqnames(gr) == as.character(seqnames(query_gr)@values)]
  
  fo <- subsetByOverlaps(x = gr, ranges = query_gr)
  
  fo <- fo[fo$type == "mRNA"]
  
  fo <- IRanges::reduce(fo)
  
  # fo <- findOverlaps(query_gr, gr, minoverlap = 1)
  
  # return(fo)
  
  seq_id <- fo %>% as_tibble() %>%
    mutate(gene_coords = paste(start, end, strand, sep = ":")) %>%
    mutate(seq_id = paste(seq_id,gene_coords, sep = "_")) %>%
    pull(seq_id)
  
  s <- start(fo)
  e <- end(fo)
  w <- width(fo)
  
  pb <- txtProgressBar(min = 0, max = length(fo), style = 3, width = 50, char = "=")
  
  out <- list()
  
  for(i in 1:length(fo)) {
    
    j <- i
    
    out[[j]] <- Biostrings::subseq(dna[keep], start = s[j], end = e[j])
    
    setTxtProgressBar(pb, j)
    
  }
  
  out <- do.call(c, out)
  
  names(out) <- seq_id
  
  return(out)
  
}


