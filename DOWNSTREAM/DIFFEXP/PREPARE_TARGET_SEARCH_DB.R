
# IDENTIFICAR PROTEINAS Y DOMINIOS PROTEICOS DE LA BIOMINERALIZACION
# 

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)


library(tidyverse)

url <- "https://raw.githubusercontent.com/RJEGR/Small-RNASeq-data-analysis/master/FUNCTIONS.R"

source(url)

wd <- "~/Documents/MIRNA_HALIOTIS/"

ref_path <- paste0(wd, "RNA_SEQ_ANALYSIS/ASSEMBLY/STRINGTIE/REFBASED_MODE/ANNOTATION/OUTPUTS/")

annot_f <- list.files(path = ref_path, pattern = "Trinotate.xls", full.names = T)

go_f <- list.files(path = ref_path, pattern = "Trinotate_report.xls.gene_ontology", full.names = T)

readMappings <-  function (file, sep = "\t", IDsep = ",")  {
  
  
  # From topGO::readMappings
  
  a <- read.delim(file = file, header = FALSE, quote = "", 
    sep = sep, colClasses = "character")
  map <- a[, 2]
  names(map) <- gsub(" ", "", a[, 1])
  return(lapply(map, function(x) gsub(" ", "", strsplit(x, 
    split = IDsep)[[1]])))
}

MAP <- readMappings(go_f)

annot <- read_tsv(annot_f, na = ".")

annot <- annot %>% drop_na(sprot_Top_BLASTP_hit)

names(annot)[1] <- "gene_id"

# annot %>% head() %>% view()

blastp_df <- split_blast(annot, hit = "BLASTP")

str(query.ids <- blastp_df %>% distinct(transcript) %>% pull(transcript)) 

# MAP <- MAP[names(MAP) %in% query.ids]

# go_id <- unlist(MAP)

# as.data.frame(go_id) %>% as_tibble(rownames = "transcript") %>% right_join(blastp_df)

# blastp_df <- blastp_df %>%
#   group_by(transcript) %>% 
#   arrange(identity) %>%
#   sample_n(1) %>%
#   ungroup()

blastp_df
# write_rds(blastp_df, file = paste0(ref_path, "blastp_df.rds"))

blastx_df <- split_blast(annot, hit = "BLASTX")

# write_rds(blastx_df, file = paste0(ref_path, "blastx_df.rds"))

go_df <- split_gene_ontology(annot, hit = "BLASTP")

pfam_df <- split_pfam(annot)

save(blastp_df, blastx_df, go_df, pfam_df, file = paste0(ref_path, "annot.Rdata"))

# EXPLORATORY   ======

go_df %>% count(ontology)

# hist(blastp_df$identity)

prot_list <- c("Perlucin", "Protein PIF", "Carbonic anhydrase", "Chitin synthase", "Insoluble matrix shell protein")

# prot_list <- str_to_lower(prot_list)

prot_list <- paste(prot_list, collapse = "|")

blastp_df %>% 
  # mutate_at(vars(all_of("name")), list(~ str_to_lower(.))) %>%
  filter_all(any_vars(grepl(prot_list, .))) %>% view() 
# mutate_at(vars(all_of(into)), list(~ str_to_sentence(.)))
# blastp_df %>% filter_all(any_vars(grepl("LOC", .))) %>%
  # distinct(gene, transcript, protein)

str(query_protein <- blastp_df %>% 
  filter_all(any_vars(grepl(prot_list, .))) %>%
  distinct(protein) %>% pull())

pfam_df %>% filter(protein %in% query_protein) %>% view()
go_df %>% filter(protein %in% query_protein) %>% view()

f <- list.files(path = ref_path, pattern = "transcripts.fa.transdecoder.rmdup.cds", full.names = T)


seqs <- Biostrings::readDNAStringSet(f)

str(keep <- names(seqs))

keep <- sapply(strsplit(keep, " "), `[`, 1)

str(keep <- unique(keep))

sum(keep <- keep %in% query_protein)

subseqs <- seqs[keep]# <--------




# blastp_df %>% distinct(name) %>% view()

blastp_df %>% distinct(uniprot) %>% view()

blastp_df %>% distinct(genus) %>% view()

blastp_df %>% 
  filter_all(any_vars(grepl("Mollusca", .))) %>%
  count(genus) %>% view()


blastp_df %>% 
  filter_all(any_vars(grepl("Mollusca", .))) %>% view()


# MAKE UTR SUBSET
# awk '{print $3}' transcripts.fa.transdecoder.gff3 | sort | uniq -c

# 108812 CDS
# 108812 exon
# 87680 five_prime_UTR
# 108812 gene
# 108812 mRNA
# 8424 start
# 100219 three_prime_UTR

# cat transcripts.fa.transdecoder.cds | seqkit replace -p "\s.+" > transcripts.fa.transdecoder.renamed.cds

# genome=transcripts.fa.transdecoder.renamed.cds

# gtf=transcripts.fa.transdecoder.gff3
# bed=transcripts.fa.transdecoder.bed


# cat $gtf | grep "three_prime_UTR" > prime_UTR.gff3

# cat prime_UTR.gff3.tmp | sed 's/$/";/' | sed 's/Parent=/Parent "/' > prime_UTR.gtf

# cat $bed | grep "prime" > prime_UTR.bed

# seqkit subseq --gtf prime_UTR.gff3 $genome --gtf-tag "Parent" --update-faidx -o prime_UTR.cds

# sequence () not found in file: transcripts.fa.transdecoder.renamed.cds

# But working:

# seqkit subseq --bed prime_UTR.gff3 $genome --chr transdecoder --update-faidx -o prime_UTR.cds

# cat prime_UTR.cds | seqkit rmdup -s -P -D duplicated_detail.txt -d duplicates.fasta -o utr_rmdup.fa


# GROUP GENES BY THEMATIC COMPONENTS ======

# EX. KEGG, BIOLOGICAL PROCESSS, CELLULAR COMPONENTS OR MOLEC. FUNCTION.

# GENE SAMPLE ABUNDANCE  ======



