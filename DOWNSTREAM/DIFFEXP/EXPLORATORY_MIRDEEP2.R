# mirdeep2 exploratory processing

# Arf format: Is a proprietary file format generated and processed by miRDeep2. It contains information of reads mapped to a reference genome. Each line in such a file contains 13 columns:
  
# read identifier
# length of read sequence
# start position in read sequence that is mapped
# end position in read sequence that is mapped
# read sequence
# identifier of the genome-part to which a read is mapped to. This is either a scaffold id or a chromosome name
# length of the genome sequence a read is mapped to
# start position in the genome where a read is mapped to
# end position in the genome where a read is mapped to
# genome sequence to which a read is mapped
# genome strand information. Plus means the read is aligned to the sense-strand of the genome. Minus means it is aligned to the antisense-strand of the genome.
# Number of mismatches in the read mapping
# Edit string that indicates matches by lowercase 'm' and mismatches by uppercase 'M'

path <- "/Users/cigom/Documents/MIRNA_HALIOTIS/CLEAN_INPUT/TEST_09022023"

fileName <- list.files(path = path, pattern = 'reads_collapsed_vs_genome.arf', full.names = T)

arf <- read_tsv(file = fileName, col_names = F)

arf[, 6] %>% filter(X6 == 'XR_006859493.2')

arf[, 6] %>% distinct(X6) %>% pull(.) -> rna_ids

# arf %>% count(X4/X7)
# which ids start with X[R|M]

sum(grepl("^X*", rna_ids)) # all the 54,083


gff_rna <- read_rds(paste0(path, "gff_rna.rds"))

gff_rna %>% count(type)

# 2 CDS        591,909
# 3 exon       665,954
# ...

gff_rna %>% filter(Parent %in% rna_ids) -> gff_rna_subset

gff_rna_subset %>% count(type)

# type       n
# <fct>  <int>
#   1 CDS   472,320
# 2 exon  520,858

# YA SABEMOS QUE LOS FRAGMENTOS DE SRNA PROVIENEN DE CDS Y EXONICOS (aparentemente), FALTA DETERMINAR LAS COORDENADAS ACORDE AL ARCHIVO GFF_rna_subset CREADO CON EL SCRIPT EXPLORATORY_GFF

arf %>% filter(X1 %in% "HR1_6762731_x220018")

