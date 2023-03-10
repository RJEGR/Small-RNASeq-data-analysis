
# GENOME FEATURES
# COMPARAR VERSION DE GENOMAS (GFF) PARA ELEGIR LA VERSION CORRECTA
# gff2bed to file conversion

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

current_ref <- '~/Documents/MIRNA_HALIOTIS/TRASH/GCF_023055435.1/ncbi_dataset/data/GCF_023055435.1/'
# prev_ref <- '~/Documents/MIRNA_HALIOTIS/GCF_003343065.1/ncbi_dataset/data/GCF_003343065.1/'

# review report at https://www.ncbi.nlm.nih.gov/genome/annotation_euk/Haliotis_rufescens/101/

file <- list.files(path = current_ref, pattern = "genomic.gff", full.names = T)

# library(Biostrings)
# library(rstatix)
library(tidyverse)

# https://edinburgh-genome-foundry.github.io/DnaFeaturesViewer/ (compare viz with this)

gff <- ape::read.gff(file)

names(gff)

head(gff)[1:8]

# attribute - A semicolon-separated list of tag-value pairs, providing additional information about each feature

# ID=exon-XM_046488722.2-3;
# Parent=rna-XM_046488722.2;
# Dbxref=GeneID:124125412,Genbank:XM_046488722.2;
# experiment=COORDINATES: polyA evidence [ECO:0006239];
# gbkey=mRNA;
# gene=LOC124125412;
# product=probable serine/threonine-protein kinase kinX;
# transcript_id=XM_046488722.2

into <- c("ID", "Parent", "Dbxref", "experiment", "gbkey", "gene", "product", "transcript_id")

nrow(gff) == length(unique(gff$seqid))

which_cols <- names(gff)[-9]

gff %>% as_tibble() %>%
  separate(attributes, into = into, sep = ";") %>%
  select(all_of(which_cols), Parent,transcript_id,product) %>%
  mutate(Parent = gsub("Parent=rna-","",Parent)) -> gff_rna
  
# gff_rna %>% mutate(transcript_id = gsub("transcript_id=","",transcript_id))

saveRDS(gff_rna, file = paste0(path, "gff_rna.rds"))

# For GCF_023055435.1 = 615 seqid (?) and 1,371,026 genomic features
# For GCF_003343065.1 = 8,371 (?) and 1,279,302 genomic features

gff %>% count(source) #%>% view()
gff %>% count(type) #%>% view()

gff %>% count(strand)

# json <- jsonlite::read_json('~/Documents/MIRNA_HALIOTIS/GCF_023055435.1/ncbi_dataset/data/assembly_data_report.jsonl')

# json$sourceDatabase

# https://cran.r-project.org/web/packages/chromoMap/vignettes/chromoMap.html

library(chromoMap)

# chromosome files

# chr_file_1 = arf[, c(6,8,9)]

chr.data <- gff_rna %>% filter(grepl("^X", Parent)) %>% select(Parent, start, end)

# chr_file_1 %>% distinct(Parent)

# chromosome name: a character representing the chromosome/contig/region name like ‘chr1’ or ‘1’ or ‘ch1’
# chromosome start: a numeric value to specify chromosome (or chromosome region) start position. If you are considering entire chromosome this value is typically 1.
# chromsome end: a numeric value specifying chromosome/contig/region end position. Again, if you are considering entire chromosome, then this value is the length of chromosome.
# centromere start (optional): centromeres will be added automatically if you provide the its start cordinates.

# annotation files

anno.data <- arf[,c(1,6,8,9, 11, 12)] 

names(anno.data) <- c("id", "Parent", "start", "end", "strand", "mismatch")

# annot.data <- anno.data %>% separate(id, into = c("Sample", "id"), sep = "_")

chr.features <- gff_rna %>% filter(grepl("^X", Parent)) %>% distinct(Parent, type, product)

anno.data <- anno.data %>% left_join(chr.features, by = "Parent")

# Element Name: a character specifying (uniquely) the elements. This can be identifiers,symbols etc.
# Chromosome Name: a character specifying the chromosome name. [NOTE: the chromosome names should be consistent in chromosome and data files.]
# Element Start: A numeric specifying element start position.
# Element End: A numeric specifying element end position.
# Data(optional): A numeric or character specifying the data value.
# Secondary Data(optional): A vector specifying the data value.useful for multi-factor scatter plots.
# Hyperlinks(optional): a character specifying the URL of the element.

write.table(chr.data, sep = "\t", file = paste0(path, "/chr.data.txt"))
write.table(anno.data, sep = "\t", file = paste0(path, "/ann.data.txt"))

ch_file <- paste0(path, "/chr.data.txt")
ann_file <- paste0(path, "/ann.data.txt")

# DIFICIL DE CARGAR, SI NO ES QUE IMPOSIBLE, PERO GRAFICOS MUY BONITOS PARA EXPLICAR LAS COORDENADAS
chromoMap(list(ch_file),list(ann_file))

# passing data.frames directly instead of files
chromoMap(list(chr.data),list(anno.data))
