# DATABASE CURATION
Genomic features and coordinates (Ex. Protein coding genes, tRNA, rRNA, ncRNAs classes, pseudogenes, etc.) may be taking from genomic GFF files. Recently (February, 2023), ENSEMBLE release the genomic features for Haliotis rufescense according to the genomic version 3055435v1rs

## ENSEMBLE:
A diferencia de NCBI, la version de ENSEMBLE tiene las anotaciones de UTR

**GTF:** Gene sets for each species. These files include annotations of both coding and non-coding genes. 
**GFF3:** GFF3 provides access to all annotated transcripts which make up an Ensembl gene set. 

> https://metazoa.ensembl.org/info/data/ftp/index.html

```bash
wget https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-56/metazoa/embl/haliotis_rufescens_gca023055435v1rs/Haliotis_rufescens_gca023055435v1rs.xgHalRufe1.0.p.56.nonchromosomal.dat.gz
```

## RNACENTRAL (snoRNA, snRNA,...)

```bash
# https://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/

wget https://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/sequences/rnacentral_active.fasta.gz
```

Son aprox 6.7 Gb de informacion comprimida (29 Gb descomp), contiene todas las secuencias de RNa de distintas bases de datos tales como: mirgenedb, mirbase, pirbase, snoRNA, rRNA, ensemble, etc ([aqui](https://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/sequences/by-database/) detalle). La anotacion se debe obtener del archivo `id_mappings`:

```bash
wget https://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/id_mapping/id_mapping.tsv.gz
# format Example
URS0000000001	ENA	GU786683.1:1..200:rRNA	77133	rRNA	
URS0000000001	ENA	GU786868.1:1..200:rRNA	77133	rRNA	
URS0000000001	ENA	GU786889.1:1..200:rRNA	77133	rRNA	
```

## MOLLUSC MIRS FROM GENOME_WIDE_ANALYSIS 
Esta base de referencia proviene del trabajo de Huang et al (2021), Tabla suplementaria Table_S1. Contiene distintas listas de miRs de 35 especies de moluscos. Un total de 5,823 secuencias unicas de miRs maduras y precursoras fueron recuperadas para su implementacion en la busqueda de miRs. Algunos de estos miRs fueron validados por steep-loop qRT-PCR usando U6-snRNA como referencia interna. El metodo de detección insilico se describe a continuacion: 

Para la preduccion de miRs maduros y precursores se consideró el siguiente criterio: 
1. Predicted mature miRs were allowed to have only 0 - 4 mismatch in sequences with know mature miRNA
2. The mismatched nucleotides were not permitted in know miRNA seed region (2 - 8 bp)
3. miRNA precursor can fold into an appropiate hairpin secondary structure that contains mature miRNA sequence whitin one arm of the hairpin and has the smallest possible folding energy.
Consultar el script

## MIRTRACE DATABASE CREATION
Revisar pagina 16 del manual

```bash
# /home/rvazquez/MIRTRACE/mirtrace/src/scripts/database_creation
# 1) export  the subdirectory /home/rvazquez/MIRTRACE/mirtrace/src/scripts/database_creation 

export PATH=$PATH:"/home/rvazquez/MIRTRACE/mirtrace/src/scripts/database_creation"

# 2) generate custom reference db than cat be used by mirtrace as follow (note: to run the script, wd must be the same dir as where the script resides)

# use python3 if neeeded

generate-mirtrace-rnatype-database.py -h

--out-dir $PWD
--species-abbrev my_abb
--species-verbosename my_sp_name
--mirna-seqs miRNA.fa
--rrna-seqs RRNA_SEQS
--trna-seqs TRNA_SEQS
--artifacts-seqs ARTIFACTS_SEQS

# emtpy fasta

# 3) RUNNING THE FILES

mirtrace qc --species my_abb --custom-db-folder $PWD ...

# aqui podriamos concatener la version mas reciente de mirbase + mirgenedb en el flag mirna-seqs, y asi sucesivamente en el flag rrna y trna incluir infomracion de rnacentral (split snorna, etc)

# revisar formatos fasta https://github.com/friedlanderlab/mirtrace/tree/master/src/lib/inputs



```

## NON REDUNDANT ANNOTATIONS
A non redundant (nr) genomic annotation mask may be build to reduce the Intra-range features (i.e overlapped annotations) and assign mapped sequences to unique/single annotations. Using proference (or hierarchical) annotation selection as follow: rRNA > tRNA, Transposable elements (TE) > protein-coding exon, other ncRNAs, introns, pseudogenes ... (Cei Abreu,2023). For this purpose, single sequence-annotation positions was obtained using the `findOverlaps` function
from `GenomicRanges` R package (Lawrence et al., 2013). In order to avoid
overlaps, we count according to the overlap of the central nucleotide of each read,
using the resize function, with parameter fix=”center” from the GenomicRanges
R package (Tesis de Isaac Martínez Ugalde).
