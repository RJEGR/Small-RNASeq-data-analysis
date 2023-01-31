Ricardo Gomez-Reyes,

Markdown for microRNA analysis

January 2023



CHECK microRNA-encoded peptides???

**Table of content**

# Main
- [1) Make a quick view of the meta data:](#1--make-a-quick-view-of-the-meta-data-)
- [2) Run Trimmomatic](#2--run-trimmomatic)
- [3) Test genome-guide assembly (hisat2)](#3--test-genome-guide-assembly--hisat2-)
- [4) StringTie](#4--stringtie)
  * [4.1) Assembly quality](#41--assembly-quality)
  * [4.2) Transrate](#42--transrate)
- [5) Annotation](#5--annotation)
- [6) Quantification](#6--quantification)



## Installing tools

### Preprocessing

**fastqc**

```bash
# FastQC v0.11.7
fastqc -v
# test
fastqc *fastq.gz
```



**multiqc**

```bash
# MultiQC v1.10.1
multiqc -v
```



**mirTrace**

miRTrace is a QC(quality control) tool for small RNA sequencing data (sRNA-Seq). Each sample is characterized by profiling sequencing quality, read length, sequencing depth and miRNA complexity and also the amounts of miRNAs versus undesirable sequences (derived from tRNAs, rRNAs and sequencing artifacts). In addition to these routine QC analyses, miRTrace can accurately and sensitively resolve taxonomic origins of small RNA-Seq data based on the composition of clade-specific miRNAs.

```bash
conda install -c bioconda mirtrace

mirtrace -v # 1.0.1

mirtrace --help

```



**cutadapt** via python

https://anaconda.org/bioconda/cutadapt

```bash
# istead of: conda create -n cutadaptenv cutadapt

conda install -c bioconda cutadapt

# Esto mantendra cutadapt automaticamente instalado:

cutadapt -v

cutadapt --help
```



**Trimgalone**

```bash
# https://github.com/FelixKrueger/TrimGalore
# Install Trim Galore

curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.6.6.tar.gz -o trim_galore.tar.gz

tar xvzf trim_galore.tar.gz
# Run Trim Galore
ln -s /Users/cigom/Documents/Tools/TrimGalore-0.6.6/trim_galore .


# lease install Cutadapt first and make sure it is in the PATH, or specify the path to the Cutadapt executable using --path_to_cutadapt /path/to/cutadapt

```



Fastq_screen (as internal control)

Check https://rjegr.github.io/Transcriptomics/markdown/Processing

### Mirdeep2

Further readings: https://github.com/rajewsky-lab/mirdeep2

```bash
# 1) Download
curl -OJX GET "https://github.com/rajewsky-lab/mirdeep2/archive/refs/heads/master.zip" -H 
"Accept: application/zip"

# 2) Quick Installing

unzip master.zip
cd master
perl install.pl

# 3) Test
cd tutorial_dir
bash run_tut.sh

```



## Downloading files

### Genome reference

Updated version from (May 2022): https://www.ncbi.nlm.nih.gov/data-hub/genome/GCF_023055435.1/

**Assembly methods** **Sequencing technology**: PacBio Sequel IIe; Dovetail's OmniC; Illumina NovaSeq

```bash
curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_023055435.1/download?include_annotation_type=GENOME_GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA&filename=GCF_023055435.1.zip" -H "Accept: application/zip"
```



Previous version (Jul 2018): https://www.ncbi.nlm.nih.gov/data-hub/genome/GCF_003343065.1/

**Assembly methods** **Sequencing technology**: PacBio; Illumina HiSeq

```bash
curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_003343065.1/download?include_annotation_type=GENOME_GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA&filename=GCF_003343065.1.zip" -H "Accept: application/zip"
```



Then unzip

```bash
unzip *.zip
```

And check

```bash
ls GCF_0*/ncbi_dataset/data/GCF*
```



### miRNA database

The miRBase v. XX entries were retrieved from http://www.mirbase.org/ftp.shtml [[DOI](doi:10.1093/nar/gky1141)], MirGeneDB 2.0 data were downloaded from http://mirgenedb.org/download [[DOI](https://doi.org/10.1093/nar/gkz885)], and previously reported bivalve miRNAs were retrieved from the electronic supplementary material files of a number of papers

```bash
mkdir "MIRBASE_"$(date +%Y%m%d)

cd MIRBASE_*

curl -OJX GET "https://www.mirbase.org/ftp/CURRENT/hairpin.fa.gz" -H "Accept: application/zip"

curl -OJX GET "https://www.mirbase.org/ftp/CURRENT/mature.fa.gz" -H "Accept: application/zip"


gunzip *gz

grep -c "^>" *.fa

```

### MirGene DB

MirGeneDB is a database of manually curated microRNA genes that have been validated and annotated as initially described in [Fromm et al. 2015 ](http://www.annualreviews.org/doi/abs/10.1146/annurev-genet-120213-092023)and [Fromm et al. 2020](https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkz885/5584683?guestAccessKey=b9fe4339-a36b-49df-96ea-9e2cd7a1c99b). MirGeneDB 2.1 includes more than 16,000 microRNA gene entries representing more than 1,500 miRNA families from 75 metazoan species and published in the [2022 NAR database issue](https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkab1101/6439665). All microRNAs can be browsed, searched and downloaded.

miRNAs obtained may BLASTed (blastn) against miRBase and MirGeneDB databases.

```bash
mkdir "MIRGENEDB"$(date +%Y%m%d)

cd MIRGENEDB*
# https://mirgenedb.org/information
curl -OJX GET "https://mirgenedb.org/static/data/ALL/ALL-pre.fas" -H "Accept: application/zip"
curl -o ALL-mat.fas -OJX GET "https://mirgenedb.org/fasta/ALL?mat=1" 
curl -o ALL.fas -OJX GET "https://mirgenedb.org/fasta/ALL?all=1"

seqkit stats *fas

```



### snc-RNA Toy-set

Dataset available from NCBI SRA database (https://www.ncbi.nlm.nih.gov/sra, Rosani et al., 2021). 

```bash
# export PATH=/~/Downloads/sratoolkit.3.0.1-mac64/bin/:$PATH

# 0)

ln -s ~/Downloads/sratoolkit.3.0.1-mac64/bin/fastq-dump .

# 1) PRJNA470756, Insight into the differential smallRNA expression of Mediterranean mussel hemocytes stimmulated in vitro and its possible relation with the immune response
esearch -db sra -query "PRJNA470756" |  efetch -format docsum | xtract -pattern Runs -ACC @acc  -element "&ACC" > SraAccList.txt

# 2)

for i in $(cat SraAccList.txt); do echo $i | xargs -n 1 -P 1 ./fastq-dump --split-files --gzip --skip-technical --outdir .; done

# 3) (optional)

ls -lth *fastq.gz

seqkit stat *fastq.gz
```



# Pipeline

### Download and backup data

```bash
# Decompress file: gzip -d gz_file -c .
# Compress file: gzip file
# also compress and backup raw data
tar -czvf usftp21.novogene.com.tar.gz usftp21.novogene.com

# backup at
/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/novogene_hr_20230105/usftp21.novogene.com.tar.gz

```



### Preprocesing

Check the quality of libraries

```bash
# 
seqkit stat 00.CleanData/HR*/*gz
# Source	File	format	type	num_seqs	sum_len	min_len	avg_len	max_len
# 00.CleanData	HR110761	FASTA	DNA	16,637,068	270,665,973	10	16.3	42
# 00.CleanData	HR110762	FASTA	DNA	17,159,495	309,349,165	10	18	42
# ....
```

Also check the quality and find adapter from raw data (`01.RawData`)

```bash
# 0) Prepare work space

mkdir RAW_INPUT && cd RAW_INPUT

ln -s /Users/cigom/Documents/MIRNA_HALIOTIS/RAW_DATA/usftp21.novogene.com/01.RawData/HR*/*gz .

# 1) fastqc *fastq.gz

mkdir -p fastqc

# usftp21.novogene.com/00.CleanData/HR*/*gz

fastqc *gz -t 1 --nogroup -o ./fastqc

# 2)
multiqc ./fastqc/*zip -o ./multiqc

```

After detect sequence adapter, procesed to cutadapter.

## Trimming adapter sequences and filtering quality reads (from raw reads)

> Recordar que la ineficiencia de la ligacion (3' o 5') resulta en una baja complejidad de la bibliotecas de secuenciacion. Para algunas bibliotecas se pierde (usando los datos filtrados de novogene) mas del 40% de las lecturas debido a errores en los adaptores. Por lo anterior, intentaremos hacer una limpieza personal de los datos crudos usando los siguientes pasos:

```bash
# Small RNA adapte sequences used for the Haliotis rufescence larvae sequencing:

RNA 5’ Adapter (RA5), part:
5'-GTTCAGAGTTCTACAGTCCGACGATC-3'
RNA 3' Adapter (RA3), part:
5'-AGATCGGAAGAGCACACGTCT-3'
```

1) Run trimGalone to confirm the presence of the adapters in the raw libraries

```bash
# ./trim_galore SRR7145597_1.fastq

mkdir filtered_data && cd filtered_data

for i in $(ls ../*.gz); do ../../trim_galore $i ; done

# Automatic adapter sequences detection
for i in $(ls ../*.gz); do ../../trim_galore $i ; done

#  choosing specific adapter sequences:
# wd: /Users/cigom/Documents/MIRNA_HALIOTIS/RAW_INPUT/filtered_and_trimmed_AGATCGGAAGAGCACACGTCT

for i in $(ls ../*.gz); do ../../trim_galore $i -a AGATCGGAAGAGCACACGTCT --length 18 --max_length 30 -o filtered_data; done
```

TrimGalone detecta automaticamente el siguiente fragment `AGATCGGAAGAGC` en 12 las bibliotecas de secuencacion, ese fragmento corresponde a un fragmento del adaptor RNA 3' Adapter (RA3):

5'- AGATCGGAAGAGC - - - - - - - - -3'# dejando 8 nucleotidos sin flanquear

5'-AGATCGGAAGAGCACACGTCT-3' 

> Por lo tanto, deciframos que el fragmento del adaptor sRNA 3' es el que se encuentra presente en casi todas las bibliotecas secuenciadas. Revisar los archivos idividuales de fastqc para identificar que todas las `Overrepresented sequences` son variaciones de este adaptor `AGATCGGAAGAGCACACGTCT`



Instad of trimgalone use directly cutadapt (better)

```bash
for i in $(ls ../*.gz); do cutadapt $i -a AGATCGGAAGAGCACACGTCT --length 18 --max_length 30 -o filtered_data; done

for i in $(ls ../*.gz); do cutadapt $i -a AGATCGGAAGAGCACACGTCT --length 18 --max_length 30 -o filtered_data; done

for i in $(ls *.fq.gz); do cutadapt -a AGATCGGAAGAGCACACGTCT -m 18 -M 30 -o filtered_and_trimmed_AGATCGGAAGAGCACACGTCT/${i%.fq.gz}_trimmed.fq.gz $i; done
 
# --untrimmed-output FILE
# --minimum-length LENGTH or -m LENGTH: Discard processed reads that are shorter than LENGTH.
# --maximum-length LENGTH or -M LENGTH: Discard processed reads that are longer than LENGTH

# Using different conf m and M

for i in $(ls *.fq.gz); do cutadapt -a AGATCGGAAGAGCACACGTCT -m 15 -M 30 -o filtered_and_trimmed_AGATCGGAAGAGCACACGTCT/${i%.fq.gz}_trimmed.fq.gz $i; done &> cutadapt.log
```

Concat summary stats as follow

```bash
# to awk dynamic regexp for slash "/" use "\/"
grep 'Command line parameters:'  m_18_M_30_filter/cutadapt.log | awk '{gsub(/[filtered_and_trimmed_AGATCGGAAGAGCACACGTCT\/,fq.gz]/, "", $11); print $11}'

# get stats
grep -A16 '=== Summary ===' m_18_M_30_filter/cutadapt.log > m_18_M_30.summary

grep 'Total reads processed:' m_18_M_30.summary | awk '{print $4}'
grep 'Reads with adapters' m_18_M_30.summary | awk '{print $4}'

grep 'Reads that were too short:' m_18_M_30.summary | awk '{print $6}'
grep 'Reads that were too long:' m_18_M_30.summary | awk '{print $6}'
grep 'Reads written (passing filters):' m_18_M_30.summary | awk '{print $5}'

grep 'Total basepairs processed:' m_18_M_30.summary | awk '{print $4}'
grep 'Total written (filtered):' m_18_M_30.summary | awk '{print $4}'
 
```



2. Also we check the presence of adapter by running **global** aligment with maft

```bash
# /Users/cigom/Documents/MIRNA_HALIOTIS/RAW_INPUT/filtered_and_trimmed_AGATCGGAAGAGCACACGTCT
cp ../HR11082.fq.gz .
gzip -d HR11082.fq.gz

# seqkit sample -p 0.1 HR110763.fq| seqkit head -n 1000 > output.fa

# /Users/cigom/Documents/Tools/mirdeep2-master/bin/fastq2fasta.pl

fastq2fasta.pl HR11082.fq > HR11082.fa

seqkit sample -p 0.1 HR11082.fa | seqkit head -n 1000 > HR11082_subset.fa

mafft --globalpair HR11082_subset.fa > HR11082_subset.afa

# clustalo -i HR110763.fa -o HR110763.afa

#

seqkit sample -p 0.1 HR11082.clean.fa | seqkit head -n 1000 > HR11082.clean_subset.fa
mafft --globalpair HR11082.clean_subset.fa > HR11082.clean_subset.afa
```

> To inhouse convertion to fastq2fasta lets: `python3 -c "from Bio import SeqIO;SeqIO.convert('file.fq', 'file.fa', 'fasta')"`



3. Then, run mitrare in order to screening  the sRNA composition from the libraries as well as sanity check the adapter trimming

`qc` Quality control mode (full set of reports). --species must be given.

`--species`  (miRBase encoding). EXAMPLE: "hsa" for Homo sapiens. "meta_species_all" for <All species concatenated together>

`w` Write QC-passed reads and unknown reads (as defined in the RNA type plot) to the output folder. Identical reads are collapsed. Entries are sorted by abundance.

```bash
# Input: Fastq files
mirtrace --help

# mirtrace qc -s meta_species_all -w *trimmed.fq.gz
# WORDIR RAW_INPUT
mkdir MIRTRACE && cd MIRTRACE
mirtrace qc -s meta_species_all -a AGATCGGAAGAGCACACGTCT -w ../*fq.gz

# Then check the mirrace-report.html

```



## Reference genome indexing

```bash
#!/bin/bash

bowtie-build rna.fna rna

```

## miRNA analysis



If test with CLEANED libraries ()

```bash
mkdir CLEAN_INPUT

cp /Users/cigom/Documents/MIRNA_HALIOTIS/RAW_DATA/usftp21.novogene.com/00.CleanData/HR*/*.gz . 

gzip -d *gz

# Parse format

sanity_check_reads_ready_file.pl HR110761.clean.fa

# You could run remove_white_space_in_id.pl file.fa > newfile 
# This will remove everything from the id line after the first whitespac

for i in $(ls *.clean.fa); do remove_white_space_in_id.pl $i > ${i%.fa}.whitespac.fa; done

# processed w/ mapper.pl
```

Elaboramos los metadatos sobre las bibliotecas de lecturas preprocesadas. **Field 1 must be contain exactly three alphabet letters**

```bash
# Create config.txt file

ls * | grep whitespac.fa | sort > f2.tmp && cut -d '.' -f 1 f2.tmp > f1.tmp

paste f1.tmp f2.tmp | column -t > config.txt && rm *tmp

#  vi config.txt
HR1  HR110761.clean.fa.gz
HR2  HR110762.clean.fa.gz
HR3  HR110763.clean.fa.gz
...
```

#### 1. Mapper module

Considerando los siguientes flags y parametros, ejecutamos el primer modulo del framework mirdeep2

**Read input file**

-d input file is a config file

-e input file is fastq format (**not compressed**)

-j  remove all entries that have a sequence that contains letters other than a,c,g,t,u,n,A,C,G,T,U,N

-m collapse reads

-p genome map to genome (must be indexed by bowtie-build).

(Optional)

-h  parse to fasta format

-q  map with one mismatch in the seed (mapping takes longer)

-k seq  clip 3' adapter sequence (if neither trimming or filtering weremade)

**Output files**

-s file         print processed reads to this file
-t file         print read mappings to this file (.arf format)

```bash
#!/bin/bash
# If neither trimming and filtering were made:
mapper.pl config.txt -d -e -j -h -l 18 -k TCGTATGCCGTCTTCTGCTTGT -m -p rna -s reads_collapsed.fa 
-t reads_collapsed_vs_genome.arf -v 2> mapper.log

# 5 files in ~ 7 Gb in size (single-end 50 bp, 26,865,638 reads)
# Hora de inicio 6:50 p.m
# Hora final 8:15 p.m

# Otherwise:

# -v              outputs progress report
mapper.pl config.txt -d -e -j -h -l 18 -m -p rna -s reads_collapsed.fa -t reads_collapsed_vs_genome.arf -v 2> mapper.log

# If fasta format
index_path=/Users/cigom/Documents/MIRNA_HALIOTIS/INDEX/GCF_023055435/

mapper.pl config.txt -d -c -j -l 18 -m -p $index_path/rna -s reads_collapsed.fa -t reads_collapsed_vs_genome.arf -v 2> mapper.log


# Hora de inicio 3:55 p.m
# Hora final 4:15 p.m

```

Check the sequences length from `reads_collapsed.fa`

```bash
# fx2tabconvert FASTA/Q to tabular format, and provide various information, like sequence length, GC content/GC skew.
# https://bioinf.shenwei.me/seqkit/usage/#fx2tab-tab2fx

seqkit fx2tab --length --name --header-line  reads_collapsed.fa

# inf not seqkit, lets:

awk '/^>/{if (l!="") print l; print; l=0; next}{l+=length($0)}END{print l}' reads_collapsed.fa | head 

seqkit fx2tab --length --name --header-line  reads_collapsed.fa | awk '{print $2}'  | sort -n | uniq -c > reads_collapsed.length_table.list
```

#### 2. Quantify

This will result in final output with name `miRNAs_expressed_all_samples*.csv`

	-p precursor.fa  miRNA precursor sequences from miRBase
	-m mature.fa     miRNA sequences from miRBase

```bash
quantifier.pl -p ../MIRBASE_20230102/hairpin.fa -m ../MIRBASE_20230102/mature.fa -r reads_collapsed.fa
```

Then create the `Differential expression` workdirectory and link the `miRNAs_expressed_all_samples*` file there wich is an tsv file.

#### 3 (Optional) miRDeep2.pl (to detect novel miRNA candidates)

Basic functions to run mirdeep2 module

> miRDeep runtime:
>
> started: 22:1:58 
> ended: 9:31:23
> total:11h:29m:25s

```bash

miRDeep2.pl reads_collapsed.fa index.fa reads_collapsed_vs_genome.arf none none none 2> report.log

# Istead of "none" file Note that there it will in practice always improve miRDeep2 performance if miRNAs from some related species is input, even if it is not closely related.

**********
The first three arguments to miRDeep2.pl must be files while arguments 4-6 can be files or must be designated as 'none'. Examples:

miRDeep2.pl reads.fa genome.fa reads_vs_genome.arf mautre_ref_miRNAs.fa mature_other_miRNAs.fa  hairpin_ref_miRNAs

or

miRDeep2.pl reads.fa genome.fa reads_vs_genome.arf none none none


reads         deep sequences in fasta format. The identifier should contain a prefix, a running
              number and a '_x' to indicate the number of reads that have this sequence.
              There should be no redundancy in the sequences.
genome        genome contigs in fasta format. The identifiers should be unique.
mappings      file_reads mapped against file_genome. The mappings should be in arf format.
              For details on the format see the documentation.
miRNAs_ref    miRBase miRNA sequences in fasta format. These should be the known mature
              sequences for the species being analyzed.
miRNAs_other  miRBase miRNA sequences in fasta format. These should be the pooled known
              mature sequences for 1-5 species closely related to the species being
              analyzed.

precursors    miRBase miRNA precursor sequences in fasta format. These should be the known precursor
              sequences for the species being analyzed.

```

mirdeep2 runs the follow steps:

```bash
#Starting miRDeep2

#testing input files
#parsing genome mappings
#excising precursors
#preparing signature
#folding precursors
#computing randfold p-values
#running miRDeep core algorithm
#running permuted controls

#doing survey of accuracy
#producing graphic results

```



Additional function for the module

```bash
# Using c. elegans as reference

extract_miRNAs.pl mature.fa cel mature > mature_cel_sp.fa

# extract_miRNAs.pl hairpin.fa ref_specie > hairpin_ref_specie.fa

miRDeep2.pl reads_collapsed.fa rna.fa reads_collapsed_vs_genome.arf mature_cel_sp.fa mature.fa precursor.fa -t C.elegans

# miRDeep2.pl reads_collapsed.fa rna.fa reads_collapsed_vs_genome.arf mature_ref_this_species.fa mature_ref_other_species.fa precursors_ref_this_species.fa 2> report.log


#../../MIRBASE_20230102/mature_cel_sp.fa ../MIRBASE_20230102/mature.fa MIRBASE_20230102/hairpin.fa 2> report.log



```



Hay que definir que estrategia de miRNAs conocidos podemos analizar para la etpa de mirdeep2.pl y quantifier.pl

```bash


#ec=`echo $?`
#if [ $ec != 0 ];then
#        echo An error occured, exit code $ec
#fi


extract_miRNAs.pl mature.fa ref_specie mature > mature_ref_specie.fa && extract_miRNAs.pl hairpin.fa ref_specie > hairpin_ref_specie.fa
```

## Differential Expression Analysis



```R
library(tidyverse)

path <- '~/Documents/MIRNA_HALIOTIS/CLEAN_INPUT/'

mtd <- read_tsv(list.files(path = path, pattern = 'METADATA', full.names = T))

fileName <- list.files(path = path, pattern = 'miRNAs_expressed_all_samples', full.names = T)

count <- read_tsv(file = fileName)

```
