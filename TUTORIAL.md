Ricardo Gomez-Reyes,

Markdown for microRNA analysis

January 2023

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



**trimmomatic**

...

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

### Preprocesing

Check the quality of libraries

```bash
seqkit stat *fastq.gz
file                format  type    num_seqs        sum_len  min_len  avg_len  max_len
SRR7145593_1.fastq  FASTQ   DNA   26,865,638  1,370,147,538       51       51       51
SRR7145594_1.fastq  FASTQ   DNA   27,317,373  1,393,186,023       51       51       51
SRR7145595_1.fastq  FASTQ   DNA   28,387,740  1,447,774,740       51       51       51
SRR7145596_1.fastq  FASTQ   DNA   26,807,946  1,367,205,246       51       51       51
SRR7145597_1.fastq  FASTQ   DNA   29,018,010  1,479,918,510       51       51       51

# After trimming and filtering

file                     format  type    num_seqs      sum_len  min_len  avg_len  max_len
SRR7145593_1_trimmed.fq  FASTQ   DNA   25,291,734  650,054,731       18     25.7       51
SRR7145594_1_trimmed.fq  FASTQ   DNA   25,787,594  683,537,844       18     26.5       51
SRR7145595_1_trimmed.fq  FASTQ   DNA   21,239,371  529,021,718       18     24.9       51
SRR7145596_1_trimmed.fq  FASTQ   DNA   21,594,261  528,528,837       18     24.5       51
SRR7145597_1_trimmed.fq  FASTQ   DNA   20,419,886  500,004,179       18     24.5       51
```

Then check the quality and find adapter

```bash
# 1) fastqc *fastq.gz

mkdir -p fastqc

fastqc *fastq.gz -t 1 --nogroup -o ./fastqc

# 2)
multiqc ./*zip -o ./multiqc

```

After detect sequence adapter, procesed to mapper.pl or cutadapter. Also check https://www.biostars.org/p/328246/ or check what is the actual Illumina small rna seq 3' sequence adapter (https://www.biostars.org/p/69245/)



## Trimming adapter sequences and filtering quality reads

> Recordar que la ineficiencia de la ligacion (3' o 5') resulta en una baja complejidad de la bibliotecas de secuenciacion.

```bash
Small RNA adapte sequences used for the Haliotis rufescence larvae sequencing:

RNA 5’ Adapter (RA5), part:
5'-GTTCAGAGTTCTACAGTCCGACGATC-3'
RNA 3' Adapter (RA3), part:
5'-AGATCGGAAGAGCACACGTCT-3'
```

Run trimGalone

```bash
# ./trim_galore SRR7145597_1.fastq

mkdir filtered_data && cd filtered_data

for i in $(ls ../*.fastq); do ../trim_galore $i ; done
```



## Reference genome indexing

```bash
#!/bin/bash

bowtie-build rna.fna rna

```

## miRNA analysis

According to mirdeep2

```bash
~/Documents/Tools/mirdeep2-master/tutorial_dir/run_tut.sh
```

Elaboramos los metadatos sobre las bibliotecas de lecturas preprocesadas. **Field 1 must be contain exactly three alphabet letters**

```bash
# Create config.txt file

ls *gz* | grep fastq | sort > f2.tmp && cut -d '.' -f 1 f2.tmp > f1.tmp

paste f1.tmp f2.tmp | column -t > config.txt && rm *tmp

#  vi config.txt
S93  SRR7145593_1.fastq
S94  SRR7145594_1.fastq
S95  SRR7145595_1.fastq
S96  SRR7145596_1.fastq
S97  SRR7145597_1.fastq
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
mapper.pl config.txt -d -e -j -h -k TCGTATGCCGTCTTCTGCTTGT -m -p rna -s reads_collapsed.fa 
-t reads_collapsed_vs_genome.arf -v 2> mapper.log

# 5 files in ~ 7 Gb in size (single-end 50 bp, 26,865,638 reads)
# Hora de inicio 6:50 p.m
# Hora final 8:15 p.m

# Otherwise:
# -v              outputs progress report
mapper.pl config.txt -d -e -j -h -m -p rna -s reads_collapsed.fa -t reads_collapsed_vs_genome.arf -v 2> mapper.log

# Hora de inicio 3:55 p.m
# Hora final 8:15 p.m

```

#### 2. Quantify

This will result in final output with name `miRNAs_expressed_all_samples*.csv`

	-p precursor.fa  miRNA precursor sequences from miRBase
	-m mature.fa     miRNA sequences from miRBase

```bash
quantifier.pl -p MIRBASE_20230102/hairpin.fa -m MIRBASE_20230102/mature.fa -r reads_collapsed.fa
```

#### 3 (Optional) miRDeep2.pl

```bash
miRDeep2.pl reads_collapsed.fa rna.fa reads_collapsed_vs_genome.arf mature_ref_this_species.fa mature_ref_other_species.fa precursors_ref_this_species.fa 2> report.log

miRDeep2.pl reads_collapsed.fa rna.fa reads_collapsed_vs_genome.arf MIRBASE_20230102/hairpin.fa MIRBASE_20230102/mature.fa 2> report.log
```



Hay que definir que estrategia de miRNAs conocidos podemos analizar para la etpa de mirdeep2.pl y quantifier.pl

```bash

miRDeep2.pl reads_collapsed.fa rna.fa reads_collapsed_vs_genome.arf mature_ref_this_species.fa mature_ref_other_species.fa precursors_ref_this_species.fa -t C.elegans

#ec=`echo $?`
#if [ $ec != 0 ];then
#        echo An error occured, exit code $ec
#fi


[rgomez@ixachi~]$: extract_miRNAs.pl mature.fa ref_specie mature >mature_ref_specie.fa
&& extract_miRNAs.pl hairpin.fa ref_specie >hairpin_ref_specie.fa
[rgomez@ixachi~]$: quantifier.pl -p precursors_ref_this_species.fa -m
mature_ref_this_species.fa -p precursos_ref_this_species.fa -r reads_collapsed.fa –t hsa -y
now
```

