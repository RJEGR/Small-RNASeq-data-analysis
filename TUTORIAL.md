Ricardo Gomez-Reyes,

Markdown for microRNA analysis

January 2023

**Table of content**

[TOC]



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
/Users/cigom/Documents/Tools/TrimGalore-0.6.6
./trim_galore

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

```bash
# 1) fastqc *fastq.gz

mkdir -p fastqc

fastqc *fastq.gz -t 1 --nogroup -o ./fastqc

# 2)
multiqc ./*zip -o ./multiqc

```

After detect sequence adapter, procesed to mapper.pl or cutadapter. Also check https://www.biostars.org/p/328246/ or check what is the actual Illumina small rna seq 3' sequence adapter (https://www.biostars.org/p/69245/)



## Cutadapt

```bash
Small RNA adapte sequences:

RNA 5’ Adapter (RA5), part:
5'-GTTCAGAGTTCTACAGTCCGACGATC-3'
RNA 3' Adapter (RA3), part:
5'-AGATCGGAAGAGCACACGTCT-3'
```





## miRNA detection

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
... 
```

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

**Output files**

-s file         print processed reads to this file
-t file         print read mappings to this file (.arf format)

```bash
#!/bin/bash

bowtie-build rna.fna rna

#ec=`echo $?`
#if [ $ec != 0 ];then
#        echo An error occured, exit code $ec
#fi

mapper.pl config.txt -d -e -j -h -k TCGTATGCCGTCTTCTGCTTGT -m -p rna -s reads_collapsed.fa –t reads_collapsed_vs_genome.arf -v 2> mapper.log

# Hora de inicio 6:50 p.m
# Hora final

```

Hay que definir que estrategia de miRNAs conocidos podemos analizar para la etpa de mirdeep2.pl y quantifier.pl

```bash

#ec=`echo $?`
#if [ $ec != 0 ];then
#       echo An error occured, exit code $ec
#fi

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

