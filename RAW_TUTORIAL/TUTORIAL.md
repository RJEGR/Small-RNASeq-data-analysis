Ricardo Gomez-Reyes,

Markdown for microRNA analysis

January 2023



CHECK microRNA-encoded peptides???
TRANSPONSABLE ELEMENTS DERIVED MIRNAS

### Server info

```bash
cat /etc/*-release | grep "DISTRIB" # see linux distro and version

#DISTRIB_ID=Ubuntu
#DISTRIB_RELEASE=20.04
#DISTRIB_CODENAME=focal
#DISTRIB_DESCRIPTION="Ubuntu 20.04.1 LTS"
# ID_LIKE=debian

free -h # check RAM memomy available

df -H # 

ssh rvazquez@200.23.162.234

```





## Installing tools

#### Python 

```bash
# Anaconda Dependencies
sudo apt-get install libgl1-mesa-glx libegl1-mesa libxrandr2 libxrandr2 libxss1 libxcursor1 libxcomposite1 libasound2 libxi6 libxtst6

# Download
wget https://repo.anaconda.com/archive/Anaconda3-2022.10-Linux-x86_64.sh

# install
bash Anaconda3-2022.10-Linux-x86_64.sh

# Include bioconda channel

conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

echo "export PATH=$PATH:/home/rvazquez/anaconda3/bin" >> $HOME/.bash_profile 


```
#### 
And also
#### 
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
# OR conda install -c "bioconda/label/cf201901" mirtrace
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



**Trimgalone (isteand of cutadapt)**

```bash
# https://github.com/FelixKrueger/TrimGalore
# Install Trim Galore

curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.6.6.tar.gz -o trim_galore.tar.gz

tar xvzf trim_galore.tar.gz
# Run Trim Galore
ln -s /Users/cigom/Documents/Tools/TrimGalore-0.6.6/trim_galore .


# lease install Cutadapt first and make sure it is in the PATH, or specify the path to the Cutadapt executable using --path_to_cutadapt /path/to/cutadapt

```



***fastq-screen***

Consider mapping to three genomes (A, B and C), the string '003' produces a file in which reads do not map to genomes A or B, but map (once or more) to genome C.  The string '--1' would generate a file in which reads uniquely map to genome C. Whether reads  map to genome A or B would be ignored.

```bash
conda install fastq-screen
conda update fastq-screen
# https://anaconda.org/bioconda/fastq-screen

fastq_screen --aligner bowtie2 --conf fastqscreen.conf --subset 0 --filter 3---- --tag --force --outdir ./fastq_screen/ 012-X04-P68-COI-ICT_S46_L001_R1_001.fastq.fastq
  
  
  for i in *qtrim.gz; do fastq_screen --aligner bowtie2 --conf fastqscreen.conf --subset 0 --tag --filter 0000330 --outdir ./fastq_screen/ $i; done
  
  # filter option is the good flag
  # 
```
Finally, modify the configure file (dowload a template [here](../examples/fastq_screen/fastq_screen.conf)) and then run the tool:


```



### Library manage

```bash
seqkit replace HR110761.clean.fa -p "^>" -r ''
# "^(\\S+)\\s?"

seqkit fx2tab *.fq.gz | head

fasta_file=HR110761.clean.fa
grep '^>' $fasta_file | sed 's/^>@//g' > ${fasta_file%.*}.ids & seqkit grep -n -f ${fasta_file%.*}.ids ${fasta_file%.clean.fa}.fq.gz -o ${fasta_file%.*}.fq


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

### (omit because doesnt work) DANSR: A tool for the detection of annotated and novel small RNAs

```bash
git clone https://github.com/ChrisMaherLab/DANSR.git

export PATH=/Users/cigom/Documents/Tools/DANSR/src/:$PATH

dansr -h

# Reference genome w/ bwa

index_path=/Users/cigom/Documents/MIRNA_HALIOTIS/INDEX/GCF_023055435_genomic

REFERENCE_FILE=$index_path/GCF_023055435.1_xgHalRufe1.0.p_genomic.fna

bwa index $REFERENCE_FILE

length=10

miRDeep2.pl reads_collapsed_l${length}.fa \
    $ref \
    reads_collapsed_l${length}_vs_${ref%.*}.arf  \
    none \
    none \
    none 2> "mirdeep2_"$length"_"$(date +%Y%m%d)".log"

# in order to create the ref gtf:
# rnacentral_active.fasta.gz: Current set of sequences that are present in at least one expert database.

curl -OJX GET "https://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/sequences/rnacentral_active.fasta.gz" # 7 Gb

bwa aln $ref RNAcentral.fa

INPUT_FILE=/Users/cigom/Documents/MIRNA_HALIOTIS/RAW_INPUT/filtered_data/HR110761_trimmed.fq.gz
PATH_TO_FEATURES=/Users/cigom/Documents/MIRNA_HALIOTIS/GENOME_20230210/ncbi_dataset/data/GCF_023055435.1

dansr \
  --input-file=$INPUT_FILE\
  --output-dir=example \
	--sample-name=example \
	--begin-no-trimming \
	--reference=$REFERENCE_FILE \
	--ref-gtf=$PATH_TO_FEATURES/genomic.gtf \
	--gtf-small=$PATH_TO_FEATURES/genomic.gtf
	
# 	--ref-gtf=$PATH_TO_REFERENCE/Homo_sapiens.GRCh38.105.gtf \
#	--gtf-small=reference/homo_sapiens.GRCh38.RNAcentral.gtf
	
	
```

### MANATEE (problemas de ejecucion debido a IntervalTree)

> corre bien con la cuenta rvazquez@200.23.162.234
>
> Pass: Temporal2022
>
> revisar manatee_20230217.log (error con el formato de annotation)

```bash

scp /Users/cigom/Documents/MIRNA_HALIOTIS/CLEAN_IN@ cPUT/HR110761.clean.fa rvazquez@200.23.162.234:/home/rvazquez/

git clone https://github.com/jehandzlik/Manatee.git

cpan install Set::IntervalTree

reference=/home/rvazquez/GENOME_20230217/ncbi_dataset/data/GCF_023055435.1/GCF_023055435.1_xgHalRufe1.0.p_genomic.fna

Manatee/bowtie-1.0.1/bowtie-build $reference ${reference%.fna}

cp config config.bkp

vi config

# configure annotation
#using GFFUtils
# https://github.com/fls-bioinformatics-core/GFFUtils
# Full documentation is available at http://gffutils.readthedocs.org/

source venv/bin/activate

gtf_extract /home/rvazquez/GENOME_20230217/ncbi_dataset/data/GCF_023055435.1/genomic.gtf --fields 'gene_name,gene_id,gene_biotype' -o genomic_subset.gtf

gtf_extract /home/rvazquez/GENOME_20230217/ncbi_dataset/data/GCF_023055435.1/genomic.gff --fields 'gene_name,gene_id,gene_biotype' --gff -o genomic_subset.gff

cut -f3 genomic_subset.gtf | sort | uniq -c

export PATH=/home/rvazquez/Manatee:$PATH

mkdir "manatee_"$(date +%Y%m%d)

# -i Path to pre-processed FASTQ or FASTA file. Valid formats: .fa, .fasta, .fastq, .fq, .fa.gz, .fasta.gz, .fastq.gz, .fq.gz.

manatee -config config -i HR110761.clean.fa -o "manatee_"$(date +%Y%m%d) 2> "manatee_"$(date +%Y%m%d)".log" &

--------------------------------MANDATORY INPUTS--------------------------------
################################################################################
Path and basename of the genome index to be searched. The basename is the name
of any of the index files up to but not including the final .1.ebwt/.rev.1.ebwt
/etc.q

index=/home/rvazquez/GENOME_20230217/ncbi_dataset/data/GCF_023055435.1/GCF_023055435.1_xgHalRufe1.0.p_genomic

############################## ANNOTATION ######################################
Path to annotation file in GFF3/GTF format. File should contain non coding
annotation.
# NOTE: Genomic annotation for ncRNAs is required as input in GTF format with the following tags in the attributes field: gene_name, gene_id, and gene_biotype.

annotation=genomic_subset.gtf

################################################################################
Path to reference FASTA/FA genome file.

genome=/home/rvazquez/GENOME_20230217/ncbi_dataset/data/GCF_023055435.1/GCF_023055435.1_xgHalRufe1.0.p_genomic.fna

# problemas or issues at:

cpan install ExtUtils::CppGuess # v ExtUtils-CppGuess-0.11
cpan install ExtUtils::MakeMaker
cpan install Set::IntervalTree


# install path 
/System/Library/Perl/5.30/ExtUtils
/Users/cigom/perl5/lib/perl5/ExtUtils/CppGuess.pm
export PERL5LIB=/Users/cigom/perl5/lib/perl5/ExtUtils/${PERL5LIB}
export PERL5LIB=/System/Library/Perl/5.30/ExtUtils${PERL5LIB}


export PERL5LIB=/Users/cigom/perl5/lib/perl5/darwin-thread-multi-2level/auto/ExtUtils:${PERL5LIB}
/Users/cigom/perl5/lib/perl5/darwin-thread-multi-2level/auto/ExtUtils
export 

export PATH=/Users/cigom/Documents/Tools/Manatee:$PATH

Manatee -h


```



### (OMIT) miRge 3.0 (problemas de instalación)

> Este paquete interesa para evaluar la entropia de isoformas de mirnas (`--isoform-entropy`) y edicion AtoI (`--AtoI`). Pero vale la pena mejor encontrar un algoritmo en R que se encargue de hacer este calculo: https://bioconductor.org/packages/release/bioc/vignettes/isomiRs/inst/doc/isomiRs.html

```bash
conda install -c bioconda mirge3
#wget https://www.python.org/ftp/python/3.7.5/python-3.7.5-macosx10.9.pkg
#sudo installer -pkg python-3.7.5-macosx10.9.pkg -target /
# Installing miRge3.0 with PyPi
python3.7 -m pip install --user cutadapt reportlab==3.5.42 biopython==1.78  scikit-learn==0.23.1  hypothesis==5.15.1 pytest==5.4.2  scipy==1.4.1  matplotlib==3.2.1  joblib==0.15.1  pandas==1.0.3 future==0.18.2

python3.7 -m pip install --user  mirge3
# upgrade
python3.7 -m pip install --user --upgrade  mirge3

export PATH=$PATH:"/Library/Frameworks/Python.framework/Versions/3.7/bin"

#
Output in format: Requested package -> Available versionsThe following specifications were found to be incompatible with your system:

  - feature:/osx-64::__osx==10.16=0
  - feature:|@/osx-64::__osx==10.16=0
  - mirge3 -> matplotlib-base -> __osx[version='>=10.11|>=10.12']

```

### (OMIT) Isomirna detection w miraaligner (from seqbuster > isomiRs fromLorena pantano) - error by using miraligner.jar (invalid or corrupt jarfile)

> https://seqcluster.readthedocs.io/mirna_annotation.html

```bash
# dependencies

# git clone https://github.com/ViennaRNA/ViennaRNA.git

curl -OJX GET "https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_5_x/ViennaRNA-2.5.1.tar.gz"

tar -xvzf ViennaRNA-2.5.1.tar.gz –C .

cd ViennaRNA*

./configure
make
sudo make install

# bedtools
# https://bedtools.readthedocs.io/en/latest/content/installation.html
curl -OJX GET "https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools-2.30.0.tar.gz"

#

curl -OJX GET "https://github.com/lpantano/seqbuster/raw/miraligner/modules/miraligner/miraligner.jar"

# desde mac darle acceso a permisos de ejecucion
#example:java -jar miraligner.jar -sub 1 -trim 3 -add 3 -s hsa -i test/test.fa -db DB -o test/out
#example: see output at miraligner/test/output.mirna & miraligner/test/output.mirna.opt


java -jar miraligner.jar -sub 1 -trim 3 -add 3 -s hsa -i ../HR110761.clean.fa -db ../MIRBASE_20230217  -o miraligner_out
```



### ilearnPlus (machine learning web method)

```bash
https://ilearnplus.erc.monash.edu/
machine-learning pipelines for computational analysis and predictions using nucleic acid and protein sequences.

https://github.com/Superzchen/iLearnPlus
```




### denovo mirna detection w brumiR (WORKING ON CLUSTER)

```bash
# Quick run
export PATH=$PATH:"~/BrumiR-3.0"

# TEST (WORKING!)
perl brumir.pl -a test/sRNA-seq.human.trim.fa.gz -p prefix


```



```bash
git clone https://github.com/camoragaq/BrumiR.git
cd ~/BrumiR-3.0

mkdir bin
cp src/* bin/

make all

# omit warning message

# add additional software:

wget https://github.com/GATB/bcalm/releases/download/v2.2.2/bcalm-binaries-v2.2.2-Linux.tar.gz

wget https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_4_x/ViennaRNA-2.4.14.tar.gz

tar -zxvf ViennaRNA-2.4.14.tar.gz

cd ViennaRNA-2.4.14.tar.gz


./configure

ln -s /usr/include/locale.h /usr/include/xlocale.h

make

sudo make install

# 

cp bcalm-binaries-v2.2.2-Linux/bin/bcalm ~/BrumiR-3.0/bin

ln -s /usr/local/bin/RNA* ~/BrumiR-3.0/bin


export PATH=$PATH:"~/BrumiR-3.0"

# TEST (WORKING!)

perl brumir.pl -a test/sRNA-seq.human.trim.fa.gz -p prefix

```







## Downloading files

### Genome reference
ENSEMBLE FEATURES: https://metazoa.ensembl.org/info/data/ftp/index.html
Annotation features: https://www.ncbi.nlm.nih.gov/genome/annotation_euk/Haliotis_rufescens/101/

Considere which file format to download:

```
**_genomic.fna.gz file:** FASTA format of the genomic sequence(s) in the assembly. Repetitive sequences in eukaryotes are masked to lower-case (see below). The FASTA title is formatted as sequence accession.version plus description. The genomic.fna.gz file includes all top-level sequences in the assembly (chromosomes, plasmids, organelles, unlocalized scaffolds, unplaced scaffolds, and any alternate loci or patch scaffolds). Scaffolds that are part of the chromosomes are not included because they are redundant with the chromosome sequences; sequences for these placed scaffolds are provided under the assembly_structure directory.
>       
> 
> *_genomic.gff.gz file
>        Annotation of the genomic sequence(s) in Generic Feature Format Version 3
>        (GFF3). Sequence identifiers are provided as accession.version.
>        Additional information about NCBI's GFF files is available at 
>        ftp://ftp.ncbi.nlm.nih.gov/genomes/README_GFF3.txt.
>        
> *_protein.faa.gz file
>        FASTA format sequences of the accessioned protein products annotated on
>        the genome assembly. The FASTA title is formatted as sequence 
>        accession.version plus description.
>        
> *_cds_from_genomic.fna.gz
>        FASTA format of the nucleotide sequences corresponding to all CDS 
>        features annotated on the assembly, based on the genome sequence. See 
>        the "Description of files" section below for details of the file format.
>        
>        
>  *_rna.fna.gz file
>        FASTA format of accessioned RNA products annotated on the genome 
>        assembly; Provided for RefSeq assemblies as relevant (Note, RNA and mRNA 
>        products are not instantiated as a separate accessioned record in GenBank
>        but are provided for some RefSeq genomes, most notably the eukaryotes.)
>        The FASTA title is provided as sequence accession.version plus 
>        description.
>   
>   *_rna_from_genomic.fna.gz
>        FASTA format of the nucleotide sequences corresponding to all RNA 
>        features annotated on the assembly, based on the genome sequence. See 
>        the "Description of files" section below for details of the file format.
>        
> Ref: https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/6454/101/GCF_023055435.1_xgHalRufe1.0.p/README.txt
> ```



Updated version from (May 2022): https://www.ncbi.nlm.nih.gov/data-hub/genome/GCF_023055435.1/

**Assembly methods** **Sequencing technology**: PacBio Sequel IIe; Dovetail's OmniC; Illumina NovaSeq

```bash
mkdir "GENOME_"$(date +%Y%m%d)
cd "GENOME_"$(date +%Y%m%d)
curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_023055435.1/download?include_annotation_type=GENOME_GFF,GENOME_GTF,GENOME_FASTA,RNA_FASTA,CDS_FASTA,PROT_FASTA&filename=GCF_023055435.1.zip" -H "Accept: application/zip"
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



And evaluate the genomic position (features)




### rfam database

```bash
mkdir "RFAM_"$(date +%Y%m%d)
cd RFAM_*
curl -OJX GET "https://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/rfam/rfam_annotations.tsv.gz"
curl -OJX GET "https://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/rfam/*.txt"
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


3. Then, run mitrace in order to screening  the sRNA composition from the libraries as well as sanity check the adapter trimming

`qc` Quality control mode (full set of reports). --species must be given.

`--species`  (miRBase encoding). EXAMPLE: "hsa" for Homo sapiens. "meta_species_all" for <All species concatenated together>

`w` Write QC-passed reads and unknown reads (as defined in the RNA type plot) to the output folder. Identical reads are collapsed. Entries are sorted by abundance.

```bash
# Input: Fastq files
mirtrace --help

# mirtrace qc -s meta_species_all -w *trimmed.fq.gz
# WORDIR RAW_INPUT
mkdir MIRTRACE && cd MIRTRACE
mirtrace qc -s meta_species_all -a AGATCGGAAGAGCACACGTCT -w --uncollapse-fasta ../*fq.gz

mirtrace qc -s meta_species_all -w --uncollapse-fasta ../*fq.gz

# Then check the mirrace-report.html


```



## Reference genome indexing

```bash
#!/bin/bash

# bowtie-build rna.fna rna
bowtie-build GCF_023055435.1_xgHalRufe1.0.p_genomic.fna GCF_023055435.1_xgHalRufe1.0.p_genomic



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

# Using complete genome as reference
# 1)
index_path=/Users/cigom/Documents/MIRNA_HALIOTIS/INDEX/GCF_023055435_genomic
ref=GCF_023055435.1_xgHalRufe1.0.p_genomic

length=15
mapper.pl config.txt -d -c -l $length -m -p $index_path/$ref -s reads_collapsed_l${length}.fa -t reads_collapsed_l${length}_vs_${ref%.*}.arf -v 2> "mapper_"$(date +%Y%m%d)".log"

# 2
index_path=/Users/cigom/Documents/MIRNA_HALIOTIS/INDEX/GCF_023055435_genomic
ref=GCF_023055435.1_xgHalRufe1.0.p_genomic

length=18
mapper.pl config.txt -d -c -l $length -m -p $index_path/$ref -s reads_collapsed_l${length}.fa -t reads_collapsed_l${length}_vs_${ref%.*}.arf -v 2> "mapper_l"$length"_"$(date +%Y%m%d)".log"

mv mapper_l18_20230210.log reads_collapsed_l18.fa reads_collapsed_l18_vs_GCF_023055435.1_xgHalRufe1.0.arf  TEST_l_18

# 3

index_path=/Users/cigom/Documents/MIRNA_HALIOTIS/INDEX/GCF_023055435_genomic
ref=GCF_023055435.1_xgHalRufe1.0.p_genomic

length=10

mapper.pl config.txt -d -c -l $length -m -p $index_path/$ref -s reads_collapsed_l${length}.fa -t reads_collapsed_l${length}_vs_${ref%.*}.arf -v 2> "mapper_l"$length"_"$(date +%Y%m%d)".log"

mv mapper_l10_20230210.log reads_collapsed_l10.fa reads_collapsed_l10_vs_GCF_023055435.1_xgHalRufe1.0.arf  TEST_l_10

tail -n14 mapper*.log | column -t


# 4) BECAUSE 18 IS THE BEST NUMBER OF DATA, LETS TURN ON -q flag in order to              allow mapping reads with one mismatch in the seed (mapping takes longer)

index_path=/Users/cigom/Documents/MIRNA_HALIOTIS/INDEX/GCF_023055435_genomic
ref=GCF_023055435.1_xgHalRufe1.0.p_genomic

length=18

mapper.pl config.txt -d -c -l $length -m -q -p $index_path/$ref -s reads_collapsed_l${length}_mismatch.fa -t reads_collapsed_l${length}_mismatch_vs_${ref%.*}.arf -v 2> "mapper_l"$length"_"$(date +%Y%m%d)".log"

tail -n14 mapper*.log | column -t

# (OMIT)
# -c              input file is fasta format
# -l discard reads shorter than int nts, default = 18
# LA VERSION DE -l 10 DEMORA MAS DE DOS DIAS EN CUANTIFICAR (quantification step)
index_path=/Users/cigom/Documents/MIRNA_HALIOTIS/INDEX/GCF_023055435/
mapper.pl config.txt -d -c -l 10 -m -p $index_path/rna -s reads_collapsed.fa -t reads_collapsed_vs_genome.arf -v 2> mapper.log


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
index_path=/Users/cigom/Documents/MIRNA_HALIOTIS/INDEX/GCF_023055435/

miRDeep2.pl reads_collapsed.fa $index_path/rna reads_collapsed_vs_genome.arf none none none 2> report.log

# Istead of "none" file. Note that there it will in practice always improve miRDeep2 performance if miRNAs from some related species is input, even if it is not closely related.

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

index_path=/Users/cigom/Documents/MIRNA_HALIOTIS/INDEX/GCF_023055435_genomic
ref=GCF_023055435.1_xgHalRufe1.0.p_genomic.fa

length=18
mautre_ref_miRNAs.fa mature_other_miRNAs.fa  hairpin_ref_miRNAs
mirbase_mat=/Users/cigom/Documents/MIRNA_HALIOTIS/MIRBASE_20230102/mature.fa
mirgene_mat=/Users/cigom/Documents/MIRNA_HALIOTIS/MIRGENEDB_20230130/ALL-mat.fas
mirgene_pre=/Users/cigom/Documents/MIRNA_HALIOTIS/MIRGENEDB_20230130/ALL-pre.fas

miRDeep2.pl reads_collapsed_l${length}.fa \
    $index_path/$ref \
    reads_collapsed_l${length}_vs_${ref%.*}.arf  \
    $mirbase_mat \
    $mirgene_mat \
    $mirgene_pre 2> "mirdeep2_"$length"_"$(date +%Y%m%d)".log"
    
#reads_collapsed_l${length}.fa
#reads_collapsed_l${length}_vs_${ref%.*}.arf
#2> "mirdeep2_"$length"_"$(date +%Y%m%d)".log"

miRDeep2.pl reads_collapsed_l${length}.fa \
    $index_path/$ref \
    reads_collapsed_l${length}_vs_${ref%.*}.arf  \
    none \ # mature_ref
    none \ # mature_other
    none 2> "mirdeep2_"$length"_"$(date +%Y%m%d)".log"
    
# or
  
miRDeep2.pl reads_collapsed_l${length}.fa $ref reads_collapsed_l${length}_vs_${ref%.*}.arf none mature.headers.fa hairpin_collapsed2dna.fa 2>report2.log


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



Hay que definir que estrategia de miRNAs conocidos podemos analizar para la etapa de mirdeep2.pl y quantifier.pl

```bash


#ec=`echo $?`
#if [ $ec != 0 ];then
#        echo An error occured, exit code $ec
#fi


extract_miRNAs.pl mature.fa ref_specie mature > mature_ref_specie.fa && extract_miRNAs.pl hairpin.fa ref_specie > hairpin_ref_specie.fa

```

### 4. MIRTOP

> Command line tool to annotate with a standard naming miRNAs e isomiRs.

```bash
https://github.com/miRTop/mirtop
```


