# Transcriptome guide-assembly

PROTOCOLO NATURE PARA UN ENSAMBLE TRANSCRIPTOMICO GUIADO: Transcript-level expression analysis of RNA-seq experiments with HISAT, StringTie and Ballgown (https://doi.org/10.1038/nprot.2016.095)

Requiered tools:
```

HISAT2
https://github.com/DaehwanKimLab/hisat2#install
STRINGTIE
https://github.com/gpertea/stringtie#obtaining-and-installing-stringtie
RSEM
https://github.com/deweylab/RSEM#-compilation--installation
BUSCO
https://anaconda.org/bioconda/busco
Esearch  (For download ncbi libraries > 5Gb. Ex. scRNA-libs)
https://www.ncbi.nlm.nih.gov/books/NBK179288/
SRA toolkit
https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit
SeqKit (Optional)
https://bioinf.shenwei.me/seqkit/download/
FASTQC
https://anaconda.org/bioconda/fastqc
MULTIQC
conda install multiqc
TRIMMOMATIC
conda install -c bioconda trimmomatic

```

Export tools:
```bash
export PATH=$PATH:"/home/rvazquez/RNA_SEQ_ANALYSIS/hisat2"
export PATH=$PATH:"/home/rvazquez/RNA_SEQ_ANALYSIS/stringtie"
export PATH=$PATH:"/home/rvazquez/anaconda3/bin/"
```

## Exploratory RNA-seq quality
```bash
mkdir -p fastqc
fastqc *.fastq -t 24 --nogroup -o ./fastqc &> fastqc.log &

#
mkdir -p multiqc
multiqc ./fastqc/*zip -o multiqc

```
https://github.com/RJEGR/Cancer_sete_T_assembly/blob/main/compile_trimmomatic.md

Identificamos que es necesario remover la presencia del adaptor universal. Hay prencencia de N en entre secuencias pero la calidad de todas las biblitecas es superior a Phred > 30.

## Trimming adapters
Using trimmomatic. Parameter description is detailed [here](http://www.usadellab.org/cms/?page=trimmomatic)

```bash
wget https://raw.githubusercontent.com/RJEGR/Cancer_sete_T_assembly/main/TruSeq3-PE-2.fa
```

```bash

for file in $(ls *_1.fastq | grep fastq)
do
withpath="${file}"
filename=${withpath##*/}
base="${filename%*_*.fastq}"
trimmomatic PE -threads 20 -phred33 \
    ${base}_1.fastq ${base}_2.fastq \
    ${base}_1.P.qtrim.fq.gz ${base}_1.UP.qtrim.fq.gz \
    ${base}_2.P.qtrim.fq.gz ${base}_2.UP.qtrim.fq.gz \
    ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 SLIDINGWINDOW:4:15 MINLEN:36 LEADING:5 TRAILING:5
done &> trimmomatic.log &

```