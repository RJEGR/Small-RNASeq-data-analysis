Transcriptome guide-assembly

PROTOCOLO NATURE PARA UN ENSAMBLE TRANSCRIPTOMICO GUIADO: Transcript-level expression analysis of RNA-seq experiments with HISAT, StringTie and Ballgown (https://doi.org/10.1038/nprot.2016.095)

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

```

Export tools:
```bash
export PATH=$PATH:"/home/rvazquez/RNA_SEQ_ANALYSIS/hisat2"
export PATH=$PATH:"/home/rvazquez/RNA_SEQ_ANALYSIS/stringtie"
export PATH=$PATH:"/home/rvazquez/anaconda3/bin/"
```

PREPROCESSING READS
```bash
mkdir -p fastqc
fastqc *.fastq -t 24 --nogroup -o ./fastqc &> fastqc.log &

#
mkdir -p multiqc
multiqc ./fastqc/*zip -o multiqc

```
https://github.com/RJEGR/Cancer_sete_T_assembly/blob/main/compile_trimmomatic.md