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
# PICARD (SAMMERGE), check java compativ.
git clone https://github.com/broadinstitute/picard.git
# git clone https://github.com/broadinstitute/picard.git --branch 2.25.0
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
After finishing
```bash
grep 'Input Read Pairs' trimmomatic.log  | awk '{print $4, $7, $12, $17, $20}' > reads_in.tmp

grep phred33 trimmomatic.log | awk '{gsub(/_/," ", $2); print $2}' | awk '{print $1}'  > samples_id.tmp

paste samples_id.tmp reads_in.tmp | column -t; rm *tmp
##
```
Good!
## Align reads to genome (HISAT2)

```bash
mkdir GENOME_INDEX;cd GENOME_INDEX


# ln -s /home/rvazquez/GENOME_20230217/ncbi_dataset/data/GCF_023055435.1/GCF_023055435.1_xgHalRufe1.0.p_genomic.newid.fa .

#genome=GCF_023055435.1_xgHalRufe1.0.p_genomic.newid.fa

cd /home/rvazquez/RNA_SEQ_ANALYSIS/ASSEMBLY

genome=/home/rvazquez/GENOME_20230217/ENSEMBLE/multi_genome.newid.fa

g_name=`basename ${genome%.fa}`

idx_file=/home/rvazquez/RNA_SEQ_ANALYSIS/GENOME_INDEX/$g_name

# hisat2-build -p 22 $genome $g_name &> hisat2_build.log &

# Run on paired-end reads

#Ex: 
# hisat2 --phred33 -p 8 -x ./GENOME_INDEX/$g_name -1 reads_f.fq -2 reads_r.fq -S output.sam


mkdir -p HISAT2_SAM_BAM_FILES

for file in $(ls *_R1.P.qtrim.fq.gz | grep gz)
do
withpath="${file}"
filename=${withpath##*/}
base="${filename%*_R*.P.qtrim.fq.gz}"
hisat2 --phred33 -p 2 -x $idx_file  \
    --rna-strandness RF \
    -1 ${base}_R1.P.qtrim.fq.gz -2 ${base}_R2.P.qtrim.fq.gz \
    --rg-id=${base} --rg SM:${base} -S HISAT2_SAM_BAM_FILES/${base}.sam
done &> hisat2.log &

# to kill the process related to it ps -o pid=pidid | xargs kill
# use ps to identify the loop process
```

## Reference guide assembly (Stringtie)

```bash
RNA_REF_GTF=/home/rvazquez/GENOME_20230217/SHORTSTACKS_REF/multi_genome.newid.gtf

RNA_ALIGN_DIR=/home/rvazquez/RNA_SEQ_ANALYSIS/ASSEMBLY/HISAT2_SAM_BAM_FILES/

WD=/home/rvazquez/RNA_SEQ_ANALYSIS/ASSEMBLY/STRINGTIE

cd $WD
ln -s $RNA_ALIGN_DIR/*.sam

# 1) sam to bam conversion and sorting

# StringTie takes as input a SAM, BAM or CRAM file sorted by coordinate (genomic location). This file should contain spliced RNA-seq read alignments such as HISAT2

for i in $(ls *sam)
do
samtools sort -@ 12 -o ${i%.sam}.sorted.bam $i
done

# 2)  run stringtie in reference-guided mode
# A reference annotation file in GTF or GFF3 format can be provided to StringTie using the -G option which can be used as 'guides' for the assembly process and help improve the transcript structure recovery for those transcripts.

# -l <label>       name prefix for output transcripts (default: MSTRG)

for i in $(ls *.sorted.bam)
do
withpath="${i}"
filename=${withpath##*/}
bs="${filename%*.sorted.bam}"
stringtie --rf -p 12 -G $RNA_REF_GTF -l $bs -o ${bs}_transcripts.gtf $i
done

# 2.1) (optional) run stringtie in de novo mode

for i in $(ls *.sorted.bam)
do
withpath="${i}"
filename=${withpath##*/}
bs="${filename%*.sorted.bam}"
stringtie --rf -p 12 -l $bs -o ${bs}_transcripts.gtf $i
done

```
## Output single transcriptome assembly

```bash
RNA_REF_GTF=/home/rvazquez/GENOME_20230217/SHORTSTACKS_REF/multi_genome.newid.gtf

RNA_REF_FASTA=/home/rvazquez/GENOME_20230217/SHORTSTACKS_REF/multi_genome.newid.fa

WD=/home/rvazquez/RNA_SEQ_ANALYSIS/ASSEMBLY/STRINGTIE

cd $WD

ls -1 *_transcripts.gtf > stringtie_gtf_list.txt

# awk '{if($3=="transcript") print}' stringtie_merged.gtf | cut -f 1,4,9 | less

# 1) Merge gtf
stringtie --rf --merge -p 24 -o transcripts.gtf -G $RNA_REF_GTF stringtie_gtf_list.txt

# 2) Generate a FASTA file with the DNA sequences for all transcripts in a GFF file

# -w flag write a fasta file with spliced exons for each transcript

gffread -w transcripts.fa -g $RNA_REF_FASTA transcripts.gtf

# 3) Optional

gffcompare -r $RNA_REF_GTF -o gffcompare transcripts.gtf
cat gffcompare.stats
```

## Output abundance matrix
```bash
# Note: 
# When the -e option is used, the reference annotation file -G is a required input and StringTie will not attempt to assemble the input read alignments but instead it will only estimate the expression levels of the "reference" transcripts provided in the -G file. With this option, no "novel" transcript assemblies (isoforms) will be produced, and read alignments not overlapping any of the given reference transcripts will be ignored, which may provide a considerable speed boost when the given set of reference transcripts is limited to a set of target genes for example.

# https://rnabio.org/module-05-isoforms/0005/05/01/Differential-Splicing/

# 1) Rerun strti using -e flag

for i in $(ls *.sorted.bam)
do
withpath="${i}"
filename=${withpath##*/}
bs="${filename%*.sorted.bam}"
stringtie --rf -p 12 -G transcripts.gtf -e -o ${bs}_transcripts.gtf $i
done

# stringtie --rf -p 4 -G transcripts.gtf -e -o ab_transcripts.gtf $RNA_ALIGN_DIR/UHR_Rep1.bam


https://github.com/RJEGR/Cancer_sete_T_assembly/blob/main/compile_trimmomatic.md#5-assembly-quality

# PREPARE DE PY
https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual

ls -l *_transcripts.gtf > samples.txt

#./sample1/sample1.gtf
#    ./sample2/sample2.gtf
#    ./sample3/sample3.gtf
    
python prepDE.py -i samples.txt

```

## notes
```bash
# 

# Merge HISAT2 BAM files
# 1) Usin samtools 
samtools merge output.bam *bam &> merge_bam.log &

#  2) 
PICARD=/home/rvazquez/RNA_SEQ_ANALYSIS/picard/build/libs/picard.jar

# java -jar $PICARD -h 
    
java -Xmx2g -jar $PICARD MergeSamFiles -OUTPUT UHR.bam -INPUT UHR_Rep1.bam -INPUT UHR_Rep2.bam -INPUT UHR_Rep3.bam

java -Xmx2g -jar $PICARD MergeSamFiles -OUTPUT output.bam \
-INPUT  SRR8956796.sorted.bam \
-INPUT  SRR8956797.sorted.bam \
-INPUT  SRR8956798.sorted.bam \
-INPUT  SRR8956799.sorted.bam \
-INPUT  SRR8956800.sorted.bam \
-INPUT  SRR8956801.sorted.bam \
-INPUT  SRR8956802.sorted.bam \
-INPUT  SRR8956803.sorted.bam \
-INPUT  SRR8956804.sorted.bam \
-INPUT  SRR8956805.sorted.bam


# (optional)

stringtie ${out_sam%.sam}.sorted.bam -o transcripts.gtf

# 2.1 Reference annotation transcripts 

# SOURCE THE GENOME GTG
ln -s /home/rvazquez/GENOME_20230217/ENSEMBLE/Haliotis_rufescens_gca023055435v1rs.xgHalRufe1.0.p.56.gtf genome.gtf


```