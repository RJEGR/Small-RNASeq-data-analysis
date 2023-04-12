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
## 0) Downloading dataset (NCBI utils)
```bash
sh -c "$(wget -q ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh -O -)"

echo "export PATH=\${PATH}:/home/rvazquez/edirect" >> ${HOME}/.bashrc

# To activate EDirect for this terminal session, please execute the following:

export PATH=${PATH}:${HOME}/edirect

# 1.2) https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/

sudo apt install ncbi-entrez-direct

wget https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/xtract.Linux.gz

ftp-cp ftp.ncbi.nlm.nih.gov /entrez/entrezdirect xtract.Linux.gz

gunzip -f xtract.Linux.gz

chmod +x xtract.Linux

mv xtract.Linux ~/.local/bin

#Add it to your path, or put it somewhere that is in your path, for example in ~/.local/bin/ so that you can get help by doing:

xtract.Linux -help

# RE-OPEN BASH SESSION

# 2) INSTALL sratoolkit FOR USE fastq-dump

wget --output-document sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz

tar -vxzf sratoolkit.tar.gz

echo "export PATH=$PATH:/home/rvazquez/RNA_SEQ_ANALYSIS/sratoolkit.3.0.1-ubuntu64/bin/" >> $HOME/.bash_profile 


which fasterq-dump # a replacement for the much older fastq-dump tool.

# 3) For best performance, obtain an API Key from NCBI, and place the following line in your .bash_profile and .zshrc configuration files:

# export NCBI_API_KEY=unique_api_key
echo "NCBI_API_KEY=74e0e4bf2d8eaa3bb742c46316dbafe12909" >> $HOME/.bash_profile

# Test

esearch -db sra -query "PRJNA488641" |  efetch -format docsum | xtract -pattern Runs -ACC @acc  -element "&ACC" > SraAccList.txt


# then download data
# Ex.
# fasterq-dump --split-files "SRR8956805" --skip-technical -p

for i in $(cat SraAccList.txt); do fasterq-dump --split-files $i --skip-technical -p; done &> fasterq_dump.log

# 

srapath SRR8956779

curl -OJX GET "https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR8956779/SRR8956779"

fasterq-dump --split-files SRR8956779


```
Additional tools
```bash
# Additional

INSTALL_DIR=/home/rvazquez/RNA_SEQ_ANALYSIS/stringtie/

cd $INSTALL_DIR

git clone https://github.com/gpertea/gffcompare
mv gffcompare gffcompare_dir; cd gffcompare_dir
make release

git clone https://github.com/gpertea/gffread
mv gffread gffread_dir; cd gffread_dir
make release

ln -s $INSTALL_DIR/gffread_dir/gffread
ln -s $INSTALL_DIR/gffcompare_dir/gffcompare
```

### Export tools:

```bash
export PATH=$PATH:"/home/rvazquez/RNA_SEQ_ANALYSIS/hisat2"
export PATH=$PATH:"/home/rvazquez/RNA_SEQ_ANALYSIS/stringtie"
export PATH=$PATH:"/home/rvazquez/anaconda3/bin/"
```

## 1) Exploratory RNA-seq quality
```bash
mkdir -p fastqc
fastqc *.fastq -t 24 --nogroup -o ./fastqc &> fastqc.log &

#
mkdir -p multiqc
multiqc ./fastqc/*zip -o multiqc

# 
export PATH=$PATH:"/home/rvazquez/seqkit_tool"
seqkit stat *fastq

```

Identificamos que es necesario remover la presencia del adaptor universal. Hay prencencia de N en entre secuencias pero la calidad de todas las biblitecas es superior a Phred > 30.

## 2) Trimming adapters
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
trimmomatic PE -threads 12 -phred33 \
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

## 3) Align reads to genome (HISAT2)

```bash
# mkdir GENOME_INDEX;cd GENOME_INDEX


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

cd /home/rvazquez/RNA_SEQ_ANALYSIS/ASSEMBLY

ln -s /home/rvazquez/RNA_SEQ_ANALYSIS/PROCESSED_LIBS/*P.qtrim.fq.gz

mkdir -p HISAT2_SAM_BAM_FILES

for file in $(ls *_1.P.qtrim.fq.gz | grep gz)
do
withpath="${file}"
filename=${withpath##*/}
base="${filename%*_*.P.qtrim.fq.gz}"
hisat2 --phred33 -p 8 -x $idx_file  \
    --rna-strandness RF \
    -1 ${base}_1.P.qtrim.fq.gz -2 ${base}_2.P.qtrim.fq.gz \
    --rg-id=${base} --rg SM:${base} -S HISAT2_SAM_BAM_FILES/${base}.sam
done 2> hisat2.log &

for i in $(ls *P.qtrim.fq.gz);do unlink $i; done

# to kill the process related to it ps -o pid=pidid | xargs kill
# use ps to identify the loop process
```

## 4) Reference guide assembly (Stringtie)

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

for i in $(ls *sam);do unlink $i; done

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

mkdir REFBASED_MODE; mv *_transcripts.gtf REFBASED_MODE

```
## 5) Denovo based assembly

```bash
# 1) (optional) run stringtie in de novo mode

mkdir DENOVO_MODE

for i in $(ls *.sorted.bam)
do
withpath="${i}"
filename=${withpath##*/}
bs="${filename%*.sorted.bam}"
stringtie --rf -p 2 -l $bs -o DENOVO_MODE/${bs}_transcripts.gtf $i
done

# 2) MERGE

ls -1 *_transcripts.gtf > stringtie_gtf_list.txt

stringtie --rf --merge -p 24 -o transcripts.gtf stringtie_gtf_list.txt

# 3) generate FASTA

RNA_REF_FASTA=/home/rvazquez/GENOME_20230217/SHORTSTACKS_REF/multi_genome.newid.fa


gffread -w transcripts.fa -g $RNA_REF_FASTA transcripts.gtf

grep "^>" -c transcripts.fa # 93,305 (OR 101,298 FROM 14 SAMPLES)

# awk '{if($3=="transcript") print}' transcripts.gtf | cut -f 1,4,9 | less

```

### 5.1) Output FASTA transcriptome assembly

Prepare guide assembly
```bash
RNA_REF_GTF=/home/rvazquez/GENOME_20230217/SHORTSTACKS_REF/multi_genome.newid.gtf

RNA_REF_FASTA=/home/rvazquez/GENOME_20230217/SHORTSTACKS_REF/multi_genome.newid.fa

WD=/home/rvazquez/RNA_SEQ_ANALYSIS/ASSEMBLY/STRINGTIE/REFBASED_MODE

cd $WD

# 1) Merge gtf

ls -1 *_transcripts.gtf > stringtie_gtf_list.txt 

stringtie --rf --merge -p 24 -o transcripts_guided.gtf -G $RNA_REF_GTF stringtie_gtf_list.txt

# 2) Get FASTA

gffread -w transcripts.fa -g $RNA_REF_FASTA transcripts.gtf

grep "^>" -c transcripts.fa # 134,293 (or 149,366 from 14 paired-end samples)

```

### GFF Compare:

```bash
REFBASED_MODE=/home/rvazquez/RNA_SEQ_ANALYSIS/ASSEMBLY/STRINGTIE/REFBASED_MODE/transcripts.gtf

DENOVO_MODE=/home/rvazquez/RNA_SEQ_ANALYSIS/ASSEMBLY/STRINGTIE/DENOVO_MODE/transcripts.gtf

RNA_REF_GTF=/home/rvazquez/GENOME_20230217/SHORTSTACKS_REF/multi_genome.newid.gtf

gffcompare -r $RNA_REF_GTF -o REFBASED.gffcompare $REFBASED_MODE

gffcompare -r $RNA_REF_GTF -o DENOVO_MODE.gffcompare $DENOVO_MODE

awk '{if($3=="transcript") print}' DENOVO_MODE.gffcompare.annotated.gtf | cut -f 1,4,9 | less

  
```

### 5.2) Output abundance matrix
For each RNA-Seq sample, run StringTie using the `-B/-b` and `-e` options in order to estimate transcript abundances and generate read coverage tables.

#### Generating guide based format
```bash
# Note: 
# When the -e option is used, the reference annotation file -G is a required input and StringTie will not attempt to assemble the input read alignments but instead it will only estimate the expression levels of the "reference" transcripts provided in the -G file. With this option, no "novel" transcript assemblies (isoforms) will be produced, and read alignments not overlapping any of the given reference transcripts will be ignored, which may provide a considerable speed boost when the given set of reference transcripts is limited to a set of target genes for example.


MODE=REFBASED_MODE

WD=/home/rvazquez/RNA_SEQ_ANALYSIS/ASSEMBLY/STRINGTIE/

cd $WD

mkdir -p QUANTIFICATION/$MODE

cd QUANTIFICATION/$MODE

ln -s $WD/*.sorted.bam .
ln -s $WD/$MODE/transcripts.gtf .


# 1) Generate *transcripts.gtf files with -e option!
# Note: Use the reference guided merged GTF in -G flag

for i in $(ls *.sorted.bam)
do
withpath="${i}"
filename=${withpath##*/}
bs="${filename%*.sorted.bam}"
stringtie --rf -p 24 -G transcripts.gtf -e -B -o ${bs}_eB_dir/${bs}_eB.gtf $i
done

# 1.1) Unlink unused files
for i in $(ls *.sorted.bam);do unlink $i; done
unlink transcripts.gtf

# https://github.com/RJEGR/Cancer_sete_T_assembly/blob/main/compile_trimmomatic.md#5-assembly-quality

# PREPARE DE PY
# https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual

# Download: https://ccb.jhu.edu/software/stringtie/dl/prepDE.py

# Generate csv sample file, Ex:
# ERR188021 <PATH_TO_ERR188021.gtf>
# ERR188023 <PATH_TO_ERR188023.gtf>
# ERR188024 <PATH_TO_ERR188024.gtf>
  
# considered subfolder ${bs}_eB_dir/${bs}_eB.gtf
for i in $(ls *_eB_dir/*_eB.gtf)
do
withpath="${i}"
fname=${withpath##*/}
bs="${fname%*_eB.gtf}"
echo "$bs" `printf "$PWD/${bs}_eB_dir/$fname"`
#echo "$bs" "$fname"
done > samples.txt

# GaD, continuar aqui 17/03/23

prepDE.py3 -i samples.txt  -v

# ..writing transcript_count_matrix.csv
# ..writing gene_count_matrix.csv
# All done.

# These count matrices (CSV files) can then be imported into R for use by DESeq2 and edgeR (using the DESeqDataSetFromMatrix and DGEList functions, respectively).

```

#### Generating De novo (Optional)
```bash
MODE=DENOVO_MODE

WD=/home/rvazquez/RNA_SEQ_ANALYSIS/ASSEMBLY/STRINGTIE/

cd $WD

mkdir -p QUANTIFICATION/$MODE

cd QUANTIFICATION/$MODE

ln -s $WD/*.sorted.bam .
ln -s $WD/$MODE/transcripts.gtf .

# 1) Generate *transcripts.gtf files with -e option!
# Note: Use the reference guided merged GTF in -G flag

# (Run as above)

```

Prepare assembly directory for file trasfer
```bash
MODE=REFBASED_MODE # DENOVO_MODE

FOLDER_PATH=/home/rvazquez/RNA_SEQ_ANALYSIS/ASSEMBLY/STRINGTIE/$MODE

tree -L 1 $FOLDER_PATH

tar -czvf ${MODE}.tar.gz $FOLDER_PATH

# 

scp rvazquez@200.23.162.234:/home/rvazquez/RNA_SEQ_ANALYSIS/ASSEMBLY/REFBASED_MODE.tar.gz 

scp rvazquez@200.23.162.234:/home/rvazquez/RNA_SEQ_ANALYSIS/ASSEMBLY/REFBASED_MODE.tar.gz .

scp rvazquez@200.23.162.234:/home/rvazquez/RNA_SEQ_ANALYSIS/ASSEMBLY/DENOVO_MODE.tar.gz .

# Extract

tar -xzvf REFBASED_MODE.tar.gz

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

### 6) Annotation
https://github.com/Trinotate/Trinotate/wiki
https://github.com/RJEGR/Transcriptomics/blob/master/markdown/trinotate.md
```bash

# you may need to install the DBD::SQLite module

perl -MCPAN -e shell
install DBD::SQLite

# Transdecoder

https://github.com/TransDecoder/TransDecoder/archive/refs/tags/TransDecoder-v5.7.0.tar.gz

tar xfv TransDecoder-v5.7.0.tar.gz
cd TransDecoder-TransDecoder-v5.7.0

chmod +x TransDecoder.LongOrfs



# sqlite
wget https://www.sqlite.org/2023/sqlite-tools-linux-x86-3410200.zip
unzip sqlite-tools-linux-x86-3410200.zip
cd sqlite-tools-linux-x86-3410200/

# HMMER/PFAM Protein Domain Identification:

wget http://eddylab.org/software/hmmer/hmmer.tar.gz
tar xvf hmmer.tar.gz
./configure
make
make install

# TmHMM
wget https://services.healthtech.dtu.dk/download/7e62eb60-c1be-404d-9d15-4bddb8600e93/tmhmm-2.0c.Linux.tar.gz

tar xvf tmhmm-2.0c.Linux.tar.gz
cd tmhmm-2.0c/bin
#Edit the header lines of the scripts `tmhmm` and `tmhmmformat.pl` to read exactly as:

#!/usr/bin/env perl

#Then, edit line 33 of tmhmm:
   #$opt_basedir = "/usr/cbs/packages/tmhmm/2.0c/tmhmm-2.0c/";
as:
   $opt_basedir = "/path/to/your/directory/containing/basedir/tmhmm-2.0c/"

#removing the comment '#' and setting the path to the directory where you installed the software.

# Emmapper: EggNOG-mapper is a tool for fast functional annotation of novel sequences.
git clone https://github.com/eggnogdb/eggnog-mapper.git

# trinotate:
wget https://github.com/Trinotate/Trinotate/archive/refs/tags/Trinotate-v4.0.0.tar.gz

tar -xvf Trinotate-v4.0.0.tar.gz

cd Trinotate-Trinotate-v4.0.0/

chmod +x Trinotate

./Trinotate --help

https://github.com/RJEGR/Transcriptomics/blob/master/markdown/trinotate.md


./util/admin/Build_Trinotate_Boilerplate_SQLite_db.pl Trinotate &> build.log &

```

#### Run
```bash
# then
WD=/home/rvazquez/Trinotate-Trinotate-v4.0.0/DATABASE
cd $WD
makeblastdb -in uniprot_sprot.pep -dbtype prot &
gunzip Pfam-A.hmm.gz &
srun $HMM/hmmpress Pfam-A.hmm &> hmmpress.log &

# 1) 
TR_PATH="/home/rvazquez/Trinotate-Trinotate-v4.0.0/"
TMHMM=$TR_PATH/"tmhmm-2.0c/bin"
SQLT=$TR_PATH"/sqlite-tools-linux-x86-3410200"
TRDC=$TR_PATH"/TransDecoder-TransDecoder-v5.7.0"
EGG=$TR_PATH"eggnog-mapper/" 

export PATH=$PATH:$TMHMM
export PATH=$PATH:$SQLT
export PATH=$PATH:$TRDC
export PATH=$PATH:$EGG
export PATH=$PATH:$TR_PATH
export PATH=$PATH:$TR_PATH/util

cp /home/rvazquez/Trinotate-Trinotate-v4.0.0/DATABASE/Trinotate.sqlite .

# optional (not working to me)
Trinotate_GTF_or_GFF3_annot_prep.pl \
    --annot Haliotis_rufescens_gca023055435v1rs.xgHalRufe1.0.p.56.gtf \
    --genome_fa transcripts.fa \
    --out_prefix out

# 2 Generate genes_trans_map

grep '^>' transcripts.fa |  sed 's/>//g' > transcript.map
cat transcript.map | awk '{gsub(/\.[0-9]$/,"",$1); print $1}' > genes.map
paste genes.map transcript.map > genes_trans_map
rm *.map


# 2.1 transdecoder first

TransDecoder.LongOrfs -t transcripts.fa --gene_trans_map genes_trans_map > transdecoder.log &

TransDecoder.Predict -t transcripts.fa --cpu 12 >> transdecoder.log & 

# 2.2) 
                
DB_DIR=/home/rvazquez/Trinotate-Trinotate-v4.0.0/DATABASE

mkdir -p $DB_DIR/EGGNOG_DATA_DIR

# run
WD=/home/rvazquez/RNA_SEQ_ANALYSIS/ASSEMBLY/STRINGTIE/REFBASED_MODE/ANNOTATION
cd $WD

Trinotate --db Trinotate.sqlite \
           --CPU 12 \
           --gene_trans_map genes_trans_map \
           --transcript_fasta transcripts.fa \
           --transdecoder_pep transcripts.fa.transdecoder.pep \
           --trinotate_data_dir $DB_DIR \
           --run "swissprot_blastp swissprot_blastx pfam tmhmmv2 EggnogMapper" &> Trinotate.log &

# /home/rvazquez/Trinotate-Trinotate-v4.0.0/PerlLib/Pipeliner.pm

# Run independently

blastp -query transcripts.fa.transdecoder.pep -db /home/rvazquez/Trinotate-Trinotate-v4.0.0/DATABASE/uniprot_sprot.pep -num_threads 12  -max_target_seqs 5 -outfmt 6 -evalue 1e-5 -out uniprot_sprot.ncbi.blastp.outfmt6 2> blastp.log &

# Pfam: identificación de dominio proteicos

hmmscan --cpu 12 -o hmmscan.out --domtblout TrinotatePFAM.out /home/rvazquez/Trinotate-Trinotate-v4.0.0/DATABASE/Pfam-A.hmm transcripts.fa.transdecoder.pep &

# Predicción de dominios transmembranal
tmhmm --short < transcripts.fa.transdecoder.pep > tmhmm.out &

# 
# transcripts.fa.transdecoder.cds

blastx -query transcripts.fa -db /home/rvazquez/Trinotate-Trinotate-v4.0.0/DATABASE/uniprot_sprot.pep -num_threads 1 -max_target_seqs 5 -outfmt 6 -evalue 1e-5 -out uniprot_sprot.ncbi.blastx.outfmt6 2> blastx.log &

# Warning: [blastx] Examining 5 or more matches is recommended
#terminate called after throwing an instance of 'std::__ios_failure'
#  what():  basic_ios::clear: iostream error
#Aborted

# evaluate /home/rvazquez/Trinotate-Trinotate-v4.0.0/util/auto_Trinotate.txt

# rnamer

RNAMMER_TRANS=/LUSTRE/bioinformatica_data/genomica_funcional/bin/Trinotate/util/rnammer_support

RNAMMER=/LUSTRE/bioinformatica_data/genomica_funcional/bin/rnammer/rnammer

srun $RNAMMER_TRANS/RnammerTranscriptome.pl --transcriptome good.Trinity.fasta --path_to_rnammer $RNAMMER &

```

#### LOAD
```bash

# THEN LOAD PROTEIN RESULTS
# Initial import of transcriptome and protein data:

# cp /home/rvazquez/Trinotate-Trinotate-v4.0.0/DATABASE/Trinotate.sqlite .
ad=/home/rvazquez/Trinotate-Trinotate-v4.0.0/util/admin/
$ad/Build_Trinotate_Boilerplate_SQLite_db.pl Trinotate
# CREATE

trinotate_data_dir=/home/rvazquez/Trinotate-Trinotate-v4.0.0/DATABASE

Trinotate --create --db myTrinotate.sqlite --trinotate_data_dir $trinotate_data_dir

sqlite=Trinotate.sqlite
gene_trans_map=genes_trans_map
transcript_fasta=transcripts.fa
transdecoder_pep=transcripts.fa.transdecoder.pep

Trinotate --db $sqlite --init --gene_trans_map $gene_trans_map --transcript_fasta $transcript_fasta --transdecoder_pep $transdecoder_pep

# Trinotate --db Trinotate.sqlite --init --gene_trans_map genes_trans_map --transcript_fasta transcripts.fa --transdecoder_pep transcripts.fa.transdecoder.pep

# Transdecoder loading protein search results:

Trinotate --db Trinotate.sqlite --LOAD_swissprot_blastp uniprot_sprot.ncbi.blastp.outfmt6

Trinotate --db Trinotate.sqlite --LOAD_pfam TrinotatePFAM.out
# Trinotate --db Trinotate.sqlite LOAD_tmhmm tmhmm.out

loaders_dir=/home/rvazquez/Trinotate-Trinotate-v4.0.0/util/trinotateSeqLoader/

$loaders_dir/Trinotate_BLAST_loader.pl --sqlite $sqlite --outfmt6 uniprot_sprot.ncbi.blastp.outfmt6 --prog blastp --dbtype Swissprot

# AND EXPORT

# Trinotate --db Trinotate.sqlite --report --incl_pep > Trinotate.xls

Trinotate --db $sqlite --report > Trinotate.xls

/home/rvazquez/Trinotate-Trinotate-v4.0.0/util/extract_GO_assignments_from_Trinotate_xls.pl --Trinotate_xls Trinotate.xls --gene > Trinotate_report.xls.gene_ontology

```