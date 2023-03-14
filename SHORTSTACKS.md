# Introduction
Mirs loci, which produce the primary stem-loop transcripts, comprise only a small percentage of the source of the entire regulatory small RNA repertory in contrast to other type pf endogenous regulatory small RNAs. Based on the small RNA alignment patterns, shortstacks identifies and annotate both, miRs and non-miRs loci, and provides detailed descriptions of the small RNA populations emanating from each locus (S. Shaid and M. Axtell, 2013).

### SHORSTACK ERROR (REVISAR)

```bash
# fail to read the header from "-" 
# CREO QUE ES ERROR DEL ARCHIVO DE REFERENCIA O FASTA, TIENE CARACETERES ESPECIALE ADEMAS DEL HEADER
# PROBLEMA EN: Your header is too big to fit in a BAM file (2325494614 bytes, while the maximum officially allowed is 2147483648). samtools version 1.9 didn't error check this very well, which is why the conversion from sam apparently worked. The latest develop branch will say:

# SOLVED: 

tail /home/rvazquez/shortStacks_20230217.log
```

Este software implementa la caja de herramientas samtools. El uso de samtools presenta errores cuando los nombres de cabecera de los archivos fastx (Ej. `> my_sequence_name`) superan un maximo de bits permitido (Ej del error: fail to read the header from ... ). La solucion se describe [aqui](https://github.com/samtools/samtools/issues/1105#issuecomment-530300080)) y suguiere que _Your header is too big to fit in a BAM file (2325494614 bytes, while the maximum officially allowed is 2147483648_.Por lo anterior, se recomienda dar formato de entrada a los inputs , tanto archivos fasta/fastq como genoma de referencia:


### 1) FORMAT DATA INPUT

```bash
export PATH=$PATH:"/home/rvazquez/seqkit_tool"

# FOR SRNA-SEQS LIBS:

# RENAME HEADERS USING seqkit replace 
# (seqkit rename only rename duplicated IDs)
# seqkit replace -p .+ -r "seq_{nr}" --nr-width 5 # use --nr-width flag to choose n digits

# 0) TEST:

# fasta_file=HR110761.clean.fa

# cat $fasta_file | seqkit replace -p .+ -r "seq_{nr}" > ${fasta_file%.fa}.newid.fa

# 1) Loop

for i in $(ls *clean.fa); do cat $i | seqkit replace -p .+ -r "seq_{nr}" > ${i%.fa}.newid.fa; done

# 2) FOR GENOMIC FILES REMOVE EXTRA ID DESCRIPTION

reference=/home/rvazquez/GENOME_20230217/ncbi_dataset/data/GCF_023055435.1/GCF_023055435.1_xgHalRufe1.0.p_genomic.fa

cat $reference |  seqkit replace -p "\s.+" > ${reference%.fa}.newid.fa
# ${fasta_file%.fa}.newid.fa

# TO CONVERT FASTA TO FASTQ
../bbmap/reformat.sh ../HR110761.clean.fa out=file.fq

exit

seqkit fx2tab HR110761.clean.newid.fa | awk '{print $0,"EEEEEEEEEEEEEEEEEEEEE"}' | seqkit tab2fx > test_tab2fx.fq

```
### 2) Filtering reads based on its biotype (mirtrace)
```bash

# test

mirtrace qc -s meta_species_all *.fq -w --uncollapse-fasta --t 20

mirtrace qc -s meta_species_all -c configfile -w --uncollapse-fasta --t 20

# SUBSET BY BIOTYPE (EX.)
# grep ">"  configfile.output/qc_passed_reads.all.uncollapsed/HR110762.clean.newid.fasta  | awk '{print $2}' |  sort | uniq -c

#   275 rnatype:artifacts
#3402010 rnatype:mirna
#1034507 rnatype:rrna
#  98267 rnatype:trna
# 4097052 rnatype:unknown

export PATH=$PATH:"/home/rvazquez/seqkit_tool"

cd qc_passed_reads.all.uncollapsed/

for i  in $(ls *fasta); do cat $i | seqkit grep -n -r -p "rnatype:mirna" -p "rnatype:unknown" >  ${i%.fasta}.subset.fasta; done

# cat HR110762.clean.newid.fasta | seqkit grep -n -r -p "rnatype:mirna" -p "rnatype:unknown" > HR110762.clean.newid.subset.fasta




```
### 3) Profilling as a function of read length 
```bash
for i in $(ls *.clean.newid.fasta); do seqkit fx2tab $i -l -g -H > ${i%.fasta}.profiling; done &

# scp -r rvazquez@200.23.162.234:/home/rvazquez/MIRTRACE/configfile.output/qc_passed_reads.all.uncollapsed/PROFILING_BY_READ_LENGTH .

```

### X) (optional) Remove duplicates preserving their reads count

```bash
# USING SEQKIT

# -s, --by-seq by seq
# -i, --ignore-case ignore case
# -D, --dup-num-file string    file to save number and list of duplicated seqs
cat $fasta_file | seqkit rmdup -s -i -o unique.clean.fa.gz -D duplicated.detail.txt

# VSEARCH 

seqs=${fasta%.fasta}.good.fasta
sorted=${seqs%.fasta}_sorted.fasta
derep=${seqs%.fasta}_derep.fasta

vsearch  -sortbylength $seqs -output $sorted
vsearch -derep_fulllength $sorted -output $derep -sizeout

# OR TALLY
# Trimmed reads between 18-38 nt were collapsed to unique sequences using tally while preserving their counts from each sample (-l 18 -u 38 -format ‘>seq%I_w%L_x%C%n%R%n’). Sequences that were only observed 1 or 2 times across the 22 samples were discarded as they most likely represent sequencing errors, and in any case contain little useful information
```



```bash
git clone https://github.com/MikeAxtell/ShortStack.git

# cd ShortStack

# Dependencies (in PATH):

#samtools (version 1.x)
#bowtie (if aligning)
#bowtie-build (if aligning and .ebwt indices not found)
#gzip (if aligning)
#RNAfold (unless running with --nohp option to disable MIRNA search)


export PATH=$PATH:"/home/rvazquez/ShortStack"
export PATH=$PATH:"/home/rvazquez/Manatee/bowtie-1.0.1"
export PATH=$PATH:"/home/rvazquez/bedtools2/bin"
export PATH=$PATH:"//usr/local/bin/" # RNAfold from Vienna 
export PATH=$PATH:"/home/rvazquez/EMBOSS-6.6.0/emboss"
export PATH=$PATH:"/home/rvazquez/seqkit_tool"


# export PATH=$PATH:"/Users/cigom/Documents/MIRNA_HALIOTIS/ShortStack"
# readfile must be in fasta (.fasta or .fa), colorspace-fasta (.csfasta), or fastq (.fastq or .fq) format, or their gzip-compressed versions (.fasta.gz, .fa.gz, .csfasta.gz, .fastq.gz, or .fq.gz) Can also be a list (seperated by spaces) of several read files.

#    --adapter [string] : sequence of 3' adapter to trim off during read-pre
#    processing. Must be at least 8 bases, with only ATCG characters. If not
#    specified, reads are assumed to be already trimmed. Ex.  --adapter CTGTAGGC

# path to reference genome in .fasta or .fa format. Required for any run.

reference=/home/rvazquez/GENOME_20230217/ncbi_dataset/data/GCF_023055435.1/GCF_023055435.1_xgHalRufe1.0.p_genomic.newid.fa

# generate fai index

samtools faidx $reference -o ${reference}.fai

# ShortStack --readfile HR110761.clean.newid.fa HR110762.clean.newid.fa --outdir test1 --genomefile ${reference}.fai --bowtie_cores 12 2> "shortStacks_"$(date +%Y%m%d)".log" &

# configure (RUNNING)
# INPUTS: ~/MIRTRACE/configfile.output/qc_passed_reads.all.uncollapsed/SUBSET_OF_KNOWN_AND_UNKNOWN_READS
#
 
ShortStack --readfile  HR110761.clean.newid.subset.fasta HR110762.clean.newid.subset.fasta HR110763.clean.newid.subset.fasta HR11081.clean.newid.subset.fasta HR11082.clean.newid.subset.fasta HR11083.clean.newid.subset.fasta HR24761.clean.newid.subset.fasta HR24762.clean.newid.subset.fasta HR24763.clean.newid.subset.fasta HR2481.clean.newid.subset.fasta  HR2482.clean.newid.subset.fasta HR2483.clean.newid.subset.fasta --outdir ShortStack_"$(date +%Y%m%d)"_test --genomefile ${reference} --bowtie_cores 20 --sort_mem 60G --mismatches 0 --dicermax 30 --mmap u --mincov 1 --pad 1 2> "ShortStack_"$(date +%Y%m%d)".log" &

#ShortStack --readfile HR110761.clean.newid.fa HR110762.clean.newid.fa --outdir ShortStack_"$(date +%Y%m%d)"_test --genomefile ${reference} --bowtie_cores 20 --sort_mem 60G --mismatches 0 --dicermax 30 --mmap u --mincov 1 --pad 1 2> "ShortStack_"$(date +%Y%m%d)".log" &



exit

```
Additional installations:

```bash
# http://emboss.open-bio.org/html/adm/ch01s01.html
#You should use the URL:

ftp://emboss.open-bio.org/pub/EMBOSS/.
#The file you need for the EMBOSS base installation is
# EMBOSS-latest.tar.gz

gunzip emboss-latest.tar.gz
tar xvf emboss-latest.tar
cd EMBOSS-*
#Compiling. If you are compiling a fresh installation:
./configure
make
#If you compile it on subsequent occasions, use the following:-
# rm config.cache;make clean

export PATH=$PATH:"/home/rvazquez/EMBOSS-6.6.0/emboss"

multi_fasta=/home/rvazquez/GENOME_20230217/multi_genome.fa
outname=`basename ${multi_fasta%.fa}`
einverted -sequence $multi_fasta -outfile nw_${outname}.inv -outseq nw_${outname}.align -gap 12 -threshold 50 -match 3 -mismatch -4
   
```

## Dataviz shortstacks outputs



