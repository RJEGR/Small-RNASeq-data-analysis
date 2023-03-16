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

# reference=/home/rvazquez/GENOME_20230217/ncbi_dataset/data/GCF_023055435.1/GCF_023055435.1_xgHalRufe1.0.p_genomic.fa

reference=/home/rvazquez/GENOME_20230217/ENSEMBLE/multi_genome.newid.fa

# cat $reference |  seqkit replace -p "\s.+" > ${reference%.fa}.newid.fa
# ${fasta_file%.fa}.newid.fa

# TO CONVERT FASTA TO FASTQ
../bbmap/reformat.sh ../HR110761.clean.fa out=file.fq

exit

seqkit fx2tab HR110761.clean.newid.fa | awk '{print $0,"EEEEEEEEEEEEEEEEEEEEE"}' | seqkit tab2fx > test_tab2fx.fq

```
### 2) Filtering reads based on its biotype (mirtrace)
Porque correr el filtrado de srna-seqs? While miRs are processed from hairpin containing precursor, tRNA are folded into characteristic L-shaped structure but both, tRNA and miRs are synthesized as a single-strand molecules than results in a rich source of rna fragments from similar size (Daniel Hasler, 2016) 
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
Hasta ahora ha dado resultados poco entendibles, dereplica una libreria de 3,098,492 (R2483.clean.newid.subset.fasta:3098492) en solamente 16 secuencias

```bash
# USING SEQKIT

# -s, --by-seq by seq
# -i, --ignore-case ignore case
# -D, --dup-num-file string    file to save number and list of duplicated seqs
cat $fasta_file | seqkit rmdup -s -i -o unique.clean.fa.gz -D duplicated.detail.txt

# TESTING BY MIRNAS READS
# VSEARCH (to profile single nucleotide difference between unique sequences) 

fasta_f=/home/rvazquez/MIRTRACE/configfile.output/qc_passed_reads.all.uncollapsed/ALL_UNCOLLAPSED_READS/HR2483.clean.newid.fasta

seqkit grep -n -r -p "rnatype:mirna" $fasta_f >  `basename ${fasta_f%.fasta}`.subset.fasta;

id=0.99
seqs=`basename ${fasta_f%.fasta}`.subset.fasta;
sorted=${seqs%.fasta}_sorted.fasta
derep=${seqs%.fasta}_derep.fasta
cntrds=${seqs%.*}_${id}_cnsns.tmp
acntrds=${cntrds%.*}_${min}_ab.fa

vsearch  -sortbylength $seqs -output $sorted --sizein

vsearch -derep_fulllength $sorted -output $derep -sizeout

vsearch  -sortbylength $udrp -output $udrp_srt
# 4.
# usearch -cluster_smallmem $udrp_srt -id  $id -centroids $cntrds -usersort -sizeout
# -consout instead of -centroids
vsearch --threads 20 -cluster_smallmem $sorted -id $id -consout $cntrds -usersort -sizeout

# save representative centroides
usearch -sortbysize $cntrds -minsize $min -output $acntrds

# Sequences distribution and size
grep '^>' $derep | cut -d"=" -f2 | sed 's/;$//g'| sort | uniq -c | sort -k2,2 -n > derep_distr.txt
grep '^>' $cntrds | cut -d"=" -f2 | sed 's/;$//g'| sort | uniq -c | sort -k2,2 -n > centroids_distr.txt


# OR TALLY
# Trimmed reads between 18-38 nt were collapsed to unique sequences using tally while preserving their counts from each sample (-l 18 -u 38 -format ‘>seq%I_w%L_x%C%n%R%n’). Sequences that were only observed 1 or 2 times across the 22 samples were discarded as they most likely represent sequencing errors, and in any case contain little useful information
```

### Running shortstacks

`-pad`: Los grupos de ARN pequeños encontrados inicialmente se fusionarán si la distancia entre ellos es menor o igual que el valor de pad. Debe ser un número entero entre 0 y 5000 0. Valor predeterminado: 75. Usamos pad=1 para recuperar formacion minima de clusteres


```bash
# VERSION 3.8
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
export PATH=$PATH:"/home/rvazquez/UCSCTOOLS/"
export PATH=$PATH:"/home/rvazquez/ShortTracks"

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

multi_fasta=/home/rvazquez/GENOME_20230217/SHORTSTACKS_REF/multi_genome.newid.fa
outname=`basename ${multi_fasta%.fa}`
einverted -sequence $multi_fasta -outfile nw_${outname}.inv -outseq nw_${outname}.align -gap 12 -threshold 50 -match 3 -mismatch -4
   
```

### Improving shorstacks (v 4.0)

El 12 de Marzo del 2023, se actualizo la version a 4.0, trayendo mejoras y facilidades en la instalacion y uso de la herramienta:

Using repeat inv file and nuclear+organelle genome (en la version superior a 3 ya no existen estos flags)

Installing newer version
```bash
# Shortstacks 4 requiered python >= 3.10.8, Therefore, lets upgrade it

sudo add-apt-repository ppa:deadsnakes/ppa
sudo apt update 
sudo apt install python3.11

# INSTALL PIP FOR THIS VERSION 
# https://pip.pypa.io/en/stable/installation/


export PATH=$PATH:"/home/rvazquez/.local/bin/"

pip3.11 install biopython
pip3.11 install tqdm
pip3.11 install tq

git clone https://github.com/MikeAxtell/ShortTracks.git

cd ShortTracks

# Replace ShortTracks header as follow:
# #!/usr/bin/env python3.11

# conda create --name strucvis strucvis

git clone https://github.com/MikeAxtell/ShortStack.git

mv ShortStack ShortStack4
cd ShortStacks4
mv ShortStack ShortStack4

vi ShortStack4
# Replace header as follow:
# #!/usr/bin/env python3.11
 
export PATH=$PATH:"/home/rvazquez/ShortStack4"
export PATH=$PATH:"/home/rvazquez/bedtools2/bin"
export PATH=$PATH:"//usr/local/bin/" # RNAfold from Vienna 
export PATH=$PATH:"/home/rvazquez/EMBOSS-6.6.0/emboss"
export PATH=$PATH:"/home/rvazquez/seqkit_tool"

export PATH=$PATH:"/home/rvazquez/bowtie-1.3.1-linux-x86_64"
# 
git clone https://github.com/MikeAxtell/strucVis.git
export PATH=$PATH:"/home/rvazquez/strucVis"

```

####  Running

```bash


export PATH=$PATH:"/home/rvazquez/UCSCTOOLS/"
export PATH=$PATH:"/home/rvazquez/ShortStack4"
export PATH=$PATH:"/home/rvazquez/bowtie-1.3.1-linux-x86_64"
export PATH=$PATH:"/usr/local/bin/"
export PATH=$PATH:"/home/rvazquez/strucVis"
export PATH=$PATH:"/home/rvazquez/UCSCTOOLS/"
export PATH=$PATH:"/home/rvazquez/ShortTracks"
export PATH=$PATH:"/home/rvazquez/bedtools2/bin"

REF=/home/rvazquez/GENOME_20230217/SHORTSTACKS_REF/multi_genome.newid.fa

# Evaluamos que no hay redundancia entre los ids y concatenamos en una sola base

MIRGENEDB=/home/rvazquez/MIRGENEDB_20230314/ALL-mat.fa
PIRBASE=/home/rvazquez/MIRTRACE/CUSTOM_DB/PIRBASE_V3/PIRBASE_rmdup.fasta

cat $MIRGENEDB $PIRBASE | seqkit rmdup -n -o knownsRNA.fa

# [INFO] 0 duplicated records removed

# generate fai index
# samtools faidx $reference -o ${reference}.fai


# TESTING VERSION 4.0

#../bbmap/reformat.sh HR110761.clean.newid.subset.fasta out=file.fq

# ShortStacks4 --genomefile $REF --knownRNAs $MIRGENEDB --dn_mirna --outdir ShortStack_"$(date +%Y%m%d)"_test --threads 20 --dicermax 30 --mmap u --mincov 0.8 --pad 1 --readfile file.fq
 
#  CONVERT ALL FASTA TO FASTQ

# for i in $(ls *newid.subset.fasta); do ../bbmap/reformat.sh $i out=${i%.fasta}.fq;done

readfile=`ls -x *.newid.subset.fq`

# --knownRNAs $MIRGENEDB

ShortStacks4 --genomefile $REF --knownRNAs knownsRNA.fa --dn_mirna --outdir ShortStack_"$(date +%Y%m%d)"_out --threads 24 --dicermax 30 --mmap u --mincov 0.8 --pad 1 --readfile $readfile &>> "ShortStack_"$(date +%Y%m%d)".log" &

# esperar 2 minutos a ver avances del archivo log

# Appending Standard Output and Standard Error: 2>&1. 
# https://www.cyberciti.biz/faq/how-to-redirect-standard-error-in-bash/
# Bash executes the redirects from left to right as follows:
# cmd >>file.txt 2>&1
# >>file.txt: Open file.txt in append mode and redirect stdout there.
# 2>&1: Redirect stderr to "where stdout is currently going". In this case, that is a file opened in append mode. In other words, the &1 reuses the file descriptor which stdout currently uses.

# to kill the process related to it ps -o pid=pidid | xargs kill
```

## Dataviz shortstacks outputs

/home/rvazquez/SHORTSTACKS/ShortStack_20230315_out


([M. Axtell Note](https://github.com/MikeAxtell/ShortStack#genome-browsers)): The output of ShortStack is designed to work with genome browsers. Specifically, the files Results.gff3, knownRNAs.gff3, the .bam files, and the .bw files can be directly visualized on either major genome browser (IGV, JBrowse). JBrowse2 has the ability to create "multi-wiggle" tracks. These tracks show multiple quantitative data tracks at once, bound to a common quantitative axis. The .bw bigwig files created by ShortStack & ShortTracks are normalized to reads-per-million, allowing direct comparisons in a multi-wiggle track. This allows visualization of size, coverage, and strandedness of the data. 

Types of bigwig files produced by ShortStack:

* readlength/stranded: A set of 8 .bw files with suffixes in the format _x_y.bw
  - x : A number ('21', '22', '23-24') or the string 'other', indicating the sizes of sRNAs tracked in that file.
  - y : Either 'p' or 'm' for the plus or minus genomic strand.
* readgroup: Only produced when there are multiple samples in the alignment. A set of n .bw files, with one per read-group. Because values are normalized to reads-per-million, these tracks are directly comparable to each other.

Load source 
```bash
export PATH=$PATH:"/Users/cigom/Documents/MIRNA_HALIOTIS/strucVis"
export PATH=$PATH:"/Users/cigom/Documents/MIRNA_HALIOTIS/ShortTracks"
export PATH=$PATH:"/Users/cigom/Documents/MIRNA_HALIOTIS/UCSCTOOLS" 
export PATH=$PATH:""
```
Then (if shortstacks is version =< 3.8)
```bash
ShortTracks --bamfile merged_alignments.bam --mode readlength --stranded 
```

Follow similar instructions than [here](https://jbrowse.org/jb2/docs/tutorials/config_gui/):

`Open JBroswer > create new assembly > show as linear > add track > multi-wiggle-track

Follow instructions to edit format: https://github.com/MikeAxtell/ShortTracks

Assembly display name: Haliotis rufescense (nuclei and organelle genome (GCA_023055435.1)
Assembly name: hr_multigenome_ensemble

consider https://gmod.github.io/JBrowseR/

```bash
printf "$PWD/%s\n" *.bam
```



