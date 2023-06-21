# DATABASE CURATION
Genomic features and coordinates (Ex. Protein coding genes, tRNA, rRNA, ncRNAs classes, pseudogenes, etc.) may be taking from genomic GFF files. Recently (February, 2023), ENSEMBLE release the genomic features for _Haliotis rufescense_ according to the genomic version 3055435v1

## ENSEMBLE:
A diferencia de NCBI, la version GFF de ENSEMBLE tiene las anotaciones de UTR.

```bash
awk '{print $3}' genome.gtf| sort | uniq -c
```

**Notar que** el ENSEMBLE tiene la version `GCA_023055435.1` del genoma de _Haliotis rufescense_  con el prefijo GCA. Lo anterior indica que esta version es una copia del ensamble de GenBank. Por otro lado, la version ` GCF_023055435`, antes descargada, indica que la copia del ensamble proviene de RefSeq. Otras diferencias importantes son las [siguientes](https://computationalbiologybytes.wordpress.com/2015/10/02/clarifying-the-differences-between-the-grc-and-refseq-versions-of-h38/):

What is the difference between GCA$ and GCF $ ?

1. The GCA indicates the GenBank copy of the assembly, and the GCF indicates the RefSeq copy of the assembly.
2. The GenBank copy is the assembly that was provided by the submitter to GenBank. The RefSeq assembly is a copy of the GenBank copy that is used as the basis for RefSeq annotation. (This is a historical precedent: RefSeq does not annotate GenBank sequences, they only annotate RefSeq sequences: http://www.ncbi.nlm.nih.gov/books/NBK50679/).
3. Ex. The GenBank copy of the human reference assembly is devoid of annotation. The RefSeq copy contains the NCBI-provided annotation.
4. There are no sequence differences between GenBank and RefSeq versions of the equivalent human assembly release (e.g. GRCh38.p4).


What is the difference between GTF and GFF format ?

**GTF:** Gene sets for each species. These files include annotations of both coding and non-coding genes. 
**GFF3:** GFF3 provides access to all annotated transcripts which make up an Ensembl gene set. 

> https://metazoa.ensembl.org/info/data/ftp/index.html

```bash
wget https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-56/metazoa/embl/haliotis_rufescens_gca023055435v1rs/Haliotis_rufescens_gca023055435v1rs.xgHalRufe1.0.p.56.nonchromosomal.dat.gz

# TOPLEVEL (*.dna.toplevel.fa.gz): These files contains all sequence regions flagged as toplevel (toplevel sequences unmasked) in an Ensembl schema. This includes chromsomes, regions not assembled into chromosomes and N padded haplotype/patch regions.

# The header line in an FASTA dump files containing DNA sequence consists of the following attributes : coord_system:version:name:start:end:strand This coordinate-system string is used in the Ensembl API to retrieve slices with the SliceAdaptor.

wget https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-56/metazoa/fasta/haliotis_rufescens_gca023055435v1rs/dna/Haliotis_rufescens_gca023055435v1rs.xgHalRufe1.0.p.dna.toplevel.fa.gz


# These files hold the transcript sequences corresponding to non-coding RNA genes (ncRNA).

wget https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-56/metazoa/fasta/haliotis_rufescens_gca023055435v1rs/ncrna/Haliotis_rufescens_gca023055435v1rs.xgHalRufe1.0.p.ncrna.fa.gz

# 
wget https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-56/metazoa/fasta/haliotis_rufescens_gca023055435v1rs/cds/Haliotis_rufescens_gca023055435v1rs.xgHalRufe1.0.p.cds.all.fa.gz



```

## Prepare query for target prediction

**cds.all.fa.gz:** : These files hold the coding sequences corresponding to Ensembl genes. Note: CDS does not contain UTR or intronic sequence (Ex. `Haliotis_rufescens_gca023055435v1rs.xgHalRufe1.0.p.cds.all.fa.gz`)

**cdna.all.fa.gz**: These files hold the cDNA sequences corresponding to Ensembl genes, excluding ncRNA genes, which are in a separate 'ncrna' Fasta file. cDNA consists of transcript sequences for actual and possible genes, including pseudogenes, NMD and the like. See the file names  explanation below for different subsets of both known and predicted transcripts (ex. `Haliotis_rufescens_gca023055435v1rs.xgHalRufe1.0.p.cdna.all.fa.gz`)

```bash
grep -c "^>" *.cdna.all.fa # 57,675
grep -c "^>" *.cds.all.fa # 55,609
grep -c "^>" *.ncrna.fa # 12,345

awk '{print $3}' *.gtf| sort | grep "utr" | uniq -c

 # 77553 five_prime_utr
 # 54432 three_prime_utr
```
Subseq UTR sequences 
```bash
export PATH=$PATH:"/home/rvazquez/seqkit_tool"

WD=/home/rvazquez/GENOME_20230217/SHORTSTACKS_REF
gtf=$WD/multi_genome.newid.gtf
genome=multi_genome.newid.fa

# cat $gtf | grep "utr" | awk '{print $3}' | sort | uniq -c
# 77553 five_prime_utr
# 54432 three_prime_utr

cat $gtf | grep "utr" > utr.gtf

seqkit subseq --gtf utr.gtf $genome --gtf-tag "gene_id" -o utr.fa

cat $gtf | grep "three_prime_utr" > three_prime_utr.gtf

seqkit subseq --gtf three_prime_utr.gtf $genome --gtf-tag "gene_id" -o three_prime_utr.fa

# REMOVING DUPLICATED sequences (-s) in --only-positive-strand (-P)

cat utr.fa | seqkit rmdup -s -P -D duplicated_detail.txt -d duplicates.fasta -o utr_rmdup.fa

# [INFO] 49035 duplicated records removed

cat three_prime_utr.fa | seqkit rmdup -s -P -D duplicated_detail.txt -d duplicates.fasta -o three_prime_utr_rmdup.fa 

# [INFO] 23561 duplicated records removed

# FOR --gtf-tag choose one of the follow:
#gene_id "GeneID_125373053"; 
#transcript_id "GeneID_125373053_t1"; 
#gene_name "Trnar-ucg-67"; 
#gene_source "RefSeq"; 
#gene_biotype "tRNA"; 
#transcript_name "Trnar-ucg-67_t1"; 
#transcript_source "RefSeq"; 
#transcript_biotype "tRNA"; 
#tag "Ensembl_canonical";

/home/rvazquez/GENOME_20230217/ENSEMBLE/utr.fa

```

## Generate multi-FASTA genome which include both, nuclear and organelle (mtDNA) genomes.
Unmasked version of genome should be used to allow the discovery of repeat-associated small RNAs, which are numerous in many species (S. Shaid and M. Axtell, 2013). We also concatenated both, mtDNA (JALGQA010000616.1) and nuclear genome (GCA_023055435.1) into one multi-FASTA file to allow comprehensive annotation/discovery of small RNAs-producing loci.  

```bash
WD=/home/rvazquez/GENOME_20230217/SHORTSTACKS_REF
cd $WD
# nuclear_file=~/GENOME_20230217/ncbi_dataset/data/GCF_023055435.1/GCF_023055435.1_xgHalRufe1.0.p_genomic.newid.fa

nuclear_file=/home/rvazquez/GENOME_20230217/ENSEMBLE/Haliotis_rufescens_gca023055435v1rs.xgHalRufe1.0.p.dna.toplevel.fa

mt_file=/home/rvazquez/GENOME_20230217/HR_mtDNA_whole_genome_shotgun_sequence_JALGQA010000616.1.fasta 

# cat $nuclear_file $mt_file > multi_genome.fa
# grep -c "^>" multi_genome.fa
# 616 scaffolds

# FOR GENOMIC FILES REMOVE EXTRA ID DESCRIPTION

cat $nuclear_file $mt_file|seqkit replace -p "\s.+" > multi_genome.newid.fa


seqkit fx2tab --length --name multi_genome.newid.fa | head

# /home/rvazquez/GENOME_20230217/ENSEMBLE/multi_genome.newid.fa

# CURATE GTF
mtGTF=/home/rvazquez/GENOME_20230217/HR_mtDNA_whole_genome_shotgun_sequence_JALGQA010000616.1.gff3

cp /home/rvazquez/GENOME_20230217/ENSEMBLE/Haliotis_rufescens_gca023055435v1rs.xgHalRufe1.0.p.56.gtf multi_genome.newid.gtf


tail -n1 $mtGTF >> multi_genome.newid.gtf

# GFF3 (gene-region info)
mtGFF=HR_mtDNA_whole_genome_shotgun_sequence_JALGQA010000616.1.gff3

cp ./home/rvazquez/GENOME_20230217/ENSEMBLE/Haliotis_rufescens_gca023055435v1rs.xgHalRufe1.0.p.56.gff3.gz multi_genome.newid.gff3.gz
gunzip multi_genome.newid.gff3.gz
tail -n1 $mtGFF >> multi_genome.newid.gff3

```
 

## RNACENTRAL (snoRNA, snRNA,...)

```bash
# https://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/

wget https://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/sequences/rnacentral_active.fasta.gz
```

Son aprox 6.7 Gb de informacion comprimida (29 Gb descomp), contiene todas las secuencias de RNa de distintas bases de datos tales como: mirgenedb, mirbase, pirbase, snoRNA, rRNA, ensemble, etc ([aqui](https://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/sequences/by-database/) detalle). La anotacion se debe obtener del archivo `id_mappings`:

```bash
wget https://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/id_mapping/id_mapping.tsv.gz
# format Example
URS0000000001	ENA	GU786683.1:1..200:rRNA	77133	rRNA	
URS0000000001	ENA	GU786868.1:1..200:rRNA	77133	rRNA	
URS0000000001	ENA	GU786889.1:1..200:rRNA	77133	rRNA	
```

## MicroRNAs DB
### miRNA database

The miRBase v. 22 entries were retrieved from http://www.mirbase.org/ftp.shtml [[DOI](doi:10.1093/nar/gky1141)], MirGeneDB 2.0 data were downloaded from http://mirgenedb.org/download [[DOI](https://doi.org/10.1093/nar/gkz885)], and previously reported bivalve miRNAs were retrieved from the electronic supplementary material files of a number of papers

```bash
mkdir "MIRBASE_"$(date +%Y%m%d)

cd MIRBASE_*

curl -OJX GET "https://www.mirbase.org/ftp/CURRENT/hairpin.fa.gz" -H "Accept: application/zip"

curl -OJX GET "https://www.mirbase.org/ftp/CURRENT/mature.fa.gz" -H "Accept: application/zip"

curl -OJX GET "https://www.mirbase.org/ftp/CURRENT/miRNA.str.zip"  -H "Accept: application/zip"

gunzip *gz

grep -c "^>" *.fa

```

### MirGene DB

MirGeneDB is a database of manually curated microRNA genes that have been validated and annotated as initially described in [Fromm et al. 2015 ](http://www.annualreviews.org/doi/abs/10.1146/annurev-genet-120213-092023)and [Fromm et al. 2020](https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkz885/5584683?guestAccessKey=b9fe4339-a36b-49df-96ea-9e2cd7a1c99b). MirGeneDB 2.1 includes more than 16,000 microRNA gene entries representing more than 1,500 miRNA families from 75 metazoan species and published in the [2022 NAR database issue](https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkab1101/6439665). All microRNAs can be browsed, searched and downloaded.

miRNAs obtained may BLASTed (blastn) against miRBase and MirGeneDB databases.

```bash
mkdir "MIRGENEDB_"$(date +%Y%m%d)

cd MIRGENEDB*

# https://mirgenedb.org/information
curl -o ALL-pre.fa-OJX GET "https://mirgenedb.org/static/data/ALL/ALL-pre.fas"
curl -o ALL-mat.fa -OJX GET "https://mirgenedb.org/fasta/ALL?mat=1" 
curl -o ALL.fa -OJX GET "https://mirgenedb.org/fasta/ALL?all=1"

seqkit stats *fa

# TEST REMOVING DUPLICATED sequences (-s) in --only-positive-strand (-P)

cat ALL-mat.fa | seqkit rmdup -s -P -D duplicated_detail.txt -d duplicates.fasta -o ALL-mat-rmdup.fasta

# 11234 duplicated records removed

#La herramienta short stacks trabaja muy bien los duplicados Ej:

# En los resultados finales, tanto Cgi-Mir-216-P1_5p y Lgi-Mir-216-P1_5p se agrupan en un mismo cluster, ademas de tomar en cuenta a Esc-Mir-216-P1_5p que difiere en una base menos al final de la lectura 3p. 

#Por lo tanto omitimos usar el formato rmdup de esta DB

>Cgi-Mir-216-P1_5p
UAAUAUCAGCUGGUAAUCCUGAG
>Lgi-Mir-216-P1_5p
UAAUAUCAGCUGGUAAUCCUGAG
---
>Esc-Mir-216-P1_5p
UAAUAUCAGCUGGUAAUCCUGA

```



## MOLLUSC MIRS FROM GENOME_WIDE_ANALYSIS 
Esta base de referencia proviene del trabajo de Huang et al (2021), Tabla suplementaria Table_S1. Contiene distintas listas de miRs de 35 especies de moluscos. Un total de 5,823 secuencias unicas de miRs maduras y precursoras fueron recuperadas para su implementacion en la busqueda de miRs. Algunos de estos miRs fueron validados por steep-loop qRT-PCR usando U6-snRNA como referencia interna. El metodo de detección insilico se describe a continuacion: 

Para la preduccion de miRs maduros y precursores se consideró el siguiente criterio: 
1. Predicted mature miRs were allowed to have only 0 - 4 mismatch in sequences with know mature miRNA
2. The mismatched nucleotides were not permitted in know miRNA seed region (2 - 8 bp)
3. miRNA precursor can fold into an appropriate hairpin secondary structure that contains mature miRNA sequence within one arm of the hairpin and has the smallest possible folding energy.
Consultar el script [`GENOME_WIDE_MIRS_MOLLUSK.R`](https://github.com/RJEGR/Small-RNASeq-data-analysis/blob/master/GENOME_WIDE_MIRS_MOLLUSK.R) para revisión del procesamiento de dichas lecturas.

## piwi interacting RNAs (piRNAs)
The use of this db is motivated by the identification of reads > 26 nt in the library from embrio-larval stages. Here, PIRBASE V3 is used as standard gold to detect piRNAs reads. Prior to the annotation the use of `shortstacks` allow to  identify piRNA candidates based several sequences features. pirBASE_v3 is downloaded manually from here: http://bigdata.ibp.ac.cn/piRBase/download.php. In addition to standard golds (n=6 species) We taking into account the presence of two mollusk > gastropods were piRNAs are identified: (Anaspidea -aca- and Biomphalaria glabrata -bgl-), one crustacean (Scylla paramamosain - spa-), c. elegans (cel) and the Zebrafish (dre).

```bash
#cat *gold*.fa aca.fa bgl.v3.0.fa dre.fa > PIRBASE.fa

#  -l, --by-length by sequence length
#  -n, --by-name  by full name instead of just id
#  -s, --by-seq  by sequence   

fasta_f=PIRBASE_V3.fasta

cat *.fa | seqkit rmdup -s -P -D duplicated_detail.txt -d duplicates.fasta -o PIRBASE_rmdup.fasta

# [INFO] 2444 duplicated records removed

/home/rvazquez/MIRTRACE/CUSTOM_DB/PIRBASE_V3/PIRBASE_rmdup.fasta
```
Evaluamos tambien que no hay duplicado a nivel de ID y secuencia entre las bases de datos de miRNA y piRNA
```bash
MIRGENEDB=/home/rvazquez/MIRGENEDB_20230314/
PIRBASE=/home/rvazquez/MIRTRACE/CUSTOM_DB/PIRBASE_V3/PIRBASE_rmdup.fasta

# 1) [INFO] 0 duplicated records removed
cat $MIRGENEDB/ALL-mat-rmdup.fasta $PIRBASE | seqkit rmdup -s -P -o fa.tmp

# 2) [INFO] 0 duplicated records removed

cat $MIRGENEDB/ALL-mat.fa $PIRBASE | seqkit rmdup -n -o fa.tmp

rm fa.tmp

```


## MIRTRACE DATABASE CREATION
Revisar pagina 16 del manual

```bash
# /home/rvazquez/MIRTRACE/mirtrace/src/scripts/database_creation
# 1) export  the subdirectory /home/rvazquez/MIRTRACE/mirtrace/src/scripts/database_creation 

export PATH=$PATH:"/home/rvazquez/MIRTRACE/mirtrace/src/scripts/database_creation"

# 2) generate custom reference db than cat be used by mirtrace as follow (note: to run the script, wd must be the same dir as where the script resides)

# use python3 if neeeded

cd /home/rvazquez/MIRTRACE/mirtrace/src/scripts/database_creation

# WORKING MODE!!
dna_sequences=/home/rvazquez/MIRTRACE/CUSTOM_DB/molluscs_mature.fa

generate-mirtrace-rnatype-database.py \
  --out-dir $PWD \
  --species-abbrev moll \
  --species-verbosename meta_mollusck_db \
  --mirna-seqs $dna_sequences

# 2)
dna_sequences=/home/rvazquez/MIRTRACE/CUSTOM_DB/rnacentral_active.fasta

# WORKING MODE!!
generate-mirtrace-rnatype-database.py \
  --out-dir $PWD \
  --species-abbrev rnac \
  --species-verbosename rna_central_db \
  --mirna-seqs $dna_sequences

# 3)

dna_sequences=/home/rvazquez/MIRTRACE/CUSTOM_DB/PIRBASE.fa

generate-mirtrace-rnatype-database.py \
  --out-dir $PWD \
  --species-abbrev pirb \
  --species-verbosename pirna_gold_db \
  --mirna-seqs $dna_sequences
#  --artifacts-seqs rnacentral_active.fasta \
#  --rrna-seqs empty.fasta \
#  --trna-seqs empty.fasta

# emtpy fasta

# 3) RUNNING THE FILES
custom_db_dir=/home/rvazquez/MIRTRACE/mirtrace/src/scripts/database_creation

# /home/rvazquez/MIRTRACE/configfile.output/qc_passed_reads.all.uncollapsed/SUBSET_OF_KNOWN_AND_UNKNOWN_READS/*clean.newid.subset.fasta

fastq_path=/home/rvazquez/MIRTRACE/configfile.output/qc_passed_reads.all.uncollapsed/SUBSET_OF_KNOWN_AND_UNKNOWN_READS/

# moll db
mirtrace qc -s moll -w --uncollapse-fasta --custom-db-folder $custom_db_dir $fastq_path/**clean.newid.subset.fasta

# rna central

mirtrace qc -s rnac -w --uncollapse-fasta --custom-db-folder $custom_db_dir $fastq_path/**clean.newid.subset.fasta

# pirna db

mirtrace qc -s pirb -w --uncollapse-fasta --custom-db-folder $custom_db_dir $fastq_path/*.clean.newid.subset.fasta


# aqui podriamos concatener la version mas reciente de mirbase + mirgenedb en el flag mirna-seqs, y asi sucesivamente en el flag rrna y trna incluir infomracion de rnacentral (split snorna, etc)

# revisar formatos fasta https://github.com/friedlanderlab/mirtrace/tree/master/src/lib/inputs



```

## Repetitive landscape annotation
Deseamos identificar elementos transponible s(Transposable elements , TE) del genoma del abulon usando RepeatMasker (Homology-based tool) o una combinacion de dos herramintas basadas en homologia (Repeatmasker, RepeatModeler) y prediccion de novo (Ej. LTR FINDER, RepeatScount).

The overall high content of repetitive elements observed in molluscan genomes makes the first step of masking extremely important. The term repeated elements designate different nucleotide structures, from several low-complexity sequences to highly complex structures such as transposable elements. Similarly to TE-derived dsRNAs, sRNAs, including microRNAs (miRNAs), short-interfering RNAs (siRNAs) and PIWI-interacting RNAs (piRNAs), play central roles in regulating TEs185 . Some challenges are shared by both sRNA-seq analysis and mRNA-seq analysis, such as mapping ambiguity or quantification186 . However, sRNA-seq analysis in the context of repeated sequences has other specificities that are detailed by Nathan R. Johnson et al (2016)


**Homology-based prediction:** Homology-based approach involves searching commonly used databases for known TEs. We used RepeatProteinMask and RepeatMasker with repbase which contains a vast amount of known TEs.

### BWA  
https://mobilednajournal.biomedcentral.com/articles/10.1186/s13100-017-0086-z#Sec2
sRNA reads of 21 nt, 22 nt and 24 nt length and mRNA reads longer than 25 nt were mapped to the maize B73 genome (RefGen_V2) and the maize TE database using bwa with zero mismatches (‘bwa aln –n 0’). Because bwa places multiply mapping reads randomly onto one mapping location under the default setting, we selected ‘bwa samse –n 100000000’ to ensure that all alignments were reported


### Repeat masker

`WD=/home/rvazquez/TRANSPOSABLE_ELEMENTS`

```bash
genome=/home/rvazquez/GENOME_20230217/SHORTSTACKS_REF/multi_genome.newid.fa

# /usr/local/RepeatMasker/Libraries
libdir=/usr/local/RepeatMasker/Libraries/
# -libdir $libdir

RepeatMasker -s -a -inv -libdir $libdir -gff -species "all" $genome &> RepeatMasker.log &

RepeatMasker -gff -species "all" $genome &> RepeatMasker.log &

# -s  Slow search; 0-5% more sensitive, 2.5 times slower than default.
# -a      shows the alignments in a .align output file; -ali(gnments) also works
# -inv    alignments are presented in the orientation of the repeat (with option -a)

#-threads 20


# -norna Because of their close similarity to SINEs and the abundance of some of their pseudogenes, RepeatMasker by default screens for matches to small pol III transcribed RNAs (mostly tRNAs and snRNAs). When you're interested in small RNA genes, you should use the -norna option that leaves these sequences unmasked, while still masking SINEs

# 
famdb.py -i $libdir/RepeatMaskerLib.h5
# -species <query species>
#Specify the species or clade of the input sequence. The species name
#must be a valid NCBI Taxonomy Database species name and be contained
# in the RepeatMasker repeat database. Some examples are:

```
RepeatMasker Installation (take a while ...)
```bash
wget http://www.repeatmasker.org/RepeatMasker/RepeatMasker-4.1.4.tar.gz

sudo cp RepeatMasker*.gz /usr/local
cd /usr/local
sudo gunzip RepeatMasker*.gz
sudo tar xvf RepeatMasker*.tar

# 2) install dependencies
version=Dfam_3.7 # Dfam release 3.7 (January 2023)

wget https://www.dfam.org/releases/${version}/families/Dfam_curatedonly.h5.gz -O Dfam37_curatedonly.h5.gz

gunzip Dfam37_curatedonly.h5.gz

mv Dfam37_curatedonly.h5 /usr/local/RepeatMasker/Libraries

cd /usr/local/RepeatMasker/Libraries

sudo mv Dfam.h5 Dfam.h5.bak

sudo ln -s /usr/local/RepeatMasker/Libraries/Dfam37_curatedonly.h5 Dfam.h5

# NOTE: This will overwrite the distributed Dfam.h5 file.

# 3) INSTALL TANDEP REPEATS FINDER

git clone https://github.com/Benson-Genomics-Lab/TRF.git
cd TRF
sudo mkdir build
sudo ../configure
sudo make
sudo make install
# To copy binary elsewhere
cp /usr/local/TRF/build/src/trf /usr/local/RepeatMasker/util/trf

# INSTALL Sequence Search Engin (HMMER by default)

conda install -c bioconda hmmer

# 4)
cd /usr/local/RepeatMasker

echo "export PATH=$PATH:/home/rvazquez/anaconda3/bin" >> $HOME/.bash_profile 

sudo chmod 777 /usr/local/RepeatMasker

perl ./configure

export PATH:$PATH: # RepeatProteinMask

echo "export PATH=$PATH:/usr/local/RepeatMasker/" >> $HOME/.bash_profile 

libdir=/usr/local/RepeatMasker/Libraries/

# run
# -gff    creates an additional General Feature Finding format output

RepeatMasker -s -libdir $libdir -gff multi_genome.newid.fa &> RepeatMasker.log2 &

# the above abort because 
# cat: write error: No space left on device
# (PRELIMINAR RESULT) zcat /home/rvazquez/TRANSPOSABLE_ELEMENTS/RM_2073304.FriMar171337092023/multi_genome.newid.fa.cat.all.gz | less

RepeatMasker -s -libdir $libdir -gff JALGQA010000001_1.fasta &> RepeatMasker.log &


# INCLUDE ALL SPECIES 
cp RepeatMasker RepeatMasker.bkp

# Fix the issue by changing RepeatMasker line 7718 from this -->
if ( $lineage eq "") {

# To this -->
if ( $lineage eq "" && $species ne 'root' && $species ne 'all') {

# Because REGEX fails to account for negative numbers also replace line 8141 with:

if ( /^STATS\s+LOCAL\s+FORWARD\s+([\-\d\.]+)\s+([\d\.]+)/ ) {

# Ref https://github.com/rmhubley/RepeatMasker/issues/195
# Note:
#Using Dfam with RepeatMasker
#============================

#RepeatMasker ships with a copy of Dfam (curated families only) in FamDB format. This can be replaced with a newer  version of Dfam, or with the full set of curated and uncurated families provided the famdb.py tool in the package is also updated.

pip3 install --user h5py
famdb.py -i dfam.h5 
#To use Dfam 3.7 with RepeatMasker 4.1.4 or earlier, first download an updated copy of the famdb.py tool from: https://github.com/Dfam-consortium/FamDB and replace the file in the RepeatMasker directory.  Then download the latest Dfam *.h5 file and rerun the RepeatMasker configure script.
```

### repeatModeler
https://github.com/Dfam-consortium/RepeatModeler
https://github.com/4ureliek/Parsing-RepeatMasker-Outputs/blob/master/parseRM.pl

```bash
# eledef
wget http://eddylab.org/software/recon/RECON1.05.tar.gz
gunzip RECON1.05.*gz; tar xvf RECON1.05.tar
cd RECON1.05/src
make
make install

cp imagespread eledef eleredef edgeredef famdef ../bin/
# /home/rvazquez/RECON1.05/bin

# REPEATCOUNT

wget 
http://www.repeatmasker.org/RepeatScout-1.0.6.tar.gz
gunzip RepeatScout.*; tar xvf RepeatScout-1.0.6.tar

RC_DIR=/usr/local/RepeatScout-1.0.5
TRF_DIR=/usr/local/RepeatMasker/util/
CDHIT_DIR=/usr/bin/
UCSCTOOLS_DIR=/home/rvazquez/UCSCTOOLS
RMBLAST_DIR=/home/rvazquez/rmblast-2.13.0/bin
NINJA_DIR=home/rvazquez/NINJA-0.95-cluster_only/NINJA
# 
wget https://www.repeatmasker.org/rmblast/rmblast-2.13.0+-x64-linux.tar.gz

gunzip rmblast-2.13.0+-x64-linux.tar.gz
tar xvf rmblast-2.13.0+-x64-linux.tar
cd rmblast-2.13.0/

# LtrHarvest - The LtrHarvest program is part of the GenomeTools suite
wget http://genometools.org/pub/genometools-1.6.0.tar.gz
gunzip genometools-1.6.0.tar.gz
tar xvf genometools-1.6.0.tar

cd genometools*
make cairo=no threads=yes
make install

# LTR_RETRIEVER

conda install -c bioconda ltr_retriever 

# ninja

git clone https://github.com/TravisWheelerLab/NINJA.git
```


## NON REDUNDANT ANNOTATIONS
A non redundant (nr) genomic annotation mask may be build to reduce the Intra-range features (i.e overlapped annotations) and assign mapped sequences to unique/single annotations. Using proference (or hierarchical) annotation selection as follow: rRNA > tRNA, Transposable elements (TE) > protein-coding exon, other ncRNAs, introns, pseudogenes ... (Cei Abreu,2023). For this purpose, single sequence-annotation positions was obtained using the `findOverlaps` function
from `GenomicRanges` R package (Lawrence et al., 2013). In order to avoid
overlaps, we count according to the overlap of the central nucleotide of each read,
using the resize function, with parameter fix=”center” from the GenomicRanges
R package (Tesis de Isaac Martínez Ugalde).
