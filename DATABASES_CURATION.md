# DATABASE CURATION
Genomic features and coordinates (Ex. Protein coding genes, tRNA, rRNA, ncRNAs classes, pseudogenes, etc.) may be taking from genomic GFF files. Recently (February, 2023), ENSEMBLE release the genomic features for Haliotis rufescense according to the genomic version 3055435v1rs

## ENSEMBLE:
A diferencia de NCBI, la version GFF de ENSEMBLE tiene las anotaciones de UTR

**GTF:** Gene sets for each species. These files include annotations of both coding and non-coding genes. 
**GFF3:** GFF3 provides access to all annotated transcripts which make up an Ensembl gene set. 

> https://metazoa.ensembl.org/info/data/ftp/index.html

```bash
wget https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-56/metazoa/embl/haliotis_rufescens_gca023055435v1rs/Haliotis_rufescens_gca023055435v1rs.xgHalRufe1.0.p.56.nonchromosomal.dat.gz
```

## Generate multi-FASTA genome which include both, nuclear and organelle (mtDNA) genomes.
Unmasked version of genome should be used to allow the discovery of repeat-associated small RNAs, which are numerous in many species (S. Shaid and M. Axtell, 2013). We also concatenated both, mtDNA (JALGQA010000616.1) and nuclear genome (GCA_023055435.1) into one multi-FASTA file to allow comprehensive annotation/discovery of small RNAs-producing loci.  

```bash
nuclear_file=~/GENOME_20230217/ncbi_dataset/data/GCF_023055435.1/GCF_023055435.1_xgHalRufe1.0.p_genomic.newid.fa
mt_file=/home/rvazquez/GENOME_20230217/HR_mtDNA_whole_genome_shotgun_sequence_JALGQA010000616.1.fasta 
cat $nuclear_file $mt_file > multi_genome.fa
# grep -c "^>" multi_genome.fa
# 616 scaffolds

# /home/rvazquez/GENOME_20230217/multi_genome.fa

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

## MOLLUSC MIRS FROM GENOME_WIDE_ANALYSIS 
Esta base de referencia proviene del trabajo de Huang et al (2021), Tabla suplementaria Table_S1. Contiene distintas listas de miRs de 35 especies de moluscos. Un total de 5,823 secuencias unicas de miRs maduras y precursoras fueron recuperadas para su implementacion en la busqueda de miRs. Algunos de estos miRs fueron validados por steep-loop qRT-PCR usando U6-snRNA como referencia interna. El metodo de detección insilico se describe a continuacion: 

Para la preduccion de miRs maduros y precursores se consideró el siguiente criterio: 
1. Predicted mature miRs were allowed to have only 0 - 4 mismatch in sequences with know mature miRNA
2. The mismatched nucleotides were not permitted in know miRNA seed region (2 - 8 bp)
3. miRNA precursor can fold into an appropriate hairpin secondary structure that contains mature miRNA sequence within one arm of the hairpin and has the smallest possible folding energy.
Consultar el script [`GENOME_WIDE_MIRS_MOLLUSK.R`](https://github.com/RJEGR/Small-RNASeq-data-analysis/blob/master/GENOME_WIDE_MIRS_MOLLUSK.R) para revisión del procesamiento de dichas lecturas.

## piwi interacting RNAs (piRNAs)
The use of this db is motivated by the identification of reads > 26 nt in the library from embrio-larval stages. Here, PIRBASE V3 is used as standard gold to detect piRNAs reads. Prior to the annotation the use of `shortstacks` allow to  identify piRNA candidates based several sequences features. pirBASE_v3 is downloaded manually from here: http://bigdata.ibp.ac.cn/piRBase/download.php. In addition to standard golds We taking into account the presence of two mollusk > gastropods were piRNAs are identified: (Anaspidea -aca- and Biomphalaria glabrata -bgl-)

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

The overall high content of repetitive elements observed in molluscan genomes makes the first step of masking extremely important. The term repeated elements designate different nucleotide structures, from several low-complexity sequences to highly complex structures such as transposable elements.

Deseamos identificar elementos transponible s(Transposable elements , TE) del genoma del abulon usando RepeatMasker (Homology-based tool) o una combinacion de dos herramintas basadas en homologia (Repeatmasker, RepeatModeler) y prediccion de novo (Ej. LTR FINDER, RepeatScount).

**Homology-based prediction:** Homology-based approach involves searching commonly used databases for known TEs. We used RepeatProteinMask and RepeatMasker with repbase which contains a vast amount of known TEs.

`WD=/home/rvazquez/TRANSPOSABLE_ELEMENTS`

```bash
genome=/home/rvazquez/GENOME_20230217/multi_genome.fa

# /usr/local/RepeatMasker/Libraries
libdir=/usr/local/RepeatMasker/Libraries/
# -libdir $libdir
RepeatMasker -s -a -inv -gff -species "Haliotis rufescens" $genome &> RepeatMasker.log &

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

# test http://www.repeatmasker.org/webrepeatmaskerhelp.html

libdir=/usr/local/RepeatMasker/Libraries/

RepeatMasker -s -libdir $libdir -gff GCF_023055435.1_xgHalRufe1.0.p_genomic.newid.fa &> RepeatMasker.log &

# Note:
#Using Dfam with RepeatMasker
#============================

#RepeatMasker ships with a copy of Dfam (curated families only) in FamDB format. This can be replaced with a newer  version of Dfam, or with the full set of curated and uncurated families provided the famdb.py tool in the package is also updated.

pip3 install --user h5py
famdb.py -i dfam.h5 
#To use Dfam 3.7 with RepeatMasker 4.1.4 or earlier, first download an updated copy of the famdb.py tool from: https://github.com/Dfam-consortium/FamDB and replace the file in the RepeatMasker directory.  Then download the latest Dfam *.h5 file and rerun the RepeatMasker configure script.
```

then repeatModeler

```bash
# eledef
wget http://eddylab.org/software/recon/RECON1.05.tar.gz
gunzip RECON1.05.*gz; tar xvf RECON1.05.tar
cd RECON1.05/src
make
make install

cp imagespread eledef eleredef edgeredef famdef ../bin/
# /home/rvazquez/RECON1.05/bin

# 

wget 
http://www.repeatmasker.org/RepeatScout-1.0.6.tar.gz
gunzip RepeatScout.*; tar xvf RepeatScout-1.0.6.tar

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
```


## NON REDUNDANT ANNOTATIONS
A non redundant (nr) genomic annotation mask may be build to reduce the Intra-range features (i.e overlapped annotations) and assign mapped sequences to unique/single annotations. Using proference (or hierarchical) annotation selection as follow: rRNA > tRNA, Transposable elements (TE) > protein-coding exon, other ncRNAs, introns, pseudogenes ... (Cei Abreu,2023). For this purpose, single sequence-annotation positions was obtained using the `findOverlaps` function
from `GenomicRanges` R package (Lawrence et al., 2013). In order to avoid
overlaps, we count according to the overlap of the central nucleotide of each read,
using the resize function, with parameter fix=”center” from the GenomicRanges
R package (Tesis de Isaac Martínez Ugalde).
