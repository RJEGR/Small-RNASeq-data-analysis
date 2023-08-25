# Curación y procesamiento de bases de datos
## Paso 1
```bash
VERSION=gca023055435v1rs
RELEASE=release-56

# A) DNA sequences
TYPE=dna
DATA=Haliotis_rufescens_${VERSION}.xgHalRufe1.0.p.dna.toplevel.fa.gz

URL=https://ftp.ebi.ac.uk/ensemblgenomes/pub/$RELEASE/metazoa/fasta/haliotis_rufescens_${VERSION}/dna/$DATA

wget $URL

# B) Sequence features (GTF)
TYPE=gtf
DATA=Haliotis_rufescens_gca023055435v1rs.xgHalRufe1.0.p.56.${TYPE}.gz

URL=https://ftp.ebi.ac.uk/ensemblgenomes/pub/$RELEASE/metazoa/${TYPE}/haliotis_rufescens_${VERSION}/$DATA

wget $URL

gunzip *.gz

```

## Paso 2
```bash
# A) Concat DNA sequences

nuclear_file=Haliotis_rufescens_gca023055435v1rs.xgHalRufe1.0.p.dna.toplevel.fa

mt_file=HR_mtDNA_whole_genome_shotgun_sequence_JALGQA010000616.1.fa 

cat $nuclear_file $mt_file | seqkit replace -p "\s.+" > multi_genome.newid.fa

# B) Concat Sequence features (GTF)

mtGTF=HR_mtDNA_whole_genome_shotgun_sequence_JALGQA010000616.1.gff3

cp Haliotis_rufescens_gca023055435v1rs.xgHalRufe1.0.p.56.gtf multi_genome.newid.gtf

tail -n1 $mtGTF >> multi_genome.newid.gtf
```
# miRs conocidos de referencia
```bash
mkdir "MIRGENEDB_"$(date +%Y%m%d)

cd MIRGENEDB*

url=https://mirgenedb.org

curl -o ALL-pre.fa -OJX GET "${url}/static/data/ALL/ALL-pre.fas"

curl -o ALL-mat.fa -OJX GET "${url}/fasta/ALL?mat=1" 

curl -o ALL.fa -OJX GET "${url}/fasta/ALL?all=1"

```
