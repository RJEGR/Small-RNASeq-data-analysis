# Quality assess and clean step
# 1) Running tools
## Fastqc & Multiqc
```bash
# A) FastQC v0.11.7

mkdir -p fastqc

fastqc *.fastq -t 24 --nogroup -o ./fastqc &> fastqc.log &

# B) MultiQC v1.10.1

multiqc ./fastqc/*zip -o multiqc
```

## Mirtrace
```bash
# 1.0.1
mirtrace qc -s meta_species_all *.fq -w --uncollapse-fasta --t 20

# Subset reads by biotype:

cd qc_passed_reads.all.uncollapsed/

for i  in $(ls *fasta); do cat $i | seqkit grep -n -r -p "rnatype:mirna" -p "rnatype:unknown" >  ${i%.fasta}.subset.fasta; done
```


# Installation
## Fastqc
```bash
# FASTQC
conda install -c bioconda fastqc
# MULTIQC
conda install multiqc
```

## Mirtrace
```bash
conda install -c bioconda mirtrace
```

## Seqkit
```bash
conda install -c bioconda seqkit
```