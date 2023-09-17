# Quality assess and clean step
# 1) Running tools
## Fastqc & Multiqc
```bash
# FastQC v0.11.7
fastqc -v
# test
fastqc *fastq.gz

# MultiQC v1.10.1
multiqc -v
```

## Mirtrace
```bash
mirtrace qc -s meta_species_all *.fq -w --uncollapse-fasta --t 20

# Subset reads by biotype:

cd qc_passed_reads.all.uncollapsed/

for i  in $(ls *fasta); do cat $i | seqkit grep -n -r -p "rnatype:mirna" -p "rnatype:unknown" >  ${i%.fasta}.subset.fasta; done
```


# Installation
## Mirtrace
```bash
conda install -c bioconda mirtrace
# OR conda install -c "bioconda/label/cf201901" mirtrace
mirtrace -v # 1.0.1

mirtrace --help
```