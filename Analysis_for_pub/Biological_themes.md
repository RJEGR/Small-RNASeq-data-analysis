# Purpose
1) Define biological-phentypic  themes

+ Calcification
+ Development
+ Growth
+ RM (Respiratory Metabolism)

2) Group transcriptome to these themes
+ Ideally, using only miRNA:mRNA transcripts (~ 170 transcripts)

3) Additionally, blast miRNA:mRNA transcripts to biomineralization-genes database

```R
# Filter miRNA:mRNA  list using XM_* ids from Haliotis_rufescens_gca023055435v1rs.xgHalRufe1.0.p.cds.all.fa

```
# Blast (diamond)
https://github.com/bbuchfink/diamond_docs

```BASH
fast=Haliotis_rufescens_gca023055435v1rs.xgHalRufe1.0.p.cds.all.fa

seqkit grep -f miRNA-mRNA.lines $fast -o miRNA-mRNA.cds

EXPORT=/LUSTRE/apps/bioinformatica/diamond_v2.1.8/
export PATH=$PATH:$EXPORT

EXPORT=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software/
export PATH=$PATH:$EXPORT

ref_db=shell_matrix_proteins.pep

diamond makedb --in $ref_db --db ${ref_db%.pep}

diamond blastx -d ${ref_db%.pep} -q miRNA-mRNA.cds -p 20 -k 1 -e 1e-5 -o ${ref_db%.pep}.diamond.blastx.outfmt6 --outfmt 6

# Second db

# seqkit replace -p "\s.+"  > Combined_Gastropoda_Bivalvia.reformat.pep

# Line 3845 has extra-space ^[:space:]>ID, remove it first

ref_db=Combined_Gastropoda_Bivalvia.pep

diamond makedb --in $ref_db --db ${ref_db%.pep}

diamond blastx -d ${ref_db%.pep} -q miRNA-mRNA.cds -p 20 -k 1 -e 1e-5 -o ${ref_db%.pep}.diamond.blastx.outfmt6 --outfmt 6

# get sequences

awk '{print $2}' shell_matrix_proteins.diamond.blastx.outfmt6 | sort | uniq > shell_matrix_proteins.diamond.blastx.outfmt6.subject

awk '{print $2}' Combined_Gastropoda_Bivalvia.diamond.blastx.outfmt6 | sort | uniq > Combined_Gastropoda_Bivalvia.diamond.blastx.outfmt6.subject

seqkit grep -f shell_matrix_proteins.diamond.blastx.outfmt6.subject shell_matrix_proteins.pep

seqkit grep -f Combined_Gastropoda_Bivalvia.diamond.blastx.outfmt6.subject Combined_Gastropoda_Bivalvia.pep


exit

```

# STRING DB
Due to low links per best hits (usually mollusk). I will use only [top organisms](https://string-db.org/cgi/about?footer_active_subpage=statistics) such as human (Model, 9606) and c. elegans (Model of larval development, 6239).


```BASH
# wget https://stringdb-downloads.org/download/protein.sequences.v12.0.fa.gz

wget https://stringdb-downloads.org/download/protein.sequences.v12.0/6239.protein.sequences.v12.0.fa.gz .

wget https://stringdb-downloads.org/download/protein.sequences.v12.0/9606.protein.sequences.v12.0.fa.gz

gunzip protein.sequences.v12.0.fa.gz

ref_db=9606.protein.sequences.v12.0.fa

diamond makedb --in $ref_db --db ${ref_db%.fa}

srun diamond blastx -d ${ref_db%.fa} -q miRNA-mRNA.cds -p 20 -k 1 -e 1e-5 -o ${ref_db%.fa}.diamond.blastx.outfmt6 --outfmt 6 &

```

# Eggnogmapper

https://github.com/eggnogdb/eggnog-mapper/wiki/

```BASH
#!/bin/sh
## Directivas
#SBATCH --job-name=emapper
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=20

EXPORT=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software/eggnog-mapper-master
export PATH=$PATH:$EXPORT

module load conda-2024
source activate base
conda activate trinotate2
# conda deactivate
# conda activate trinotate_test

emapper.py -i miRNA-mRNA.cds --itype CDS --translate --cpu 20 --data_dir /LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Trinotate/TRINOTATE_DB//EGGNOG_DATA_DIR -o eggnog_mapper --override

exit


```


