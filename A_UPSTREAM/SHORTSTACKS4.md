# small RNA assembly, annotation and quantification
Ricardo Gómez Reyes
# 1) Running tools
## Shortstacks
```bash
MIRGENEDB=/home/rvazquez/MIRGENEDB_20230314/ALL-mat.fa
PIRBASE=/home/rvazquez/MIRTRACE/CUSTOM_DB/PIRBASE_V3/PIRBASE_rmdup.fasta

cat $MIRGENEDB $PIRBASE | seqkit rmdup -n -o knownsRNA.fa

for i in $(ls *newid.subset.fasta); do ../bbmap/reformat.sh $i out=${i%.fasta}.fq;done

readfile=`ls -x *.newid.subset.fq`

ShortStacks4 --genomefile $REF --known_miRNAs knownsRNA.fa --dn_mirna --outdir ShortStack_"$(date +%Y%m%d)"_out --threads 24 --dicermax 30 --mmap u --mincov 0.8 --pad 1 --readfile $readfile &>> "ShortStack_"$(date +%Y%m%d)".log" &
```

# 2) Installation
```bash
conda create --name ShortStack4 shortstack 
conda activate ShortStack4
```