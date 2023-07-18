# MIR Functional regulatory prediction
Ricardo Gomez-Reyes
# 1) Running tools
## RNAhybrid
```bash
# Usage: RNAhybrid [options] [target sequence] [query sequence].

target=three_prime_utr.fa
query=mature_star_mir.fa
output=${query%.*}_vs_${target%.*}

RNAhybrid -t $target -m 20000 -q $query -n 50 -f 2,8 -s 3utr_human > ${output}_RNAhybrid.out &

```
In order to output a flat text and tabular file derived from RNAhybrid results(*RNAhybrid.out.tsv):

```bash

grep -E "target:|miRNA :|mfe|p-value|length|position" mir_vs_utr_rmdup_RNAhybrid.out > mir_vs_utr_rmdup_RNAhybrid.out.Lines

FILE=mir_vs_utr_rmdup_RNAhybrid.out.Lines

grep -E "target:" $FILE | sed 's/target: //g' > TARGET.tmp
grep -E "miRNA :" $FILE | sed 's/miRNA : //g' > QUERY.tmp
grep -E "mfe" $FILE | sed 's/mfe: //g' | sed 's/kcal\/mol//g' > MFE.tmp
grep -E "p-value" $FILE | sed 's/p-value: //g' > PVAL.tmp
grep -E "position" $FILE | sed 's/position //g' > POS.tmp

grep -A1 "target:" $FILE | grep "length:" | sed 's/length: //g' > TLEN.tmp

grep -A1 "miRNA :" $FILE | grep "length:" | sed 's/length: //g' > QLEN.tmp
 
paste TLEN.tmp QLEN.tmp > LEN.tmp

paste TARGET.tmp QUERY.tmp MFE.tmp PVAL.tmp POS.tmp LEN.tmp | tr -s '[:blank:]' > mir_vs_utr_rmdup_RNAhybrid.out.tsv

rm *.tmp

```

## TargetScan
Using perl script (`targetscan_70.pl`) to identify conserved and non conserved targets sites using a custom set of data. 
```bash

# ./targetscan_70.pl miRNA_file UTR_file PredictedTargetsOutputFile

target=utr_rmdup_ts.txt
query=mature_star_mir.txt

prefix=${target%.*}

output=${query%.*}_vs_$prefix

./targetscan_70.pl $query $target ${output}_targetscan.out &> targetscan.log &



```

# Install tools:
## RNAhybrid
```bash
wget https://bibiserv.cebitec.uni-bielefeld.de/applications/rnahybrid/resources/downloads/RNAhybrid-2.1.2.tar.gz

gunzip RNAhybrid-2.1.2.tar.gz
tar xvf RNAhybrid-2.1.2.tar
cd RNAhybrid-2.1.2
./configure
make
make install

chmod +x RNAhybrid

export PATH=$PATH:$WD

```
## TargetScan
```bash
wget
https://www.targetscan.org/vert_80/vert_80_data_download/targetscan_70.zip

unzip targetscan_70.zip
```
Code from below is recommended by the nf-core framework for community-curated bioinformatics pipelines (Ewels, et al., 2020, doi: 10.1038/s41587-020-0439-x)

### Prepare TargetScan data input
```bash
# After rnahybrid, > 90K utr (~ 35K unique) sequences were significant targeted by mirs:

file=mir_vs_utr_rmdup_RNAhybrid.out.psig.ids 

seqkit replace -p "\s.+" utr_rmdup.fa | seqkit grep -f $file > ${file%.ids}.fa

target=mir_vs_utr_rmdup_RNAhybrid.out.psig.fa 

query=mature_star_mir.fa

prefix=${target%.*}

# 1) Generate aligned UTR or gene (with gaps from alignment)
# ex:
# BMP8B	9606	GUCCACCCGCCCGGC
# BMP8B	9615	-GUG--CUGCCCACC

# Muscle for low number of seqs
# Mafft for high number of seqs 

mafft --thread 12 $target > ${prefix}.aln &

# format for targetscan

align=${prefix}.aln

awk -f linearizefasta.awk < $align > align.tmp

cat align.tmp | cut -f2 > seq.tmp

cat $align | grep ">" | sed 's/>//g' | sed 's/ /;/' > id.tmp

paste id.tmp seq.tmp | awk -v OFS="\t" '{print $1, "0000", $2}' > ${align%.*}_ts.txt


# 2) - miRNA_file    => miRNA families by species
# ex:
# let-7/98	GAGGUAG	9606;10090;10116
# miR-127/127-3p	GAGGUAG	9606;10090

## Convert to TargetScan
# https://github.com/nf-core/circrna/blob/master/bin/targetscan_format.sh

## Subseq window sequences between 2 to 7 (i.e. seed sequence)

grep -v ">" $query > mature_sequence
## Extract seed sequence (7bp after 1st)
awk '{print substr($1, 2, 7)}' mature_sequence > seed_sequence
## Isolate ID (awk last field (NF))
grep ">" $query | awk -F ' ' '{print $NF}' | sed 's/>//g' > miR_ID
## Combine
paste miR_ID seed_sequence > targetscan_tmp.txt
## Correct delimiter, add dummy species
awk -v OFS="\t" '{print $1, $2, "0000"}' targetscan_tmp.txt > mature_star_mir.txt

```

## BASH scripts
```bash
vi linearizefasta.awk

# Paste:

/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;}
     {printf("%s",$0);}
END  {printf("\n");}
```

## Subseq UTR
