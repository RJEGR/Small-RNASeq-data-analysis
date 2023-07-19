
# microRNA functional annotation
The mayority of miRs are conserved in both sequence and function across Metazoa. Best know function for miRs is the repression of protein-coding gene expression which is involved in 

## DATABASES
http://www.mirtar2go.org/downloadPage.html
https://mirtarbase.cuhk.edu.cn/~miRTarBase/miRTarBase_2022/php/index.php

## RNAhybrid

https://bibiserv.cebitec.uni-bielefeld.de/rnahybrid?id=rnahybrid_manual_manual

```bash
wget https://bibiserv.cebitec.uni-bielefeld.de/applications/rnahybrid/resources/downloads/RNAhybrid-2.1.2.tar.gz

gunzip RNAhybrid-2.1.2.tar.gz

tar xvf RNAhybrid-2.1.2.tar

cd RNAhybrid-2.1.2

./configure
make
make install

chmod +x /home/rvazquez/RNAhybrid-2.1.2/src/RNAhybrid
export PATH=$PATH:"/home/rvazquez/RNAhybrid-2.1.2/src/"

# Usage: RNAhybrid [options] [target sequence] [query sequence].


```

The use

```bash
export PATH=$PATH:"/home/rvazquez/RNAhybrid-2.1.2/src/"

RNAhybrid -t genome_or_transcriptome.fa -q miRNAs.fasta

 -m <max targetlength>
  -n <max query length>

```

## miRanda
https://cbio.mskcc.org/miRNA2003/miranda.html (old version)
upgrade: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2945792/
```bash
wget https://cbio.mskcc.org/miRNA2003/src1.9/binaries/miRanda-1.9-i686-linux-gnu.tar.gz

gunzip miRanda-1.9-i686-linux-gnu.tar.gz; tar xvf miRanda-1.9-i686-linux-gnu.tar

chmod +x /home/rvazquez/miRanda-1.9-i686-linux-gnu/bin/miranda

# No such file direct error
# https://www.baeldung.com/linux/no-such-file-or-directory-error

# evaluar que se tengan las dependencias de miranda (TRUE)

objdump -p miranda | grep NEEDED

# NEEDED               libm.so.6
# NEEDED               libc.so.6

ls /lib/x86_64-linux-gnu/lib*.so.6
# al parecer el binario is "not a dynamic executable linux"
# it is just a problem w/ a 32-bit binary system. then:

sudo apt-get install gcc-multilib

uname -a

# x86_64 GNU/Linux
# Then

export PATH=$PATH:"/home/rvazquez/miRanda-1.9-i686-linux-gnu/bin"

miranda --h


```

## SubmiRine
```bash
# https://research.nhgri.nih.gov/software/SubmiRine/download.shtml

wget https://research.nhgri.nih.gov/software/SubmiRine/downloads/SubmiRine.tar.gz
gunzip SubmiRine.tar.gz
tar xvf SubmiRine.tar

export PATH=$PATH:"/home/rvazquez/SubmiRine/src"
# aux  MiRNA.py  SubmiRine_Compare.R  SubmiRine_Search.py  Target.py  utils.py  UTR.py

#  Python version 2.7
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda install -c anaconda numpy

SubmiRine_Search.py 
# CHECK PY VERSION
```

## miRNAconsTargets
miRNAconsTarget: consensus target prediction on user provided input data based on Miranda , PITA and TargetSpy for animal microRNAs.

```bash
git clone https://github.com/bioinfoUGR/sRNAtoolbox.git

export PATH=$PATH:"/home/rvazquez/SRNA_TOOLBOX/sRNAtoolbox/exec"

miRNAconsTargets.jar --help

target=utr.fa
query=mir.fasta

miRNAconsTargets.jar $query $target $PWD 12 TS:MIRANDA &> miRNAconsTargets.log &

# Program string: each program has a tag: TS (targetSpy); PITA (pita); MIRANDA (miranda); simple seed method (SEED), For example, TS:MIRANDA will lauch the target prediction using TargetSpy and MIRANDA
# TEST

# miRNAconsTargets   microRNA_file   utr_file    output_directory(ABSOLUTE PATH MUST BE GIVEN)    p(threads)    program_string    program_parameters


# https://arn.ugr.es/srnatoolbox/static/sRNAtoolbox_manual.pdf

```
## DMISO (ERROR complejo de compilar tensorflow)
https://www.nature.com/articles/s41598-022-14890-8
https://github.com/amlantalukder/DMISO

```bash
#  (1) Install Python 3
#   (2) Install numpy, itertools, sklearn packages for Python 3.

wget http://hulab.ucf.edu/research/projects/DMISO/DMISO.zip
mv tool DMISO 
cd DMISO;chmod +x dmiso.py

# add header line in dmiso.py as follow
# #!/usr/bin/env python3

# INSTALL PIP FOR THIS VERSION 
# https://pip.pypa.io/en/stable/installation/


export PATH=$PATH:"/home/rvazquez/.local/bin/"

pip3 install numpy more-itertools scikit-learn

# Keras and TensorFlow are open source Python libraries for working with neural networks, creating machine learning models and performing deep learning. Because Keras is a high level API for TensorFlow, they are installed together.

pip3 install keras tensorflow-cpu

pip3 install -U tensorflow # tensorflow is not available for py 3.11

pip3 freeze | grep -E "keras|tensorflow|numpy|itertools|learn"

# Test
./dmiso.py -p examples/test_pairs.txt -o examples/test_output.txt
# OR
./dmiso.py -m examples/test_miRNAs.fa -t examples/test_mRNAs.fa

```
Tensorflow installation:
`2023-03-17 14:58:14.090618: I tensorflow/core/platform/cpu_feature_guard.cc:193] This TensorFlow binary is optimized with oneAPI Deep Neural Network Library (oneDNN) to use the following CPU instructions in performance-critical operations:  FMA
To enable them in other operations, rebuild TensorFlow with the appropriate compiler flags`

```bash
 # bazel
 
 sudo apt install apt-transport-https curl gnupg -y

curl -fsSL https://bazel.build/bazel-release.pub.gpg | gpg --dearmor >bazel-archive-keyring.gpg

sudo mv bazel-archive-keyring.gpg /usr/share/keyrings

echo "deb [arch=amd64 signed-by=/usr/share/keyrings/bazel-archive-keyring.gpg] https://storage.googleapis.com/bazel-apt stable jdk1.8" | sudo tee /etc/apt/sources.list.d/bazel.list

sudo apt update && sudo apt install bazel-5.3.0

sudo apt update && sudo apt full-upgrade
```
Then configure and build tensorflow
```bash

# Clone the TensorFlow repository
 git clone https://github.com/tensorflow/tensorflow.git
 
 # Navigate to the TensorFlow repository
 cd tensorflow
 
 git checkout r2.11
 
 # Configure the build
./configure
 
 # Use bazel build to create the TensorFlow 2.x package-builder with CPU-only support:
 bazel build [--config=option] //tensorflow/tools/pip_package:build_pip_package
 
 bazel build --config=opt //tensorflow/tools/pip_package:build_pip_package
 
bazel build [--config=option] //tensorflow/tools/pip_package:build_pip_package
```

## miRAW (hard to train the model )
```bash
git clone https://app86@bitbucket.org/bipous/miraw_dl4mirna_bin
# git clone git@bitbucket.org:bipous/miraw_dl4mirna_binaries.git

chmod +x miraw_dl4mirna_binaries/miRAW.jar

export PATH=$PATH:"/home/rvazquez/MIRS_FUNCTIONAL_ANNOT/miRAW/miraw_dl4mirna_binaries"



miRAW.jar --help

Options:
	GenePrediction: predict miRNA-gene interactions
	ClassifySites:
	TrainModel:
	CreateDataSets:
	CrossValidation:
	
miRAW.jar GenePrediction help

miRAW GeneClassification evaluates a set of 3'UTR gene and microRNAs transcripts to predict their interactions and target sites.

Usage: miRaw GeneClassification option configurationFile
	option: predict - predicts the target sites and saves them to a file.
	option: evaluate - predicts the gene-miRNA interactions and validates them against a Unified file.

	Configuration must have the following parametes
		UnifiedFile: location csv file containing the geneName, geneId, miRRNAname, 3UTR transcritpt, mature miRNA transcript, validation (only for evaluation)
		DLModel: location of the Deep Learning Model to use
		ExperimentFolder: folder where the experiment results will be saved
		ExperimentName: name of the experiment
		
#

git clone https://bitbucket.org/bipous/miraw_-dl4mirnatargeting_code.git

miraw_-dl4mirnatargeting_code/TestData/ConfigFile/genePrediction.properties

########################################
###########GenePrediction###############
########################################
ExperimentName=FullCycle
ExperimentFolder=/Users/apla/Documents/MirnaTargetDatasets/BioinformaticsPaper/GenePredResults
UnifiedFile=/Users/apla/Documents/MirnaTargetDatasets/BioinformaticsPaper/DataSet/Original/tarBamirseValidTranscript.csv
DLModel=/Users/apla/Documents/MirnaTargetDatasets/BioinformaticsPaper/CVResults/BioinformaticsCV/Set1/lastModel.bin

```

## targetscan 7
https://elifesciences.org/articles/05005#abstract


## TargetNet (well mantained)
https://github.com/mswzeus/TargetNet

Datasets: The complete TargetNet algorithm was evaluated with two types of experimentally verified miRNA–mRNA pair datasets, (i) miRAW and (ii) log-fold change (LFC) test datasets.

**miRNA–mRNA pair datasets:**
+ **Train validated:** First, we used miRAW test datasets with binary labels indicating functional and non-functional targets. They originated from DIANA-TarBase (Vlachos et al., 2015) and MirTarBase (Chou et al., 2016) databases. The miRAW test datasets can help us evaluate the functional miRNA target classification performance of TargetNet.

+ **Train validated**: Second, we used eleven microarray and two RNA-seq LFC test datasets with real-valued labels indicating the level of functionality of miRNA targets. In each microarray and RNA-seq dataset, a miRNA was individually transfected into HeLa and HEK329 cells, respectively. Then, the log-fold change of mRNA expression was measured. More negative labels indicate more functional miRNA–mRNA pairs, which strongly down-regulate the targeted genes. The LFC test datasets can help us to evaluate how well TargetNet distinguishes high-functional miRNA targets.

**miRNA–CTS pair datasets**
**- Predicted model:** The prediction model of TargetNet was trained with miRAW miRNA–CTS pair datasets. 


4.1.1 miRNA–mRNA pair datasets

```bash
git clone https://github.com/mswzeus/TargetNet.git

cd TargetNet

conda env create -f TargetNet.yaml  # take more than 1 hour

conda create -n TargetNet python=3.8

source activate TargetNet

conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --add channels pytorch
# conda list

# -c CHANNEL

conda install pytorch==1.10.2 torchvision==0.11.3 cudatoolkit=11.3 -c pytorch -c conda-forge

pip install -r requirements.txt

# 

export PATH=$PATH:"/home/rvazquez/TargetNet"

python3.8 train_model.py --help

CUDA_VISIBLE_DEVICES=0

# training the model
python3.8 train_model.py --data-config config/data/miRAW_train.json --model-config config/model/TargetNet.json --run-config config/run/run.json --output-path results/TargetNet_training/

#
CUDA_VISIBLE_DEVICES=0
python3.8 evaluate_model.py --data-config config/data/miRAW_eval.json --model-config config/model/TargetNet.json --run-config config/run/run.json --checkpoint pretrained_models/TargetNet.pt --output-path results/TargetNet-evaluation/

# For using other datasets, modify the data paths specified in the miRAW_eval.json data-config file.

# config/data/miRAW_train.json

{
  "train_path": "data/miRAW_Train_Validation.txt",
  "val_path": "data/miRAW_Train_Validation.txt"
}

# mirna_id        mirna_seq       mrna_id mrna_seq        label   split

# config/data/miRAW_eval.json
# Your imput samples to deep search test

{
"test0_path": "data/miRAW_Test0.txt",
"test1_path": "data/miRAW_Test1.txt",
"test2_path": "data/miRAW_Test2.txt",
"test3_path": "data/miRAW_Test3.txt",
"test4_path": "data/miRAW_Test4.txt",
"test5_path": "data/miRAW_Test5.txt",
"test6_path": "data/miRAW_Test6.txt",
"test7_path": "data/miRAW_Test7.txt",
"test8_path": "data/miRAW_Test8.txt",
"test9_path": "data/miRAW_Test9.txt"
}

# RESULTS
results/TargetNet-evaluation/
```




## Others
https://bitbucket.org/account/user/bipous/projects/MIRAW.
https://bitbucket.org/leslielab/chimiric/src/master/

###  RUNNING TOOLS

```bash
# vi linearizefasta.awk
/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;}
     {printf("%s",$0);}
END  {printf("\n");}

# 

awk -f linearizefasta.awk < three_prime_utr_rmdup.aln
```

```bash
export PATH=$PATH:"/home/rvazquez/seqkit_tool"

export PATH=$PATH:"/home/rvazquez/miRanda-1.9-i686-linux-gnu/bin"

export PATH=$PATH:"/home/rvazquez/RNAhybrid-2.1.2/src/"

ln -s /home/rvazquez/SHORTSTACKS/knownsRNA.fa
ln -s /home/rvazquez/SHORTSTACKS/ShortStack_20230315_out/mir.fasta/mir.fasta

# ln -s /home/rvazquez/GENOME_20230217/ENSEMBLE/utr.fa .
ln -s /home/rvazquez/GENOME_20230217/SHORTSTACKS_REF/utr_rmdup.fa .

ln -s /home/rvazquez/GENOME_20230217/SHORTSTACKS_REF/three_prime_utr_rmdup.fa .

# Prepare mir.fasta
# mir.fasta: This is a FASTA formatted file containing hairpin, mature miRNA, and miRNA* sequences derived from ShortStack's identification of MIRNA loci.

seqkit grep -r -p "mature|star" mir.fasta > mature_star_mir.fa

target=three_prime_utr.fa

query=mature_star_mir.fa

output=${query%.*}_vs_${target%.*}

# miranda $query $target > miranda.out &
miranda $query $target -quiet -out ${output}_miranda.out

# free(): invalid next size (fast)
# [1]+  Aborted                 (core dumped) miranda $query $target > miranda.out

# Segmentation fault (core dumped)
# sudo apt-get clean all

# 2) (RNAHYB Tiene sezgos debido a sus calibraciones en el flag -s, ademas no sabemos los rangos del flag -d que omite a -s)

# RNAhybrid [options] [target sequence] [query sequence].

RNAhybrid -t $target -m 20000 -q $query -n 50 -f 2,8 -s 3utr_human > ${output}_RNAhybrid.out &

# Create table based on output:
# grep -E "target:|miRNA :|mfe|p-value|length|position"

grep -E "target:|miRNA :|mfe|p-value|length|position" mir_vs_utr_rmdup_RNAhybrid.out > mir_vs_utr_rmdup_RNAhybrid.out.Lines

FILE=mir_vs_utr_rmdup_RNAhybrid.out.Lines

grep -E "target:" $FILE | sed 's/target: //g' > TARGET.tmp
grep -E "miRNA :" $FILE | sed 's/miRNA : //g' > QUERY.tmp
grep -E "mfe" $FILE |sed 's/mfe: //g'  | sed 's/kcal\/mol//g' > MFE.tmp
grep -E "p-value" $FILE | sed 's/p-value: //g' > PVAL.tmp
# grep -E "length" $FILE | sed 's/length: //g' > LEN.tmp
grep -E "position" $FILE | sed 's/position //g' > POS.tmp

# how to split LEN-2-rows in two cols (LEN.tmp)
#

grep -A1 "target:" $FILE | grep "length:" | sed 's/length: //g' > TLEN.tmp

grep -A1 "miRNA :" $FILE | grep "length:" | sed 's/length: //g' > QLEN.tmp
 
paste TLEN.tmp QLEN.tmp > LEN.tmp

paste TARGET.tmp QUERY.tmp MFE.tmp PVAL.tmp POS.tmp LEN.tmp | tr -s '[:blank:]' > mir_vs_utr_rmdup_RNAhybrid.out.tsv

rm *.tmp

scp rvazquez@200.23.162.234:/home/rvazquez/MIRS_FUNCTIONAL_ANNOT/mir_vs_utr_rmdup_RNAhybrid.out.tsv .

# evalute if significance:

# grep "p-value:" mir_vs_utr_rmdup_RNAhybrid.out  | cut -f2 -d ":" | sort -n | head

  
# -s (3utr_fly|3utr_worm|3utr_human)

# The helix constraint format (-f) is "from,to", eg. -f 2,7 forces structures to have a helix from position 2 to 7 with respect to the query.

# https://manpages.ubuntu.com/manpages/trusty/man1/RNAhybrid.1.html

# (-d <xi>,<theta>) <xi> and <theta> are the position and shape parameters, respectively, of the extreme value distribution assumed for p-value calculation. If omitted, they are estimated from the maximal duplex energy of the query. In that case, a data set name has to be given with the -s flag.


```
# TEST TARGET SCAN
```bash

# TARGETSCAN
# After rnahybrid, > 90K utr sequences were significant targeted by mirs:



file=mir_vs_utr_rmdup_RNAhybrid.out.psig.ids # 37 663 unique 

seqkit replace -p "\s.+" utr_rmdup.fa | seqkit grep -f $file > ${file%.ids}.fa

grep -c "^>" ${file%.ids}.fa

target=mir_vs_utr_rmdup_RNAhybrid.out.psig.fa # utr_rmdup.fa
query=mature_star_mir.fa

prefix=${target%.*}
output=${query%.*}_vs_$prefix

# 1) align UTR
# contains three fields (tab-delimited):

# a. Gene/UTR ID or name
# b. Species ID for this gene/UTR (must match ID in miRNA file)
# c. Aligned UTR or gene (with gaps from alignment)

# Generate aligned UTR or gene (with gaps from alignment)
# ex:
# BMP8B	9606	GUCCACCCGCCCGGC
# BMP8B	9615	-GUG--CUGCCCACC

# Muscle for low number of seqs
# Mafft for high n 

# cat $target | seqkit sample -n 1000 -o ${prefix}.subset.fa
# --globalpair # flag hard to run
# fast way to align:

mafft --thread 20 $target > ${prefix}.aln &

# ps -u rvazquez
# # to kill the process related to ps -o pid=pidid | xargs kill

##format for targetscan

align=${prefix}.aln

awk -f linearizefasta.awk < $align > align.tmp
cat align.tmp | cut -f2 > seq.tmp

cat $align | grep ">" | sed 's/>//g' | sed 's/ /;/' > id.tmp

paste id.tmp seq.tmp | awk -v OFS="\t" '{print $1, "0000", $2}' > ${align%.*}_ts.txt


# 2) - miRNA_file    => miRNA families by species
# contains three fields (tab-delimited):
# a. miRNA family ID/name
# b. seed region (7mer) for this miRNA
# c. semicolon-delimited list of species IDs in which this miRNA has been annotated
# ex:
# let-7/98	GAGGUAG	9606;10090;10116
# miR-127/127-3p	GAGGUAG	9606;10090

## Convert to TargetScan
# https://github.com/nf-core/circrna/blob/master/bin/targetscan_format.sh

## Isolate the sequences
query=mature_star_mir.fa

grep -v ">" $query > mature_sequence
## Extract seed sequence (7bp after 1st)
awk '{print substr($1, 2, 7)}' mature_sequence > seed_sequence
## Isolate ID (awk last field (NF))
grep ">" $query | awk -F ' ' '{print $NF}' | sed 's/>//g' > miR_ID
## Combine
paste miR_ID seed_sequence > targetscan_tmp.txt
## Correct delimiter, add dummy species
awk -v OFS="\t" '{print $1, $2, "0000"}' targetscan_tmp.txt > mature.txt

# ./targetscan_70.pl miRNA_file UTR_file PredictedTargetsOutputFile
./targetscan_70.pl mature.txt ${align%.*}_ts.txt ${output}_targetscan.out &> targetscan.log &

# Processing JALGQA010000008.1_47463137-47463197:+;LOC125376846
#  ....

scp rvazquez@200.23.162.234:/home/rvazquez/MIRS_FUNCTIONAL_ANNOT/TARGETSCAN/mature_star_mir_vs_mir_vs_utr_rmdup_RNAhybrid.out.psig_targetscan.out .
# /home/rvazquez/MIRS_FUNCTIONAL_ANNOT/TARGETSCAN

# RUNNING /home/rvazquez/MIRS_FUNCTIONAL_ANNOT/BLAST_FROM_PREDICTED_TARGETS

# creating a diamond-formatted database file
# https://github.com/bbuchfink/diamond/wiki

mkdir DATABASE
cd DATABASE
ln -s /home/rvazquez/Trinotate-Trinotate-v4.0.0/DATABASE/uniprot_sprot.pep .

cd ..

query=mir_vs_utr_rmdup_RNAhybrid.out.psig.gene_features_known_annot.fa

../diamond makedb --in uniprot_sprot.pep -d uniprot_sprot

# running a search in blastx mode (less than 1 min)

./diamond blastx -d DATABASE/uniprot_sprot -q $query -o uniprot_sprot.ncbi.blastx.tsv

scp rvazquez@200.23.162.234:/home/rvazquez/MIRS_FUNCTIONAL_ANNOT/BLAST_FROM_PREDICTED_TARGETS/uniprot_sprot.ncbi.blastx.tsv .

# UNIPROT TO GENE ONTOLOGY
# GET GENE ONTOLOGY OBO FORMAT
http://current.geneontology.org/ontology/

https://github.com/Trinotate/Trinotate/blob/master/util/admin/Build_Trinotate_Boilerplate_SQLite_db.pl

# blastx -query $query -db /home/rvazquez/Trinotate-Trinotate-v4.0.0/DATABASE/uniprot_sprot.pep -num_threads 1 -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -out uniprot_sprot.ncbi.blastx.outfmt6 2> blastx.log &