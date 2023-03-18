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

chmod + /home/rvazquez/RNAhybrid-2.1.2/src/RNAhybrid
export PATH=$PATH:"/home/rvazquez/RNAhybrid-2.1.2/src/"

# Usage: RNAhybrid [options] [target sequence] [query sequence].


```

The use

```bash
export PATH=$PATH:"/home/rvazquez/RNAhybrid-2.1.2/src/"

RNAhybrid -t genome_or_transcriptome.fa -q miRNAs.fasta


```

## miRanda
https://cbio.mskcc.org/miRNA2003/miranda.html
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

## TargetNet (well mantained)
https://github.com/mswzeus/TargetNet

```bash
git clone https://github.com/mswzeus/TargetNet.git
cd TargetNet
conda env create -f TargetNet.yaml  # take more than 1 hour

conda create -n TargetNet python=3.8

conda activate TargetNet

conda install pytorch==1.10.2 torchvision==0.11.3 cudatoolkit=11.3 -c pytorch -c conda-forge

pip install -r requirements.txt

```

## Others
https://bitbucket.org/account/user/bipous/projects/MIRAW.
https://bitbucket.org/leslielab/chimiric/src/master/

###  RUNNING TOOLS


```bash

export PATH=$PATH:"/home/rvazquez/miRanda-1.9-i686-linux-gnu/bin"

export PATH=$PATH:"/home/rvazquez/RNAhybrid-2.1.2/src/"

ln -s /home/rvazquez/SHORTSTACKS/knownsRNA.fa
ln -s /home/rvazquez/SHORTSTACKS/ShortStack_20230314_test1/mir.fasta


target=/home/rvazquez/GENOME_20230217/ENSEMBLE/utr.fa
query=mir.fasta



miranda $query $target > miranda.out &

# free(): invalid next size (fast)
# [1]+  Aborted                 (core dumped) miranda $query $target > miranda.out

# Segmentation fault (core dumped)
# sudo apt-get clean all

# 2) (RNAHYB Tiene sezgos debido a sus calibraciones en el flag -s, ademas no sabemos los rangos del flag -d que omite a -s)

# RNAhybrid [options] [target sequence] [query sequence].
RNAhybrid -t $target -m 20000 -q $query -n 50 -f 2,8 -s 3utr_human > RNAhybrid.out &

# seqkit stats $target # 16,632 max len
# seqkit stats $query # 49 max len
# -m <max targetlength>
# -n <max query length>
  
# -s (3utr_fly|3utr_worm|3utr_human)

# The helix constraint format (-f) is "from,to", eg. -f 2,7 forces structures to have a helix from position 2 to 7 with respect to the query.

# https://manpages.ubuntu.com/manpages/trusty/man1/RNAhybrid.1.html

# (-d <xi>,<theta>) <xi> and <theta> are the position and shape parameters, respectively, of the extreme value distribution assumed for p-value calculation. If omitted, they are estimated from the maximal duplex energy of the query. In that case, a data set name has to be given with the -s flag.

```