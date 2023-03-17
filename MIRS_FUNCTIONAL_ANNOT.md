# microRNA functional annotation
The mayority of miRs are conserved in both sequence and function across Metazoa. Best know function for miRs is the repression of protein-coding gene expression which is involved in 

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

export PATH=$PATH:"/home/rvazquez/MIRS_FUNCTIONAL_ANNOT/SubmiRine/src"

# CHECK PY VERSION
```


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