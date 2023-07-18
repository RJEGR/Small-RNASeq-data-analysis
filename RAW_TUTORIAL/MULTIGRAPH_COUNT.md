# MGCount 

```bash
pip3 install git+https://github.com/hitaandrea/MGcount.git
python3 -m mgcount 
```

Esta herramienta de quiere de formato GTF muy detallado que se genera con diferentes fuentes... algo dificil por ahora, pero la herramienta es util para visualizacion de datos finales.

MGcount starts with a set of genomic alignments of RNA-seq reads (one BAM file per sample/cell) and a set of RNA feature annotations stored in a single gene transfer format (GTF) file. 

The scope of the MGcount quantification is bounded by the features annotated in the reference GTF file. To maximize the scope of the analysis, we combined annotations from DASHR, RNAcentral, miRbase and Ensembl in a single GTF file. The MGcount repository provides integrated GTF annotations for human, arabidopsis, mouse and nematode, and the corresponding R scripts used for their generation. These can be used as a template to integrate annotations from other species: https://github.com/hitaandrea/MGcount/tree/main/R



```bash
# featureCounts not found. Please, check the program is available on the path specified by --featureCounts_path argument (default: /user/bin/featureCounts)

wget -o subread-2.0.4-source.tar.gz  https://sourceforge.net/projects/subread/files/subread-2.0.4/subread-2.0.4-source.tar.gz/download


make -f Makefile.MacOS

export PATH=$PATH:"/Users/cigom/Documents/MIRNA_HALIOTIS/MGCounts/subread-2.0.4-source/bin"

# test (WORK)
python3 -m mgcount -T 2 --gtf annotations_gtf/Homo_sapiens.GRCh38.custom.gtf --outdir outputs --bam_infiles input_bamfilenames.txt --featureCounts_path /Users/cigom/Documents/MIRNA_HALIOTIS/MGCounts/subread-2.0.4-source/bin/featureCounts

# OWN DATA

printf "$PWD/%s\n" *.bam > input_bamfilenames.txt


python3 -m mgcount -T 2 --gtf Haliotis_rufescens_gca023055435v1rs.xgHalRufe1.0.p.56.gtf --outdir outputs --bam_infiles input_bamfilenames.txt --featureCounts_path /Users/cigom/Documents/MIRNA_HALIOTIS/MGCounts/subread-2.0.4-source/bin/featureCounts

```