# Comparing to MAGs from NCBI


```
/work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/compare
```


#unzip files
```
for i in 4484-113/*.fna.gz; do gunzip $i; done
```


#change file extensions
```for i in 4484-113/*.fna; do sample=$(basename $i .fna).fa; mv $i 4484-113/$sample; done```
