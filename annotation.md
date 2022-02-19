# Annotation




## CANT-HYD HMMs


Test for loop to run CANT-HYD HMMs on Purple Haze 0-4 cm bins. 
```
for i in bins/bin_PurpleHaze_04.*/genes.faa; do  hmmsearch --tblout hmmsearch_$(basename $(dirname $i)).tblout ../CANT-HYD.hmm $i > $(basename $(dirname $i)).out;done

```
