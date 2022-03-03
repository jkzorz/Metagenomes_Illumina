# Annotation

## DRAM annotation 

Copy all "good" (>75% complete, <5% contamination) bins into a new directory. 

```
cat good_bins.list | while read line; do cp $line good_bins/; done
```

Run DRAM script on all good bins: 

```
#!/bin/bash
###### Reserve computing resources ######
#SBATCH --mail-user=jacqueline.zorz@ucalgary.ca
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=500GB
#SBATCH --time=24:00:00
#SBATCH --partition=bigmem

#set environmental parameters


cd /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/annotation
source /home/jacqueline.zorz/software/miniconda3/etc/profile.d/conda.sh 
conda activate DRAM


####command

DRAM.py annotate -i '/work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/binning/drep_out/dereplicated_genomes/good_bins/*.fa' -o /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/annotation/dram/good_bin_annotations

```

DRAM annotation of all bins won't finish in less than 24 hours (max time allowed with bigmem server cluster). Thus will need to use workaround to get through all bins in "good" bins directory: 

From: https://github.com/WrightonLabCSU/DRAM/issues/54
"Easiest way is to figure out which MAGs that DRAM completed for, move the output somewhere else, move the completed inputs somewhere else, then re-run DRAM on the MAGs that DRAM had not started/completed. The final DRAM outputs are just a cat of the tRNA, rRNA and annotation.tsv files from each individual MAG"

Can then use concatenated tRNA, rRNA, and annotation.tsv files from each MAG as input for DRAM.py distill 

## Barrnap 
Use barrnap to grab rRNA genes from bins 

```
for i in *.fa; do barrnap $i --outseq barrnap_16S/'rrna_'$i; done
```




## CANT-HYD HMMs


Test for loop to run CANT-HYD HMMs on Purple Haze 0-4 cm bins. 
```
for i in bins/bin_PurpleHaze_04.*/genes.faa; do  hmmsearch --tblout hmmsearch_$(basename $(dirname $i)).tblout ../CANT-HYD.hmm $i > $(basename $(dirname $i)).out;done

```
