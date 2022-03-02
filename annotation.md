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



## CANT-HYD HMMs


Test for loop to run CANT-HYD HMMs on Purple Haze 0-4 cm bins. 
```
for i in bins/bin_PurpleHaze_04.*/genes.faa; do  hmmsearch --tblout hmmsearch_$(basename $(dirname $i)).tblout ../CANT-HYD.hmm $i > $(basename $(dirname $i)).out;done

```
