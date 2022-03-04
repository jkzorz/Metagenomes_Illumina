# Incorporating Nanopore long reads into Illumina metagenomes

Using long Nanopore reads to try to improve assembly of Illumina metagenomes. 

Nanopore does it's own quality control for the reads (failed vs passed folders). Ignore the failed reads and concatenate all passed reads into one large file. Use porechop to remove adapters.

```
#start interactive session to run commands on arc server 
salloc --mem=50G -c 20 -N 1 -n 1  -t 05:00:00

#concatenate all fastq files that passed QC
cat fastq_pass/* > seqs.fastq

#activate porechop environment
conda activate porechop

#run porechop to trim adapters
porechop -i seqs.fastq -o seqs_trimmed.fastq -t 8

```


## Assembly using nanopore and Illumina reads

Unicycler? (https://github.com/rrwick/Unicycler#installation) - only for single genomes (not metagenomes)

Try using metaspades because it has a "nanopore" option. But they don't guarantee good results: https://github.com/ablab/spades#meta

```
#!/bin/bash
###### Reserve computing resources ######
#SBATCH --mail-user=jacqueline.zorz@ucalgary.ca
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=180GB
#SBATCH --time=24:00:00
#SBATCH --partition=bigmem,cpu2019,cpu2021


###### Set environment variables ######
echo "Starting run at : 'date'"
source /home/jacqueline.zorz/software/miniconda3/etc/profile.d/conda.sh 
conda activate spades

cd /work/ebg_lab/gm/gapp/jzorz/2015-nanopore/

#command

metaspades.py -1 /work/ebg_lab/gm/gapp/jzorz/2015-nanopore/JZ-Condor-BG15-Background-23-24-28_Li32330_S106_R1_QC.fastq -2 /work/ebg_lab/gm/gapp/jzorz/2015-nanopore/JZ-Condor-BG15-Background-23-24-28_Li32330_S106_R2_QC.fastq --nanopore /work/ebg_lab/gm/gapp/jzorz/2015-nanopore/2015_metagenome_test_trimmed.fastq -m 180 -t 40 -o metaspades_hybrid_assembly_BG15_2428


```




## Long read assembly

For assembly of long reads independently: 

Start by trimming adapters: 
```
#trim adapters
cat *.fastq.gz > 2B3_D53_2428_seqs.fastq.gz
conda activate porechop
porechop -i 2B3_D53_2428_seqs.fastq.gz -o 2B3_D53_2428_seqs_trimmed.fastq.gz -t 20

```
Perform assembly of long reads with flye: 

```
flye --nano-raw 2B3_D53_2428_seqs_trimmed.fastq.gz --meta --genome-size 4m --out-dir flye_assembly -i 0 -t 20
```



