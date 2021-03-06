# Assembly


Comparative review on metagenome assembly: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0169662

### Megahit 
- Better for microdiversity 
- Faster/less RAM
- Produces shorter contigs
- Might struggle with larger datasets in assembling abundant species
- biased towards relatively low coverage genomes
- "divide and conquer" approach might be best option for large datasets? 


### MetaSpades
- Produces longer contigs
- Captures high diversity even at high complexity and low read coverage
- Takes much more memory and time 


Skip read merging with bbmerge? 

2A-1 ("The Hole") 12-16cm sample: Megahit assembly took ~5-7 hours, Metaspades assembly took 19 hours. 

## Assembly with Megahit

Assemble each sample separately with Megahit 

**Megahit version: v1.2.9**

```
#!/bin/bash
###### Reserve computing resources ######
#SBATCH --mail-user=jacqueline.zorz@ucalgary.ca
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=180GB
#SBATCH --time=168:00:00
#SBATCH --partition=cpu2019,cpu2021


###### Set environment variables ######
echo "Starting run at : 'date'"
source /home/jacqueline.zorz/software/miniconda3/etc/profile.d/conda.sh 
conda activate megahit
cd /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/megahit/


#do the assembly using quality controlled reads
for f in /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/bbduk/cat_qc/*R1_QC.fastq; 

do 	
	megahit -1 $f -2 $(dirname $f)/$(basename $f R1_QC.fastq)R2_QC.fastq -t 40 -o megahit_$(basename $f _R1_QC.fastq) --min-contig-len 500;
done
```

## Megahit with subsampled reads

Using megahit to assemble reads subsampled to 10% of original metagenome: ~10M paired reads. 

```
#!/bin/bash
###### Reserve computing resources ######
#SBATCH --mail-user=jacqueline.zorz@ucalgary.ca
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=150GB
#SBATCH --time=24:00:00
#SBATCH --partition=bigmem,cpu2019,cpu2021


###### Set environment variables ######
echo "Starting run at : 'date'"
source /home/jacqueline.zorz/software/miniconda3/etc/profile.d/conda.sh 
conda activate megahit


#do the assembly using quality controlled reads
megahit -1 /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/bbduk/JZ-Condor-2B1-PurplePatch-A54-0-4_Li32225_S1_R1_QC_subsample10pc.fastq -2 /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/bbduk/JZ-Condor-2B1-PurplePatch-A54-0-4_Li32225_S1_R2_QC_subsample10pc.fastq -t 25 -o megahit_JZ-Condor-2B1-PurplePatch-A54-0-4_subsample10pc --min-contig-len 500
```

## Quality of Assembly with MetaQuast

Example: 

```
metaquast.py --fast -m 500 -t 20 -o metaquast -l megahit final.contigs.fa
```


## Stats.sh bbtools

Use stats.sh to get quick overview of assembly statistics 
```
conda activate bbtools

stats.sh in=/work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/megahit/megahit_JZ-Condor-2B1-PurplePatch-A54-0-4_norm100/final.contigs.fa
```


