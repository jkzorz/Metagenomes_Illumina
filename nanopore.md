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

Four rounds of mapping and polishing with bwa and racon 

**Having issues with bwa mem mapping on larger nanopore files... May need to substitute bwa mem for minimap2 here...**


```
#!/bin/bash
###### Reserve computing resources ######
#SBATCH --mail-user=jacqueline.zorz@ucalgary.ca
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=50GB
#SBATCH --time=3:00:00
#SBATCH --partition=cpu2019,apophis-bf,pawson-bf,razi-bf

ASSEMBLY=${1?Error: missing assembly location}
SEQS=${2?Error: missing fastq sequence location}

###### Set environment variables ######
source /home/jacqueline.zorz/software/miniconda3/etc/profile.d/conda.sh 
conda activate racon

###### Set environment variables ######
echo "Starting run at : 'date'"
#step 1
/work/ebg_lab/gm/gapp/jzorz/bwa/bwa index $ASSEMBLY


/work/ebg_lab/gm/gapp/jzorz/bwa/bwa mem -t 10 -x ont2d $ASSEMBLY $SEQS -o assembly.mapping.sam


/work/ebg_lab/gm/gapp/jzorz/racon/build/bin/racon -t 10 $SEQS assembly.mapping.sam $ASSEMBLY > racon1.fasta

sed -n '/^>/,$p' racon1.fasta | sed 's/\s.*$//g' > racon1.mod.fasta

#step 2
/work/ebg_lab/gm/gapp/jzorz/bwa/bwa index racon1.mod.fasta

/work/ebg_lab/gm/gapp/jzorz/bwa/bwa mem -t 10 -x ont2d racon1.mod.fasta $SEQS -o racon1.mapping.sam

/work/ebg_lab/gm/gapp/jzorz/racon/build/bin/racon -t 10 $SEQS racon1.mapping.sam racon1.mod.fasta > racon2.fasta

sed -n '/^>/,$p' racon2.fasta | sed 's/\s.*$//g' > racon2.mod.fasta

#step 3
/work/ebg_lab/gm/gapp/jzorz/bwa/bwa index racon2.mod.fasta

/work/ebg_lab/gm/gapp/jzorz/bwa/bwa mem -t 10 -x ont2d racon2.mod.fasta $SEQS -o racon2.mapping.sam

/work/ebg_lab/gm/gapp/jzorz/racon/build/bin/racon -t 10 $SEQS racon2.mapping.sam racon2.mod.fasta > racon3.fasta

sed -n '/^>/,$p' racon3.fasta | sed 's/\s.*$//g' > racon3.mod.fasta

#step 4 
/work/ebg_lab/gm/gapp/jzorz/bwa/bwa index racon3.mod.fasta

/work/ebg_lab/gm/gapp/jzorz/bwa/bwa mem -t 10 -x ont2d racon3.mod.fasta $SEQS -o racon3.mapping.sam

/work/ebg_lab/gm/gapp/jzorz/racon/build/bin/racon -t 10 $SEQS racon3.mapping.sam racon3.mod.fasta > racon4.fasta

sed -n '/^>/,$p' racon4.fasta | sed 's/\s.*$//g' > racon4.mod.fasta

##
echo "Job finished with exit code $? at: 'date'"
##
```


**Instead of using bwa mem like above, manually ran all polishing steps using minimap2 and racon**

```
conda activate minimap 
#step 1 - mapping
minimap2 -ax map-ont flye_assembly/assembly.fasta Nanopore_2A2_seqs_trimmed.fastq > map1.sam
#step 1 - polish 
/work/ebg_lab/gm/gapp/jzorz/racon/build/bin/racon -t 20 Nanopore_2A2_seqs_trimmed.fastq map1.sam flye_assembly/assembly.fasta > racon1.fasta
sed -n '/^>/,$p' racon1.fasta | sed 's/\s.*$//g' > racon1.mod.fasta

#step 2 - mapping
conda activate minimap
minimap2 -ax map-ont racon1.mod.fasta Nanopore_2A2_seqs_trimmed.fastq > racon1.map.sam
#step 2 - polish 
/work/ebg_lab/gm/gapp/jzorz/racon/build/bin/racon -t 20 Nanopore_2A2_seqs_trimmed.fastq racon1.map.sam racon1.mod.fasta > racon2.fasta
sed -n '/^>/,$p' racon2.fasta | sed 's/\s.*$//g' > racon2.mod.fasta

#step 3 - mapping
minimap2 -ax map-ont racon2.mod.fasta Nanopore_2A2_seqs_trimmed.fastq > racon2.map.sam
#step 3 - polish
/work/ebg_lab/gm/gapp/jzorz/racon/build/bin/racon -t 20 Nanopore_2A2_seqs_trimmed.fastq racon2.map.sam racon2.mod.fasta > racon3.fasta
sed -n '/^>/,$p' racon3.fasta | sed 's/\s.*$//g' > racon3.mod.fasta

#step 4 - mapping 
minimap2 -ax map-ont racon3.mod.fasta Nanopore_2A2_seqs_trimmed.fastq > racon3.map.sam
#step 4 - polish
/work/ebg_lab/gm/gapp/jzorz/racon/build/bin/racon -t 20 Nanopore_2A2_seqs_trimmed.fastq racon3.map.sam racon3.mod.fasta > racon4.fasta
sed -n '/^>/,$p' racon4.fasta | sed 's/\s.*$//g' > racon4.mod.fasta
```




Final round of polishing with medaka: 
Takes a long time with large files - over 1 day. 

```
#!/bin/bash
###### Reserve computing resources ######
#SBATCH --mail-user=jacqueline.zorz@ucalgary.ca
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=180GB
#SBATCH --time=100:00:00
#SBATCH --partition=cpu2021,cpu2019


###### Set environment variables ######
source /home/jacqueline.zorz/software/miniconda3/etc/profile.d/conda.sh 
conda activate medaka


medaka_consensus -i Nanopore_2A2_seqs_trimmed.fastq -d racon4.mod.fasta -t 25 -m r941_min_high_g303 -o Medaka_polish4

```

Map reads to polished contigs with minimap

```
minimap2 -ax map-ont Medaka_polish4/consensus.fasta Nanopore_2A2_seqs_trimmed.fastq > 2A2_polished_seqs.sam


###Output of 2B3_D53 run:

2428_seqs_trimmed.fastq.gz > 2B3_D53_2428.sam
[M::mm_idx_gen::0.124*0.88] collected minimizers
[M::mm_idx_gen::0.146*1.19] sorted minimizers
[M::main::0.146*1.19] loaded/built the index for 190 target sequence(s)
[M::mm_mapopt_update::0.155*1.18] mid_occ = 23
[M::mm_idx_stat] kmer size: 15; skip: 10; is_hpc: 0; #seq: 190
[M::mm_idx_stat::0.161*1.17] distinct minimizers: 611337 (97.89% are singletons); average occurrences: 1.052; average spacing: 5.334; total length: 3428912
[M::worker_pipeline::16.099*2.01] mapped 303815 sequences
[M::main] Version: 2.22-r1101
[M::main] CMD: minimap2 -ax map-ont Medaka_polish/consensus.fasta 2B3_D53_2428_seqs_trimmed.fastq.gz
[M::main] Real time: 16.110 sec; CPU: 32.312 sec; Peak RSS: 1.160 GB

###Output of 2A2_D52 run:
[M::mm_idx_gen::4.506*1.33] collected minimizers
[M::mm_idx_gen::5.373*1.60] sorted minimizers
[M::main::5.375*1.60] loaded/built the index for 20653 target sequence(s)
[M::mm_mapopt_update::5.732*1.56] mid_occ = 34
[M::mm_idx_stat] kmer size: 15; skip: 10; is_hpc: 0; #seq: 20653
[M::mm_idx_stat::5.927*1.54] distinct minimizers: 24662668 (77.10% are singletons); average occurrences: 1.427; average spacing: 5.338; total length: 187914556
[M::main] Version: 2.22-r1101
[M::main] CMD: minimap2 -ax map-ont Medaka_polish4/consensus.fasta Nanopore_2A2_seqs_trimmed.fastq
[M::main] Real time: 1698.012 sec; CPU: 5103.984 sec; Peak RSS: 4.051 GB


```

Sort and index sam/bam files and create depth file 

```
conda activate samtools 

samtools view -b 2A2_polished_seqs.sam -o 2A2_polished_seqs.bam

samtools sort -o 2A2_polished_seqs_sort.bam 2A2_polished_seqs.bam

samtools index 2A2_polished_seqs_sort.bam

#metabat
conda activate metabat
#create depth file
jgi_summarize_bam_contig_depths --outputDepth depth.txt *sort.bam

```

## Binning 
Try different binning methods and compare results:
- medaka consensus fasta without a depth file
- medaka consensus fasta with a depth file
- short read polished medaka consensus fasta without a depth file 


### Racon polish only
```
metabat -i racon4.mod.fasta -o metabat_racon_polish_only/bin_test -v
```

24 bins created. Lots of contamination in large bins (>1800%). Two bins with >75% completeness, <5% contamination, one more bin had 88% completeness, 8% contamination. 



### Medaka polish 
```
metabat -i Medaka_polish4/consensus.fasta -o metabat_medaka_polish/bin_test2 -v
```

28 bins created. 

## Polishing Nanopore assembly with Illumina short reads 




