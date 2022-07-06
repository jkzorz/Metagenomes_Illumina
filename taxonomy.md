## Phyloflash

Use Phyloflash to assign taxonomy to full length 16S rRNA sequences 

```
conda activate phyloflash

phyloFlash.pl -lib Phylo_test_2AT_2428 -poscov -treemap -log -read1 /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/bbduk/cat_qc/JZ-Condor-2AT-175NW-E46-24-28_Li32317_S93_R1_QC.fastq -read2 /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/bbduk/cat_qc/JZ-Condor-2AT-175NW-E46-24-28_Li32317_S93_R2_QC.fastq -readlength 150

```





## Kraken 

Try assigning taxonomy to contigs/reads with Kraken 



```
#!/bin/bash
###### Reserve computing resources ######
#SBATCH --mail-user=jacqueline.zorz@ucalgary.ca
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=180GB
#SBATCH --time=5:00:00
#SBATCH --partition=cpu2019,cpu2021,cpu2021-bf24,bigmem,cpu2019-bf05

###### Set environment variables ######
echo "Starting run at : 'date'"
source /home/jacqueline.zorz/software/miniconda3/etc/profile.d/conda.sh 
conda activate kraken2

cd /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/taxonomy



#kraken2 --paired --db /work/ebg_lab/referenceDatabases/KrakenGTDB --threads 40 --classified-out PurplePatch2428_pair_gtdb.kraken.classified#.fasta --unclassified-out PurplePatch2428_pair_gtdb#.kraken.unclassified.fasta --output PurplePatch2428_pair_gtdb.kraken.output --report PurplePatch2428_pair_gtdb.kraken.report /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/bbduk/cat_qc/JZ-Condor-2B1-PurplePatch-A54-24-28_Li32230_S6_R1_QC.fastq /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/bbduk/cat_qc/JZ-Condor-2B1-PurplePatch-A54-24-28_Li32230_S6_R2_QC.fastq

#kraken2 --paired --db /work/ebg_lab/referenceDatabases/KrakenGTDB --threads 40 --classified-out BG15_04_pair.kraken.classified#.fasta --unclassified-out BG15_04_pair#.kraken.unclassified.fasta --output BG15_04_pair.kraken.output --report BG15_04_pair.kraken.report /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/bbduk/cat_qc/JZ-Condor-BG15-Background-23-0-4_Li32305_S81_R1_QC.fastq /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/bbduk/cat_qc/JZ-Condor-BG15-Background-23-0-4_Li32305_S81_R2_QC.fastq

#kraken2 --db /work/ebg_lab/referenceDatabases/KrakenGTDB --threads 40 --classified-out PurplePatch2428_contig_gtdb.kraken.classified.fasta --unclassified-out PurplePatch2428_contig_gtdb.kraken.unclassified.fasta --output PurplePatch2428_contig_gtdb.kraken.output --report PurplePatch2428_contig_gtdb.kraken.report /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/megahit/megahit_hc_positive/megahit_JZ-Condor-2B1-PurplePatch-A54-24-28_Li32230_S6/final.contigs.fa


#kraken2 --db /work/ebg_lab/referenceDatabases/KrakenGTDB --threads 40 --classified-out PurplePatch2428_contig_gtdb_unbinned.kraken.classified.fasta --unclassified-out PurplePatch2428_contig_gtdb_unbinned.kraken.unclassified.fasta --output PurplePatch2428_contig_gtdb_unbinned.kraken.output --report PurplePatch2428_contig_gtdb_unbinned.kraken.report /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/binning/PurplePatch_2428_depth/unbinned.fna


kraken2 --db /work/ebg_lab/referenceDatabases/KrakenGTDB --threads 40 --classified-out 2B3_D53_2428_nanopore_gtdb.kraken.classified.fasta --unclassified-out 2B3_D53_2428_nanopore.kraken.unclassified.fasta --output 2B3_D53_2428_nanopore.kraken.output --report 2B3_D53_2428_nanopore.kraken.report /work/ebg_lab/gm/gapp/jzorz/2B3_2428_D53_fastq_pass/2B3_D53_2428_seqs_trimmed.fastq.gz

kraken2 --db /work/ebg_lab/referenceDatabases/KrakenGTDB --threads 40 --classified-out 2B3_D53_2428_nanopore_contig_gtdb.kraken.classified.fasta --unclassified-out 2B3_D53_2428_nanopore_contig.kraken.unclassified.fasta --output 2B3_D53_2428_nanopore_contig.kraken.output --report 2B3_D53_2428_nanopore_contig.kraken.report /work/ebg_lab/gm/gapp/jzorz/2B3_2428_D53_fastq_pass/Medaka_polish/consensus.fasta
```

