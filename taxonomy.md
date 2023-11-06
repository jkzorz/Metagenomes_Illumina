## GToTree

Uses a set of single copy marker genes to produce phylogenetic tree. **[More info here]**(https://github.com/AstrobioMike/GToTree/wiki/example-usage#alteromonas-example). The default run ended up removing ~300 MAGs because they didn't have enough marker genes. 

```
conda activate gtotree
cd /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/dereplicated_genomes

#shows list of in-built HMMs
gtt-hmms

#create list of MAG file paths
ls *.fa > Final_MAG_List.txt

#basic run 
GToTree -f Final_MAG_List.txt -H Bacteria_and_Archaea -j 15 -o gtotree_test

#bacteria only
GToTree -f Final_MAG_List.txt -H Bacteria -j 15 -o gtotree_bacteria_only_test

#archaea only
GToTree -f Final_MAG_List.txt -H Archaea -j 15 -o gtotree_archaea_only_test
```

Can also just perform the multiple sequence alignment of marker genes by using -N flag. MSA can then be imported into another program (e.g., IQTree). -x Flag overrides the default super5 muscle sequence alignment algorithm used with many genomes. 


## Phyloflash

Use Phyloflash to assign taxonomy to full length 16S rRNA sequences 

```
conda activate phyloflash

phyloFlash.pl -lib Phylo_test_2AT_2428 -poscov -treemap -log -read1 /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/bbduk/cat_qc/JZ-Condor-2AT-175NW-E46-24-28_Li32317_S93_R1_QC.fastq -read2 /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/bbduk/cat_qc/JZ-Condor-2AT-175NW-E46-24-28_Li32317_S93_R2_QC.fastq -readlength 150

```


## FastTree 

Use FastTree to make tree of gtdbtk.ar53.user_msa.fasta and gtdbtk.bac120.user_msa.fasta files. (Needed to gunzip file first) 

```
#archaea tree
FastTree < gtdbtk.ar53.user_msa.fasta > gtdbtk.ar53.user_msa.tre

#bacteria tree
FastTree < gtdbtk.bac120.user_msa.fasta > gtdbtk.bac120.user_msa.tre

```

## FastANI

Use FastANI to determine relatedness between genomes 

```
conda activate fastani

#create list of MAGs to compare
ls > MAG_list.txt #need to delete the MAG_list.txt line from file

#run fastANI - all MAGs against all MAGs
fastANI --ql MAG_list.txt --rl MAG_list.txt -o MAG_FastANI.tsv

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

## inStrain
Trying inStrain to use with reference Atribacteria genome from gtdb (https://gtdb.ecogenomic.org/genome?gid=GCA_001773955.1). Same genus as the poor quality MAGs recovered, but different species. Installation was slow with Conda, so used Mamba instead. 

Use Hole 24-28cm sample (high Atribacteria content according to 16S). First need to map reads to gtdb Atribacteira genome contig file.

```
conda activate bbtools

bbmap.sh ref=gtdb_atribacteria_GCA_001773955.1_ASM177395v1_genomic.fna in=/work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/bbduk/cat_qc/JZ-Condor-2A1-TheHole-C54-24-28_Li32297_S73_R1_QC.fastq in2=/work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/bbduk/cat_qc/JZ-Condor-2A1-TheHole-C54-24-28_Li32297_S73_R2_QC.fastq outm=Hole2428_instrain.sam covstats=covstats_instrain_test.txt scafstats=scafstats_instrain_test.txt threads=40 minid=0.
 ```
 
 convert sam file to sorted bam file 
 
 ```
 conda activate samtools
 
 samtools view -b Hole2428_instrain.sam | samtools sort -o Hole2428_instrain_sort.bam
 ```
 
 Run inStrain profile function. Have to specify the use full header parameter because bbmap doesn't shorten contig headers in mapping files. 

```
 conda activate instrain 
 
 inStrain profile Hole2428_instrain_sort.bam gtdb_atribacteria_GCA_001773955.1_ASM177395v1_genomic.fna -o gtdb_atribacteria_instrain_test --use_full_fasta_header
```
 


