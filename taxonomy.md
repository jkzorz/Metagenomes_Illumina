## GToTree

Uses a set of single copy marker genes to produce phylogenetic tree. **[(More info here)](https://github.com/AstrobioMike/GToTree/wiki/example-usage#alteromonas-example)**. The default run ended up removing ~300 MAGs because they didn't have enough marker genes. 

```
conda activate gtotree
#version:  GToTree v1.8.3
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

#alignment only - archaea
GToTree -f archaea_list.txt -H Archaea -j 10 -N -n 3 -o gtotree_alignment_archaea

#alignment only - bacteria
GToTree -f bacteria_list.txt -H Bacteria -j 10 -N -n 3 -o gtotree_alignment_bacteria

#with parameters for SCG best hit, IQ-Tree, and keeping MAGs with >30% of SCGs
#bacteria
GToTree -f bacteria_list.txt -H Bacteria -j 10 -n 3 -T IQ-TREE -B -G 0.3 -o gtotree_tree_bacteria2

#archaea
GToTree -f archaea_list.txt -H Archaea -j 10 -n 3 -T IQ-TREE -B -G 0.3 -o gtotree_tree_archaea2


```

Can also just perform the multiple sequence alignment of marker genes by using -N flag. MSA can then be imported into another program (e.g., IQTree). -x Flag overrides the default super5 muscle sequence alignment algorithm used with many genomes. 

**Also ran with drep 98% MAGs**. Will probably use this tree for manuscript: 
```
conda activate gtotree
#version:  GToTree v1.8.3
cd /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/dereplicated_genomes98

#create list of MAG file paths
ls *.fa > all_MAG_list98.txt

#bacteria and archaea
GToTree -f all_MAG_list98.txt -H Bacteria_and_Archaea -j 10 -n 3 -T IQ-TREE -B -G 0.3 -o gtotree_tree_all_bac_arc_drep98

#bac and arch - Hug et al.
GToTree -f all_MAG_list98.txt -H Universal_Hug_et_al -j 10 -n 3 -T IQ-TREE -B -G 0.3 -o gtotree_tree_all_hug_drep98

```
**Pull representative genomes from GTDB to make tree for novel phylum**

```
cd /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/dereplicated_genomes98

#get bacterial gtdb representatives
gtt-get-accessions-from-GTDB -t Bacteria --GTDB-representatives-only

#subset to representatives from each class
gtt-subset-GTDB-accessions -i GTDB-Bacteria-domain-GTDB-rep-metadata.tsv --get-only-individuals-for-the-rank class

```

**Run GToTree with PVC clade MAGs (dereplicated at 90%) and representative GTDB genomes to make tree for novel phylum** 
```
cd /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/dereplicated_genomes98
#slurm script: gtotree_newphylum.slurm

#gototree command: best-hit mode, needs 50% SCG to be included, GTDB taxonomy at Phylum level 
GToTree -f PVC_new_phylum_MAG_list_90.txt -a subset-accessions.txt -H Bacteria -D -L Phylum -j 10 -n 3 -T IQ-TREE -B -G 0.5 -o gtotree_new_phylum_GTDB-drep90
```

## Phyloflash

Use Phyloflash to assign taxonomy to full length 16S rRNA sequences 

```
conda activate phyloflash

phyloFlash.pl -lib Phylo_test_2AT_2428 -poscov -treemap -log -read1 /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/bbduk/cat_qc/JZ-Condor-2AT-175NW-E46-24-28_Li32317_S93_R1_QC.fastq -read2 /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/bbduk/cat_qc/JZ-Condor-2AT-175NW-E46-24-28_Li32317_S93_R2_QC.fastq -readlength 150

```

For loop to run on all samples:
```
#!/bin/bash
###### Reserve computing resources ######
#SBATCH --mail-user=jacqueline.zorz@ucalgary.ca
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=180GB
#SBATCH --time=72:00:00
#SBATCH --partition=cpu2019,cpu2021,cpu2023

###### Set environment variables ######
echo "Starting run at : 'date'"
source /home/jacqueline.zorz/software/miniconda3/etc/profile.d/conda.sh 
conda activate phyloflash

cd /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/taxonomy/PhyloFlash


for x in /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/bbduk/cat_qc/*R1_QC.fastq.gz; 

do 
	R1=$x; 
	R2=$(dirname $x)/$(basename $x R1_QC.fastq.gz)R2_QC.fastq.gz; 
	reads=$(basename $x _R1_QC.fastq.gz); 
	read_short=${reads%_Li*}; #eg JZ-Condor-2AT-700NW-B7-0-4
	read_short2=${read_short:10};  #eg 2AT-700NW-B7-0-4

		
	phyloFlash.pl -lib Phylo_${read_short2} -poscov -treemap -log -read1 $R1 -read2 $R2 -readlength 150 -dbhome /home/jacqueline.zorz/software/miniconda3/envs/phyloflash/lib/phyloFlash/138.1/ -taxlevel 6; done
```

**Collect all assembled Eukaryotic rRNA sequences**

```
#search for Eukaryota in sequence taxonomy "libNAME.phyloFlash.extractedSSUclassifications.csv"
grep "Eukaryota" *extracted* > Eukaryota_rRNA_list.txt

#remove extra bits to get the header found in the fasta files 
awk -F: '{split($2, a, ","); print a[1] "_"}' Eukaryota_rRNA_list.txt > Eukaryota_rRNA_list_short.txt

#use grep to collect Eukaryotic sequences
for i in *.fasta; do grep -A1 -f Eukaryota_rRNA_list_short.txt $i >> Eukaryota_seqs.fasta; done

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

## SingleM

Use SingleM to profile taxonomy of reads in metagenomes. It also provides the microbial content of a sample and can create krona plots. 

```
conda activate singleM

```

For loop: 

```
#!/bin/bash
###### Reserve computing resources ######
#SBATCH --mail-user=jacqueline.zorz@ucalgary.ca
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=180GB
#SBATCH --time=72:00:00
#SBATCH --partition=cpu2019,cpu2021,cpu2023

###### Set environment variables ######
echo "Starting run at : 'date'"
source /home/jacqueline.zorz/software/miniconda3/etc/profile.d/conda.sh 
conda activate singleM

cd /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/taxonomy/singleM


for x in /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/bbduk/cat_qc/*R1_QC.fastq.gz; 

do 
	R1=$x; 
	R2=$(dirname $x)/$(basename $x R1_QC.fastq.gz)R2_QC.fastq.gz; 
	reads=$(basename $x _R1_QC.fastq.gz); 
	read_short=${reads%_Li*}; #eg JZ-Condor-2AT-700NW-B7-0-4
	read_short2=${read_short:10}.tbl;  #eg 2AT-700NW-B7-0-4

		
	singlem pipe -1 $R1 -2 $R2 -p $read_short2 --threads 25 --taxonomic-profile-krona ${read_short2}_krona --metapackage ~/singleM_db/S4.3.0.GTDB_r220.metapackage_20240523.smpkg.zb 

#then run microbial fraction on each profile
	singlem microbial_fraction --forward $R1 --reverse $R2 -p $read_short2 > ${read_short2}_mf.tsv

done

```

Create summary profiles at each taxonomic level of all output tables

```
for i in *.tbl; do singlem summarise --input-taxonomic-profile $i --output-species-by-site-relative-abundance-prefix $(basename $i .tbl); done
```

## Whokaryote

Whokaryote (https://github.com/LottePronk/whokaryote) is used to predict if a contig is prokaryotic or eukaryotic in origin. Running on samples with a high predicted eukaryotic content with hopes to identify contigs that could be used to create a eukaryotic bin.

```
#activate whokaryote environment
conda activate whokaryote

#runs fairly quickly - example with one assembly
cd /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/taxonomy/whokaryote

whokaryote.py --contigs 175NW_contigs/final.contigs.fa --outdir whokaryote_175NW_Illumina_contigs

```

## Eukaryote and Cyanobacteria ASVs

Using phytoref database (http://phytoref.sb-roscoff.fr/downloads) to refine taxonomic classification of eukaryotic and cyanobacteria ASV sequences.

```
#pull all eukaryotic/cyano sequences from ASV fasta file
conda activate seqkit 

grep -f euk_asv_IDs.txt ~/Documents/University\ of\ Calgary/PostDoc/Atlantic\ Condor\ 2021/UofC_Analysis/Dec13/Condor_ASVseqs.fasta > eukaryotic_ASV_sequences.fasta

```

BLAST sequences against phytoref database (Cyanobacteria and PhytoRef databases concatenated)

```
cat phytoref_cyano_sequences.txt phytoref_euk_sequences.txt > phytoref_full_database.fasta

#make blast db of phytoref sequences
makeblastdb -in phytoref_full_database.fasta -out phytoref_full_database.db -dbtype nucl

#blast using evalue and perc_identity cutoff 
lastn -query eukaryotic_ASV_sequences.fasta -db phytoref_full_database.db -outfmt 6 -out blast_phytoref_asvs_results.tbl -evalue 1e-20 -perc_identity 95


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
 


