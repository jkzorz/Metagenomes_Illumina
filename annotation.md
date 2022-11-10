# Annotation

## DRAM annotation 

**DRAM version: 1.3**

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

Example command for concatenating files from all dram output folders: 
```
for i in good_bin_annotations*/annotations.tsv; do cat $i >> all_annotations.tsv; done
for i in good_bin_annotations*/trnas.tsv; do cat $i >> all_trnas.tsv; done
for i in good_bin_annotations*/rrnas.tsv; do cat $i >> all_rrnas.tsv; done

#distill
DRAM.py distill -i all_annotations.tsv -o genome_summaries --trna_path all_trnas.tsv --rrna_path all_rrnas.tsv
```

**Das Tool and dRep bins**

1750 bins, divided into ~100 bin folders to run DRAM separately: 
Made individual dram folders (dram1, dram2, dram3, etc) manually

*copy_bins_command.sh*
```
dir=1
counter=1

for file in /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/binning/dastool/drep_dastool_out/dereplicated_genomes/*.fa
do
   cp $file dram$dir/
   ((counter++))
   (( $counter%100 == 1 )) && ((dir++))
done

```




## Barrnap 
Use barrnap to grab rRNA genes from bins 

**barrnap version: 0.9**

```
for i in *.fa; do barrnap $i --outseq barrnap_16S/'rrna_'$i; done
```

das tool and drep bins:

```
for i in /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/binning/dastool/drep_dastool_out/dereplicated_genomes/*.fa; do barrnap $i --outseq /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/annotation/barrnap_drep_dastool/'rrna_'$(basename $i); done
```

Collect all MAG 16S sequences and append MAG name to fasta header 
```
for i in barrnap_drep_dastool/rrna_*.fa; do sample="$(basename $i .fa)"; echo -n $(grep ">16S" $i) >> test.fa && echo ":${sample}" >> test.fa; grep "16S" $i -A 1 | tail -n 1 >> test.fa ; done

#get rid of MAGs without 16S
sed 's/^:.*//' test.fa > test2.fa

#get rid of empty lines
sed '/^$/d' test2.fa > 16S_MAG_sequences_4db.fa
```

### Blast 16S Miseq sequences against MAGs

```
makeblastdb -in 16S_MAG_sequences_4db.fa -out 16S_MAG_sequences_4db.db -dbtype nucl

blastn -query may17_ASVseqs.fasta -db 16S_MAG_sequences_4db.db -outfmt 6 -out blast_16S_results.tbl -max_target_seqs 1 -perc_identity 97 -word_size 200

#allow for multiple hits of ASV to MAG 16S
#blastn -query may17_ASVseqs.fasta -db 16S_MAG_sequences_4db.db -outfmt 6 -out blast_16S_results2.tbl -max_target_seqs 5 -perc_identity 99.5 -word_size 200

#word_size doesn't mean length of alignment, but length of "word" used for blast search 
blastn -query may17_ASVseqs.fasta -db 16S_MAG_sequences_4db.db -outfmt 6 -out blast_MAG16S_May17_results_99.tbl -max_target_seqs 100 -perc_identity 99


#allow for shorter matches because some of the 16S sequences aren't completely full length
#blastn -query may17_ASVseqs.fasta -db 16S_MAG_sequences_4db.db -outfmt 6 -out blast_16S_results3.tbl -max_target_seqs 5 -perc_identity 97

```




## CANT-HYD HMMs


Test for loop to run CANT-HYD HMMs on Purple Haze 0-4 cm bins. 
```
for i in bins/bin_PurpleHaze_04.*/genes.faa; do  hmmsearch --tblout hmmsearch_$(basename $(dirname $i)).tblout ../CANT-HYD.hmm $i > $(basename $(dirname $i)).out;done

```

Try on all good MAGs

```
conda activate hmmer
cd /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/binning/drep_out/dereplicated_genomes/good_bins

#copy over the hmms
cp /work/ebg_lab/gm/CANT-HYD/Final_important_files_JZ/CANT-HYD.hmm .

#copy over prodigal predicted genes for the good bins 
for i in *.fa; do cp -R ../../bins/$(basename $i .fa) . ; done

#cant hyd for loop - use noise cutoff
for i in predicted_genes/bin*/genes.faa; do  hmmsearch --cut_nc --tblout hmmsearch_$(basename $(dirname $i)).tblout CANT-HYD.hmm $i > $(basename $(dirname $i)).out;done

```

All das tool and drep MAGs: 

```
for i in /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/binning/dastool/drep_dastool_out/dereplicated_genomes/*.fa; do hmmsearch --cut_nc --tblout hmmsearch_$(basename $i .fa).tblout CANT-HYD.hmm /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/binning/dastool/drep_dastool_out/bins/$(basename $i .fa)/genes.faa > $(basename $i).out;done
```



## MetaErg 2.0

docker link: https://hub.docker.com/r/kinestetika/metaerg
github: https://github.com/kinestetika/MetaErg

Docker is not for HPC clusters, singularity is generally used instead. But singularity can be used to pull docker images and convert to singularity images (sif files)
```
module load singularity
singularity pull docker://kinestetika/metaerg
```

Can use 'singularity run' to enter singularity container for metaerg, and run metaerg scripts: 

```
singularity run /work/ebg_lab/referenceDatabases/metaerg_latest.sif
```

Use this to download metaerg databases. For some reason, I only seem to be able to access my home directory from the singularity container? 
```
metaerg --download_database --database_dir metaerg_test/
```



