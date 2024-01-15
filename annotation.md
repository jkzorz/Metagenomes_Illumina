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
DRAM.py distill -i all_annotations.tsv -o genome_summaries --trna_path all_trnas.tsv --rrna_path all_rrnas.tsv --genomes_per_product 250
```

DRAM distill heatmap information: https://github.com/WrightonLabCSU/DRAM/blob/master/data/function_heatmap_form.tsv


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


## Kofamscan

Use kofamscan to annotate proteins with kofam HMMs

**Downloading kofamscan**
```
#create kofamscan conda environment
conda create -n kofamscan
conda activate kofamscan
#install dependencies with conda
conda install ruby
conda install hmmer
conda install parallel

#create kofam directory in home 
mkdir ~/kofamscan
cd ~/kofamscan

#download profiles, ko lists, and kofam-scan 
wget ftp://ftp.genome.jp/pub/tools/kofamscan/kofamscan.tar.gz
wget ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz
gunzip ko_list.gz
tar xf profiles.tar.gz
wget https://www.genome.jp/ftp/tools/kofam_scan/kofam_scan-1.3.0.tar.gz
tar xf kofam_scan-1.3.0.tar.gz

#update config file with locations of profile and ko list (optional, can also be provided as parameters in command):
profile: /home/jacqueline.zorz/kofamscan/profiles
ko_list: /home/jacqueline.zorz/kofamscan/ko_list
```
**Running kofamscan:**

```
 ~/kofamscan/kofam_scan-1.3.0/exec_annotation -o kofam_test.txt --cpu 20 --tmp-dir ko-fam_temp -E 0.0001 -f detail-tsv ../checkm2_output_mags/protein_files/concoct_2B1-TinyBubbles-C18-0-4_106_sub.faa
```

**kofamscan for loop with slurm:**

```
#!/bin/bash
###### Reserve computing resources ######
#SBATCH --mail-user=jacqueline.zorz@ucalgary.ca
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=180GB
#SBATCH --time=24:00:00
#SBATCH --partition=bigmem,cpu2019,cpu2021

#set environmental parameters

source /home/jacqueline.zorz/software/miniconda3/etc/profile.d/conda.sh 
conda activate kofamscan

cd /work/ebg_lab/gm/gapp/taylor/kofam_hmms

for i in /work/ebg_lab/gm/gapp/taylor/checkm2_output_mags/protein_files/*.faa; do mag=$(basename $i .faa); ~/kofamscan/kofam_scan-1.3.0/exec_annotation -o kofam_$mag.txt --cpu 20 --tmp-dir ko-fam_temp -E 0.0001 -f detail-tsv $i; rm ko-fam_temp/ -R; done
```

**kofamscan on all final MAGs: detail-tsv version**
```
#!/bin/bash
###### Reserve computing resources ######
#SBATCH --mail-user=jacqueline.zorz@ucalgary.ca
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=180GB
#SBATCH --time=168:00:00
#SBATCH --partition=cpu2019,cpu2021

#set environmental parameters

source /home/jacqueline.zorz/software/miniconda3/etc/profile.d/conda.sh 
conda activate kofamscan

cd /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/binning/dastool/drep_dastool_out2/dereplicated_genomes/Good_Bins_CheckM2/drep_Final_Nanopore_output/dereplicated_genomes/kofamscan/

for i in /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/binning/dastool/drep_dastool_out2/dereplicated_genomes/Good_Bins_CheckM2/drep_Final_Nanopore_output/dereplicated_genomes/checkm2_output_all/protein_files/*.faa; do mag=$(basename $i .faa); ~/kofamscan/kofam_scan-1.3.0/exec_annotation -o kofam_$mag.txt --cpu 20 --tmp-dir ko-fam_temp -E 0.0001 -f detail-tsv $i; rm ko-fam_temp/ -R; done
```
Results now located at: ```/work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/dereplicated_genomes/kofamscan```

**kofamscan on all final MAGs: mapper version for input to kegg mapper**
```
#!/bin/bash
###### Reserve computing resources ######
#SBATCH --mail-user=jacqueline.zorz@ucalgary.ca
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=180GB
#SBATCH --time=168:00:00
#SBATCH --partition=cpu2019,cpu2021

#set environmental parameters

source /home/jacqueline.zorz/software/miniconda3/etc/profile.d/conda.sh 
conda activate kofamscan

cd /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/binning/dastool/drep_dastool_out2/dereplicated_genomes/Good_Bins_CheckM2/drep_Final_Nanopore_output/dereplicated_genomes/kofamscan_mapper/

for i in /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/binning/dastool/drep_dastool_out2/dereplicated_genomes/Good_Bins_CheckM2/drep_Final_Nanopore_output/dereplicated_genomes/checkm2_output_all/protein_files/*.faa; do mag=$(basename $i .faa); ~/kofamscan/kofam_scan-1.3.0/exec_annotation -o kofam_mapper_$mag.txt --cpu 20 --tmp-dir ko-fam_temp -E 0.0001 -f mapper $i; rm ko-fam_temp/ -R; done

```

Results now located at: ```/work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/dereplicated_genomes/kofamscan_mapper```

## Genomad

Use Genomad to identify plasmids and viruses in MAGs

```
conda activate genomad

cd /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/dereplicated_genomes98

```
Bash script: 

```
#!/bin/bash
###### Reserve computing resources ######
#SBATCH --mail-user=jacqueline.zorz@ucalgary.ca
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=180GB
#SBATCH --time=48:00:00
#SBATCH --partition=bigmem,cpu2019,cpu2021,cpu2021-bf24


###### Set environment variables ######
echo "Starting run at : 'date'"
source /home/jacqueline.zorz/software/miniconda3/etc/profile.d/conda.sh 
conda activate genomad

cd /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/dereplicated_genomes98/genomad


for i in /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/dereplicated_genomes98/*.fa; do mag=$(basename $i); genomad end-to-end --cleanup --splits 8 --enable-score-calibration $i $mag.genomad_default /work/ebg_lab/referenceDatabases/genomad/genomad_db/; done

```
Score calibration not needed if sample has fewer than 1000 sequences (https://portal.nersc.gov/genomad/score_calibration.html) 

```
#collect all virus and plasmid sequences:
for i in *genomad_default/*summary/*_virus_proteins.faa; do cat $i >> all_virus_proteins.faa; done
for i in *genomad_default/*summary/*_plasmid_proteins.faa; do cat $i >> all_plasmid_proteins.faa; done

#collect just gene names
grep ">" all_virus_proteins.faa > virus_proteins_list.txt
grep ">" all_plasmid_proteins.faa > plasmid_proteins_list.txt

#add MAG names to genomad protein names:
cd /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/dereplicated_genomes98/genomad/
for i in *genomad_default/*_annotate/*_proteins.faa; do mag=$(basename $i _proteins.faa); sed "s/>/>${mag}:/g" $i > ${mag}_2.faa; done

#collect all gene names with MAG name in header and put in new file:
for i in *_2.faa; do grep ">" --no-group-separator $i >> MAG_contig_gene.txt;done
#remove everything after the space to make file smaller
sed 's/ .*//g' MAG_contig_gene.txt > MAG_contig_gene2.txt
```

## dbCAN3 

The latest dbCAN version, with cgc (cazyme gene cluster) clusters and substrate specific dbcan database: (https://dbcan.readthedocs.io/en/latest/index.html)

```
#!/bin/bash
###### Reserve computing resources ######
#SBATCH --mail-user=jacqueline.zorz@ucalgary.ca
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=180GB
#SBATCH --time=60:00:00
#SBATCH --partition=cpu2019,cpu2021

###### Set environment variables ######
echo "Starting run at : 'date'"
source /home/jacqueline.zorz/software/miniconda3/etc/profile.d/conda.sh 
conda activate dbcan

cd /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/dereplicated_genomes98/dbcan


for i in /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/dereplicated_genomes98/*.fa; do mag=$(basename $i); run_dbcan $i prok -c cluster --out_dir dbcan_$mag --db_dir /work/ebg_lab/referenceDatabases/dbcan_db/; done
```


## Barrnap 
Use barrnap to grab rRNA genes from bins 

**barrnap version: 0.9**

```
#e.g.
for i in *.fa; do barrnap $i --outseq barrnap_16S/'rrna_'$i; done
```

das tool and drep bins (update to checkm2 bins):

```
conda activate barrnap
cd /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/annotation/barrnap_drep_dastool

for i in /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/binning/dastool/drep_dastool_out2/dereplicated_genomes/*.fa; do barrnap $i --outseq /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/annotation/barrnap_drep_dastool/'rrna_'$(basename $i); done
```

Collect all MAG 16S sequences and append MAG name to fasta header 
```
for i in barrnap_drep_dastool/rrna_*.fa; do sample="$(basename $i .fa)"; echo -n $(grep ">16S" $i) >> test.fa && echo ":${sample}" >> test.fa; grep "16S" $i -A 1 | tail -n 1 >> test.fa ; done

#get rid of MAGs without 16S
sed 's/^:.*//' test.fa > test2.fa

#get rid of empty lines
sed '/^$/d' test2.fa > 16S_MAG_sequences_4db.fa
```

**Python script** (python_16S_parse2.py) to collect 16S rRNA genes from barrnap output and put in new files. 
```
#! /usr/bin/python

import os, re
from itertools import islice 

word = '16S'

for file in os.listdir("/work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/annotation/barrnap_drep_dastool"):
	if file.endswith(".fa"):
		#search for 16S
		fasta = open(file, 'r').read().splitlines()
			#read lines in list
		for index, line in enumerate(fasta):
			if re.search(word, line):
				#create new file 
				new = ["header", file]
				new2 = "_".join(new)
				with open(new2, 'a') as f:
					temp = fasta[index:index+2]
					print('\n'.join(temp), file=f)
```


**Python script** (change_header.py) to add MAG name to fasta header
```
#! /usr/bin/python

import os, re
from itertools import islice 

word = "16S"

for file in os.listdir("/work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/annotation/barrnap_drep_dastool"):
	if file.startswith("header"):
		input = open(file, 'r')
		file2 = re.sub('.fa', '', file)
		file3 = re.sub('header_rrna_', '', file2)
		temp = ['temp', file]
		temp2 = open('_'.join(temp), 'w')
		for line in input:
			text = ['16S',file3]
			if word in line: 
				temp2.write(re.sub(word, '_'.join(text), line))
			else: 
				temp2.write(line)
			
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

## Prodigal gene prediction 

Use Prodigal to predict genes in MAGs

```
cd /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/dereplicated_genomes98

for i in *.fa; do prodigal -a drep98_proteins/${i}a -d drep98_proteins/${i}sta -i $i -o drep98_proteins/$(basename $i .fa).out; done
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

#cant hyd for loop - noise cutoff and drep 98 MAGs
cd /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/dereplicated_genomes98/drep98_proteins
for i in genes_protein/*.faa; do  hmmsearch --cut_nc --tblout hmmsearch_$(basename $i).tblout CANT-HYD.hmm $i > $(basename $i).out;done

```

All das tool and drep MAGs: 

```
for i in /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/binning/dastool/drep_dastool_out/dereplicated_genomes/*.fa; do hmmsearch --cut_nc --tblout hmmsearch_$(basename $i .fa).tblout CANT-HYD.hmm /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/binning/dastool/drep_dastool_out/bins/$(basename $i .fa)/genes.faa > $(basename $i).out;done
```

## HMM search 
Using TIGRfam HMM to search for reductive dehalogenase in MAGs
```
conda activate hmmer

for i in /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/binning/dastool/drep_dastool_out2/checkm2_output/protein_files/*.faa; do hmmsearch --cut_nc --tblout rdhA_hmmsearch_$(basename $i .faa).tblout TIGR02486.1.HMM $i > $(basename $i .faa).out;done

```

### Python scripts for managing HMM output
From all files, grab rows with contig info
```
#! /usr/bin/python

import os, re
from itertools import islice 

word = 'k141_*'

for file in os.listdir("/work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/annotation/Jayne_rdhA"):
	if file.endswith(".tblout"):
		#search for rows with content
		fasta = open(file, 'r').readlines()
		#read lines in list
		for line in fasta:
			if re.search(word, line):
				#create new file 
				new = ["summary", file]
				new2 = "_".join(new)
				with open(new2, 'a') as f:
					print(line, file=f)
```

From summary files, add file name to final column, overwrite summary file

```
#! /usr/bin/python

import os, re, pandas
from itertools import islice

word = 'k141_*'

for file in os.listdir("/work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/annotation/Jayne_rdhA"):
        if file.startswith("summary"):
                df = pandas.read_csv(file, header=None)
                df[len(df.columns)] = file
                df.to_csv(file, header=None, index=None)

```

### Pfam HMMs

Need to first use hmmfetch to grab specific HMM from all Pfam HMMs (Pfam-A.HMM file)

```
hmmfetch /work/ebg_lab/referenceDatabases/Pfam/Pfam-A.hmm PF13486.8 > /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/annotation/Jayne_rdhA/PF13486.8.hmm 
```

Search genes against Pfam HMM
```
for i in /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/binning/dastool/drep_dastool_out2/checkm2_output/protein_files/*.faa; do hmmsearch --cut_nc --tblout Pfam_hmmsearch_$(basename $i .faa).tblout PF13486.8.hmm $i > $(basename $i .faa).out;done
```
Then use the python scripts above to manipulate data. 

### Pull gene sequences with HMM match

Pull gene sequences that match TIGRfam or Pfam HMM 
```
cd /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/annotation/Jayne_rdhA/

#change , to space in csv file for easier use with awk 
sed 's/,/ /g' PF13486/MAG_info_PF13486.csv > MAG_info_PF13486_update.csv
sed 's/,/ /g' TIGR02486/MAG_info_rdhA_TIGR02486.csv > MAG_info_rdhA_TIGR02486_update.csv

#make list of genes from summarized output
awk '{ print $NF"_"$1" " }' MAG_info_rdhA_TIGR02486_update.csv > TIGR02486_genes2.list
awk '{ print $NF"_"$1" " }' MAG_info_PF13486_update.csv > PF13486_genes2.list

#fix list to match contig headers
sed 's/summary_rdhA_hmmsearch_//' TIGR02486_genes2.list > TIGR02486_genes3.list
sed 's/.tblout//' TIGR02486_genes3.list > TIGR02486_genes4.list
sed 's/summary_Pfam_hmmsearch_//' PF13486_genes2.list > PF13486_genes3.list
sed 's/.tblout//' PF13486_genes3.list > PF13486_genes4.list

#use grep to find genes in DRAM files 
for i in /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/annotation/dram_drep_dastool/dram*/genes.faa; do grep -Ff TIGR02486_genes.list $i -A1 --no-group-separator >> TIGR02486_genes.fasta; done

for i in /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/annotation/dram_drep_dastool/dram*/genes.faa; do grep -Ff PF13486_genes.list $i -A1 --no-group-separator >> PF13486_genes.fasta; done
```

DRAM files don't contain all the genes from some of the MAGs from the updated checkm2 list. Try to pull genes from updated checkm2 list. First need to convert files to single line fasta, and put MAG names in gene headers. 

**Note: Duplicate contig/gene names in some MAGs... 

```
#change headers
for i in /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/binning/dastool/drep_dastool_out2/checkm2_output/protein_files/*.faa; do sample=$(basename $i .faa); sed "s/>/>${sample}_/" $i > /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/binning/dastool/drep_dastool_out2/checkm2_output/protein_files/headers/${sample}.faa; done

#convert to single line
for i in /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/binning/dastool/drep_dastool_out2/checkm2_output/protein_files/headers/*.faa; do awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' $i > /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/binning/dastool/drep_dastool_out2/checkm2_output/protein_files/headers/singleline/$(basename $i); done

#repeat search - pfam
for i in /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/binning/dastool/drep_dastool_out2/checkm2_output/protein_files/headers/*.faa; do grep -Ff PF13486_genes4.list $i -A1 --no-group-separator >> PF13486_genes_checkm2.fasta; done

#repeat search - tigrfam
for i in /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/binning/dastool/drep_dastool_out2/checkm2_output/protein_files/headers/*.faa; do grep -Ff TIGR02486_genes4.list $i -A1 --no-group-separator >> TIGR02486_genes_checkm2.fasta; done


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
singularity run ~/metaerg_latest.sif
```

Use this to download metaerg databases. For some reason, I only seem to be able to access my home directory from the singularity container? 
```
mkdir metaerg_database
metaerg --download_database --database_dir metaerg_database/
metaerg --create_database S --database_dir  metaerg_database/
```

Running Metaerg on one genome sample: 
```
metaerg --contig_file bin_Hole_1216_85.fa --database_dir ../metaerg_database/ --path_to_signalp ../metaerg_database/signalp-5.0b.Linux.tar.gz --path_to_tmhmm ../metaerg_database/tmhmm-2.0c.Linux.tar.gz
```

Running Metaerg on multiple genomes: 
```
 metaerg --contig_file ~/metaerg_abyssubacteria/ --database_dir ../metaerg_database/ --path_to_signalp ../metaerg_database/signalp-6.0g.fast.tar.gz --path_to_tmhmm ../metaerg_database/tmhmm-2.0c.Linux.tar.gz --file_extension .fa --rename_contigs --rename_genomes
```

Running a Singularity command in batch on slurm: 
```
#!/bin/bash
###### Reserve computing resources ######
#SBATCH --mail-user=jacqueline.zorz@ucalgary.ca
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=100GB
#SBATCH --time=8:00:00
#SBATCH --partition=cpu2019,cpu2021,cpu2021-bf24,bigmem


#might not be necessary... 
source /home/jacqueline.zorz/software/miniconda3/etc/profile.d/conda.sh 
conda activate metaerg2

cd ~/metaerg_abyssubacteria

singularity run ../metaerg_latest.sif metaerg --contig_file ~/metaerg_abyssubacteria/ --database_dir ~/metaerg_database/ --path_to_signalp ~/metaerg_database/signalp-6.0g.fast.tar.gz --path_to_tmhmm ~/metaerg_database/tmhmm-2.0c.Linux.tar.gz --file_extension .fa --rename_contigs --rename_genomes
```

