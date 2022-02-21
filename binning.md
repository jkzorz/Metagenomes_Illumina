# Binning 

Binning with metabat, and using jgi_summarize_contigs to create depth file.
  
Example code with one sample at a time (as they finish mapping stage): 

```
#!/bin/bash
###### Reserve computing resources ######
#SBATCH --mail-user=jacqueline.zorz@ucalgary.ca
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=100GB
#SBATCH --time=3:00:00
#SBATCH --partition=cpu2019,cpu2021,cpu2021-bf24,cpu2019-bf05,cpu2017-bf05

###### Set environment variables ######
echo "Starting run at : 'date'"
source /home/jacqueline.zorz/software/miniconda3/etc/profile.d/conda.sh 
conda activate metabat

cd /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/binning/


summarize depths to make depth file 
jgi_summarize_bam_contig_depths --outputDepth /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/binning/PurpleHaze_04_depth/PurpleHaze_04_depth.txt /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/mapping/mapping_hc_positive/*contig_2A2-PurpleHaze-D52-0-4_sorted.bam


#binning 
metabat -i /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/megahit/megahit_hc_positive/megahit_JZ-Condor-2A2-PurpleHaze-D52-0-4_Li32235_S11/final.contigs.fa -a /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/binning/PurpleHaze_04_depth/PurpleHaze_04_depth.txt -o bin_PurpleHaze_04 -v
```


## CheckM

Example checkM script 
```
#!/bin/bash
###### Reserve computing resources ######
#SBATCH --mail-user=jacqueline.zorz@ucalgary.ca
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=100GB
#SBATCH --time=3:00:00
#SBATCH --partition=cpu2019,cpu2021,cpu2021-bf24,cpu2019-bf05,cpu2017-bf05

###### Set environment variables ######
echo "Starting run at : 'date'"
source /home/jacqueline.zorz/software/miniconda3/etc/profile.d/conda.sh 
conda activate checkm

cd /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/binning/

checkm lineage_wf -f /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/binning/PurpleHaze_04_depth/CheckM_PurpleHaze_04.txt -t 10 -x fa /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/binning/PurpleHaze_04_depth/ /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/binning/PurpleHaze_04_depth/
```


## Gtdbtk
Example Gtdbtk script 

```
#!/bin/bash
###### Reserve computing resources ######
#SBATCH --mail-user=jacqueline.zorz@ucalgary.ca
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=250GB
#SBATCH --time=4:00:00
#SBATCH --partition=bigmem

###### Set environment variables ######
echo "Starting run at : 'date'"
source /home/jacqueline.zorz/software/miniconda3/etc/profile.d/conda.sh 
conda activate gtdbtk

###### Run your script ######

gtdbtk classify_wf --genome_dir /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/binning/PurpleHaze_04_depth/ -x fa --out_dir /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/binning/PurpleHaze_04_depth/gtdbtk_bins --cpus 20


#gtdbtk classify_wf --genome_dir /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/binning/Test_Hole_1216_depth/ -x fa --out_dir /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/binning/Test_Hole_1216_depth/gtdbtk_bins --cpus 20


#gtdbtk classify_wf --genome_dir /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/binning/Test_Hole_2428_depth/ -x fa --out_dir /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/binning/Test_Hole_2428_depth/gtdbtk_bins --cpus 20


##
echo "Job finished with exit code $? at: 'date'"
##

```

## barrnap 

Use barrnap to grab rRNA genes from bins 
```
for i in *.fa; do barrnap $i --outseq barrnap/'rrna_'$i; done
```



## dRep 

Test with bins from 3 samples. Need to put locations of each bin separately into a text file 

```
#e.g. make list of bin locations 
ls Test_Hole_04_depth/*.fa > bin_locations.txt
ls Test_Hole_1216_depth/*.fa >> bin_locations.txt
ls Test_Hole_2428_depth/*.fa >> bin_locations.txt

conda activate drep 

dRep dereplicate drep_test/ -g bin_locations.txt
```

didn't work because of issues with checkm and drep python versions? Can alternatively supply checkM info in a csv:
--genomeInfo GENOMEINFO
                        location of .csv file containing quality information on the genomes. Must contain:
                       "genome"(basename of .fasta file of that genome), "completeness"(0-100 value for completeness of
                        the genome), "contamination"(0-100 value of the contamination of the genome) (default: None)


Ended up copying bin files over to new directory. Will delete later (~4GB)

```
 for i in *_depth/bin*.fa;do cp $i bins_hc_positive/; done
```

