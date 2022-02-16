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
