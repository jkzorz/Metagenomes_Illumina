# BBduk QC 

Had java issues with original bbtools that I had previously downloaded. Created new conda environment and re-downloaded bbtools re online instructions (https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/installation-guide/)

```
#new conda environment
conda create -n bbduk
conda activate bbduk 

#downloaded manually from here: https://sourceforge.net/projects/bbmap/
#transferred file over to home directory on server
#ran unzip command
tar -xvzf BBMap_(version).tar.gz
```

Created for loop to loop through all JZ* samples in metagenome directory and ran in slurm batch script: 
```
#!/bin/bash
###### Reserve computing resources ######
#SBATCH --mail-user=jacqueline.zorz@ucalgary.ca
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=100GB
#SBATCH --time=12:00:00
#SBATCH --partition=bigmem,cpu2019,cpu2021,cpu2021-bf24

###### Set environment variables ######
echo "Starting run at : 'date'"
source /home/jacqueline.zorz/software/miniconda3/etc/profile.d/conda.sh 
conda activate bbduk

cd /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/bbduk/

###### Run your script ######

#step 1
for f in /work/ebg_lab/gm/Metagenomes_2022/GAPP2/JZ*/*R1_001.fastq.gz;  
do 
	R1=$f
	R2=$(dirname $f)/$(basename $f R1_001.fastq.gz)R2_001.fastq.gz
	sample=$(basename $(dirname $f))
	base=$(basename $f _R1_001.fastq.gz)
	
	#step 1
	~/bbmap/bbduk.sh in=$R1 in2=$R2 out=${sample}_${base}_R1_lastbase_rm.fastq out2=${sample}_${base}_R2_lastbase_rm.fastq ftm=5;
	
	#step 2
	~/bbmap/bbduk.sh in=${sample}_${base}_R1_lastbase_rm.fastq in2=${sample}_${base}_R2_lastbase_rm.fastq out=${sample}_${base}_R1_adapter_rm.fastq out2=${sample}_${base}_R2_adapter_rm.fastq ref=/home/jacqueline.zorz/software/miniconda3/envs/bbtools/bbtools/lib/resources/adapters.fa tbo tpe k=23 mink=11 hdist=1 ktrim=r;
	
	#step3
	~/bbmap/bbduk.sh in=${sample}_${base}_R1_adapter_rm.fastq in2=${sample}_${base}_R2_adapter_rm.fastq out=${sample}_${base}_R1_dec.fastq out2=${sample}_${base}_R2_dec.fastq outm=${sample}_${base}_contaminants.fq ref=/home/jacqueline.zorz/software/miniconda3/envs/bbtools/bbtools/lib/resources/phix_adapters.fa.gz k=31 hdist=1 stats=${sample}_${base}_stats.txt;
	
	#step4
	~/bbmap/bbduk.sh in=${sample}_${base}_R1_dec.fastq in2=${sample}_${base}_R2_dec.fastq out=${sample}_${base}_R1_qc.fastq out2=${sample}_${base}_R2_qc.fastq qtrim=rl trimq=15 minlength=30 entropy=0.5;

	#remove inbetween files 
	rm ${sample}_${base}_R1_lastbase_rm.fastq
	rm ${sample}_${base}_R2_lastbase_rm.fastq
	rm ${sample}_${base}_R1_adapter_rm.fastq
	rm ${sample}_${base}_R2_adapter_rm.fastq
	rm ${sample}_${base}_R1_dec.fastq
	rm ${sample}_${base}_R2_dec.fastq

done

###
```

## Concatenate resulting 4 qc files into one large file per sample (R1, R2) 

```
#!/bin/bash
###### Reserve computing resources ######
#SBATCH --mail-user=jacqueline.zorz@ucalgary.ca
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=100GB
#SBATCH --time=12:00:00
#SBATCH --partition=bigmem,cpu2019,cpu2021,cpu2021-bf24


#for concatenating qc'ed reads - R1
for i in JZ*L001_R1_qc.fastq; do cat $i $(basename $i L001_R1_qc.fastq)L00{2,3,4}_R1_qc.fastq  > cat_qc/$(basename $i L001_R1_qc.fastq)R1_QC.fastq; done

#for concatenating qc'ed reads - R2
for i in JZ*L001_R2_qc.fastq; do cat $i $(basename $i L001_R2_qc.fastq)L00{2,3,4}_R2_qc.fastq  > cat_qc/$(basename $i L001_R2_qc.fastq)R2_QC.fastq; done
```



## Trying different methods to improve downstream assembly and binning

Trying different methods like subsampling and normalization of read coverage to improve downstream assembly and binning. 

### Subsampling

Subsampled reads to 10% of original
```
conda activate bbtools

#subsample
reformat.sh in=cat_qc/JZ-Condor-2B1-PurplePatch-A54-0-4_Li32225_S1_R1_QC.fastq in2=cat_qc/JZ-Condor-2B1-PurplePatch-A54-0-4_Li32225_S1_R2_QC.fastq out=JZ-Condor-2B1-PurplePatch-A54-0-4_Li32225_S1_R1_QC_subsample10pc.fastq out2=JZ-Condor-2B1-PurplePatch-A54-0-4_Li32225_S1_R2_QC_subsample10pc.fastq samplerate=0.1

```
Output: 
```
Input is being processed as paired
Input:                          225676410 reads                 31544289481 bases
Processed:                      22569544 reads                  3154567908 bases
Output:                         22569544 reads (10.00%)         3154567908 bases (10.00%)

Time:                           302.843 seconds.
Reads Processed:      22569k    74.53k reads/sec
Bases Processed:       3154m    10.42m bases/sec

```







## Version: 

```
(bbtools) bash-4.4$ bbduk.sh -version
java -Djava.library.path=/home/jacqueline.zorz/software/miniconda3/envs/bbtools/bbtools/lib/jni/ -ea -Xmx41364m -Xms41364m -cp /home/jacqueline.zorz/software/miniconda3/envs/bbtools/bbtools/lib/current/ jgi.BBDukF -version
BBMap version 37.62
```

