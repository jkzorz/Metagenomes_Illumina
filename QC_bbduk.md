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

Subsampled reads to 10% of original using reformat.sh from bbtools
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

Seems that for both megahit and metaspades assemblies, the subsampling approach did not help (at least at 10% of reads). 

### Normalization

Used bbnorm.sh to normalize read depth...

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
conda activate bbtools

cd /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/bbduk/

#bbnorm 
#e.g. bbnorm.sh in=reads.fq out=normalized.fq target=100 min=5
bbnorm.sh -Xmx90g in=cat_qc/JZ-Condor-2B1-PurplePatch-A54-0-4_Li32225_S1_R1_QC.fastq in2=cat_qc/JZ-Condor-2B1-PurplePatch-A54-0-4_Li32225_S1_R2_QC.fastq out=JZ-Condor-2B1-PurplePatch-A54-0-4_Li32225_S1_R1_QC_norm100.fastq out2=JZ-Condor-2B1-PurplePatch-A54-0-4_Li32225_S1_R2_QC_norm100.fastq target=100 min=5
```

Output: 
```
Starting run at : 'date'
java -ea -Xmx90g -Xms90g -cp /home/jacqueline.zorz/software/miniconda3/envs/bbtools/bbtools/lib/current/ jgi.KmerNormalize bits=32 -Xmx90g in=cat_qc/JZ-Condor-2B1-PurplePatch-A54-0-4_Li32225_S1_R1_QC.fastq in2=cat_qc/JZ-Condor-2B1-PurplePatch-A54-0-4_Li32225_S1_R2_QC.fastq out=JZ-Condor-2B1-PurplePatch-A54-0-4_Li32225_S1_R1_QC_norm100.fastq out2=JZ-Condor-2B1-PurplePatch-A54-0-4_Li32225_S1_R2_QC_norm100.fastq target=100 min=5
Executing jgi.KmerNormalize [bits=32, -Xmx90g, in=cat_qc/JZ-Condor-2B1-PurplePatch-A54-0-4_Li32225_S1_R1_QC.fastq, in2=cat_qc/JZ-Condor-2B1-PurplePatch-A54-0-4_Li32225_S1_R2_QC.fastq, out=JZ-Condor-2B1-PurplePatch-A54-0-4_Li32225_S1_R1_QC_norm100.fastq, out2=JZ-Condor-2B1-PurplePatch-A54-0-4_Li32225_S1_R2_QC_norm100.fastq, target=100, min=5]

BBNorm version 37.62

   ***********   Pass 1   **********   


Settings:
threads:          	48
k:                	31
deterministic:    	true
toss error reads: 	false
passes:           	1
bits per cell:    	16
cells:            	33.77B
hashes:           	3
base min quality: 	5
kmer min prob:    	0.5

target depth:     	400
min depth:        	3
max depth:        	500
min good kmers:   	15
depth percentile: 	64.8
ignore dupe kmers:	true
fix spikes:       	false

Made hash table:  	hashes = 3   	 mem = 62.84 GB   	cells = 33.73B   	used = 62.439%
Warning:  This table is somewhat full, which may reduce accuracy.  Ideal load is under 60% used.
For better accuracy, use the 'prefilter' flag; run on a node with more memory; quality-trim or error-correct reads; or increase the values of the minprob flag to reduce spurious kmers.  In practice you should still get good normalization results even with loads over 90%, but the histogram and statistics will be off.

Estimated unique kmers:     	11011026132

Table creation time:		450.588 seconds.
Started output threads.
Table read time: 		629.586 seconds.   	50103.21 kb/sec
Total reads in:  		225676410	60.554% Kept
Total bases in:  		31544289481	60.597% Kept
Error reads in:  		110650129	49.030%
Error pairs in:  		63714725 	56.466%
Error type 1:    		93542624 	41.450%
Error type 2:    		16690362 	7.396%
Error type 3:    		3120408  	1.383%
Total kmers counted:          	24773733315
Total unique kmer count:      	11989703856
Includes forward kmers only.
The unique kmer estimate can be more accurate than the unique count, if the tables are very full.
The most accurate value is the greater of the two.

Percent unique:               	48.40%
Depth average:                	2.07	(unique kmers)
Depth median:                 	1	(unique kmers)
Depth standard deviation:     	4.22	(unique kmers)
Corrected depth average:      	1.68	

Depth average:                	10.68	(all kmers)
Depth median:                 	3	(all kmers)
Depth standard deviation:     	72.79	(all kmers)

Approx. read depth median:    	3.82

   ***********   Pass 2   **********   


Settings:
threads:          	48
k:                	31
deterministic:    	true
toss error reads: 	false
passes:           	1
bits per cell:    	16
cells:            	33.77B
hashes:           	3
base min quality: 	5
kmer min prob:    	0.5

target depth:     	100
min depth:        	5
max depth:        	100
min good kmers:   	15
depth percentile: 	54.0
ignore dupe kmers:	true
fix spikes:       	false

Made hash table:  	hashes = 3   	 mem = 62.84 GB   	cells = 33.73B   	used = 30.826%

Estimated unique kmers:     	4144231361

Table creation time:		301.637 seconds.
Started output threads.
Table read time: 		318.300 seconds.   	60053.32 kb/sec
Total reads in:  		136655458	61.691% Kept
Total bases in:  		19114983525	61.491% Kept
Error reads in:  		39561734 	28.950%
Error pairs in:  		26946040 	39.436%
Error type 1:    		22102320 	16.174%
Error type 2:    		17123870 	12.531%
Error type 3:    		1718251  	1.257%
Total kmers counted:          	15015168646
Total unique kmer count:      	4919973474
Includes forward kmers only.
The unique kmer estimate can be more accurate than the unique count, if the tables are very full.
The most accurate value is the greater of the two.

Percent unique:               	32.77%
Depth average:                	3.05	(unique kmers)
Depth median:                 	2	(unique kmers)
Depth standard deviation:     	5.69	(unique kmers)
Corrected depth average:      	2.89	

Depth average:                	13.66	(all kmers)
Depth median:                 	5	(all kmers)
Depth standard deviation:     	45.31	(all kmers)

Approx. read depth median:    	6.37

Removing temp files.

Total time:      		1700.925 seconds.   	29783.36 kb/sec
```


## Version: 

```
(bbtools) bash-4.4$ bbduk.sh -version
java -Djava.library.path=/home/jacqueline.zorz/software/miniconda3/envs/bbtools/bbtools/lib/jni/ -ea -Xmx41364m -Xms41364m -cp /home/jacqueline.zorz/software/miniconda3/envs/bbtools/bbtools/lib/current/ jgi.BBDukF -version
BBMap version 37.62
```

