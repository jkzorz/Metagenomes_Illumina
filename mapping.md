# Mapping reads to contigs with bbmap 

For now, splitting samples into hydrocarbon positive (2A1, 2A2, 2B1) and hydrocarbon negative (BG15, 2AT). Mapping reads from all hydrocarbon positive samples to all hydrocarbon positive contigs, and reads from all hydrocarbon negative samples to all hydrocarbon negative contigs.  

There are 19 hydrocarbon positive samples, so 19 reads x 19 contigs = 361 resulting bam files... 

Chose minid=0.9 to speed up process. Added 'if' statements so that command could be run multiple times without overwriting existing files. 

**bbmap (bbtools) version: 37.62**

```
#!/bin/bash
###### Reserve computing resources ######
#SBATCH --mail-user=jacqueline.zorz@ucalgary.ca
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=300GB
#SBATCH --time=24:00:00
#SBATCH --partition=bigmem


###### Set environment variables ######
echo "Starting run at : 'date'"
source /home/jacqueline.zorz/software/miniconda3/etc/profile.d/conda.sh 
conda activate bbtools

cd /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/mapping/mapping_hc_positive


for f in /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/megahit/megahit_hc_positive/megahit*/final.contigs.fa;  

do
	contig=$f
	contig_name=$(basename $(dirname $f))  #megahit_JZ-Condor-2A1-TheHole-C54-24-28_Li32297_S73
	contig_short=${contig_name:18} #2B1-TinyBubbles-C18-24-30_Li32239_S15
	contig_short2=${contig_short%_Li*} #2B1-TinyBubbles-C18-24-30

	
	for x in /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/bbduk/cat_qc/*R1_QC.fastq;
	do 
		R1=$x
		R2=$(dirname $x)/$(basename $x R1_QC.fastq)R2_QC.fastq
		reads=$(basename $x _R1_QC.fastq) #eg JZ-Condor-2B1-TinyBubbles-C18-0-4_Li32293_S69
		read_short=${reads%_Li*}
		read_short2=${read_short:10}
		read_names="reads_${read_short2}"
		sam="${read_names}_contig_${contig_short2}.sam"
		bam="${sam::-3}bam"
		sorted_bam="${bam::-4}_sorted.bam"
		covstats="${sam::-4}_covstats.txt"
		scafstats="${sam::-4}_scafstats.txt"

		if [ -e "$sam" ]; then
			continue    # the output file already exists, so skip re-creating it
		fi

		if [ -e "$covstats" ]; then
			continue    # the output file already exists, so skip re-creating it
		fi

		#step1 - mapping and sam file creation 
		bbmap.sh ref=$contig in=$R1 in2=$R2 outm=$sam covstats=$covstats scafstats=$scafstats threads=40 minid=0.9


		##Step2 sam to bam
		samtools view -b $sam | samtools sort -o $sorted_bam


		###Step3 indexing sorted bam file
		samtools index $sorted_bam
		
		##Step4 remove sam file 
		rm $sam;
		
	done
done 

```

## Create a file with the total number of reads mapped to each sample's contigs

Do for negative and positive hc directories. 
```
cd /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/mapping/mapping_hc_negative
for file in *covstats.txt; do echo "$file: $(awk 'NR>1 {sum += $5} END {print sum}' "$file")"; done > ../hc_negative_covstats_sums.txt

cd /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/mapping/mapping_hc_positive
for file in *covstats.txt; do echo "$file: $(awk 'NR>1 {sum += $5} END {print sum}' "$file")"; done > ../hc_positive_covstats_sums.txt
```
