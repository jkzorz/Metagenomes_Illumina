# Binning 

Binning with metabat, and using jgi_summarize_contigs to create depth file.
  
Example code with one sample at a time (as they finish mapping stage): 

**Metabat version: 2.12.1**

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

**CheckM version: v1.1.3 **
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

Example checkm coverage: 
```
checkm coverage -x fa . coverage_PurplePatch_04.tsv ../../mapping/mapping_hc_positive/reads_2A1-TheHole-C54-0-4_contig_2B1-PurplePatch-A54-0-4_sorted.bam 
```

Example checkm unbinned: 
```
checkm unbinned . /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/megahit/megahit_hc_positive/megahit_JZ-Condor-2B1-PurplePatch-A54-24-28_Li32230_S6/final.contigs.fa unbinned.fna unbinned_stats.tsv -x fa
```


## Gtdbtk
Example Gtdbtk script 

**Gtdbtk version: v1.6.0, release 207**

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

Gtdbtk script used for dereplicated bins. Updated to gtdbtk v2.0.0 and db release 207 

```
#!/bin/bash
###### Reserve computing resources ######
#SBATCH --mail-user=jacqueline.zorz@ucalgary.ca
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=500GB
#SBATCH --time=12:00:00
#SBATCH --partition=bigmem

###### Set environment variables ######
echo "Starting run at : 'date'"
source /home/jacqueline.zorz/software/miniconda3/etc/profile.d/conda.sh 
conda activate gtdbtk2

###### Run your script ######

gtdbtk classify_wf --genome_dir /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/binning/drep_out/dereplicated_genomes/ -x fa --out_dir /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/binning/drep_out/gtdbtk_bins_update207 --cpus 40


```




## dRep 

Test with bins from 3 samples. Need to put locations of each bin separately into a text file 

**dRep version: v3.2.2 **

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
 for i in *_depth/bin*.fa;do cp $i bins_all/; done
```

Created a consolidated checkM file so dRep doesn't have to run CheckM again: 

```
for i in *_depth/CheckM*.txt; do echo $i; awk '$1 ~/^bin/ {print $1 ".fa," $13 "," $14}' $i >> checkm_consolidated_drep.csv; done

#add headers
sed -i '1s/^/genome,completeness,contamination\n/' checkm_consolidated_drep.csv

```

dRep script: 

```
#!/bin/bash
###### Reserve computing resources ######
#SBATCH --mail-user=jacqueline.zorz@ucalgary.ca
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=180GB
#SBATCH --time=20:00:00
#SBATCH --partition=bigmem,cpu2019,cpu2021

###### Set environment variables ######
echo "Starting run at : 'date'"
source /home/jacqueline.zorz/software/miniconda3/etc/profile.d/conda.sh 
conda activate drep

cd /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/binning/


#run drep 
dRep dereplicate drep_out/ -p 25 -g /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/binning/bins_all/*.fa --genomeInfo /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/binning/checkm_consolidated_drep.csv
```

## CheckM coverage of good bins 

Using checkM coverage to estimate abundance of bins in community. (First moved all .bam and .bam.bai files to one folder)

```
#!/bin/bash
###### Reserve computing resources ######
#SBATCH --mail-user=jacqueline.zorz@ucalgary.ca
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=180GB
#SBATCH --time=96:00:00
#SBATCH --partition=cpu2019,cpu2021

###### Set environment variables ######
echo "Starting run at : 'date'"
source /home/jacqueline.zorz/software/miniconda3/etc/profile.d/conda.sh 
conda activate checkm

cd /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/binning/drep_out/

#script 

checkm coverage -t 20 -x fa /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/binning/drep_out/dereplicated_genomes/good_bins/ coverage.goodbins.tsv /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/mapping/mapping_bam/*.bam 

```

## DAS Tool for bin refinement

**Dastool version: 1.1.4**

May want to go back and try DAS Tool before dREP. 

First need fasta to contig file. This command didn't originally execute from dastool installation so I copied the code and put it in a new file directly. 
Need all bin files in same folder

```
sh Fasta_to_Contig2Bin.sh -i dastool_2B1-PurplePatch-A54-12-16/ -e fa > 2B1-PurplePatch-A54-12-16_contigs2bin.tsv
```

Next run DasTool command 
```
DAS_Tool -i 2B1-PurplePatch-A54-12-16_contigs2bin.tsv -c /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/binning/concoct/concoct_nospace_contigs/2B1-PurplePatch-A54-12-16_nospace_final.contigs.fa -o dastool_2B1-PurplePatch-A54-12-16/ --write_bins --write_bin_evals --write_unbinned -t 20
```




## Calculating coverage of MAGs

Issues: huge files, lots of unassembled and unbinned reads, many samples. Would like to avoid re-mapping reads to bins  

Potential solution: 
https://bitbucket.org/berkeleylab/metabat/issues/111/jgi_summarize_bam_contig_depths-to
You can take the weighted average of the set of contigs for each MAG. Maybe Import the depths.txt to a spreadsheet and sum (length*avg_coverage) / sum(length)

How JGI_summarize_bam_contigs_depth works: https://bitbucket.org/berkeleylab/metabat/issues/48/jgi_summarize_bam_contig_depths-coverage
1) The edges of a contig are generally excluded from the coverage counts up to a default of 75 bases or the average read length (--includeEdgeBases, --maxEdgeBases). This is because, generally mappers have a difficult time aligning a partial read to a contig when it would extend off edge and the coverage ramps up from 0 to the true coverage in this region 2) reads that map imperfectly are excluded when the %ID of the mapping drops below a threshold (--percentIdentity=97). MetaBAT is designed to resolve strain variation and mapping reads with low %ID indicate that the read actually came from a different strain/species. 3) clips/insertions/deletions/mismatches are excluded from the coverage count -- only the read bases that exactly match the reference are counted as coverage. This generally has a small effect, except in the case of long reads from PacBio and Nanopore.

I think the above responses will just give average coverage of the MAG in the sample. Marc suggested: **(sumproduct of depth and contig length)/total nucleotides in samples after qc**

This will give total number of nucleotides mapped to the MAG's contigs (sum(lengthxcontig depth)), as a fraction of all the nucleotides in the sample (total nucleotides in sample after QC) 

## Testing vamb

vamb: https://github.com/RasmussenLab/vamb

Need sample name in contig with separator "-" . Need to change both contig and depth file 
```
sed 's/>/>PurplePatch2428-/' /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/megahit/megahit_hc_positive/megahit_JZ-Condor-2B1-PurplePatch-A54-24-28_Li32230_S6/final.contigs.fa > /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/megahit/megahit_hc_positive/megahit_JZ-Condor-2B1-PurplePatch-A54-24-28_Li32230_S6/header_final.contigs.fa

sed 's/k141/PurplePatch2428-k141/g' /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/binning/PurplePatch_2428_depth/PurplePatch_2428_depth.txt > /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/binning/PurplePatch_2428_depth/header_PurplePatch_2428_depth.txt
```

vamb script: 

```
#!/bin/bash
###### Reserve computing resources ######
#SBATCH --mail-user=jacqueline.zorz@ucalgary.ca
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=180GB
#SBATCH --time=5:00:00
#SBATCH --partition=cpu2019,cpu2021,cpu2021-bf24,bigmem,cpu2019-bf05



###### Set environment variables ######
echo "Starting run at : 'date'"
source /home/jacqueline.zorz/software/miniconda3/etc/profile.d/conda.sh 
conda activate vamb

cd /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/binning/

#vamb --outdir path/to/outdir --fasta /path/to/catalogue.fna.gz --bamfiles /path/to/bam/*.bam -o C --minfasta 200000

vamb --outdir /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/binning/vamb_PurplePatch_2428 --fasta /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/megahit/megahit_hc_positive/megahit_JZ-Condor-2B1-PurplePatch-A54-24-28_Li32230_S6/header_final.contigs.fa --jgi /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/binning/PurplePatch_2428_depth/PurplePatch_2428_depth.txt -o - --minfasta 100000

```

## Testing CONCOCT

Need to remove spaces from headers of contig file? Having issues with contig headers. Removed spaces from headers in contig.fa file, but headers didn't match bam files then and I believe this is what was causing errors with the concoct_coverage_table.py script. According to CASI II, CONCOCT is supposed to be good for metagenomes with high strain diversity (test run with one sample managed to get low quality atribacteria bin). Will try using CONCOCT on metagenome data, but will need to pre-process contigs and sorted bam files so that they do not contain spaces. 

**Concoct version: 1.1.0**

### Pre-processing contig and sam files 

Removing everything after the first space from all contig files. Run from: /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/binning/concoct/concoct_nospace_contigs
```
for f in /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/megahit/megahit_hc_positive/megahit*/final.contigs.fa;

do 
	contig=$f
	contig_name=$(basename $(dirname $f))  #megahit_JZ-Condor-2A1-TheHole-C54-24-28_Li32297_S73
	contig_short=${contig_name:18} #2B1-TinyBubbles-C18-24-30_Li32239_S15
	contig_short2=${contig_short%_Li*} #2B1-TinyBubbles-C18-24-30

sed 's/\s.*$//' $f > ${contig_short2}_nospace_final.contigs.fa;

done

for f in /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/megahit/megahit_hc_negative/megahit*/final.contigs.fa;

do 
	contig=$f
	contig_name=$(basename $(dirname $f))  #megahit_JZ-Condor-2A1-TheHole-C54-24-28_Li32297_S73
	contig_short=${contig_name:18} #2B1-TinyBubbles-C18-24-30_Li32239_S15
	contig_short2=${contig_short%_Li*} #2B1-TinyBubbles-C18-24-30

sed 's/\s.*$//' $f > ${contig_short2}_nospace_final.contigs.fa;

done
```

**Bam headers**

Remove everything after first space, but keep "LN" at the end in bam headers. Need to first convert bam to sam, change headers, then copy headers to new bam file. Then need to re-index bam file.

```
conda activate samtools


cd /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/binning/concoct/concoct_nospace_bam

for BAM in /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/mapping/mapping_bam/*.bam;
do
	short_name=$(basename ${BAM})
	new_name="header_nospace_${short_name}"

	samtools view -H $BAM > test_header.sam
	sed "s/\sflag.*LN/\tLN/" test_header.sam > header_corrected.sam
	samtools reheader header_corrected.sam $BAM > $new_name
	samtools index $new_name
	rm test_header.sam
	rm header_corrected.sam
done

```


**Concoct for loop**

For loop to run concoct on each sample. Takes a really long time (~5 hours per sample). Split contig files into 3 groups and ran 3 different scripts to speed up the process. 
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
#SBATCH --partition=cpu2019,cpu2021

###### Set environment variables ######
echo "Starting run at : 'date'"
source /home/jacqueline.zorz/software/miniconda3/etc/profile.d/conda.sh 
conda activate concoct


cd /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/binning/concoct/

mkdir /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/binning/concoct/concoct_nospace_contigs/concoct_nospace_contigs_10K/

for i in /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/binning/concoct/concoct_nospace_contigs/*.fa;

do

#step 1
cut_up_fasta.py $i -c 10000 -o 0 --merge_last -b $(basename $i .fa)_10K.bed > /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/binning/concoct/concoct_nospace_contigs/concoct_nospace_contigs_10K/$(basename $i .fa)_10K.fa



#step 2
concoct_coverage_table.py $(basename $i .fa)_10K.bed /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/binning/concoct/concoct_nospace_bam/header_nospace_reads*$(basename $i _nospace_final.contigs.fa)_sorted.bam > $(basename $i _nospace_final.contigs.fa)_coverage_table.tsv

#step 3
concoct --composition_file /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/binning/concoct/concoct_nospace_contigs/concoct_nospace_contigs_10K/$(basename $i .fa)_10K.fa --coverage_file $(basename $i _nospace_final.contigs.fa)_coverage_table.tsv -b $(basename $i _nospace_final.contigs.fa)_concoct_output/ -t 25

#step 4
merge_cutup_clustering.py $(basename $i _nospace_final.contigs.fa)_concoct_output/$(basename $i _nospace_final.contigs.fa)_clustering_gt1000.csv > $(basename $i _nospace_final.contigs.fa)_concoct_output/$(basename $i _nospace_final.contigs.fa)_clustering_merged.csv

#step 5
#step 5
mkdir $(basename $i _nospace_final.contigs.fa)_concoct_output/fasta_bins

#step 6
extract_fasta_bins.py $i $(basename $i _nospace_final.contigs.fa)_concoct_output/$(basename $i _nospace_final.contigs.fa)_clustering_merged.csv --output_path $(basename $i _nospace_final.contigs.fa)_concoct_output/fasta_bins;

done
```

## Metabinner

Can use coverage file generated by CONCOCT. 

**Coverage profiles**

Seems you can't define an output name or location. Defaults to contig location. May be an issue in a for loop 
Also need full path to python script. 
```
python ~/software/miniconda3/envs/metabinner/bin/scripts/gen_kmer.py /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/binning/concoct/concoct_nospace_contigs/2B1-PurplePatch-A54-12-16_nospace_final.contigs.fa 1000 4
```

Filter out short contigs. Also can't define name or location of output... 
```
python ~/software/miniconda3/envs/metabinner/bin/scripts/Filter_tooshort.py /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/binning/concoct/concoct_nospace_contigs/2B1-PurplePatch-A54-12-16_nospace_final.contigs.fa 1000
```

**Running metabinner**

Need to be really particular with including all parameters and including full paths. 
-s size of metagenome (small, large, huge)
-p path to metabinner scripts (don't include final /) 

```
bash ~/software/miniconda3/envs/metabinner/bin/run_metabinner.sh -a 2B1-PurplePatch-A54-12-16_nospace_final.contigs_1000.fa -o /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/binning/metabinner/test -d /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/binning/concoct/2B1-PurplePatch-A54-12-16_coverage_table.tsv -k /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/binning/metabinner/kmer_4_f1000.csv -t 10 -s huge -p ~/software/miniconda3/envs/metabinner/bin
```
