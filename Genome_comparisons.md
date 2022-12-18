# Comparing to MAGs from NCBI



**unzip files**
```
/work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/compare

for i in */*.fna.gz; do gunzip $i; done
```


**change file extensions** 
File extenstions need to be the same for checkm and gtdbtk

```
for i in */*.fna; do sample=$(basename $i .fna).fa; phylum=$(dirname $i); mv $i $phylum/$sample; done
```

## CheckM2 on all genomes
```
#!/bin/bash
###### Reserve computing resources ######
#SBATCH --mail-user=jacqueline.zorz@ucalgary.ca
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=100GB
#SBATCH --time=24:00:00
#SBATCH --partition=cpu2019,cpu2021,cpu2021-bf24

###### Set environment variables ######
echo "Starting run at : 'date'"
source /home/jacqueline.zorz/software/miniconda3/etc/profile.d/conda.sh 
conda activate checkm2

for i in /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/compare/*/

do ~/checkm2/bin/checkm2 predict --threads 30 --input $i --output-directory $i/checkm2_compare_output -x fa

done
```



## GTDBTK on all genomes
```
#!/bin/bash
###### Reserve computing resources ######
#SBATCH --mail-user=jacqueline.zorz@ucalgary.ca
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=350GB
#SBATCH --time=16:00:00
#SBATCH --partition=bigmem

###### Set environment variables ######
echo "Starting run at : 'date'"
source /home/jacqueline.zorz/software/miniconda3/etc/profile.d/conda.sh 
conda activate gtdbtk2

###### Run your script ######

cd /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/compare

for i in /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/compare/*/checkm2*

do 
sample=$(dirname $i)
gtdbtk classify_wf --genome_dir $sample -x fa --out_dir $sample/gtdbtk_bins/ --cpus 20
done
```

## DRAM on all genomes
```
#!/bin/bash
###### Reserve computing resources ######
#SBATCH --mail-user=jacqueline.zorz@ucalgary.ca
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=800GB
#SBATCH --time=24:00:00
#SBATCH --partition=bigmem

#set environmental parameters
source /home/jacqueline.zorz/software/miniconda3/etc/profile.d/conda.sh 
conda activate DRAM

cd /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/compare

####command

DRAM.py annotate -i '/work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/compare/*/*.fa' -o /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/compare/dram --threads 20

```

Distill results: 
```
DRAM.py distill -i all_compare_annotations.tsv -o genome_summaries --trna_path all_compare_trnas.tsv --rrna_path all_compare_rrnas.tsv
```

## Barrnap 

Use Barrnap to retrieve rrna sequences:

```
conda activate barrnap

cd /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/compare

for i in /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/compare/*/*.fa; do sample=rrna_$(basename $i); phy=$(basename $(dirname $i)); mkdir /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/compare/$phy/barrnap_$phy; barrnap $i --outseq $phy/barrnap_$phy/$sample; done

```



