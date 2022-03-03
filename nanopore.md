# Incorporating Nanopore long reads into Illumina metagenomes

Using long Nanopore reads to try to improve assembly of Illumina metagenomes. 

Nanopore does it's own quality control for the reads (failed vs passed folders). Ignore the failed reads and concatenate all passed reads into one large file. Use porechop to remove adapters.

```
#start interactive session to run commands on arc server 
salloc --mem=50G -c 20 -N 1 -n 1  -t 05:00:00

#concatenate all fastq files that passed QC
cat fastq_pass/* > seqs.fastq

#activate porechop environment
conda activate porechop

#run porechop to trim adapters
porechop -i seqs.fastq -o seqs_trimmed.fastq -t 8

```


## Assembly usingn nanopore and Illumina reads

Unicycler? (https://github.com/rrwick/Unicycler#installation)







