# Pairwise comparison with Mash

Downloaded mash from: https://github.com/marbl/Mash/releases
No dependencies, so just downloaded and transferred over to home drive on server. Unzipped using ``` tar -xf```


```
#create sketch 
~/mash-Linux64-v2.3/mash sketch -o test /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/megahit/megahit_hc_negative/megahit_JZ-Condor-2AT-175NW-E46-0-4_Li32302_S78/final.contigs.fa /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/megahit/megahit_hc_negative/megahit_JZ-Condor-2AT-175NW-E46-12-16_Li32320_S96/final.contigs.fa

#find distances between contigs 
 ~/mash-Linux64-v2.3/mash dist test.msh /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/megahit/megahit_hc_negative/megahit_JZ-Condor-2AT-700NW-B7-0-4_Li32307_S83/final.contigs.fa > out_test.txt

```

Headers of tab-delimited output are: Reference-ID, Query-ID, Mash-distance, P-value, and Matching-hashes

### Contigs

First ran Mash on contig files. Needed to gather locations of all megahit contig files: 

```
#list contig file locations 
for i in /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/megahit/megahit_hc_positive/megahit_JZ-Condor*/final.contigs.fa; do printf $i' '; done > contigs_list.txt
for i in /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/megahit/megahit_hc_negative/megahit_JZ-Condor*/final.contigs.fa; do printf $i' '; done >> contigs_list.txt

```

Then ran mash sketch and dist by copying and pasting locations from contigs_list.txt to end of following commands: 

```
#sketch
~/mash-Linux64-v2.3/mash sketch -o contigs

#dist
~/mash-Linux64-v2.3/mash dist contigs.msh
```




