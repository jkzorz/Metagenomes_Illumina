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
~/mash-Linux64-v2.3/mash dist contigs.msh ... > contigs.mash
```

Clean up results file to remove paths and file suffixes

```
sed 's/\/work\/ebg_lab\/gm\/gapp\/jzorz\/Metagenomes_Illumina\/megahit\/megahit_hc_positive\/megahit_JZ-Condor-//g' contigs_mash.txt > contigs_mash2.txt
sed 's/\/work\/ebg_lab\/gm\/gapp\/jzorz\/Metagenomes_Illumina\/megahit\/megahit_hc_negative\/megahit_JZ-Condor-//g' contigs_mash2.txt > contigs_mash3.txt

sed 's/\/final\.contigs\.fa//g' contigs_mash3.txt > contigs_mash4.txt
```

### NMDS in R 

Brought contigs_mash4.txt to local computer. Deleted p-value and Matching-hashes columns and saved as csv. Ran the following commands in R to make an NMDS plot of the data 

```
library(tidyverse)
library(vegan)
setwd("~/University of Calgary/PostDoc/Metagenomes")
set.seed(123)

mash = read.csv("contigs_mash4.csv", header = FALSE)

#turn back into wide matrix
mash2 = mash %>% pivot_wider(values_from = V3, names_from = V2)

#turn into distance matrix object
mash3 = as.matrix(mash2[,-1])
dish = as.dist(mash3)

#nmds of dish matrix 
nmds = metaMDS(dish)

#extract scores
data.scores = as.data.frame(scores(nmds))

#add columns
data.scores$sample = row.names(data.scores)
data.scores2 =  data.scores %>% separate(sample, c("Site", "Subsite", "Core", "Depth1", "Depth2"), sep = "-")
data.scores2$Depth2 = gsub(pattern = "_.*", replacement = "", data.scores2$Depth2)

#plot sites
gg = ggplot(data.scores2, aes(x = NMDS1, y = NMDS2)) + geom_point(aes(colour = Site, size = as.numeric(Depth1))) + theme(panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey40"), legend.key = element_blank()) + labs(size = "Depth (cmbsf)") + scale_colour_manual(values = c("#ED315D", "#F78C6B", "#049F76", "#FCC088", "#83D483")) + scale_radius(range = c(2,6), breaks = c(0,12,24))

#plot subsites
gg = ggplot(data.scores2, aes(x = NMDS1, y = NMDS2)) + geom_point(aes(colour = Subsite, size = as.numeric(Depth1))) + theme(panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey40"), legend.key = element_blank()) + labs(size = "Depth (cmbsf)") + scale_radius(range = c(2.5,6), breaks = c(0,12,24)) + scale_colour_manual(values = c("#2B193D", "#585285", "#4E61A6", "#8993BD", "#ced0db", "#C5979D", "#293E3C", "#547856", "#71A352", "#bfd1b4"))+ guides(colour = guide_legend(override.aes = list(size=4)))

```


