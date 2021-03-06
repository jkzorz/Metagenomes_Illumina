## Assembly 
The assembly with MegaHit has generally been poor for individual samples (e.g. N50 <1000 bp, < 10% of reads assembled? - need to check this stat), and there are very few abundant species making high quality bins. The metagenomes from this study (https://www.nature.com/articles/s41396-019-0448-2#Sec2) were sequenced using a MiSeq sequencer. Only 5-15M reads obtained, but they managed to get 3 high quality Atribacteria bins. They also had a very high percentage of reads that assembled. They used metaspades and Maxbin for assembly and binning, respectively. 


According to Kraken. 5% of the raw reads from PurplePatch 24-28cm were Atribacteria (16% of assigned reads, 6 million reads), and 5.98% of the contigs from the same sample (12% of assigned contigs, 23K contigs) were assigned to Atribacteria. Of the unbinned contigs from the same sample 6% of contigs were Atribacteria (22k contigs, 13% of assigned contigs, 2 Atribacteria bins were actually made - 40 and 26% complete). Archaea actually don't seem to be assembling well as they are about 50% of raw reads, and only ~14% of contig sequences) 


## Mapping

Samples that are same depth (i.e. 0-4cm) have much higher mapping percentage than samples from same core but different depths (in some instances). 

E.g. 

The hole (2A-1) 0-4 cm vs The hole (2A-1) 12-16 cm: 5.5% of reads mapped

The hole (2A-1) 0-4 cm vs Purple haze (2A-2) 0-4 cm: 35% of reads mapped. 

The hole and purple haze were ~ 1 km apart. Therefore it seems that samples that are the same depth but 1 km apart are more similar than samples that are 8 cm different in depth but in the exact same location


## Binning

Depth files really help binning success.

e.g. 

Sample 2A-1 "The Hole" 12-16 cm
No depth file: 39 bins, 3 "good" bins (>90% completeness, <5% contamination)
Depth file from 37 samples mapping: 131 bins, 19 "good" bins (>90% completeness, <5% contamination)


Most contigs aren't binning still... 
Example from Purple Patch 0-4 cm 

**Unbinned:** 
1749417 contigs

**Total:**
1852709 contigs 

Only ~6% of contigs (103,292) binned. 


### Storage space

```
(base) bash-4.4$ cd /work/ebg_lab/eb/sodalakes/
(base) bash-4.4$ du -sh
1.9T    .

(base) bash-4.4$ cd ../autofermentation/
(base) bash-4.4$ du -sh
1.8T    .

(base) bash-4.4$ cd ../../gm/gapp/jzorz/
(base) bash-4.4$ du -sh
7.1T    .
```
