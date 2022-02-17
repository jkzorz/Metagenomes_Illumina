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
