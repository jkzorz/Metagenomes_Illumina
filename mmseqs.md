# MMSEQs

Use mmseqs to cluster MAG proteins to look for instances of HGT. 

```
#first needed to create concatenated faa file with MAG names in gene headers
cd /work/ebg_lab/gm/gapp/jzorz/Metagenomes_Illumina/dereplicated_genomes98/drep98_proteins/

#add MAG name to headers
for i in genes_protein/*.faa; do mag=$(basename $i .faa); sed "s/>/>${mag}:/g" $i > genes_protein_header/${mag}_header.faa;done

#condatenate files
cat genes_protein_header/*.faa > all_mag_proteins_drep98_header.faa

```

Run mmseqs to cluster all proteins with >80% identity over 80% of the protein
```
conda activate mmseqs

#run mmseqs with easy-cluster
#coverage mode 0: coverage of query and target
mmseqs easy-cluster all_mag_proteins_drep98_header.faa all_mag_proteins_drep98_cov80_pc80 tmp --min-seq-id 0.8 -c 0.8 --cov-mode 0


```


