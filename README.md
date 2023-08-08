# Variable genomic trajectories of speciation across animals #


### Analytical approach ###
For each pair of genomes (congeneric species), microsynteny and macrosynteny are both analysed.

The pipeline processor [run_synteny_analysis.py](https://github.com/PalMuc/speciation_synteny/blob/main/02-processing_scripts/run_synteny_analysis.py) is coded in Python, and run simply as:

`../run_synteny_analysis.py -i ../species_pair_list.tab`

For each species pair, for example the tuna, this begins with the scaffolds, proteins, and GFF downloaded from [NCBI](https://www.ncbi.nlm.nih.gov/genome/?term=txid8240[orgn]):

```
GCF_910596095.1_fThuMac1.1_genomic.fna.gz
GCF_910596095.1_fThuMac1.1_genomic.gff.gz
GCF_910596095.1_fThuMac1.1_protein.faa.gz
GCF_914725855.1_fThuAlb1.1_genomic.fna.gz
GCF_914725855.1_fThuAlb1.1_genomic.gff.gz
GCF_914725855.1_fThuAlb1.1_protein.faa.gz
```

and this generates the following files for each species:

* [get_genbank_longest_isoforms.py](https://github.com/PalMuc/speciation_synteny/blob/main/02-processing_scripts/get_genbank_longest_isoforms.py) filtered proteins with isoforms removed `.x.faa`, like: `GCF_910596095.1_fThuMac1.1_protein.x.faa` and `GCF_914725855.1_fThuAlb1.1_protein.x.faa`
* [get_genbank_longest_isoforms.py](https://github.com/PalMuc/speciation_synteny/blob/main/02-processing_scripts/get_genbank_longest_isoforms.py) filtered GFFs corresponding to the proteins `.x.gff`, like: `GCF_910596095.1_fThuMac1.1_genomic.x.gff` , `GCF_914725855.1_fThuAlb1.1_genomic.x.gff`
* [DIAMOND](https://github.com/bbuchfink/diamond) results `fThuAlb1_vs_fThuMac1.blastp.tab` and `fThuAlb1_vs_fThuMac1.renamed.blastp.tab`
* [scaffold_synteny.py](https://github.com/wrf/genomeGTFtools/blob/master/scaffold_synteny.py) results `fThuAlb1_vs_fThuMac1.scaffold_synteny.tab` and `fThuAlb1_vs_fThuMac1.scaffold_synteny.pdf`
* [microsynteny.py](https://github.com/wrf/genomeGTFtools/blob/master/microsynteny.py) results `fThuAlb1_vs_fThuMac1.microsynteny.tab` and `fThuAlb1_vs_fThuMac1.microsynteny.pdf`
* [fastarenamer.py](https://bitbucket.org/wrf/sequences/src/master/fastarenamer.py) renamed versions of proteins for clustering `.x.n.faa`, like: `GCF_910596095.1_fThuMac1.1_protein.x.n.faa` , `GCF_914725855.1_fThuAlb1.1_protein.x.n.faa`
* [makehomologs.py](https://bitbucket.org/wrf/sequences/src/master/makehomologs.py) clustering outputs `fasta_clusters.H.thunnus_clusters_v1.tab` `clusters_thunnus_clusters_v1.tar.gz` and log `thunnus_clusters_v1.2023-08-02-010624.mh.log`
* [alignment_conserved_site_to_dots.py](https://bitbucket.org/wrf/sequences/src/master/alignment_conserved_site_to_dots.py) accumulated tabular output `fThuAlb1_vs_fThuMac1.homologs_identity.tab`

Subsequent processing occurs using [several R scripts](https://github.com/PalMuc/speciation_synteny/tree/main/03-graphic_scripts), for analysis and plotting.


### Full citation ###


<a rel="license" href="http://creativecommons.org/licenses/by-sa/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-sa/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-sa/4.0/">Creative Commons Attribution-ShareAlike 4.0 International License</a>.
