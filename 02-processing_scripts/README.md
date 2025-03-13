## processing scripts ##
Data were downloaded from NCBI in 2023. Because multiple assembly versions exist for some species, refer to the data table [species_pair_list.tab](https://github.com/PalMuc/congeneric_synteny/blob/main/02-processing_scripts/species_pair_list.tab) for the precise `GCF_` accessions used. **NOTE: some of these accessions may not longer be the primary assembly for this species.**

Some additional code is cloned from [WRF's genomeGTFtools repo](https://github.com/wrf/genomeGTFtools) and [WRF's sequence processing repo](https://bitbucket.org/wrf/sequences/). The versions used in the paper are included in this folder.

The pipeline processor [run_synteny_analysis.py](https://github.com/PalMuc/speciation_synteny/blob/main/02-processing_scripts/run_synteny_analysis.py) is coded in Python, and using the above table [species_pair_list.tab](https://github.com/PalMuc/congeneric_synteny/blob/main/02-processing_scripts/species_pair_list.tab) that indicates the input files, is run simply as:

`run_synteny_analysis.py -i species_pair_list.tab`

For each pair of genomes (congeneric species), microsynteny and macrosynteny are both analysed.
