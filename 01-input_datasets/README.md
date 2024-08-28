# datasets for synteny analysis #

Using additional code at [WRF's genomeGTFtools](https://github.com/wrf/genomeGTFtools), [WRF's sequence processing](https://bitbucket.org/wrf/sequences/)

## *Acropora* ##
Using [*Acropora hyacinthus*](https://marinegenomics.oist.jp/ahya/viewer/download?project_id=91) and [*Acropora millepora*](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_013753865.1/) by [Shinzato et al 2021](https://doi.org/10.1093/molbev/msaa216)

```
~/gffread-0.12.7.Linux_x86_64/gffread -g ahya.fasta -w ahya.nucl.fa ahya.gff
grep -E -f  ../../excluder_list.txt ahya.nucl.fa | cut -d ' ' -f 1 > ahya.nucl.x.names
excludeAinB.py ahya.nucl.x.names ahya.nucl.fa > ahya.nucl.x.fa
# Found 4067 terms to exclude
# Reading sequences from ahya.nucl.fa
# Excluded 4067 names out of 27215 sequences from A in B
# Wrote 23148 sequences
~/minimap2-2.23_x64-linux/minimap2 -a -x splice --secondary=no GCA_020536085.1_Ahyacinthus.chrsV1_genomic.fna ahya.nucl.x.fa | ~/samtools-1.14/samtools sort - -o ahya.nucl.vs_chrsV1.bam
~/git/pinfish/spliced_bam2gff/spliced_bam2gff -M -s ahya.nucl.vs_chrsV1.bam > ahya.nucl.vs_chrsV1.gtf
~/git/genomeGTFtools/misc/stringtie_gtf_to_gff3.py ahya.nucl.vs_chrsV1.gtf -r | sed s/exon/CDS/g > ahya.nucl.vs_chrsV1.gff
~/gffread-0.12.7.Linux_x86_64/gffread -g GCA_020536085.1_Ahyacinthus.chrsV1_genomic.fna -x ahya.vs_chrsV1.nucl.fa ahya.nucl.vs_chrsV1.gff
prottrans.py -n -r ahya.vs_chrsV1.nucl.fa | remove_identical_seqs.py - > ahya.vs_chrsV1.prot.fa
# Read 24624 sequences, wrote 24149
# Total counted bases: 8818454
grep -E -f  ../../excluder_list.txt ahya.vs_chrsV1.prot.fa | cut -d ' ' -f 1 > ahya.vs_chrsV1.prot.x.names
excludeAinB.py ahya.vs_chrsV1.prot.x.names ahya.vs_chrsV1.prot.fa > ahya.vs_chrsV1.prot.x.fa
# Found 2909 terms to exclude
# Reading sequences from ahya.vs_chrsV1.prot.fa
# Excluded 2966 names out of 26944 sequences from A in B
# Wrote 23978 sequences
getgtfgenes.py ahya.vs_chrsV1.prot.x.fa ahya.nucl.vs_chrsV1.gff -f > ahya.nucl.vs_chrsV1.x.gff
# Reading queries from ahya.vs_chrsV1.prot.x.fa  12:01:07 2023
# Read 184190 lines for 22664 entries  12:01:07 2023
# Searching features from ahya.nucl.vs_chrsV1.gff  12:01:07 2023
# Counted 277897 lines for 29031 transcripts  12:01:08 2023
# Wrote 211019 lines, with 24704 top-level features

~/project/speciation-synteny/get_genbank_longest_isoforms.py -g GCF_013753865.1_Amil_v2.1_genomic.gff.gz -p GCF_013753865.1_Amil_v2.1_protein.faa.gz -f GCF_013753865.1_Amil_v2.1_genomic.x.gff -o GCF_013753865.1_Amil_v2.1_protein.x.faa
# Parsing protein fasta GCF_013753865.1_Amil_v2.1_protein.faa.gz as gzipped   15:32:15 2023
# Read 41860 proteins   15:32:16 2023
# Calculating protein lengths   15:32:16 2023
# Got 41860 protein lengths   15:32:16 2023
# Parsing protein fasta GCF_013753865.1_Amil_v2.1_genomic.gff.gz as gzipped   15:32:16 2023
# Found 30136 genes for 41860 proteins   15:32:21 2023
# Names parsed as gene-LOC122956530  XP_044172110.1
# Wrote 30136 proteins as CDS   15:32:21 2023
# Wrote 30136 proteins   15:32:21 2023

~/diamond-v2.0.13/diamond-v2.0.13 makedb --in GCF_013753865.1_Amil_v2.1_protein.x.faa -d GCF_013753865.1_Amil_v2.1_protein.x.faa
~/diamond-v2.0.13/diamond-v2.0.13 blastp -q ahya.vs_chrsV1.prot.x.fa -d GCF_013753865.1_Amil_v2.1_protein.x.faa -o AhyaV1_vs_Amilv2.blastp.tab
#Reported 170796 pairwise alignments, 170796 HSPs.
#20912 queries aligned.
~/git/genomeGTFtools/scaffold_synteny.py -b AhyaV1_vs_Amilv2.blastp.tab -q ahya.nucl.vs_chrsV1.x.gff -d GCF_013753865.1_Amil_v2.1_genomic.gff.gz -f GCA_020536085.1_Ahyacinthus.chrsV1_genomic.fna.gz -F GCF_013753865.1_Amil_v2.1_genomic.fna.gz --ignore-gene-features --genbank-gff > AhyaV1_vs_Amilv2.scaffold_synteny.tab
# Parsing genomic contigs GCA_020536085.1_Ahyacinthus.chrsV1_genomic.fna.gz as gzipped  11:57:14 2023
# Found 907 contigs  11:57:18 2023
# Sorting contigs by length, keeping up to 99999Mbp
# Kept 672 contigs, for 445017284 bases, last contig was 9979bp long  11:57:18 2023
# Parsing genomic contigs GCF_013753865.1_Amil_v2.1_genomic.fna.gz as gzipped  11:57:18 2023
# Found 854 contigs  11:57:24 2023
# Sorting contigs by length, keeping up to 99999Mbp
# Kept 672 contigs, for 474498957 bases, last contig was 9974bp long  11:57:24 2023
# Parsing loci from ahya.nucl.vs_chrsV1.x.gff  11:57:24 2023
# Found 24149 genes  11:57:24 2023
# GFF names parsed as ahya_s0119.g34.t1b from ahya_s0119.g34.t1b
# 24149 RNA to protein IDs parsed as ahya_s0119.g34.t1b from None
# Parsing loci from GCF_013753865.1_Amil_v2.1_genomic.gff.gz as gzipped  11:57:24 2023
# Found 87852 genes  11:57:28 2023
# GFF names parsed as rna-XM_044316175.1 from rna-XM_044316175.1
# 41860 RNA to protein IDs parsed as rna-XM_044316175.1 from XP_044172110.1
# Parsing tabular blast output AhyaV1_vs_Amilv2.blastp.tab  11:57:28 2023
# Found blast hits for 20891 query sequences, removed 1768 hits by evalue  11:57:29 2023
# Removed 0 queries and 0 subjects with 250 or more hits
# Blast names parsed as ahya_s0119.g34.t1b from ahya_s0119.g34.t1b, and XP_029212286.2 from XP_029212286.2
# Kept 20891 blast hits
# Determining match positions  11:57:29 2023
# Wrote match positions for 20865 genes
# Wrote target positions for 16557 genes
Rscript ~/git/genomeGTFtools/synteny_2d_plot.R AhyaV1_vs_Amilv2.scaffold_synteny.tab A.hyacinthus A.millepora 240 

~/git/genomeGTFtools/microsynteny.py -b AhyaV1_vs_Amilv2.blastp.tab -q ahya.nucl.vs_chrsV1.x.gff -d GCF_013753865.1_Amil_v2.1_genomic.gff.gz --genbank-gff -z 100000 -T > AhyaV1_vs_Amilv2.microsynteny.tab
Rscript ~/git/genomeGTFtools/synteny_block_length_plot.R AhyaV1_vs_Amilv2.microsynteny.tab 0.95

fastarenamer.py -p Ahya -j '|' ahya.vs_chrsV1.prot.x.fa > ahya.vs_chrsV1.prot.x.n.fa
fastarenamer.py -p Amil -j '|' GCF_013753865.1_Amil_v2.1_protein.x.faa > GCF_013753865.1_Amil_v2.1_protein.x.n.faa
~/diamond-v2.0.13/diamond-v2.0.13 makedb --in GCF_013753865.1_Amil_v2.1_protein.x.n.faa -d GCF_013753865.1_Amil_v2.1_protein.x.n.faa
~/diamond-v2.0.13/diamond-v2.0.13 blastp -q ahya.vs_chrsV1.prot.x.n.fa -d GCF_013753865.1_Amil_v2.1_protein.x.n.faa -o AhyaV1_vs_Amilv2.blastp.n.tab
#Reported 170796 pairwise alignments, 170796 HSPs.
#20912 queries aligned.
makehomologs.py -i AhyaV1_vs_Amilv2.blastp.n.tab -f ahya.vs_chrsV1.prot.x.n.fa GCF_013753865.1_Amil_v2.1_protein.x.n.faa -p 234 -o acropora_clusters_v1 -z 2 -s 2 -M 10 -H 9 -c

for FILE in clusters_acropora_clusters_v1/*.fasta ; do BASE="${FILE%.fasta}" ; mafft $FILE > $BASE.aln ; done

for FILE in clusters_acropora_clusters_v1/homologs_acropora_clusters_v1_*.aln ; do alignment_conserved_site_to_dots.py -t -n -a $FILE >> AhyaV1_vs_Amilv2.homologs_identity.tab ; done
```


## Placozoa ##
Using [*Hoilungia hongkongensis*](https://bitbucket.org/molpalmuc/hoilungia-genome/src/master/) and [*Trichoplax adhaerens* redo](https://bitbucket.org/wrf/genome-reannotations/src/master/jbrowse-tracks/trichoplax/) by [Eitel et al 2018](https://doi.org/10.1371/journal.pbio.2005359)

```
for FILE in 6554-orthogs_proteins_alignments/*.fasta ; do alignment_conserved_site_to_dots.py -n -t -a $FILE >> hoilungia_v_trichoplax.homologs_identity.tab ; done

```

## Hexactinellida ##
Using [*Aphrocallistes vastus*](https://github.com/PalMuc/Aphrocallistes_vastus_genome) by [Francis et al 2023](https://doi.org/10.1098/rsos.230423) and [*Oopsacas minuta*](https://www.ncbi.nlm.nih.gov/datasets/taxonomy/111878/) by [Santini et al 2023](https://doi.org/10.1186/s12915-023-01619-w)

```
zcat Avas.v1.29_annotations.prot.vs_oopsacas_gb.tab.gz | sed s/LOD99/"Omin|LOD99"/g | sed s/Avas/"Avas|Avas"/g > Avas.v1.29_annotations.prot.vs_oopsacas_gb.tab

cat ~/git/Aphrocallistes_vastus_genome/synteny/vs_oopsacas/avas_omin_proteins_combined.fasta | sed s/O_minuta/Omin/g | sed s/Avas/"Avas|Avas"/g > avas_omin_proteins_combined.fasta

makehomologs.py -i Avas.v1.29_annotations.prot.vs_oopsacas_gb.tab -f avas_omin_proteins_combined.fasta -p 234 -o hexact_clusters_v1 -z 2 -s 2 -M10 -H 9 -c`
for FILE in clusters_hexact_clusters_v1/*.fasta ; do BASE="${FILE%.fasta}" ; mafft --anysymbol $FILE > $BASE.aln ; done

for FILE in clusters_hexact_clusters_v1/homologs*.aln ; do alignment_conserved_site_to_dots.py -t -n -a $FILE >> avas_omin_proteins_combined.homologs_identity.tab ; done
```

## downloads from NCBI ##

* [Octopus sinensis (East Asian common octopus)](https://www.ncbi.nlm.nih.gov/genome/?term=txid2607531[Organism:noexp])
* [Genome assembly ASM119413v2 reference Octopus bimaculoides](https://www.ncbi.nlm.nih.gov/data-hub/genome/GCF_001194135.2/)
* [Crassostrea gigas (Pacific oyster)](https://www.ncbi.nlm.nih.gov/genome/10758)
* [Crassostrea virginica (eastern oyster)](https://www.ncbi.nlm.nih.gov/genome/398)
* [Crassostrea angulata (assembly ASM2561291v2)](https://www.ncbi.nlm.nih.gov/genome/12241)
* [Daphnia pulex (common water flea)](https://www.ncbi.nlm.nih.gov/genome/288)
* [Daphnia magna (crustaceans)](https://www.ncbi.nlm.nih.gov/assembly/GCF_020631705.1)
* [Vespa crabro (European hornet)](https://www.ncbi.nlm.nih.gov/genome/?term=txid7445[orgn]&shouldredirect=false)
* [Vespa velutina](https://www.ncbi.nlm.nih.gov/genome/?term=txid202808[orgn]&shouldredirect=false)
* [Anastrepha ludens (Mexican fruit fly)](https://www.ncbi.nlm.nih.gov/genome/?term=txid28586[orgn]&shouldredirect=false)
* [Anastrepha obliqua](https://www.ncbi.nlm.nih.gov/genome/?term=txid95512[orgn]&shouldredirect=false)
* [Culex pipiens (northern house mosquito)](https://www.ncbi.nlm.nih.gov/genome/?term=txid42434[orgn]&shouldredirect=false)
* [Culex quinquefasciatus (southern house mosquito)](https://www.ncbi.nlm.nih.gov/genome/?term=txid7176[orgn]&shouldredirect=false)
* [Bombus pyrosoma](https://www.ncbi.nlm.nih.gov/genome/?term=txid396416%5Borgn%5D)
* [Bombus terrestris (buff-tailed bumblebee)](https://www.ncbi.nlm.nih.gov/genome/?term=txid30195%5Borgn%5D)
* [Cervus elaphus (red deer)](https://www.ncbi.nlm.nih.gov/genome/10790)
* [Cervus canadensis](https://www.ncbi.nlm.nih.gov/genome/34916)
* [Perca flavescens (yellow perch)](https://www.ncbi.nlm.nih.gov/genome/?term=txid8167[orgn])
* [Epinephelus lanceolatus (giant grouper)](https://www.ncbi.nlm.nih.gov/genome/?term=txid310571[orgn])
* [Epinephelus fuscoguttatus (brown-marbled grouper)](https://www.ncbi.nlm.nih.gov/genome/?term=txid293821[orgn])
* [Epinephelus moara (kelp grouper)](https://www.ncbi.nlm.nih.gov/genome/?term=txid300413[orgn])
* [Bufo bufo (common toad)](https://www.ncbi.nlm.nih.gov/genome/?term=txid8384[orgn])
* [Bufo gargarizans (Asiatic toad)](https://www.ncbi.nlm.nih.gov/genome/?term=txid30331[orgn])
* [Oreochromis aureus (blue tilapia)](https://www.ncbi.nlm.nih.gov/genome/?term=txid47969[orgn])
* [Oreochromis niloticus (Nile tilapia)](https://www.ncbi.nlm.nih.gov/genome/?term=txid8128[orgn])
* [Thunnus maccoyii (southern bluefin tuna)](https://www.ncbi.nlm.nih.gov/genome/?term=txid8240[orgn])
* [Thunnus albacares (yellowfin tuna)](https://www.ncbi.nlm.nih.gov/genome/?term=txid8236[orgn])
* [Mauremys reevesii (Reeves's turtle)](https://www.ncbi.nlm.nih.gov/genome/?term=txid260615[orgn])
* [Mauremys mutica (yellowpond turtle)](https://www.ncbi.nlm.nih.gov/genome/?term=txid74926[orgn])
* [Lytechinus pictus (painted urchin)](https://www.ncbi.nlm.nih.gov/genome/?term=txid7653[orgn])
* [Lytechinus variegatus (green sea urchin)](https://www.ncbi.nlm.nih.gov/genome/3495)


```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/013/358/895/GCF_013358895.1_ZZ_aureus/GCF_013358895.1_ZZ_aureus_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/013/358/895/GCF_013358895.1_ZZ_aureus/GCF_013358895.1_ZZ_aureus_protein.faa.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/013/358/895/GCF_013358895.1_ZZ_aureus/GCF_013358895.1_ZZ_aureus_genomic.gff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/858/045/GCF_001858045.2_O_niloticus_UMD_NMBU/GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/858/045/GCF_001858045.2_O_niloticus_UMD_NMBU/GCF_001858045.2_O_niloticus_UMD_NMBU_protein.faa.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/858/045/GCF_001858045.2_O_niloticus_UMD_NMBU/GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.gff.gz
```

```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/027/943/255/GCF_027943255.1_idAnaObli1_1.0/GCF_027943255.1_idAnaObli1_1.0_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/027/943/255/GCF_027943255.1_idAnaObli1_1.0/GCF_027943255.1_idAnaObli1_1.0_protein.faa.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/027/943/255/GCF_027943255.1_idAnaObli1_1.0/GCF_027943255.1_idAnaObli1_1.0_genomic.gff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/028/408/465/GCF_028408465.1_idAnaLude1.1/GCF_028408465.1_idAnaLude1.1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/028/408/465/GCF_028408465.1_idAnaLude1.1/GCF_028408465.1_idAnaLude1.1_protein.faa.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/028/408/465/GCF_028408465.1_idAnaLude1.1/GCF_028408465.1_idAnaLude1.1_genomic.gff.gz

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/015/732/765/GCF_015732765.1_VPISU_Cqui_1.0_pri_paternal/GCF_015732765.1_VPISU_Cqui_1.0_pri_paternal_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/015/732/765/GCF_015732765.1_VPISU_Cqui_1.0_pri_paternal/GCF_015732765.1_VPISU_Cqui_1.0_pri_paternal_protein.faa.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/015/732/765/GCF_015732765.1_VPISU_Cqui_1.0_pri_paternal/GCF_015732765.1_VPISU_Cqui_1.0_pri_paternal_genomic.gff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/801/865/GCF_016801865.2_TS_CPP_V2/GCF_016801865.2_TS_CPP_V2_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/801/865/GCF_016801865.2_TS_CPP_V2/GCF_016801865.2_TS_CPP_V2_protein.faa.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/801/865/GCF_016801865.2_TS_CPP_V2/GCF_016801865.2_TS_CPP_V2_genomic.gff.gz

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/014/825/855/GCF_014825855.1_ASM1482585v1/GCF_014825855.1_ASM1482585v1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/014/825/855/GCF_014825855.1_ASM1482585v1/GCF_014825855.1_ASM1482585v1_protein.faa.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/014/825/855/GCF_014825855.1_ASM1482585v1/GCF_014825855.1_ASM1482585v1_genomic.gff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/910/591/885/GCF_910591885.1_iyBomTerr1.2/GCF_910591885.1_iyBomTerr1.2_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/910/591/885/GCF_910591885.1_iyBomTerr1.2/GCF_910591885.1_iyBomTerr1.2_protein.faa.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/910/591/885/GCF_910591885.1_iyBomTerr1.2/GCF_910591885.1_iyBomTerr1.2_genomic.gff.gz
```

```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/015/342/785/GCF_015342785.2_UCSD_Lpic_2.1/GCF_015342785.2_UCSD_Lpic_2.1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/015/342/785/GCF_015342785.2_UCSD_Lpic_2.1/GCF_015342785.2_UCSD_Lpic_2.1_protein.faa.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/015/342/785/GCF_015342785.2_UCSD_Lpic_2.1/GCF_015342785.2_UCSD_Lpic_2.1_genomic.gff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/018/143/015/GCF_018143015.1_Lvar_3.0/GCF_018143015.1_Lvar_3.0_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/018/143/015/GCF_018143015.1_Lvar_3.0/GCF_018143015.1_Lvar_3.0_protein.faa.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/018/143/015/GCF_018143015.1_Lvar_3.0/GCF_018143015.1_Lvar_3.0_genomic.gff.gz
```



