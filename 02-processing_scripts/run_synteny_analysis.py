#!/usr/bin/env python
#
# run_synteny_analysis.py created 2023-06-23

"""run_synteny_analysis.py  v1.0 last modified 2023-06-28

run_synteny_analysis.py -i species_pair_list.tab

    species pairs are in a 12 column tabular format:

group_genus   species1   species2   spec1_code   spec2_code   scaf1   gff1   prot1   scaf2   gff2   prot2   color

perca   Perca.flavescens   Perca fluviatilis   PFLA1   Pfluv1   GCF_004354835.1_PFLA_1.0_genomic.fna.gz   GCF_004354835.1_PFLA_1.0_genomic.gff.gz   GCF_004354835.1_PFLA_1.0_protein.faa.gz   GCF_010015445.1_GENO_Pfluv_1.0_genomic.fna.gz   GCF_010015445.1_GENO_Pfluv_1.0_genomic.gff.gz   GCF_010015445.1_GENO_Pfluv_1.0_protein.faa.gz   0.51
"""

import sys
import os
import argparse
import subprocess
import time
import gzip
import glob
from collections import defaultdict
from Bio import SeqIO


def logreport(loginput, logout):
	# if input is string, then output like this:
	# This is a sample string for stderr Wed Oct 23 2013 16:31
	# or:
	# Wed Oct 23 2013 16:31
	# This is the sample string written to log
	if type(loginput) is str:
		print("{}\n{}".format(loginput, time.asctime()),  file=sys.stderr )
		print("#TIME-{}\n{}".format(time.asctime(), loginput), file=logout )
	# otherwise if input is a list, then assume it is a list of arguments and print as such
	elif type(loginput) is list:
		print("Making system call:\n{}".format(" ".join(loginput)),  file=sys.stderr )
		print("#TIME-{}\n{}".format(time.asctime(), " ".join(loginput)), file=logout )


def clean_name(namestring):
	for char in "0123456789":
		namestring = namestring.replace(char,"")
	return namestring

def get_z_value(scaffolds1, scaffolds2):
	"""read all scaffolds and return the calculated value for option z"""
	total_bp = 0
	for seqfile in [scaffolds1, scaffolds2]:
		for seqrec in SeqIO.parse(gzip.open(seqfile,'rt'),'fasta'):
			total_bp += len(seqrec.seq)
	# return average for both genomes
	return int( total_bp/5000/2 )

def zcat_proteins(zipped_proteins, is_dummy_run, logout):
	"""unzip protein gzip and return file name of unzipped"""
	# zcat GCF_018143015.1_Lvar_3.0_protein.faa.gz > GCF_018143015.1_Lvar_3.0_protein.faa
	unzipped_prots = zipped_proteins.rsplit('.',1)[0]
	zcat_Args = ['zcat', zipped_proteins]
	logreport(zcat_Args,logout)
	subprocess.run(zcat_Args, stdout=open(unzipped_prots, 'wb'), check=True)
	return os.path.abspath(unzipped_prots)


def filter_prots_and_gff(unzipped_prots, genome_gff, is_dummy_run, logout):
	"""filter a single protein from both fasta and GFF files, and return the two filenames"""
	excluded_proteins = "{}.x.faa".format( unzipped_prots.rsplit('.',2)[0] ) # should remove .faa and add .x.faa
	excluded_gff = "{}.x.gff".format( genome_gff.rsplit('.',2)[0] ) # should remove .gff.gz and add .x.gff
	ggli_Args = [ os.path.expanduser('~/project/speciation-synteny/get_genbank_longest_isoforms.py'), '-g', genome_gff, '-p', unzipped_prots, '-f', excluded_gff , '-o', excluded_proteins ]
	logreport(ggli_Args,logout)
	if is_dummy_run:
		pass
	else:
		subprocess.run(ggli_Args, check=True)
	return os.path.abspath(excluded_proteins) , os.path.abspath(excluded_gff)


def grep_excluded_isoform_names(exclude_list, unzipped_prots, is_dummy_run, logout):
	"""grep non-first isoform names from the protein list, and then generate a protein file without those names"""
	# grep -E -f  ../../excluder_list GCF_018143015.1_Lvar_3.0_protein.faa | cut -d ' ' -f 1 > GCF_018143015.1_Lvar_3.0_protein.x.names
	# excludeAinB.py GCF_018143015.1_Lvar_3.0_protein.x.names GCF_018143015.1_Lvar_3.0_protein.faa -f > GCF_018143015.1_Lvar_3.0_protein.x.faa
	exclude_names = "{}.names".format( unzipped_prots.rsplit('.',1)[0] ) # should remove .faa and add .names
	with open(exclude_names, 'w') as output:
		grep_process = subprocess.Popen(['grep', '-E', '-f', exclude_list, unzipped_prots], stdout=subprocess.PIPE)
		cut_process = subprocess.Popen(['cut', '-f', '1', '-d' ,' '], stdin=grep_process.stdout, stdout=output)
		cut_process.communicate()
	excluded_proteins = "{}.x.faa".format( unzipped_prots.rsplit('.',1)[0] ) # should remove .faa and add .x.faa
	exAB_Args = ['excludeAinB.py', exclude_names, unzipped_prots , '-f' ]
	#rmX_Args = ['remove_identical_seqs.py', '-' ]
	logreport(exAB_Args,logout)
	#logreport(rmX_Args,logout)
	subprocess.run(exAB_Args, stdout=open(excluded_proteins, 'wb'), check=True)
	#with open(excluded_proteins, 'w') as output:
	#	exAB_process = subprocess.Popen(exAB_Args, stdout=subprocess.PIPE)
	#	rmX_process = subprocess.Popen(rmX_Args, stdin=exAB_process.stdout, stdout=output)
	#	rmX_process.communicate()
	return excluded_proteins


def exclude_isoforms_from_gff(excluded_proteins, genome_gff, is_dummy_run, logout):
	"""exclude alternate isoforms from gff, and return the name of the gff without those entries"""
	# getgtfgenes.py GCF_018143015.1_Lvar_3.0_protein.x.faa GCF_018143015.1_Lvar_3.0_genomic.gff -f --genbank-gff > GCF_018143015.1_Lvar_3.0_genomic.x.gff
	excluded_gff = "{}.x.gff".format( genome_gff.rsplit('.',2)[0] ) # should remove .gff.gz and add .x.gff
	ggg_Args = ['getgtfgenes.py', excluded_proteins, genome_gff, '-f', '--genbank-gff' ]
	logreport(ggg_Args,logout)
	subprocess.run(ggg_Args, stdout=open(excluded_gff, 'wb'), check=True)
	return os.path.abspath(excluded_gff)


def diamond_makedb(protein_fasta, is_dummy_run, logout):
	"""run diamond makedb"""
	# ~/diamond-v2.0.13/diamond-v2.0.13 makedb --in GCF_015342785.2_UCSD_Lpic_2.1_protein.x.faa -d GCF_015342785.2_UCSD_Lpic_2.1_protein.x.faa
	diamond_Args = [ os.path.expanduser('~/diamond-v2.0.13/diamond-v2.0.13'), 'makedb', '--in', protein_fasta, '-d', protein_fasta]
	logreport(diamond_Args,logout)
	try:
		subprocess.run(diamond_Args, check=True)
		return True
	except subprocess.CalledProcessError:
		return None


def diamond_blastp(sp1_proteins, sp2_proteins, output_name, is_dummy_run, logout):
	"""run diamond blastp, and return the name of the output table"""
	# ~/diamond-v2.0.13/diamond-v2.0.13 blastp -q GCF_018143015.1_Lvar_3.0_protein.x.faa -d GCF_015342785.2_UCSD_Lpic_2.1_protein.x.faa -o lvar3.0_vs_lpic2.1.blastp.tab
	output_table = "{}.blastp.tab".format( output_name )
	diamond_Args = [ os.path.expanduser('~/diamond-v2.0.13/diamond-v2.0.13'), 'blastp', '-q', sp1_proteins, '-d', sp2_proteins, '-o', output_table]
	logreport(diamond_Args,logout)
	if is_dummy_run:
		pass
	else:
		subprocess.run(diamond_Args, check=True)
	return os.path.abspath(output_table)


def run_scaffold_synteny(blast_table, genome1_fasta, genome2_fasta, genome1_gff, genome2_gff, spec_name1, spec_name2, color_index, is_dummy_run, logout):
	"""run scaffold synteny and the following Rscript and return the synteny tabular file"""
	# ~/git/genomeGTFtools/scaffold_synteny.py -b lvar3.0_vs_lpic2.1.blastp.tab -q GCF_018143015.1_Lvar_3.0_genomic.gff.gz -f GCF_018143015.1_Lvar_3.0_genomic.fna.gz -d GCF_015342785.2_UCSD_Lpic_2.1_genomic.gff.gz -F GCF_015342785.2_UCSD_Lpic_2.1_genomic.fna.gz -l 860 -L 990 --ignore-gene-features --genbank-gff > lvar3.0_vs_lpic2.1.blastp.scaffold_synteny.tab
	# Rscript ~/git/genomeGTFtools/synteny_2d_plot.R lvar3.0_vs_lpic2.1.blastp.scaffold_synteny.tab L.variegatus L.pictus 100
	scaffold_synteny_tab = "{}.scaffold_synteny.tab".format( blast_table.rsplit('.',2)[0] )
	scaffold_synteny_args = [ os.path.expanduser('~/git/genomeGTFtools/scaffold_synteny.py'), '-b', blast_table, '-q', genome1_gff, '-f', genome1_fasta, '-d', genome2_gff, '-F', genome2_fasta, '--ignore-gene-features', '--genbank-gff'] # '-l', length1, '-L', length2, 
	logreport(scaffold_synteny_args,logout)
	if is_dummy_run:
		pass
	else:
		subprocess.run(scaffold_synteny_args, stdout=open(scaffold_synteny_tab, 'wb'), check=True)
		sp_Args = ['Rscript', os.path.expanduser('~/git/genomeGTFtools/synteny_2d_plot.R'), scaffold_synteny_tab, spec_name1, spec_name2, color_index]
		logreport(sp_Args,logout)
		subprocess.run(sp_Args, check=True)
	return scaffold_synteny_tab


def run_microsynteny(blast_table, genome1_gff, genome2_gff, color_index, is_dummy_run, z_value, logout):
	"""run microsynteny and return file name of the output tabular file"""
	# ~/git/genomeGTFtools/microsynteny.py -b lvar3.0_vs_lpic2.1.blastp.tab -q GCF_018143015.1_Lvar_3.0_genomic.x.gff -d GCF_015342785.2_UCSD_Lpic_2.1_genomic.x.gff --genbank-gff > lvar3.0_vs_lpic2.1.blastp.microsynteny.tab
	microsynteny_tab = "{}.microsynteny.tab".format( blast_table.rsplit('.',2)[0] )
	microsynteny_args = [os.path.expanduser('~/git/genomeGTFtools/microsynteny.py'), '-b', blast_table, '-q', genome1_gff, '-d', genome2_gff, '--genbank-gff', '-z', str(z_value), "-T" , "-s", "4" ]
	logreport(microsynteny_args,logout)
	if is_dummy_run:
		pass
	else:
		subprocess.run(microsynteny_args, stdout=open(microsynteny_tab, 'wb'), check=True)
		sp_Args = ['Rscript', os.path.expanduser('~/git/genomeGTFtools/synteny_block_length_plot.R'), microsynteny_tab, color_index]
		logreport(sp_Args,logout)
		subprocess.run(sp_Args, check=True)
	return microsynteny_tab


def run_fastarenamer(filtered_proteins, species_tag, is_dummy_run, logout):
	"""run fastarenamer and return filename of """
	# fastarenamer.py -p Lvar -j '|' GCF_018143015.1_Lvar_3.0_protein.x.faa > GCF_018143015.1_Lvar_3.0_protein.x.n.faa
	renamed_proteins = "{}.n.faa".format( filtered_proteins.rsplit('.',1)[0] )
	fr_Args = ['fastarenamer.py', '-p', species_tag, '-j', '|', filtered_proteins]
	logreport(fr_Args,logout)
	subprocess.run(fr_Args, stdout=open(renamed_proteins, 'wb'), check=True)
	return os.path.abspath(renamed_proteins)


def run_makehomologs(blast_table, protein1_fasta, protein2_fasta, genus_name, is_dummy_run, logout):
	"""run makehomologs"""
	# makehomologs.py -i lvar3.0_vs_lpic2.1.blastp.n.tab -f GCF_018143015.1_Lvar_3.0_protein.x.n.faa GCF_015342785.2_UCSD_Lpic_2.1_protein.x.n.faa -p 234 -o lytechinus_clusters_v1 -z 2 -s 2 -M 10 -H 9 -c
	# Rscript ~/git/supermatrix/plot_homolog_output_logs.R lytechinus_clusters_v1.2023-06-16-134251.mh.log fasta_clusters.H.lytechinus_clusters_v1.tab
	cluster_name = "{}_clusters_v1".format(genus_name)
	mh_Args = ['makehomologs.py', '-i', blast_table, '-f', protein1_fasta, protein2_fasta, '-p', '234', '-o', cluster_name, "-z", "2", "-s", "2", "-M" "10", "-H", "9", "-c"]
	logreport(mh_Args,logout)
	subprocess.run(mh_Args)
	cluster_folder = "clusters_{}".format(cluster_name)
	return os.path.abspath(cluster_folder)


def run_mafft_on_clusters(cluster_folder, blast_table, is_dummy_run, logout):
	"""run MAFFT on each fasta file; --anysymbol is needed as some proteins will have U or X"""
	# for FILE in clusters_lytechinus_clusters_v1/*.fasta ; do BASE="${FILE%.fasta}" ; mafft $FILE > $BASE.aln ; done
	# for FILE in clusters_lytechinus_clusters_v1/homologs_lytechinus_clusters_v1_*.aln ; do alignment_conserved_site_to_dots.py -t -n -a $FILE >> lvar3.0_vs_lpic2.1.homologs_identity.tab ; done
	file_list = glob.glob( os.path.join(cluster_folder , "homologs*.fasta") )
	concatenated_lines = []
	alignment_counter = 0
	logreport("# INFORMATION: starting MAFFT commands", logout)
	for fastafile in file_list:
		homolog_alignment = "{}.aln".format( fastafile.rsplit('.',1)[0] )
		mafft_Args = ['mafft', "--anysymbol", fastafile ]
		subprocess.run(mafft_Args, stdout=open(homolog_alignment, 'wb'), check=True)
		acsd_Args = ['alignment_conserved_site_to_dots.py', '-t', '-n', '-a', homolog_alignment ]
		acsd_process = subprocess.Popen(acsd_Args, stdout=subprocess.PIPE )
		acsd_output, _ = acsd_process.communicate()
		concatenated_lines.append(acsd_output.decode().strip())
		alignment_counter += 1
	# only report last command as an example
	logreport(mafft_Args,logout)
	logreport("# INFORMATION: aligned {} pairs with MAFFT".format( alignment_counter ), logout)

	identity_counter = 0
	cluster_identity_file = "{}.homologs_identity.tab".format( blast_table.rsplit('.',2)[0] )
	with open(cluster_identity_file,'w') as cif:
		for line in concatenated_lines:
			print(line, file=cif)
			identity_counter += 1
	logreport("# INFORMATION: Pairwise identities for {} pairs".format( identity_counter ), logout)
	return cluster_identity_file


def main(argv, wayout):
	if not len(argv):
		argv.append("-h")
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('-i','--input-table', help="text table listing all inputs")
	parser.add_argument('-l','--log-dir', help="optional directory of log files", default='./')
	parser.add_argument('-p','--pretend', help="skip actual running of programs", action="store_true")
	parser.add_argument('--genbank-gff', help="use presets when proteins and GFF files are from GenBank", action="store_true")
	args = parser.parse_args(argv)

	print("# Starting process:  {}".format( time.asctime() ), file=sys.stderr)
	startclock=time.asctime()
	starttime= time.time()

	# filename processing for output log
	log_file = "{}.synteny_analysis.log".format(time.strftime("%Y-%m-%d-%H%M%S"))
	log_path = os.path.join(args.log_dir, log_file)
	print("# Log files being written to {}:  {}".format(log_path , time.asctime() ), file=sys.stderr)
	logout = open(log_path, 'w')

	if args.pretend:
		logreport("# INFORMATION: running in PRETEND mode -p", logout)

	pair_counter = 0
	genus_tracker = defaultdict(int)
	for line in open(args.input_table):
		# each line is a pair of two species
		line = line.strip()
		if line and line[0]!="#": # skip empty or comment lines
			lsplits = line.strip().split('\t')
			# group_genus    species1    species2    spec1_code    spec2_code    scaf1    gff1    prot1    scaf2    gff2    prot2    color
			genus_name = lsplits[0]
			if genus_name =="group_genus": # meaning first line, so skip
				continue

			pair_counter += 1
			genus_tracker[genus_name] += 1
			# in case species are used twice
			if genus_tracker.get(genus_name,0) > 1: # meaning first should be species, second should be species2
				genus_code = "{}{}".format( genus_name, genus_tracker.get(genus_name) )
			else:
				genus_code = genus_name
			# this should always turn species2 into species
			genus_folder = clean_name(genus_name)

			species1_name = lsplits[1]
			species2_name = lsplits[2]
			spec1_code = lsplits[3]
			spec2_code = lsplits[4]
			scaffolds1 = os.path.join(genus_folder, lsplits[5])
			scaffolds2 = os.path.join(genus_folder, lsplits[8])
			gff1 = os.path.join(genus_folder, lsplits[6])
			gff2 = os.path.join(genus_folder, lsplits[9])
			prots1 = os.path.join(genus_folder, lsplits[7])
			prots2 = os.path.join(genus_folder, lsplits[10])
			color_index = lsplits[11]

			# assume all files are direct downloads from NCBI, so .gz
			# unzip both sets
			#unzip_prots1 = zcat_proteins(prots1, args.pretend, logout)
			#unzip_prots2 = zcat_proteins(prots2, args.pretend, logout)

			# remove isoforms from protein sets
			#filtered_prots1 = grep_excluded_isoform_names(args.exclude, unzip_prots1, args.pretend, logout)
			#filtered_prots2 = grep_excluded_isoform_names(args.exclude, unzip_prots2, args.pretend, logout)
			filtered_prots1, filtered_gff1 = filter_prots_and_gff(prots1, gff1, args.pretend, logout)
			filtered_prots2, filtered_gff2 = filter_prots_and_gff(prots2, gff2, args.pretend, logout)

			# make diamond blast db and run
			a_to_b_name = "{}/{}_vs_{}".format( genus_folder, spec1_code, spec2_code )
			is_db_made = diamond_makedb(filtered_prots2, args.pretend, logout)
			blast_table = diamond_blastp(filtered_prots1, filtered_prots2, a_to_b_name, args.pretend, logout)

			# make scaffold synteny table and dot plot
			ss_color = str(int(float(color_index)*255))
			syn_table = run_scaffold_synteny(blast_table, scaffolds1, scaffolds2, gff1, gff2, species1_name, species2_name, ss_color, args.pretend, logout)

			# make microsynteny table and plot
			#filtered_gff1 = exclude_isoforms_from_gff(filtered_prots1, gff1, args.pretend, logout)
			#filtered_gff2 = exclude_isoforms_from_gff(filtered_prots2, gff2, args.pretend, logout)
			z_value = get_z_value(scaffolds1, scaffolds2)
			run_microsynteny(blast_table, filtered_gff1, filtered_gff2, color_index, args.pretend, z_value, logout)

			# rename prots by code for rapid indexing in makehomologs
			rename_prots1 = run_fastarenamer(filtered_prots1, spec1_code, args.pretend, logout)
			rename_prots2 = run_fastarenamer(filtered_prots2, spec2_code, args.pretend, logout)
			a_to_b_rename = "{}.renamed".format( a_to_b_name )
			is_db_made = diamond_makedb(rename_prots2, args.pretend, logout)
			blast_table_rn = diamond_blastp(rename_prots1, rename_prots2, a_to_b_rename, args.pretend, logout)
			cluster_folder = run_makehomologs(blast_table_rn, rename_prots1, rename_prots2, genus_code, args.pretend, logout)

			# run mafft on all clusters
			homolog_identity_tab = run_mafft_on_clusters(cluster_folder, blast_table, args.pretend, logout)

	logreport("#TIME-Process completed in %.1f minutes" % ((time.time()-starttime)/60), logout )
	print( "Done:  {}".format(time.asctime()) , file=sys.stderr)
	logout.close()

if __name__ == "__main__":
	main(sys.argv[1:],sys.stdout)
