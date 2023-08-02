#!/usr/bin/env python
#
# get_genbank_longest_isoforms.py v1.0 2023-06-16

"""get_genbank_longest_isoforms.py  last modified 2023-08-01

get_genbank_longest_isoforms.py -g GCF_013753865.1_Amil_v2.1_genomic.gff.gz -p GCF_013753865.1_Amil_v2.1_protein.faa.gz -f GCF_013753865.1_Amil_v2.1_genomic.x.gff -o GCF_013753865.1_Amil_v2.1_protein.x.faa
"""

import sys
import argparse
import time
import gzip
from collections import defaultdict
from Bio import SeqIO

def make_exclude_dict(excludefile):
	'''read file of list of contigs, and return a dict where keys are contig names to exclude'''
	sys.stderr.write("# Reading exclusion list {}  {}\n".format(excludefile, time.asctime() ) )
	exclusion_dict = {}
	for term in open(excludefile,'r'):
		term = term.strip()
		if term[0] == ">":
			term = term[1:]
		exclusion_dict[term] = True
	sys.stderr.write("# Found {} contigs to exclude  {}\n".format( len(exclusion_dict), time.asctime() ) )
	return exclusion_dict

def read_protein_fasta(seqfile):
	"""read protein fasta file, return SeqIO dict"""
	if seqfile.rsplit('.',1)[-1]=="gz": # autodetect gzip format
		opentype = gzip.open
		sys.stderr.write("# Parsing protein fasta {} as gzipped  {}\n".format(seqfile, time.asctime() ) )
	else: # otherwise assume normal open for fasta format
		opentype = open
		sys.stderr.write("# Parsing protein fasta {}  {}\n".format(seqfile, time.asctime() ) )
	seqdict = SeqIO.to_dict(SeqIO.parse(opentype(seqfile,'rt'), "fasta"))
	sys.stderr.write("# Read {} proteins  {}\n".format( len(seqdict), time.asctime() ) )
	return seqdict

def get_protein_lengths(seqdict):
	"""from the dictionary of sequences, return a dict where key is seq ID and value is integer length"""
	lengthdict = {}
	sys.stderr.write("# Calculating protein lengths  {}\n".format( time.asctime() ) )
	for seqrec in seqdict.values():
		lengthdict[seqrec.id] = len(seqrec.seq)
	sys.stderr.write("# Got {} protein lengths  {}\n".format( len(lengthdict), time.asctime() ) )
	return lengthdict

def parse_gtf(gtffile):
	"""read gtf file and return, two dicts"""
	#NC_052623.1	Gnomon	gene	145196	298441	.	+	.	ID=gene-DCHS1;Dbxref=GeneID:120379788;Name=DCHS1;gbkey=Gene;gene=DCHS1;gene_biotype=protein_coding
	#NC_052623.1	Gnomon	mRNA	145196	298441	.	+	.	ID=rna-XM_039497313.1;Parent=gene-DCHS1;Dbxref=GeneID:120379788,Genbank:XM_039497313.1;Name=XM_039497313.1;gbkey=mRNA;gene=DCHS1;model_evidence=Supporting evidence includes similarity to: 3 Proteins%2C and 100%25 coverage of the annotated genomic feature by RNAseq alignments;product=dachsous cadherin-related 1%2C transcript variant X2;transcript_id=XM_039497313.1
	#NC_052623.1	Gnomon	exon	145196	145404	.	+	.	ID=exon-XM_039497313.1-1;Parent=rna-XM_039497313.1;Dbxref=GeneID:120379788,Genbank:XM_039497313.1;gbkey=mRNA;gene=DCHS1;product=dachsous cadherin-related 1%2C transcript variant X2;transcript_id=XM_039497313.1
	#NC_052623.1	Gnomon	CDS	159356	161107	.	+	0	ID=cds-XP_039353247.1;Parent=rna-XM_039497313.1;Dbxref=GeneID:120379788,Genbank:XP_039353247.1;Name=XP_039353247.1;gbkey=CDS;gene=DCHS1;product=protocadherin-16;protein_id=XP_039353247.1

	mrna_to_gene = {} # key is mRNA ID, value is gene ID
	gene_to_proteins = defaultdict(list) # key is gene ID, value is list of protein isoforms
	prot_to_lines = defaultdict(list) # key is prot ID, value is list of GFF lines

	used_prots = {} # key is prot, value is True

	if gtffile.rsplit('.',1)[-1]=="gz": # autodetect gzip format
		opentype = gzip.open
		sys.stderr.write("# Parsing protein fasta {} as gzipped  {}\n".format(gtffile, time.asctime() ) )
	else: # otherwise assume normal open for fasta format
		opentype = open
		sys.stderr.write("# Parsing protein fasta {}  {}\n".format(gtffile, time.asctime() ) )
	for line in opentype(gtffile, 'rt'):
		line = line.strip()
		if line and not line[0]=="#": # ignore empty lines and comments
			lsplits = line.split("\t")
			scaffold = lsplits[0]
			feature = lsplits[2]
			attributes = lsplits[8]
			attrd = dict([(field.strip().split("=",1)) for field in attributes.split(";") if field.count("=")])
			if feature=="mRNA":
				rna_id = attrd.get("ID")
				parent_gene = attrd.get("Parent")
				mrna_to_gene[rna_id] = parent_gene
			elif feature=="CDS":
				protein_id = attrd.get("protein_id")
				if used_prots.get(protein_id, False): #TODO why is this here
					pass
				parent_rna = attrd.get("Parent")
				parent_gene = mrna_to_gene.get(parent_rna)
				gene_to_proteins[parent_gene].append( protein_id )
				prot_to_lines[protein_id].append( line )
				used_prots[protein_id] = True
	sys.stderr.write("# Found {} genes for {} proteins  {}\n".format( len(gene_to_proteins), len(prot_to_lines), time.asctime() ) )
	sys.stderr.write("# Names parsed as {}  {}\n".format( parent_gene, protein_id ) )
	return gene_to_proteins, prot_to_lines

def get_longest_protein(gene_to_proteins, prot_to_lines, prot_dict, prot_lengths, outgff, outprots):
	write_counter = 0
	longest_prots = {}
	with open(outgff,'wt') as outgff:
		for gene, protset in gene_to_proteins.items():
			if gene is not None: # unsure why this happens
				prot_sd = {}
				for prot in protset:
					prot_sd[prot] = int(prot_lengths.get(prot,1))
				longest_prot_id = sorted( prot_sd.items(), key=lambda x: x[1], reverse=True)[0][0]
				longest_prots[longest_prot_id] = True
				#sys.stderr.write("# {}  {}  {}\n".format( gene, prot, longest_prot_id ) )
				write_counter += 1
				outgff.write( "{}\n".format( "\n".join( prot_to_lines.get(longest_prot_id) ) ) )
	sys.stderr.write("# Wrote {} proteins as CDS  {}\n".format( write_counter, time.asctime() ) )

	write_counter = 0
	with open(outprots,'wt') as outpro:
		for longest_prot_id in longest_prots:
			write_counter += 1
			outpro.write( prot_dict.get(longest_prot_id).format("fasta") )
	sys.stderr.write("# Wrote {} proteins  {}\n".format( write_counter, time.asctime() ) )
	return True

def main(argv, wayout):
	if not len(argv):
		argv.append("-h")
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('-g','--gff', help="GFF annotation from GenBank, can be .gz", required=True)
	parser.add_argument('-p','--proteins', help="fasta file of proteins, can be .gz", required=True)
	parser.add_argument('-f','--output-gff', help="filename of output GFF containing only the mRNA/exon/CDS of the longest protein per locus" )
	parser.add_argument('-o','--output-proteins', help="filename of output fasta of the longest protein per locus" )
	parser.add_argument('-E','--exclude', help="file of list of bad contigs, one per line")
	parser.add_argument('--print-no-match', help="print lines for all queries, including those without blast matches", action="store_true")
	args = parser.parse_args(argv)

	exclusiondict = make_exclude_dict(args.exclude) if args.exclude else {}

	prot_dict = read_protein_fasta(args.proteins)
	prot_lengths = get_protein_lengths(prot_dict)
	gene_to_proteins, prot_to_lines = parse_gtf(args.gff)
	get_longest_protein(gene_to_proteins, prot_to_lines, prot_dict, prot_lengths, args.output_gff, args.output_proteins)


if __name__ == "__main__":
	main(sys.argv[1:],sys.stdout)
