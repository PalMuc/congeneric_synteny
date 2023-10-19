#!/usr/bin/env python
#
# v3.5 allow gzip input 2023-10-07
# v3.4 python3 update 2019-10-12
# v3.3 allow standard input 2019-03-05
# v3.2 added X removal from seqnamer.py 2014-06-10
# v3.1 fixed bug with calling string split 0 2014-05-16
# v3.0 removed biopython dependence 2014-04-24
# v2.2 added exception for split for multiple naming types in one fasta 2014-03-13
# v2.1 added replace as in seqnamer.py 2013-04-08
# v2.0 adapted from renamefasta and copyblastserver 2013-04-01
# v1.1 checks version of python and generates string accordingly 12/12/12
# v1.0 takes a term and adds to the ID of a fasta sequence, with the number 01/08/2012

'''script to rename fasta files    last modified 2023-10-07
fastarenamer.py -a "_S1" seqsfile.fa > renamedseqs.fa

    to append name and sequence numbering
fastarenamer.py -p NAME -n seqsfile.fa > renamedseqs.fa
    so:
>Locus_1_Transcript_1/15
TTGGATATCACCAATAGTTAGGGCGACGAAGCTGTCATTATTACGCGTGTAATATCCGAG

    will print as:
>NAME_000001_Locus_1_Transcript_1/15
TTGGATATCACCAATAGTTAGGGCGACGAAGCTGTCATTATTACGCGTGTAATATCCGAG

    for replacement with -r, all characters in the string are replaced
>seq1|ABCD.1234|GENE
    -r ".|" becomes:
>seq1_ABCD_2134_GENE

    can take stdin, if not using numbering mode (-n), e.g. from gzip
gzip -dc seqsfile.fa.gz | fastarenamer.py -p Genus_species - | gzip > seqsfile.renamed.fa.gz
'''

import sys
import argparse
import os
import gzip

def seqcount_to_fstring(fastafile):
	fastacount = 0
	with open(fastafile,'r') as ff:
		for line in ff:
			if line[0]==">":
				fastacount += 1
	# literally counting the number of digits in the length
	reclog = str(len(str(fastacount)))
	# extract version information from sys, and change string accordingly
	# for python 2.7 or greater, use "_{:06}_" instead of "_{0:06}_", as field identifiers are not required
	if sys.version_info[1] >= 7:
		return "{:0"+reclog+"}"
	else:
		return "{0:0"+reclog+"}"

def main(argv, wayout):
	if not len(argv):
		argv.append("-h")
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('input_file', type=argparse.FileType('rU'), default = '-', help="fasta format file, can be .gz, can be stdin as -")
	parser.add_argument('-a','--append', help="string to append to each seq ID", default="")
	parser.add_argument('-p','--prepend', help="string to prepend to each seq ID", default="")
	parser.add_argument('-n','--number', action='store_true', help="prepend a unique number to each ID, cannot use with stdin")
	parser.add_argument('-d','--delimiter', help="delimiter for split, default ' '", default=" ")
	parser.add_argument('-j','--join', help="character to join append or prepend, default '_'", default="_")
	parser.add_argument('-r','--replace', metavar="CHARS", help="replace characters in the header with '_'")
	parser.add_argument('--replace-char', metavar="CHAR", help="character to substitute for -r, default is '_'", default="_")
	parser.add_argument('-s','--split', type=int, help="split piece using delimiter -d, ex: 0", default=None)
	parser.add_argument('-S','--swissprot', action='store_true', help="defaults for SwissProt IDs")
	parser.add_argument('-T','--trinity', action='store_true', help="defaults for TRINITY v2.4 IDs")
	parser.add_argument('-v','--verbose', action='store_true', help="verbose output")
	parser.add_argument('-x','--xreplace', action='store_true', help="replace nucleotide sequence X with N")
	args = parser.parse_args(argv)

	# if numbering sequences with padded zeroes, read the entire file once
	if args.number:
		fstring = seqcount_to_fstring(args.input_file)

	count = 0
	if args.swissprot: # reassign some of the namespace values
		args.delimiter = "|"
		args.split = 2

	if args.input_file.name.rsplit('.',1)[-1]=="gz":
		input_handler = gzip.open(args.input_file.name,'rt')
	else:
		input_handler = args.input_file
	for line in input_handler:
		line = line.strip()
		if line: # remove blank lines
			if line[0]==">": # if fasta header
				count+=1
				# take everything but the first character
				fastaline = line[1:]
				# changed to nonetype detection, to allow for splitting index 0
				if not args.split is None:
					# try to split anyway
					try:
						fastaline = fastaline.split(args.delimiter)[args.split]
					except IndexError: # if it does not work, then ignore and continue
						pass
					# because velvet transcripts are in the form of:
					# Locus_2478_Transcript_1/4_Confidence_0.444_Length_2024
					# if the "length_xx" part of the name should be removed, use:
					#seq.id = outdir_name+fstring.format(count)+seq.id.rsplit("_",2)[0]
				if args.swissprot: # from >tr|E1BRG3|E1BRG3_CHICK
				# convert >E1BRG3_CHICK Uncharacterized protein OS=Gallus gallus GN=KDM1B PE=4 SV=2
				# into >E1BRG3_CHICK
					fastaline = fastaline.split()[0]
				if args.trinity:
					fastaline = fastaline.split()[0].replace("TRINITY_","")
				if args.replace:
					for n in args.replace:
						fastaline = fastaline.replace(n,args.replace_char)
				if args.number:
					outitemlist = [args.prepend, fstring.format(count), fastaline, args.append]
				else:
					outitemlist = [args.prepend, fastaline, args.append]
				modline = ">" + (args.join).join([x for x in outitemlist if x])
				wayout.write(modline + os.linesep)
			else: # otherwise print line
				if args.xreplace:
					wayout.write(line.replace("X","N") + os.linesep)
				else:
					wayout.write(line + os.linesep)

if __name__ == "__main__":
	main(sys.argv[1:], sys.stdout)

