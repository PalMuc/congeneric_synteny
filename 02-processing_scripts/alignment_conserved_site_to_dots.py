#!/usr/bin/env python
#
# alignment_dots.py v1 2017-08-15
# renamed alignment_conserved_site_to_dots.py v1.1 2023-01-31 python3 update
# v1.2 2023-06-21 add no print mode
# v1.21 2023-07-26 specify which sequence had gaps
# v1.3 2023-10-09 allow binary output mode for indel ASR

'''alignment_conserved_site_to_dots.py  last modified 2023-10-09
within an alignment, convert conserved letters to dots

alignment_conserved_site_to_dots.py -a mox_all.aln -s DOPO_HUMAN > mox_all_dot_version.aln

    use -t to print tabular stats for each sequence pairwise with taxon of -s
    the 11 fields are:
target  taxon  gaps  gap%  identities  id%  differences diff%  id/diff  q_gaps  s_gaps
'''

import sys
import time
import argparse
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO

def dotted_alignment(alignment_file, alignformat, target_seqid, do_binary, do_tabular, skip_print, wayout):
	sys.stderr.write( "# Reading alignment from {}\n".format( alignment_file ) )
	alignment = AlignIO.read( alignment_file, alignformat )

	al_length = alignment.get_alignment_length()
	num_taxa = len(alignment)

	sys.stderr.write( "# Alignment contains {} taxa for {} sites, including gaps\n".format( num_taxa, al_length ) )

	targetseq = None
	newalignment = None
	for i,seqrec in enumerate(alignment):
		if target_seqid is None and i==0:
			target_seqid = str(seqrec.id)
			sys.stderr.write( "# Defaulting to first sequence {}\n".format( target_seqid ) )
		if seqrec.id==target_seqid:
			targetseq = seqrec.seq
			newalignment = alignment[i:i+1]
	if targetseq is None:
		sys.stderr.write( "# ERROR: CANNOT FIND SEQUENCE {}, CHECK OPTION -s OR ALIGNMENT\n".format( target_seqid ) )
		return None

	for seqrec in alignment:
		nullcount = 0 # should equal t_gaps + s_gaps + any X
		t_gapcount = 0
		s_gapcount = 0
		identical = 0
		if do_binary: # meaning 01 mode
			binstring = "".join( ['1' if s == '-' else '0' for s in seqrec.seq ] )
			newalignment.append(  SeqRecord( seq=Seq(binstring),id=seqrec.id, description=seqrec.description )  )
		else: # meaning dot mode
			if seqrec.id==target_seqid:
				continue # skip target sequence, since used as base
			else:
				dottedstring = ""
				for tl,sl in zip( targetseq, seqrec.seq ): # iterate through the two sequences
					if tl=="-":
						t_gapcount+= 1
					if sl=="-":
						s_gapcount+= 1
						# add check for both being gaps
						if tl=="-" and al_length==2: # meaning only 2 sequences, and both are gaps
							sys.stderr.write( "# WARNING: alignment {} of length {} has gaps in both sequences {} {} \n".format( alignment_file, al_length, seqrec.id, target_seqid ) )
					if tl=="-" or sl=="-" or tl=="X" or sl=="X": # if either one is gaps, always display gaps
						dottedstring += sl
						nullcount += 1
					elif tl==sl:
						dottedstring += "."
						identical += 1
					else: # implying not the same and neither is gap
						dottedstring += sl
				newalignment.append(  SeqRecord( seq=Seq(dottedstring),id=seqrec.id, description=seqrec.description )  )
			differences = al_length - nullcount - identical
			try: # make string of ratio
				ident_to_diff_ratio = "{:.2f}".format( identical*1.0/differences )
			except ZeroDivisionError: # so that 9999 appears as integer and not 9999.00
				ident_to_diff_ratio = "9999"
			if do_tabular: # meant to import as a table in something else
				outline = "{}\t{}\t{}\t{:.1}\t{}\t{:.1}\t{}\t{:.1}\t{}\t{}\t{}\n".format( target_seqid, seqrec.id, nullcount, nullcount*1.0/al_length, identical, identical*1.0/al_length, differences, differences*1.0/al_length, ident_to_diff_ratio, t_gapcount, s_gapcount )
				if skip_print:
					wayout.write( outline )
				else:
					sys.stderr.write( outline )
			else: # include tags for each value
				sys.stderr.write( "{}: gaps:{} ({:.1%}) ; identites:{} ({:.1%}) ; differences:{} ({:.1%})\n".format( seqrec.id, nullcount, nullcount*1.0/al_length, identical, identical*1.0/al_length, differences, differences*1.0/al_length ) )
	if not skip_print:
		AlignIO.write([newalignment[int(do_binary):]], wayout, alignformat)

def main(argv, wayout):
	if not len(argv):
		argv.append('-h')
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument("-a","--alignment", help="multiple sequence alignment", required=True)
	parser.add_argument("-f","--format", default="fasta", help="alignment format [fasta]")
	parser.add_argument("-b","--binary", help="convert alignment to 0 and 1, for sites and gaps, where gaps are 1", action="store_true")
	parser.add_argument("-n","--no-print", help="only calculate, do not print alignment", action="store_true")
	parser.add_argument("-s","--sequence", help="sequence ID for reference sequence, otherwise defaults to first sequence")
	parser.add_argument("-t","--tabular", help="report stats as tabular, cannot use -b", action="store_true")
	args = parser.parse_args(argv)

	if args.tabular and args.binary:
		sys.stderr.write( "# WARNING: binary mode -b cannot be used with tabular mode -t, ignoring tabular\n" )
	dotted_alignment( args.alignment, args.format, args.sequence, args.binary, args.tabular, args.no_print, wayout )

if __name__ == "__main__":
	main(sys.argv[1:], sys.stdout)
