### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
					python extract_myb_domain.py
					--in <INPUT_FASTA_FILE>
					--out <OUTPUT_FASTA_FILE>
					
					bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""

import os, sys, re

# --- end of imports --- #


def load_multiple_fasta_file( filename ):
	"""! @brief load all sequences from a multiple fasta file """
	
	sequences = {}
	
	with open( filename, "r" ) as f:
	 	header = f.readline().strip()[1:]
	 	if '\t' in header:
			header = header.split('\t')[0]
		line = f.readline()
		seq = ""
		while line:
			if line[0] == '>':
				try:
					sequences[ header ] 
					print "ERROR: duplicated FASTA header: " + header
				except:
					sequences.update( { header: seq } )
				header = line.strip()[1:]
				if '\t' in header:
					header = header.split('\t')[0]
				seq = ""
			else:
				seq += line.strip()
			line = f.readline()
		try:
			sequences[ header ] 
			print "ERROR: duplicated FASTA header: " + header
		except:
			sequences.update( { header: seq } )
	return sequences


def main( arguments ):
	"""! @brief run everything """
	
	myb_prot_seq_input_file = arguments[ arguments.index( '--in' )+1 ]
	output_file = arguments[ arguments.index( '--out' )+1 ]
	
	seqs = load_multiple_fasta_file( myb_prot_seq_input_file  )
	print len( seqs.keys() )
	
	
	with open(  output_file, "w" ) as out:
		for key in sorted( seqs.keys() ):
			seq = seqs[ key ]
			try:
				match = re.findall( "\w{5}W\w{85,100}W\w{7}", seq )[0]
				out.write( '>' + key + '\n' + match + '\n' )
			except:
				print key


if '--in' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
