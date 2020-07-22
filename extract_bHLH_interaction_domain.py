### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
					python extract_bHLH_interaction_domain.py
					--in <INPUT_FILE>
					--out <OUTPUT_FILE>
					"""

import re, sys

# --- end of imports --- #


def load_multiple_fasta_file( filename ):
	"""! @brief load all sequences from a multiple fasta file """
	
	sequences = {}
	
	with open( filename, "r" ) as f:
	 	header = f.readline().strip()[1:].split('\t')[0]
		line = f.readline()
		seq = ""
		while line:
			if line[0] == '>':
				try:
					sequences[ header ] 
					print header
				except:
					sequences.update( { header: seq } )
				header = line.strip()[1:].split('\t')[0]
				seq = ""
			else:
				seq += line.strip()
			line = f.readline()
		try:
			sequences[ header ] 
			print header
		except:
			sequences.update( { header: seq } )
	return sequences


def main( arguments ):
	
	myb_prot_seq_input_file = arguments[ arguments.index('--in')+1 ]
	output_file = arguments[ arguments.index('--out')+1 ]
	
	seqs = load_multiple_fasta_file( myb_prot_seq_input_file  )
	
	
	with open(  output_file, "w" ) as out:
		for key in sorted( seqs.keys() ):
			seq = seqs[ key ]
			try:
				match = re.findall( "[DE]{1}L\w{2}[RK]{1}\w{3}L\w{6}L\w{3}R", seq )[0]
				out.write( '>' + key + '\n' + match + '\n' )
			except:
				print key
				
	
if '--in' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
