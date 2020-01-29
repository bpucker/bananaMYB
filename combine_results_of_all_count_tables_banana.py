### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.2 ###

## reference: https://doi.org/10.3389/fmolb.2018.00062 ###


__usage__ = """
	python combine_results_of_all_count_tables.py
	--in <FULL_PATH_TO_INPUT_DIR>
	--gff <FULL_PATH_TO_GFF3_FILE>
	--prefix <FULL_PATH_AND_PREFIX_FOR_OUTPUT>
				"""

import glob, re, sys
from operator import itemgetter

# --- end of imports --- #

def load_raw_counts( count_table, raw_counts ):
	"""! @brief load all data from given count table """
	
	tissue = count_table.split('/')[-1].split('.')[0]
	
	with open( count_table, "r" ) as f:
		f.readline()
		f.readline()
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			try:
				try:
					raw_counts[ parts[0] ].update( { tissue: int( parts[6] ) } )
				except IndexError:
					print "ERROR:" + count_table
					return raw_counts
			except KeyError:
				try:
					raw_counts.update( { parts[0]: { tissue: int( parts[6] ) } } )
				except IndexError:
					print "ERROR:" + count_table
					return raw_counts
			line = f.readline()
	return raw_counts


def load_TPMs( count_table, TPMs ):
	"""! @brief load all data from given count table """
	
	tissue = count_table.split('/')[-1].split('.')[0]
	
	lib_size_counter = 0
	with open( count_table, "r" ) as f:
		f.readline()
		f.readline()
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			try:
				lib_size_counter += int( parts[ 6 ] )
			except IndexError:
				print "ERROR:" + count_table
				return TPMs, False
			line = f.readline()
	if lib_size_counter > 0:
		factor = 1000000.0  / lib_size_counter
		with open( count_table, "r" ) as f:
			f.readline()
			f.readline()
			line = f.readline()
			while line:
				parts = line.strip().split('\t')
				try:
					TPMs[ parts[0] ].update( { tissue: int( parts[ 6 ] )*factor } )
				except KeyError:
					TPMs.update( { parts[0]: { tissue: int( parts[ 6 ] )*factor } } )
				line = f.readline()
		return TPMs, True
	else:
		return TPMs, False


def get_collapsed_exon_length( exons ):
	"""! @brief collapse overlapping exon fragments to get effective exon length """
	
	pos_to_sort = []
	for exon in exons:
		pos_to_sort.append( { 'start': exon[0], 'end': exon[1] } )
	sorted_positions = sorted( pos_to_sort, key=itemgetter('start', 'end') )
	final_exon_positions = []
	for element in sorted_positions:
		if len( final_exon_positions ) == 0:
			final_exon_positions.append( element )
		elif element['start'] < final_exon_positions[ -1 ]['end'] and element['end'] > final_exon_positions[ -1 ]['end']:
			final_exon_positions[-1] = { 'start': final_exon_positions[ -1 ]['start'], 'end': element['end'] }
		elif element['start'] > final_exon_positions[ -1 ]['end']:
			final_exon_positions.append( element )
	total_length = 0
	for exon in final_exon_positions:
		total_length += exon['end']-exon['start']
	return total_length


def get_gene_lengths( gff3_file ):
	"""! @brief get gene lengths """
	
	rna_to_gene = {}
	with open( gff3_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				if parts[2] in [ "mRNA", "transcript", "lnc_RNA", "tRNA", "rRNA", "primary_transcript", "miRNA", "snoRNA", "snRNA", "antisense_RNA", "" ]:
					try:
						ID = re.findall( "Ma\d+_t\d+", parts[-1] )[0]
						try:
							parent = re.findall( "Ma\d+_g\d+", parts[-1] )[0]
						except IndexError:
							try:
								parent = re.findall( "Ma\d+_t\d+", parts[-1] )[0]
							except:
								pass	#print line
						rna_to_gene.update( { ID: parent } )
					except IndexError:
						try:
							ID = re.findall( "Ma\d+_t\d+", parts[-1] )[0]
							try:
								parent = re.findall( "Ma\d+_g\d+", parts[-1] )[0]
							except IndexError:
								try:
									parent = re.findall( "Ma\d+_t\d+", parts[-1] )[0]
								except:
									print line
						except IndexError:
							pass	#print line
			line = f.readline()
	
	exons_per_gene = {}
	with open( gff3_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				if parts[2] in  [ "exon", "CDS", "five_prime_UTR", "three_prime_UTR" ]:
					start, end = int( parts[3] ), int( parts[4] )
					try:
						parent = re.findall( "Ma\d+_t\d+", parts[-1] )[0]
						try:
							exons_per_gene[ rna_to_gene[ parent ] ].append( ( start, end ) )
						except KeyError:
							exons_per_gene.update( { rna_to_gene[ parent ]: [ (start, end) ] } )
					except IndexError:
						try:
							parent = re.findall( "Ma\d+_g\d+", parts[-1] )[0]
							try:
								exons_per_gene[ parent ].append( ( start, end ) )
							except KeyError:
								exons_per_gene.update( { parent: [ (start, end) ] } )
						except IndexError:
							pass	#print line
			line = f.readline()
	
	# --- calculate exon length per gene by collapsing overlaps --- #
	length_per_gene = {}
	for key in exons_per_gene.keys():
		exons = exons_per_gene[ key ]
		length_per_gene.update( { key: get_collapsed_exon_length( exons ) } )
	return length_per_gene


def main( arguments ):
	"""! @brief run all parts """
	
	input_dir= arguments[ arguments.index( '--in' )+1 ]
	gff3_file = arguments[ arguments.index( '--gff' )+1 ]
	prefix = arguments[ arguments.index( '--prefix' )+1 ]
	
	input_files = glob.glob( input_dir + "*.count_table" ) + glob.glob( input_dir + "*/*.count_table" )+ glob.glob( input_dir + "*.countTable.txt" ) + glob.glob( input_dir + "*/*.countTable.txt" )
	
	print "number of identified count tables: " + str( len( input_files ) )
	
	gene_lengths = get_gene_lengths( gff3_file )
	
	raw_counts = {}
	TPMs = {}
	samples = []
	for each in input_files:
		raw_counts = load_raw_counts( each, raw_counts )
		TPMs, status = load_TPMs( each, TPMs )
		if status:
			samples.append( each.split('/')[-1].split('.')[0] )
	
	sample_names_file = prefix + "_SAMPLE_NAMES.txt"
	with open( sample_names_file, "w" ) as out:
		out.write( "\n".join( samples ) )
	
	output_raw_counts_file = prefix + "_raw_counts.txt"
	output_TPM_file = prefix + "_TPMs.txt"
	output_FPKM_file = prefix + "_FPKMs.txt"
	
	with open( output_raw_counts_file, "w" ) as out:
		out.write( "\t".join( [ "gene" ] + samples ) + '\n' )
		for gene in sorted( raw_counts.keys() ):
			new_line = [ gene ]
			for each in samples:
				try:
					new_line.append( raw_counts[ gene ][ each ] )
				except KeyError:
					new_line.append( 0 )
			out.write( "\t".join( map( str, new_line ) ) + '\n' )
	
	with open( output_TPM_file, "w" ) as out:
		out.write( "\t".join( [ "gene" ] + samples ) + '\n' )
		for gene in sorted( TPMs.keys() ):
			new_line = [ gene ]
			for each in samples:
				new_line.append( TPMs[ gene ][ each ] )
			out.write( "\t".join( map( str, new_line ) ) + '\n' )
	
	with open( output_FPKM_file, "w" ) as out:
		fin_samples = []
		for each in samples:
			fin_samples.append( each.replace( "_pass_1", "" ).replace( "_pass", "" ) )
		out.write( "\t".join( [ "gene" ] + fin_samples ) + '\n' )
		for gene in sorted( TPMs.keys() ):
			values = TPMs[ gene ]
			FPKMs = []
			try:
				gene_len = gene_lengths[ gene ]
				for each in samples:
					FPKMs.append( ( 1000.0*values[ each ] ) / gene_len )
				out.write( "\t".join( map( str, [ gene ] + FPKMs ) ) + '\n' )
			except KeyError:
				print gene
	
	
	# --- file validation --- #
	raw_counter =[]
	with open( output_raw_counts_file, "r" ) as f:
		line = f.readline()
		while line:
			raw_counter.append( line.count('\t') )
			line = f.readline()
	print list( set( raw_counter ) )
	
	tpm_counter =[]
	with open( output_raw_counts_file, "r" ) as f:
		line = f.readline()
		while line:
			tpm_counter.append( line.count('\t') )
			line = f.readline()
	print list( set( tpm_counter ) )
	
	fpkm_counter =[]
	with open( output_raw_counts_file, "r" ) as f:
		line = f.readline()
		while line:
			fpkm_counter.append( line.count('\t') )
			line = f.readline()
	print list( set( fpkm_counter ) )


if __name__ == '__main__':
	
	if '--in' in sys.argv and '--gff' in sys.argv and '--prefix' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )
