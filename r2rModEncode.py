

# merge read-read alignment + with read-genome alignment 
# for base modifications

# iterate through genome aligned read + and then use read-read alignment index
# for fast random access via read name (won't work vice versa because chroms are indexed)

#check if MM/ML exists and update accordingly





# Overview
# Grab Reads
# Grab Aligned Pairs for Read2Read Alignment
# For Modified Base of interest [i.e. m6(A), m5(C)]:
#		Get forward positions of modified base
# 		Use get aligned pairs to convert into ref-read coordinates	
#		Convert into MM tag format



import pysam 
import sys
import numpy as np 
import gzip
import subprocess
import shlex
import argparse
import tqdm

def inputArgs():
	'''Parse in arguments. '''
	
	parser = argparse.ArgumentParser()

	parser.add_argument('-r','--r2r', 
		type = str, 
		help = 'read-2-read bam alignment with mods in MM/ML field added by modkit repair')

	parser.add_argument('-o','--output',
		type=str,
		default = 'aligned_reads.r2r_mods.bam',
		help = 'output bam name')

	parser.add_argument('-a','--alignment',
		type=str,
		help = 'ref-reads aligned to ref genome bam file')

	parser.add_argument('-b','--bases',
		type=str,
		help = 'modification base, options = [ A,C ]')

	parser.add_argument('--replace', action='store_true')
	args = parser.parse_args()

	return args.r2r, args.output, args.alignment, args.bases, args.replace

def coordinateConversion_MMTag(sequence, base, modification_coords):
	'''Sequence is array of bases. Base is a base 
		for conversion to the coordinate system 
		used in the MM tag.'''

	mask = sequence == base
 	# find all masks = boolean Array 
	
	coords = modification_coords[ sequence[modification_coords] == base] 
	
	# when working with double stranded data we can only use modifications
	# on the specifc base that the modification would fall on on that strand
	# i.e. As  on + , Ts on -, we only want the mods that = that base of interest

	MM_coords = ','.join(list((np.diff(np.cumsum(mask)[coords]) - 1).astype(str)))

	return MM_coords 


def processAlignments(r2r_bam, alignment_bam, output_bam, base,replaceMM):

	if base == "A":
		tag = 'A+a'
		modkey = ('A', 0, 'a')
		byteBase = bytes(base, encoding='utf-8')
	elif base == "C":
		tag = 'C+m'
		modkey = ('C', 0, 'm')
		byteBase = bytes(base, encoding='utf-8')
	else:
		sys.stderr.write('Base not supported, only A/C supported currently\n')
		exit()
	for ref_read in tqdm.tqdm(alignment_bam):
		
		for r2r_read in r2r_bam.fetch(str(ref_read.query_name), until_eof=True):
			if r2r_read.is_mapped:
				if not r2r_read.is_secondary and not r2r_read.is_supplementary:
					
					mods = r2r_read.modified_bases_forward
					
					forward_mods = np.array(mods[modkey]).astype(int).T #0th is position, 1th is quality
					
					aligned_pairs = np.array(r2r_read.get_aligned_pairs(with_seq=False,matches_only=True)).astype(int).T #0th is read, 1th is ref_read, 2nd is ref base

					idx = np.isin(aligned_pairs[0],forward_mods[0],assume_unique=True)

					positions = aligned_pairs[1][idx] # new positions from sup model

					quals = forward_mods[1][np.isin(forward_mods[0],aligned_pairs[0],assume_unique=True)]

					forseq = ref_read.get_forward_sequence()

					first_occurence = forseq.index(base)

					adjust_positions = np.concatenate([[first_occurence],positions]) 
					# need to get the position of the first T because everything is offset from there 
					# when calculating MM tag
					
					sequence = np.frombuffer(bytes(forseq, "utf-8"), dtype="S1")

					MM_coords = coordinateConversion_MMTag(sequence,byteBase,adjust_positions)
					if len(MM_coords) == 0:
						continue
					MM_tag = tag +"," + MM_coords
					ML_tag = [int(ml) for ml in list(quals)]
					
					if replaceMM or not ref_read.has_tag("MM"):
					
						ref_read.set_tag("MM", MM_tag, replace=True)
						ref_read.set_tag("ML", ML_tag, replace=True)

					else:

						old_mm = ref_read.get_tag("MM")
						old_ml = ref_read.get_tag("ML")

						
						newMM = old_mm + ";" + MM_tag
						newML = list(old_ml) + ML_tag
						
						ref_read.set_tag("MM", newMM, replace=True)
						ref_read.set_tag("ML", newML, replace=True)

					output_bam.write(ref_read)
					break 
				

def main():

	#get input args
	r2r, output, alignment, base,replaceMM = inputArgs()


	alignment_bam = pysam.AlignmentFile(alignment,'rb',check_sq=False) 
	
	header = alignment_bam.header.to_dict()

	output_bam = pysam.AlignmentFile(output, "wb", header = header)

	r2r_bam = pysam.AlignmentFile(r2r,'rb',check_sq=False) 

	processAlignments(r2r_bam, alignment_bam, output_bam, base,replaceMM)




	
if __name__=="__main__":
	main()

