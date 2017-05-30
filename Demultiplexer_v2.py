#!/usr/bin/python
'''
This Python scripts is part of the Arabidopsis thaliana iCLIP Analysis Pipeline.
It takes a fastq-file (Illumina/Solexa) and demultiplexes the reads corresponding to
the given barcodes. Reads, that can't be sorted to the barcodes, are written into an extra file.

Contact:
Katja Meyer (katja.meyer@uni-bielefeld.de)           
Claus Weinholdt (claus.weinholdt@informatik.uni-halle.de)   
Kevin Lamkiewicz (klamkiewicz@techfak.uni-bielefeld.de) 
'''

# ---------------%Imports%--------------- #

# import for the arguments of the command line
import argparse as args
# used to exit the program with proper feedback
import sys
# easy reading for the barcode files
import csv
# special dictonary to store key:list-values
from collections import defaultdict
# chunk generation for FastQ reading
import itertools
# reg expression 
import re
# file handle 
import os 

from Bio.Seq import Seq

# ---------------%Imports%--------------- #

# -----------%Global Variables%---------- #

# -----------%Global Variables%---------- #

# -------------%File-Reading%------------ #

def readFile():
	'''
	This method reads the fastq and barcode files. Furthermore for every barcode given, a writing buffer
	is created. For each read in the fastq file, the position 3 to 6 is analysed (in our case, the barcode of the samples
	were in this range). If the current read doesn't have a valid barcode, an extra file with not sortable reads is created.
	'''

	myBarcodes = []

	countBarcode = {}
	mapForWrite = {}
	try:
		with open(arguments.barcode_file, 'r') as csvfile:
			barCodeReader = csv.reader(csvfile, delimiter=';')
			for index,row in enumerate(barCodeReader):
				#print(index)
				if len(row) == 2:
					row[1] = str(Seq(row[1][2:6]).reverse_complement()) # cut out the barcode and make the reverse_complement
					countBarcode[row[1]] = 0
					myBarcodes.append(row)
					print (row)
		myBarcodes = dict(myBarcodes)
	except IOError:
		print ('Barcode-File not found. Check your path: {}').format(arguments.barcode_file) 
		sys.exit(1)

	for key,value in myBarcodes.items():
		filePath = '{}/{}{}.fastq'.format(arguments.output_path,arguments.pre_fastq_file, key) #key is the current barcode name
		mapForWrite[value] = open(filePath, 'w') # for each barcode, another output file is created
	notSort = open('{}/{}NotSortable.fastq'.format(arguments.output_path,arguments.pre_fastq_file) , 'w') # for all reads without a valid barcode

	with open(arguments.fastq_file, 'r') as fastQfile:
		while True:
			next_four_lines = list(itertools.islice(fastQfile, 4))
			if not next_four_lines:
				break;
			written=1				
			for i in range(3,-1,-1): #3, 2, 1, 0]
				if i == 3: # check only at position 4 
					barcode = next_four_lines[1][i:i+4] #next_four_lines[1][3:7]
					if barcode in mapForWrite:
						countBarcode[barcode] += 1

						## cut barcode form sequence 
						if arguments.cut_barcode == "TRUE":
							next_four_lines[1]=next_four_lines[1][i:len(next_four_lines[1])] 
							next_four_lines[3]=next_four_lines[3][i:len(next_four_lines[3])]

						mapForWrite[barcode].write("".join(next_four_lines))
						written=0
						break			
			if written:
				notSort.write("".join(next_four_lines))
	
	for key,value in countBarcode.items():	
		print([key,value])

# -------------%File-Reading%------------ #

# -------------%Arg-Parsing%------------- #

parser = args.ArgumentParser()
parser.add_argument("fastq_file", type=str, help="Path to the fastq-file.")
parser.add_argument("barcode_file", type=str, help="Path to the file containing the Barcodes.\n Usage: Identifier;Barcode\n")
parser.add_argument("output_path", type=str, help="Path to directory where the results are saved in")
parser.add_argument("pre_fastq_file", type=str, help="Prefix for the fastq-file")
parser.add_argument("cut_barcode", type=str, help="If TRUE barcode is removed from sequence")

arguments = parser.parse_args()

if __name__ == '__main__':
	readFile()

# -------------%Arg-Parsing%------------- #
