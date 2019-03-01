"""
	Program: Smartpore Pipeline
	Version: 1.0
	Author: Mohamad al Kadi
	Revision: Nicolas Jung
	Description: A pipeline finding the genes position from Nanopore reads and the genome they belong to.
"""

#Native imports
import sys, os, argparse, subprocess
#Installed imports
from Bio import SeqIO
#Additional imports
import conf
from log import logger

#Global options
log=logger()
config=conf.getconf()

def main(argv):
	global log
	global config
	## Arguments definition ##
	parser = argparse.ArgumentParser("A pipeline finding the genes position from Nanopore reads and the genome they belong to.")
	parser.add_argument('-i', help='Input file. This should be the path to a fastq file.', required=True)
	parser.add_argument('-o', help='Output path. This should be the path to output our data.', required=True)
	parser.add_argument('-v', help='Verbose. Display additional information', action='store_true', required=False)
	parser.add_argument('-s', help='Runs silently, cancels verbose.', action='store_true', required=False)
	args = parser.parse_args(argv)


	## Arguments check ##
	#Required
	inputPath = args.i
	if not os.path.isfile(inputPath):
		log.tell("The input file %s doesn't exist"%(args.i))
		return
	#Check the file format
	try:
		inputRecords = list(SeqIO.parse(inputPath, "fastq"))
	except ValueError:
		log.tell("The input file must be in the fastq format.")
		return
	if len(inputRecords) == 0:
		log.tell("The input file must be in the fastq format.")
		return

	outputPath = args.o
	if not os.path.isdir(outputPath):
		log.tell("The output directory %s doesn't exist"%(args.o))
		return

	#Optional
	#More output should happen if the script runson verbose mode
	log.verbose = args.v
	#No output should happen if the script runs silently
	log.silent = args.s

	## Smartpore Pipeline ##

	#Find the strands

	#Map the reads to the genome

	#Get the bed file from our first mapping

	#Separate the bedfile between positive and negative reads

	#Generate a fasta file from each bed file and the full genome

	#Map the newly created fasta file

	#Get a sorted bam file from our mapping

	#Get the coverage for each position from our map results

	#Determining genes start and end position from the coverage_v1

	#Determining genes end positions from 3' coverage

	#Intersecting genes start and end positions from both analysis to refine our analysis

	#Combining the results



if __name__ == "__main__":
	main(sys.argv[1:])
