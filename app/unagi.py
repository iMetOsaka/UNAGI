"""
	Program: UNAGI Pipeline
	Version: 1.0.1
	Author: Mohamad al Kadi
	Revision: Nicolas Jung
	Description: A pipeline finding the genes position from Nanopore reads and the genome they belong to.
"""

#Native imports
import sys, os, argparse, subprocess, gzip
#Installed imports
from Bio import SeqIO
#Additional imports
from log import logger
import conf

#Global options
log=logger()
config=conf.getconf()
version="v1.0.1"

def main(argv):
	global log
	global config
	global version
	## Arguments definition ##
	parser = argparse.ArgumentParser("UNAGI pipeline.")
	parser.add_argument('-v', '--version', action='version', version="%(prog)s ("+version+")")
	parser.add_argument('-i', '--input', help='Input file. This should be the path to a fastq file.', required=True)
	parser.add_argument('-o', '--output', help='Output path. This should be the path to output the result files to.', required=True)
	parser.add_argument('-g', '--genome', help='Genome File. This should be a valid genome file for minimap.', required=True)
	parser.add_argument('-s', '--stranded', help='Stranded. Should be selected if the input file contains reads that are already stranded.', action='store_true', required=False)
	parser.add_argument('-V', '--verbose', help='Verbose. Displays additional information while running.', action='store_true', required=False)
	parser.add_argument('-S', '--silent', help='Runs silently (no console output), cancels verbose.', action='store_true', required=False)
	args = parser.parse_args(argv)


	#Config file check
	if config == None:
		log.tell("The configuration file is not ini formatted")
		return

	## Arguments check ##
	#Required
	#Input
	inputFile = args.input
	if not os.path.isfile(inputFile):
		log.tell("The input file %s doesn't exist"%(args.input))
		return
	#Check the file format
	try:
		log.tell("Reading the input file")
		inputRecords = list(SeqIO.parse(inputFile, "fastq"))
	except ValueError:
		try:
			with gzip.open(inputFile, "rt") as inputArchive:
				inputRecords = list(SeqIO.parse(inputArchive, "fastq"))
		except ValueError:
			log.tell("The input file must be in the fastq or fastq.gz format.")
			return

	if len(inputRecords) == 0:
		log.tell("The input file must be in the fastq or fastq.gz format.")
		return

	#Output
	outputPath = args.output
	if not os.path.isdir(outputPath):
		log.tell("The output directory %s doesn't exist, creating it"%(outputPath))
		os.mkdir(outputPath)

	transitionnalOutputPath=os.path.join(outputPath,"transitionnal")
	if not os.path.isdir(transitionnalOutputPath):
		os.mkdir(transitionnalOutputPath)


	#Genome
	genomeFile = args.genome
	if not os.path.isfile(genomeFile):
		log.tell("The genome file %s doesn't exist"%(args.genome))
		return

	#Optional
	#If the -s option is selected, the stranded file is the input file. Otherwise, it has to be generated
	stranding = not args.stranded
	if stranding:
		strandedFile = os.path.join(transitionnalOutputPath,config["stranded_file"])
	else:
		strandedRecords = inputRecords
		strandedFile = inputFile

	#More output should happen if the script runson verbose mode
	log.verbose = args.verbose
	#No output should happen if the script runs silently
	log.silent = args.silent

	## UNAGI Pipeline ##

	if(stranding):
		log.tell("Finding the different strands in the input reads")
		#Find the strands if needed
		strandedRecords = findStrands(inputRecords)
		#Writing the stranded reads to a fastq file
		with gzip.open(strandedFile, "wt") as strandedGzFile:
			SeqIO.write(strandedRecords, strandedGzFile, "fastq")
		log.tell("A total of %i records out of %i (%i%%) were successfully stranded"%(len(strandedRecords), len(inputRecords), round(len(strandedRecords)*100/len(inputRecords))))

	#Clearing memory of huge variables:
	inputRecords=None
	del inputRecords
	strandedRecords=None
	del strandedRecords

	#Map the reads to the genome
	log.tell("Mapping the reads to the genome")
	minimap(strandedFile, genomeFile, os.path.join(transitionnalOutputPath,config["raw_mapped_sam_file"]))

	#Get a bam file from the results and sort it
	log.tell("Sorting the mapped reads")
	samToBam(os.path.join(transitionnalOutputPath,config["raw_mapped_sam_file"]),os.path.join(transitionnalOutputPath,config["raw_mapped_bam_file"]))
	sortBam(os.path.join(transitionnalOutputPath,config["raw_mapped_bam_file"]),os.path.join(transitionnalOutputPath,config["raw_sorted_bam_file"]))

	#Get the bed file from our first mapping and separating it between positive and negative reads
	log.tell("Generating the positive and negative bed files")
	bamToBed(os.path.join(transitionnalOutputPath,config["raw_sorted_bam_file"]),os.path.join(transitionnalOutputPath,config["raw_sorted_bed_file"]))
	bamToBed(os.path.join(transitionnalOutputPath,config["raw_sorted_bam_file"]),os.path.join(transitionnalOutputPath,config["raw_sorted_split_bed_file"]),split=True)
	separateBedFile(os.path.join(transitionnalOutputPath,config["raw_sorted_bed_file"]),os.path.join(transitionnalOutputPath,config["positive_bed_file"]),os.path.join(transitionnalOutputPath,config["negative_bed_file"]))

	#Generate a fasta file from each bed file and the full genome
	log.tell("Generating fasta files from our bed files and the genome")
	bedToFasta(os.path.join(transitionnalOutputPath,config["positive_bed_file"]),genomeFile,os.path.join(transitionnalOutputPath,config["positive_refined_fasta_file"]))
	bedToFasta(os.path.join(transitionnalOutputPath,config["negative_bed_file"]),genomeFile,os.path.join(transitionnalOutputPath,config["negative_refined_fasta_file"]))

	#Map the newly created fasta file
	log.tell("Mapping the new fasta files to the genome")
	minimap(os.path.join(transitionnalOutputPath,config["positive_refined_fasta_file"]), genomeFile, os.path.join(transitionnalOutputPath,config["positive_mapped_sam_file"]))
	minimap(os.path.join(transitionnalOutputPath,config["negative_refined_fasta_file"]), genomeFile, os.path.join(transitionnalOutputPath,config["negative_mapped_sam_file"]))

	#Get a sorted bam file from our mapping
	log.tell("Sorting the new mapped reads")
	samToBam(os.path.join(transitionnalOutputPath,config["positive_mapped_sam_file"]),os.path.join(transitionnalOutputPath,config["positive_mapped_bam_file"]))
	samToBam(os.path.join(transitionnalOutputPath,config["negative_mapped_sam_file"]),os.path.join(transitionnalOutputPath,config["negative_mapped_bam_file"]))
	sortBam(os.path.join(transitionnalOutputPath,config["positive_mapped_bam_file"]),os.path.join(transitionnalOutputPath,config["positive_sorted_bam_file"]))
	sortBam(os.path.join(transitionnalOutputPath,config["negative_mapped_bam_file"]),os.path.join(transitionnalOutputPath,config["negative_sorted_bam_file"]))

	#Get the coverage for each position from our map results
	log.tell("Generating the genome coverage for each position")
	genomecov(os.path.join(transitionnalOutputPath,config["positive_sorted_bam_file"]),genomeFile,os.path.join(transitionnalOutputPath,config["positive_coverage_file"]))
	genomecov(os.path.join(transitionnalOutputPath,config["negative_sorted_bam_file"]),genomeFile,os.path.join(transitionnalOutputPath,config["negative_coverage_file"]))
	combineCoverage(os.path.join(transitionnalOutputPath,config["positive_coverage_file"]),os.path.join(transitionnalOutputPath,config["negative_coverage_file"]),os.path.join(transitionnalOutputPath,config["total_coverage_file"]))

	#Determine genes start and end positions from the coverage, store the lists in variable for later use
	log.tell("Determining genes start and end positions from the coverage")
	positiveChrList = findCoverageBreaks(os.path.join(transitionnalOutputPath,config["positive_coverage_file"]),os.path.join(transitionnalOutputPath,config["positive_genelist_incomplete_file"]),"+")
	negativeChrList = findCoverageBreaks(os.path.join(transitionnalOutputPath,config["negative_coverage_file"]),os.path.join(transitionnalOutputPath,config["negative_genelist_incomplete_file"]),"-")

	#Determine genes end positions from 3' coverage
	log.tell("Determining genes end positions from 3' coverage")
	genomecov3dash(os.path.join(transitionnalOutputPath,config["positive_sorted_bam_file"]),os.path.join(transitionnalOutputPath,config["positive_3dash_file"]))
	genomecov3dash(os.path.join(transitionnalOutputPath,config["negative_sorted_bam_file"]),os.path.join(transitionnalOutputPath,config["negative_3dash_file"]))

	#Determine genes start positions from 5' coverage
	log.tell("Determining genes start positions from 5' coverage")
	genomecov5dash(os.path.join(transitionnalOutputPath,config["positive_sorted_bam_file"]),os.path.join(transitionnalOutputPath,config["positive_5dash_file"]))
	genomecov5dash(os.path.join(transitionnalOutputPath,config["negative_sorted_bam_file"]),os.path.join(transitionnalOutputPath,config["negative_5dash_file"]))

	#Intersect genes start and end positions from the coverage analysis and cuts to refine our analysis
	log.tell("Intersecting genes start and end positions from the coverage analysis and cuts")
	positiveChrList = insertDashCuts(os.path.join(transitionnalOutputPath,config["positive_3dash_file"]), os.path.join(transitionnalOutputPath,config["positive_5dash_file"]), positiveChrList, os.path.join(transitionnalOutputPath,config["positive_genelist_file"]), "+")
	negativeChrList = insertDashCuts(os.path.join(transitionnalOutputPath,config["negative_3dash_file"]), os.path.join(transitionnalOutputPath,config["negative_5dash_file"]), negativeChrList, os.path.join(transitionnalOutputPath,config["negative_genelist_file"]), "-")

	#Combining the results
	log.tell("Combining the results in a single file")
	combineChrLists(positiveChrList,negativeChrList,os.path.join(outputPath,config["full_genelist_file"]))

	#Finding splicing isoforms
	log.tell("Finding splicing isoforms")
	findSpliceSites(os.path.join(transitionnalOutputPath,config["raw_sorted_split_bed_file"]),os.path.join(transitionnalOutputPath,config["raw_splice_sites_file"]))
	filterByCoverage(os.path.join(transitionnalOutputPath,config["raw_splice_sites_file"]),os.path.join(transitionnalOutputPath,config["total_coverage_file"]),os.path.join(transitionnalOutputPath,config["coverage_filtered_splice_sites_file"]))
	filterByGenome(os.path.join(transitionnalOutputPath,config["coverage_filtered_splice_sites_file"]),genomeFile,os.path.join(transitionnalOutputPath,config["genome_filtered_splice_sites_file"]))
	log.tell("Selecting the relevant isoforms")
	findBestIsoform(os.path.join(transitionnalOutputPath,config["genome_filtered_splice_sites_file"]),os.path.join(transitionnalOutputPath,config["best_filtered_splice_sites_file"]))
	#Translating to bed format for the final results
	convertIsoformstoBed(os.path.join(transitionnalOutputPath,config["genome_filtered_splice_sites_file"]),os.path.join(outputPath,config["final_splice_sites_file"]))
	convertIsoformstoBed(os.path.join(transitionnalOutputPath,config["best_filtered_splice_sites_file"]),os.path.join(outputPath,config["final_unique_splice_sites_file"]))
	# --> UTR Isoforms

	#End message
	log.tell("The full gene list has been generated at: %s"%(os.path.join(outputPath,config["full_genelist_file"])))

## Custom Functions ##

#Find he forwards and backwards strands from a set of unstranded records
def findStrands(records):
	global config
	global log

	#Find the forwards reads
	log.write("Finding the forwards reads.")
	forwards=[]
	for record in records:
		if config["strand_forwards_identifier"] in record.seq[int(config["strand_start_index"]):int(config["strand_end_index"])]:
			forwards.append(record)
	log.write("%i forwards reads found"%(len(forwards)))

	#Find the backwards reads
	log.write("Finding the backwards reads.")
	backwards=[]
	for record in records:
		 if config["strand_backwards_identifier"] in record.seq[int(config["strand_start_index"]):int(config["strand_end_index"])]:
 			backwards.append(record)
	log.write("%i backwards reads found"%(len(backwards)))

	#Getting the reverse complement for the backwards reads
	log.write("Getting the reverse complement for the backwards reads.")
	for record in backwards:
		record.seq = record.seq.reverse_complement()

	return forwards + backwards

#Sorts the entries in a bed file into postive an negative ones
def separateBedFile(sourceFile,positiveFile,negativeFile):
	global config
	global log
	#counters to give a feedack on the number of entries in each file at the end
	pos=0
	neg=0
	#The sign is in 6th position so for each line, we write it into the postive or negative file accordingly
	with open(positiveFile, "w") as positive, open(negativeFile, "w") as negative, open(sourceFile, "r") as source:
		for line in source:
			parts=line.strip().split()
			if len(parts) >= 6 and parts[5] == "+":
				positive.write(line)
				pos=pos+1
			elif len(parts) >= 6 and parts[5] == "-":
				negative.write(line)
				neg=neg+1
	log.write("Found %i positive entries and %i negative ones.")

#Creates a table of start and stop position of genes for each chromosomes based on a coverage threshold
#If the strand is specified, the strands of all entries will be set to it in the resulting bed file
def findCoverageBreaks(sourceFile, outputFile, strand="."):
	global config
	global log
	chromosomes=dict() 	#The list of all chromosomes with genes position for each one
	currentGenes=list() #The list of gene positions for the currently parsed chromosome
	currentGene={"start":None,"end":None,"strand":strand} #The current gene for which we search the position
	currentChr=None	#The name of the currently parsed chromosome
	geneNumber=0	#The number of genes found for information
	with open(sourceFile, "r") as coverage:
		for line in coverage:
			chr, pos, cov=line.strip().split("\t")
			#If we reach a new chromosome, we save the current genes positions in an array and begin a new one
			if chr != currentChr:
				if currentChr is not None: chromosomes[currentChr]=list(currentGenes)
				currentChr=chr
				currentGenes=list()
				currentGene={"start":None,"end":None,"strand":strand}

			#In any case try to find the start of a new gene if we don't have one
			if currentGene["start"] is None:
				if int(cov) > int(config["min_coverage_threshold"]): currentGene["start"]=int(pos)
			#Or the end of the current gene if we already have a start
			else:
				#When we find the end of a gene, we register it in our list and start lookin for a new one
				if int(cov) < int(config["min_coverage_threshold"]):
					currentGene["end"]=int(pos)
					currentGenes.append(dict(currentGene))
					currentGene={"start":None,"end":None,"strand":strand}
					geneNumber = geneNumber+1

	#Add the last chromsome list
	if currentChr is not None: chromosomes[currentChr]=list(currentGenes)

	#Write the output as a bed file for later use
	with open(outputFile, "w") as bedGenes:
		for chr, genelist in chromosomes.items():
			for gene in genelist:
				bedGenes.write("%s\t%i\t%i\t.\t.\t%s\n"%(chr,gene["start"],gene["end"],gene["strand"]))

	log.write("Coverage Breaks: Found %i chromosomes and %i genes total."%(len(chromosomes),geneNumber))
	return chromosomes

#Mixes a list of 3' or 5' results into a list of chromosomes start and end sites as  created by findCoverageBreaks
def insertDashCuts(sourceFile3, sourceFile5, chromosomeList, outputFile, strand="."):
	for dash in ["3","5"]:
		#Read the 3' or 5' coverage file o find clusters of start sites or end sites
		currentDashSite={"chr":None,"position":None,"coverage":None}	#The line read currently
		lastDashSite=None		#The last line that had a coverage over theset threshold
		currentCluster=list()	#The current cluster of results
		currentChr=None			#The current chromosome name
		currentChrList=list()	#The list 3' and 5' sites for the current chromosome
		dashSiteList=dict()		#The full list of relevant 3' or 5' sites
		if dash == "3":
			sourceFile=sourceFile3
		elif dash == "5":
			sourceFile=sourceFile5
		with open(sourceFile, "r") as dashList:
			for line in dashList:
				currentDashSite["chr"],  currentDashSite["position"],  currentDashSite["coverage"] = line.strip().split("\t")
				#First check if we are entering a new chromosome
				if currentDashSite["chr"] != currentChr:
					#if so, store the list of sites for the previous chromosome, then empty the list and reset the variables
					if currentChr is not None and currentCluster != list():
						currentChrList.append(clusterMean(currentCluster))
						dashSiteList[currentChr]=list(currentChrList)
					currentCluster=list()
					currentChr=currentDashSite["chr"]
					currentChrList=list()
					lastDashSite=None

				#lines with a coverage less than a threshold (min_3dash_threshold) are ignored
				if int(currentDashSite["coverage"]) > int(config["min_"+dash+"dash_threshold"]):
					#Peaks within a distance of each other (min_3dash_peak_separation) are clustered together
					if lastDashSite is None or int(currentDashSite["position"])-int(lastDashSite["position"]) < int(config["min_"+dash+"dash_peak_separation"]):
						currentCluster.append(dict(currentDashSite))
				#If the peak is far enough from the last cluster, register the mean position for the peak and reset the cluster with our new result
				elif currentCluster != list():
					#The mean position is stored in the current chromosome list
					currentChrList.append(clusterMean(currentCluster))
					currentCluster=list()
					currentCluster.append(dict(currentDashSite))
					#Once the site has been saved, keep it to compare it to the next
					lastDashSite=dict(currentDashSite)

		#Add the last chromosome list
		if currentChr is not None and currentCluster != list():
			currentChrList.append(clusterMean(currentCluster))
			dashSiteList[currentChr]=list(currentChrList)
		#The dashSiteList now contains every cut site given by our 3' or 5' file, organized by chomosome name

		#Copy the dash site list to the correct variable
		if dash == "3":
			dashSiteList3=dict(dashSiteList)
		elif dash == "5":
			dashSiteList5=dict(dashSiteList)

	#Compare the dashSiteList with our chromosomeList chromosome by chromosome
	newChrList=dict()
	geneNumber=0
	for chr, geneList in chromosomeList.items():
		#Compare each cut site with each gene site
		newGeneList=list()
		for gene in geneList:
			if chr in dashSiteList3.keys():
				for cut3 in dashSiteList3[chr]:
					#If the 3 dash cut is within the start and end of a gene and separated enough from the end of it,
					#cut the gene accordingly and look for a 5' site to start the next gene
					if gene["start"] < cut3 < gene["end"]:
						# if dash == "3" and gene["end"]-cut > int(config["min_3dash_end_separation"]):
						if gene["end"]-cut3 > int(config["min_3dash_end_separation"]):
							#Create a new gene ending at the corresponding 3' position
							newGeneList.append({"start":gene["start"],"end":cut3,"strand":strand})
							geneNumber = geneNumber+1
							gene["start"]=cut3
							#Look for a 5' site to start the next gene
							if chr in dashSiteList5.keys():
								for cut5 in dashSiteList5[chr]:
									#The 5' site must be between a reasonable length from the 3' cut and the end of the gene
									if cut3-int(config["max_5dash_offset"]) < cut5 <gene["end"]:
										#Reduce the current gene to start at the 5' just found
										gene["start"]=cut5
										break

			#The renaining gene i added to the list at the end
			newGeneList.append(gene)
			geneNumber = geneNumber+1
		#The newGeneList is added to our newChrList
		newChrList[chr]=list(newGeneList)

	#Write the output as a bed file
	with open(outputFile, "w") as bedGenes:
		for chr, genelist in newChrList.items():
			for gene in genelist:
				bedGenes.write("%s\t%i\t%i\t.\t.\t%s\n"%(chr,gene["start"],gene["end"],gene["strand"]))

	log.write("Coverage Breaks + Cuts: Found %i chromosomes and %i genes total."%(len(newChrList),geneNumber))
	return newChrList

#function used to calculate the mean position of a cluster of dash sites
def clusterMean(cluster):
	totalSites=0
	totalPosition=0
	#The mean is calculated proportionnaly to the number of sites or each position
	for site in cluster:
		totalSites = totalSites+int(site["coverage"])
		totalPosition = totalPosition+(int(site["position"])*int(site["coverage"]))

	if totalSites == 0:
		return 0
	return int(round(totalPosition/totalSites))

#Combines two chromosomes list together respecting the ordering
def combineChrLists(chrList1, chrList2, outputFile):
	mergedChrList=dict()
	geneNumber=0
	for chr, geneList1 in chrList1.items():
		mergedGeneList=list()
		for gene in geneList1:
			#chrList2[chr][0]["start"] is the start position of the first gene left in chrList2 for the current chromosome
			#While there are genes in chrList2 that start before the current gene, add them to the mergedGeneList and remove them from chrList2
			while len(chrList2[chr]) > 0 and chrList2[chr][0]["start"] < gene["start"]:
				mergedGeneList.append(chrList2[chr].pop(0))
				geneNumber = geneNumber+1
			#When no genes are left that start before the current one, add it to the mergedGeneList
			mergedGeneList.append(gene)
			geneNumber = geneNumber+1
		#When all genes from the current chromosome have been added for chrList1, add any remaining genes from chrList2 to the mergedGeneList
		for gene in chrList2[chr]:
			mergedGeneList.append(gene)
			geneNumber = geneNumber+1

		#Finaly, add the mergedGeneList to the mergedChrList
		mergedChrList[chr]=mergedGeneList

	#The mergedChrList should now contain every gene for every chromosome sorted by chromosome and position

	#Write the output as a bed file
	with open(outputFile, "w") as bedGenes:
		for chr, genelist in mergedChrList.items():
			for gene in genelist:
				bedGenes.write("%s\t%i\t%i\t.\t.\t%s\n"%(chr,gene["start"],gene["end"],gene["strand"]))

	log.write("Combining: Found %i chromosomes and %i genes total."%(len(mergedChrList),geneNumber))

# Find splice sites in a bed file created with the -split option
def findSpliceSites(sourceFile, outputFile):
	#read the file and compare each line with its neighbours
	with open(sourceFile, "r") as bedFile:
		firstline=True
		lastline=None
		exons=[]
		splicedTranscirpts=[]
		for line in bedFile:
			#Seach for duplicated instances of reads in the split file (split sites) by comparing their names
			currentline=line.strip().split("\t")
			if lastline != None and lastline[3] == currentline[3]:
				if firstline:
					exons.append(lastline)
					firstline=False
				exons.append(currentline)
			else:
				#If we found a duplicate, then we add it to the list along with its splice sites
				if not firstline:
					splicedTranscirpts.append(getSplicedTranscript(exons))
				#When a new read name is found, reinitialize the variables
				exons=[]
				firstline=True
			lastline=currentline

		#If the last lines were a duplicate, then we add it to the list along with its splice sites
		if not firstline:
			splicedTranscirpts.append(getSplicedTranscript(exons))

	#sort our duplicate list by key, then by chromosome name for later manipulation
	splicedTranscirpts.sort(key=lambda transcript: transcript[3])
	splicedTranscirpts.sort(key=lambda transcript: transcript[0])

	#Keep only the transcripts with at least 3
	isoforms=[]
	currentIsoform=[]
	firstline=True
	lastSpliceSites=None
	for transcript in splicedTranscirpts:
		spliceSites = " ".join(transcript[3])
		if spliceSites == lastSpliceSites:
			currentIsoform.append(transcript)
		else:
			if firstline:
				firstline=False
			else:
				#Add the isoform to the list if there are at least 3 instances of it
				if len(currentIsoform) > 2:
					isoforms.append(getIsoform(currentIsoform))

			currentIsoform=[transcript]
		lastSpliceSites=spliceSites


	with open(outputFile, "w") as outFile:
		for transcript in isoforms:
			outFile.write("%s\t%i\t%i\t%i\t%s\n"%(transcript[0],transcript[1],transcript[2],transcript[3],"\t".join(transcript[4])))

def getIsoform(transcripts):
	if len(transcripts)<=0:
		return["",0,0,0,[]]
	start=0
	end=0
	for transcript in transcripts:
		#the first part of each transcript should be the same chromosome name
		chr = transcript[0]
		#the final start and end positions will be the mean of all start and end positions, respectively
		start+=int(transcript[1])
		end+=int(transcript[2])
		#the last part of each transcript should be the same splice sites
		spliceSites=transcript[3]

	return [chr,start/len(transcripts),end/len(transcripts),len(transcripts)-1,spliceSites]

def getSplicedTranscript(exons):
	start=None
	lastexon=None
	spliceSites=[]
	for exon in exons:
		#the first part of each read should be the same chromosome name
		chr = exon[0]
		#the end position of the previous exon is the start of a splice site
		if lastexon is not None:
			spliceSites.append(lastexon[2])
		#the start position of the first exon should be the start of the transcript
		if start is None:
			start = exon[1]
		#the start position of the other exons are the end of splice sites
		else:
			spliceSites.append(exon[1])
		lastexon=exon
	#the remaining end position at the end of the loop is the end of the transcript
	end=lastexon[2]
	return [chr,start,end,spliceSites]

#Filters isoforms, removing the ones that are not represented enough compared to the coverage
def filterByCoverage(rawIsoformFile,coverageFile,filteredIsoformFile):
	global config

	coverageMap = {}
	currentChr=None
	#Read the fill coverage file for quick access
	with open(coverageFile,'r') as coverage:
		for line in coverage:
			chr, position, value=line.strip().split("\t")
			if chr != currentChr:
				coverageMap[chr]=[None] #the position 0 does not exist in the coverage file, so we fill it with None
				currentChr=chr
			coverageMap[chr].append(int(value))

	#Check the coverage at every splice position
	with open(rawIsoformFile,'r') as rawIsomophs:
		with open(filteredIsoformFile,'w') as filteredIsoforms:
			for line in rawIsomophs:
				parts = line.strip().split("\t")
				chr=parts[0]
				readCount=int(parts[3])
				spliceSites=parts[4:]
				valid=False
				for spliceSite in spliceSites:
					valid = valid or (readCount != 0 and coverageMap[chr][int(spliceSite)] / readCount < int(config["max_coverage_to_splice_ratio"]))
				if valid:
					filteredIsoforms.write(line)

#Filters isoforms, removing the ones that are not represented enough compared to the coverage
def filterByGenome(rawIsoformFile,genomeFile,filteredIsoformFile):
	global config

	genomeMap = {}
	currentChr=None
	#Read the fill genome file for quick access
	with open(genomeFile,'r') as genome:
		for line in genome:
			if line[0]==">":
				parts=line.strip().split(" ")
				chr=parts[0][1:]
				if chr != currentChr:
					genomeMap[chr]=""
				currentChr=chr
			else:
				genomeMap[currentChr]+=line.strip().upper()

	#Check the genome for poly-A and poly-T content within each exon and delete the exon if it is found
	with open(rawIsoformFile,'r') as rawIsomophs:
		with open(filteredIsoformFile,'w') as filteredIsoforms:
			for line in rawIsomophs:
				parts = line.strip().split("\t")
				chr=parts[0]
				start=parts[1]
				end=parts[2]
				readCount=parts[3]
				spliceSites=parts[4:]
				#Translate the splice sites into exon addresses
				exons=[]
				exon=[start,None]
				for spliceSite in spliceSites:
					if exon is None:
						exon=[spliceSite,None]
					else:
						exon[1]=spliceSite
						exons.append(exon)
						exon=None
				exon[1]=end
				exons.append(exon)

				#Look for poly-A, poly-T or suspiciously high A or T count in the exon
				spliceSites=[]
				start=None
				for	exon in exons:
					sequence=genomeMap[chr][int(exon[0]):int(exon[1])]
					polyA="A" * int(config["max_polyA_length"])
					polyT="T" * int(config["max_polyA_length"])
					if len(sequence) < int(config["min_exon_length"]) and len(sequence) > 0 and (polyA in sequence or polyT in sequence or sequence.count('A')/len(sequence) > float(config["max_single_base_to_length_ratio"]) or sequence.count('T')/len(sequence) > float(config["max_single_base_to_length_ratio"])):
						pass
					else:
						if start is None:
							start=exon[0]
						else:
							spliceSites.append(exon[0])
						spliceSites.append(exon[1])
				if len(spliceSites) > 2:
					end=spliceSites.pop()
					filteredIsoforms.write("\t".join((chr,start,end,readCount,"\t".join(spliceSites)))+"\n")

#Creates a file containing only the most frequent isomoprhs from a file containing a list of them
def findBestIsoform(sourceFile, outputFile):
	global config

	#read the file and compare each line with its neighbours
	with open(sourceFile, "r") as inFile:
		with open(outputFile, "w") as outFile:
			lastline=None
			isoforms=[]
			for line in inFile:
				#Seach for instances sharing approximately the same splice site
				currentline=line.strip().split("\t")
				if lastline is not None and (abs(int(lastline[4]) - int(currentline[4])) < int(config["max_splice_difference"]) or abs(int(lastline[5]) - int(currentline[5])) <  int(config["max_splice_difference"])):
					isoforms.append(currentline)
				else:
					#Once a new transcript is hit, write the best isoform of the list to the output file
					if lastline != None:
						outFile.write("\t".join(getBest(isoforms))+"\n")
					#Then reinitialize with the curent line
					isoforms=[currentline]
				lastline=currentline

			#Write the output for the last batch of isoforms
			outFile.write("\t".join(getBest(isoforms))+"\n")

#Returns the best isoform of a list
def getBest(isoforms):
	best=None
	bestScore=0
	for isoform in (isoforms):
		score=isoform[3]
		if best is None or score > bestScore:
			best=isoform
			bestScore=score
	return best

#Writes a bed file from a file containing isoforms
def convertIsoformstoBed(sourceFile, outputFile):
	with open(sourceFile, "r") as inFile:
		with open(outputFile, "w") as outFile:
			for line in inFile:
				parts=line.strip().split("\t")
				chr=parts[0]
				start=parts[1]
				end=parts[2]
				readCount=parts[3]
				spliceSites=parts[4:]
				#Translate the splice sites into exons
				exons=[]
				exon=[int(start),None]
				for spliceSite in spliceSites:
					if exon is None:
						exon=[int(spliceSite),None]
					else:
						exon[1]=int(spliceSite)
						exons.append(exon)
						exon=None
				exon[1]=int(end)
				exons.append(exon)
				#Create the blocksizes and blockstarts fields from the exons
				blockSizes=""
				blockStarts=""
				for exon in exons:
					blockSizes+="%i,"%(exon[1]-exon[0])
					blockStarts+="%i,"%(exon[0]-int(start))
				outFile.write("%s\t%s\t%s\t\t500\t.\t%s\t%s\t0\t%s\t%s\t%s\n"%(chr,start,end,start,end,len(exons),blockSizes,blockStarts))






#Combines two coverage files for the same chromosomes and adds their coverage together
def combineCoverage(positiveCoverageFile,negativeCoverageFile,allCoverageFile):
	with open(positiveCoverageFile, "r") as positive:
		with open(negativeCoverageFile, "r") as negative:
			with open(allCoverageFile, "w") as all:
				totalCoverage={}
				#Count all coverage in the positive file and organize it in memory
				for posline in positive:
					posparts=posline.split("\t")
					if posparts[0] in totalCoverage.keys():
						if posparts[1] in totalCoverage[posparts[0]].keys():
							totalCoverage[posparts[0]][posparts[1]] += int(posparts[2])
						else:
							totalCoverage[posparts[0]][posparts[1]] = int(posparts[2])
					else:
						totalCoverage[posparts[0]] = {}
						totalCoverage[posparts[0]][posparts[1]] = int(posparts[2])

				#Add up all coverage in the negative file
				for negline in negative:
					negparts=negline.split("\t")
					if negparts[0] in totalCoverage.keys():
						if negparts[1] in totalCoverage[negparts[0]].keys():
							totalCoverage[negparts[0]][negparts[1]] += int(negparts[2])
						else:
							totalCoverage[negparts[0]][negparts[1]] = int(negparts[2])
					else:
						totalCoverage[negparts[0]] = {}
						totalCoverage[negparts[0]][negparts[1]] = int(negparts[2])

				#Write the total coverage in the shared file
				for chr in totalCoverage:
					for position in totalCoverage[chr]:
						all.write("%s\t%s\t%i\n"%(chr,position,totalCoverage[chr][position]))




## Command line functions ##

#These functions all call a tool through command line (subprocess) and wait for it to finish
#Options for these commands are set i the config file (conf.ini)

def minimap(sourceFile, genomeFile, outputFile):
	global config
	global log
	command=config["minimap_path"]+" "+config["minimap_options"]+" "+genomeFile+" "+sourceFile+" > "+outputFile
	log.write("Running minimap with the following command: "+command)
	process = subprocess.Popen(command, shell=True)
	process.wait()
	log.write("Finished running minimap")

def samToBam(sourceFile, outputFile):
	global config
	global log
	command=config["samtools_path"]+" "+config["samtobam_options"]+" "+sourceFile+" > "+outputFile
	log.write("Running samtools with the following command: "+command)
	process = subprocess.Popen(command, shell=True)
	process.wait()
	log.write("Finished running samtools")

def bamToBed(sourceFile, outputFile, split=False):
	global config
	global log
	if split:
		command=config["bedtools_path"]+" "+config["bamtobed_split_options"]+" "+sourceFile+" > "+outputFile
	else:
		command=config["bedtools_path"]+" "+config["bamtobed_options"]+" "+sourceFile+" > "+outputFile
	log.write("Running bedtools with the following command: "+command)
	process = subprocess.Popen(command, shell=True)
	process.wait()
	log.write("Finished running bedtools")

def bedToFasta(sourceFile, genomeFile, outputFile):
	global config
	global log
	command=config["bedtools_path"]+" "+config["getfasta_options"]+" -fi "+genomeFile+" -bed "+sourceFile+" -fo "+outputFile
	log.write("Running bedtools with the following command: "+command)
	process = subprocess.Popen(command, shell=True)
	process.wait()
	log.write("Finished running bedtools")

def sortBam(sourceFile, outputFile):
	global config
	global log
	command=config["samtools_path"]+" "+config["sortbam_options"]+" "+sourceFile+" -o "+outputFile+" > "+outputFile
	log.write("Running samtools with the following command: "+command)
	process = subprocess.Popen(command, shell=True)
	process.wait()
	log.write("Finished running samtools")

def genomecov(sourceFile, genomeFile, outputFile):
	global config
	global log
	command=config["bedtools_path"]+" "+config["genomecov_options"]+" -ibam "+sourceFile+" -g "+genomeFile+" > "+outputFile
	log.write("Running bedtools with the following command: "+command)
	process = subprocess.Popen(command, shell=True)
	process.wait()
	log.write("Finished running bedtools")

def genomecov3dash(sourceFile, outputFile):
	global config
	global log
	command=config["bedtools_path"]+" "+config["genomecov3dash_options"]+" -ibam "+sourceFile+" > "+outputFile
	log.write("Running bedtools with the following command: "+command)
	process = subprocess.Popen(command, shell=True)
	process.wait()
	log.write("Finished running bedtools")

def genomecov5dash(sourceFile, outputFile):
	global config
	global log
	command=config["bedtools_path"]+" "+config["genomecov5dash_options"]+" -ibam "+sourceFile+" > "+outputFile
	log.write("Running bedtools with the following command: "+command)
	process = subprocess.Popen(command, shell=True)
	process.wait()
	log.write("Finished running bedtools")



if __name__ == "__main__":
	main(sys.argv[1:])
