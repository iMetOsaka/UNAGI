# General

UNAGI (UNAnnotated Gene Identifier) is a pipeline allowing for the detection and annotation of transcripts from long-read sequences (usually from Oxford Nanopore Technologies Sequencers) with the help of a genome file from the studied organism.

It will typically be used as follows
On Linux:
```
$ ./unagi -i [path_to_fastq_file] -o [path_to_desired_output] -g [path_to_genome_file]
```

On Windows:
```
> unagi.bat -i [path_to_fastq_file] -o [path_to_desired_output] -g [path_to_genome_file]
```

The results are classified between transitional files in their own folder and final files.
Final files are the main results. They are presented in the bed format.
Transitional files are the files generated throughout the pipeline run. They can be in different formats and will be discussed in details in the "Detailed Configuration" part of this document.

# Installation

UNAGI requires python and BioPython to run. It was tested using python3 but should run with python 2.7 as well.

## Test if python is present on your system
You can check if python is already on your system by typing the following command in your command line
```
$ python3 --version
```
If you get a version number, you can skip to installing BioPython
Otherwise you will have to install python using one of the following methods depending on your OS.

### On Ubuntu
install python through apt-get:
```
$ sudo apt-get update
$ sudo apt-get install python3.6
```
### On CentOS/RHEL
install python through yum:
```
$ sudo yum upgrade
$ sudo yum install python3
```
### On Windows
Download the installer from [the official python download page](https://www.python.org/downloads/windows/) and run it.
Be sure to check the option "Add python.exe to Path" during the installation.

## Install BioPython

Once python is installed, the BioPython module is easily installed using pip from the command line:
```
$ pip install biopython
```

## Check UNAGI

UNAGI should now be ready to run. Check it by going to the unagi directory and type the following command:
```
cd /path/to/unagi/directory
./unagi -v
```

Alternatively, if you want to run a specific version of python, you can use:
```
cd /path/to/unagi/directory
/my/preferred/python app/unagi.py -v
```

# Options

**Required**
```
-i  *or* --input
```
Input file. This should be the path to a fastq file.

```
-o  *or* --output
```
Output path. This should be the path to output the result files to.

```
-g  *or* --genome
```
Genome File. This should be a valid genome file for minimap.

**Optional**
```
-v  *or* --version
```
Displays the program version and exits.

```
-s *or* --stranded
```
Stranded. Should be selected if the input file contains reads that are already stranded.
The input file should then be the fastq file containing the stranded reads

```
-V *or* --verbose
```
Verbose. Displays additional information while running.

```
-S  *or* --silent
```
Runs silently (no console output), cancels verbose.

# Detailed Configuration

More granular configuration can be made through the conf.ini file located in the app folder.

###Path to Utilities
minimap_path=./tools/minimap2/minimap2
samtools_path=./tools/samtools/samtools
bedtools_path=./tools/bedtools/bin/bedtools

###Research options
strand_forwards_identifier=GGG
strand_backwards_identifier=TTT
strand_start_index=24
strand_end_index=29
min_coverage_threshold=6
min_3dash_threshold=10
min_3dash_peak_separation=30
min_3dash_end_separation=100
min_5dash_threshold=10
min_5dash_peak_separation=30
min_5dash_start_separation=100
max_5dash_offset=200
max_coverage_to_splice_ratio=100
min_exon_length=100
max_polyA_length=8
max_single_base_to_length_ratio=0.7
max_splice_difference=10

###Command options
minimap_options=-ax splice -G 5k --secondary=no
samtobam_options=view -Sb
bamtobed_options=bamtobed -i
bamtobed_split_options=bamtobed -split -i
getfasta_options=getfasta
sortbam_options=sort
genomecov_options=genomecov -d
genomecov3dash_options=genomecov -dz -3
genomecov5dash_options=genomecov -dz -5


###Transitional output names
stranded_file=stranded.fastq
raw_mapped_sam_file=raw_mapped.sam
raw_mapped_bam_file=raw_mapped.bam
raw_sorted_bam_file=raw_sorted.bam
raw_sorted_bed_file=raw_sorted.bed
raw_sorted_split_bed_file=raw_sorted_split.bed
positive_bed_file=positive.bed
negative_bed_file=negative.bed
positive_refined_fasta_file=positive.fasta
negative_refined_fasta_file=negative.fasta
positive_mapped_sam_file=positive_mapped.sam
negative_mapped_sam_file=negative_mapped.sam
positive_mapped_bam_file=positive_mapped.bam
negative_mapped_bam_file=negative_mapped.bam
positive_sorted_bam_file=positive_sorted.bam
negative_sorted_bam_file=negative_sorted.bam
positive_coverage_file=positive_coverage.txt
negative_coverage_file=negative_coverage.txt
total_coverage_file=total_coverage.txt
positive_genelist_incomplete_file=positive_genelist_incomplete.bed
negative_genelist_incomplete_file=negative_genelist_incomplete.bed
positive_3dash_file=positive_3dash_positions.txt
negative_3dash_file=negative_3dash_positions.txt
positive_5dash_file=positive_5dash_positions.txt
negative_5dash_file=negative_5dash_positions.txt
positive_genelist_file=positive_genelist.bed
negative_genelist_file=negative_genelist.bed
raw_splice_sites_file=raw_splice_sites.txt
coverage_filtered_splice_sites_file=coverage_filtered_splice_sites.txt
genome_filtered_splice_sites_file=genome_filtered_splice_sites.txt
best_filtered_splice_sites_file=best_filtered_splice_sites.txt

###Final output names
full_genelist_file=Full_Genelist.bed
final_splice_sites_file=Splicing_Isoforms.bed
final_unique_splice_sites_file=Splicing_Isoforms_Unique.bed
