
from Bio import SeqIO
input_file = "trimmed_barcode01.fastq"
forward = (rec for rec in  SeqIO.parse(input_file, "fastq")  if 'GGG' in rec.seq[24:29])
count = SeqIO.write(forward, "f.fastq", "fastq")
print("Saved %i reads" % count)
complement = (rec for rec in  SeqIO.parse(input_file, "fastq")  if 'TTT' in rec.seq[24:29])
count = SeqIO.write(complement, "complement.fastq", "fastq")
print("Saved %i reads" % count)
#implement seqtk alternative same package
#implement cat for python
#minimap > samtools samtobam > samtools bamsorted > bamtobed
#seperate + and -
#bed to fasta
#minimap > samtools samtobam > samtools bamsorted
#bedtools genomecov > awk 1&3 print
#python coverage_v1.py coverage.txt (covert bed) >
#bedtools genomcov -3 cov > awk '$3 > 10 { print }' _58pl_3.bed > _58pl_3_10.bed
#
