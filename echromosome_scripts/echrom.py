#!/lustre/projects/staton/software/ActivePython-2.7.8/bin/python

from Bio import SeqIO
import sys

# intialize variables
cellP = 0

# load reference record
seq_record = SeqIO.index("Ppersica_139_v1.hardmasked.fa", "fasta")

# load and cycle through layout file for echromosomes
for line in open('echromsLayout.tsv', 'r'):
	cell = line.split("\t")

	# first pass set chromosomes
	if cellP == 0:
		cellP = cell[0]
		sys.stdout.write(">echromosome_")
		print(cellP)

	# did we hit a new chromosome?
	if cell[0] != cellP:
		cellP = cell[0]
		print('')
		sys.stdout.write(">echromosome_")
		print(cellP)
	chrom  = cell[2]
	start  = int(cell[3])
	end    = int(cell[4])
	orient = int(cell[5])

	# check orientation and print if 0 forward all else reverse_complement
	if orient == 0: 
		sys.stdout.write(str(seq_record[chrom].seq[start:end]))
	else:
		sys.stdout.write(str(seq_record[chrom].seq[start:end].reverse_complement()))

# print newline
print('')
