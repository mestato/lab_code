#!/lustre/projects/staton/software/ActivePython-2.7.8/bin/python

from Bio import SeqIO
handle = open("test.fasta", "rU")
for record in SeqIO.parse(handle, "fasta") :
    print record.id
handle.close()
