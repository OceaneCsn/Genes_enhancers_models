import sys
from Bio import SeqIO

coding_tss = open("Coding_tss.txt", 'r')
tsss = coding_tss.readlines()
tsss = [tss.split('\n')[0] for tss in tsss]
coding_tss.close()

genes = {}

# stores the sequence only if it corresponds to a most expressed coding tss
for record in SeqIO.parse("fantom_TSS.fa", "fasta"):
	if record.id in tsss:
		genes[record.id] = record.seq
		

for ID in list(genes.keys()):
	print('>'+ ID +'\n'+str(genes[ID])+'\n')