import pandas as pd
import sys
from Bio import SeqIO
cell_type = sys.argv[1]
pairs_folder = sys.argv[2]

####### writes a fasta for all the sequences used for the dataset of one cell type
# uses the file generated by EG_pairs_celltype.py for the variable pairs


enhancers = {}

#dataset with enhancers and gene ids

print('EG_pairs_'+cell_type+'.txt')
pairs = pd.read_csv(pairs_folder+'/EG_pairs_'+cell_type+'.txt', sep = ';')

#fasta generated
fasta_enhencers = open('Reg_gene_expr_dexter/enhancers_'+cell_type+'.fa', "w")
fasta_enhencers_pos = open('Reg_gene_expr_dexter/enhancers_pos'+cell_type+'.fa', "w")
fasta_enhencers_neg = open('Reg_gene_expr_dexter/enhancers_neg'+cell_type+'.fa', "w")

print("processing enhancers")
for record in SeqIO.parse("Data/fantom_enhancers.fa", "fasta"):
	enhancers[record.id] = record.seq

print("creating enhancers fasta")


gene_expr_pos = open('Reg_gene_expr_dexter/'+'expression_pos_'+cell_type+'.csv', "w")
gene_expr_neg = open('Reg_gene_expr_dexter/'+'expression_neg_'+cell_type+'.csv', "w")


gene_expr_pos.write("gene"+'\t'+"col1"+'\n')
gene_expr_neg.write("gene"+'\t'+"col1"+'\n')


for i in range(len(pairs)):
	e = pairs["Enhancer"].iloc[i]
	p = pairs["Gene"].iloc[i]

	fasta_enhencers.write('>'+ '/'.join([e, p]) +'\n'+str(enhancers[e])+'\n')

	if int(pairs["Interaction"].iloc[i])== 1:
		gene_expr_pos.write('/'.join([e, p])+'\t'+str(pairs["gene_expression"].iloc[i])+'\n')
		fasta_enhencers_pos.write('>'+ '/'.join([e, p]) +'\n'+str(enhancers[e])+'\n')
	else:
		gene_expr_neg.write('/'.join([e, p])+'\t'+str(pairs["gene_expression"].iloc[i])+'\n')
		fasta_enhencers_neg.write('>'+ '/'.join([e, p]) +'\n'+str(enhancers[e])+'\n')

fasta_enhencers.close()
gene_expr_pos.close()
gene_expr_neg.close()