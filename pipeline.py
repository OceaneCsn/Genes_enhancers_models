import os 
import sys 

ct = sys.argv[1]



print("Building a balanced dataset of EP pairs...")

os.system("Rscript build_pairs_ct.R "+ct+" 00")

print("Creating fastas for each pairs...")

os.system("python3 balanced_fasta_cell_type.py "+ct+ " EG_pairs 00")


fasta_enhancers = "Fastas/enhancers_"+ct+".fa"
fasta_genes = "Fastas/genes_"+ct+"_"+".fa"

compo_enh = "Compo_nucleo_enhancers/enhancers_"+ct+"_"+nb_nucl+".txt"
compo_proms = "Compo_nucleo_promoters/genes_"+ct+"_"+nb_nucl+".xt"

print("Computing nucleotide composition for enhancers...")

os.system("python3 Quadrinucleotides.py "+ fasta_enhancers + " "+nb_nucl+" > " + compo_enh)

print("Computing nucleotide composition for promoters...")

os.system("python3 Quadrinucleotides.py "+ fasta_genes + " "+nb_nucl+" > " + compo_proms)

os.system("python3 dexter_ct.py "+ct)