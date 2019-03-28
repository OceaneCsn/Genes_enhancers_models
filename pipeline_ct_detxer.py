import os 
import sys 

if len(sys.argv) == 2:
	ct = sys.argv[1]
	nb_nucl = '2'
	mode = "all"
	pairs_folder = 'EG_pairs'

if len(sys.argv) == 3:
	ct = sys.argv[1]
	nb_nucl = sys.argv[2]
	mode = "all"
	pairs_folder = 'EG_pairs'

if len(sys.argv) == 4:
	ct = sys.argv[1]
	nb_nucl = sys.argv[2]
	mode = sys.argv[3]
	pairs_folder = 'EG_pairs'

if len(sys.argv) == 5:
	ct = sys.argv[1]
	nb_nucl = sys.argv[2]
	mode = sys.argv[3]
	pairs_folder = sys.argv[4]


if mode != "only_compo":


	print("Building a balanced dataset of EP pairs...")

	os.system("Rscript build_pairs_ct.R "+ct)

print("Creating fastas for each pairs...")

os.system("python3 balanced_fasta_cell_type.py "+ct+ " "+pairs_folder)


fasta_enhancers = "Fastas/enhancers_"+ct+".fa"
fasta_genes = "Fastas/genes_"+ct+".fa"

compo_enh = "Compo_nucleo_enhancers/enhancers_"+ct+"_"+nb_nucl+".txt"
compo_proms = "Compo_nucleo_promoters/genes_"+ct+"_"+nb_nucl+".txt"

print("Computing nucleotide composition for enhancers...")


os.system("python3 Quadrinucleotides.py "+ fasta_enhancers + " "+nb_nucl+" > " + compo_enh)

print("Computing nucleotide composition for promoters...")

os.system("python3 Quadrinucleotides.py "+ fasta_genes + " "+nb_nucl+" > " + compo_proms)


if mode != "no_model" and mode != "only_compo":

	print("Creating a LASSO model...")

	os.system("Rscript lasso_roc_ct.R "+ct+ " "+nb_nucl)
	#os.system("Rscript dexter_ct.R "+ct+ " "+nb_nucl)

	print("DONE! See the results in the Figures folder")