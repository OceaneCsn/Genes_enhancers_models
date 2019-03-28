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

print(mode)

if mode != "only_compo":


	print("Building a balanced dataset of EP pairs...")

	os.system("Rscript build_pairs_ct.R "+ct)

	print("Creating fastas for each pairs...")

	os.system("python3 balanced_fasta_cell_type.py "+ct+ " "+pairs_folder)


os.system("python3 dexter_ct.py "+ct)

if mode != "no_model" and mode != "only_compo":

	print("Creating a LASSO model...")

	os.system("Rscript dexter_ct.R "+ct+ " "+nb_nucl)

	print("DONE! See the results in the Figures folder")