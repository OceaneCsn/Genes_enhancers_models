import os


#os.system("Rscript build_pairs_ct.R none plot")


cts = []
files = []
for r,d,f in os.walk("EG_pairs"):
	files = f
	for file in f:
		cts.append(file.split('.')[0].split('_')[-1])

cts_file = open("cell_types_plt.txt", "w")

for ct in cts:
	print("Compositions (quadri and dexter) for cell type : ",ct)
	cts_file.write(ct+'\t')
	os.system("python3 pipeline_ct.py "+ct+" 2 only_compo EG_pairs")
	os.system("python3 pipeline_ct_dexter.py "+ct+" 2 only_compo EG_pairs")

#os.system("Rscript Plot_auc_length.R")
