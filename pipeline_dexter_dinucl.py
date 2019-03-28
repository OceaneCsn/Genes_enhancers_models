import os

cts = []
files = []
for r,d,f in os.walk("Length_AUC_plot"):
	files = f
	for file in f:
		cts.append(file.split('.')[0].split('_')[-1])



for ct in cts:
	print("####### cell type : ",ct)
	os.system("python3 pipeline_ct_dexter.py "+ct+" 2 only_compo EG_pairs")


#os.system("Rscript Plot_dexter_dinucl.R")
