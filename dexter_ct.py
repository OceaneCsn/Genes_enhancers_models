import os 
import sys

ct = sys.argv[1]

#creates the directories if needed
if not os.path.isdir("./Dexter_results/Enhancers"+ct):
	os.mkdir("./Dexter_results/Enhancers/"+ct)

if not os.path.isdir("./Indexes/Enhancers/"+ct):
	os.mkdir("./Indexes/Enhancers/"+ct)

if not os.path.isdir("./Dexter_results/Promoters/"+ct):
	os.mkdir("./Dexter_results/Promoters/"+ct)

if not os.path.isdir("./Indexes/Promoters/"+ct):
	os.mkdir("./Indexes/Promoters/"+ct)


# Indexing for enhancers
os.system("motifIndex index ./Indexes/Enhancers/index_"+ct+" ./Fastas/enhancers_"+ct+".fa")

# Dexter variables for enhancers
os.system("python3 ~/DExTER/Main.py -fasta_file ./Fastas/enhancers_"+ct+".fa -logistic -index_file ./Indexes/Enhancers/index_"+ct+".ssa -expression_file ./Classif_dexter/Enhancers/classif_"+ct+".csv -target_condition col1 -experience_directory ./Dexter_results/Enhancers/"+ct+" -alignement_point_index 0 -nb_regions 10")

# Indexing for promoters
os.system("motifIndex index ./Indexes/Promoters/index_"+ct+" ./Fastas/genes_"+ct+".fa")

# Dexter variables for promoters
os.system("python3 ~/DExTER/Main.py -fasta_file ./Fastas/genes_"+ct+".fa -logistic -index_file ./Indexes/Promoters/index_"+ct+".ssa -expression_file ./Classif_dexter/Promoters/classif_"+ct+".csv -target_condition col1 -experience_directory ./Dexter_results/Promoters/"+ct+" -alignement_point_index 0 -nb_regions 10")