import os
import sys

#ct = sys.argv[1]
cts_file = open("cell_types.txt", "r")
cts = [ct.split('\n')[0] for ct in cts_file.readlines()]

for ct in cts[:-1]:


	os.system('Rscript build_pairs_ct.R '+ct+' EG_pairs')

	os.system('python3 fasta_enhancers_gene_expr.py '+ct+' EG_pairs')

	#pos results
	if not os.path.isdir("./Dexter_results/Regression_gene_expr/"+ct+"_pos"):
		os.mkdir("./Dexter_results/Regression_gene_expr/"+ct+"_pos")


	#neg results
	if not os.path.isdir("./Dexter_results/Regression_gene_expr/"+ct+"_neg"):
		os.mkdir("./Dexter_results/Regression_gene_expr/"+ct+"_neg")


	#neg indexes
	if not os.path.isdir("./Indexes/Regression_gene_expr/"+ct+"_neg"):
		os.mkdir("./Indexes/Regression_gene_expr/"+ct+"_neg")


	#pos indices
	if not os.path.isdir("./Indexes/Regression_gene_expr/"+ct+"_pos"):
		os.mkdir("./Indexes/Regression_gene_expr/"+ct+"_pos")


	#pos indexing and domains search
	os.system("motifIndex index ./Indexes/Regression_gene_expr/index_"+ct+"_pos ./Reg_gene_expr_dexter/enhancers_pos"+ct+".fa")
	os.system("python3 ~/DExTER/Main.py -fasta_file ./Reg_gene_expr_dexter/enhancers_pos"+ct+".fa -index_file ./Indexes/Regression_gene_expr/index_"+ct+"_pos.ssa -expression_file ./Reg_gene_expr_dexter/expression_pos_"+ct+".csv -target_condition col1 -experience_directory ./Dexter_results/Regression_gene_expr/"+ct+"_pos -alignement_point_index 0 -nb_regions 10")

	#neg indexing and domain search
	os.system("motifIndex index ./Indexes/Regression_gene_expr/index_"+ct+"_neg ./Reg_gene_expr_dexter/enhancers_neg"+ct+".fa")
	os.system("python3 ~/DExTER/Main.py -fasta_file ./Reg_gene_expr_dexter/enhancers_neg"+ct+".fa -index_file ./Indexes/Regression_gene_expr/index_"+ct+"_neg.ssa -expression_file ./Reg_gene_expr_dexter/expression_neg_"+ct+".csv -target_condition col1 -experience_directory ./Dexter_results/Regression_gene_expr/"+ct+"_neg -alignement_point_index 0 -nb_regions 10")


#os.system("Rscript Dexter_Lasso_comparison.R")