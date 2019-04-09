import os
import sys

ct = sys.argv[1]


os.system('Rscript build_pairs_ct.R '+ct+' EG_pairs')

os.system('python3 fasta_enhancers_gene_expr.py '+ct+' EG_pairs')


if not os.path.isdir("./Dexter_results/Regression_gene_expr/"+ct):
	os.mkdir("./Dexter_results/Regression_gene_expr/"+ct)

if not os.path.isdir("./Indexes/Regression_gene_expr/"+ct):
	os.mkdir("./Indexes/Regression_gene_expr/"+ct)

os.system("motifIndex index ./Indexes/Regression_gene_expr/index_"+ct+" ./Reg_gene_expr_dexter/enhancers_"+ct+".fa")


os.system("python3 ~/DExTER/Main.py -fasta_file ./Reg_gene_expr_dexter/enhancers_"+ct+".fa -index_file ./Indexes/Regression_gene_expr/index_"+ct+".ssa -expression_file ./Reg_gene_expr_dexter/expression_"+ct+".csv -target_condition col1 -experience_directory ./Dexter_results/Regression_gene_expr/"+ct+" -alignement_point_index 0 -nb_regions 10")


os.system("Rscript plot_hist_error_reg_ct.R "+ct)