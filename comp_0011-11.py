import os

cts = ["CNhs11252", "CNhs11322","CNhs10870", "CNhs11872", "CNhs13637", "CNhs10882", "CNhs11761", "CNhs13195", "CNhs10727","CNhs10723", "CNhs14045", "CNhs12845", "CNhs12516", "CNhs13479", "CNhs12806", "CNhs10854", "CNhs11934", "CNhs10731","CNhs12191" ,"CNhs11253", "CNhs10853", "CNhs11790" ,"CNhs12321" ,"CNhs11083", "CNhs12311" ,"CNhs11732", "CNhs11247","CNhs11945" ,"CNhs10868", "CNhs12839", "CNhs12075", "CNhs12368" ,"CNhs10728", "CNhs10855", "CNhs11882" ,"CNhs11864"]

for ct in cts:
	print("python3 pipeline_ct.py "+ct+" 2 all EG_pairs && python3 pipeline_ct.py "+ct+" 2 all EG_pairs_00 00")
	os.system("python3 pipeline_ct.py "+ct+" 2 all EG_pairs && python3 pipeline_ct.py "+ct+" 2 all EG_pairs_00 00")
	os.system("Rscript lasso_1100_11.R "+ct)