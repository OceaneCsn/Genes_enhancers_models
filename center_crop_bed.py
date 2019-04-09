import sys

bed = sys.argv[1]

enhancers = open(bed, "r")


for i, line in enumerate(enhancers.readlines()):
	l = line.split('\t')
	center = int((int(l[1])+ int(l[2]))/2)
	print(l[0]+'\t'+str(center-501)+'\t'+str(center+500)+'\t'+l[3].split('\n')[0])