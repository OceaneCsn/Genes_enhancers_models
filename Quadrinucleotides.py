import pandas as pd
import sys
import re
fastafile = sys.argv[1]
#nombre de nucleotides a prendre en compte dans le calcul des frequences
nbnucl = int(sys.argv[2])
                
fasta = open(fastafile, 'r')
seq = fasta.readlines()
dinuc = ["GpC", "TpC", "CpA", "CpT", "ApC", "ApT", "TpA", "CpG", "ApA", "CpC"]
trinuc = ['TTT', 'AAA', 'ATT', 'AAT', 'TCT', 'AGA', 'TTA', 'TAA', 'CTT', 'AAG', 'TTC', 'GAA', 'TAT', 'ATA', 'TGT', 'CTG', 'ACA', 'CAG', 'TGA', 'TCA', 'CAT', 'ATG', 'TTG', 'CAA', 'CCT', 'AGG', 'TGG', 'CCA', 'AGT', 'ACT', 'GAG', 'CTC', 'GGA', 'TCC', 'GTT', 'AAC', 'TGC', 'GCA', 'GTG', 'CAC', 'GCT', 'AGC', 'GGG', 'CCC', 'GAT', 'ATC', 'TAG', 'CTA', 'GTA', 'TAC', 'GCC', 'GGC', 'GGT', 'ACC', 'GTC', 'GAC', 'CGG', 'CCG', 'CGT', 'ACG', 'GCG', 'CGC', 'TCG', 'CGA']
quadri = list()

for tri in trinuc:
        for n in ['A', 'C', 'G', 'T']:
                quadri.append(tri+n)

#retourne la sequence complementaire du parametre
def compl(seq):
        nuc = list(seq)
        comp = []
        cn = ''
        for n in nuc:
                if(n == 'A'):
                        cn = 'T'
                if(n == 'T'):
                        cn = 'A'
                if(n == 'G'):
                        cn = 'C'
                if(n == 'C'):
                        cn = 'G'
                comp.append(cn)
        return ''.join(comp)

def reverseCompl(mot):
        return ''.join(list(reversed(list(compl(mot)))))

#pour les mots pairs, retourne les mots qui sont leur reverse compl
def soucis(seq):
        soucis = []
        for mot in seq:
                if(mot == reverseCompl(mot)):
                        soucis.append(mot)
        return soucis


soucisq = soucis(quadri)
#pas de probleme pour les trinuc
splitdinuc = [''.join(d.split('p')) for d in dinuc]
soucisd = soucis(splitdinuc)

#on supprime les inverses des complementaires des sequences : aura 32 variables trinucleotides plutot que 64
usefultrinuc = []
uselesstrinuc = []
for t in trinuc:
        if(t not in uselesstrinuc):
                usefultrinuc.append(t)
                uselesstrinuc.append(reverseCompl(t))

#pareil pour les quadriusefulq = []
uselessq = []
usefulq = []
for q in quadri:
        if(q not in uselessq):
                usefulq.append(q)
                if(q not in soucisq):
                        uselessq.append(reverseCompl(q))

#remplissage du dictionnaire des frequences
data = {}
for i in range(0,len(seq)):
        if seq[i].startswith('>'):
                if i!= 0:
                        if nbnucl>=1:
                        #nucleotides
                                d["AT"] = (s.count('A') + s.count('T'))/float(len(s))
                                d["GC"] = (s.count('G') + s.count('C'))/float(len(s))
                        
                        if(nbnucl>=2):
                        #dinucleotides
                                for din in dinuc:
                                        dn = din.split('p')
                                        dn = ''.join(dn)
                                        if(dn not in soucisd):
                                                add = s.count(reverseCompl(dn))
                                        else:
                                                add = 0
                                        d[din] = (s.count(dn)+add) /float(len(s))
                                        
                        #trinucleotides
                        if(nbnucl>=3):
                                for tn in usefultrinuc:
                                        d[tn] = (s.count(tn)+s.count(reverseCompl(tn))) /float(len(s))
                                
                        if(nbnucl>=4):
                        #quadrinucleotides
                                for q in usefulq:
                                        add=0
                                        if(q not in soucisq):
                                                add = s.count(reverseCompl(q))
                                        d[q] = (s.count(q)+add) /float(len(s))
                #nom de la sequence comme cle du dictionnaire
                data[seq[i][1:-1]] = {}
                d = data[seq[i][1:-1]]
                s = ''
        else:
                #calcul des taux
                s += seq[i].upper()[0:-1]

                
                                
df = pd.DataFrame.from_dict(data, orient='columns')
df.to_csv(sys.stdout, sep = ';')