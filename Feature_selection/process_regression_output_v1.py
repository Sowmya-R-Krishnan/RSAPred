#Program to process regression output and extract the best models

import sys
import csv
import pandas as pd

data = sys.argv[1]
pass_models = sys.argv[2]
pass_stringent = sys.argv[3]

rna_features = ["A", "G", "C", "U", "AA", "AG", "AC", "AU", "GA", "GG", "GC", "GU", "CA", "CG", "CC", "CU", "UA", "UG", "UC", "UU", "AAA", "AAG", "AAC", "AAU", "AGA", "AGC", "AGU", "ACA", "ACC", "ACU", "AUG", "AUC", "AUU", "GAA", "GAG", "GAC", "GAU", "GGA", "GGG", "GGC", "GGU", "GCA", "GCG", "GCC", "GCU", "GUA", "GUG", "GUC", "GUU", "CAA", "CAG", "CAC", "CAU", "CGA", "CGG", "CGC", "CGU", "CCA", "CCG", "CCC", "CCU", "CUA", "CUG", "CUC", "UAC", "UAU", "UGA", "UGG", "UGC", "UGU", "UCA", "UCG", "UCC", "UCU", "UUG", "UUC", "UUU", "AAAG", "AAAU", "AAGU", "AAUG", "AAUC", "AGCA", "AGCG", "AGUG", "ACAC", "ACAU", "ACCC", "ACUA", "ACUG", "AUGG", "AUGU", "AUCU", "AUUC", "AUUU", "GAAA", "GAAU", "GAGC", "GACA", "GACU", "GAUG", "GGAA", "GGAG", "GGGA", "GGGC", "GGGU", "GGCA", "GGCC", "GGUU", "GCAA", "GCAC", "GCGA", "GCGU", "GCUG", "GUGG", "GUGC", "GUGU", "GUCA", "GUCC", "GUUG", "CAAA", "CAAG", "CAAC", "CAGU", "CACC", "CACU", "CAUG", "CAUU", "CGAC", "CGCU", "CCAA", "CCAG", "CUGA", "CUGC", "CUGU", "CUCA", "UAUU", "UGAA", "UGAG", "UGAC", "UGAU", "UGGA", "UGGG", "UGGC", "UGCG", "UGCU", "UGUG", "UGUC", "UGUU", "UCAA", "UCAC", "UCAU", "UCCA", "UCUG", "UCUC", "UUGA", "UUGG", "UUCU", "UUUG", "DNC_AA", "DNC_AG", "DNC_AC", "DNC_AU", "DNC_GA", "DNC_GG", "DNC_GC", "DNC_GU", "DNC_CA", "DNC_CG", "DNC_CC", "DNC_CU", "DNC_UA", "DNC_UG", "DNC_UC", "DNC_UU", "DNC_Feat_17", "DNC_Feat_18", "DNC_Feat_19", "DNC_Feat_20", "DNC_Feat_21", "DNC_Feat_22", "DNC_Feat_23", "DNC_Feat_24", "A(((", "A((.", "A(..", "A(.(", "A.((", "A..(", "A...", "G(((", "G((.", "G(..", "G(.(", "G.((", "G..(", "G...", "C(((", "C((.", "C(.(", "C.((", "C.(.", "C..(", "C...", "U(((", "U((.", "U(..", "U(.(", "U.((", "U..(", "U...", "A,A", "A,C", "A,U", "A,A-U", "A,U-A", "A,C-G", "G,A", "G,G", "G,U", "G,A-U", "G,G-C", "G,C-G", "C,A", "C,C", "C,U", "C,U-A", "C,G-C", "C,C-G", "U,A", "U,G", "U,C", "U,U", "U,G-C", "U,G-U", "A-U,A", "A-U,C", "A-U,A-U", "A-U,U-A", "A-U,G-C", "A-U,C-G", "A-U,U-G", "U-A,A", "U-A,G", "U-A,U", "U-A,U-A", "U-A,G-CU-A,C-G", "U-A,G-U", "U-A,U-G", "G-C,A", "G-C,C", "G-C,U", "G-C,A-U", "G-C,U-A", "G-C,G-C", "G-C,C-GG-C,G-U", "G-C,U-G", "C-G,A", "C-G,G", "C-G,A-U", "C-G,U-A", "C-G,G-C", "C-G,C-G", "C-G,U-G", "G-U,A-UG-U,U-A", "G-U,G-C", "G-U,C-G", "G-U,G-U", "U-G,G", "U-G,A-U", "U-G,C-G", "U-G,G-U", "U-G,U-G", "PS_Feat_101", "PS_Feat_102", "PS_Feat_103", "PS_Feat_104", "PS_Feat_105", "PS_Feat_106", "PS_Feat_107", "PS_Feat_108"]

p1 = []
p2 = []
with open(data) as output:
	for line in output.readlines():
		line = line.strip()
		contents = line.split('\t')
		corr_coeff = float(contents[-1])

		if(corr_coeff>=0.60):
			rna_check = 0
			p1.append(contents)
			for feat in contents[:-1]:
				if(feat in rna_features):
					rna_check = rna_check + 1
			if(rna_check>0 and rna_check<=(len(contents)-2)):
				p2.append(contents)

print(len(p1), len(p2))

with open(pass_stringent, 'w') as f2:
	for model in p2:
		print('\t'.join(model), file=f2)

print("Regression output processed.")

