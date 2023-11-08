#Calculate parallel pseudo-dinucleotide composition features for RNA sequence - based on repRNA method

import sys
import math
import pandas as pd
import itertools
import numpy as np
from collections import Counter

data = sys.argv[1]
prop_data = sys.argv[2]
outpath = sys.argv[3]

df = pd.read_csv(data, sep='\t', header=0)
prop_df = pd.read_csv(prop_data, sep='\t', header=0)

seqs = list(df['Target_RNA_sequence'])
names = list(df['Target_RNA_name'])
ids = list(df['Target_RNA_ID'])

def calc_correlation_vector(seq, lambda_val, prop_df):
	mu_val = len(list(prop_df.columns))-1
	
	corr_vector = []
	seq_len = len(seq)
	
	#TODO: In case of confusion, refer repRNA paper for details
	for i in range(1, lambda_val+1):
		corr_sum = 0

		for j in range(1, seq_len-(i+1)):
			#TODO: If confused, to avoid index errors use try-catch block
			dinucl1 = seq[j-1]+seq[j+1-1]  #Adding -1 to all indices to account for string index starting from 0
			dinucl2 = seq[j+i-1]+seq[j+i+1-1]  #Adding -1 to all indices to account for string index starting from 0

			prop_vec1 = prop_df[prop_df['Dinucleotide']==dinucl1].values[0]
			prop_vec2 = prop_df[prop_df['Dinucleotide']==dinucl2].values[0]
				
			corr_val = 0
			for u in range(1, mu_val+1):  #Starts from 1 to avoid the first column (dinucleotide type)
				corr_val = corr_val + (float(prop_vec1[u]) - float(prop_vec2[u]))**2

			corr_sum = corr_sum + corr_val/mu_val

		corr_vector.append(corr_sum/(seq_len-(i+1)))  #TODO: Never forget to bracket the subtraction term in denominator!

	return corr_vector

#Based on http://www.hgmd.cf.ac.uk/docs/nuc_lett.html
seq_repl = {"X":"A", "N":"A", " ":"", "R":"G", "Y":"C", "K":"G", "M":"A", "S":"G", "W":"A", "B":"G", "D":"G", "H":"A", "V":"G"}

seqs_new = []
names_new = []
ids_new = []
for i, seq in enumerate(seqs):
	try:
		seq = seq.replace("T", "U")
		seq_chars = list(seq)
		chars = Counter(seq_chars)
		for char in chars.keys():
			if(char not in ["A", "G", "C", "U"]):
				seq = seq.replace(char, seq_repl[char])

		seqs_new.append(seq)
		names_new.append(names[i])
		ids_new.append(ids[i])
	except:
		#print(names[i], ids[i])
		continue

bases_prime = ["A", "G", "C", "U"]
di_strings = [p for p in itertools.product(bases_prime, repeat=2)]

for i, di_str in enumerate(di_strings):
	di_strings[i] = ''.join(di_str)

#Fixed parameters: lambda_val < minimum_query_length (depends on the database of sequences used as input)
lambda_val = 8
w_val = 0.5

with open(outpath+"pPseudoDNC_features_v1.out", 'w') as f:
	header_str = "Target_RNA_ID\tTarget_RNA_name\t"+"\tDNC_".join(di_strings)+"\t"
	for i in range(17, 17+lambda_val):
		header_str = header_str+"DNC_Feat_"+str(i)+"\t"

	print(header_str, file=f)

	for i, seq in enumerate(seqs_new):
		dinucleotides = []
		feature_vec = [0]*(16+lambda_val)  #16 = No. of possible dinucleotides
		dinucl_norm_freq = []

		corr_vector = calc_correlation_vector(seq, lambda_val, prop_df)

		for j in range(len(seq)):
			if(j+1<len(seq)):
				dinucleotides.append(seq[j]+seq[j+1])

		dinucl_freq = Counter(dinucleotides)
		for s, di_str in enumerate(di_strings):
			dinucl_norm_freq.append(dinucl_freq[di_str]/len(seq))

		for s, di_str in enumerate(di_strings):
			feature_vec[s] = np.round(dinucl_norm_freq[s]/(sum(dinucl_norm_freq)+(w_val*sum(corr_vector))), 5)
		
		for s in range(16, 16+lambda_val):
			feature_vec[s] = np.round((w_val*corr_vector[s-16])/(sum(dinucl_norm_freq)+(w_val*sum(corr_vector))), 5)
			
		feat_str = ""
		feat_str = feat_str+ids_new[i]+"\t"+names_new[i]+"\t"

		for feat in feature_vec:
			feat_str = feat_str+str(feat)+"\t"

		print(feat_str, file=f)

print("pPseudoDNC features calculated successfully")
































