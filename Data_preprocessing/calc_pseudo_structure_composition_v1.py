#Calculate pseudo structure status composition features for RNA sequence - based on repRNA method

import sys
import math
import RNA
import pandas as pd
import itertools
import numpy as np
from collections import Counter

data = sys.argv[1]
outpath = sys.argv[2]

df = pd.read_csv(data, sep='\t', header=0)

seqs = list(df['Target_RNA_sequence'])
names = list(df['Target_RNA_name'])
ids = list(df['Target_RNA_ID'])

def calc_correlation_vector(seq, comp_seq, lambda_val, free_energy):	
	corr_vector = []
	seq_len = len(seq)
	
	#TODO: In case of confusion, refer repRNA paper for details
	for i in range(1, lambda_val+1):
		corr_sum = 0

		for j in range(1, seq_len-i+1):
			#TODO: If confused, to avoid index errors use try-catch block
			bp1 = comp_seq[j-1]  #Adding -1 to the indices to account for string index starting from 0 in python
			bp2 = comp_seq[j+i-1]  #Adding -1 to the indices to account for string index starting from 0 in python

			f1 = free_energy[bp1]
			f2 = free_energy[bp2]
				
			corr_val = (f1-f2)**2
			corr_sum = corr_sum + corr_val

		corr_vector.append(corr_sum/(seq_len-i))  #TODO: Never forget to bracket the subtraction term in denominator!

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

#Fixed parameters
n_val = 2
lambda_val = 8
w_val = 0.5

nfeatures = 10**n_val
struct_status = ["A", "G", "C", "U", "A-U", "U-A", "G-C", "C-G", "G-U", "U-G"]
free_energy = {"A": 0, "G": 0, "C": 0, "U": 0, "A-U": -2, "U-A": -2, "G-C": -3, "C-G": -3, "G-U": -1, "U-G": -1}  #Values based on https://doi.org/10.1371/journal.pone.0121501 ~ -1 kcal/mol per H-bond
feat_desc = [p for p in itertools.product(struct_status, repeat=n_val)]  #Enumerate all possible feature combinations

for k, feat in enumerate(feat_desc):
	feat_desc[k] = ','.join(feat)

with open(outpath+"Pseudo_structure_status_composition_v1.out", 'w') as f:
	header_str = "Target_RNA_ID\tTarget_RNA_name\t"+"\t".join(feat_desc)+"\t"
	for i in range(nfeatures+1, nfeatures+1+lambda_val):
		header_str = header_str+"PS_Feat_"+str(i)+"\t"

	print(header_str, file=f)

	for i, seq in enumerate(seqs_new):
		feature_vec = [0]*(nfeatures+lambda_val)
		(ss, mfe) = RNA.fold(seq)

		basepairs = []
		bracket_stack = []
		for j, char in enumerate(ss):
			if(char=="("):
				bracket_stack.append(j)  #Equivalent to stack push step
			if(char==")"):
				id1 = bracket_stack.pop()  #Returns the last element in the stack/list and removes it from the list
				id2 = j
				basepairs.append((id1, id2))

		#TODO: Refer Fig. 3 of https://doi.org/10.1371/journal.pone.0121501 article for information on comp_seq
		comp_seq = ""
		for j, base in enumerate(seq):
			index_found = 0
			for j1, j2 in basepairs:
				if(j1==j and index_found==0):
					comp_seq = comp_seq+seq[j2]
					index_found = 1
					break
				elif(j2==j and index_found==0):
					comp_seq = comp_seq+seq[j1]
					index_found = 1
					break

			if(index_found==0):
				comp_seq = comp_seq+"."

		#print(seq)
		#print(comp_seq)

		composition_vec = []
		for b1, b2 in zip(seq, comp_seq):
			if(b2=="."):
				composition_vec.append(b1)  #Unpaired base
			else:
				composition_vec.append(b1+"-"+b2)  #Base pair

		#Currently works only for n_val = 2
		#TODO: Generalize the seq_desc enumeration code for any n_val value
		seq_desc = []
		for j in range(len(composition_vec)):
			if(j+1<len(composition_vec)):
				seq_desc.append(composition_vec[j]+","+composition_vec[j+1])

		seq_desc_freq = Counter(seq_desc)

		feat_norm_freq_vec = []
		for s, feat in enumerate(feat_desc):
			try:
				feat_norm_freq_vec.append(seq_desc_freq[feat]/len(seq_desc))  #Normalize by number of description vectors possible
			except:
				feat_norm_freq_vec.append(0)
				continue

		corr_vector = calc_correlation_vector(seq, composition_vec, lambda_val, free_energy)
		for s, feat in enumerate(feat_desc):
			feature_vec[s] = np.round(feat_norm_freq_vec[s]/(sum(feat_norm_freq_vec)+(w_val*sum(corr_vector))), 5)
			
		for s in range(nfeatures, nfeatures+lambda_val):
			feature_vec[s] = np.round((w_val*corr_vector[s-nfeatures])/(sum(feat_norm_freq_vec)+(w_val*sum(corr_vector))), 5)

		feat_str = ""
		feat_str = feat_str+ids_new[i]+"\t"+names_new[i]+"\t"

		for feat in feature_vec:
			feat_str = feat_str+str(feat)+"\t"

		print(str(feat_str), file=f)

print("Pseudo structure status composition features calculated successfully")































