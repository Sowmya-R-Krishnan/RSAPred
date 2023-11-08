#Calculate triplet structure composition features for RNA sequence - based on repRNA method

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
tri_strings = [p for p in itertools.product(bases_prime, repeat=3)]

for i, tri_str in enumerate(tri_strings):
	tri_strings[i] = ''.join(tri_str)

ss_combos = ["(((", "((.", "(..", "(.(", ".((", ".(.", "..(", "..."]


with open(outpath+"Triplet_features_v1.out", 'w') as f:
	header_str = "Target_RNA_ID\tTarget_RNA_name\t"
	for base in bases_prime:
		for combo in ss_combos:
			header_str = header_str+base+combo+"\t"

	print(header_str, file=f)

	for i, seq in enumerate(seqs_new):
		feature_vec = [0]*(4*len(ss_combos))  #8 per nucleotide

		(ss, mfe) = RNA.fold(seq)  #Predicting the secondary structure using ViennaRNA package
		ss = ss.replace(")", "(")  #Replacing closing brackets with opening brackets to match possible ss_combos

		trinucleotides = []
		tri_ss = []
		tri_ss_only = []
		for j in range(len(seq)):
			if(j+2<len(seq)):
				trinucleotides.append(seq[j]+seq[j+1]+seq[j+2])
				tri_ss.append((seq[j+1], ss[j]+ss[j+1]+ss[j+2]))
				tri_ss_only.append(ss[j]+ss[j+1]+ss[j+2])

		trinucl_freq = Counter(tri_ss_only)
		tri_ss_freq = Counter(tri_ss)

		index = 0
		for base in bases_prime:
			for combo in ss_combos:
				try:
					feature_vec[index] = np.round(tri_ss_freq[(base, combo)]/sum(tri_ss_freq.values()), 5)  #Normalize by the total number of possible triplets in the sequence
					index = index + 1
				except:
					index = index + 1
					continue

		feat_str = ""
		feat_str = feat_str+ids_new[i]+"\t"+names_new[i]+"\t"

		for feat in feature_vec:
			feat_str = feat_str+str(feat)+"\t"

		print(feat_str, file=f)

print("Triplet sequence-structure composition features calculated successfully")
































