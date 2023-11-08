#Calculate nucleotide composition features of input RNA sequence upto length (K) = 6

import sys
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

#Calculate the composition vector of a list of strings/substrings in a given sequence/string
def calc_composition_vector(sequence, string_list):
	comp_vec = [0]*len(string_list)
	for i, string in enumerate(string_list):
		string = ''.join(list(string))
		if(sequence.count(string)>0):
			comp_vec[i] = comp_vec[i] + sequence.count(string)

	for i, count in enumerate(comp_vec):
		comp_vec[i] = np.round(count/len(sequence), 5)
	
	return comp_vec		

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
		print(names[i], ids[i])


#Based on Cartesian Product - https://stackoverflow.com/questions/3099987/generating-permutations-with-repetitions
bases_prime = ["A", "G", "C", "U"]
di_strings = [p for p in itertools.product(bases_prime, repeat=2)]
tri_strings = [p for p in itertools.product(bases_prime, repeat=3)]
tetra_strings = [p for p in itertools.product(bases_prime, repeat=4)]
penta_strings = [p for p in itertools.product(bases_prime, repeat=5)]
hexa_strings = [p for p in itertools.product(bases_prime, repeat=6)]

#Dimensionality of the corresponding composition vectors
print("No. of dinucleotides possible: "+str(len(di_strings)))
print("No. of trinucleotides possible: "+str(len(tri_strings)))
print("No. of tetranucleotides possible: "+str(len(tetra_strings)))
print("No. of pentanucleotides possible: "+str(len(penta_strings)))
print("No. of hexanucleotides possible: "+str(len(hexa_strings)))


with open(outpath+"Mononucleotide_v1.out", "w") as f:
	comb_strs = []
	for string in bases_prime:
		string = ''.join(list(string))
		comb_strs.append(string)

	print("Target_RNA_ID\tTarget_RNA_name\t"+"\t".join(comb_strs), file=f)
	for i, seq in enumerate(seqs_new):
		mono_vec = calc_composition_vector(seq, bases_prime)
		print(ids_new[i]+"\t"+names_new[i]+"\t"+"\t".join([str(i) for i in mono_vec]), file=f)
	f.close()
		
with open(outpath+"Dinucleotide_v1.out", "w") as f:
	comb_strs = []
	for string in di_strings:
		string = ''.join(list(string))
		comb_strs.append(string)

	print("Target_RNA_ID\tTarget_RNA_name\t"+"\t".join(comb_strs), file=f)
	for i, seq in enumerate(seqs_new):
		di_vec = calc_composition_vector(seq, di_strings)
		print(ids_new[i]+"\t"+names_new[i]+"\t"+"\t".join([str(i) for i in di_vec]), file=f)
	f.close()

with open(outpath+"Trinucleotide_v1.out", "w") as f:
	comb_strs = []
	for string in tri_strings:
		string = ''.join(list(string))
		comb_strs.append(string)

	print("Target_RNA_ID\tTarget_RNA_name\t"+"\t".join(comb_strs), file=f)
	for i, seq in enumerate(seqs_new):
		tri_vec = calc_composition_vector(seq, tri_strings)
		print(ids_new[i]+"\t"+names_new[i]+"\t"+"\t".join([str(i) for i in tri_vec]), file=f)
	f.close()

with open(outpath+"Tetranucleotide_v1.out", "w") as f:
	comb_strs = []
	for string in tetra_strings:
		string = ''.join(list(string))
		comb_strs.append(string)

	print("Target_RNA_ID\tTarget_RNA_name\t"+"\t".join(comb_strs), file=f)
	for i, seq in enumerate(seqs_new):
		tetra_vec = calc_composition_vector(seq, tetra_strings)
		print(ids_new[i]+"\t"+names_new[i]+"\t"+"\t".join([str(i) for i in tetra_vec]), file=f)
	f.close()

with open(outpath+"Pentanucleotide_v1.out", "w") as f:
	comb_strs = []
	for string in penta_strings:
		string = ''.join(list(string))
		comb_strs.append(string)

	print("Target_RNA_ID\tTarget_RNA_name\t"+"\t".join(comb_strs), file=f)
	for i, seq in enumerate(seqs_new):
		penta_vec = calc_composition_vector(seq, penta_strings)
		print(ids_new[i]+"\t"+names_new[i]+"\t"+"\t".join([str(i) for i in penta_vec]), file=f)
	f.close()

with open(outpath+"Hexanucleotide_v1.out", "w") as f:
	comb_strs = []
	for string in hexa_strings:
		string = ''.join(list(string))
		comb_strs.append(string)

	print("Target_RNA_ID\tTarget_RNA_name\t"+"\t".join(comb_strs), file=f)
	for i, seq in enumerate(seqs_new):
		hexa_vec = calc_composition_vector(seq, hexa_strings)
		print(ids_new[i]+"\t"+names_new[i]+"\t"+"\t".join([str(i) for i in hexa_vec]), file=f)
	f.close()

print("Nucleotide composition features calculated successfully")

















