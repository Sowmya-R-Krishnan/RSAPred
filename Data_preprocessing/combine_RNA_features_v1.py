#Python program to combine multiple RNA feature calculation results into a single feature file for each RNA sequence

import sys
import csv
import pandas as pd

f1_data = sys.argv[1]  #Mononucleotide composition
f2_data = sys.argv[2]  #Dinucleotide composition
f3_data = sys.argv[3]  #Trinucleotide composition
f4_data = sys.argv[4]  #Tetranucleotide composition
f5_data = sys.argv[5]  #pPseudoDNC composition
f6_data = sys.argv[6]  #Triplet structure composition
f7_data = sys.argv[7]  #Pseudo structure composition
outfile = sys.argv[8]  #Output file with all features

mono_df = pd.read_csv(f1_data, sep="\t", header=0)
di_df = pd.read_csv(f2_data, sep="\t", header=0)
tri_df = pd.read_csv(f3_data, sep="\t", header=0)
tetra_df = pd.read_csv(f4_data, sep="\t", header=0)
ppseudoDNC_df = pd.read_csv(f5_data, sep="\t", header=0)
triplet_df = pd.read_csv(f6_data, sep="\t", header=0)
pseudo_struct_df = pd.read_csv(f7_data, sep="\t", header=0)

#Drop the entire column in the dataframe if at least one of the values is NaN
mono_df.dropna(axis=1, how='any', thresh=None, subset=None, inplace=True)
di_df.dropna(axis=1, how='any', thresh=None, subset=None, inplace=True)
tri_df.dropna(axis=1, how='any', thresh=None, subset=None, inplace=True)
tetra_df.dropna(axis=1, how='any', thresh=None, subset=None, inplace=True)
ppseudoDNC_df.dropna(axis=1, how='any', thresh=None, subset=None, inplace=True)
triplet_df.dropna(axis=1, how='any', thresh=None, subset=None, inplace=True)
pseudo_struct_df.dropna(axis=1, how='any', thresh=None, subset=None, inplace=True)

#print(mono_df)

merge_df = pd.merge(mono_df, di_df,  how='left', on=['Target_RNA_ID','Target_RNA_name'])
merge_df = pd.merge(merge_df, tri_df,  how='left', on=['Target_RNA_ID','Target_RNA_name'])
merge_df = pd.merge(merge_df, tetra_df,  how='left', on=['Target_RNA_ID','Target_RNA_name'])
merge_df = pd.merge(merge_df, ppseudoDNC_df,  how='left', on=['Target_RNA_ID','Target_RNA_name'])
merge_df = pd.merge(merge_df, triplet_df,  how='left', on=['Target_RNA_ID','Target_RNA_name'])
merge_df = pd.merge(merge_df, pseudo_struct_df,  how='left', on=['Target_RNA_ID','Target_RNA_name'])
#print(merge_df)

expected_columns = len(mono_df.columns)+len(di_df.columns)-2+len(tri_df.columns)-2+len(tetra_df.columns)-2+len(ppseudoDNC_df.columns)-2+len(triplet_df.columns)-2+len(pseudo_struct_df.columns)-2
observed_columns = len(merge_df.columns)

#Check if number of columns of final df exactly as expected - TODO: Check after each merge to pinpoint the df with problem??
if(expected_columns==observed_columns):
	print("Feature merge successful!")
	merge_df.to_csv(outfile, sep='\t', index=False, header=True)
else:
	print("Feature merge failed. Kindly check dataframe headers for any potential overlaps and re-try.")

































