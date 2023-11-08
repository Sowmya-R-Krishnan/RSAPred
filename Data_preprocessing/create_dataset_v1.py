#Mapping the RNA and small molecule features to create the final dataset for model training

import sys
import csv
import pandas as pd

rna_feat = sys.argv[1]
mol_feat = sys.argv[2]
dataset_raw = sys.argv[3]
outfile = sys.argv[4]

rna_df = pd.read_csv(rna_feat, sep="\t", header=0)
rna_df.dropna(axis=1, how='any', thresh=None, subset=None, inplace=True)

mol_df = pd.read_csv(mol_feat, sep=",", header=0)
mol_df.dropna(axis=1, how='any', thresh=None, subset=None, inplace=True)

data_df = pd.read_csv(dataset_raw, sep="\t", header=0)

entries = list(data_df['Entry_ID'])
rna_ids = list(data_df['Target_RNA_ID'])
mol_ids = list(data_df['Molecule_ID'])

colnames = []
colnames.extend(list(rna_df.columns))
colnames.extend(list(mol_df.columns))
colnames.extend(["pKd"])
colnames.remove("name")
colnames.insert(2, "name")
final_df = pd.DataFrame(columns = colnames)
print(len(rna_df.index), len(mol_df.index), len(data_df.index))

for entry, rid, mid in zip(entries, rna_ids, mol_ids):
	try:
		data_point = pd.DataFrame(columns = colnames)
		rna_feat = rna_df.loc[rna_df['Target_RNA_ID'] == rid].to_dict('list')
		mol_feat = mol_df.loc[mol_df['name'] == mid].to_dict('list')
		kd_val = {'pKd': [float(data_df.loc[data_df['Entry_ID'] == entry]['pKd'])]}

		data_point = rna_feat
		data_point.update(mol_feat)
		data_point.update(kd_val)
		row = pd.DataFrame.from_dict(data_point)
		#print(data_point)
		final_df = pd.concat([final_df, row], ignore_index=True)  #Append the row to the final dataframe to create the dataset
		#print(final_df)
	except:
		print("Failed: ", entry, rid, mid)
		continue

final_df.to_csv(outfile, sep="\t", index=False, header=True)

	
	
































