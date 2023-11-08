#Program to compare every successive set of feature combination results and extract the best models

import sys
import csv
import pandas as pd

data1 = sys.argv[1]
data2 = sys.argv[2]
least_featnum = sys.argv[3]
outfile = sys.argv[4]

least_featnum = int(least_featnum)

df1 = pd.read_csv(data1, sep='\t', header=None, on_bad_lines='skip')
df2 = pd.read_csv(data2, sep='\t', header=None, on_bad_lines='skip')

min_model_metric = max(list(df1[least_featnum]))  #TODO: Change column header according to number of columns - Last column has model metric
best_combo = df1.loc[df1[least_featnum]==min_model_metric]

print("Intial best model:", len(best_combo.index))
print(best_combo)
print("")

best_models = df2.loc[df2[least_featnum+1]>min_model_metric]
print("Final best models:", len(best_models.index))
print(best_models)

best_models.to_csv(outfile, sep='\t', index=False, header=False)





























