#Program to extract feature statistics from successive regression results and fix a subset of features with good results in at least 5 models
#TODO: Probe different cut-offs to shortlist features

import sys
import csv
import pandas as pd
from collections import Counter

data = sys.argv[1]
nfeat = sys.argv[2]
outfile1 = sys.argv[3]  #Feature stats file
outfile2 = sys.argv[4]  #Pass features file for next combination of regression

df = pd.read_csv(data, sep='\t', header=None, on_bad_lines='skip')

feats_all = []
for i in range(int(nfeat)):
	feat_list = list(df[i])
	feats_all.extend(feat_list)

feat_freq = Counter(feats_all)
feat_freq_df = pd.DataFrame(feat_freq.items(), columns=['Feature', 'Frequency'])
feat_freq_df.sort_values(by=['Frequency'], ascending=False, inplace=True)
feat_freq_df.to_csv(outfile1, sep='\t', header=True, index=False)

feat_pass_list = list(feat_freq_df['Feature'][0:200])  #Top 200 features for model
print("No. of features which appear in at least 5 models: "+str(len(feat_pass_list)))

with open(outfile2, 'w') as f:
	for feat in feat_pass_list:
		print(feat, file=f)































