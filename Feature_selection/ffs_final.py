#Progressively add non-redundant features to the regression model

import sys
import csv
import pickle
import pandas as pd
import numpy as np
from itertools import combinations
from sklearn import linear_model
from scipy import stats
from scipy.stats import pearsonr
from tqdm import tqdm
from multiprocessing import Pool
import itertools
from collections import Counter

data = sys.argv[1]
two_combos = sys.argv[2]
nomit = sys.argv[3]
n_feat = sys.argv[4]
pass_feat = sys.argv[5]
outpath = sys.argv[6]

n_feat = int(n_feat)
nomit = int(nomit)

df = pd.read_csv(data, sep='\t', header=0)

feats_final = []
if(pass_feat=="all"):
	feats_final = list(df.columns)[nomit:-1]
else:
	with open(pass_feat) as feats:
		for feat in feats.readlines():
			feat = feat.strip()
			feats_final.append(feat)

print("Total no. of features = "+str(len(feats_final)))

with open(two_combos, 'rb') as f:
	pair_corrs = pickle.load(f)

#Using set() instead of list() to speed up search process - Proven to be blazing fast (https://stackoverflow.com/questions/5993621/fastest-way-to-search-a-list-in-python)
pair_corrs = set(pair_corrs)

def corr_function(data, combination):
	corr = data[combination[0]].corr(data[combination[1]])
	return corr

def regression_model(data, combination, y):
	X = data[list(combination)]
	model = linear_model.LinearRegression()
	model.fit(X, y)
	predicted = model.predict(X)
	pearson_corr, p_value = pearsonr(y, predicted)
	return(combination, pearson_corr, model.coef_)

def find_pass_combos(combination, pair_corrs, data, y):
	feat_pairs = set(list(combinations(combination, 2)))
	feat_pass = [True for pair in feat_pairs if(pair in pair_corrs)]
	if(len(feat_pass)==len(feat_pairs)):
		return regression_model(data, combination, y)

y = df["pKd"]

#-------------------------------------------------------------------------------------------------------------------------------------
models = []
for n in range(n_feat, n_feat+1):
	models = []
	feat_combinations = set(list(combinations(feats_final, n)))
	print("No. of possible "+str(n)+" feature combinations: "+str(len(feat_combinations)))

	if(n>2):
		feat_combinations_pass = []
		for combination in feat_combinations:
			feat_pairs = set(list(combinations(combination, 2)))
			feat_pass = [True for pair in feat_pairs if(pair in pair_corrs)]
			if(len(feat_pass)==len(feat_pairs)):
				feat_combinations_pass.append(combination)
	else:
		feat_combinations_pass = [combination for combination in feat_combinations if(combination in pair_corrs)]

	print("Feasible "+str(n)+" feature combinations enumerated: "+str(len(feat_combinations_pass)))

	try:
		pool = Pool(4)
		results = pool.starmap(regression_model, zip(itertools.repeat(df), feat_combinations_pass, itertools.repeat(y)))
	finally:
		pool.close()
		pool.join()

	models = [(list(combo), results[i][0]) for i, combo in enumerate(feat_combinations_pass)]
	
	with open(outpath+"Aptamers_best_"+str(n)+"_feature_combos.log", 'w') as f:
		for feat_combo, corr in models:
			print('\t'.join(feat_combo)+"\t"+str(corr), file=f)

	print(str(n)+" feature combinations done.")

			
		


























