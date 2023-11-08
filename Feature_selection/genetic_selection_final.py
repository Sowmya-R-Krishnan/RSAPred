#Perform backward feature elimination using a MLR model to compare with forward feature selection results

import sys
import csv
import pandas as pd
import numpy as np
import pickle
from scipy import stats
from scipy.stats import pearsonr
from sklearn import linear_model, metrics
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from genetic_selection import GeneticSelectionCV
from sklearn.metrics import make_scorer
from matplotlib import pyplot as plt
import seaborn as sns

data = sys.argv[1]
nomit = sys.argv[2]
outfile = sys.argv[3]

nomit = int(nomit)
df = pd.read_csv(data, sep='\t', header=0)

feats_final = []
feats_final = list(df.columns)[nomit:-1]
print("Total no. of features = "+str(len(feats_final)))

X_dataset = df[feats_final].to_numpy()
Y_dataset = df['pKd'].to_numpy()

#Pearson correlation scoring metric for CV
def pearson_func(data1, data2):
	pearson_corr, p_value = pearsonr(data1, data2)
	return pearson_corr

pearson_r_scorer = make_scorer(pearson_func)

mlr = linear_model.LinearRegression()
selector = GeneticSelectionCV(mlr, cv=10, verbose=1, scoring=pearson_r_scorer, max_features=55, n_jobs=-1)
selector.fit(X_dataset, Y_dataset)

print("Optimal number of features: "+str(selector.n_features_))

with open(outfile, "wb") as f:
	pickle.dump(selector, file=f)
print("Genetic selection results pickled successfully.")

final_feat_set = []
for i, val in enumerate(selector.support_):
	if(val):
		final_feat_set.append(feats_final[i])

print(final_feat_set)

mlr = linear_model.LinearRegression()
X_new = df[final_feat_set].to_numpy()
Y_new = df['pKd'].to_numpy()
mlr.fit(X_new, Y_new)
Y_pred = mlr.predict(X_new)
pearson_corr, p_value = pearsonr(Y_new, Y_pred)
mae = mean_absolute_error(Y_new, Y_pred)
mse = mean_squared_error(Y_new, Y_pred)
r2 = r2_score(Y_new, Y_pred)

print("MAE: "+str(mae))
print("MSE: "+str(mse))
print("R2 score: "+str(r2))
print("Pearson correlation: "+str(pearson_corr))


