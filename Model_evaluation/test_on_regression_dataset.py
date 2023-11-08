#Program to test the performance of a model on an external test dataset - regression dataset

import sys
import csv
import pickle
import pandas as pd
import numpy as np
from itertools import combinations
from sklearn import linear_model, metrics, ensemble
from scipy import stats
from scipy.stats import pearsonr
import itertools
from sklearn.model_selection import cross_val_score
from sklearn.metrics import make_scorer
from sklearn.metrics import mean_absolute_error

data = sys.argv[1]
dataset = sys.argv[2]
best_models = sys.argv[3]
outfile = sys.argv[4]

best_df = pd.read_csv(best_models, sep='\t', header=None)
print("No. of best models = "+str(len(best_df.index)))
print("No. of features in best models = "+str(len(best_df.columns)-1))

last_column = len(best_df.columns)-1
model_params = np.asarray(list(best_df[last_column]))
top_model = np.argmax(model_params)
top_model_specs = list(best_df.iloc[top_model])

top_model_feat = top_model_specs[:-1]
print("Top model features: "+','.join(top_model_feat))
print("Correlation of top model: "+str(top_model_specs[-1]))

final_df = pd.read_csv(data, sep="\t", header=0)

model_dataset = pd.read_csv(dataset, sep='\t', header=0)
X_dataset = model_dataset[top_model_feat].to_numpy()
Y_dataset = model_dataset['pKd'].to_numpy()
model = linear_model.LinearRegression()
model.fit(X_dataset, Y_dataset)

Y_test = final_df['pKd'].to_numpy()
test_X = final_df[top_model_feat].to_numpy()
Y_pred = model.predict(test_X)
print(pearsonr(Y_test, Y_pred), mean_absolute_error(Y_test, Y_pred))

out = open(outfile, "w")
for val1, val2 in zip(Y_test, Y_pred):
	print(val1, val2, file=out)






















