#Program to calculate model metrics for the best feature combination with 10-fold cross-validation

import sys
import csv
import pickle
import pandas as pd
import numpy as np
from itertools import combinations
from sklearn import linear_model, metrics
from scipy import stats
from scipy.stats import pearsonr
from tqdm import tqdm
from multiprocessing import Pool
import itertools
from sklearn.model_selection import cross_val_score
from sklearn.metrics import make_scorer
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score

data = sys.argv[1]
best_models = sys.argv[2]
nfeat = sys.argv[3]
nomit = sys.argv[4]
outpath = sys.argv[5]

nfeat = int(nfeat)
nomit = int(nomit)
df = pd.read_csv(data, sep='\t', header=0)

feats_final = []
feats_final = list(df.columns)[nomit:-1]
print("Total no. of features = "+str(len(feats_final)))

best_df = pd.read_csv(best_models, sep='\t', header=None)
print("No. of best models = "+str(len(best_df.index)))
print("No. of features in best models = "+str(len(best_df.columns)-1))

last_column = len(best_df.columns)-1
model_params = list(best_df[last_column])
top_model_r = max(model_params)

best_df.rename(columns={nfeat: "r_val"}, inplace=True)
top_models_df = best_df[best_df['r_val']==top_model_r]

for m_index, top_model in top_models_df.iterrows():
	top_model_specs = list(top_model)

	top_model_feat = top_model_specs[:-1]
	print("Top model features: "+','.join(top_model_feat))
	print("Correlation of top model: "+str(top_model_specs[-1]))

	#Creating a dataset out of best model features for cross-validation
	X_dataset = df[top_model_feat].to_numpy()
	Y_dataset = df['pKd'].to_numpy()

	mae_val = []
	mse_val = []

	Y_pred_loo = []
	Y_test_loo = []

	loo_cv = LeaveOneOut()
	cv_epoch = 1
	for train_index, test_index in loo_cv.split(X_dataset, Y_dataset):
		X_train, X_test = X_dataset[train_index], X_dataset[test_index]
		Y_train, Y_test = Y_dataset[train_index], Y_dataset[test_index]

		model = linear_model.LinearRegression()
		model.fit(X_train, Y_train)
		Y_pred = model.predict(X_test)

		Y_pred_loo.append(Y_pred[0])
		Y_test_loo.append(Y_test[0])
		
		mae = mean_absolute_error(Y_test, Y_pred)
		mse = mean_squared_error(Y_test, Y_pred)
		mae_val.append(mae)
		mse_val.append(mse)
		cv_epoch = cv_epoch + 1

	with open(outpath+"Results_JK_"+str(m_index)+".csv", "w") as f:
		for val1, val2 in zip(Y_test_loo, Y_pred_loo):
			print(str(val1)+"\t"+str(val2), file=f)

	pearson_corr, p_value = pearsonr(Y_test_loo, Y_pred_loo)
	r2 = r2_score(Y_test_loo, Y_pred_loo)

	print(" ")
	print("Model",m_index,"LOO CV results")
	print("-------------------")
	print("Average MAE: "+str(np.mean(mae_val)))
	print("Average MSE: "+str(np.mean(mse_val)))
	print("R2 score: "+str(r2))
	print("Pearson correlation: "+str(pearson_corr))
	print("-------------------")
	print(" ")



















