#Program to calculate model metrics for the best feature combination with 10-fold cross-validation

import sys
import csv
import pickle
import pandas as pd
import numpy as np
from itertools import combinations
from sklearn import linear_model, metrics
from sklearn.model_selection import StratifiedKFold
from sklearn.utils import shuffle
from scipy import stats
from scipy.stats import pearsonr
from tqdm import tqdm
from multiprocessing import Pool
import itertools
from sklearn.model_selection import cross_val_score
from sklearn.metrics import make_scorer
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.preprocessing import KBinsDiscretizer
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score

data = sys.argv[1]
best_models = sys.argv[2]
nomit = sys.argv[3]
outpath = sys.argv[4]

#------------------------------------------------------------------------------------------------------------------------------------
#Adapted from https://github.com/scikit-learn/scikit-learn/issues/4757#issuecomment-694924478
class regressor_stratified_cv:
    def __init__(self,n_splits=10,n_repeats=2,group_count=10,random_state=0,strategy='quantile'):
        self.group_count=group_count
        self.strategy=strategy
        self.cvkwargs=dict(n_splits=n_splits,n_repeats=n_repeats,random_state=random_state)  #Added shuffle here
        self.cv=RepeatedStratifiedKFold(**self.cvkwargs)
        self.discretizer=KBinsDiscretizer(n_bins=self.group_count,encode='ordinal',strategy=self.strategy)  
            
    def split(self,X,y,groups=None):
        kgroups=self.discretizer.fit_transform(y[:,None])[:,0]
        return self.cv.split(X,kgroups,groups)
    
    def get_n_splits(self,X,y,groups=None):
        return self.cv.get_n_splits(X,y,groups)
#------------------------------------------------------------------------------------------------------------------------------------

nomit = int(nomit)
df = pd.read_csv(data, sep='\t', header=0)
#print("Total no. of features = "+str(len(df.columns)-nomit-1))

feats_final = []
feats_final = list(df.columns)[nomit:-1]
print("Total no. of features = "+str(len(feats_final)))

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

#Creating a dataset out of best model features for cross-validation
df = shuffle(df)
X_dataset = df[top_model_feat].to_numpy()
Y_dataset = df['pKd'].to_numpy()

skf = regressor_stratified_cv(n_splits=10, n_repeats=1, random_state=42, group_count=5, strategy='uniform')

mae_val = []
mse_val = []
r2_val = []
pearson_rval = []

for train_index, test_index in skf.split(X_dataset, Y_dataset):
	X_train, X_test = X_dataset[train_index], X_dataset[test_index]
	Y_train, Y_test = Y_dataset[train_index], Y_dataset[test_index]

	model = linear_model.LinearRegression()
	model.fit(X_train, Y_train)
	Y_pred = model.predict(X_test)
	
	mae = mean_absolute_error(Y_test, Y_pred)
	mse = mean_squared_error(Y_test, Y_pred)
	r2 = r2_score(Y_test, Y_pred)
	pearson_corr, p_value = pearsonr(Y_test, Y_pred)
	mae_val.append(mae)
	mse_val.append(mse)
	r2_val.append(r2)
	pearson_rval.append(pearson_corr)

print(" ")
print("10-fold CV results")
print("-------------------")
print("Average MAE: "+str(np.mean(mae_val)))
print("Average MSE: "+str(np.mean(mse_val)))
print("Average R2 score: "+str(np.mean(r2_val)))
print("Average Pearson correlation: "+str(np.mean(pearson_rval)))



















