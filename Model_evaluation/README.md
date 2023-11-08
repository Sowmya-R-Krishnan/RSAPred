# Model evaluation with external test datasets - Regression and classification data

## Source of the datasets
* Regression dataset - HIV-1 TAR RNA QSAR study (Cai et al., 2022)
* Classification dataset - ROBIN repository (Yazdani et al., 2022)

## Description
* All the datasets provided have been pre-processed as per the pipeline explained in Dataset_preprocessing README file and provided here. They can be directly used for model evaluation.

## Prerequisites
* Anaconda or Miniconda with Python 3.8.
* A conda environment with the following libraries:
	* Python (>v3.8)
	* pandas
	* numpy
	* scipy
	* scikit-learn
	* pickle
	* matplotlib
	* seaborn
	* Generic libraries: sys, re, csv, os, collections...

## Sample commands (Follow the same order to reproduce the files provided in sample_output folder)
```
* python test_on_regression_dataset.py "./data/QSAR_dataset_final.csv" "../Feature_selection/data/Final_sample_dataset_v1.csv" "best_models.log" "test_results.csv"
* python test_on_classification_dataset.py <positive dataset> <negative dataset> "../Feature_selection/data/Final_sample_dataset_v1.csv" "best_models.log"
```

## Miscellaneous
* The complete dataset can be downloaded for each RNA subtype from the <a href="https://web.iitm.ac.in/bioinfo2/R_SIM/" target="_blank">R-SIM database</a>
* All programs take inputs as command-line arguments only
* The same pipeline is followed to prepare all custom test datasets used in the study, in the format necessary for model evaluation
* The header fields or columns provided in the sample dataset must be present (column order can be random) in your custom dataset file for the programs to work in the intended fashion
