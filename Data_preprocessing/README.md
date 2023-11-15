# Dataset preprocessing stages

## Source of the features
* RNA sequence-based features from repRNA package custom implementation
* Small molecule structure-based features from Mordred package

## Description
* RNA sequence-based features are obtained using 4 custom python scripts (/data/rna_data.csv) and combined to get the final set of RNA features
* OpenBabel is used to obtain the 3D structures with hydrogens in SDF format for the unique small molecules in the dataset (/data/mol_data.csv)
* Mordred CLI is used to obtain the small molecule features using the sample command provided below
* Finally, RNA sequence-based features and small molecule features are combined together to obtain the training dataset

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
	* mordred
	* openbabel
        * viennarna
	* Generic libraries: sys, re, csv, os, collections...

## Sample commands (Follow the same order to reproduce the files provided in sample_output folder)
```
* python calc_kmer_composition_features_v1.py "./data/rna_data.csv" "./sample_output/"
* python calc_ppseudoDNC_features_v1.py "./data/rna_data.csv" "./data/Physicochemical_indices_RNA.csv" "./sample_output/"
* python calc_triplet_str_composition_features_v1.py "./data/rna_data.csv" "./sample_output/"
* python calc_pseudo_structure_composition_v1.py "./data/rna_data.csv" "./sample_output/"
* obabel -ismi ./data/mol_data.smi -osdf -O ./sample_output/mol_data.sdf --gen3d -h
* python -m mordred ./sample_output/mol_data.sdf -t sdf -o ./sample_output/Mol_features_v1.csv -3
* python combine_RNA_features_v1.py "./sample_output/Mononucleotide_v1.out" "./sample_output/Dinucleotide_v1.out" "./sample_output/Trinucleotide_v1.out" "./sample_output/Tetranucleotide_v1.out" "./sample_output/pPseudoDNC_features_v1.out" "./sample_output/Triplet_features_v1.out" "./sample_output/Pseudo_structure_status_composition_v1.out" "./sample_output/RNA_features_v1.csv"
* python create_dataset_v1.py "./sample_output/RNA_features_v1.csv" "./sample_output/Mol_features_v1.csv" "./data/sample_data.csv" "./sample_output/Final_sample_dataset_v1.csv"
```

## Miscellaneous
* The complete dataset can be downloaded for each RNA subtype from the <a href="https://web.iitm.ac.in/bioinfo2/R_SIM/" target="_blank">R-SIM database</a>
* All programs take inputs as command-line arguments only
* The same pipeline is followed to prepare all custom test datasets used in the study, in the format necessary for model evaluation
* The header fields or columns provided in the sample dataset must be present (column order can be random) in your custom dataset file for the programs to work in the intended fashion
