#!/usr/bin/env python
"""
PASNet: Pathway-Associated Sparse Deep Neural Network for Genomic Prediction

This script trains and evaluates a PASNet model for binary classification tasks
in bioinformatics, particularly for genomic prediction with pathway information.

FIXED VERSION: Resolves numpy.object_ dtype conversion error
"""

import numpy as np
import pandas as pd
import math
import copy
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

import torch
import torch.nn.functional as F
import torch.nn as nn
import torch.optim as optim

from sklearn.metrics import precision_score, recall_score, confusion_matrix
from sklearn.metrics import roc_curve, auc as sklearn_auc
from sklearn.metrics import roc_auc_score, f1_score
from sklearn.metrics import roc_curve, auc
import shap
import pickle

# Set random seed for reproducibility
torch.manual_seed(0)

def vectorized_label(target, n_class):
	'''convert target(y) to be one-hot encoding format(dummy variable)
	'''
	TARGET = np.array(target).reshape(-1)

	return np.eye(n_class)[TARGET]


def load_data(path, dtype):
	'''Load data, and then covert it to a Pytorch tensor.
	Input:
		path: path to input dataset (which is expected to be a csv file).
		dtype: define the data type of tensor (i.e. dtype=torch.FloatTensor)
	Output:
		X: a Pytorch tensor of 'x'.
		Y: a Pytorch tensor of 'y'(one-hot encoding).
	'''
	data = pd.read_excel(path)
	data = data.drop('Unnamed: 0', axis=1)

	x = data.drop(["Status"], axis = 1).values
	y = data.loc[:, ["Status"]].values

	X = torch.from_numpy(x).type(dtype)
	Y = torch.from_numpy(vectorized_label(y, 2)).type(dtype)
	###if gpu is being used
	if torch.cuda.is_available():
		X = X.cuda()
		Y = Y.cuda()
	###
	return(X, Y)


def load_pathway(path, dtype):
	'''Load a bi-adjacency matrix of pathways, and then covert it to a Pytorch tensor.
	Input:
		path: path to input dataset (which is expected to be a csv file).
		dtype: define the data type of tensor (i.e. dtype=torch.FloatTensor)
	Output:
		PATHWAY_MASK: a Pytorch tensor of the bi-adjacency matrix of pathways.
	'''
	# pathway_mask = pd.read_csv(path, index_col = 0).as_matrix()
	pathway_mask = pd.read_excel(path, index_col = 0)

	PATHWAY_MASK = torch.from_numpy(pathway_mask.values.astype(np.float32)).type(dtype)
	###if gpu is being used
	if torch.cuda.is_available():
		PATHWAY_MASK = PATHWAY_MASK.cuda()
	###
	return(PATHWAY_MASK)


