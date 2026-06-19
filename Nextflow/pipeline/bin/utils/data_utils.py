import numpy as np
import pandas as pd
import torch

def vectorized_label(target, n_class):
    TARGET = np.array(target).reshape(-1)
    return np.eye(n_class)[TARGET]

def load_data(path, dtype, comparison):
    data = pd.read_excel(path)
    if 'Unnamed: 0' in data.columns:
        data = data.drop('Unnamed: 0', axis=1)
    x = data.drop([comparison], axis=1).values
    y = data.loc[:, [comparison]].values
    X = torch.from_numpy(x).type(dtype)
    Y = torch.from_numpy(vectorized_label(y, 2)).type(dtype)
    if torch.cuda.is_available():
        X = X.cuda()
        Y = Y.cuda()
    return X, Y

def load_pathway(path, dtype):
    pathway_mask = pd.read_excel(path, index_col=0)
    PATHWAY_MASK = torch.from_numpy(pathway_mask.values.astype(np.float32)).type(dtype)
    if torch.cuda.is_available():
        PATHWAY_MASK = PATHWAY_MASK.cuda()
    return PATHWAY_MASK
