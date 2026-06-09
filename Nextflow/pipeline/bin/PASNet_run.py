#!/usr/bin/env python
"""
PASNet: Pathway-Associated Sparse Deep Neural Network for Genomic Prediction

This script trains and evaluates a PASNet model for binary classification tasks
in bioinformatics, particularly for genomic prediction with pathway information.
"""

import argparse
import numpy as np
import pandas as pd
import math
import copy
import os
from scipy.interpolate import interp1d
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend to avoid display/permission errors
import matplotlib.pyplot as plt

import torch
import torch.nn.functional as F
import torch.nn as nn
import torch.optim as optim

from sklearn.metrics import precision_score, recall_score, confusion_matrix
from sklearn.metrics import roc_curve, auc as sklearn_auc
from sklearn.metrics import roc_auc_score, f1_score
import shap
import pickle


def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='PASNet Training')
    parser.add_argument('--train-data',    required=True,  help='Train_Data Input')
    parser.add_argument('--val-data',      required=True,  help='Evaluation_Data Input')
    parser.add_argument('--pathway-mask',  required=True,  help='Pathway mask matrix file')
    parser.add_argument('--outdir',        required=True,  help='Output directory')
    parser.add_argument('--comparison',    required=True,  help='Comparison column name')
    parser.add_argument('--val-data-grid', required=False, help='Evaluation_Griding_Data Input')

    return parser.parse_args()


# Set random seed for reproducibility
torch.manual_seed(0)

# Global dtype
dtype = torch.FloatTensor


def vectorized_label(target, n_class):
	'''convert target(y) to be one-hot encoding format(dummy variable)
	'''
	TARGET = np.array(target).reshape(-1)

	return np.eye(n_class)[TARGET]


def load_data(path, dtype, comparison):
    """Load data and convert it to a Pytorch tensor.

    Input:
        path:       path to input dataset (expected to be an xlsx file).
        dtype:      data type of tensor (e.g. torch.FloatTensor)
        comparison: column name used as the label / target variable
    Output:
        X: a Pytorch tensor of 'x'.
        Y: a Pytorch tensor of 'y' (one-hot encoding).
    """
    data = pd.read_excel(path)
    # Drop the auto-index column if present
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
    """Load a bi-adjacency matrix of pathways and convert it to a Pytorch tensor.

    Input:
        path:  path to input dataset (expected to be an xlsx file).
        dtype: data type of tensor (e.g. torch.FloatTensor)
    Output:
        PATHWAY_MASK: a Pytorch tensor of the bi-adjacency matrix of pathways.
    """
    pathway_mask = pd.read_excel(path, index_col=0)
    PATHWAY_MASK = torch.from_numpy(pathway_mask.values.astype(np.float32)).type(dtype)

    if torch.cuda.is_available():
        PATHWAY_MASK = PATHWAY_MASK.cuda()

    return PATHWAY_MASK


def auc(y_true, y_pred):
    """Compute ROC-AUC score."""
    if torch.cuda.is_available():
        y_true = y_true.cpu().detach()
        y_pred = y_pred.cpu().detach()

    score = roc_auc_score(y_true.detach().numpy(), y_pred.detach().numpy())
    return score


def calc_auc(y_true, y_pred):
    """Compute ROC-AUC score (alias used during evaluation loop)."""
    return auc(y_true, y_pred)


def f1(y_true, y_pred):
    """Compute F1 score."""
    y    = torch.argmax(y_true, dim=1)
    pred = torch.argmax(y_pred, dim=1)

    if torch.cuda.is_available():
        y    = y.cpu().detach()
        pred = pred.cpu().detach()

    score = f1_score(y.detach().numpy(), pred.detach().numpy())
    return score


def bce_for_one_class(predict, target, lts=False):
    """Calculate cross entropy in average for samples belonging to the same class.

    lts = False: non-LTS samples are obtained.
    """
    lts_idx = torch.argmax(target, dim=1)
    idx = 1 if lts else 0
    y    = target[lts_idx == idx]
    pred = predict[lts_idx == idx]
    cost = F.binary_cross_entropy(pred, y)
    return cost


def binary_cross_entropy_for_imbalance(predict, target):
    """Calculate cross entropy for imbalanced data in binary classification."""
    total_cost = (bce_for_one_class(predict, target, lts=True)
                  + bce_for_one_class(predict, target, lts=False))
    return total_cost


class PASNet(nn.Module):
    def __init__(self, In_Nodes, Pathway_Nodes, Hidden_Nodes, Out_Nodes, Pathway_Mask):
        super(PASNet, self).__init__()
        self.sigmoid      = nn.Sigmoid()
        self.softmax      = nn.Softmax(dim=1)
        self.pathway_mask = Pathway_Mask

        # gene layer --> pathway layer
        self.sc1 = nn.Linear(In_Nodes,      Pathway_Nodes)
        # pathway layer --> hidden layer
        self.sc2 = nn.Linear(Pathway_Nodes, Hidden_Nodes)
        # hidden layer --> output layer
        self.sc3 = nn.Linear(Hidden_Nodes,  Out_Nodes)

        # randomly select a small sub-network
        self.do_m1 = torch.ones(Pathway_Nodes)
        self.do_m2 = torch.ones(Hidden_Nodes)

        if torch.cuda.is_available():
            self.do_m1 = self.do_m1.cuda()
            self.do_m2 = self.do_m2.cuda()

    def forward(self, x):
        # Force the connections between gene layer and pathway layer
        self.sc1.weight.data = self.sc1.weight.data.mul(self.pathway_mask)
        x = self.sigmoid(self.sc1(x))
        if self.training:
            x = x.mul(self.do_m1)
        x = self.sigmoid(self.sc2(x))
        if self.training:
            x = x.mul(self.do_m2)
        x = self.softmax(self.sc3(x))
        return x


def dropout_mask(n_node, drop_p):
    """Construct a binary matrix to randomly drop nodes in a layer."""
    keep_p = 1.0 - drop_p
    mask   = torch.Tensor(np.random.binomial(1, keep_p, size=n_node))
    if torch.cuda.is_available():
        mask = mask.cuda()
    return mask


def s_mask(sparse_level, param_matrix, nonzero_param_1D, dtype):
    """Construct a binary sparsity mask for weights between two consecutive layers."""
    non_neg_param_1D = torch.abs(nonzero_param_1D)
    num_param        = nonzero_param_1D.size(0)
    top_k            = math.ceil(num_param * (100 - sparse_level) * 0.01)

    sorted_non_neg_param_1D, _ = torch.topk(non_neg_param_1D, top_k)
    param_mask = torch.abs(param_matrix) > sorted_non_neg_param_1D.min()
    param_mask = param_mask.type(dtype)

    if torch.cuda.is_available():
        param_mask = param_mask.cuda()

    return param_mask


def trainPASNet(train_x, train_y, eval_x, eval_y, pathway_mask,
                In_Nodes, Pathway_Nodes, Hidden_Nodes, Out_Nodes,
                Learning_Rate, L2_Lambda, nEpochs, Dropout_Rates,
                optimizer="Adam"):
    """Train a PASNet model and return predictions, losses, and the trained model."""

    net = PASNet(In_Nodes, Pathway_Nodes, Hidden_Nodes, Out_Nodes, pathway_mask)
    if torch.cuda.is_available():
        net.cuda()

    if optimizer == "SGD":
        opt = optim.SGD(net.parameters(),  lr=Learning_Rate, weight_decay=L2_Lambda)
    else:
        opt = optim.Adam(net.parameters(), lr=Learning_Rate, weight_decay=L2_Lambda)

    for epoch in range(nEpochs):
        net.train()
        opt.zero_grad()

        # Randomize dropout masks
        net.do_m1 = dropout_mask(Pathway_Nodes, Dropout_Rates[0])
        net.do_m2 = dropout_mask(Hidden_Nodes,  Dropout_Rates[1])

        pred = net(train_x)
        loss = binary_cross_entropy_for_imbalance(pred, train_y)
        loss.backward()
        opt.step()

        # Force pathway mask
        net.sc1.weight.data = net.sc1.weight.data.mul(net.pathway_mask)

        # Obtain sub-network connections
        do_m1_grad      = copy.deepcopy(net.sc2.weight._grad.data)
        do_m2_grad      = copy.deepcopy(net.sc3.weight._grad.data)
        do_m1_grad_mask = torch.where(do_m1_grad == 0, do_m1_grad, torch.ones_like(do_m1_grad))
        do_m2_grad_mask = torch.where(do_m2_grad == 0, do_m2_grad, torch.ones_like(do_m2_grad))

        net_sc2_weight = copy.deepcopy(net.sc2.weight.data)
        net_sc3_weight = copy.deepcopy(net.sc3.weight.data)

        net_state_dict = net.state_dict()

        # Sparse Coding
        copy_net        = copy.deepcopy(net)
        copy_state_dict = copy_net.state_dict()

        for name, param in copy_state_dict.items():
            if "weight" not in name:
                continue
            if "sc1" in name:
                continue

            if "sc2" in name:
                active_param = net_sc2_weight.mul(do_m1_grad_mask)
            if "sc3" in name:
                active_param = net_sc3_weight.mul(do_m2_grad_mask)

            nonzero_param_1d = active_param[active_param != 0]
            if nonzero_param_1d.size(0) == 0:
                break

            copy_param_1d = copy.deepcopy(nonzero_param_1d)
            S_set         = torch.arange(100, -1, -10)[1:]
            copy_param    = copy.deepcopy(active_param)
            S_loss        = []

            for S in S_set:
                param_mask        = s_mask(sparse_level=S.item(), param_matrix=copy_param,
                                           nonzero_param_1D=copy_param_1d, dtype=dtype)
                transformed_param = copy_param.mul(param_mask)
                copy_state_dict[name].copy_(transformed_param)
                copy_net.train()
                y_tmp    = copy_net(train_x)
                loss_tmp = binary_cross_entropy_for_imbalance(y_tmp, train_y)
                S_loss.append(loss_tmp)

            S_loss        = [tensor.detach() for tensor in S_loss]
            interp_S_loss = interp1d(S_set, S_loss, kind='cubic')
            interp_S_set  = torch.linspace(min(S_set), max(S_set), steps=100)
            interp_loss   = interp_S_loss(interp_S_set)
            optimal_S     = interp_S_set[np.argmin(interp_loss)]

            optimal_param_mask = s_mask(sparse_level=optimal_S.item(), param_matrix=copy_param,
                                        nonzero_param_1D=copy_param_1d, dtype=dtype)

            if "sc2" in name:
                final_optimal_param_mask = torch.where(
                    do_m1_grad_mask == 0, torch.ones_like(do_m1_grad_mask), optimal_param_mask)
                optimal_transformed_param = net_sc2_weight.mul(final_optimal_param_mask)
            if "sc3" in name:
                final_optimal_param_mask = torch.where(
                    do_m2_grad_mask == 0, torch.ones_like(do_m2_grad_mask), optimal_param_mask)
                optimal_transformed_param = net_sc3_weight.mul(final_optimal_param_mask)

            copy_state_dict[name].copy_(optimal_transformed_param)
            net_state_dict[name].copy_(optimal_transformed_param)

        if epoch == nEpochs - 1:
            net.train()
            train_pred = net(train_x)
            train_loss = binary_cross_entropy_for_imbalance(train_pred, train_y).view(1,)

            net.eval()
            eval_pred = net(eval_x)
            eval_loss = binary_cross_entropy_for_imbalance(eval_pred, eval_y).view(1,)

    # Save weight/bias matrices
    pd.DataFrame(net_state_dict['sc3.weight']).to_excel('sc3_weights.xlsx')
    pd.DataFrame(net_state_dict['sc3.bias']).to_excel('sc3_biases.xlsx')
    pd.DataFrame(net_state_dict['sc2.weight']).to_excel('sc2_weights.xlsx')
    pd.DataFrame(net_state_dict['sc2.bias']).to_excel('sc2_biases.xlsx')
    pd.DataFrame(net_state_dict['sc1.weight']).to_excel('sc1_weights.xlsx')
    pd.DataFrame(net_state_dict['sc1.bias']).to_excel('sc1_biases.xlsx')

    return train_pred, eval_pred, train_loss, eval_loss, net


def main():
    args   = parse_args()
    outdir = args.outdir

    # ── PASNet Settings ────────────────────────────────────────────────────────
    pathway_mask  = load_pathway(args.pathway_mask, dtype)
    In_Nodes      = pathway_mask.shape[1]
    Pathway_Nodes = pathway_mask.shape[0]
    Hidden_Nodes  = pathway_mask.shape[0]
    Out_Nodes     = 2

    # ── Initial Settings for Empirical Search ─────────────────────────────────
    Learning_Rates = [0.05, 0.01, 0.007, 0.005, 0.001, 0.0007, 0.0005, 0.0001]
    L2_Lambdas     = [3e-4, 5e-4, 7e-4, 1e-3, 3e-3, 5e-3]
    Dropout_Rates  = [0.8, 0.7]
    nEpochs        = 1

    # ── Load data for empirical / grid search ─────────────────────────────────
    x_train, y_train = load_data(args.train_data,    dtype, args.comparison)
    x_valid, y_valid = load_data(args.val_data_grid, dtype, args.comparison)

    opt_l2   = 0
    opt_lr   = 0
    opt_loss = torch.Tensor([float("Inf")])
    if torch.cuda.is_available():
        opt_loss = opt_loss.cuda()

    # ── Grid search for optimal hyperparameters ────────────────────────────────
    for lr in Learning_Rates:
        for l2 in L2_Lambdas:
            pred_tr, pred_val, loss_tr, loss_val, model = trainPASNet(
                x_train, y_train, x_valid, y_valid, pathway_mask,
                In_Nodes, Pathway_Nodes, Hidden_Nodes, Out_Nodes,
                lr, l2, nEpochs, Dropout_Rates, optimizer="Adam")

            if loss_val < opt_loss:
                opt_l2   = l2
                opt_lr   = lr
                opt_loss = loss_val

            print("L2: ", l2, "LR: ", lr, "Loss in Validation: ", loss_val)
            print("Optimal L2: ", opt_l2, "Optimal LR: ", opt_lr)
            auc_te = auc(y_valid, pred_val)
            f1_te  = f1(y_valid, pred_val)
            print("AUC in Test: ", auc_te, "F1 in Test: ", f1_te)

    # ── Retrain with optimal hyperparameters ──────────────────────────────────
    pathway_mask  = load_pathway(args.pathway_mask, dtype)
    In_Nodes      = pathway_mask.shape[1]
    Pathway_Nodes = pathway_mask.shape[0]
    Hidden_Nodes  = pathway_mask.shape[0]
    Out_Nodes     = 2

    nEpochs       = 1
    Dropout_Rates = [0.8, 0.7]

    N = 1  # number of repeated times
    K = 1  # number of folds

    test_auc = []
    test_f1  = []

    for replicate in range(N):
        for fold in range(K):
            print("replicate: ", replicate, "fold: ", fold)

            # FIX: pass args.comparison to every load_data() call
            x_train, y_train = load_data(args.train_data, dtype, args.comparison)
            x_eval,  y_eval  = load_data(args.val_data,   dtype, args.comparison)

            pred_train, pred_eval, loss_train, loss_eval, model = trainPASNet(
                x_train, y_train, x_eval, y_eval, pathway_mask,
                In_Nodes, Pathway_Nodes, Hidden_Nodes, Out_Nodes,
                opt_lr, opt_l2, nEpochs, Dropout_Rates, optimizer="Adam")

            if torch.cuda.is_available():
                pred_eval = pred_eval.cpu().detach()

            np.savetxt(
                "PASNet_pred_" + str(replicate) + "_" + str(fold) + ".txt",
                pred_eval.detach().numpy(), delimiter=","
            )

            auc_te = calc_auc(y_eval, pred_eval)
            f1_te  = f1(y_eval, pred_eval)
            print("AUC in Test: ", auc_te, "F1 in Test: ", f1_te)
            test_auc.append(auc_te)
            test_f1.append(f1_te)

    np.savetxt("PASNet_AUC.txt", test_auc, delimiter=",")
    np.savetxt("PASNet_F1.txt",  test_f1,  delimiter=",")

    # ── Save the model ─────────────────────────────────────────────────────────
    model_output_dir = os.path.join(outdir, "Output/GO")
    os.makedirs(model_output_dir, exist_ok=True)
    with open(os.path.join(model_output_dir, 'model.pkl'), 'wb') as file:
        pickle.dump(model, file)

    print("\n=== PASNet training completed successfully ===")


if __name__ == "__main__":
    main()
