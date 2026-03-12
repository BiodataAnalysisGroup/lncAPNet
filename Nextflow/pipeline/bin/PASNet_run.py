#!/usr/bin/env python
"""
PASNet: Pathway-Associated Sparse Deep Neural Network for Genomic Prediction

This script trains and evaluates a PASNet model for binary classification tasks
in bioinformatics, particularly for genomic prediction with pathway information.
"""

import argparse
import sys
import os
import numpy as np
import pandas as pd
import math
import copy
from scipy.interpolate import interp1d
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
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

# Set random seed for reproducibility
torch.manual_seed(0)
np.random.seed(0)


def vectorized_label(target, n_class):
    """Convert target(y) to one-hot encoding format (dummy variable)"""
    TARGET = np.array(target).reshape(-1)
    return np.eye(n_class)[TARGET]


def load_data(path, dtype):
    """Load data and convert it to a PyTorch tensor.
    
    Args:
        path: Path to input dataset (csv or excel file)
        dtype: Data type of tensor (e.g., torch.FloatTensor)
    
    Returns:
        X: PyTorch tensor of features
        Y: PyTorch tensor of labels (one-hot encoding)
    """
    # Detect file type and load accordingly
    if path.endswith('.csv'):
        data = pd.read_csv(path)
    elif path.endswith('.xlsx') or path.endswith('.xls'):
        data = pd.read_excel(path)
    else:
        raise ValueError(f"Unsupported file format: {path}")
    
    # Drop unnamed columns if present
    if 'Unnamed: 0' in data.columns:
        data = data.drop('Unnamed: 0', axis=1)
    
    x = data.drop(["Status"], axis=1).values
    y = data.loc[:, ["Status"]].values
    
    X = torch.from_numpy(x).type(dtype)
    Y = torch.from_numpy(vectorized_label(y, 2)).type(dtype)
    
    if torch.cuda.is_available():
        X = X.cuda()
        Y = Y.cuda()
    
    return X, Y


def load_pathway(path, dtype):
    """Load a bi-adjacency matrix of pathways and convert to PyTorch tensor.
    
    Args:
        path: Path to pathway matrix file (csv or excel)
        dtype: Data type of tensor
    
    Returns:
        PATHWAY_MASK: PyTorch tensor of the bi-adjacency matrix
    """
    if path.endswith('.csv'):
        pathway_mask = pd.read_csv(path, index_col=0)
    elif path.endswith('.xlsx') or path.endswith('.xls'):
        pathway_mask = pd.read_excel(path, index_col=0)
    else:
        raise ValueError(f"Unsupported file format: {path}")
    
    PATHWAY_MASK = torch.from_numpy(pathway_mask.values.astype(np.float32)).type(dtype)
    
    if torch.cuda.is_available():
        PATHWAY_MASK = PATHWAY_MASK.cuda()
    
    return PATHWAY_MASK


def auc(y_true, y_pred):
    """Calculate AUC score"""
    if torch.cuda.is_available():
        y_true = y_true.cpu().detach()
        y_pred = y_pred.cpu().detach()
    
    auc_score = roc_auc_score(y_true.detach().numpy(), y_pred.detach().numpy())
    return auc_score


def f1(y_true, y_pred):
    """Calculate F1 score"""
    # Convert one-hot encoding to integer
    y = torch.argmax(y_true, dim=1)
    # Estimated targets (either 0 or 1)
    pred = torch.argmax(y_pred, dim=1)
    
    if torch.cuda.is_available():
        y = y.cpu().detach()
        pred = pred.cpu().detach()
    
    f1_score_val = f1_score(y.detach().numpy(), pred.detach().numpy())
    return f1_score_val


def bce_for_one_class(predict, target, lts=False):
    """Calculate cross entropy in average for samples belonging to the same class.
    
    Args:
        predict: Predicted values
        target: True labels
        lts: If False, non-severe samples; if True, severe samples
    """
    lts_idx = torch.argmax(target, dim=1)
    idx = 1 if lts else 0
    y = target[lts_idx == idx]
    pred = predict[lts_idx == idx]
    
    cost = F.binary_cross_entropy(pred, y)
    return cost


def binary_cross_entropy_for_imbalance(predict, target):
    """Calculate cross entropy for imbalanced data in binary classification"""
    total_cost = bce_for_one_class(predict, target, lts=True) + \
                 bce_for_one_class(predict, target, lts=False)
    return total_cost


def dropout_mask(n_node, drop_p):
    """Construct a binary matrix to randomly drop nodes in a layer.
    
    Args:
        n_node: Number of nodes in the layer
        drop_p: Probability that a node is dropped
    
    Returns:
        mask: Binary matrix (1=keep, 0=drop)
    """
    keep_p = 1.0 - drop_p
    mask = torch.Tensor(np.random.binomial(1, keep_p, size=n_node))
    
    if torch.cuda.is_available():
        mask = mask.cuda()
    
    return mask


def s_mask(sparse_level, param_matrix, nonzero_param_1D, dtype):
    """Construct a binary matrix w.r.t. a sparsity level of weights.
    
    Args:
        sparse_level: Percentage value in [0, 100) for proportion of weights to drop
        param_matrix: Weight matrix for entire network
        nonzero_param_1D: 1D of non-zero param_matrix
        dtype: Data type of tensor
    
    Returns:
        param_mask: Binary matrix (1=keep, 0=drop)
    """
    # Take absolute values
    non_neg_param_1D = torch.abs(nonzero_param_1D)
    num_param = nonzero_param_1D.size(0)
    top_k = math.ceil(num_param * (100 - sparse_level) * 0.01)
    
    sorted_non_neg_param_1D, indices = torch.topk(non_neg_param_1D, top_k)
    param_mask = torch.abs(param_matrix) > sorted_non_neg_param_1D.min()
    param_mask = param_mask.type(dtype)
    
    if torch.cuda.is_available():
        param_mask = param_mask.cuda()
    
    return param_mask


class PASNet(nn.Module):
    """Pathway-Associated Sparse Neural Network"""
    
    def __init__(self, In_Nodes, Pathway_Nodes, Hidden_Nodes, Out_Nodes, Pathway_Mask):
        super(PASNet, self).__init__()
        self.sigmoid = nn.Sigmoid()
        self.softmax = nn.Softmax(dim=1)
        self.pathway_mask = Pathway_Mask
        
        # Gene layer --> pathway layer
        self.sc1 = nn.Linear(In_Nodes, Pathway_Nodes)
        # Pathway layer --> hidden layer
        self.sc2 = nn.Linear(Pathway_Nodes, Hidden_Nodes)
        # Hidden layer --> output layer
        self.sc3 = nn.Linear(Hidden_Nodes, Out_Nodes)
        
        # Dropout masks for sub-network
        self.do_m1 = torch.ones(Pathway_Nodes)
        self.do_m2 = torch.ones(Hidden_Nodes)
        
        if torch.cuda.is_available():
            self.do_m1 = self.do_m1.cuda()
            self.do_m2 = self.do_m2.cuda()
    
    def forward(self, x):
        # Force connections between gene layer and pathway layer
        self.sc1.weight.data = self.sc1.weight.data.mul(self.pathway_mask)
        x = self.sigmoid(self.sc1(x))
        
        if self.training:
            x = x.mul(self.do_m1)
        
        x = self.sigmoid(self.sc2(x))
        
        if self.training:
            x = x.mul(self.do_m2)
        
        x = self.softmax(self.sc3(x))
        return x


def trainPASNet(train_x, train_y, eval_x, eval_y, pathway_mask,
                In_Nodes, Pathway_Nodes, Hidden_Nodes, Out_Nodes,
                Learning_Rate, L2_Lambda, nEpochs, Dropout_Rates,
                optimizer_name="Adam", dtype=torch.FloatTensor):
    """Train PASNet model"""
    
    net = PASNet(In_Nodes, Pathway_Nodes, Hidden_Nodes, Out_Nodes, pathway_mask)
    
    if torch.cuda.is_available():
        net.cuda()
    
    # Initialize optimizer
    if optimizer_name == "SGD":
        opt = optim.SGD(net.parameters(), lr=Learning_Rate, weight_decay=L2_Lambda)
    else:
        opt = optim.Adam(net.parameters(), lr=Learning_Rate, weight_decay=L2_Lambda)
    
    for epoch in range(nEpochs):
        net.train()
        opt.zero_grad()
        
        # Randomize dropout masks
        net.do_m1 = dropout_mask(Pathway_Nodes, Dropout_Rates[0])
        net.do_m2 = dropout_mask(Hidden_Nodes, Dropout_Rates[1])
        
        pred = net(train_x)
        loss = binary_cross_entropy_for_imbalance(pred, train_y)
        loss.backward()
        opt.step()
        
        net.sc1.weight.data = net.sc1.weight.data.mul(net.pathway_mask)
        
        # Obtain sub-network's connections
        do_m1_grad = copy.deepcopy(net.sc2.weight._grad.data)
        do_m2_grad = copy.deepcopy(net.sc3.weight._grad.data)
        do_m1_grad_mask = torch.where(do_m1_grad == 0, do_m1_grad, torch.ones_like(do_m1_grad))
        do_m2_grad_mask = torch.where(do_m2_grad == 0, do_m2_grad, torch.ones_like(do_m2_grad))
        
        # Copy weights
        net_sc2_weight = copy.deepcopy(net.sc2.weight.data)
        net_sc3_weight = copy.deepcopy(net.sc3.weight.data)
        
        net_state_dict = net.state_dict()
        
        # Sparse Coding
        copy_net = copy.deepcopy(net)
        copy_state_dict = copy_net.state_dict()
        
        for name, param in copy_state_dict.items():
            if not "weight" in name:
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
            S_set = torch.arange(100, -1, -10)[1:]
            copy_param = copy.deepcopy(active_param)
            S_loss = []
            
            for S in S_set:
                param_mask = s_mask(sparse_level=S.item(), param_matrix=copy_param,
                                   nonzero_param_1D=copy_param_1d, dtype=dtype)
                transformed_param = copy_param.mul(param_mask)
                copy_state_dict[name].copy_(transformed_param)
                copy_net.train()
                y_tmp = copy_net(train_x)
                loss_tmp = binary_cross_entropy_for_imbalance(y_tmp, train_y)
                S_loss.append(loss_tmp)
            
            # Apply cubic interpolation
            S_loss = [tensor.detach() for tensor in S_loss]
            interp_S_loss = interp1d(S_set, S_loss, kind='cubic')
            interp_S_set = torch.linspace(min(S_set), max(S_set), steps=100)
            interp_loss = interp_S_loss(interp_S_set)
            optimal_S = interp_S_set[np.argmin(interp_loss)]
            optimal_param_mask = s_mask(sparse_level=optimal_S.item(), param_matrix=copy_param,
                                       nonzero_param_1D=copy_param_1d, dtype=dtype)
            
            if "sc2" in name:
                final_optimal_param_mask = torch.where(do_m1_grad_mask == 0,
                                                       torch.ones_like(do_m1_grad_mask),
                                                       optimal_param_mask)
                optimal_transformed_param = net_sc2_weight.mul(final_optimal_param_mask)
            if "sc3" in name:
                final_optimal_param_mask = torch.where(do_m2_grad_mask == 0,
                                                       torch.ones_like(do_m2_grad_mask),
                                                       optimal_param_mask)
                optimal_transformed_param = net_sc3_weight.mul(final_optimal_param_mask)
            
            copy_state_dict[name].copy_(optimal_transformed_param)
            net_state_dict[name].copy_(optimal_transformed_param)
        
        # Print progress every 1000 epochs
        if (epoch + 1) % 1000 == 0:
            net.eval()
            with torch.no_grad():
                eval_pred = net(eval_x)
                eval_loss = binary_cross_entropy_for_imbalance(eval_pred, eval_y)
                eval_auc = auc(eval_y, eval_pred)
                eval_f1 = f1(eval_y, eval_pred)
            print(f"Epoch {epoch+1}/{nEpochs}: Eval Loss={eval_loss:.4f}, AUC={eval_auc:.4f}, F1={eval_f1:.4f}")
    
    # Final evaluation
    net.train()
    train_pred = net(train_x)
    train_loss = binary_cross_entropy_for_imbalance(train_pred, train_y).view(1,)
    
    net.eval()
    eval_pred = net(eval_x)
    eval_loss = binary_cross_entropy_for_imbalance(eval_pred, eval_y).view(1,)
    
    return train_pred, eval_pred, train_loss, eval_loss, net


def save_model_weights(model, output_dir):
    """Save model weights to Excel files"""
    os.makedirs(output_dir, exist_ok=True)
    
    state_dict = model.state_dict()
    
    # Save sc3 weights and biases
    pd.DataFrame(state_dict['sc3.weight'].cpu().numpy()).to_excel(
        os.path.join(output_dir, 'sc3_weights.xlsx'))
    pd.DataFrame(state_dict['sc3.bias'].cpu().numpy()).to_excel(
        os.path.join(output_dir, 'sc3_biases.xlsx'))
    
    # Save sc2 weights and biases
    pd.DataFrame(state_dict['sc2.weight'].cpu().numpy()).to_excel(
        os.path.join(output_dir, 'sc2_weights.xlsx'))
    pd.DataFrame(state_dict['sc2.bias'].cpu().numpy()).to_excel(
        os.path.join(output_dir, 'sc2_biases.xlsx'))
    
    # Save sc1 weights and biases
    pd.DataFrame(state_dict['sc1.weight'].cpu().numpy()).to_excel(
        os.path.join(output_dir, 'sc1_weights.xlsx'))
    pd.DataFrame(state_dict['sc1.bias'].cpu().numpy()).to_excel(
        os.path.join(output_dir, 'sc1_biases.xlsx'))
    
    print(f"Model weights saved to {output_dir}")


def plot_roc_curve(y_true, y_pred, output_path):
    """Plot and save ROC curve"""
    # Convert one-hot y_true to class indices
    if len(y_true.shape) > 1 and y_true.shape[1] > 1:
        y_true_idx = y_true.argmax(dim=1).cpu().numpy()
    else:
        y_true_idx = y_true.cpu().numpy()
    
    # Convert logits to softmax probabilities for class 1
    probs = torch.softmax(y_pred, dim=1)[:, 1].cpu().numpy()
    
    # Compute ROC curve and AUC
    fpr, tpr, _ = roc_curve(y_true_idx, probs)
    roc_auc = sklearn_auc(fpr, tpr)
    
    # Plot ROC curve
    plt.figure()
    plt.plot(fpr, tpr, label=f"ROC curve (AUC = {roc_auc:.2f})")
    plt.plot([0, 1], [0, 1], 'k--', label="Random")
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title("ROC Curve")
    plt.legend(loc="lower right")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
    
    # Save fpr and tpr data
    roc_data_path = output_path.replace('.png', '_data.txt')
    np.savetxt(roc_data_path, np.column_stack((fpr, tpr)),
               fmt="%.6f", header="FPR\tTPR", delimiter="\t")
    
    print(f"ROC curve saved to {output_path}")
    print(f"ROC data saved to {roc_data_path}")
    
    return roc_auc


def generate_shap_analysis(model, x_train, x_eval, feature_names, output_dir, max_display=50):
    """Generate SHAP analysis and feature importance"""
    os.makedirs(output_dir, exist_ok=True)
    
    print("Generating SHAP values (this may take a while)...")
    explainer = shap.DeepExplainer(model, x_train)
    shap_values = explainer.shap_values(x_eval)
    
    # Calculate average absolute SHAP values for feature importance
    vals = np.abs(shap_values).mean(axis=1).mean(axis=0)
    feature_importance = pd.DataFrame(
        list(zip(feature_names, vals)),
        columns=['Feature', 'Importance']
    )
    
    # Sort features by importance
    feature_importance.sort_values(by='Importance', ascending=False, inplace=True)
    
    # Save feature importance table
    importance_path = os.path.join(output_dir, 'feature_importance.csv')
    feature_importance.to_csv(importance_path, index=False)
    print(f"Feature importance saved to {importance_path}")
    
    # Generate and save summary plot
    shap_plot_path = os.path.join(output_dir, 'shap_summary.png')
    fig = shap.summary_plot(shap_values, x_eval.cpu().numpy(),
                           feature_names=feature_names,
                           max_display=max_display, show=False)
    plt.savefig(shap_plot_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"SHAP summary plot saved to {shap_plot_path}")


def main():
    parser = argparse.ArgumentParser(
        description='Train and evaluate PASNet model for genomic prediction',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Input/Output arguments
    parser.add_argument('--train-data', type=str, required=True,
                       help='Path to training data file (csv or excel)')
    parser.add_argument('--val-data', type=str, required=True,
                       help='Path to validation data file (csv or excel)')
    parser.add_argument('--pathway-mask', type=str, required=True,
                       help='Path to pathway mask file (csv or excel)')
    parser.add_argument('--output-dir', type=str, default='./output',
                       help='Directory to save output files')
    
    # Network architecture arguments
    parser.add_argument('--in-nodes', type=int, required=True,
                       help='Number of input nodes (genes)')
    parser.add_argument('--pathway-nodes', type=int, required=True,
                       help='Number of pathway nodes')
    parser.add_argument('--hidden-nodes', type=int, required=True,
                       help='Number of hidden nodes')
    parser.add_argument('--out-nodes', type=int, default=2,
                       help='Number of output nodes (classes)')
    
    # Training hyperparameters
    parser.add_argument('--learning-rate', type=float, default=0.01,
                       help='Learning rate')
    parser.add_argument('--l2-lambda', type=float, default=0.0003,
                       help='L2 regularization parameter')
    parser.add_argument('--epochs', type=int, default=300,
                       help='Number of training epochs')
    parser.add_argument('--dropout-rates', type=float, nargs=2, default=[0.8, 0.7],
                       help='Dropout rates for pathway and hidden layers')
    parser.add_argument('--optimizer', type=str, default='Adam', choices=['Adam', 'SGD'],
                       help='Optimizer to use')
    
    # Experiment settings
    parser.add_argument('--n-replicates', type=int, default=1,
                       help='Number of replicates')
    parser.add_argument('--n-folds', type=int, default=1,
                       help='Number of folds')
    
    # Analysis options
    parser.add_argument('--save-weights', action='store_true',
                       help='Save model weights to Excel files')
    parser.add_argument('--save-predictions', action='store_true',
                       help='Save predictions to text files')
    parser.add_argument('--plot-roc', action='store_true',
                       help='Generate ROC curve plot')
    parser.add_argument('--shap-analysis', action='store_true',
                       help='Perform SHAP analysis for feature importance')
    parser.add_argument('--shap-max-display', type=int, default=50,
                       help='Maximum features to display in SHAP plot')
    
    # Device arguments
    parser.add_argument('--no-cuda', action='store_true',
                       help='Disable CUDA even if available')
    
    args = parser.parse_args()
    
    # Set device
    if args.no_cuda:
        torch.cuda.is_available = lambda: False
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Set data type
    dtype = torch.FloatTensor
    
    # Load data
    print("Loading data...")
    pathway_mask = load_pathway(args.pathway_mask, dtype)
    
    # Store results
    test_auc = []
    test_f1 = []
    
    for replicate in range(args.n_replicates):
        for fold in range(args.n_folds):
            print(f"\n{'='*60}")
            print(f"Replicate: {replicate+1}/{args.n_replicates}, Fold: {fold+1}/{args.n_folds}")
            print(f"{'='*60}")
            
            x_train, y_train = load_data(args.train_data, dtype)
            x_eval, y_eval = load_data(args.val_data, dtype)
            
            print(f"Training data shape: {x_train.shape}")
            print(f"Validation data shape: {x_eval.shape}")
            print(f"Pathway mask shape: {pathway_mask.shape}")
            
            # Train model
            print(f"\nTraining PASNet for {args.epochs} epochs...")
            pred_train, pred_eval, loss_train, loss_eval, model = trainPASNet(
                x_train, y_train, x_eval, y_eval, pathway_mask,
                args.in_nodes, args.pathway_nodes, args.hidden_nodes, args.out_nodes,
                args.learning_rate, args.l2_lambda, args.epochs, args.dropout_rates,
                args.optimizer, dtype
            )
            
            # Calculate metrics
            if torch.cuda.is_available():
                pred_eval = pred_eval.cpu().detach()
            
            auc_te = auc(y_eval, pred_eval)
            f1_te = f1(y_eval, pred_eval)
            
            print(f"\n{'='*60}")
            print(f"Results - Replicate {replicate+1}, Fold {fold+1}")
            print(f"{'='*60}")
            print(f"Test AUC: {auc_te:.4f}")
            print(f"Test F1:  {f1_te:.4f}")
            
            test_auc.append(auc_te)
            test_f1.append(f1_te)
            
            # Save predictions
            if args.save_predictions:
                pred_path = os.path.join(args.output_dir,
                                        f'predictions_rep{replicate}_fold{fold}.txt')
                np.savetxt(pred_path, pred_eval.detach().numpy(), delimiter=',')
                print(f"Predictions saved to {pred_path}")
            
            # Save model weights
            if args.save_weights:
                weights_dir = os.path.join(args.output_dir, f'weights_rep{replicate}_fold{fold}')
                save_model_weights(model, weights_dir)
            
            # Save model
            model_path = os.path.join(args.output_dir, f'model_rep{replicate}_fold{fold}.pkl')
            with open(model_path, 'wb') as f:
                pickle.dump(model, f)
            print(f"Model saved to {model_path}")
            
            # Plot ROC curve
            if args.plot_roc:
                roc_path = os.path.join(args.output_dir,
                                       f'roc_curve_rep{replicate}_fold{fold}.png')
                roc_auc = plot_roc_curve(y_eval, pred_eval, roc_path)
            
            # SHAP analysis
            if args.shap_analysis:
                # Load feature names
                if args.train_data.endswith('.csv'):
                    train_df = pd.read_csv(args.train_data)
                else:
                    train_df = pd.read_excel(args.train_data)
                
                if 'Unnamed: 0' in train_df.columns:
                    train_df = train_df.drop('Unnamed: 0', axis=1)
                train_df = train_df.drop(['Status'], axis=1)
                feature_names = train_df.columns.tolist()
                
                shap_dir = os.path.join(args.output_dir, f'shap_rep{replicate}_fold{fold}')
                generate_shap_analysis(model, x_train, x_eval, feature_names,
                                     shap_dir, args.shap_max_display)
    
    # Save summary results
    auc_path = os.path.join(args.output_dir, 'summary_AUC.txt')
    f1_path = os.path.join(args.output_dir, 'summary_F1.txt')
    np.savetxt(auc_path, test_auc, delimiter=',')
    np.savetxt(f1_path, test_f1, delimiter=',')
    
    print(f"\n{'='*60}")
    print("Summary Statistics")
    print(f"{'='*60}")
    print(f"Mean AUC: {np.mean(test_auc):.4f} ± {np.std(test_auc):.4f}")
    print(f"Mean F1:  {np.mean(test_f1):.4f} ± {np.std(test_f1):.4f}")
    print(f"\nResults saved to {args.output_dir}")


if __name__ == '__main__':
    main()