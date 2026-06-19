import math
import copy
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from .metrics import binary_cross_entropy_for_imbalance # Internal cross-module import

class PASNet(nn.Module):
    def __init__(self, In_Nodes, Pathway_Nodes, Hidden_Nodes, Out_Nodes, Pathway_Mask):
        super(PASNet, self).__init__()
        self.sigmoid      = nn.Sigmoid()
        self.softmax      = nn.Softmax(dim=1)
        self.pathway_mask = Pathway_Mask

        self.sc1 = nn.Linear(In_Nodes,      Pathway_Nodes)
        self.sc2 = nn.Linear(Pathway_Nodes, Hidden_Nodes)
        self.sc3 = nn.Linear(Hidden_Nodes,  Out_Nodes)

        self.do_m1 = torch.ones(Pathway_Nodes)
        self.do_m2 = torch.ones(Hidden_Nodes)

        if torch.cuda.is_available():
            self.do_m1 = self.do_m1.cuda()
            self.do_m2 = self.do_m2.cuda()

    def forward(self, x):
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
    keep_p = 1.0 - drop_p
    mask   = torch.Tensor(np.random.binomial(1, keep_p, size=n_node))
    if torch.cuda.is_available():
        mask = mask.cuda()
    return mask

def s_mask(sparse_level, param_matrix, nonzero_param_1D, dtype):
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
