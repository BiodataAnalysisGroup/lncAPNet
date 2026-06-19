import torch
import torch.nn.functional as F
from sklearn.metrics import roc_auc_score, f1_score

def auc(y_true, y_pred):
    if torch.cuda.is_available():
        y_true = y_true.cpu().detach()
        y_pred = y_pred.cpu().detach()
    score = roc_auc_score(y_true.detach().numpy(), y_pred.detach().numpy())
    return score

def calc_auc(y_true, y_pred):
    return auc(y_true, y_pred)

def f1(y_true, y_pred):
    y = torch.argmax(y_true, dim=1)
    pred = torch.argmax(y_pred, dim=1)
    if torch.cuda.is_available():
        y = y.cpu().detach()
        pred = pred.cpu().detach()
    score = f1_score(y.detach().numpy(), pred.detach().numpy())
    return score

def bce_for_one_class(predict, target, lts=False):
    lts_idx = torch.argmax(target, dim=1)
    idx = 1 if lts else 0
    y = target[lts_idx == idx]
    pred = predict[lts_idx == idx]
    cost = F.binary_cross_entropy(pred, y)
    return cost

def binary_cross_entropy_for_imbalance(predict, target):
    total_cost = (bce_for_one_class(predict, target, lts=True)
                  + bce_for_one_class(predict, target, lts=False))
    return total_cost
