import torch, sys, os
from sklearn.metrics import average_precision_score
from sklearn import metrics
import pickle, itertools, math, copy
import numpy as np

from operator import itemgetter, add
from sklearn.metrics import confusion_matrix
import torch.nn.functional as F

import torch.nn as nn
from tqdm import tqdm
from torch.optim.lr_scheduler import ExponentialLR

from sklearn.metrics import balanced_accuracy_score, roc_auc_score, precision_recall_curve, \
    average_precision_score, f1_score, auc, recall_score, precision_score

from torch.autograd import Variable
from sklearn.preprocessing import StandardScaler
std_scaler = StandardScaler()

# Setting up the device
os.environ["CUDA_DEVICE_ORDER"]="PCI_BUS_ID"   
os.environ["CUDA_VISIBLE_DEVICES"]="0"
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")


class CNN_model(nn.Module):
    def __init__(self, rn, pns, dropout, pockn):
        super(CNN_model, self).__init__()
        self.rn = rn  # Input feature size
        self.pns = pns  # Pooling size
        self.pockn = pockn  # Number of protein pockets
        # Define APM processing pipelines for proteins and ligands
        self.ptn_apm = nn.Sequential(
            nn.Conv1d(rn, rn // 2, 2, 1),
            nn.Conv1d(rn // 2, rn // 4, 2, 1),
            nn.BatchNorm1d(rn // 4),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.AdaptiveMaxPool1d(pns)
        )
        self.lgd_apm = nn.Sequential(
            nn.Conv1d(rn, rn // 2, 2, 1),
            nn.Conv1d(rn // 2, rn // 4, 2, 1),
            nn.BatchNorm1d(rn // 4),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.AdaptiveMaxPool1d(pns)
        )
        self.bilstm = nn.LSTM(rn // 4 * pns, rn // 4 * pns, bidirectional=True, dropout=dropout)
        self.fc_in = nn.Linear(rn // 4 * pns * pockn * 2, rn // 4 * pns * pockn)
        self.fc_out = nn.Linear(rn // 4 * pns * pockn, 1)
        self.attention = MultiHeadAttention(rn // 4 * pns * 2, rn // 4 * pns * 2, 2)

    def forward(self, g):
        feature_ptn = self.ptn_apm(g[1])  # Process protein APM
        feature_lgd = self.lgd_apm(g[0])  # Process ligand APM
        sequence = torch.cat((feature_lgd, feature_ptn), dim=0).view(1, -1, self.rn // 4 * self.pns)

        mask = torch.eye(self.pockn, dtype=torch.uint8).view(1, self.pockn, self.pockn).cuda()
        mask[0, sequence.size()[1]:self.pockn, :] = 0
        mask[0, :, sequence.size()[1]:self.pockn] = 0
        mask[0, :, sequence.size()[1] - 1] = 1
        mask[0, sequence.size()[1] - 1, :] = 1
        mask[0,  sequence.size()[1] - 1,  sequence.size()[1] - 1] = 0
        sequence = F.pad(input=sequence, pad=(0, 0, 0, self.pockn - sequence.size()[1]), mode='constant', value=0)
        sequence = sequence.permute(1, 0, 2)
        
        h_0 = Variable(torch.zeros(2, 1, int(self.rn/4)*self.pns).cuda())
        c_0 = Variable(torch.zeros(2, 1, int(self.rn/4)*self.pns).cuda())
        
        output, _ = self.bilstm(sequence, (h_0, c_0))
        
        output = output.permute(1, 0,  2)
        out,att = self.attention(output, mask=mask, return_attention=True)  # Apply multi-head attention
        out = F.relu(self.fc_in(out.view(-1, out.size()[1] * out.size()[2])))  # Fully connected layer
        out = torch.sigmoid(self.fc_out(out))  # Sigmoid for binary classification
        return out


# Compute scaled dot-product attention
def scaled_dot_product(q, k, v, mask=None):
    d_k = q.size()[-1]
    attn_logits = torch.matmul(q, k.transpose(-2, -1))
    attn_logits = attn_logits / math.sqrt(d_k)
    if mask is not None:
        attn_logits = attn_logits.masked_fill(mask == 0, -9e15)
    attention = F.softmax(attn_logits, dim=-1)
    values = torch.matmul(attention, v)
    return values, attention


# Multi-head attention mechanism
class MultiHeadAttention(nn.Module):
    def __init__(self, input_dim, embed_dim, num_heads):
        super().__init__()
        assert embed_dim % num_heads == 0, "Embedding dimension must be 0 modulo number of heads."
        self.embed_dim = embed_dim
        self.num_heads = num_heads
        self.head_dim = embed_dim // num_heads
        # Stack all weight matrices 1...h together for efficiency
        # Note that in many implementations you see "bias=False" which is optional
        self.qkv_proj = nn.Linear(input_dim, 3 * embed_dim)
        self.o_proj = nn.Linear(embed_dim, embed_dim)
        self._reset_parameters()
    def _reset_parameters(self):
        # Original Transformer initialization, see PyTorch documentation
        nn.init.xavier_uniform_(self.qkv_proj.weight)
        self.qkv_proj.bias.data.fill_(0)
        nn.init.xavier_uniform_(self.o_proj.weight)
        self.o_proj.bias.data.fill_(0)
    def forward(self, x, mask=None, return_attention=False):
        batch_size, seq_length, embed_dim = x.size()
        qkv = self.qkv_proj(x)
        # Separate Q, K, V from linear output
        qkv = qkv.reshape(batch_size, seq_length, self.num_heads, 3 * self.head_dim)
        qkv = qkv.permute(0, 2, 1, 3)  # [Batch, Head, SeqLen, Dims]
        q, k, v = qkv.chunk(3, dim=-1)
        # Determine value outputs
        values, attention = scaled_dot_product(q, k, v, mask=mask)
        values = values.permute(0, 2, 1, 3)  # [Batch, SeqLen, Head, Dims]
        values = values.reshape(batch_size, seq_length, embed_dim)
        o = self.o_proj(values)
        if return_attention:
            return o, attention
        else:
            return o



def fwd_pass(X, y, model_t, optimizer, criterion, train=False):
    if train:
        optimizer.zero_grad()
    out = []
    for item in X:
        x = [0, 0]
        x[0] = item[0].to(device)
        x[1] = item[1].to(device)
        out.append(model_t(x))
        del x
    out = torch.stack(out, 0).view(-1, 1).to(device)
    y = torch.Tensor(y).view(-1, 1).to(device)
    loss = criterion(out, y)
    matches = [torch.round(i) == torch.round(j) for i, j in zip(out, y)]
    acc = matches.count(True) / len(matches)
    if train:
        loss.backward()
        optimizer.step()
    return acc, loss, out


def get_roce(predList, targetList, roceRate):
    p = sum(targetList)
    n = len(targetList) - p
    predList = [[index, x] for index, x in enumerate(predList)]
    predList = sorted(predList, key=lambda x:x[1], reverse=True)
    tp1 = 0
    fp1 = 0
    maxIndexs = []
    for x in predList:
        if(targetList[x[0]] == 1):
            tp1 += 1
        else:
            fp1 += 1
            if(fp1>((roceRate*n)/100)):
                break
    roce = (tp1*n)/(p*fp1)
    return roce


def test_func(model_f, y_label, X_test_f):
    y_pred = []
    y_label = torch.Tensor(y_label)
    print("Testing:")
    print("-------------------")
    with tqdm(range(0, len(X_test_f), 1)) as tepoch:
        for i in tepoch:
            with torch.no_grad():
                x = [0, 0]
                x[0] = X_test_f[i][0].to(device)
                x[1] = X_test_f[i][1].to(device)
                y_pred.append(model_f(x).cpu())
    y_pred = torch.cat(y_pred, dim=0)
    y_pred_c = [round(i.item()) for i in y_pred]
    auroc = str(roc_auc_score(y_label, y_pred))
    print("AUROC: " + auroc, end=" ")
    return auroc



def setting_dataset(cmpdf, ptndf, label_df, tag, pockn, pair_typeN):
    subsets = label_df[label_df.iloc[:,3]==tag]
    subsets = subsets[subsets.iloc[:,0].isin(cmpdf.index)&subsets.iloc[:,1].isin(ptndf.index)]
    apm = np.array(cmpdf.loc[subsets.iloc[:,0],:])        
    apm = torch.tensor(apm, dtype = torch.float).view(len(subsets),1,pair_typeN,10)
    ppock = []
    for pds in list(subsets.iloc[:,1]):
        pcks = ptndf.loc[pds,:]
        if len(pcks)>(pockn-1):
            pcks = pcks.sample(frac = 1, random_state = 7)[0:(pockn-1)]
        ppock.append(torch.tensor(np.array(pcks).reshape([-1,pair_typeN,10]), dtype = torch.float))
    X = [[apm[i], ppock[i]] for i in range(len(apm))]
    y = list(subsets.iloc[:,2])
    return X, y
