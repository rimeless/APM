
import sys
import torch
src_dir = '/spstorage/USERS/gina/source' 
sys.path.append(src_dir)
import os
os.environ["CUDA_LAUNCH_BLOCKING"]="1"   
os.environ["CUDA_DEVICE_ORDER"]="PCI_BUS_ID"   
os.environ["CUDA_VISIBLE_DEVICES"]="5"
device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

from sklearn.metrics import average_precision_score
from sklearn import metrics
import pickle

from operator import itemgetter, add
import itertools
from sklearn.metrics import confusion_matrix
import torch.nn.functional as F


from APM_models import *
from models import *

from sklearn.preprocessing import StandardScaler
std_scaler = StandardScaler()

import math

import torch.nn as nn
from tqdm import tqdm


from torch.optim.lr_scheduler import ExponentialLR

from sklearn.metrics import balanced_accuracy_score, roc_auc_score, precision_recall_curve, \
    average_precision_score, f1_score, auc, recall_score, precision_score


from torch.autograd import Variable


def fwd_pass(X, y, model_t, train=False):
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
    # roce1 = get_roce(y_pred, y_label, 0.5)
    # roce2 = get_roce(y_pred, y_label, 1)
    # roce3 = get_roce(y_pred, y_label, 2)
    # roce4 = get_roce(y_pred, y_label, 5)
    auroc = str(roc_auc_score(y_label, y_pred))
    print("AUROC: " + auroc, end=" ")
    # print("PRAUC: " + str(average_precision_score(y_label, y_pred)), end=" ")
    # print("F1 Score: " + str(f1_score(y_label, y_pred_c)), end=" ")
    # print("Precision Score:" + str(precision_score(y_label, y_pred_c)), end=" ")
    # print("Recall Score:" + str(recall_score(y_label, y_pred_c)), end=" ")
    # print("Balanced Accuracy Score " + str(balanced_accuracy_score(y_label, y_pred_c)), end=" ")
    # print("0.5 re Score " + str(roce1), end=" ")
    # print("1 re Score " + str(roce2), end=" ")
    # print("2 re Score " + str(roce3), end=" ")
    # print("5 re Score " + str(roce4), end=" ")
    # print("-------------------")
    return auroc




def setting_dataset(lgd, ptn, opt, pockn, df, nm, rn, std_scaler):
    subsets = df[df.tag==opt[0]]
    if len(opt)>1:
        subsets = subsets[subsets.type==opt[1]]
    if len(opt)>2:
        subsets = subsets[subsets.rigidity==opt[2]]
    if nm:
        if opt[0]=='train':
            apm = std_scaler.fit_transform(lgd.loc[subsets.cid,:]) 
        else:
            apm = std_scaler.transform(lgd.loc[subsets.cid,:]) 
    else:
        apm = np.array(lgd.loc[subsets.cid,:])        
    apm = torch.tensor(apm, dtype = torch.float).view(len(subsets),1,rn,10)
    ppock = []
    for pds in list(subsets.pdb_id):
        pcks = pock.loc[pds,:]
        if (len(pcks)>(pockn-1))&(len(pcks)!=rn*10):
            pcks = pcks.sample(frac = 1, random_state = 7)[0:(pockn-1)]
        ppock.append(torch.tensor(np.array(pcks).reshape([-1,rn,10]), dtype = torch.float))
    X = [[apm[i], ppock[i]] for i in range(len(apm))]
    y = list(subsets.label)
    return X, y, std_scaler


class CNN_model(nn.Module):
    def __init__(self, rn, pns, dropout, pockn):
        super(CNN_model, self).__init__()
        self.rn = rn
        self.pns = pns 
        self.pockn = pockn
        self.ptn_apm = nn.Sequential(
            nn.Conv1d(in_channels = rn, out_channels = int(rn/2), kernel_size = 2, stride = 1, padding = 0, dilation =1, groups= 1, bias = True, padding_mode= 'zeros'),
            nn.Conv1d(in_channels = int(rn/2), out_channels = int(rn/4), kernel_size = 2, stride = 1, padding = 0, dilation =1, groups= 1, bias = True, padding_mode= 'zeros'),
            nn.BatchNorm1d(int(rn/4), eps=1e-05, momentum=0.1, affine=True, track_running_stats=True),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.AdaptiveMaxPool1d(pns)
        )
        self.lgd_apm = nn.Sequential(
            nn.Conv1d(in_channels = rn, out_channels = int(rn/2), kernel_size = 2, stride = 1, padding = 0, dilation =1, groups= 1, bias = True, padding_mode= 'zeros'),
            nn.Conv1d(in_channels = int(rn/2), out_channels = int(rn/4), kernel_size = 2, stride = 1, padding = 0, dilation =1, groups= 1, bias = True, padding_mode= 'zeros'),
            nn.BatchNorm1d(int(rn/4), eps=1e-05, momentum=0.1, affine=True, track_running_stats=True),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.AdaptiveMaxPool1d(pns)
        )
        self.bilstm = nn.LSTM(int(rn/4)*pns, int(rn/4)*pns, num_layers=1, bidirectional=True, dropout=dropout)
        self.fc_in = nn.Linear(int(rn/4)*pns*pockn*2, int(rn/4)*pns*pockn) #1922
        self.fc_out = nn.Linear(int(rn/4)*pns*pockn, 1)
        self.attention = MultiHeadAttention(int(rn/4)*pns*2, int(rn/4)*pns*2, 2)
    #    self.W_s1 = nn.Linear(60, 45) #62
    #    self.W_s2 = nn.Linear(45, 30)
    #def attention_net(self, lstm_output):
    #    attn_weight_matrix = self.W_s2(torch.tanh(self.W_s1(lstm_output)))
    #    attn_weight_matrix = attn_weight_matrix.permute(0, 2, 1)
    #    attn_weight_matrix = F.softmax(attn_weight_matrix, dim=2)
    #    return attn_weight_matrix
    def forward(self, g):
        feature_ptn = g[1]
        feature_lgd = g[0]
        feature_ptn = self.ptn_apm(feature_ptn)
        feature_lgd = self.lgd_apm(feature_lgd)
        
        sequence = torch.cat((feature_lgd, feature_ptn), dim=0).view(1, -1, int(self.rn/4)*self.pns)
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
        out,att = self.attention(output, mask=mask, return_attention=True)
        #attn_weight_matrix = self.attention_net(output)
        #out = torch.bmm(attn_weight_matrix, output)
        out = F.relu(self.fc_in(out.view(-1, out.size()[1]*out.size()[2])))
        out = torch.sigmoid(self.fc_out(out))
        return out,att




class CNN_model_new(nn.Module):
    def __init__(self, rn, pns, dropout, pockn):
        super(CNN_model, self).__init__()
        self.rn = rn
        self.pns = pns 
        self.pockn = pockn
        self.ptn_apm = nn.Sequential(
            nn.Conv1d(in_channels = rn, out_channels = int(rn/2), kernel_size = 2, stride = 1, padding = 0, dilation =1, groups= 1, bias = True, padding_mode= 'zeros'),
            nn.Conv1d(in_channels = int(rn/2), out_channels = int(rn/4), kernel_size = 2, stride = 1, padding = 0, dilation =1, groups= 1, bias = True, padding_mode= 'zeros'),
            nn.BatchNorm1d(int(rn/4), eps=1e-05, momentum=0.1, affine=True, track_running_stats=True),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.AdaptiveMaxPool1d(pns)
        )
        self.lgd_apm = nn.Sequential(
            nn.Conv1d(in_channels = rn, out_channels = int(rn/2), kernel_size = 2, stride = 1, padding = 0, dilation =1, groups= 1, bias = True, padding_mode= 'zeros'),
            nn.Conv1d(in_channels = int(rn/2), out_channels = int(rn/4), kernel_size = 2, stride = 1, padding = 0, dilation =1, groups= 1, bias = True, padding_mode= 'zeros'),
            nn.BatchNorm1d(int(rn/4), eps=1e-05, momentum=0.1, affine=True, track_running_stats=True),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.AdaptiveMaxPool1d(pns)
        )
        self.bilstm = nn.LSTM(int(rn/4)*pns, int(rn/4)*pns, num_layers=1, bidirectional=True, dropout=dropout)
        self.fc_in = nn.Linear(int(rn/4)*pns*pockn*2, int(rn/4)*pns*pockn) #1922
        self.fc_out = nn.Linear(int(rn/4)*pns*pockn, 1)
        self.attention = MultiHeadAttention(int(rn/4)*pns*2, int(rn/4)*pns*2, 2)
    #    self.W_s1 = nn.Linear(60, 45) #62
    #    self.W_s2 = nn.Linear(45, 30)
    #def attention_net(self, lstm_output):
    #    attn_weight_matrix = self.W_s2(torch.tanh(self.W_s1(lstm_output)))
    #    attn_weight_matrix = attn_weight_matrix.permute(0, 2, 1)
    #    attn_weight_matrix = F.softmax(attn_weight_matrix, dim=2)
    #    return attn_weight_matrix
    def forward(self, g):
        feature_ptn = g[1]
        feature_lgd = g[0]
        feature_ptn = self.ptn_apm(feature_ptn)
        feature_lgd = self.lgd_apm(feature_lgd)
        
        feature_ptn

        h_0 = Variable(torch.zeros(2, 1, int(self.rn/4)*self.pns).cuda())
        c_0 = Variable(torch.zeros(2, 1, int(self.rn/4)*self.pns).cuda())
        
        output, _ = self.bilstm(sequence, (h_0, c_0))


        sequence = torch.cat((feature_lgd, feature_ptn), dim=0).view(1, -1, int(self.rn/4)*self.pns)
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
        out,att = self.attention(output, mask=mask, return_attention=True)
        #attn_weight_matrix = self.attention_net(output)
        #out = torch.bmm(attn_weight_matrix, output)
        out = F.relu(self.fc_in(out.view(-1, out.size()[1]*out.size()[2])))
        out = torch.sigmoid(self.fc_out(out))
        return out,att


dtip = '/spstorage/USERS/gina/Project/AP/DB/InterpretableDTIP'
# allsets = pd.read_csv(dtip+'/data/data_dataset_3D.csv', index_col = 0)
# allsets = pd.read_csv(dtip+'/data/data_dataset_3D_target_split.csv', index_col = 0)

# allsets = pd.read_csv(dtip+'/human/human_dataset_3D_cid_split.csv', index_col = 0)

pns = 7
pkn = ['0','11','6']
ltsn = ['0','1','2']
lrt = 1e-3
# CNN bindingdb basic 


# 
pocket_r = [[[[],[]],[[],[]]] , [[[],[]],[[],[]]]]
pocket_r1 = [[[],[]],[[],[]]]
for i in range(len(X_train)):
	pid = trdf.iloc[i,3]
	actidx = list(apock.loc[pid,:][550]).index(1)
	if actidx!=29:
		mr = model([X_train[i][0].to(device), X_train[i][1].to(device)])
		pidx = X_train[i][1].shape[0]
		for j in range(2):
			att = [ab.item() for ab in mr[1][0][j][X_train[i][1].shape[0]]]
			if att.index(max(att)) == actidx:
				pocket_r[j][y_train[i]][round(mr[0].item())].append(i)


len(pocket_r[0][0][0])
len(pocket_r[0][0][1])
len(pocket_r[0][1][0])
len(pocket_r[0][1][1])

len(pocket_r[1][0][0])
len(pocket_r[1][0][1])
len(pocket_r[1][1][0])
len(pocket_r[1][1][1])

### resolution 등의 cut - off 를 설정 해서 잘라 봤을 때 performance가 오르지? 
rridx = pd.read_csv(dtip+'/pdbval/rrridx.csv', index_col=0)
rridx = pd.read_csv(dtip+'/pdbval/eeidx.csv', index_col=0)
rridx = pd.read_csv(dtip+'/pdbval/gbb.csv', index_col=0)

cmp2cid = pd.read_csv(DIR_AP+'/DB/InterpretableDTIP/pdbval/cmp2cid2.txt', header=None, sep='\t')
cmp2cid.columns = ['cmp','cid']
/data/AP1_no_both_random.pkl
lts = pd.read_pickle(dtip+'/data/AP1_both_random.pkl')

lts = pd.read_pickle(dtip+'/pdbval/AP1_renew222.pkl')
lts = pd.read_pickle(dtip+'/data/AP_z.pkl')
lts2 = pd.read_pickle(dtip+'/data/AP1_z_all.pkl')
lts = pd.concat([lts, lts2[~lts2.index.isin(lts.index)]])
lts2 = pd.read_pickle(dtip+'/data/AP1_covid_z_both.pkl')
lts = pd.concat([lts, lts2[~lts2.index.isin(lts.index)]])
lts = lts[lts.index.isin(list(set(itertools.chain(*[ff.PUBCHEM_CID for ff in fadx2]))))]
lts2 = pd.read_pickle(dtip+'/covid/AP1_no_both.pkl')
lts = pd.concat([lts, lts2[~lts2.index.isin(lts.index)]])

lts2 = pd.read_pickle(dtip+'/data/AP1_no_both_random2.pkl')


lts = pd.concat([lts, pd.read_pickle(dtip+'/data/AP1_no_both_random.pkl')])
lts = pd.concat([lts, pd.read_pickle(dtip+'/data/AP1_no_both_random2.pkl')])
lts = lts[lts]

lts2 = pd.read_pickle(dtip+'/data/AP1_covid_z_both.pkl')
lts = pd.concat([lts, lts2[~lts2.index.isin(lts.index)]])


lts = pd.read_pickle(dtip+'/data/AP_no_both.pkl')

lts = pd.read_pickle(dtip+'/pdbval/AP1_renew_pdb.pkl')
lts = pd.read_pickle(dtip+'/pdbval/AP1_renew.pkl')
lts = pd.concat([lts, pd.read_pickle(dtip+'/pdbval/AP1_renew2.pkl')])
lts = pd.read_pickle(dtip+'/pdbval/AP1_renew_model.pkl')
lts = pd.read_pickle(dtip+'/pdbval/AP1_renew_ideal.pkl')
rridx = pd.merge(rridx , cmp2cid)
pock = pd.read_pickle('{}/pdbval/pdb_sub_pockets_{}_20.pkl'.format(dtip,pkn[idx]))
ppock = pd.read_pickle('{}/pdbval/sub_pockets_11_20_renew2_29pocks_kcluster.pkl'.format(dtip))

pock = pd.read_pickle('{}/pdbval/sub_pockets_11_20_29pocks_kcluster_all.pkl'.format(dtip))
pock = pd.read_pickle('{}/pdbval/sub_pockets_11_20_29pocks_kcluster_betalactam_seprate.pkl'.format(dtip))
pock = pd.read_pickle(f'{dtip}/covid/sub_pockets_11_20_29pocks_kcluster_bigger.pkl')
pock = pd.read_pickle(f'{dtip}/covid/sub_pockets_11_20_renew2_29pocks_kcluster_both.pkl')
pock = pd.read_pickle(f'{dtip}/pdbval/sub_pockets_11_20_renew2_29pocks_kcluster_both.pkl')
pock = pd.read_pickle(f'{dtip}/pdbval/sub_pockets_11_20_renew2_29pocks_kcluster_no_both.pkl')



lts = pd.read_pickle(chemm+'/AP1_no_both.pkl')


pock = pd.read_pickle('{}/sub_pockets_11_20_renew2_29pocks_kcluster_no_both.pkl'.format(dtip))
# pock['PDB'] = [cd.upper()[0:4] for cd in list(pock.index)]




pock = pock[~pock.index.isin(['3g33D','6cnkE','3e00A'])]
pock.index= [cd[0:4].upper() for cd in list(pock.index)]
pock.columns = list(range(550))
pock2 = pd.read_pickle(f'{chemm}/sub_pockets_11_20_renew2_29pocks_kcluster_no_both.pkl')
pock2=pock2[~pock2.index.isin(pock.index)]

pock = pd.concat([pock,pock2])


pock = pd.read_pickle('{}/pdbval/sub_pockets_11_20_29pocks_kcluster.pkl'.format(dtip))

dropout=0.3
rn=55
pockn = 30
model = CNN_model(rn, pns, dropout, pockn)
model.to(device)
MODEL_NAME = 'ap1_lr0.001_dr0.3_30pock_7pns_False_final'
# MODEL_NAME = 'ap1_lr0.001_dr0.3_30pock_7pns_False_rp0'
model.load_state_dict(torch.load("{}/data/{}.pth".format(dtip,MODEL_NAME)))




dropout=0.1
rn=66
pockn = 30
model = CNN_model(rn, pns, dropout, pockn)
model.to(device)
MODEL_NAME = 'ap1_lr0.001_dr0.3_30pock_7pns_False_final'
MODEL_NAME = 'lumped_ori_ap1_lr0.0001_dr0.1_30pock_7pns_False_512batch_sample_rp1'
MODEL_NAME = 'lumped_ori_ap1_lr0.0001_dr0.1_30pock_7pns_False_512batch_sample_rp4'

model.load_state_dict(torch.load("{}/data/{}.pth".format(dtip,MODEL_NAME)))

### covid
pock = pd.read_pickle(f'{dtip}/covid/sub_pockets_11_20_29pocks_kcluster.pkl')
pock = pd.read_pickle(f'{dtip}/covid/sub_pockets_11_20_side_29pocks_kcluster.pkl')
pock = pd.read_pickle(f'{dtip}/covid/sub_pockets_11_20_renew2_29pocks_kcluster.pkl')
pock = pd.read_pickle(f'{dtip}/covid/sub_pockets_11_20_renew2_29pocks_kcluster_both.pkl')
pock = pd.read_pickle(f'{dtip}/covid/sub_pockets_11_20_renew2_29pocks_kcluster_only.pkl')

lts = pd.read_pickle(dtip+'/covid/AP_jcim_no_both.pkl')

lts = pd.read_pickle(dtip+'/data/AP1_renew2.pkl')
lts2 = pd.read_pickle(dtip+'/data/AP1_renew2_screen.pkl')
lts3 = pd.read_pickle(dtip+'/data/AP1_renew2_highes.pkl')
lts = pd.concat([lts,lts2,lts3])
pock = pd.read_pickle(f'{dtip}/covid/sub_pockets_11_20_renew2_29pocks_kcluster_bigger.pkl')
pock = pd.read_pickle(f'{dtip}/covid/sub_pockets_11_20_renew2_29pocks_kcluster_lgd.pkl')
pock = pd.read_pickle(f'{dtip}/covid/sub_pockets_11_10_29pocks_kcluster.pkl')

pock = pd.read_pickle(f'{dtip}/covid/sub_pockets_11_20_29pocks_kcluster_bigger.pkl')
pock2 = pd.read_pickle(f'{dtip}/covid/sub_pockets_11_20_29pocks_kcluster.pkl')

pdbs = list(set(pock.index))

lts3 = pd.read_pickle(dtip+'/covid/AP1.pkl')
lts2 = pd.read_pickle(dtip+'/covid/AP1_2.pkl')

lts = pd.concat([lts,lts2[~lts2.index.isin(lts.index)], lts3[~lts3.index.isin(lts.index)]])


apm = np.array(lts) 
apm = torch.tensor(apm, dtype = torch.float).view(len(lts),1,rn,10)

pdbs = ['7aqj', '7a1u','7ax6', '7awu','6ynq', '7aha', '7baj', '6yb7','7ap6', '6y2e']
pdbs = ['7abu','7amj','7aol','7aqe','7aku','7ay7', '7aga', '6yvf','7ans','7ar5']
pdbs = ['7axm', '7axo', '7awr','7avd', '7aww','7aws','7adw', '7arf','7af0', '8ecg']
pdbs = []
pdbs = ['','']

pdbs = ['7aqj', '7a1u','7ax6', '7awu','6ynq', '7aha', '7baj', '6yb7','7ap6', '6y2e','7abu','7amj','7aol','7aqe','7aku','7ay7', '7aga', '6yvf','7ans','7ar5','7axm', '7axo', '7awr','7avd', '7aww','7aws','7adw', '7arf','7af0']


done_pdb = []
# done_pdb = done_pdb[0:10]
yydf = pd.DataFrame()
attdf = pd.DataFrame()
attdf_2 = pd.DataFrame()
attdf01 = pd.DataFrame()
attdf_201 = pd.DataFrame()
for pds in pdbs:
    done_pdb.append(pds)
    print(pds)
    pcks = pock.loc[pds,:]
    pcks = pcks.iloc[:,0:550]
    if (len(pcks)>(pockn-2))&(len(pcks)!=rn*10):
        pcks = pcks[0:(pockn-2)]
    # ppock = torch.tensor(np.array(pcks).reshape([-1,rn,10]), dtype = torch.float)
    ppock = torch.tensor(np.array(pcks, dtype = float).reshape([-1,rn,10]), dtype = torch.float)
    X = [[apm[i], ppock] for i in range(len(lts))]
    yy = []
    atts = []
    atts2 = []
    atts01 = []
    atts02 =[]
    for idx, X_test in enumerate(X):
        x = [0, 0]
        x[0] = X_test[0].to(device)
        x[1] = X_test[1].to(device)
        res = model(x)
        att = res[1]
        pnn = X_test[1].shape[0]
        # wt = [ww.item() for ww in att[0][0][0][0:X_test[1].shape[0]].cpu()]
        atts01.append(np.argmax(np.array([att[0][0][i][0:pnn][0].cpu().item() for i in range(pnn)]))-1)
        atts02.append(np.argmax(np.array([att[0][1][i][0:pnn][0].cpu().item() for i in range(pnn)]))-1)
        atts.append(np.argmax(np.array([ww.item() for ww in att[0][0][0][0:pnn].cpu()]))-1)
        # wt = [ww.item() for ww in att[0][1][0][0:pnn].cpu()]
        atts2.append(np.argmax(np.array([ww.item() for ww in att[0][1][0][0:pnn].cpu()]))-1)
        yy.append(res[0].item())
        if len(yy) % 100000 == 0:
            print(len(yy))
    yydf2 = pd.DataFrame(yy).T   
    attdf2 = pd.DataFrame(atts).T  
    attdf22 = pd.DataFrame(atts2).T   
    attdf201 = pd.DataFrame(atts01).T  
    attdf2201 = pd.DataFrame(atts02).T   
    yydf2.columns = list(lts.index)
    attdf2.columns = list(lts.index)   
    attdf22.columns = list(lts.index)   
    attdf201.columns = list(lts.index)   
    attdf2201.columns = list(lts.index)  
    yydf = pd.concat([yydf, yydf2])
    attdf = pd.concat([attdf, attdf2])
    attdf_2 = pd.concat([attdf_2, attdf22])
    attdf01 = pd.concat([attdf01, attdf201])
    attdf_201 = pd.concat([attdf_201, attdf2201])
    yydf.index = done_pdb
    attdf.index = done_pdb
    attdf_2.index = done_pdb
    attdf01.index = done_pdb
    attdf_201.index = done_pdb
    yydf.columns = list(lts.index)
    attdf.columns = list(lts.index)
    attdf_2.columns = list(lts.index)
    attdf01.columns = list(lts.index)
    attdf_201.columns = list(lts.index)
    yydf.to_pickle(f'{dtip}/covid/yydf_rp_wt_z2.pkl')
    attdf.to_pickle(f'{dtip}/covid/attdf_rp_wt_z2.pkl')
    attdf_2.to_pickle(f'{dtip}/covid/attdf2_rp_wt_z2.pkl')  
    attdf01.to_pickle(f'{dtip}/covid/attdf01_rp_wt_z2.pkl')
    attdf_201.to_pickle(f'{dtip}/covid/attdf201_rp_wt_z2.pkl')   


cds = list(set(lts.index)-set(clts.index))
rcds = random.sample(cds,100000)

lts = pd.concat([lts[lts.index.isin(rcds)],clts])
apm = np.array(lts) 
apm = torch.tensor(apm, dtype = torch.float).view(len(lts),1,rn,10)


apm = np.array(ldf.loc[list(adf.pdb),:]) 
apm = torch.tensor(apm, dtype = torch.float).view(len(ldf.loc[list(adf.pdb),:]),1,rn,10)


rapm = np.array(rldf.loc[list(adf.pdb),:]) 
rapm = torch.tensor(rapm, dtype = torch.float).view(len(rldf.loc[list(adf.pdb),:]),1,rn,10)


zapm = np.array(zdf.loc[list(adf.iloc[:,2]),:]) 
zapm = torch.tensor(zapm, dtype = torch.float).view(len(zdf.loc[list(adf.iloc[:,2]),:]),1,rn,10)

ppocks = []
for pds in list(adf.pdb):
    print(pds)
    pcks = pock.loc[pds,:]
    pcks = pcks.iloc[:,0:550]
    if (len(pcks)>(pockn-2))&(len(pcks)!=rn*10):
        pcks = pcks[0:(pockn-2)]
    # ppock = torch.tensor(np.array(pcks).reshape([-1,rn,10]), dtype = torch.float)
    ppock = torch.tensor(np.array(pcks, dtype = float).reshape([-1,rn,10]), dtype = torch.float)
    ppocks.append(ppock)


rep_pres = []
for ix in range(1000):
    X = [[zapm[i], ppocks[i]] for i in range(len(zapm))]
    pres = []
    for idx, X_test in enumerate(X):
        x = [0, 0]
        x[0] = X_test[0].to(device)
        x[1] = X_test[1].to(device)
        res = model(x)
        pres.append(res[0].cpu().item())
    rep_pres.append(pres)



rep_apres = []
for ix in range(1000):
    X = [[apm[i], ppocks[i]] for i in range(len(zapm))]
    pres = []
    for idx, X_test in enumerate(X):
        x = [0, 0]
        x[0] = X_test[0].to(device)
        x[1] = X_test[1].to(device)
        res = model(x)
        pres.append(res[0].cpu().item())
    rep_apres.append(pres)



rep_rapres = []
for ix in range(1000):
    X = [[rapm[i], ppocks[i]] for i in range(len(zapm))]
    pres = []
    for idx, X_test in enumerate(X):
        x = [0, 0]
        x[0] = X_test[0].to(device)
        x[1] = X_test[1].to(device)
        res = model(x)
        pres.append(res[0].cpu().item())
    rep_rapres.append(pres)


rep_rdf = pd.DataFrame(rep_rapres)
rep_adf = pd.DataFrame(rep_apres)
rep_zdf = pd.DataFrame(rep_pres)
rep_rdf.columns = list(adf.pdb)
rep_adf.columns = list(adf.pdb)
rep_zdf.columns = list(adf.pdb)
rep_rdft = rep_rdf.T
rep_adft = rep_adf.T
rep_zdft = rep_zdf.T

rep_rdft = rep_rdft[rep_rdft.index.isin(iidf.pdb)]
rep_adft = rep_adft[rep_adft.index.isin(iidf.pdb)]
rep_zdft = rep_zdft[rep_zdft.index.isin(iidf.pdb)]



dres = pd.read_csv(DIR_AP+'/DB/gen_dir/result1',sep='\t',header=None)
dres.columns = ['pdb','affinity']
dres['pdb'] = [a.split('_')[0] for a in list(dres.pdb)]
dres['affinity'] = [-float(a) for a in list(dres.affinity)]
rrocdf = rocdf[rocdf.subset>300]
ddres = dres[dres.pdb.isin(rrocdf.pdb)]
rrocdf = rrocdf[rrocdf.pdb.isin(ddres.pdb)]
daff = list(ddres.loc[list(rrocdf.pdb),:].affinity)
rroc = list(rrocdf.roc)


scipy.stats.pearsonr(daff,rroc)




apmdf3 = pd.DataFrame(press)
rrs = [0.6194885361552027, 0.888888888888889, 0.24489795918367352, 0.7957142857142858, 0.6962900316126546, 0.9708960104643558, 0.6785714285714286, 0.3780487804878049, 0.3909090909090909, 0.5413770130923417, 0.6028368794326241, 0.8712089447938505, 0.725, 0.74057792626458, 0.4418103448275862, 0.40162271805273836, 0.8184984520123839, 0.665340712556209]
[apmdf3.iloc[i,i] for i in range(18)]
scipy.stats.pearsonr(rrs,[apmdf3.iloc[i,i] for i in range(18)])


done_pdb = []
# done_pdb = done_pdb[0:10]
yys =[]

w1s = []
bw1s = []
w2s = []
bw2s = []
for pds in pdbs:
    done_pdb.append(pds)
    print(pds)
    pcks = pock.loc[pds,:]
    pcks = pcks.iloc[:,0:550]
    if (len(pcks)>(pockn-2))&(len(pcks)!=rn*10):
        pcks = pcks[0:(pockn-2)]
    # ppock = torch.tensor(np.array(pcks).reshape([-1,rn,10]), dtype = torch.float)
    ppock = torch.tensor(np.array(pcks, dtype = float).reshape([-1,rn,10]), dtype = torch.float)
    X = [[apm[i], ppock] for i in range(len(lts2))]
    yy = []
    w_mh1 = [[],[],[]]
    w_mh2 = [[],[],[]]
    bw_mh1 = [[],[],[]]
    bw_mh2 =[[],[],[]]
    for idx, X_test in enumerate(X):
        x = [0, 0]
        x[0] = X_test[0].to(device)
        x[1] = X_test[1].to(device)
        res = model(x)
        att = res[1]
        pnn = X_test[1].shape[0]
        # wt = [ww.item() for ww in att[0][0][0][0:X_test[1].shape[0]].cpu()]
        wws = att[0][0][pnn][0:pnn].cpu().detach().numpy()
        rrr=ss.rankdata(-wws)
        for ridx in range(1,4):
            if ridx in rrr:
                bw_mh1[ridx-1].append(np.where(rrr==ridx)[0][0])
                w_mh1[ridx-1].append(wws[np.where(rrr==ridx)[0][0]])
            else:
                bw_mh1[ridx-1].append(100)
                w_mh1[ridx-1].append(100)
        wws = att[0][1][pnn][0:pnn].cpu().detach().numpy()
        rrr=ss.rankdata(-wws)
        for ridx in range(1,4):
            if ridx in rrr:
                bw_mh2[ridx-1].append(np.where(rrr==ridx)[0][0])
                w_mh2[ridx-1].append(wws[np.where(rrr==ridx)[0][0]])
            else:
                bw_mh2[ridx-1].append(100)
                w_mh2[ridx-1].append(100)
        yy.append(res[0].item())
        if len(yy) % 100000 == 0:
            print(len(yy))
    yys.append(yy)
    w1s.append(w_mh1)
    bw1s.append(bw_mh1)
    w2s.append(w_mh2)
    bw2s.append(bw_mh2)

for pidx, pds in enumerate(pdbs):
    done_pdb.append(pds)
    print(pds)
    pcks = pock.loc[pds,:]
    pcks = pcks.iloc[:,0:550]
    if (len(pcks)>(pockn-2))&(len(pcks)!=rn*10):
        pcks = pcks[0:(pockn-2)]
    # ppock = torch.tensor(np.array(pcks).reshape([-1,rn,10]), dtype = torch.float)
    ppock = torch.tensor(np.array(pcks, dtype = float).reshape([-1,rn,10]), dtype = torch.float)
    X = [[apm[i], ppock] for i in range(len(lts))]
    yy = []
    w_mh1 = [[],[],[]]
    w_mh2 = [[],[],[]]
    bw_mh1 = [[],[],[]]
    bw_mh2 =[[],[],[]]
    for idx, X_test in enumerate(X):
        x = [0, 0]
        x[0] = X_test[0].to(device)
        x[1] = X_test[1].to(device)
        res = model(x)
        att = res[1]
        pnn = X_test[1].shape[0]
        # wt = [ww.item() for ww in att[0][0][0][0:X_test[1].shape[0]].cpu()]
        wws = att[0][0][pnn][0:pnn].cpu().detach().numpy()
        rrr=ss.rankdata(-wws)
        for ridx in range(1,4):
            if ridx in rrr:
                bw_mh1[ridx-1].append(np.where(rrr==ridx)[0][0])
                w_mh1[ridx-1].append(wws[np.where(rrr==ridx)[0][0]])
            else:
                bw_mh1[ridx-1].append(100)
                w_mh1[ridx-1].append(100)
        wws = att[0][1][pnn][0:pnn].cpu().detach().numpy()
        rrr=ss.rankdata(-wws)
        for ridx in range(1,4):
            if ridx in rrr:
                bw_mh2[ridx-1].append(np.where(rrr==ridx)[0][0])
                w_mh2[ridx-1].append(wws[np.where(rrr==ridx)[0][0]])
            else:
                bw_mh2[ridx-1].append(100)
                w_mh2[ridx-1].append(100)
        yy.append(res[0].item())
        if len(yy) % 100000 == 0:
            print(len(yy))
    yys[pidx].extend(yy)
    w1s[pidx][0].extend(w_mh1[0])
    w1s[pidx][1].extend(w_mh1[1])
    w1s[pidx][2].extend(w_mh1[2])
    bw1s[pidx][0].extend(bw_mh1[0])
    bw1s[pidx][1].extend(bw_mh1[1])
    bw1s[pidx][2].extend(bw_mh1[2])
    w2s[pidx][0].extend(w_mh2[0])
    w2s[pidx][1].extend(w_mh2[1])
    w2s[pidx][2].extend(w_mh2[2])
    bw2s[pidx][0].extend(bw_mh2[0])
    bw2s[pidx][1].extend(bw_mh2[1])
    bw2s[pidx][2].extend(bw_mh2[2])






for ppp in range(29):
    bw2s[ppp][0].extend(bw2s[ppp+29][0])
    bw2s[ppp][1].extend(bw2s[ppp+29][1])
    bw2s[ppp][2].extend(bw2s[ppp+29][2])




yydf2 = pd.DataFrame(yys)
yydf2.index = pdbs
yydf2.columns = list(lts2.index)

# yydf2 = pd.concat([yydf2,ryydf2],axis=1)
# yydf2 = pd.concat([yydf2,yyy.T[~yyy.T.index.isin(yydf2.columns)].T],axis=1)

25219879
17337753
8035385
acts = list(id2inhit[id2inhit.inhibit<50].cid)
iacts = list(id2inhit[id2inhit.inhibit>50].cid)

yydf2.T.sort_values('7ay7').loc[:,'7ay7']
yydf2.loc['7ay7',acts],yydf2.loc['7ay7',acts]
print(roc_auc_score([1]*len(acts) + [0]*len(iacts) , list(yydf2.loc['7ay7',acts])+list(yydf2.loc['7ay7',iacts])))
print(roc_auc_score([1]*len(acts) + [0]*len(iacts) , list(yydf2.loc['7ay7',acts])+list(yydf2.loc['7ay7',iacts])))

fprs, tprs, thresholds  = roc_curve([1]*len(acts) + [0]*len(iacts), list(yydf2.loc['7ay7',acts])+list(yydf2.loc['7ay7',iacts]))

plt.plot(fprs, tprs)

pid = '7ay7'
for aidx,adx1 in enumerate(xadx2):
    if xaids[aidx] in ['1409608','1409613','1409614']:
        adx = adx1.loc[:,['PUBCHEM_CID','Activity Comment']].iloc[2:,:]
        adx['act'] = ['Active' if cm=='ACTIVE' else 'Inactive' if cm[0:8]=='INACTIVE' else cm for cm in list(adx['Activity Comment'])]
        adx = adx.dropna()
        adx = adx.loc[:,['PUBCHEM_CID','act']]
    else:
        adx = adx1.loc[:,['PUBCHEM_CID','PUBCHEM_ACTIVITY_OUTCOME']].iloc[2:,:]
        adx = adx.dropna()
    adx.columns = ['cid','act']
    adx = adx[~adx.cid.duplicated()]
    apx = pd.merge(adx,cid2p)
    apx = apx[(~np.isnan(apx.prt))&(~apx.cid.isin(pcds))&(apx.prt.isin(pcds))]
    if len(apx)!=0:
        apx = apx.loc[:,['prt','act']]
        apx.columns = ['cid','act']
        adx = pd.concat([adx,apx])
    adx = adx[adx.cid.isin(pcds)]
    prnks = []
    # aqq = yydf.loc[pid,:]
    aqq = dbns.loc[:,'0']
    arnk = [];acrnk = [];icrnk = []
    scrnk = list(aqq.sort_values(ascending=False).index)
    if 'Active' in list(adx.act):
        act = list(adx[adx.act=='Active'].cid)
        if 'Inactive' in list(adx.act):
            iact = list(adx[adx.act=='Inactive'].cid)
        else:
            iact = list(adx[adx.act!='Active'].cid)
        if xaids[aidx] in ['1409579','1409595','1409599','1409585','1409613']:
            scs = roc_auc_score([1]*len(act) + [0]*len(iact) , list(aqq[act])+list(aqq[iact]))
            aas.append(roc_auc_score([1]*len(act) + [0]*len(iact), list(yydf.loc['7ay7',act])+list(yydf.loc['7ay7',iact])))
            bbs.append(roc_auc_score([1]*len(act) + [0]*len(iact), list(dbns.loc[act,'0'])+list(dbns.loc[iact,'0'])))
            fprs, tprs, thresholds  = roc_curve([1]*len(act) + [0]*len(iact), list(yydf.loc['7ay7',act])+list(yydf.loc['7ay7',iact]))
            fprs, tprs, thresholds  = roc_curve([1]*len(act) + [0]*len(iact), list(dbns.loc[act,'0'])+list(dbns.loc[iact,'0']))
            plt.plot(fprs, tprs, label = f'AID_{xaids[aidx]}')


plt.figure(figsize=(4, 7))

plt.boxplot(aas, positions = [0.89],
                patch_artist = True,
                boxprops=dict(facecolor='palevioletred', color='k'),widths = 0.22)

plt.boxplot(bbs, positions = [1.11],
                patch_artist = True,
                boxprops=dict(facecolor='dodgerblue', color='k'),widths = 0.22)


plt.plot([], c='palevioletred', label='Atom-pair')
plt.plot([], c='dodgerblue', label='DBN(ECFP+CC)')
plt.legend()
# plt.legend(loc = 'upper left')
plt.xticks([1], [1])
plt.ylabel('ROC_AUC_SCORE')
# 
plt.ylim([0.3,1])
plt.title('COVID19 Mpro inhibitor assay')
plt.savefig(f'{dtip}/covid/assay_boxplot.png', dpi = 100)

# plt.show()


fprs, tprs, thresholds  = roc_curve([1]*len(acts) + [0]*len(iacts), list(dbns.loc[acts,'0'])+list(dbns.loc[iacts,'0']))

fprs, tprs, thresholds  = roc_curve([1]*len(acts) + [0]*len(iacts), list(yydf2.loc['7ay7',acts])+list(yydf2.loc['7ay7',iacts]))
            
plt.plot(fprs, tprs, label = f'JCIM_ASSAY')


plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('ROC AUC CURVE for inhibitor screening of COVID19 Mpro assay_DBN')
plt.plot([0, 1], [0, 1], linestyle='--', color='red')#, label='Random Classifier') 
plt.legend()
plt.show()

with open(dtip+'/covid/yys_best1_0_2_5.pkl', 'rb') as f:
    [yys5, w1s5,bw1s5,w2s5,bw2s5,hits,act1,iact1] = pickle.load(f)


yydf2 = pd.DataFrame(yys5)
yydf2.index = ['7ay7']
yydf2.columns = list(lts2.index)
rocs = []
rocs2 = []
fprs, tprs, thresholds  = roc_curve([1]*len(act1) + [0]*len(iact1) , list(yydf2.loc['7ay7',act1])+list(yydf2.loc['7ay7',iact1]))
plt.plot(fprs, tprs, label = f'ff')

for sidx in range(2,4):
    ssubaddp = subaddp[subaddp.source==sidx]
    for act_c in [3]:
        act = hits[sidx-2]
        iact = list(set(ssubaddp[ssubaddp.iloc[:,act_c]>100].cid))
        if (len(act)>1)& (len(iact)>1):
            fprs, tprs, thresholds  = roc_curve([1]*len(act) + [0]*len(iact) , list(yydf2.loc['7ay7',act])+list(yydf2.loc['7ay7',iact]))
            plt.plot(fprs, tprs, label = f'dd')




if (np.mean(itemgetter(*[0,2,5,8,11,14])(rocs))>best_rocc)|(np.mean(itemgetter(*[1,2,5,8,11,14])(rocs))>best_rocc):




yys2 = yys.copy()
w1s2 = w1s.copy()
bw1s2 = bw1s.copy()
w2s2 = w2s.copy()
bw2s2 = bw2s.copy()


yydf = pd.DataFrame(yys)
yydf.index = pdbs
yydf.columns = pre_cds + list(lts.index)
yydf.to_pickle(f'{dtip}/covid/yydf_no_both_random.pkl')


for i in range(3):
    w1df = pd.DataFrame([ws[i] for ws in w1s])
    w1df.index = pdbs
    w1df.columns = pre_cds + list(lts.index)
    w1df.to_pickle(f'{dtip}/covid/w1df_{i}_no_both_random.pkl')



for i in range(3):
    w2df = pd.DataFrame([ws[i] for ws in w2s])
    w2df.index = pdbs
    w2df.columns = pre_cds + list(lts.index)
    w2df.to_pickle(f'{dtip}/covid/w2df_{i}_no_both_random.pkl')


for i in range(3):
    bw1df = pd.DataFrame([ws[i] for ws in bw1s])
    bw1df.index = pdbs
    bw1df.columns = pre_cds + list(lts.index)
    bw1df.to_pickle(f'{dtip}/covid/bw1df_{i}_no_both_random.pkl')


for i in range(3):
    bw2df = pd.DataFrame([ws[i] for ws in bw2s])
    bw2df.index = pdbs
    bw2df.columns = pre_cds + list(lts.index)
    bw2df.to_pickle(f'{dtip}/covid/bw2df_{i}_no_both_random.pkl')


yydf = pd.concat([yy1,yydf])

yydf = yydf.T[~np.isnan(yydf.T.index)].T
yydf = yydf.T[~yydf.T.index.duplicated()].T




attdf.to_pickle(f'{dtip}/covid/attdf_no_both.pkl')
attdf_2.to_pickle(f'{dtip}/covid/attdf2_no_both.pkl')  
attdf01.to_pickle(f'{dtip}/covid/attdf01_no_both.pkl')
attdf_201.to_pickle(f'{dtip}/covid/attdf201_no_both.pkl')   





yydf = pd.concat([pd.read_pickle(f'{dtip}/covid/yydf_rp_wt_10_{i}.pkl') for i in range(1,7)])
attdf = pd.concat([pd.read_pickle(f'{dtip}/covid/attdf_rp_wt_10_{i}.pkl') for i in range(1,7)])
attdf_2 = pd.concat([pd.read_pickle(f'{dtip}/covid/attdf2_rp_wt_10_{i}.pkl') for i in range(1,7)])
attdf01 = pd.concat([pd.read_pickle(f'{dtip}/covid/attdf01_rp_wt_10_{i}.pkl') for i in range(1,7)])
attdf_201 = pd.concat([pd.read_pickle(f'{dtip}/covid/attdf201_rp_wt_10_{i}.pkl') for i in range(1,7)])
attdf = attdf[~attdf.index.duplicated()]
attdf.to_pickle(f'{dtip}/covid/attdf_rp_wt_10.pkl')



done_pdb = []
# done_pdb = done_pdb[0:10]
yydf = pd.DataFrame()
attdf = pd.DataFrame()
attdf_2 = pd.DataFrame()
attdf01 = pd.DataFrame()
attdf_201 = pd.DataFrame()
for pds in pdbs:
    done_pdb.append(pds)
    print(pds)
    pcks = pock.loc[pds,:]
    pcks = pcks.iloc[:,0:660]
    if (len(pcks)>(pockn-1))&(len(pcks)!=rn*10):
        pcks = pcks[0:(pockn-1)]
    # ppock = torch.tensor(np.array(pcks).reshape([-1,rn,10]), dtype = torch.float)
    ppock = torch.tensor(np.array(pcks, dtype = float).reshape([-1,rn,10]), dtype = torch.float)
    X = [[apm[i], ppock] for i in range(len(lts))]
    yy = []
    atts = []
    atts2 = []
    atts01 = []
    atts02 =[]
    for idx, X_test in enumerate(X):
        x = [0, 0]
        x[0] = X_test[0].to(device)
        x[1] = X_test[1].to(device)
        res = model(x)
        att = res[1]
        pnn = X_test[1].shape[0]
        # wt = [ww.item() for ww in att[0][0][0][0:X_test[1].shape[0]].cpu()]
        atts01.append(np.argmax(np.array([att[0][0][i][0:pnn][0].cpu().item() for i in range(pnn)]))-1)
        atts02.append(np.argmax(np.array([att[0][1][i][0:pnn][0].cpu().item() for i in range(pnn)]))-1)
        atts.append(np.argmax(np.array([ww.item() for ww in att[0][0][0][0:pnn].cpu()]))-1)
        # wt = [ww.item() for ww in att[0][1][0][0:pnn].cpu()]
        atts2.append(np.argmax(np.array([ww.item() for ww in att[0][1][0][0:pnn].cpu()]))-1)
        yy.append(res[0].item())
        if len(yy) % 100000 == 0:
            print(len(yy))
    yydf2 = pd.DataFrame(yy).T   
    attdf2 = pd.DataFrame(atts).T  
    attdf22 = pd.DataFrame(atts2).T   
    attdf201 = pd.DataFrame(atts01).T  
    attdf2201 = pd.DataFrame(atts02).T   
    yydf2.columns = list(lts.index)
    attdf2.columns = list(lts.index)   
    attdf22.columns = list(lts.index)   
    attdf201.columns = list(lts.index)   
    attdf2201.columns = list(lts.index)  
    yydf = pd.concat([yydf, yydf2])
    attdf = pd.concat([attdf, attdf2])
    attdf_2 = pd.concat([attdf_2, attdf22])
    attdf01 = pd.concat([attdf01, attdf201])
    attdf_201 = pd.concat([attdf_201, attdf2201])
    yydf.index = done_pdb
    attdf.index = done_pdb
    attdf_2.index = done_pdb
    attdf01.index = done_pdb
    attdf_201.index = done_pdb
    yydf.columns = list(lts.index)
    attdf.columns = list(lts.index)
    attdf_2.columns = list(lts.index)
    attdf01.columns = list(lts.index)
    attdf_201.columns = list(lts.index)
    yydf.to_pickle(f'{dtip}/covid/yydf_renew2_rp_wt_att1.pkl')
    attdf.to_pickle(f'{dtip}/covid/attdf_renew2_rp_wt_att1.pkl')
    attdf_2.to_pickle(f'{dtip}/covid/attdf2_renew2_rp_wt_att1.pkl')  
    attdf01.to_pickle(f'{dtip}/covid/attdf01_renew2_rp_wt_att1.pkl')
    attdf_201.to_pickle(f'{dtip}/covid/attdf201_renew2_rp_wt_att1.pkl')   




yydf.to_pickle(f'{dtip}/covid/yydf_renew2_rp_4.pkl')
attdf.to_pickle(f'{dtip}/covid/attdf_renew2_rp_4.pkl')
attdf_2.to_pickle(f'{dtip}/covid/attdf2_renew2_rp_4.pkl')       

done_pdb = []
# done_pdb = done_pdb[0:10]
yydf = pd.DataFrame()
attdf = pd.DataFrame()
attdf_2 = pd.DataFrame()
for pds in pdbs:
    done_pdb.append(pds)
    print(pds)
    pcks = pock.loc[pds,0:549]
    if (len(pcks)>(pockn-1))&(len(pcks)!=rn*10):
        pcks = pcks[0:(pockn-1)]
    # ppock = torch.tensor(np.array(pcks).reshape([-1,rn,10]), dtype = torch.float)
    ppock = torch.tensor(np.array(pcks, dtype = float).reshape([-1,rn,10]), dtype = torch.float)
    X = [[apm[i], ppock] for i in range(len(lts))]
    yy = []
    atts = []
    atts2 = []
    for idx, X_test in enumerate(X):
        x = [0, 0]
        x[0] = X_test[0].to(device)
        x[1] = X_test[1].to(device)
        res = model(x)
        att = res[1]
        wt = [ww.item() for ww in att[0][0][X_test[1].shape[0]][0:X_test[1].shape[0]].cpu()]
        atts.append(wt.index(max(wt)))
        wt = [ww.item() for ww in att[0][1][X_test[1].shape[0]][0:X_test[1].shape[0]].cpu()]
        atts2.append(wt.index(max(wt)))
        yy.append(res[0].item())
        if len(yy) % 100000 == 0:
            print(len(yy))
    yydf2 = pd.DataFrame(yy).T   
    attdf2 = pd.DataFrame(atts).T  
    attdf22 = pd.DataFrame(atts2).T          
    yydf2.columns = list(lts.index)
    attdf2.columns = list(lts.index)   
    attdf22.columns = list(lts.index)   
    yydf = pd.concat([yydf, yydf2])
    attdf = pd.concat([attdf, attdf2])
    attdf_2 = pd.concat([attdf_2, attdf22])
    yydf.index = done_pdb
    attdf.index = done_pdb
    attdf_2.index = done_pdb
    yydf.columns = list(lts.index)
    attdf.columns = list(lts.index)
    attdf_2.columns = list(lts.index)
    yydf.to_pickle(f'{dtip}/covid/yydf_rp_wt_10_9.pkl')
    attdf.to_pickle(f'{dtip}/covid/attdf_rp_wt_10_9.pkl')
    attdf_2.to_pickle(f'{dtip}/covid/attdf_2_rp_wt_10_9.pkl')       


pock = pd.read_pickle(f'{dtip}/covid/sub_pockets_11_20_29pocks_kcluster_bigger.pkl')

done_pdb = []
# done_pdb = done_pdb[0:10]
yydf = pd.DataFrame()
attdf = pd.DataFrame()
attdf_2 = pd.DataFrame()
for pds in pdbs:
    done_pdb.append(pds)
    print(pds)
    pcks = pock.loc[pds,0:549]
    if (len(pcks)>(pockn-1))&(len(pcks)!=rn*10):
        pcks = pcks[0:(pockn-1)]
    # ppock = torch.tensor(np.array(pcks).reshape([-1,rn,10]), dtype = torch.float)
    ppock = torch.tensor(np.array(pcks, dtype = float).reshape([-1,rn,10]), dtype = torch.float)
    X = [[apm[i], ppock] for i in range(len(lts))]
    yy = []
    atts = []
    atts2 = []
    for idx, X_test in enumerate(X):
        x = [0, 0]
        x[0] = X_test[0].to(device)
        x[1] = X_test[1].to(device)
        res = model(x)
        att = res[1]
        wt = [ww.item() for ww in att[0][0][X_test[1].shape[0]][0:X_test[1].shape[0]].cpu()]
        atts.append(wt.index(max(wt)))
        wt = [ww.item() for ww in att[0][1][X_test[1].shape[0]][0:X_test[1].shape[0]].cpu()]
        atts2.append(wt.index(max(wt)))
        yy.append(res[0].item())
        if len(yy) % 100000 == 0:
            print(len(yy))
    yydf2 = pd.DataFrame(yy).T   
    attdf2 = pd.DataFrame(atts).T  
    attdf22 = pd.DataFrame(atts2).T          
    yydf2.columns = list(lts.index)
    attdf2.columns = list(lts.index)   
    attdf22.columns = list(lts.index)   
    yydf = pd.concat([yydf, yydf2])
    attdf = pd.concat([attdf, attdf2])
    attdf_2 = pd.concat([attdf_2, attdf22])
    yydf.index = done_pdb
    attdf.index = done_pdb
    attdf_2.index = done_pdb
    yydf.columns = list(lts.index)
    attdf.columns = list(lts.index)
    attdf_2.columns = list(lts.index)
    yydf.to_pickle(f'{dtip}/covid/yydf_2.pkl')
    attdf.to_pickle(f'{dtip}/covid/attdf_2.pkl')
    attdf_2.to_pickle(f'{dtip}/covid/attdf2_2.pkl')       




###




rridx = rridx[rridx.PDB.isin(pock.index)]

rridx.columns = ['PDB', 'resolution', 'release year', 'K', 'Kd/Ki', 'reference','cmp', 'cid']
idx =1
## beta lactam 
lts = pd.read_pickle(dtip+'/pdbval/AP1_beta.pkl')
lts = pd.concat([lts, pd.read_pickle(dtip+'/pdbval/AP1_beta2.pkl')])
lts = pd.concat([lts, pd.read_pickle(dtip+'/pdbval/AP1_beta3.pkl')])
lts = pd.concat([lts, pd.read_pickle('{}/data/AP{}.pkl'.format(dtip,'1_renew'))])
chemm = '/spstorage/USERS/gina/Project/AP/DB/ChEMBL'
lts = pd.concat([lts, pd.read_pickle('{}/AP_{}.pkl'.format(chemm,ltsn[idx]))])
lts = pd.concat([lts, pd.read_pickle('{}/AP_{}_neg.pkl'.format(chemm,ltsn[idx]))])
lts2 = pd.read_pickle('{}/AP_{}_renew2.pkl'.format(chemm,ltsn[idx]))
lts2 = lts2[~lts2.index.isin(lts.index)]
lts = pd.concat([lts, lts2])
lts2 = pd.read_pickle('{}/human/AP1_renew.pkl'.format(dtip))
lts2 = lts2[~lts2.index.isin(lts.index)]
lts = pd.concat([lts, lts2])
lts = lts[~lts.index.duplicated()]
lts = lts[(lts.T.sum()!=0)]
lts = lts.sort_index()

pock = pd.read_pickle('{}/pdbval/sub_pockets_11_20_29pocks_kcluster_betalactam_seprate.pkl'.format(dtip))

dstss = list(pock.iloc[:,551])
dstss.sort()

apm = np.array(lts.loc[rridx.cid,:]) 
apm = torch.tensor(apm, dtype = torch.float).view(len(rridx),1,rn,10)
# pdbs = list(set(pock.index))
ppock = {}
for pds in list(set(rridx.PDB)):
    pcks = pock.loc[pds,0:549]
    ppock[pds] = torch.tensor(np.array(pcks).reshape([-1,rn,10]), dtype = torch.float)


done_pdb = []
# done_pdb = done_pdb[0:10]
yys =[]

w1s = []
bw1s = []
w2s = []
bw2s = []

X = [[apm[i], ppock[pds]] for i,pds in enumerate(rridx.PDB)]
yy = []
w_mh1 = [[],[],[]]
w_mh2 = [[],[],[]]
bw_mh1 = [[],[],[]]
bw_mh2 =[[],[],[]]
for idx, X_test in enumerate(X):
    x = [0, 0]
    x[0] = X_test[0].to(device)
    x[1] = X_test[1].to(device)
    res = model(x)
    att = res[1]
    pnn = X_test[1].shape[0]
    # wt = [ww.item() for ww in att[0][0][0][0:X_test[1].shape[0]].cpu()]
    wws = att[0][0][pnn][0:pnn].cpu().detach().numpy()
    rrr=ss.rankdata(-wws)
    for ridx in range(1,4):
        if ridx in rrr:
            bw_mh1[ridx-1].append(np.where(rrr==ridx)[0][0])
            w_mh1[ridx-1].append(wws[np.where(rrr==ridx)[0][0]])
        else:
            bw_mh1[ridx-1].append(100)
            w_mh1[ridx-1].append(100)
    wws = att[0][1][pnn][0:pnn].cpu().detach().numpy()
    rrr=ss.rankdata(-wws)
    for ridx in range(1,4):
        if ridx in rrr:
            bw_mh2[ridx-1].append(np.where(rrr==ridx)[0][0])
            w_mh2[ridx-1].append(wws[np.where(rrr==ridx)[0][0]])
        else:
            bw_mh2[ridx-1].append(100)
            w_mh2[ridx-1].append(100)
    yy.append(res[0].item())
    if len(yy) % 100000 == 0:
        print(len(yy))



pdbss = list(rridx.PDB)
for btps in ['both']:
    pock = pd.read_pickle(f'{dtip}/pdbval/sub_pockets_11_20_renew2_29pocks_kcluster_{btps}.pkl')
    for ltps in ['','_model','_ideal']:
        best_scss1=1; best_scss2=1;  best_cscss1=1; best_cscss2=1; best_pwtss = 0.1
        if ltps=='':
            lts = pd.read_pickle(f'{dtip}/pdbval/AP1_{btps}{ltps}.pkl')
            apm = np.array(lts.loc[rridx.cid,:]) 
        elif ltps=='_pdb':
            lts = pd.read_pickle(f'{dtip}/pdbval/AP1_renew2_{btps}{ltps}.pkl')
            apm = np.array(lts.loc[rridx.PDB,:]) 
        else:
            lts = pd.read_pickle(f'{dtip}/pdbval/AP1_renew2_{btps}{ltps}.pkl')
            apm = np.array(lts.loc[rridx.cmp,:]) 
        apm = torch.tensor(apm, dtype = torch.float).view(len(rridx),1,rn,10)
        ppock = {}
        for pds in pdbss:
            pcks = pock.loc[pds,0:549]
            ppock[pds] = torch.tensor(np.array(pcks).reshape([-1,rn,10]), dtype = torch.float)
        X = [[apm[i], ppock[pds]] for i,pds in enumerate(rridx.PDB)]
        lpocks = [len(ppock[pds]) for pds in pdbss]
        rpocki = dict(zip(pdbss,[np.argmin(pock.loc[pds,:][550]) if pock.loc[pds,:][550][np.argmin(pock.loc[pds,:][550])]<=20 else '-' for pds in pdbss]))
        for epc in range(100):
            rr1 = [];  ww1 = [];rr2 = [];   ww2 = [] ; yy =[]
            for idx, X_test in enumerate(X):
                x = [0, 0]
                x[0] = X_test[0].to(device)
                x[1] = X_test[1].to(device)
                res = model(x)
                att = res[1]
                pnn = X_test[1].shape[0]
                wws = att[0][0][pnn][0:pnn].cpu().detach().numpy()
                rrr=ss.rankdata(-wws)
                ww1.extend(wws)
                rr1.extend(rrr)
                wws = att[0][1][pnn][0:pnn].cpu().detach().numpy()
                rrr=ss.rankdata(-wws)
                ww2.extend(wws)
                rr2.extend(rrr)
                yy.append(res[0].item())
                if len(yy) % 100000 == 0:
                    print(len(yy))
            rpock = pock.loc[pdbss,:].copy()
            rpock['rank1'] = rr1
            rpock['wht1'] = ww1
            rpock['rank2'] = rr2
            rpock['wht2'] = ww2
            best_pids = [] ; best_pwts = [] ; b1 = [];  b2 = []; cb1 = []; cb2 = []
            for pid in pdbss:
                prpock = rpock.loc[pid,:]
                minidx = np.argmin(prpock[550])
                if rpocki[pid]=='-':
                    aaxx = 0; aaxx2 = 0 ;cb1.append(5) ;cb2.append(5)
                    for rk in range(1,4):
                        if prpock.iloc[minidx,551]==rk:
                            b1.append(rk); aaxx = 1
                            if (rk==1)&(len(prpock)==29):
                                if prpock.iloc[minidx,552]>0.10:
                                    best_pwts.append(prpock.iloc[minidx,552])
                                    best_pids.append(pid)
                        if prpock.iloc[minidx,553]==rk:
                            if (rk==1)&(len(prpock)==29):
                                if prpock.iloc[minidx,554]>0.10:
                                    best_pwts.append(prpock.iloc[minidx,554])
                                    best_pids.append(pid)
                            b2.append(rk); aaxx2 = 1
                    if aaxx==0: b1.append(4)
                    if aaxx2==0:b2.append(4)
                else:
                    aaxx = 0; aaxx2 = 0 ;
                    for rk in range(1,4):
                        if prpock.iloc[minidx,551]==rk:
                            b1.append(rk); cb1.append(rk); aaxx = 1
                            if (rk==1)&(len(prpock)==29):
                                if prpock.iloc[minidx,552]>0.10:
                                    best_pwts.append(prpock.iloc[minidx,552])
                                    best_pids.append(pid)
                        if prpock.iloc[minidx,553]==rk:
                            if (rk==1)&(len(prpock)==29):
                                if prpock.iloc[minidx,554]>0.10:
                                    best_pwts.append(prpock.iloc[minidx,554])
                                    best_pids.append(pid)
                            b2.append(rk);cb2.append(rk); aaxx2 = 1
                    if aaxx==0: b1.append(4);cb1.append(4)
                    if aaxx2==0:b2.append(4);cb2.append(4)
            bb2 = [[b2[i] for i,lp in enumerate(lpocks) if lp==pn] for pn in range(2,30)]
            bef2 = [(bbi.count(1)/len(bbi))/(1/(i+2)) for i,bbi in enumerate(bb2)]
            br21= [(bbi.count(1)/len(bbi)) for i,bbi in enumerate(bb2)]
            br22= [(bbi.count(2)/len(bbi)) for i,bbi in enumerate(bb2)]
            br23= [(bbi.count(3)/len(bbi)) for i,bbi in enumerate(bb2)]
            br2a= [br21[i] +br22[i]+ br23[i] for i,bbi in enumerate(bb2)]
            cbb2 = [[cb2[i] for i,lp in enumerate(lpocks) if (lp==pn)&(cb2[i]!=5)] for pn in range(2,30)]
            cbef2 = [(bbi.count(1)/len(bbi))/(1/(i+2)) for i,bbi in enumerate(cbb2)]
            cbr21= [(bbi.count(1)/len(bbi)) for i,bbi in enumerate(cbb2)]
            cbr22= [(bbi.count(2)/len(bbi)) for i,bbi in enumerate(cbb2)]
            cbr23= [(bbi.count(3)/len(bbi)) for i,bbi in enumerate(cbb2)]
            cbr2a= [cbr21[i] +cbr22[i]+ cbr23[i] for i,bbi in enumerate(cbb2)]
            bb1 = [[b1[i] for i,lp in enumerate(lpocks) if lp==pn] for pn in range(2,30)]
            bef1 = [(bbi.count(1)/len(bbi))/(1/(i+2)) for i,bbi in enumerate(bb1)]
            br11= [(bbi.count(1)/len(bbi)) for i,bbi in enumerate(bb1)]
            br12= [(bbi.count(2)/len(bbi)) for i,bbi in enumerate(bb1)]
            br13= [(bbi.count(3)/len(bbi)) for i,bbi in enumerate(bb1)]
            br1a= [br11[i] +br12[i]+ br13[i] for i,bbi in enumerate(bb1)]
            cbb1 = [[cb1[i] for i,lp in enumerate(lpocks) if (lp==pn)&(cb1[i]!=5)] for pn in range(2,30)]
            cbef1 = [(bbi.count(1)/len(bbi))/(1/(i+2)) for i,bbi in enumerate(cbb1)]
            cbr11= [(bbi.count(1)/len(bbi)) for i,bbi in enumerate(cbb1)]
            cbr12= [(bbi.count(2)/len(bbi)) for i,bbi in enumerate(cbb1)]
            cbr13= [(bbi.count(3)/len(bbi)) for i,bbi in enumerate(cbb1)]
            cbr1a= [cbr21[i] +cbr22[i]+ cbr23[i] for i,bbi in enumerate(cbb1)]
            scss1 = np.median(bef1)
            scss2 = np.median(bef2)
            cscss1 = np.median(cbef1)
            cscss2 = np.median(cbef2)
            if scss1>best_scss1:
                best_scss1 = scss1
                with open(f'{dtip}/pdbval/best_ef1_{btps}{ltps}.pkl','wb') as f:
                    pickle.dump([rpock, yy, best_pids, best_pwts],f)
            if scss2>best_scss2:
                best_scss2 = scss2
                with open(f'{dtip}/pdbval/best_ef2_{btps}{ltps}.pkl','wb') as f:
                    pickle.dump([rpock, yy, best_pids, best_pwts],f)
            if cscss1>best_cscss1:
                best_cscss1 = cscss1
                with open(f'{dtip}/pdbval/best_cef1_{btps}{ltps}.pkl','wb') as f:
                    pickle.dump([rpock, yy, best_pids, best_pwts],f)
            if cscss2>best_cscss2:
                best_cscss2 = cscss2
                with open(f'{dtip}/pdbval/best_cef2_{btps}{ltps}.pkl','wb') as f:
                    pickle.dump([rpock, yy, best_pids, best_pwts],f)
            if len(best_pwts)>=1:
                if max(best_pwts)>best_pwtss:
                    best_pwtss = max(best_pwts)
                    with open(f'{dtip}/pdbval/best_wt_{btps}{ltps}.pkl','wb') as f:
                        pickle.dump([rpock, yy, best_pids, best_pwts],f)



with open(f'{dtip}/pdbval/best_ef2_{btps}{ltps}2.pkl','rb') as f:
    [rpock, yy, best_pids, best_pwts] = pickle.load(f)



plt.bar(range(2,30), br21, color='b')
plt.bar(range(2,30), br22, color='r', bottom=br21) 
plt.bar(range(2,30), br23, color='c', bottom=[br21[i] +br22[i] for i in range(len(br21))]) 
plt.plot(range(2,30),[1 if ii<=3 else 3/ii for ii in range(2,30)], color = 'c')
plt.plot(range(2,30),[2/ii for ii in range(2,30)],color='r')
plt.plot(range(2,30),[1/ii for ii in range(2,30)],color ='b')


plt.show()

plt.bar(range(2,30), cbef2, color='b')
plt.axhline(1, c='red')
plt.show()



wws = [0.435018, 0.022736, 0.016995, 0.014793, 0.011895, 0.004475]
colors = ['red','hotpink','cornflowerblue','lightseagreen','darkviolet','forestgreen']
for ii,ww in enumerate(wws):
    plt.bar(ii, ww, color=colors[ii])
    plt.text(ii-0.33, ww+0.005, round(ww,3))


plt.title('predicted weight on potential binding site')
plt.xlabel('Binding Site')
plt.ylabel('weight')
plt.xticks(range(0,6),range(1,7))
plt.show()

with open(f'{dtip}/target_ef_no_both2.pkl', 'rb') as f:
    resnb = pickle.load(f)

resnb = dict(zip(cdss,ress))
# with open(f'{dtip}/target_ef_ecfp4_random.pkl', 'wb') as f:
#     pickle.dump(resnb, f)


# with open(f'{dtip}/target_ef_no_both_short.pkl', 'wb') as f:
#     pickle.dump(resnbb, f)


with open(f'{dtip}/target_ef_no_both_short.pkl', 'rb') as f:
    resnb = pickle.load(f)


with open(f'{dtip}/target_ef_no_both_550K_long.pkl', 'rb') as f:
    resns = pickle.load(f)


resnb.update(resns)

# with open(f'{dtip}/target_ef_no_both.pkl', 'rb') as f:
#     resnb3 = pickle.load(f)



with open(f'{dtip}/target_ef_ecfp4.pkl', 'rb') as f:
    resef = pickle.load(f)



# with open(f'{dtip}/target_ef_graph.pkl', 'rb') as f:
#     resgr = pickle.load(f)


with open(f'{dtip}/target_ef_graph_short.pkl', 'rb') as f:
    resgrst = pickle.load(f)

# with open(f'{dtip}/target_ef_ECFP6_2048.pkl', 'rb') as f:
#     res6 = pickle.load(f)



# with open(f'{dtip}/target_ef_ECFP4_1024.pkl', 'rb') as f:
#     res4 = pickle.load(f)



resec6 = {}
for rp in range(11):
    with open(f'{dtip}/target_ef_ECFP6_2048_{rp}.pkl', 'rb') as f:
        ress = pickle.load(f)
    resec6.update(ress)



# with open(f'{dtip}/target_ef_graph_short.pkl', 'rb') as f:
#     resgrst = pickle.load(f)

cofs = [100,200,500,1000,1500,2000,5000,10000,50000,100000]
cofs = [100,500,1000,10000]




cds = intersect(list(resnb.keys()),list(resef.keys()))
cds = intersect(cds, list(resec6.keys()))



ap = [[] for i in range(24)]
fp = [[] for i in range(24)]
fp6 = [[] for i in range(24)]
grst = [[] for i in range(24)]
clens = [[] for i in range(24)]
sts = [[] for i in range(24)]
for cd in cds:
    for ii in range(24):
        ij = ii
        if (ii//6==1)|(ii//6==2):ij = ii+6
        elif ii//6==3:ij = ii +24
        if (resnb[cd][ii]!='-')&(resgrst[cd][ii]!='-')&(resec6[cd][ii]!='-')&(resef[cd][ij]!='-'):
            tds = intersect(list(resnb[cd][ii].keys()), list(resgrst[cd][ii].keys()))
            tds = intersect(tds, list(resec6[cd][ii].keys()))
            tds = intersect(tds, list(resef[cd][ij].keys()))
            for td in tds:
                sts[ii].append((cd,td))
                # clens[ii].append(len(aaadf[aaadf.tid==td]))
                ap[ii].append(resnb[cd][ii][td])
                fp6[ii].append(resec6[cd][ii][td])
                fp[ii].append(resef[cd][ij][td])
                grst[ii].append(resgrst[cd][ii][td])


ap[20] = [20 if aa>20 else aa for aa in ap[20]]
fp[20] = [20 if aa>20 else aa for aa in fp[20]]
grst[20] = [20 if aa>20 else aa for aa in grst[20]]
fp6[20] = [20 if aa>20 else aa for aa in fp6[20]]
# fp61 = [20 if aa>20 else aa for aa in fp6[20]]

ap[14] = [20 if aa>20 else aa for aa in ap[14]]
fp[14] = [20 if aa>20 else aa for aa in fp[14]]
grst[14] = [20 if aa>20 else aa for aa in grst[14]]
fp6[14] = [20 if aa>20 else aa for aa in fp6[14]]


ap[8] = [20 if aa>20 else aa for aa in ap[8]]
fp[8] = [20 if aa>20 else aa for aa in fp[8]]
grst[8] = [20 if aa>20 else aa for aa in grst[8]]
fp6[8] = [20 if aa>20 else aa for aa in fp6[8]]

plt.figure(figsize=(8,6))

plt.boxplot(list(itemgetter(*[3,2,1,0])([aa for ii,aa in enumerate(ap) if ii%6==2])), positions = [a+0.3 for a in range(1,5)],
                patch_artist = True,
                vert=False,
                capprops=dict(color='k'),
                whiskerprops=dict(color='k'),
                flierprops=dict(color='k', markeredgecolor='k'),
                boxprops=dict(facecolor='palevioletred', color='k'),
                medianprops=dict(color='k'),
                widths = 0.15)


plt.boxplot(list(itemgetter(*[3,2,1,0])([aa for ii,aa in enumerate(grst) if ii%6==2])), positions = [a+0.1 for a in range(1,5)],
                patch_artist = True,
                vert=False,
                capprops=dict(color='k'),
                whiskerprops=dict(color='k'),
                flierprops=dict(color='k', markeredgecolor='k'),
                boxprops=dict(facecolor='dodgerblue', color='k'),
                medianprops=dict(color='k'),
                widths = 0.15)



plt.boxplot(list(itemgetter(*[3,2,1,0])([aa for ii,aa in enumerate(fp) if ii%6==2])), positions = [a-0.1 for a in range(1,5)],
                patch_artist = True,
                vert=False,
                capprops=dict(color='k'),
                whiskerprops=dict(color='k'),
                flierprops=dict(color='k', markeredgecolor='k'),
                boxprops=dict(facecolor='#47A992', color='k'),
                medianprops=dict(color='k'),
                widths = 0.15)



plt.boxplot(list(itemgetter(*[3,2,1,0])([aa for ii,aa in enumerate(fp6) if ii%6==2])), positions = [a-0.3 for a in range(1,5)],
                patch_artist = True,
                vert=False,
                capprops=dict(color='k'),
                whiskerprops=dict(color='k'),
                flierprops=dict(color='k', markeredgecolor='k'),
                boxprops=dict(facecolor='#7A3E3E', color='k'),
                medianprops=dict(color='k'),
                widths = 0.15)


plt.xscale('log')
# plt.title('Enrichment score ')
plt.xlabel('Enrichment score of similarity')
plt.ylabel('cut-off')

plt.plot([], c='palevioletred', label='Atom-pair')
plt.plot([], c='dodgerblue', label='graph2vec')
plt.plot([], c='#47A992', label='ECFP4')
plt.plot([], c='#7A3E3E', label='ECFP6')
plt.legend(loc = 'lower right')
plt.yticks(range(1,5), [10000,1000,500,100])
# plt.yticks(range(1,11), reversed(cofs))
plt.show()





plt.scatter(ap1,fp1,s=5)
plt.xlabel('AP')
plt.ylabel('FP')
plt.show()

plt.scatter(ap1,clens[2],s=5)
plt.xlabel('AP')
plt.ylabel('tg_lens')
plt.show()


plt.figure(figsize=(20,8))
for i in range(len(cofs)):
    plt.subplot(2,5,i+1)
    ap1 = [math.log(aa,2) for aa in ap[i*6+2]]
    fp1 = [math.log(aa,2) for aa in fp[i*6+2]]
    plt.hist((ap1,fp1), bins = 30, histtype = 'bar')
    plt.title(cofs[i])


plt.show()


ap = [[] for i in range(60)]
fp = [[] for i in range(60)]
clens = [[] for i in range(60)]
sts = [[] for i in range(60)]
for cd in cds:
    for ii in range(60):
        if (resnb[cd][ii]!='-')& (resef[cd][ii]!='-'):
            tds = intersect(list(resnb[cd][ii].keys()),list(resef[cd][ii].keys()))
            for td in tds:
                sts[ii].append((cd,td))
                clens[ii].append(len(aaadf[aaadf.tid==td]))
                ap[ii].append(resnb[cd][ii][td])
                fp[ii].append(resef[cd][ii][td])


ap1 = [20 if aa>20 else aa for aa in ap[38]]
fp1 = [20 if aa>20 else aa for aa in fp[38]]

plt.scatter(ap1,fp1,s=5)
plt.xlabel('AP')
plt.ylabel('FP')
plt.show()

plt.scatter(ap1,clens[2],s=5)
plt.xlabel('AP')
plt.ylabel('tg_lens')
plt.show()


plt.figure(figsize=(20,8))
for i in range(len(cofs)):
    plt.subplot(2,5,i+1)
    ap1 = [math.log(aa,2) for aa in ap[i*6+2]]
    fp1 = [math.log(aa,2) for aa in fp[i*6+2]]
    plt.hist((ap1,fp1), bins = 30, histtype = 'bar')
    plt.title(cofs[i])


plt.show()



plt.figure(figsize=(20,8))
for i in range(len(cofs)):
    plt.subplot(2,5,i+1)
    ap1 = [50 if aa>50 else aa for aa in ap[i*6+2]]
    fp1 = [50 if aa>50 else aa for aa in fp[i*6+2]]
    plt.scatter(ap1,fp1,s=3)
    plt.xlabel('AP')
    plt.ylabel('FP')
    plt.title(cofs[i])

plt.show()



plt.figure(figsize=(20,8))
for i in range(len(cofs)):
    plt.subplot(2,5,i+1)
    ap1 = [math.log(aa,2) for aa in ap[i*6+5]]
    fp1 = [math.log(aa,2) for aa in fp[i*6+5]]
    plt.hist((ap1,fp1), bins = 30, histtype = 'bar')
    plt.title(cofs[i])


plt.show()

ap1 = [math.log(aa,2) for aa in ap[18]]
fp1 = [math.log(aa,2) for aa in fp[18]]
plt.hist((ap1,fp1), histtype = 'bar')

plt.show()


                # patch_artist = True,
                # boxprops=dict(facecolor=colors[act_c], color='k'),widths = 0.2)

plt.boxplot(ap[12], positions = [0.75],
                patch_artist = True,
                boxprops=dict(facecolor='palevioletred', color='k'),widths = 0.25)


plt.boxplot(fp[12], positions = [1.25],
                patch_artist = True,
                boxprops=dict(facecolor='dodgerblue', color='k'),widths = 0.25)


plt.boxplot(ap[14], positions = [0.85],
                patch_artist = True,
                boxprops=dict(facecolor='palevioletred', color='k'),widths = 0.25)


plt.boxplot(fp[14], positions = [1.15],
                patch_artist = True,
                boxprops=dict(facecolor='dodgerblue', color='k'),widths = 0.25)




itertools.chain(*[0,2,3,7])([aa for ii,aa in enumerate(ap) if ii%6==2])



plt.boxplot(list(itemgetter(*[0,2,3,7])([aa for ii,aa in enumerate(ap) if ii%6==2])), positions = [a-0.15 for a in range(1,5)],
                patch_artist = True,
                vert=True,
                capprops=dict(color='palevioletred'),
                whiskerprops=dict(color='palevioletred'),
                flierprops=dict(color='palevioletred', markeredgecolor='palevioletred'),
                boxprops=dict(facecolor='white', color='palevioletred'),
                medianprops=dict(color='palevioletred'),
                widths = 0.25)


plt.boxplot(list(itemgetter(*[0,2,3,7])([aa for ii,aa in enumerate(fp) if ii%6==2])), positions = [a+0.15 for a in range(1,5)],
                patch_artist = True,
                vert=True,
                capprops=dict(color='dodgerblue'),
                whiskerprops=dict(color='dodgerblue'),
                flierprops=dict(color='dodgerblue', markeredgecolor='dodgerblue'),
                boxprops=dict(facecolor='white', color='dodgerblue'),
                medianprops=dict(color='dodgerblue'),
                widths = 0.25)


plt.yscale('log')
# plt.title('Enrichment score ')
plt.ylabel('Enrichment score of similarity')
plt.xlabel('cut-off')

plt.plot([], c='palevioletred', label='Atom-pair')
plt.plot([], c='dodgerblue', label='ECFP4')
plt.legend()
# plt.legend(loc = 'upper left')
plt.xticks(range(1,5), [100,500,1000,10000])
# plt.savefig(f'{dtip}/EF_nb_ef.png', dpi = 100)
plt.show()

plt.figure(figsize=(8,6))

plt.boxplot(list(itemgetter(*reversed(range(10)))([aa for ii,aa in enumerate(ap) if ii%6==1])), positions = [a+0.15 for a in range(1,11)],
                patch_artist = True,
                vert=False,
                capprops=dict(color='palevioletred'),
                whiskerprops=dict(color='palevioletred'),
                flierprops=dict(color='palevioletred', markeredgecolor='palevioletred'),
                boxprops=dict(facecolor='white', color='palevioletred'),
                medianprops=dict(color='palevioletred'),
                widths = 0.25)


plt.boxplot(list(itemgetter(*reversed(range(10)))([aa for ii,aa in enumerate(fp) if ii%6==1])), positions = [a-0.15 for a in range(1,11)],
                patch_artist = True,
                vert=False,
                capprops=dict(color='dodgerblue'),
                whiskerprops=dict(color='dodgerblue'),
                flierprops=dict(color='dodgerblue', markeredgecolor='dodgerblue'),
                boxprops=dict(facecolor='white', color='dodgerblue'),
                medianprops=dict(color='dodgerblue'),
                widths = 0.25)
# palevioletred

plt.figure(figsize=(8,6))

plt.boxplot(list(itemgetter(*[7,3,2,0])([aa for ii,aa in enumerate(ap) if ii%6==2])), positions = [a+0.15 for a in range(1,5)],
                patch_artist = True,
                vert=False,
                capprops=dict(color='k'),
                whiskerprops=dict(color='k'),
                flierprops=dict(color='k', markeredgecolor='k'),
                boxprops=dict(facecolor='palevioletred', color='k'),
                medianprops=dict(color='k'),
                widths = 0.25)


plt.boxplot(list(itemgetter(*[7,3,2,0])([aa for ii,aa in enumerate(fp) if ii%6==2])), positions = [a-0.15 for a in range(1,5)],
                patch_artist = True,
                vert=False,
                capprops=dict(color='k'),
                whiskerprops=dict(color='k'),
                flierprops=dict(color='k', markeredgecolor='k'),
                boxprops=dict(facecolor='dodgerblue', color='k'),
                medianprops=dict(color='k'),
                widths = 0.25)


plt.xscale('log')
# plt.title('Enrichment score ')
plt.xlabel('Enrichment score of similarity')
plt.ylabel('cut-off')

plt.plot([], c='palevioletred', label='Atom-pair')
plt.plot([], c='dodgerblue', label='ECFP4')
plt.legend(loc = 'upper left')
# plt.yticks(range(1,5), [10000,1000,500,100])
plt.yticks(range(1,11), reversed(cofs))
plt.show()
# plt.savefig(f'{dtip}/EF_nb_ef.png', dpi = 100)


plt.scatter(ap[12],fp[12])



                best_scss = scss
                plt.bar(range(1,30), tips_sum_by_day_male, color='b', alpha=alpha)
                plt.bar(range(1,30), tips_sum_by_day_female, color='r', alpha=alpha, bottom=tips_sum_by_day_male) 
                plt.bar(range(1,30), tips_sum_by_day_female, color='r', alpha=alpha, bottom=tips_sum_by_day_male) 
                plt.title('1 / 2 / 3 percentile rank', fontsize=20)
                plt.ylabel('percentile', fontsize=18)
                plt.xlabel('Day', fontsize=18)
                plt.xticks(range(1,30), range(1,30), fontsize=15)
                plt.legend((p1[0], p2[0]), ('Male', 'Female'), fontsize=15)
                plt.savefig(f'{dtip}/pdbval/best_{btps}{ltps}.png')
                plt.bar(range(1,30), tips_sum_by_day_male, color='b', alpha=alpha)
                plt.bar(range(1,30), tips_sum_by_day_female, color='r', alpha=alpha) 
                plt.bar(range(1,30), tips_sum_by_day_female, color='r', alpha=alpha) 
                plt.title('Enrichment score for 1st / 2nd / 3rd', fontsize=20)
                plt.ylabel('EF', fontsize=18)
                plt.xlabel('Day', fontsize=18)
                plt.xticks(range(1,30), range(1,30), fontsize=15)
                plt.legend((p1[0], p2[0]), ('Male', 'Female'), fontsize=15)
                plt.savefig(f'{dtip}/pdbval/best_{btps}{ltps}.png')
                with open(f'{dtip}/pdbval/best_{btps}{ltps}.pkl','wb') as f:
                    pickle.dump([rpock, yy, best_pids, best_pwts],f)
        bb1 = [[b1[i] for i,lp in enumerate(lpocks) if lp==pn] for pn in range(2,30)]
        br1 = [(bbi.count(1)/len(bbi))/(1/(i+2)) for i,bbi in enumerate(bb1)]
        br13 = [(len([ii for ii in bbi if ii<=3])/len(bbi))/(3/(i+2)) for i,bbi in enumerate(bb1)]


bb1 = [[r1[i] for i,lp in enumerate(lpocks) if (lp==pn)&(r1[i]!=5)] for pn in range(2,30)]
br1 = [(bbi.count(1)/len(bbi))/(1/(i+2)) for i,bbi in enumerate(bb1)]
br13 = [(len([ii for ii in bbi if ii<=3])/len(bbi))/(3/(i+2)) for i,bbi in enumerate(bb1)]



b2 = [bw_mh2[0][i]==cpocki[i] for i in range(len(rridx))]

yys.append(yy)
w1s.append(w_mh1)
bw1s.append(bw_mh1)
w2s.append(w_mh2)
bw2s.append(bw_mh2)


dtip = '/spstorage/USERS/gina/Project/AP/DB/InterpretableDTIP'


with open(f'{dtip}/pdbval/best_{btps}{ltps}.pkl','rb') as f:
    aa = pickle.load(f)


rpock, yy, best_pids, best_pwts = aa

pdbss = list(set(rpock.index))
lpocks = [len(rpock.loc[pid,:]) for pid in pdbs]
ppock = [pid for pid in pdbs if len(rpock.loc[pid,:])==29]

aa = rpock[rpock.index.isin(ppock)]

max(aa[aa.rank1==1].wht1)
max(aa[aa.rank2==1].wht2)


pdbs = []
for dd in dstss:
    pds = list(pock[pock.iloc[:,551]==dd].index)[0]
    if pds not in pdbs:
        pdbs.append(pds)



done_pdb = []
# done_pdb = done_pdb[0:10]
yydf = pd.DataFrame()
attdf = pd.DataFrame()
attdf_2 = pd.DataFrame()
for pds in pdbs[200:300]+pdbs[500:600]:
    # if pds not in done_pdb + xx:
    done_pdb.append(pds)
    print(pds)
    pcks = pock.loc[pds,0:549]
    if (len(pcks)>(pockn-1))&(len(pcks)!=rn*10):
        pcks = pcks[0:(pockn-1)]
    ppock = torch.tensor(np.array(pcks).reshape([-1,rn,10]), dtype = torch.float)
    X = [[apm[i], ppock] for i in range(len(lts))]
    yy = []
    atts = []
    atts2 = []
    for idx, X_test in enumerate(X):
        x = [0, 0]
        x[0] = X_test[0].to(device)
        x[1] = X_test[1].to(device)
        res = model(x)
        att = res[1]
        wt = [ww.item() for ww in att[0][0][X_test[1].shape[0]][0:X_test[1].shape[0]].cpu()]
        atts.append(wt.index(max(wt)))
        wt = [ww.item() for ww in att[0][1][X_test[1].shape[0]][0:X_test[1].shape[0]].cpu()]
        atts2.append(wt.index(max(wt)))
        yy.append(res[0].item())
        if len(yy) % 100000 == 0:
            print(len(yy))
    yydf2 = pd.DataFrame(yy).T   
    attdf2 = pd.DataFrame(atts).T  
    attdf22 = pd.DataFrame(atts2).T          
    yydf2.columns = list(lts.index)
    attdf2.columns = list(lts.index)   
    attdf22.columns = list(lts.index)   
    yydf = pd.concat([yydf, yydf2])
    attdf = pd.concat([attdf, attdf2])
    attdf_2 = pd.concat([attdf_2, attdf22])
    yydf.index = done_pdb
    attdf.index = done_pdb
    attdf_2.index = done_pdb
    yydf.columns = list(lts.index)
    attdf.columns = list(lts.index)
    attdf_2.columns = list(lts.index)
    yydf.to_pickle(f'{dtip}/yydf_all3.pkl')
    attdf.to_pickle(f'{dtip}/attdf_all3.pkl')
    attdf_2.to_pickle(f'{dtip}/attdf2_all3.pkl')       



beta_int[beta_int.UniProtKB.isin(df.Entry)]
beta_int[beta_int.cid.isin(lts.index)]



lts = lts[~lts.index.duplicated()]


df = df[~df.Entry.duplicated()]

beta_int = pd.merge(beta_int,df.loc[:,['Entry','PDB']].drop_duplicates(), left_on = 'UniProtKB',right_on = 'Entry')
beta_int=beta_int[beta_int.cid.isin(lts.index)]
beta_int = beta_int[(beta_int.standard_value>=100000)|(beta_int.standard_value<=1000)]
beta_int = beta_int[beta_int.PDB.isin(pock.index)]

# bbb = bb[bb.PDB.isin(pdb_list.PDB)]
apm = np.array(lts.loc[beta_int.cid,:]) 


apm = torch.tensor(apm, dtype = torch.float).view(len(beta_int),1,rn,10)
ppock = {}
for pds in list(set(beta_int.PDB)):
    pcks = pock.loc[pds,0:549]
    if (len(pcks)>(pockn-1))&(len(pcks)!=rn*10):
        pcks = pcks[0:(pockn-1)]
    ppock[pds] = torch.tensor(np.array(pcks).reshape([-1,rn,10]), dtype = torch.float)

# apm = np.array(lts) 

# apm = torch.tensor(apm, dtype = torch.float).view(len(lts),1,rn,10)

yy = []
X = [[apm[i], ppock[pid]] for i, pid in enumerate(list(beta_int.PDB))]
ll = [1 if ab<=1000 else 0 for ab in list(beta_int.standard_value)]
atts = []
for idx, X_test in enumerate(X):
    x = [0, 0]
    x[0] = X_test[0].to(device)
    x[1] = X_test[1].to(device)
    res = model(x)
    att = res[1]
    wt = [ww.item() for ww in att[0][0][X_test[1].shape[0]][0:X_test[1].shape[0]].cpu()]
    atts.append(wt.index(max(wt)))
    yy.append(res[0].item())
    if len(yy)%100000==0:print(len(yy))




beta_int['prob'] = yy
beta_int['att_pock_idx'] = atts

beta_int = beta_int.loc[:,['PDB','UniProtKB','cid','standard_value','prob','att_pock_idx', 'Organism','Protein names','standard_relation','standard_units','standard_flag','standard_type','activity_comment','assay_id','Entrez','etid']]

beta_int.to_pickle(dtip+'/pdbval/beta_int_with_prob2.pkl')

aucs_int10000_2 = roc_auc_score(ll,yy)




## all pdb (gpdb)
# cid / pdb / ideal / model
llss = []
aucss = []
pockets =[]
ellss =[]
lgd_dsts =[]
for sidx in [4]:
    for eidx in [8]:
        if (sidx==5) & (eidx==5):
            continue
        bb= rridx[(rridx['K']>=eidx)|(rridx['K']<=sidx)]
        for ap_idx, ap_type in enumerate([('','cid')]):
            # lts = pd.read_pickle(f'{dtip}/pdbval/AP1_renew{ap_type[0]}.pkl')
            apm = np.array(lts.loc[bb[ap_type[1]],:])        
            apm = torch.tensor(apm, dtype = torch.float).view(len(bb),1,rn,10)
            ppock = []
            for pds in list(bb.PDB):
                pcks = pock.loc[pds,0:549]
                if (len(pcks)>(pockn-1))&(len(pcks)!=rn*10):
                    pcks = pcks[0:(pockn-1)]
                ppock.append(torch.tensor(np.array(pcks).reshape([-1,rn,10]), dtype = torch.float))
            X = [[apm[i], ppock[i]] for i in range(len(apm))]
            kk = list(bb['K'])
            pids = list(bb.PDB)
            for rpp in range(5):
yy= []
ll = []
pocket_sub = [[[],[]],[[],[]]] ; lgd_dst = [[[],[]],[[],[]]]
for idx, X_test in enumerate(X):
    x = [0, 0]
    x[0] = X_test[0].to(device)
    x[1] = X_test[1].to(device)
    res = model(x)
    # att = res[1]
    yy.append(res[0].item())
    if kk[idx]>=eidx:
        rsp = 1
        ll.append(1)
    else:
        ll.append(0)
        rsp = 0
    for i in range(2):
        wt = [ww.item() for ww in att[0][i][X_test[1].shape[0]][0:X_test[1].shape[0]].cpu()]
        if list(pock.index).count(pids[idx])>1:
            lgdst = list(pock.loc[pids[idx],550])
        else:
            lgdst = [pock.loc[pids[idx],550]]
        if wt.index(max(wt)) == lgdst.index(min(lgdst)):
            if pids[idx] not in pocket_sub[rsp][round(res[0].item())]:
                pocket_sub[rsp][round(res[0].item())].append(pids[idx])
        lgd_dst[rsp][round(res[0].item())].append(lgdst[wt.index(max(wt))])
                print((sidx,eidx))
                aucs = roc_auc_score(ll,yy)
                print(aucs)
                pockets.append(pocket_sub)
                aucss.append(aucs)
                llss.append(sum([len(pocket_sub[i][j]) for i in range(2) for j in range(2)]))
                lgd_dsts.append(lgd_dst)
                

## all pdb (gpdb)
# cid / pdb / ideal / model
llss = []
aucss = []
pockets =[]
ellss =[]
lgd_dsts =[]
yys= []
sidx = 4 
eidx = 8
rridx.columns = ['PDB', 'resolution', 'release year', 'K', 'Kd/Ki', 'reference', 'cmp', 'cid']
rridx = rridx[rridx.PDB.isin(pock.index)]
rridx= rridx[(rridx['K']>=eidx)|(rridx['K']<=sidx)]


apm = np.array(lts.loc[rridx.cid,:])        
apm = torch.tensor(apm, dtype = torch.float).view(len(rridx),1,rn,10)
ppock = []
for pds in list(rridx.PDB):
    pcks = pock.loc[pds,0:659]
    if (len(pcks)>(pockn-1))&(len(pcks)!=rn*10):
        pcks = pcks[0:(pockn-1)]
    ppock.append(torch.tensor(np.array(pcks).reshape([-1,rn,10]), dtype = torch.float))



X = [[apm[i], ppock[i]] for i in range(len(apm))]
kk = list(rridx['K'])
ll = [1 if kkk>=8 else 0 for kkk in kk]
pids = list(rridx.PDB)
pock.columns = list(range(660))+['dst']

for rpp in range(100):
    yy= []
    dsts = []
    # pocket_sub = [[[],[]],[[],[]]] ; lgd_dst = [[[],[]],[[],[]]]
    porank = [[],[]]; poatt = [[],[]]; pocklen = []; bestpock_dst = [[],[]]; 
    for idx, X_test in enumerate(X):
        x = [0, 0]
        x[0] = X_test[0].to(device)
        x[1] = X_test[1].to(device)
        res = model(x)
        att = res[1]
        yy.append(res[0].item())
        rsp = ll[idx]
        if list(pock.index).count(pids[idx])>1:
            lgdst = list(pock.loc[pids[idx],:].dst)
        else:
            lgdst = [pock.loc[pids[idx],'dst']]
        dsts.append(min(lgdst))
        pocklen.append(len(lgdst))
        att_idx = lgdst.index(min(lgdst))
        for i in range(2):
            wt = [ww.item() for ww in att[0][i][X_test[1].shape[0]][0:X_test[1].shape[0]].cpu()]    
            att_v = wt[att_idx]
            bestpock_dst[i].append(lgdst[wt.index(max(wt))])
            wt.sort(reverse=True)
            porank[i].append(wt.index(att_v))
            poatt[i].append(att_v)
    aucs = roc_auc_score(ll,yy)
    print(aucs)
    poranks.append(porank)
    aucss.append(aucs)
    poatts.append(poatt)
    yys.append(yy)
    pocklens.append(pocklen)
    dstss.append(dsts)
    bestpock_dsts.append(bestpock_dst)

apm = np.array(lts.loc[rridx.cid,:])        
apm = torch.tensor(apm, dtype = torch.float).view(len(rridx),1,rn,10)
ppock = []
for pds in list(rridx.PDB):
    pcks = pock.loc[pds,0:659]
    if (len(pcks)>(pockn-1))&(len(pcks)!=rn*10):
        pcks = pcks[0:(pockn-1)]
    ppock.append(torch.tensor(np.array(pcks).reshape([-1,rn,10]), dtype = torch.float))



X = [[apm[i], ppock[i]] for i in range(len(apm))]
kk = list(rridx['K'])
ll = [1 if kkk>=8 else 0 for kkk in kk]
pids = list(rridx.PDB)
pock.columns = list(range(660))+['dst']

for rpp in range(100):
    yy= []
    dsts = []
    # pocket_sub = [[[],[]],[[],[]]] ; lgd_dst = [[[],[]],[[],[]]]
    porank = [[],[]]; poatt = [[],[]]; pocklen = []; bestpock_dst = [[],[]]; 
    for idx, X_test in enumerate(X):
        x = [0, 0]
        x[0] = X_test[0].to(device)
        x[1] = X_test[1].to(device)
        res = model(x)
        att = res[1]
        yy.append(res[0].item())
        rsp = ll[idx]
        if list(pock.index).count(pids[idx])>1:
            lgdst = list(pock.loc[pids[idx],:].dst)
        else:
            lgdst = [pock.loc[pids[idx],'dst']]
        dsts.append(min(lgdst))
        pocklen.append(len(lgdst))
        att_idx = lgdst.index(min(lgdst))
        for i in range(2):
            wt = [ww.item() for ww in att[0][i][X_test[1].shape[0]][0:X_test[1].shape[0]].cpu()]    
            att_v = wt[att_idx]
            bestpock_dst[i].append(lgdst[wt.index(max(wt))])
            wt.sort(reverse=True)
            porank[i].append(wt.index(att_v))
            poatt[i].append(att_v)
    aucs = roc_auc_score(ll,yy)
    print(aucs)
    poranks.append(porank)
    aucss.append(aucs)
    poatts.append(poatt)
    yys.append(yy)
    pocklens.append(pocklen)
    dstss.append(dsts)
    bestpock_dsts.append(bestpock_dst)


dstss = []
poatts = []
pocklens = []
yys = []
aucss = []
poranks = []
bestpock_dsts = []
for rpp in range(100):
    yy= []
    dsts = []
    # pocket_sub = [[[],[]],[[],[]]] ; lgd_dst = [[[],[]],[[],[]]]
    porank = [[],[]]; poatt = [[],[]]; pocklen = []; bestpock_dst = [[],[]]; 
    for idx, X_test in enumerate(X):
        x = [0, 0]
        x[0] = X_test[0].to(device)
        x[1] = X_test[1].to(device)
        res = model(x)
        att = res[1]
        yy.append(res[0].item())
        rsp = ll[idx]
        if list(pock.index).count(pids[idx])>1:
            lgdst = list(pock.loc[pids[idx],:].dst)
        else:
            lgdst = [pock.loc[pids[idx],'dst']]
        dsts.append(min(lgdst))
        pocklen.append(len(lgdst))
        att_idx = lgdst.index(min(lgdst))
        for i in range(2):
            wt = [ww.item() for ww in att[0][i][X_test[1].shape[0]][0:X_test[1].shape[0]].cpu()]    
            att_v = wt[att_idx]
            bestpock_dst[i].append(lgdst[wt.index(max(wt))])
            wt.sort(reverse=True)
            porank[i].append(wt.index(att_v))
            poatt[i].append(att_v)
    aucs = roc_auc_score(ll,yy)
    print(aucs)
    poranks.append(porank)
    aucss.append(aucs)
    poatts.append(poatt)
    yys.append(yy)
    pocklens.append(pocklen)
    dstss.append(dsts)
    bestpock_dsts.append(bestpock_dst)

    
import matplotlib.pyplot as plt


from scipy.stats.stats import pearsonr   

pearsonr(yys[6],poranks[6][0])
pearsonr(yys[6],poranks[6][1])



acts = [idx for idx,pdd in enumerate(list(rridx.PDB)) if pdd in list(bb.PDB)]

acts = [idx for idx,pdd in enumerate(list(rridx.PDB)) if (pdd in list(bb.PDB)) & (dstss[6][idx] <=10)]

acts = [idx for idx,lb in enumerate(ll) if lb == 1]
decs = [idx for idx,lb in enumerate(ll) if lb == 0]

acts = [idx for idx,lb in enumerate(yys[6]) if lb >0.5]

acts = [idx for idx,lb in enumerate(yys[6]) if (lb >=0.9) & (ll[idx]==1)]


decs = [idx for idx,lb in enumerate(yys[6]) if lb <=0.5]
acts = [idx for idx,lb in enumerate(dstss[60]) if lb <=10]
acts = [idx for idx,lb in enumerate(yys[6]) if (lb >=0.5) & (ll[idx]==1) & (dstss[6][idx] <=10)]
np.mean(list(itemgetter(*acts)(poranks[60][0])))
roc_auc_score(list(itemgetter(*acts)(ll)),list(itemgetter(*acts)(yys[6])))
# list(itemgetter(*acts)(yys[6]))
pearsonr(dstss[6],poatts[6][1])

pearsonr(list(itemgetter(*acts)(yys[6])),list(itemgetter(*acts)(poatts[6][0])))
pearsonr(list(itemgetter(*acts)(yys[6])),list(itemgetter(*acts)(poatts[6][1])))


pearsonr(list(itemgetter(*acts)(yys[60])),list(itemgetter(*acts)(poatts[60][0])))
pearsonr(list(itemgetter(*acts)(yys[60])),list(itemgetter(*acts)(poatts[60][1])))


pearsonr(list(itemgetter(*acts)(yys[60])),list(itemgetter(*acts)(poranks[60][0])))
pearsonr(list(itemgetter(*acts)(yys[60])),list(itemgetter(*acts)(poranks[60][1])))



pearsonr(list(itemgetter(*acts)(yys[6])),list(itemgetter(*acts)(poranks[6][0])))
pearsonr(list(itemgetter(*acts)(yys[6])),list(itemgetter(*acts)(poranks[6][1])))


plt.scatter(yys[6],poranks[6][0])
plt.show()

acts = [idx for idx,lb in enumerate(yys[10]) if(dstss[10][idx] ==bestpock_dsts[0][0][idx])]
acts = [idx for idx,lb in enumerate(yys[10]) if (bestpock_dsts[0][0][idx] <=10)]
acts = [idx for idx,lb in enumerate(yys[10]) if (dstss[10][idx] <=10)]

acts = [idx for idx,lb in enumerate(yys[10]) if (lb >=0.5) & (ll[idx]==1) & (dstss[10][idx] <=10)]
np.mean(list(itemgetter(*acts)(poranks[10][0])))

roc_auc_score(list(itemgetter(*acts)(ll)),list(itemgetter(*acts)(yys[10])))
pearsonr(list(itemgetter(*acts)(ll)),list(itemgetter(*acts)(yys[10])))

acts = [idx for idx,lb in enumerate(yys[6]) if (dstss[6][idx] <=16.556128286175618)]
roc_auc_score(list(itemgetter(*acts)(ll)),list(itemgetter(*acts)(yys[6])))

for ii in range(6,20):
    acts = [idx for idx,lb in enumerate(yys[10]) if (bestpock_dsts[0][0][idx] <=ii)]
    if len(set(list(itemgetter(*acts)(ll))))!=1:
        print(ii)
        print(roc_auc_score(list(itemgetter(*acts)(ll)),list(itemgetter(*acts)(yys[10]))))


ii =27
acts = [idx for idx,lb in enumerate(yys[60]) if (bestpock_dsts[60][0][idx] <=ii)|(bestpock_dsts[60][1][idx] <=ii)]
roc_auc_score(list(itemgetter(*acts)(ll)),list(itemgetter(*acts)(yys[60])))
len(acts)
acts = [idx for idx,lb in enumerate(yys[60]) if (bestpock_dsts[60][0][idx] >ii)&(bestpock_dsts[60][1][idx] >ii)]
len(acts)
roc_auc_score(list(itemgetter(*acts)(ll)),list(itemgetter(*acts)(yys[60])))

for ii in range(6,20):
    acts = [idx for idx,lb in enumerate(yys[60]) if (dstss[60][idx] <=ii)]
    if len(set(list(itemgetter(*acts)(ll))))!=1:
        print(ii)
        len(acts)
        print(roc_auc_score(list(itemgetter(*acts)(ll)),list(itemgetter(*acts)(yys[60]))))

aaa=[[],[]]
for ii in range(27,51):
    acts = [idx for idx,lb in enumerate(yys[60]) if (bestpock_dsts[60][0][idx] <=ii)|(bestpock_dsts[60][1][idx] <=ii)]
    if len(acts)>1:
        if len(set(list(itemgetter(*acts)(ll))))!=1:
            print(ii)
            len(acts)
            aaa[0].append(roc_auc_score(list(itemgetter(*acts)(ll)),list(itemgetter(*acts)(yys[60]))))
            print(roc_auc_score(list(itemgetter(*acts)(ll)),list(itemgetter(*acts)(yys[60]))))
            acts = [idx for idx,lb in enumerate(yys[60]) if (bestpock_dsts[60][0][idx] >ii)&(bestpock_dsts[60][1][idx] >ii)]
            print(roc_auc_score(list(itemgetter(*acts)(ll)),list(itemgetter(*acts)(yys[60]))))
            aaa[1].append(roc_auc_score(list(itemgetter(*acts)(ll)),list(itemgetter(*acts)(yys[60]))))
            


iiss = []
for iii in range(200):
    iis = []
    for ii in range(6,50):
        acts = [idx for idx,lb in enumerate(yys[iii]) if (bestpock_dsts[iii][0][idx] <=ii)]
        if len(acts)>1:
            if len(set(list(itemgetter(*acts)(ll))))!=1:
                iis.append(roc_auc_score(list(itemgetter(*acts)(ll)),list(itemgetter(*acts)(yys[iii]))))
    iiss.append(iis)
  
            # [iii, ii, roc_auc_score(list(itemgetter(*acts)(ll)),list(itemgetter(*acts)(yys[iii])))]
biiss = []
for iii in range(200):
    iis = []
    for ii in range(6,50):
        acts = [idx for idx,lb in enumerate(yys[iii]) if (bestpock_dsts[iii][1][idx] <=ii)]
        if (len(acts)>620)&(len(acts)<680):
            iis.append(roc_auc_score(list(itemgetter(*acts)(ll)),list(itemgetter(*acts)(yys[iii]))))
            acts = [idx for idx,lb in enumerate(yys[iii]) if (bestpock_dsts[iii][1][idx] >ii)]
            iis.append(roc_auc_score(list(itemgetter(*acts)(ll)),list(itemgetter(*acts)(yys[iii]))))
    biiss.append(iis)




for iii in range(200):
    if aucss[iii]>0.76:
        print(iii)
        print(biiss[iii])



iiss = []
for iii in range(200):
    iis = []
    for ii in range(6,50):
        acts = [idx for idx,lb in enumerate(yys[iii]) if (bestpock_dsts[iii][0][idx] <=ii)]
        if len(acts)>1:
            if len(set(list(itemgetter(*acts)(ll))))!=1:
                iis.append(roc_auc_score(list(itemgetter(*acts)(ll)),list(itemgetter(*acts)(yys[iii]))))
    iiss.append(iis)
            # [iii, ii, roc_auc_score(list(itemgetter(*acts)(ll)),list(itemgetter(*acts)(yys[iii])))]

for iii in range(200):
    if aucss[iii]>0.76:
        print(iiss[iii])



for ii in range(6,20):
    acts = [idx for idx,lb in enumerate(yys[6]) if (dstss[6][idx] <=ii)]
    if len(set(list(itemgetter(*acts)(ll))))!=1:
        print(ii)
        print(roc_auc_score(list(itemgetter(*acts)(ll)),list(itemgetter(*acts)(yys[6]))))



for ii in range(3,50):
    acts = [idx for idx,pdd in enumerate(list(rridx.PDB)) if (pdd in list(bb.PDB)) & (dstss[9][idx] <=ii)]
    if len(set(list(itemgetter(*acts)(ll))))!=1:
        print(ii)
        print(roc_auc_score(list(itemgetter(*acts)(ll)),list(itemgetter(*acts)(yys[9]))))



aaa = [[0.7812744145277157,0.7741299489144317, 0.7777281043234245, 0.7737228860239913, 0.7701044406987917, 0.7672438330170778, 0.7698017310789049, 0.7678645954092278, 0.7718622080866978, 0.770740263559595, 0.7703775074604895, 0.7677291869128604, 0.7708304757154436, 0.7728787359329231, 0.7725529123768695, 0.7714699074074074, 0.7721779687513445, 0.7725186925870541, 0.7700116659936286, 0.7705769401035902, 0.7698605315176807, 0.768818036003615, 0.7672614995623845, 0.767568658041262], [0.7485528489952562,0.7563578207749988, 0.746348374069792, 0.7495555698344363, 0.7525554588951717, 0.757995143947277, 0.7524668593640984, 0.7544347826086956, 0.7354583727124323, 0.7389277389277389, 0.7361845972957084, 0.7494835680751174, 0.7143589743589743, 0.6866096866096866, 0.6726190476190476, 0.6789473684210525, 0.6535812672176309, 0.6390168970814132, 0.672008547008547, 0.6497326203208557, 0.6466973886328724, 0.6638655462184875, 0.7333333333333334, 0.7306666666666667]]

aaa = [itemgetter(*[0,5,10,15,20])(aaa[0]),itemgetter(*[0,5,10,15,20])(aaa[1])]
ticks = [ib for ib in range(len(aaa[0]))]
poss = []
for ibb in [-0.15,0.15]:
    bb= [a+ibb for a in list(range(len(aaa[0])))]
    poss.append(bb)



cols= ['green','purple','dodgerblue', 'lightpink','orange','grey','red']
flierprops=dict(markersize=1)

label = list(range(27,51,5))
index = np.arange(len(label))

for ibb in range(len(poss)):
    plt.bar(poss[ibb], aaa[ibb], width = 0.3)
    for idx in range(len(poss[ibb])):
        plt.text([abb-0.165 for abb in poss[ibb]][idx], [ab + 0.01 for ab in aaa[ibb]][idx], [round(ab,2) for ab in aaa[ibb]][idx])
    plt.title('PDBbind', fontsize=10)
    plt.xlabel('distance cut-off between exact ligand and best attention pocket from model', fontsize=8)
    plt.ylabel('AUROC', fontsize=8)
    plt.xticks(index, label, fontsize=7)


plt.legend(['short dst','far dst'])
plt.show()



# cid / pdb / ideal / model
llss3 = [[],[],[],[]]
aucss3 = [[],[],[],[]]
pockets3 =[[],[],[],[]]
ellss3 =[]
lgd_dsts3 =[[],[],[],[]]
for sidx in [4,5]:
    for eidx in range(5,10):
        if (sidx==5) & (eidx==5):
            continue
        bb= rridx[(rridx['K']>=eidx)|(rridx['K']<=sidx)]
        for ap_idx, ap_type in enumerate([('','cid'),('_pdb','PDB'),('_model','cmp'),('_ideal','cmp')]):
            lts = pd.read_pickle(f'{dtip}/pdbval/AP1_renew{ap_type[0]}.pkl')
            apm = np.array(lts.loc[bb[ap_type[1]],:])        
            apm = torch.tensor(apm, dtype = torch.float).view(len(bb),1,rn,10)
            ppock = []
            for pds in list(bb.PDB):
                pcks = pock.loc[pds,0:549]
                if (len(pcks)>(pockn-1))&(len(pcks)!=rn*10):
                    pcks = pcks[0:(pockn-1)]
                ppock.append(torch.tensor(np.array(pcks).reshape([-1,rn,10]), dtype = torch.float))
            X = [[apm[i], ppock[i]] for i in range(len(apm))]
            kk = list(bb['K'])
            pids = list(bb.PDB)
            for rpp in range(5):
                yy= []
                ll = []
                pocket_sub = [[[],[]],[[],[]]] ; lgd_dst = [[[],[]],[[],[]]]
                for idx, X_test in enumerate(X):
                    x = [0, 0]
                    x[0] = X_test[0].to(device)
                    x[1] = X_test[1].to(device)
                    res = model(x)
                    att = res[1]
                    yy.append(res[0].item())
                    if kk[idx]>=eidx:
                        rsp = 1
                        ll.append(1)
                    else:
                        ll.append(0)
                        rsp = 0
                    for i in range(2):
                        wt = [ww.item() for ww in att[0][i][X_test[1].shape[0]][0:X_test[1].shape[0]].cpu()]
                        lgdst = list(pock.loc[pids[idx],550])
                        if wt.index(max(wt)) == lgdst.index(min(lgdst)):
                            if pids[idx] not in pocket_sub[rsp][round(res[0].item())]:
                                pocket_sub[rsp][round(res[0].item())].append(pids[idx])
                        lgd_dst[rsp][round(res[0].item())].append(lgdst[wt.index(max(wt))])
                print((sidx,eidx))
                auc = roc_auc_score(ll,yy)
                print(auc)
                pockets3[ap_idx].append(pocket_sub)
                aucss3[ap_idx].append(auc)
                llss3[ap_idx].append(sum([len(pocket_sub[i][j]) for i in range(2) for j in range(2)]))
                lgd_dsts3[ap_idx].append(lgd_dst)
                





with open('{}/pdbval/apm_perf_except_with_dst.pkl'.format(dtip), 'wb') \
        as f_dbn:
    pickle.dump([aucss3, llss3, pockets3,lgd_dsts3], f_dbn)






pocket_iii = [[[],[]],[[],[]],[[],[]]] 

kk = list(rr['-logKd/Ki'])
pids = list(rr.PDB)

lls_i = []
aucs_i = []
pocket_iii =[]
for ii in range(1000):
    yy= []
    ll = []
	pocket_sub = [[[],[]],[[],[]],[[],[]]] 
	for idx, X_test in enumerate(X):
		x = [0, 0]
		x[0] = X_test[0].to(device)
		x[1] = X_test[1].to(device)
		res = model(x)
		att = res[1]
		yy.append(res[0].item())
		if kk[idx]>=9:
			rsp = 2
			ll.append(1)
		else:
			ll.append(0)
			if kk[idx]<=4:
				rsp = 0
			else:
				rsp = 1
		for i in range(2):
			wt = [ww.item() for ww in att[0][i][X_test[1].shape[0]][0:X_test[1].shape[0]].cpu()]
			lgdst = list(pock.loc[pids[idx],550])
			if wt.index(max(wt)) == lgdst.index(min(lgdst)):
				if pids[idx] not in pocket_sub[rsp][round(res[0].item())]:
					pocket_sub[rsp][round(res[0].item())].append(pids[idx])
	auc = roc_auc_score(ll,yy)
	pocket_iii.append(pocket_sub)
	aucs_i.append(auc)
	lls_i.append(sum([len(pocket_sub[i][j]) for i in range(3) for j in range(2)]))





[len(pocket_rrr[i][j]) for i in range(3) for j in range(2)]
sum([len(pocket_rrr[i][j]) for i in range(3) for j in range(2)])

[len(set(pocket_rr[0][i][j]+pocket_rr[1][i][j])) for i in range(3) for j in range(2)]


list(itertools.chain(*list(itertools.chain(*list(itertools.chain(*pocket_rr))))))
[np.mean([min(pock.loc[ppp,550]) for ppp in pp[i][j]]) for pp in pocket_rr for i in range(3) for j in range(2)]
[[min(pock.loc[ppp,550]) for ppp in pp[i][j]] for pp in pocket_rr for i in range(3) for j in range(2)]

[[len(pock.loc[ppp,550]) for ppp in pp[i][j]] for pp in pocket_rr for i in range(3) for j in range(2)]


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
    # roce1 = get_roce(y_pred, y_label, 0.5)
    # roce2 = get_roce(y_pred, y_label, 1)
    # roce3 = get_roce(y_pred, y_label, 2)
    # roce4 = get_roce(y_pred, y_label, 5)
    auroc = str(roc_auc_score(y_label, y_pred))
    print("AUROC: " + auroc, end=" ")



cmp2cid = pd.read_csv(dtip+'/pdbval/cmp2cid2.txt' , sep = '\t', header = None)

cmp2cid = cmp2cid.dropna()

rridx = rridx[rridx.cmp.isin(cmp2cid.iloc[:,0])]
rridx['resolution'] = [float(aa) for aa in list(rridx.resolution)]
rridx[rridx.resolution<2]



## APM human 


yy = []
for idx, X_test in enumerate(X):
    x = [0, 0]
    x[0] = X_test[0].to(device)
    x[1] = X_test[1].to(device)
    res = model(x)
    att = res[1]
    yy.append(res[0].item())



auc = roc_auc_score(y,yy)
print(auc)


### AttentionSite

rridx = pd.read_csv(dtip+'/pdbval/subdf.csv', index_col = 0)

rridx = pd.read_csv(dtip+'/pdbval/rrridx.csv', index_col = 0)

rridx = pd.read_csv(dtip+'/pdbval/eeidx.csv', index_col = 0)



with open(f"/spstorage/USERS/gina/Project/AP/DB/AttentionSite/dataset/pdbbind_pdb2g_all.pkl", 'rb') as fp:
    p_graphs = pickle.load(fp)

model = DTITAG2()

model.load_state_dict(torch.load("/spstorage/USERS/gina/Project/AP/DB/AttentionSite/model_pkl/bindingdb_real.pkl"))
model.to(device)

pp = pd.read_csv(dtip+'/data/pdb_ptn_loc_refined.csv', index_col =0)
pp = pp[pp.iloc[:,0]=='original']

lls_att3 = []
aucs_att3 = []
pocket_att3 =[]
ellss3 =[]
lgd_dsts2 = []
for sidx in [4,5]:
    for eidx in range(5,10):
        if (sidx==5) & (eidx==5):
            continue
        bb= rridx[(rridx['K']>=eidx)|(rridx['K']<=sidx)]
        ll = list(bb['K'])
        negative_train =[]
        ffls = []
        pdbs = []
        for i,pps in enumerate(list(bb.PDB)):
            # print(f"{i}")
            try:
                constructed_graphs = p_graphs[pps]
                llm = Chem.SDMolSupplier('{}/DB/refined-set/{}/{}_ligand.sdf'.format(DIR_AP,pps,pps), removeHs=True)
                smile = Chem.rdmolfiles.MolToSmiles(llm[0])
                if llm is None:
                    llm = Chem.rdmolfiles.MolFromMol2File(DIR_AP+'/DB/refined-set/'+pps+'/'+pps+'_ligand.mol2', removeHs=True)
                    smile = Chem.rdmolfiles.MolToSmiles(llm)
                g = smiles_to_bigraph(smile, node_featurizer=node_featurizer)
                g = dgl.add_self_loop(g)
                negative_train.append(((constructed_graphs, g), ll[i]))
                pdbs.append(pps)
            except Exception as e:
                # print(e)
                ffls.append([i,pps])
                continue
        fflss = []
        for i,pps in ffls:
            # print(f"{i}")
            try:
                constructed_graphs = p_graphs[pps]
                llm = Chem.rdmolfiles.MolFromMol2File(DIR_AP+'/DB/refined-set/'+pps+'/'+pps+'_ligand.mol2', removeHs=True)
                smile = Chem.rdmolfiles.MolToSmiles(llm)
                g = smiles_to_bigraph(smile, node_featurizer=node_featurizer)
                g = dgl.add_self_loop(g)
                negative_train.append(((constructed_graphs, g), ll[i]))
                pdbs.append(pps)
            except Exception as e:
                # print(e)
                fflss.append(pps)
                continue
        X = [i[0] for i in negative_train]
        y = [float(i[1]) for i in negative_train]
        yy= []
        ll = []
        pocket_sub = [[[],[]],[[],[]]] 
        lgd_dst = [[[],[]],[[],[]]] 
        print((sidx,eidx))
        for idx, X_test in enumerate(X):
            subdf = pp[pp.iloc[:,1]==pdbs[idx]]
            if len(subdf)>139:
                ellss3.append(idx)
                continue
            x = [0, 0]
            x[0] = X_test[0].to(device)
            x[1] = X_test[1].to(device)
            res = model(x)
            att = res[1]
            yy.append(res[0].item())
            if y[idx]>=eidx:
                rsp = 1
                ll.append(1)
            else:
                ll.append(0)
                rsp = 0
            for i in range(2):
                wt = [ww.item() for ww in att[0][i][len(subdf)][0:len(subdf)].cpu()]
                lgdst = list(subdf.iloc[:,2])
                lgd_dst[rsp][round(res[0].item())].append(lgdst[wt.index(max(wt))])
                if wt.index(max(wt)) == lgdst.index(min(lgdst)):
                    if pdbs[idx] not in pocket_sub[rsp][round(res[0].item())]:
                        pocket_sub[rsp][round(res[0].item())].append(pdbs[idx])
        auc = roc_auc_score(ll,yy)
        print(auc)
        pocket_att3.append(pocket_sub)
        aucs_att3.append(auc)
        lls_att3.append(sum([len(pocket_sub[i][j]) for i in range(2) for j in range(2)]))
        lgd_dsts2.append(lgd_dst)





with open('{}/pdbval/attn_perf_except_with_dst.pkl'.format(dtip), 'wb') \
        as f_dbn:
    pickle.dump([aucs_att3, lls_att3, pocket_att3, ellss3,lgd_dsts2], f_dbn)





[len(pocket_sub[i][j]) for i in range(3) for j in range(2)]
sum([len(pocket_sub[i][j]) for i in range(3) for j in range(2)])

[len(set(pocket_rr[0][i][j]+pocket_rr[1][i][j])) for i in range(3) for j in range(2)]



pdbs = list(bb.PDB)
lls_att = []
aucs_att = []
pocket_att =[]

yy= []
ll = []
pocket_sub = [[[],[]],[[],[]],[[],[]]] 

for idx, X_test in enumerate(X):
	subdf = pp[pp.iloc[:,1]==pdbs[idx]]
	if len(subdf)>139:
		continue
	x = [0, 0]
	x[0] = X_test[0].to(device)
	x[1] = X_test[1].to(device)
	res = model(x)
	att = res[1]
	yy.append(res[0].item())
	if y[idx]>=6:
		rsp = 1
		ll.append(1)
	else:
		ll.append(0)
		rsp = 0
	for i in range(2):
		wt = [ww.item() for ww in att[0][i][len(subdf)][0:len(subdf)].cpu()]
		lgdst = list(subdf.iloc[:,2])
		if wt.index(max(wt)) == lgdst.index(min(lgdst)):
			if pdbs[idx] not in pocket_sub[rsp][round(res[0].item())]:
				pocket_sub[rsp][round(res[0].item())].append(pdbs[idx])
auc = roc_auc_score(ll,yy)
pocket_att.append(pocket_sub)
aucs_att.append(auc)
lls_att.append(sum([len(pocket_sub[i][j]) for i in range(3) for j in range(2)]))



with open('{}/pdbval/attn_perf.pkl'.format(dtip), 'wb') \
        as f_dbn:
    pickle.dump([aucs_att, lls_att, pocket_att, ellss], f_dbn)



with open('{}/pdbval/drugvqa_perf.pkl'.format(dtip), 'wb') \
        as f_dbn:
    pickle.dump([testAucs, other_perfs], f_dbn)


## drugVQA


for sidx in [4,5]:
	for eidx in range(5,10):
		if (sidx==5) & (eidx==5):
			continue
		bb = eeidx[(eeidx['K']>=eidx)|(eeidx['K']<=sidx)]
		bb['act'] = [1 if aa>=eidx else 0 for aa in list(bb.K)]
		raw_all = [' '.join([str(bb.iloc[i,3]),str(bb.iloc[i,0]),str(bb.iloc[i,7])]) for i in range(len(bb))]
		tt = '\n'.join(raw_all)
		with open(f"/spstorage/USERS/gina/Project/AP/DB/drugVQA/data/pdbbind/dataPre/{sidx}_{eidx}_except", 'w') as fp:
		    fp.write(tt)




with open('{}/pdbval/drugvqa_perf_except.pkl'.format(dtip), 'wb') \
        as f_dbn:
    pickle.dump([testAucs, other_perfs], f_dbn)




with open('{}/ChEMBL/drugvqa_perf_except.pkl'.format(dtip), 'wb') \
        as f_dbn:
    pickle.dump([testAucs, other_perfs], f_dbn)




## DBN


allsets = pd.read_csv(dtip+'/data/data_dataset_3D.csv', index_col = 0)
# cid_info =  pd.read_csv(dtip+'/cid_list.csv')
# allsets['rigid'] = ['rigid' if cid_info[cid_info.cid==cd]['rotbonds'].item()<6 else 'flex' for cd in list(allsets.cid)]


with open('{}/DBN_best_model_bindingDB_final.pkl'.format(base_path), 'rb') \
        as f_dbn:
    dbn = pickle.load(f_dbn)




alts = pd.read_csv('{}/BindingDB_lgd_df.csv'.format(base_path), index_col = 0)
apock = pd.read_csv('{}/BindingDB_ptn_df.csv'.format(base_path), index_col = 0)


pock = pd.read_csv('/spstorage/USERS/gina/Project/AP/DB/DeepDTIs_DBN/data/pdbbind_ptn_df_all.csv', index_col = 0)
lts = pd.read_csv('/spstorage/USERS/gina/Project/AP/DB/DeepDTIs_DBN/data/pdbbind_lgd_df_all.csv', index_col = 0)



rridx = pd.read_csv(dtip+'/pdbval/rrridx.csv', index_col = 0)

aa=pd.read_csv('/spstorage/USERS/gina/Project/AP/DB/pdb_ligand/INDEX_refined_name.2020', sep='\t')
bb = [aaa.split('  ')[0:3] for aaa in list(aa.iloc[:,0])]
bb = pd.DataFrame(bb)
bb.columns = ['PDB','year','uniprot']

rridx = pd.merge(rridx, bb.loc[:,['PDB','uniprot']])
# rridx = rridx[~rridx.uniprot.isin(['P66992', 'P0A5R0', 'B2DJD9', 'Q6G8R1', 'B7GVP4', '------', 'A5MTN0', 'U6NCW5', 'D5CV28', 'P0A4X6', 'D0ZP76'])]



perfss =[]
for sidx in [4,5]:
	for eidx in range(5,10):
		if (sidx==5) & (eidx==5):
			continue
		bb= rridx[(rridx['K']>=eidx)|(rridx['K']<=sidx)]
		bb['label'] = [1 if aa>=eidx else 0 for aa in list(bb['K'])]
		X = pd.concat([pd.concat([lts.loc[bb.PDB, :],alts.loc[allsets.cid, :]]).reset_index(drop=T), pd.concat([pock.loc[bb.uniprot, :],apock.loc[allsets.uniprot_id, :]]).reset_index(drop=T)], axis=1)
		XX = X.to_numpy()
		YY = bb.label.to_numpy()
		X = XX.astype(np.float32)
		Y = YY
		min_max_scaler = MinMaxScaler()  ## min max scaler
		min_max_scaler.fit(X)
		X = min_max_scaler.transform(X)
		test = X[0:len(bb)]
		test_y_pred_proba = dbn.predict_proba(test)[:, 1]
		test_y_pred = dbn.predict(test)
		acc = accuracy_score(YY, test_y_pred)
		auc = roc_auc_score(YY, test_y_pred_proba)
		tpr = list(test_y_pred[YY == 1]).count(1) / float(list(YY).count(1))
		tnr = list(test_y_pred[YY == 0]).count(0) / float(list(YY).count(0))
		perfss.append([acc,auc,tpr,tnr])





rridx = pd.read_csv(dtip+'/pdbval/eeidx.csv', index_col = 0)

aa=pd.read_csv('/spstorage/USERS/gina/Project/AP/DB/pdb_ligand/INDEX_refined_name.2020', sep='\t')
bb = [aaa.split('  ')[0:3] for aaa in list(aa.iloc[:,0])]
bb = pd.DataFrame(bb)
bb.columns = ['PDB','year','uniprot']
rridx = pd.merge(rridx, bb.loc[:,['PDB','uniprot']])


perfss2 =[]
for sidx in [4,5]:
	for eidx in range(5,10):
		if (sidx==5) & (eidx==5):
			continue
		bb= rridx[(rridx['K']>=eidx)|(rridx['K']<=sidx)]
		bb['label'] = [1 if aa>=eidx else 0 for aa in list(bb['K'])]
		X = pd.concat([pd.concat([lts.loc[bb.PDB, :],alts.loc[allsets.cid, :]]).reset_index(drop=T), pd.concat([pock.loc[bb.uniprot, :],apock.loc[allsets.uniprot_id, :]]).reset_index(drop=T)], axis=1)
		XX = X.to_numpy()
		YY = bb.label.to_numpy()
		X = XX.astype(np.float32)
		Y = YY
		min_max_scaler = MinMaxScaler()  ## min max scaler
		min_max_scaler.fit(X)
		X = min_max_scaler.transform(X)
		test = X[0:len(bb)]
		test_y_pred_proba = dbn.predict_proba(test)[:, 1]
		test_y_pred = dbn.predict(test)
		acc = accuracy_score(YY, test_y_pred)
		auc = roc_auc_score(YY, test_y_pred_proba)
		tpr = list(test_y_pred[YY == 1]).count(1) / float(list(YY).count(1))
		tnr = list(test_y_pred[YY == 0]).count(0) / float(list(YY).count(0))
		perfss2.append([acc,auc,tpr,tnr])



with open('{}/pdbval/dbn_perf_except.pkl'.format(dtip), 'wb') \
        as f_dbn:
    pickle.dump(perfss2, f_dbn)



## GNN 


humans = '\n'.join([' '.join([str(aa) for aa in allsets.iloc[i,[5,6,2,3]]]) for i in range(len(allsets))]+[' '.join([str(aa) for aa in allsets[allsets.tag=='test'].iloc[i,[5,6,2,4]]]) for i in range(len(allsets[allsets.tag=='test']))])

with open(base_path + '/original/human_random.txt', 'w') as f:
    f.write(humans)



rridx = pd.read_csv(dtip+'/pdbval/eeidx.csv', index_col = 0)
ptn_df = pd.read_csv(dtip+'/pdbval/uni2seq_all.csv', index_col=0)

aa=pd.read_csv('/spstorage/USERS/gina/Project/AP/DB/pdb_ligand/INDEX_refined_name.2020', sep='\t')
bb = [aaa.split('  ')[0:3] for aaa in list(aa.iloc[:,0])]
bb = pd.DataFrame(bb)
bb.columns = ['PDB','year','uniprot']


allsets = pd.merge(rridx, bb)
allsets = pd.merge(allsets, ptn_df)



for sidx in [4,5]:
	for eidx in range(5,10):
		if (sidx==5) & (eidx==5):
			continue
		cc= allsets[(allsets['K']>=eidx)|(allsets['K']<=sidx)]
		cc['label'] = [1 if aa>=eidx else 0 for aa in list(cc['K'])]
		pdbbind = '\n'.join([' '.join([str(aa) for aa in cc.loc[i,['smiles','sequence','label']]]) for i in list(cc.index)])
		with open(f'{base_path}/original/pdbbind_{sidx}_{eidx}_except.txt', 'w') as f:
		    f.write(pdbbind)




fingerprint_dict = load_pickle('../dataset/human/bindingdb/radius2_ngram3/fingerprint_dict.pickle')
word_dict = load_pickle('../dataset/human/bindingdb/radius2_ngram3/word_dict.pickle')

bkeys = list(fingerprint_dict.keys())
wkeys = list(word_dict.keys())
fii = 20841
wii = 8299
for sidx in [4,5]:
    for eidx in range(5,10):
        if (sidx==5) & (eidx==5):
            continue
        print((sidx,eidx))
DATASET = f'pdbbind_{sidx}_{eidx}'
dir_input = (f'../dataset/{DATASET}/input/radius{radius}_ngram{ngram}/')
ff = load_pickle(dir_input + 'fingerprint_dict.pickle')
ww = load_pickle(dir_input + 'word_dict.pickle')    
for fk in list(ff.keys()):
    if fk not in bkeys:
        fingerprint_dict[fk] = fii
        fii +=1
for wk in list(ww.keys()):
    if wk not in wkeys:
        word_dict[wk] = wii
        wii +=1
bkeys = list(fingerprint_dict.keys())
wkeys = list(word_dict.keys())
        # DATASET = f'pdbbind_{sidx}_{eidx}_except'
        # dir_input = (f'../dataset/{DATASET}/input/radius{radius}_ngram{ngram}/')
        # ff = load_pickle(dir_input + 'fingerprint_dict.pickle')
        # ww = load_pickle(dir_input + 'word_dict.pickle')    
        # for fk in list(ff.keys()):
        #     if fk not in bkeys:
        #         fingerprint_dict[fk] = fii
        #         fii +=1
        # for wk in list(ww.keys()):
        #     if wk not in bkeys:
        #         word_dict[wk] = wii
        #         wii +=1
        # bkeys = list(fingerprint_dict.keys())
        # wkeys = list(word_dict.keys())




DATASET=human_random
DATASET=chembl
DATASET=pdbbind_4_5_except
DATASET=pdbbind_5_9

radius=2
ngram=3

python preprocess_pdb.py $DATASET $radius $ngram


DATASET=pdbbind_5_9_except

radius=2
ngram=3

python preprocess_pdb.py $DATASET $radius $ngram


# python preprocess_human.py $DATASET $radius $ngram
# python preprocess_bindingdb.py $DATASET $radius $ngram
# python preprocess_chembl.py $DATASET $radius $ngram



DATASET=pdbbind_4_5
# DATASET=celegans
# DATASET=yourdata

# radius=1
radius=2
# radius=3

# ngram=2
ngram=3

dim=10
layer_gnn=3
side=5
window=11
layer_cnn=3
layer_output=3
lr=1e-3
lr_decay=0.5
decay_interval=10
weight_decay=1e-6
iteration=100

setting=$DATASET--radius$radius--ngram$ngram--dim$dim--layer_gnn$layer_gnn--window$window--layer_cnn$layer_cnn--layer_output$layer_output--lr$lr--lr_decay$lr_decay--decay_interval$decay_interval--weight_decay$weight_decay--iteration$iteration
python run_pdb.py $DATASET $radius $ngram $dim $layer_gnn $window $layer_cnn $layer_output $lr $lr_decay $decay_interval $weight_decay $iteration $setting
python run_human_split.py $DATASET $radius $ngram $dim $layer_gnn $window $layer_cnn $layer_output $lr $lr_decay $decay_interval $weight_decay $iteration $setting
python run_bindingdb_split.py $DATASET $radius $ngram $dim $layer_gnn $window $layer_cnn $layer_output $lr $lr_decay $decay_interval $weight_decay $iteration $setting





###
allsets = pd.read_csv(dtip+'/data/data_dataset_3D.csv', index_col = 0)


lts = pd.read_pickle('{}/data/AP{}.pkl'.format(dtip,ltsn[idx]))

lts = pd.read_pickle('{}/data/AP1_z_no_both.pkl'.format(dtip))

allsets = allsets.sample(frac = 1, random_state = 7)    
rn = int(len(lts.columns)/10)
norm = False
pockn = 30
pock = pd.read_pickle('{}/sub_pockets_{}_20_{}pocks_kcluster.pkl'.format(dtip,pkn[idx], str(pockn)))
pock = pd.read_pickle('{}/sub_pockets_11_20_renew2_29pocks_kcluster_no_both.pkl'.format(dtip))
# X_train,y_train, std_scaler = setting_dataset(lts, pock, ['train'],pockn,allsets, norm, rn, std_scaler)
# X_dev,y_dev, std_scaler = setting_dataset(lts, pock, ['dev'],pockn,allsets, norm, rn, std_scaler)
X_seen_test,y_seen_test, std_scaler = setting_dataset(lts, pock, ['test','seen'],pockn,allsets, norm, rn, std_scaler)
X_unseen_test,y_unseen_test, std_scaler = setting_dataset(lts, pock, ['test','unseen'],pockn,allsets, norm, rn, std_scaler)
dropout = 0.3

best_unseen = 0.9
# ysts
for epc in range(10000):
    yy_seen = []
    yy_unseen = []
    print(epc)
    for idx, X_test in enumerate(X_seen_test):
        x = [0, 0]
        x[0] = X_test[0].to(device)
        x[1] = X_test[1].to(device)
        res = model(x)
        yy_seen.append(res[0].item())
    for idx, X_test in enumerate(X_unseen_test):
        x = [0, 0]
        x[0] = X_test[0].to(device)
        x[1] = X_test[1].to(device)
        res = model(x)
        yy_unseen.append(res[0].item())
    seen_rscr = roc_auc_score(y_seen_test, yy_seen)
    unseen_rscr = roc_auc_score(y_unseen_test, yy_unseen)
    if unseen_rscr>best_unseen:
        best_unseen = unseen_rscr
        ysts2 = yy_seen.copy()
        yusts2 = yy_unseen.copy()




## best pdbval

bb = rridx[(rridx.K>=8)|(rridx.K<=4)]
pdbss = list(bb.PDB)
btps ='no_both'
pock = pd.read_pickle(f'{dtip}/pdbval/sub_pockets_11_20_renew2_29pocks_kcluster_{btps}.pkl')
lts = pd.read_pickle(f'{dtip}/pdbval/AP1_{btps}.pkl')
apm = np.array(lts.loc[bb.cid,:]) 
apm = torch.tensor(apm, dtype = torch.float).view(len(bb),1,rn,10)
ppock = {}
for pds in pdbss:
    pcks = pock.loc[pds,0:549]
    ppock[pds] = torch.tensor(np.array(pcks).reshape([-1,rn,10]), dtype = torch.float)


best_pdbroc = 0.76
X = [[apm[i], ppock[pds]] for i,pds in enumerate(bb.PDB)]
lpocks = [len(ppock[pds]) for pds in pdbss]
rpocki = dict(zip(pdbss,[np.argmin(pock.loc[pds,:][550]) if pock.loc[pds,:][550][np.argmin(pock.loc[pds,:][550])]<=20 else '-' for pds in pdbss]))
for epc in range(10000):
    print(epc)
    rr1 = [];  ww1 = [];rr2 = [];   ww2 = [] ; yy =[]
    for idx, X_test in enumerate(X):
        x = [0, 0]
        x[0] = X_test[0].to(device)
        x[1] = X_test[1].to(device)
        res = model(x)
        att = res[1]
        pnn = X_test[1].shape[0]
        wws = att[0][0][pnn][0:pnn].cpu().detach().numpy()
        rrr=ss.rankdata(-wws)
        ww1.extend(wws)
        rr1.extend(rrr)
        wws = att[0][1][pnn][0:pnn].cpu().detach().numpy()
        rrr=ss.rankdata(-wws)
        ww2.extend(wws)
        rr2.extend(rrr)
        yy.append(res[0].item())
        if len(yy) % 100000 == 0:
            print(len(yy))
    pdbroc = roc_auc_score([0]*len(bb[bb.K<=4]) + [1]*len(bb[bb.K>=8]),yy)
    rpock = pock.loc[pdbss,:].copy()
    rpock['rank1'] = rr1
    rpock['wht1'] = ww1
    rpock['rank2'] = rr2
    rpock['wht2'] = ww2
    if pdbroc>best_pdbroc:
        best_pdbroc = pdbroc
        print(pdbroc)
        with open(f'{dtip}/pdbval/best_pdbroc.pkl','wb') as f:
            pickle.dump([rpock, yy, pdbroc],f)
 


# best chemble

# btps ='no_both'
# pock = pd.read_pickle(f'{dtip}/pdbval/sub_pockets_11_20_renew2_29pocks_kcluster_{btps}.pkl')
# lts = pd.read_pickle(f'{dtip}/pdbval/AP1_{btps}.pkl')





pdbss = list(set(pock.index))

ppock = {}
for pds in pdbss:
    pcks = pock.loc[pds,0:549]
    ppock[pds] = torch.tensor(np.array(pcks).reshape([-1,rn,10]), dtype = torch.float)


# best_pdbroc = 0.76


apm = np.array(lts.loc[cc.cid,:]) 
apm = torch.tensor(apm, dtype = torch.float).view(len(cc),1,rn,10)
cX = [[apm[i], ppock[pds]] for i,pds in enumerate(cc.To)]

uu = uu[uu.To.isin(pock.index)]
apm = np.array(lts.loc[uu.cid,:]) 
apm = torch.tensor(apm, dtype = torch.float).view(len(uu),1,rn,10)
uX = [[apm[i], ppock[pds]] for i,pds in enumerate(uu.To)]

best_ccroc = 0.80
best_uuroc = 0.73
for epc in range(10000):
    print(epc)
    yy =[]
    for idx, X_test in enumerate(cX):
        x = [0, 0]
        x[0] = X_test[0].to(device)
        x[1] = X_test[1].to(device)
        res = model(x)
        yy.append(res[0].item())
    ccroc = roc_auc_score(list(cc.label),yy)
    if ccroc>best_ccroc:
        best_ccroc = ccroc
        print(ccroc)
        with open(f'{chemm}/best_ccroc.pkl','wb') as f:
            pickle.dump([cc, yy, ccroc],f)
    yy =[]
    for idx, X_test in enumerate(uX):
        x = [0, 0]
        x[0] = X_test[0].to(device)
        x[1] = X_test[1].to(device)
        res = model(x)
        yy.append(res[0].item())
    uuroc = roc_auc_score(list(uu.label),yy)
    if uuroc>best_uuroc:
        best_uuroc = uuroc
        print(uuroc)
        with open(f'{chemm}/best_uuroc.pkl','wb') as f:
            pickle.dump([uu, yy, uuroc],f)





### chem old 

chemm = '/spstorage/USERS/gina/Project/AP/DB/ChEMBL'
uu=pd.read_csv(chemm+'/filtered_unseen_chembl.csv', index_col=0)
cc=pd.read_csv(chemm+'/filtered_seen_chembl.csv', index_col=0)


idx=1
lts = pd.read_pickle(chemm+'/AP1_no_both.pkl')

# pock = pd.read_pickle('{}/sub_pockets_{}_20_renew2_29pocks_kcluster.pkl'.format(chemm,pkn[idx]))
# pock.index = [abbb.split('_')[0] for abbb in list(pock.index)]
# # pock = pd.concat([pock[~pock.index.isin(pock[pock.iloc[:,550]==1].index)],pock[pock.iloc[:,550]==1]])
# apock = pd.read_pickle('{}/sub_pockets_11_20_renew2_29pocks_kcluster_no_both.pkl'.format(dtip))
# apock.index = [abbb[0:4].upper() for abbb in list(apock.index)]
# apock.columns = list(range(550))
# pock = pd.concat([pock, apock])

# ppock = pd.read_pickle('{}/pdbval/sub_pockets_11_20_renew2_29pocks_kcluster.pkl'.format(dtip))
# ppock.index = [abbb[0:4].upper() for abbb in list(ppock.index)]
# ppock = ppock.iloc[:,0:550]
# ppock.columns = list(range(550))
# ppock = ppock[~ppock.index.isin(pock.index)]
# pock = pd.concat([pock,ppock])

# list(set(uu.To)-set(pock.index))
# ppock = {}
# for pds in list(set(uu.To)):
#     pcks = pock.loc[pds,:]
#     if (len(pcks)>(pockn-1))&(len(pcks)!=rn*10):
#         pcks = pcks.sample(frac = 1, random_state = 7)[0:(pockn-1)]
#     ppock[pds]= torch.tensor(np.array(pcks).reshape([-1,rn,10]), dtype = torch.float)


pock = pd.read_pickle('{}/sub_pockets_{}_20_{}pocks_kcluster_less_than100.pkl'.format(chemm,pkn[idx], str(pockn)))
pock.index = [abbb.split('_')[0] for abbb in list(pock.index)]
# pock = pd.concat([pock[~pock.index.isin(pock[pock.iloc[:,550]==1].index)],pock[pock.iloc[:,550]==1]])
apock = pd.read_pickle('{}/sub_pockets_11_20_{}pocks_kcluster.pkl'.format(dtip, str(pockn)))
apock.index = [abbb[0:4].upper() for abbb in list(apock.index)]
apock.columns = list(range(550))
pock = pd.concat([pock, apock])
ppock = {}
for pds in list(set(list(cc.To)+list(uu.To))):
    pcks = pock.loc[pds,:]
    if (len(pcks)>(pockn-1))&(len(pcks)!=rn*10):
        pcks = pcks.sample(frac = 1, random_state = 7)[0:(pockn-1)]
    ppock[pds]= torch.tensor(np.array(pcks).reshape([-1,rn,10]), dtype = torch.float)
    



apm = np.array(lts.loc[cc.cid,:]) 
apm = torch.tensor(apm, dtype = torch.float).view(len(cc),1,rn,10)
cX = [[apm[i], ppock[pds]] for i,pds in enumerate(cc.To)]

uu = uu[uu.To.isin(pock.index)]
apm = np.array(lts.loc[uu.cid,:]) 
apm = torch.tensor(apm, dtype = torch.float).view(len(uu),1,rn,10)
uX = [[apm[i], ppock[pds]] for i,pds in enumerate(uu.To)]


best_ccroc = 0.80
best_uuroc = 0.73
for epc in range(10000):
    print(epc)
    yy =[]
    for idx, X_test in enumerate(cX):
        x = [0, 0]
        x[0] = X_test[0].to(device)
        x[1] = X_test[1].to(device)
        res = model(x)
        yy.append(res[0].item())
    ccroc = roc_auc_score(list(cc.label),yy)
    if ccroc>best_ccroc:
        best_ccroc = ccroc
        print(ccroc)
        with open(f'{chemm}/best_ccroc_old.pkl','wb') as f:
            pickle.dump([cc, yy, ccroc],f)
    yy =[]
    for idx, X_test in enumerate(uX):
        x = [0, 0]
        x[0] = X_test[0].to(device)
        x[1] = X_test[1].to(device)
        res = model(x)
        yy.append(res[0].item())
    uuroc = roc_auc_score(list(uu.label),yy)
    if uuroc>best_uuroc:
        best_uuroc = uuroc
        print(uuroc)
        with open(f'{chemm}/best_uuroc_old.pkl','wb') as f:
            pickle.dump([uu, yy, uuroc],f)


