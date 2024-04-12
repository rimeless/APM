import os, sys, datetime, pickle, math, argparse
import pandas as pd
import numpy as np
import itertools
from itertools import combinations_with_replacement, product
from operator import itemgetter, add

from rdkit import Chem
from rdkit import RDConfig
from rdkit.Chem import AllChem
from rdkit.Chem import ChemicalFeatures
from sklearn import metrics
from scipy import stats
import deepchem

pk = deepchem.dock.ConvexHullPocketFinder()



factory = ChemicalFeatures.BuildFeatureFactory('data/BaseFeatures_h.fdef')
fn = factory.GetFeatureFamilies() 


rdkit_morgan_radius=2
radius =2
gap = 0.6
inv = 1.2

##################

def dst(pr):
    xd = abs(pr[0][0]-pr[1][0])
    yd = abs(pr[0][1]-pr[1][1])
    zd = abs(pr[0][2]-pr[1][2])
    dst = math.sqrt(xd**2+yd**2+zd**2)
    return dst


def filter_feats_both(feat):
  dns = [f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')]
  feat = [f for f in feat if not (f.GetFamily()=='Acceptor')&(f.GetAtomIds() in dns)]
  return feat

def filter_feats_renew(feat):
  arh = [f.GetAtomIds() for f in feat if (f.GetFamily()=='Aromatic')]
  arh = list(itertools.chain(*arh))
  ions = [f.GetAtomIds() for f in feat if (f.GetFamily()=='NegIonizable')|(f.GetFamily()=='PosIonizable')]
  ions = list(itertools.chain(*ions))
  feat = [f for f in feat if not (f.GetFamily()=='Hydrophobe')&(any([ff in arh for ff in f.GetAtomIds()]))]
  feat = [f for f in feat if not (f.GetFamily()=='Donor')&(len(set(f.GetAtomIds()).intersection(ions))!=0)]
  feat = [f for f in feat if not (f.GetFamily()=='Acceptor')&(len(set(f.GetAtomIds()).intersection(ions))!=0)]
  feat = [f for f in feat if not (f.GetFamily()=='LumpedHydrophobe')&(any([ff in arh for ff in f.GetAtomIds()]))]
  return feat


sym2feat = {'C':7,'O':10,'N':10, 'S':10}
def feature_AH(row, feats, atms, mps, idx, col, inv, das, fn):
  # basic charac. extra atom to atom charac. 7 C hydrophobic / 10 O  N  S hydrophilic / 11 etc
  ft_idx = list(set(list(itertools.chain(*[f.GetAtomIds() for f in feats[idx]]))))
  uu_idx = [i for i in range(len(mps[idx])) if i not in ft_idx]
  use_pos = np.array([ft.GetPos() for ft in feats[idx]])
  use_fam = [fn.index(ft.GetFamily()) if not ft.GetAtomIds() in das[idx] else len(fn) for ft in feats[idx]]
  if (len(uu_idx) != 0) & (len(use_pos) !=0):
    use_idx = np.array(itemgetter(*uu_idx)(mps[idx]))
    if len(uu_idx)>1:
      use_atms = [sym2feat[at.GetSymbol()] if at.GetSymbol() in ['C','O','N','S'] else len(fn)+2 for at in itemgetter(*uu_idx)(atms[idx])]
      poses = np.concatenate((use_pos,use_idx))
    else:
      atsym = itemgetter(*uu_idx)(atms[idx]).GetSymbol()
      poses = np.concatenate((use_pos,np.array([use_idx])))
      if atsym in ['C','O','N','S']:
        use_atms = [sym2feat[atsym]]
      else:
        use_atms = [len(fn)+2]
    use_types =  use_fam + use_atms
  else:
    poses=use_pos
    use_types = use_fam
  ft=itertools.combinations(range(len(poses)), 2)
  lt=[0]*len(row)*len(col)
  for f in ft:
    dst = np.linalg.norm(poses[f[0]]-poses[f[1]]) 
    if dst!=0:
      cv= math.floor(math.log(dst/col[0], inv))+1
      if (cv>0) & (cv<len(col)-1):n1=col[cv-1];n2=col[cv+1]
      else:
        if cv<=0:n1=0;n2=col[0]
        else:n1=col[-1];n2=20     
      b=(stats.norm.cdf(n2,dst,0.5)-stats.norm.cdf(n1,dst,0.5))
      b1=stats.norm.cdf(n1,dst,0.5)
      b2=1-(stats.norm.cdf(n2,dst,0.5))
      Tpair=str(min(itemgetter(*f)(use_types)))+str(max(itemgetter(*f)(use_types)))
      if (dst!=0) & (Tpair in row) & (dst<20):
        tidx = row.index(Tpair)
        if dst<col[0]:cv=0
        elif dst>=col[-1]:cv=len(col)-1 
        if (cv>0) & (cv<len(col)-1):
          lt[tidx*len(col)+cv] += b
          lt[tidx*len(col)+cv-1] += b1
          lt[tidx*len(col)+cv+1] += b2
        else:
          if cv==0:
            lt[tidx*len(col)] += (b+b1)
            lt[tidx*len(col)+1] += b2
          else:
            lt[tidx*len(col)+cv] += (b+b2)
            lt[tidx*len(col)+cv-1] += b1
  return lt


def gen_featdf(mm, use_feats, mn):
  feats = [factory.GetFeaturesForMol(m) for m in mm]
  feats = [[f for f in feat if (f.GetFamily() in use_feats)] for feat in feats]
  feats = [filter_feats_renew(feat) for feat in feats]
  das = [[f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')] for feat in feats]
  das = [list(set([x for x in da if da.count(x) > 1])) for da in das] 
  feats = [filter_feats_both(feat) for feat in feats]
  podfs = []
  for idx in range(len(feats)):
	  use_pos = np.array([ft.GetPos() for ft in feats[idx]])
	  use_fam = [ft.GetFamily() if not ft.GetAtomIds() in das[idx] else 'D-A' for ft in feats[idx]]
	  podf = pd.DataFrame(use_pos)
	  podf['feat'] = use_fam
	  podf['cmp'] = mn[idx]
	  podfs.append(podf)
  podfs = pd.concat(podfs)
  return podfs

def gen_frag_df(frag_df, mm, mn):
    frag_mols = [Chem.MolFromSmiles(fg) for fg in list(frag_df.smiles)]
    frgs = pd.DataFrame([[1 if m.HasSubstructMatch(fg) else 0 for fg in frag_mols] for m in mm])
    mtcs = frgs.sum()[frgs.sum()==max(frgs.sum())].sort_index()
    frag_df = pd.concat([frag_df, frgs.T], axis=1)
    sub_bb = frag_df.loc[list(mtcs.index),:]
    sub_bb = sub_bb.sort_values('mw', ascending=False)
    sub_bb = sub_bb.sort_values('atom_num', ascending=False)
    sub_bb['Fragment'] = [f'static/fragments/efid_{i}.png' for i in list(sub_bb.efid)]
    sub_bb['from'] = [','.join([f'{str(mn[i])}' for i,x in enumerate(rcmp) if x==1]) for rcmp in sub_bb.iloc[:,6:].values]
    sub_bb = sub_bb.loc[:,['Fragment','from','smiles','name','class']]
    return sub_bb

def gen_AP(mm, use_feats, rns, cnrs, inv, fn):
  fns = [i for i, f in enumerate(fn) if f in use_feats]
  mn =  [m.GetProp('_Name') for m in mm]
  feats = [factory.GetFeaturesForMol(m) for m in mm]
  mps = [m.GetConformer().GetPositions() for m in mm] 
  atms = [m.GetAtoms() for m in mm]
  feats = [filter_feats_renew(feat) for feat in feats]
  fns = fns + [len(fn)]
  das = [[f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')] for feat in feats]
  das = [list(set([x for x in da if da.count(x) > 1])) for da in das] 
  feats = [filter_feats_both(feat) for feat in feats]
  fns = fns + list(range(len(fn)+1,len(fn)+3))
  rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
  vt = [feature_AH(rns, feats, atms, mps, idx, cnrs, inv, das, fn) for idx in range(len(mps))]
  vdf = pd.DataFrame(vt)
  vdf.index = mn
  return vdf


def gen_FP(mm):
  smiles_list = [Chem.MolToSmiles(m) for m in mm]
  mn =  [m.GetProp('_Name') for m in mm]
  fp = [list(Chem.rdMolDescriptors.GetMorganFingerprintAsBitVect(m, radius)) for m in mm if m is not None]
  fps = pd.DataFrame(fp)
  fps.index = mn
  return fps, smiles_list

def gen_pockets(pdb_file, pk, rns, cnrs, inv, fn, pockn):
    pockets = pk.find_pockets(pdb_file)
    pnm = os.path.splitext(os.path.basename(pdb_file))[0]
    # pocket df
    subcents = []
    pckts = []
    for i, pocket in enumerate(pockets):
        x_min = pocket.x_range[0]
        x_max = pocket.x_range[1]
        y_min = pocket.y_range[0]
        y_max = pocket.y_range[1]
        z_min = pocket.z_range[0]
        z_max = pocket.z_range[1]
        xg = (x_max - x_min)//10
        xrgs = [[sid, sid+10] for sid in range(int(x_min),int(x_max),round((x_max-10-x_min)/xg))[0:int(xg+1)]]
        yg = (y_max - y_min)//10
        yrgs = [[sid, sid+10] for sid in range(int(y_min),int(y_max),round((y_max-10-y_min)/yg))[0:int(yg+1)]]
        zg = (z_max - z_min)//10
        zrgs = [[sid, sid+10] for sid in range(int(z_min),int(z_max),round((z_max-10-z_min)/zg))[0:int(zg+1)]]
        for spck in list(itertools.product(*[xrgs,yrgs,zrgs])):
            cents = [np.mean(sp) for sp in spck]
            if not any([dst([cents, sbcts])<5 for sbcts in subcents]):
                pckts.append(list(itertools.chain(*spck)))
                subcents.append(cents)    
    subdf = pd.DataFrame(pckts)     
    ppm = Chem.rdmolfiles.MolFromPDBFile(pdb_file, removeHs=True)
    ppf = factory.GetFeaturesForMol(ppm)
    ppuse_pos = np.array([ft.GetPos() for ft in ppf])
    vts = []
    for ilen in range(len(subdf)):
        x_min, x_max, y_min, y_max, z_min, z_max= subdf.iloc[ilen, 0:6]
        feat = []
        for i, f in enumerate(ppf):
            if x_min < ppuse_pos[i][0] < x_max and y_min < ppuse_pos[i][1] < y_max and z_min < ppuse_pos[i][2] < z_max:
                feat.append(f)
        das = [[]]
        if len(feat)>5:
            feat = filter_feats_renew(feat)
            das = [f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')]
            das = list(set([x for x in das if das.count(x) > 1]))
            # feat = filter_feats_both(feat)
            mps = ppm.GetConformer().GetPositions()                         
            atms = [[am for i,am in enumerate(ppm.GetAtoms()) if x_min < mps[i][0] < x_max and y_min < mps[i][1] < y_max and z_min < mps[i][2] < z_max]]
            mps = [pm for pm in ppm.GetConformer().GetPositions() if x_min < pm[0] < x_max and y_min < pm[1] < y_max and z_min < pm[2] < z_max]
            mps = [np.stack(mps)]
            vt = feature_AH(rns, [feat], atms, mps, 0, cnrs, inv, [das], fn)
            vts.append(vt)
    pock = pd.DataFrame(vts).drop_duplicates()
    pock = pock[(pock.T.sum()>10) & (pock.T.sum()<100)]
    nsubss = [] 
    for pid in list(set(pock.index)):
        subpock = pock.loc[pid, :]
        if len(subpock.shape)==1:
            nsubdf = pd.DataFrame(subpock).T
        elif len(subpock)>pockn:
            k_means = KMeans(init="k-means++", n_clusters=pockn, n_init=20)#, n_jobs = 50)
            k_means.fit(subpock)
            nsubs = []
            for i in range(pockn):
                iidx = [ii for ii, pp in enumerate(k_means.labels_) if pp == i]
                ssub = subpock.iloc[iidx,:]
                dsts = [np.linalg.norm(ssub.iloc[a,:].to_numpy()-k_means.cluster_centers_[0]) for a in range(len(ssub))]
                nsubs.append(ssub.iloc[dsts.index(min(dsts)),:])
            nsubdf = pd.DataFrame(nsubs)
        else:
            nsubdf = subpock
        nsubss.append(nsubdf)
    pcks = pd.concat(nsubss)
    pcks.index = pnm
    return pcks



def generate_APM(input_type, input_path, dstbin, use_feats, pockn):
  rto = 20/(1.2**dstbin)
  cnrs = [rto*(1.2**g) for g in range(dstbin)]
  feat_types = use_feats.split(',')
  rn = int((len(feat_types) + 3)*(len(feat_types) +4)/2)
  fns = [fn.index(f) for f in feat_types] + [9, 10, 11]
  rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
  if input_type=='compound':
    mm = Chem.SDMolSupplier(input_path, removeHs=True) 
    mm = [m for m in mm if m is not None]
    apm = gen_AP(mm, feat_types, rns, cnrs, inv, fn)
  else:
    apm = gen_pockets(input_path, pk, rns, cnrs, inv, fn, pockn)  
  return apm
