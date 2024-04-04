
import sys
src_dir = '/spstorage/USERS/gina/source' 
sys.path.append(src_dir)
from sklearn.cluster import AgglomerativeClustering
from basic import *
from APM import *
from operator import itemgetter, add
import glob
from collections import OrderedDict

gap = 0.6



def filter_feats_both(feat):
	dns = [f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')]
	feat = [f for f in feat if not (f.GetFamily()=='Acceptor')&(f.GetAtomIds() in dns)]
	return feat



def filter_feats_dup(feat):
	dns = [f.GetAtomIds() for f in feat if (f.GetFamily()=='Acceptor')]
	feat = [f for f in feat if not (f.GetFamily()=='Donor')&(f.GetAtomIds() in dns)]
	arm = [f.GetAtomIds() for f in feat if (f.GetFamily()=='Aromatic')]
	arm = [a for b in arm for a in b]
	feat = [f for f in feat if not (f.GetFamily() in ['Donor','PosIonizable'])&(len(intersect(set(f.GetAtomIds()),arm))!=0)]
	return feat

# def filter_feats_hyp(feat):
# 	arh = [f.GetAtomIds() for f in feat if (f.GetFamily()=='Aromatic')]
# 	feat = [f for f in feat if not (f.GetFamily()=='Hydrophobe')&(f.GetAtomIds() in arh)]
# 	return feat


def filter_feats_renew(feat):
	arh = [f.GetAtomIds() for f in feat if (f.GetFamily()=='Aromatic')]
	hydp = [f.GetAtomIds() for f in feat if (f.GetFamily()=='Hydrophobe')]
	ions = [f.GetAtomIds() for f in feat if (f.GetFamily()=='NegIonizable')|(f.GetFamily()=='PosIonizable')]
	ions = list(itertools.chain(*ions))
	feat = [f for f in feat if not (f.GetFamily()=='Hydrophobe')&(f.GetAtomIds() in arh)]
	feat = [f for f in feat if not (f.GetFamily()=='Donor')&(len(intersect(set(f.GetAtomIds()),ions))!=0)]
	feat = [f for f in feat if not (f.GetFamily()=='Acceptor')&(len(intersect(set(f.GetAtomIds()),ions))!=0)]
	return feat


def filter_feats_renew2(feat):
	arh = [f.GetAtomIds() for f in feat if (f.GetFamily()=='Aromatic')]
	arh = list(itertools.chain(*arh))
	ions = [f.GetAtomIds() for f in feat if (f.GetFamily()=='NegIonizable')|(f.GetFamily()=='PosIonizable')]
	ions = list(itertools.chain(*ions))
	feat = [f for f in feat if not (f.GetFamily()=='Hydrophobe')&(any([ff in arh for ff in f.GetAtomIds()]))]
	feat = [f for f in feat if not (f.GetFamily()=='Donor')&(len(intersect(set(f.GetAtomIds()),ions))!=0)]
	feat = [f for f in feat if not (f.GetFamily()=='Acceptor')&(len(intersect(set(f.GetAtomIds()),ions))!=0)]
	feat = [f for f in feat if not (f.GetFamily()=='LumpedHydrophobe')&(any([ff in arh for ff in f.GetAtomIds()]))]
	return feat


# def filter_feats_renew3(feat):
# 	arh = [f.GetAtomIds() for f in feat if (f.GetFamily()=='Aromatic')]
# 	arh = list(itertools.chain(*arh))
# 	ions = [f.GetAtomIds() for f in feat if (f.GetFamily()=='NegIonizable')|(f.GetFamily()=='PosIonizable')]
# 	ions = list(itertools.chain(*ions))
# 	dns = [f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')]
# 	feat = [f for f in feat if not (f.GetFamily()=='Acceptor')&(f.GetAtomIds() in dns)]
# 	feat = [f for f in feat if not (f.GetFamily()=='Hydrophobe')&(any([ff in arh for ff in f.GetAtomIds()]))]
# 	feat = [f for f in feat if not (f.GetFamily()=='Donor')&(len(intersect(set(f.GetAtomIds()),ions))!=0)]
# 	feat = [f for f in feat if not (f.GetFamily()=='Acceptor')&(len(intersect(set(f.GetAtomIds()),ions))!=0)]
# 	feat = [f for f in feat if not (f.GetFamily()=='LumpedHydrophobe')&(any([ff in arh for ff in f.GetAtomIds()]))]
# 	return feat



def filter_feats_renew3(feat):
	arh = [f.GetAtomIds() for f in feat if (f.GetFamily()=='Aromatic')]
	arh = list(itertools.chain(*arh))
	feat = [f for f in feat if not (f.GetFamily()!='Aromatic')&(any([ff in arh for ff in f.GetAtomIds()]))]
	halo = 	[f.GetAtomIds() for f in feat if (f.GetFamily()=='Halogen')]
	halo = list(itertools.chain(*halo))
	feat = [f for f in feat if not (f.GetFamily()=='Acceptor')&(f.GetAtomIds() in halo)]
	ions = [f.GetAtomIds() for f in feat if (f.GetFamily()=='NegIonizable')|(f.GetFamily()=='PosIonizable')]
	ions = list(itertools.chain(*ions))
	feat = [f for f in feat if not (f.GetFamily() in ['Donor', 'Acceptor', 'LumpedHydrophobe'])&(len(intersect(set(f.GetAtomIds()),ions))!=0)]
	lump = 	[f.GetAtomIds() for f in feat if (f.GetFamily()=='LumpedHydrophobe')]
	lump = list(itertools.chain(*lump))	
	feat = [f for f in feat if not (f.GetFamily()=='Hydrophobe')&(any([ff in lump for ff in f.GetAtomIds()]))]
	return feat



def filter_feats_renew_posarm_all(feat):
	arh = [f.GetAtomIds() for f in feat if (f.GetFamily()=='Aromatic')]
	arh = list(itertools.chain(*arh))
	feat = [f for f in feat if not (f.GetFamily() in ['Hydrophobe','Donor', 'Acceptor', 'LumpedHydrophobe'])&(any([ff in arh for ff in f.GetAtomIds()]))]
	halo = 	[f.GetAtomIds() for f in feat if (f.GetFamily()=='Halogen')]
	halo = list(itertools.chain(*halo))
	feat = [f for f in feat if not (f.GetFamily()=='Acceptor')&(f.GetAtomIds() in halo)]
	ions = [f.GetAtomIds() for f in feat if (f.GetFamily()=='NegIonizable')|(f.GetFamily()=='PosIonizable')]
	ions = list(itertools.chain(*ions))
	feat = [f for f in feat if not (f.GetFamily() in ['Donor', 'Acceptor', 'LumpedHydrophobe'])&(len(intersect(set(f.GetAtomIds()),ions))!=0)]
	lump = 	[f.GetAtomIds() for f in feat if (f.GetFamily()=='LumpedHydrophobe')]
	lump = list(itertools.chain(*lump))	
	feat = [f for f in feat if not (f.GetFamily()=='Hydrophobe')&(any([ff in lump for ff in f.GetAtomIds()]))]
	return feat

def pairLr_dist(fts, row, col, inv, das, fn):
	ft=itertools.combinations(fts, 2)
	lt=[0]*len(row)*len(col)
	for f in ft:
		dst=cal_dst(f)	
		if dst!=0:
			cv= math.floor(math.log(dst/col[0], inv))+1
			if (cv>0) & (cv<len(col)-1):n1=col[cv-1];n2=col[cv+1]
			else:
				if cv<0:n1=0;n2=col[0]
				else:n1=col[-1];n2=20			
			b=(stats.norm.cdf(n2,dst,0.5)-stats.norm.cdf(n1,dst,0.5))
			b1=stats.norm.cdf(n1,dst,0.5)
			b2=1-(stats.norm.cdf(n2,dst,0.5))
			aType=fn.index(f[0].GetFamily())
			if f[0].GetAtomIds() in das:aType=9
			bType=fn.index(f[1].GetFamily())
			if f[1].GetAtomIds() in das:bType=9
			Tpair=str(min(aType,bType))+str(max(aType,bType))
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


# def bulk_cids_ratio(sdfs):
# 	mm = Chem.SDMolSupplier(sdfs, removeHs=True) 
# 	mn = [m.GetPropsAsDict()['PUBCHEM_COMPOUND_CID'] for m in mm]
# 	das = [f.GetAtomIds() for f in feats if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')]
# 	das = unique([x for x in das if das.count(x) > 1])
# 	feats = [filter_feats_hyp_acc(m) for m in mm]
# 	vt = [pairLr_dist(f, rn6, cn_r, inv, das) for f in feats]
# 	vdf = pd.DataFrame(vt)
# 	vdf.index = mn
# 	return vdf

# def MMFFLr_dist_gap_feats(feats, mmffp, mps, idx, col, gap, inv, das):
# 	row = [str(x)+str(y) for x,y in combinations_with_replacement([0,1,2,3,5,6,7,8,9,10,11,12,13],2)]
# 	use_pos = np.array([ft.GetPos() for ft in feats[idx]])
# 	use_fam = [fn.index(ft.GetFamily()) for ft in feats[idx]]
# 	mmffps = [m for m in mmffp[idx][1::] if len(mps[idx])>=int(m.split(' ')[0])]
# 	use_idx = np.array([mps[idx][int(m.split(' ')[0])-1] for m in mmffps])# if bool(mps[idx][int(m.split(' ')[0])-1] in use_pos)])# not in f_idx]# if int(m.split(' ')[0])<mhac[idx]]
# 	use_mmff = [int(float(m.split(' ')[1])//gap + 11) for m in mmffp[idx][1::]]# if bool(mps[idx][int(m.split(' ')[0])-1] in use_pos)]# if int(m.split(' ')[0])<mhac[idx]]
# 	if len(use_idx) != 0:poses = np.concatenate((use_pos,use_idx))
# 	else:poses=use_pos
# 	use_types = use_fam + use_mmff
# 	ft=itertools.combinations(range(len(poses)), 2)
# 	lt=[0]*len(row)*len(col)
# 	for f in ft:
# 		dst = np.linalg.norm(poses[f[0]]-poses[f[1]])	
# 		if dst!=0:
# 			cv= math.floor(math.log(dst/col[0], inv))+1
# 			if (cv>0) & (cv<len(col)-1):n1=col[cv-1];n2=col[cv+1]
# 			else:
# 				if cv<0:n1=0;n2=col[0]
# 				else:n1=col[-1];n2=20			
# 			b=(stats.norm.cdf(n2,dst,0.5)-stats.norm.cdf(n1,dst,0.5))
# 			b1=stats.norm.cdf(n1,dst,0.5)
# 			b2=1-(stats.norm.cdf(n2,dst,0.5))
# 			Tpair=str(min(itemgetter(*f)(use_types)))+str(max(itemgetter(*f)(use_types)))
# 			if (dst!=0) & (Tpair in row) & (dst<20):
# 				tidx = row.index(Tpair)
# 				if dst<col[0]:cv=0
# 				elif dst>=col[-1]:cv=len(col)-1	
# 				if (cv>0) & (cv<len(col)-1):
# 					lt[tidx*len(col)+cv] += b
# 					lt[tidx*len(col)+cv-1] += b1
# 					lt[tidx*len(col)+cv+1] += b2
# 				else:
# 					if cv==0:
# 						lt[tidx*len(col)] += (b+b1)
# 						lt[tidx*len(col)+1] += b2
# 					else:
# 						lt[tidx*len(col)+cv] += (b+b2)
# 						lt[tidx*len(col)+cv-1] += b1
# 	return lt


gap=0.6
def feature_PP(row,feats, mmffp, mps, idx, col, gap, inv, das, fn):
	use_pos = np.array([ft.GetPos() for ft in feats[idx]])
	use_fam = [fn.index(ft.GetFamily()) if not ft.GetAtomIds() in das[idx] else len(fn) for ft in feats[idx]]
	mmffps = [m for m in mmffp[idx][1::] if len(mps[idx])>=int(m.split(' ')[0])]
	use_idx = np.array([mps[idx][int(m.split(' ')[0])-1] for m in mmffps if not bool(mps[idx][int(m.split(' ')[0])-1] in use_pos)])# not in f_idx]# if int(m.split(' ')[0])<mhac[idx]]
	use_mmff = [int(float(m.split(' ')[1])//gap + len(fn)+3) for m in mmffps if not bool(mps[idx][int(m.split(' ')[0])-1] in use_pos)]# if int(m.split(' ')[0])<mhac[idx]]
	if (len(use_idx) != 0) & (len(use_pos) !=0):poses = np.concatenate((use_pos,use_idx))
	else:poses=use_pos
	use_types = use_fam + use_mmff
	ft=itertools.combinations(range(len(poses)), 2)
	lt=[0]*len(row)*len(col)
	for f in ft:
		dst = np.linalg.norm(poses[f[0]]-poses[f[1]])	
		if dst!=0:
			cv= math.floor(math.log(dst/col[0], inv))+1
			if (cv>0) & (cv<len(col)-1):n1=col[cv-1];n2=col[cv+1]
			else:
				if cv<0:n1=0;n2=col[0]
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

sym2feat = {'C':10,'O':11,'N':12, 'S':13}
def feature_AA(row, feats, atms, mps, idx, col, inv, das, fn):
	# basic charac. extra atom to atom charac. 10 C / 11 O / 12 N / 13 S / 14 etc
	ft_idx = unique(list(itertools.chain(*[f.GetAtomIds() for f in feats[idx]])))
	uu_idx = [i for i in range(len(mps[idx])) if i not in ft_idx]
	use_pos = np.array([ft.GetPos() for ft in feats[idx]])
	use_fam = [fn.index(ft.GetFamily()) if not ft.GetAtomIds() in das[idx] else len(fn) for ft in feats[idx]]
	if (len(uu_idx) != 0) & (len(use_pos) !=0):
		use_idx = np.array(itemgetter(*uu_idx)(mps[idx]))
		if len(uu_idx)>1:
			use_atms = [sym2feat[at.GetSymbol()] if at.GetSymbol() in ['C','O','N','S'] else len(fn)+5 for at in itemgetter(*uu_idx)(atms[idx])]
			poses = np.concatenate((use_pos,use_idx))
		else:
			atsym = itemgetter(*uu_idx)(atms[idx]).GetSymbol()
			poses = np.concatenate((use_pos,np.array([use_idx])))
			if atsym in ['C','O','N','S']:
				use_atms = [sym2feat[atsym]]
			else:
				use_atms = [len(fn)+5]
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
				if cv<0:n1=0;n2=col[0]
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


sym2feat = {'C':10,'O':11,'N':12, 'S':13}
def feature_AP(row, feats, atms, mmffp, gap, mps, idx, col, inv, das, fn):
	# basic charac. extra atom to atom charac. 10 C / 11 O / 12 N / 13 S / 14 etc
	ft_idx = unique(list(itertools.chain(*[f.GetAtomIds() for f in feats[idx]])))
	uu_idx = [i for i in range(len(mps[idx])) if i not in ft_idx]
	use_pos = np.array([ft.GetPos() for ft in feats[idx]])
	use_fam = [fn.index(ft.GetFamily()) if not ft.GetAtomIds() in das[idx] else len(fn) for ft in feats[idx]]
	if (len(uu_idx) != 0) & (len(use_pos) !=0):
		use_idx = np.array(itemgetter(*uu_idx)(mps[idx]))
		if len(uu_idx)>1:
			use_atms = [sym2feat[at.GetSymbol()] for at in itemgetter(*uu_idx)(atms[idx]) if at.GetSymbol() in ['C','O','N','S']]
			use_idx = [use_idx[i] for i,at in enumerate(itemgetter(*uu_idx)(atms[idx])) if at.GetSymbol() in ['C','O','N','S']]
			if len(use_atms)>0:
				poses = np.concatenate((use_pos,use_idx))
			else:
				poses = use_pos
		else:
			atsym = itemgetter(*uu_idx)(atms[idx]).GetSymbol()
			if atsym in ['C','O','N','S']:
				use_atms = [sym2feat[atsym]]
				poses = np.concatenate((use_pos,np.array([use_idx])))
			else:
				use_atms = []
				poses = use_pos
		use_types =  use_fam + use_atms
	else:
		poses=use_pos
		use_types = use_fam
	mmffps = [m for m in mmffp[idx][1::] if len(mps[idx])>=int(m.split(' ')[0])]
	use_idx = np.array([mps[idx][int(m.split(' ')[0])-1] for m in mmffps if not bool(mps[idx][int(m.split(' ')[0])-1] in poses)])# not in f_idx]# if int(m.split(' ')[0])<mhac[idx]]
	use_mmff = [int(float(m.split(' ')[1])//gap + len(fn)+7) for m in mmffps if not bool(mps[idx][int(m.split(' ')[0])-1] in poses)]# if int(m.split(' ')[0])<mhac[idx]]
	if (len(use_idx) != 0) & (len(use_pos) !=0):poses = np.concatenate((poses,use_idx))
	use_types = use_types + use_mmff
	ft=itertools.combinations(range(len(poses)), 2)
	lt=[0]*len(row)*len(col)
	for f in ft:
		dst = np.linalg.norm(poses[f[0]]-poses[f[1]])	
		if dst!=0:
			cv= math.floor(math.log(dst/col[0], inv))+1
			if (cv>0) & (cv<len(col)-1):n1=col[cv-1];n2=col[cv+1]
			else:
				if cv<0:n1=0;n2=col[0]
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



sym2feat3 = {'P':10, 'S':10}
def feature_AH2(row, feats, atms, mps, idx, col, inv, das, fn):
	# basic charac. extra atom to atom charac. 7 C hydrophobic / 10 O  N  S hydrophilic / 11 etc
	ft_idx = unique(list(itertools.chain(*[f.GetAtomIds() for f in feats[idx]])))
	uu_idx = [i for i in range(len(mps[idx])) if i not in ft_idx]
	use_pos = np.array([ft.GetPos() for ft in feats[idx]])
	use_fam = [fn.index(ft.GetFamily()) if not ft.GetAtomIds() in das[idx] else len(fn) for ft in feats[idx]]
	if (len(uu_idx) != 0) & (len(use_pos) !=0):
		use_idx = np.array(itemgetter(*uu_idx)(mps[idx]))
		use_atms = []
		if len(uu_idx)>1:
			use_uidx = []
			for uidx, at in enumerate(itemgetter(*uu_idx)(atms[idx])):
				if at.GetSymbol() in ['P','S']:
					use_atms.append(10)
					use_uidx.append(uidx)
				elif len(at.GetBonds())==1:
					use_atms.append(11)
					use_uidx.append(uidx)
			if len(use_uidx)>1:
				use_idx = np.array(itemgetter(*use_uidx)(use_idx))
				#use_atms = [10 if at.GetSymbol() in ['P','S'] else len(fn)+2 for at in itemgetter(*uu_idx)(atms[idx])]
				poses = np.concatenate((use_pos,use_idx))
			elif len(use_uidx)==1:
				use_idx = np.array(itemgetter(*use_uidx)(use_idx))
				poses = np.concatenate((use_pos,np.array([use_idx])))
			else:
				poses=use_pos
				use_types = use_fam
		else:
			at = itemgetter(*uu_idx)(atms[idx])
			if at.GetSymbol() in ['P','S']:
				use_atms = [10]
				poses = np.concatenate((use_pos,np.array([use_idx])))
			elif len(at.GetBonds())==1:
				use_atms = [11]
				poses = np.concatenate((use_pos,np.array([use_idx])))
			else:
				poses=use_pos
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
				#if cv<0:n1=0;n2=col[0]
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


sym2feat2 = {'C':7,'O':10,'N':10, 'S':10}
def feature_AH(row, feats, atms, mps, idx, col, inv, das, fn):
	# basic charac. extra atom to atom charac. 7 C hydrophobic / 10 O  N  S hydrophilic / 11 etc
	ft_idx = unique(list(itertools.chain(*[f.GetAtomIds() for f in feats[idx]])))
	uu_idx = [i for i in range(len(mps[idx])) if i not in ft_idx]
	use_pos = np.array([ft.GetPos() for ft in feats[idx]])
	use_fam = [fn.index(ft.GetFamily()) if not ft.GetAtomIds() in das[idx] else len(fn) for ft in feats[idx]]
	if (len(uu_idx) != 0) & (len(use_pos) !=0):
		use_idx = np.array(itemgetter(*uu_idx)(mps[idx]))
		if len(uu_idx)>1:
			use_atms = [sym2feat2[at.GetSymbol()] if at.GetSymbol() in ['C','O','N','S'] else len(fn)+2 for at in itemgetter(*uu_idx)(atms[idx])]
			poses = np.concatenate((use_pos,use_idx))
		else:
			atsym = itemgetter(*uu_idx)(atms[idx]).GetSymbol()
			poses = np.concatenate((use_pos,np.array([use_idx])))
			if atsym in ['C','O','N','S']:
				use_atms = [sym2feat2[atsym]]
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
				#if cv<0:n1=0;n2=col[0]
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


sym2feat2 = {'C':7,'O':10,'N':10, 'S':10}
def feature_AHP(row, feats, atms, mps, idx, col, inv, das, posarm, fn):
	# basic charac. extra atom to atom charac. 7 C hydrophobic / 10 O  N  S hydrophilic / 11 etc / 12 negative aromatic
	ft_idx = unique(list(itertools.chain(*[f.GetAtomIds() for f in feats[idx]])))
	uu_idx = [i for i in range(len(mps[idx])) if i not in ft_idx]
	use_pos = np.array([ft.GetPos() for ft in feats[idx]])
	# use_fam = [fn.index(ft.GetFamily()) if not ft.GetAtomIds() in das[idx] else len(fn) if not any([x in posarm[idx] for x in ft.GetAtomIds()]) else 12 for ft in feats[idx]]
	use_fam = [len(fn) if ft.GetAtomIds() in das[idx] else fn.index(ft.GetFamily()) if not any([x in posarm[idx] for x in ft.GetAtomIds()]) else 12 for ft in feats[idx]]
	# use_fam = [fn.index(ft.GetFamily()) for ft in use_fam]
	if (len(uu_idx) != 0) & (len(use_pos) !=0):
		use_idx = np.array(itemgetter(*uu_idx)(mps[idx]))
		if len(uu_idx)>1:
			use_atms = [sym2feat2[at.GetSymbol()] if at.GetSymbol() in ['C','O','N','S'] else len(fn)+2 for at in itemgetter(*uu_idx)(atms[idx])]
			poses = np.concatenate((use_pos,use_idx))
		else:
			atsym = itemgetter(*uu_idx)(atms[idx]).GetSymbol()
			poses = np.concatenate((use_pos,np.array([use_idx])))
			if atsym in ['C','O','N','S']:
				use_atms = [sym2feat2[atsym]]
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
				if cv<0:n1=0;n2=col[0]
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


def bulk_cids(ssets):
	sdf, fact, rns, cnrs, fn = ssets
	mm = Chem.SDMolSupplier(sdf, removeHs=True) 
	mn = [m.GetPropsAsDict()['PUBCHEM_COMPOUND_CID'] for m in mm if m is not None]
	feats = [factory.GetFeaturesForMol(m) for m in mm if m is not None]
	das = [[]]*len(feats)
	if fact > 0:
		feats = [filter_feats_renew(feat) for feat in feats]
		if fact > 1:
			das = [[f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')] for feat in feats]
			das = [unique([x for x in da if da.count(x) > 1]) for da in das]
			if fact > 4:				
				if fact > 5:
					atms = [m.GetAtoms() for m in mm if m is not None]
					mps = [m.GetConformer().GetPositions() for m in mm if m is not None]
					if fact==7:
						mpd = [m.GetPropsAsDict() for m in mm if m is not None]
						mmffp = [m['PUBCHEM_MMFF94_PARTIAL_CHARGES'] if 'PUBCHEM_MMFF94_PARTIAL_CHARGES' in m.keys() else '' for m in mpd]
						mmffp = [m.split('\n') if (m!=0) & (m!='') else '' for m in mmffp]
						vt = [feature_PP(feats, mmffp, mps, idx, cnrs, gap, inv, das, fn) for idx in range(len(mpd))]				
					else:
						vt = [feature_AA(rns, feats, atms, mps, idx, cnrs, inv, das, fn) for idx in range(len(mpd))]
				else:
					mpd = [m.GetPropsAsDict() for m in mm if m is not None]
					mmffp = [m['PUBCHEM_MMFF94_PARTIAL_CHARGES'] if 'PUBCHEM_MMFF94_PARTIAL_CHARGES' in m.keys() else '' for m in mpd]
					mmffp = [m.split('\n') if (m!=0) & (m!='') else '' for m in mmffp]		
					vt = [feature_AP(rns, feats, atms, mmffp, gap,mps, idx, cnrs, inv, das, fn) for idx in range(len(mpd))]
			else:
				vt = [pairLr_dist(f, rns, cnrs, inv, das[i], fn) for i, f in enumerate(feats)]
		else:
			vt = [pairLr_dist(f, rns, cnrs, inv, das[i], fn) for i, f in enumerate(feats)]
	else:
		vt = [pairLr_dist(f, rns, cnrs, inv, das[i], fn) for i, f in enumerate(feats)]
	vdf = pd.DataFrame(vt)
	vdf.index = mn
	return vdf

dstbin = 10
f_list=['','','','_all'] +['_h']*8
st = 20/(inv**dstbin) # 16 bin
cnrs = [st*(inv**g) for g in range(dstbin)]
fnss = [[0,1,2,3,5,6],[0,1,2,3,5,6],[0,1,2,3,5,6,9],[0,1,2,3,4,5,7,8,10],[0,1,2,3,4,6,7,9],[0,1,2,3,4,6,7,9,10,11,12,13,14],[0,1,2,3,4,6,7,9,10,11,12,13,14],[0,1,2,3,4,6,7,9,10,11,12,13,14,15,16,17,18],[0,1,2,3,4,6,7,9,10,11,12,13,14],[0,1,2,3,4,6,7,9,10,11,12,13,14],[0,1,2,3,4,6,7,9,10,11,12,13,14,15,16,17,18],[0,1,2,3,4,6,7,9,10,11]]




cid_list = 
use_feats = ['Donor', 'Acceptor', 'NegIonizable', 'PosIonizable', 'Halogen', 'Aromatic', 'Hydrophobe']
use_feats = ['Donor', 'Acceptor', 'NegIonizable', 'PosIonizable', 'Aromatic', 'Hydrophobe']
use_feats = ['Donor', 'Acceptor', 'NegIonizable', 'PosIonizable', 'Halogen', 'Aromatic', 'Hydrophobe', 'LumpedHydrophobe']
# D+A type use / AH / AA / PP / AP
tps = [1,0,1,0,0] # 2
tps = [1,1,0,0,0] # 1
# tps = [1,0,0,0,0] # 0
tps = [0,0,0,0,0] # 0

str_dir = '/spstorage/DB/PUBCHEM/3D_1conf'
str_dir_files = glob.glob(f'{str_dir}/*')

cidss = [[] for i in range(6702)]
for cid in cid_list:
	if cid < 167550000:
		cidss[int(cid // 25000)].append(cid)


sdfs = ['{}/{}_{}.sdf'.format(str_dir,str(i * 25000 +1).zfill(8), str((i + 1) * 25000).zfill(8)) for i,cd in enumerate(cidss)
 if len(cd)!=0]
cds = [cd for cd in cidss if len(cd)!=0]

infos = [[cds[i], sdfs[i], tps] for i in range(len(cds)) if sdfs[i] in str_dir_files]


start_time = time.time()
ltiss3 = mclapply(infos,AP_from_cid, ncpu= 55) # 100 / 10 855s   / 700 / 20 --> 57m | 3466s
print("---{}s seconds---".format(time.time()-start_time))

lltss2 = pd.concat(ltiss)
lltss2.to_pickle(dtip+'/data/AP1_no_both_left.pkl')



gap = 0.6
def AP_from_cid(info):
	cids, sdf, tps = info
	factory = ChemicalFeatures.BuildFeatureFactory('{}/etc/BaseFeatures_h.fdef'.format(DIR_AP))
	fn = factory.GetFeatureFamilies()	
	fns = [i for i, f in enumerate(fn) if f in use_feats]
	mm = Chem.SDMolSupplier(sdf, removeHs=True) 
	mm = [m for m in mm if m is not None]
	mm = [m for m in mm if m.GetPropsAsDict()['PUBCHEM_COMPOUND_CID'] in cids]
	mn =  [m.GetPropsAsDict()['PUBCHEM_COMPOUND_CID'] for m in mm]
	feats = [factory.GetFeaturesForMol(m) for m in mm]
	if sum(tps[1:5])>0:
		mps = [m.GetConformer().GetPositions() for m in mm]	
		if tps[3]!=1:
			atms = [m.GetAtoms() for m in mm]
		if sum(tps[3:5])>0:
			mpd = [m.GetPropsAsDict() for m in mm]
			mmffp = [m['PUBCHEM_MMFF94_PARTIAL_CHARGES'] if 'PUBCHEM_MMFF94_PARTIAL_CHARGES' in m.keys() else '' for m in mpd]
			mmffp = [m.split('\n') if (m!=0) & (m!='') else '' for m in mmffp]		
	# featur n selection
	if tps[0]==1:
		feats = [filter_feats_renew_posarm_all(feat) for feat in feats]
		fns = fns + [len(fn)]
		das = [[f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')] for feat in feats]
		das = [unique([x for x in da if da.count(x) > 1]) for da in das]	
		posarm = [[f.GetAtomIds() for f in feat if (f.GetFamily()=='PosIonizable')|(f.GetFamily()=='Aromatic')] for feat in feats]
		posarm = [[a for b in posa for a in b] for posa in posarm]
		posarm = [unique([x for x in da if da.count(x) > 1]) for da in posarm]	
		feats = [filter_feats_dup(feat) for feat in feats]
	if tps[1]==1:
		fns = fns + list(range(len(fn)+1,len(fn)+4))
		rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
		vt = [feature_AHP(rns, feats, atms, mps, idx, cnrs, inv, das, posarm,fn) for idx in range(len(mps))]
	elif sum(tps[2:4])>0:
		fns = fns + list(range(len(fn)+1,len(fn)+6))
		rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
		if tps[2]==1:
			vt = [feature_AA(rns, feats, atms, mps, idx, cnrs, inv, das, fn) for idx in range(len(mps))]
		else:
			vt = [feature_PP(rns,feats, mmffp, mps, idx, cnrs, gap, inv, das, fn) for idx in range(len(mpd))]
	elif tps[4]==1:
		fns = fns + list(range(len(fn)+1,len(fn)+10))
		rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
		vt = [feature_AP(rns, feats, atms, mmffp, gap,mps, idx, cnrs, inv, das, fn) for idx in range(len(mpd))]			
	else:
		rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
		vt = [pairLr_dist(f, rns, cnrs, inv, das[i], fn) for i, f in enumerate(feats)]
	vdf = pd.DataFrame(vt)
	vdf.index = mn
	return vdf

gap = 0.6
def AP_from_cid(info):
	cids, sdf, tps = info
	factory = ChemicalFeatures.BuildFeatureFactory('{}/etc/BaseFeatures_h.fdef'.format(DIR_AP))
	fn = factory.GetFeatureFamilies()	
	fns = [i for i, f in enumerate(fn) if f in use_feats]
	mm = Chem.SDMolSupplier(sdf, removeHs=True) 
	mm = [m for m in mm if m is not None]
	mm = [m for m in mm if m.GetPropsAsDict()['PUBCHEM_COMPOUND_CID'] in cids]
	mn =  [m.GetPropsAsDict()['PUBCHEM_COMPOUND_CID'] for m in mm]
	feats = [factory.GetFeaturesForMol(m) for m in mm]
	if sum(tps[1:5])>0:
		mps = [m.GetConformer().GetPositions() for m in mm]	
		if tps[3]!=1:
			atms = [m.GetAtoms() for m in mm]
		if sum(tps[3:5])>0:
			mpd = [m.GetPropsAsDict() for m in mm]
			mmffp = [m['PUBCHEM_MMFF94_PARTIAL_CHARGES'] if 'PUBCHEM_MMFF94_PARTIAL_CHARGES' in m.keys() else '' for m in mpd]
			mmffp = [m.split('\n') if (m!=0) & (m!='') else '' for m in mmffp]		
	# featur n selection
	if tps[0]==1:
		feats = [filter_feats_renew2(feat) for feat in feats]
		fns = fns + [len(fn)]
		das = [[f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')] for feat in feats]
		das = [unique([x for x in da if da.count(x) > 1]) for da in das]	
		feats = [filter_feats_both(feat) for feat in feats]
	if tps[1]==1:
		fns = fns + list(range(len(fn)+1,len(fn)+3))
		rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
		vt = [feature_AH(rns, feats, atms, mps, idx, cnrs, inv, das, fn) for idx in range(len(mps))]
	elif sum(tps[2:4])>0:
		fns = fns + list(range(len(fn)+1,len(fn)+6))
		rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
		if tps[2]==1:
			vt = [feature_AA(rns, feats, atms, mps, idx, cnrs, inv, das, fn) for idx in range(len(mps))]
		else:
			vt = [feature_PP(rns,feats, mmffp, mps, idx, cnrs, gap, inv, das, fn) for idx in range(len(mpd))]
	elif tps[4]==1:
		fns = fns + list(range(len(fn)+1,len(fn)+10))
		rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
		vt = [feature_AP(rns, feats, atms, mmffp, gap,mps, idx, cnrs, inv, das, fn) for idx in range(len(mpd))]			
	else:
		rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
		vt = [pairLr_dist(f, rns, cnrs, inv, das[i], fn) for i, f in enumerate(feats)]
	vdf = pd.DataFrame(vt)
	vdf.index = mn
	return vdf



sym2feat2 = {'C':7,'O':10,'N':10, 'S':10}
def feature_dsts(row, feats, atms, mps, idx, col, inv, das, fn):
	# basic charac. extra atom to atom charac. 7 C hydrophobic / 10 O  N  S hydrophilic / 11 etc
	ft_idx = unique(list(itertools.chain(*[f.GetAtomIds() for f in feats[idx]])))
	uu_idx = [i for i in range(len(mps[idx])) if i not in ft_idx]
	use_pos = np.array([ft.GetPos() for ft in feats[idx]])
	use_fam = [fn.index(ft.GetFamily()) if not ft.GetAtomIds() in das[idx] else len(fn) for ft in feats[idx]]
	if (len(uu_idx) != 0) & (len(use_pos) !=0):
		use_idx = np.array(itemgetter(*uu_idx)(mps[idx]))
		if len(uu_idx)>1:
			use_atms = [sym2feat2[at.GetSymbol()] if at.GetSymbol() in ['C','O','N','S'] else len(fn)+2 for at in itemgetter(*uu_idx)(atms[idx])]
			poses = np.concatenate((use_pos,use_idx))
		else:
			atsym = itemgetter(*uu_idx)(atms[idx]).GetSymbol()
			poses = np.concatenate((use_pos,np.array([use_idx])))
			if atsym in ['C','O','N','S']:
				use_atms = [sym2feat2[atsym]]
			else:
				use_atms = [len(fn)+2]
		use_types =  use_fam + use_atms
	else:
		poses=use_pos
		use_types = use_fam
	ft=itertools.combinations(range(len(poses)), 2)
	lt=[0]*len(row)*len(col)
	dst = [np.linalg.norm(poses[f[0]]-poses[f[1]]) for f in ft]
	return dst


gap = 0.6
def dsts_from_cid(info):
	cids, sdf, tps = info
	factory = ChemicalFeatures.BuildFeatureFactory('{}/etc/BaseFeatures_h.fdef'.format(DIR_AP))
	fn = factory.GetFeatureFamilies()	
	fns = [i for i, f in enumerate(fn) if f in use_feats]
	mm = Chem.SDMolSupplier(sdf, removeHs=True) 
	mm = [m for m in mm if m is not None]
	mm = [m for m in mm if m.GetPropsAsDict()['PUBCHEM_COMPOUND_CID'] in cids]
	mn =  [m.GetPropsAsDict()['PUBCHEM_COMPOUND_CID'] for m in mm]
	feats = [factory.GetFeaturesForMol(m) for m in mm]
	if sum(tps[1:5])>0:
		mps = [m.GetConformer().GetPositions() for m in mm]	
		if tps[3]!=1:
			atms = [m.GetAtoms() for m in mm]
		if sum(tps[3:5])>0:
			mpd = [m.GetPropsAsDict() for m in mm]
			mmffp = [m['PUBCHEM_MMFF94_PARTIAL_CHARGES'] if 'PUBCHEM_MMFF94_PARTIAL_CHARGES' in m.keys() else '' for m in mpd]
			mmffp = [m.split('\n') if (m!=0) & (m!='') else '' for m in mmffp]		
	# featur n selection
	if tps[0]==1:
		feats = [filter_feats_renew2(feat) for feat in feats]
		fns = fns + [len(fn)]
		das = [[f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')] for feat in feats]
		das = [unique([x for x in da if da.count(x) > 1]) for da in das]	
		feats = [filter_feats_both(feat) for feat in feats]
	if tps[1]==1:
		fns = fns + list(range(len(fn)+1,len(fn)+3))
		rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
		dsts = [feature_dsts(rns, feats, atms, mps, idx, cnrs, inv, das, fn) for idx in range(len(mps))]
		dsts = list(itertools.chain(*dsts))
	return tuple(dsts)

aa = dsts_from_cid(infos[-10])

dstss = mclapply(infos, dsts_from_cid, ncpu = 220)

with open(dtip+'/dstss.pkl','wb') as f:
    pickle.dump([dstss, dstsss],f)


dstss_all = list(itertools.chain(*dstsss))
# density plot 
import seaborn as sns
fig = plt.figure(figsize = (5,5))
ax = fig.add_subplot(1,1,1)
# cols= ['#F45050','#FF8400', '#3C486B','#87CBB9']
cols= ['#C3E2C2','#CD8D7A', '#3081D0','#DBCC95']

sns.kdeplot(dstss_all, shade=True, color = '#12372A', bw=0.4)
plt.tick_params(
    axis='y',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False)

plt.xlim([-1.5,20])
ax.set_xlabel('Distance (Å)')
# ax.set_xticks([0]+[round(a,1) for a in cnrs]+[20])
plt.savefig(f'{dtip}/figures/pair_dst_dense_0.4.png', dpi=300, bbox_inches = 'tight', transparent = True)


for i in range(2):
	plt.subplot(2,1,i+1)
	for ibb in range(4):
		sns.kdeplot(bbb[ibb][i], shade=True, color = cols[ibb], bw=0.1)
	handles = [Patch(facecolor=c, edgecolor='black', label=nms[i]) for i,c in enumerate(cols)]
	plt.legend(handles=handles, fontsize = '8.5', loc = 'upper right', ncol=1, framealpha=0)
	plt.tick_params(
	    axis='y',          # changes apply to the x-axis
	    which='both',      # both major and minor ticks are affected
	    bottom=False,      # ticks along the bottom edge are off
	    top=False,         # ticks along the top edge are off
	    labelbottom=False)
	plt.ylabel(['BindingDB','ChEMBL'][i])

# plt.legend(handles=handles, fontsize = '8.5', loc = 'lower center', ncol=4, framealpha=0)

plt.savefig(f'{dtip}/figures/single_multi_dens.png', dpi=300, bbox_inches = 'tight', transparent = True)





infos = [[f'{natdir}/out.sdf',tps]]

vt = []
mnss = []
for i,idx in enumerate(range(len(mps))):
	try:
		vt.append(feature_AH(rns, feats, atms, mps, idx, cnrs, inv, das, fn))
		mnss.append(mns[i])
	except:
		print(mns[i])


vt2 = []
mnss2 = []
feats2 = [filter_feats_both(feat) for feat in feats]
for i,idx in enumerate(range(len(mps))):
	try:
		vt2.append(feature_AH(rns, feats2, atms, mps, idx, cnrs, inv, das, fn))
		mnss2.append(mns[i])
	except:
		print(mns[i])


fps= []
mnssf = []
for i,m in enumerate(mm):
	try:
		if m is not None:
			fps.append(list(Chem.rdMolDescriptors.GetMorganFingerprintAsBitVect(m, radius)))
			mnssf.append(mns[i])
	except:
		print(mns[i])


external_lts = pd.DataFrame(vt)
external_bts = pd.DataFrame(vt2)
external_fps = pd.DataFrame(fps)
external_lts.index = mns
external_bts.index = mns
external_fps.index = mns
external_lts.to_pickle(f'{natdir}/AP1_no_both_external_smiles.pkl')
external_bts.to_pickle(f'{natdir}/AP1_both_external_smiles.pkl')
external_fps.to_pickle(f'{natdir}/fps_external_smiles.pkl')




gap = 0.6
def AP_from_file(info):
	sdf, tps = info
	factory = ChemicalFeatures.BuildFeatureFactory('{}/etc/BaseFeatures_h.fdef'.format(DIR_AP))
	fn = factory.GetFeatureFamilies()	
	fns = [i for i, f in enumerate(fn) if f in use_feats]
	mm = Chem.SDMolSupplier(sdf, removeHs=True) 
	mm = [m for m in mm if m is not None]
	mns = [m.GetProp('_Name') for m in mm]
	feats = [factory.GetFeaturesForMol(m) for m in mm]
	if sum(tps[1:5])>0:
		mps = [m.GetConformer().GetPositions() for m in mm]	
		if tps[3]!=1:
			atms = [m.GetAtoms() for m in mm]
		if sum(tps[3:5])>0:
			mpd = [m.GetPropsAsDict() for m in mm]
			mmffp = [m['PUBCHEM_MMFF94_PARTIAL_CHARGES'] if 'PUBCHEM_MMFF94_PARTIAL_CHARGES' in m.keys() else '' for m in mpd]
			mmffp = [m.split('\n') if (m!=0) & (m!='') else '' for m in mmffp]		
	# featur n selection
	if tps[0]==1:
		feats = [filter_feats_renew2(feat) for feat in feats]
		fns = fns + [len(fn)]
		das = [[f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')] for feat in feats]
		das = [unique([x for x in da if da.count(x) > 1]) for da in das]	
		# feats = [filter_feats_both(feat) for feat in feats]
	if tps[1]==1:
		fns = fns + list(range(len(fn)+1,len(fn)+3))
		rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
		vt = [feature_AH(rns, feats, atms, mps, idx, cnrs, inv, das, fn) for idx in range(len(mps))]
	elif sum(tps[2:4])>0:
		fns = fns + list(range(len(fn)+1,len(fn)+6))
		rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
		if tps[2]==1:
			vt = [feature_AA(rns, feats, atms, mps, idx, cnrs, inv, das, fn) for idx in range(len(mps))]
		else:
			vt = [feature_PP(rns,feats, mmffp, mps, idx, cnrs, gap, inv, das, fn) for idx in range(len(mpd))]
	elif tps[4]==1:
		fns = fns + list(range(len(fn)+1,len(fn)+10))
		rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
		vt = [feature_AP(rns, feats, atms, mmffp, gap,mps, idx, cnrs, inv, das, fn) for idx in range(len(mpd))]			
	else:
		rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
		vt = [pairLr_dist(f, rns, cnrs, inv, das[i], fn) for i, f in enumerate(feats)]
	vdf = pd.DataFrame(vt)
	vdf.index = mns
	return vdf


gap = 0.6
def AP_from_dud(info):
	cids, sdf, tps = info
	factory = ChemicalFeatures.BuildFeatureFactory('{}/etc/BaseFeatures_h.fdef'.format(DIR_AP))
	fn = factory.GetFeatureFamilies()	
	fns = [i for i, f in enumerate(fn) if f in use_feats]
	mm = Chem.SDMolSupplier(sdf, removeHs=True) 
	mm = [m for m in mm if m is not None]
	mm = [m for m in mm if m.GetProp('_Name') in cids]
	mn =  [m.GetProp('_Name') for m in mm]
	feats = [factory.GetFeaturesForMol(m) for m in mm]
	if sum(tps[1:5])>0:
		mps = [m.GetConformer().GetPositions() for m in mm]	
		if tps[3]!=1:
			atms = [m.GetAtoms() for m in mm]
		if sum(tps[3:5])>0:
			mpd = [m.GetPropsAsDict() for m in mm]
			mmffp = [m['PUBCHEM_MMFF94_PARTIAL_CHARGES'] if 'PUBCHEM_MMFF94_PARTIAL_CHARGES' in m.keys() else '' for m in mpd]
			mmffp = [m.split('\n') if (m!=0) & (m!='') else '' for m in mmffp]		
	# featur n selection
	if tps[0]==1:
		feats = [filter_feats_renew2(feat) for feat in feats]
		fns = fns + [len(fn)]
		das = [[f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')] for feat in feats]
		das = [unique([x for x in da if da.count(x) > 1]) for da in das]	
		# feats = [filter_feats_both(feat) for feat in feats]
	if tps[1]==1:
		fns = fns + list(range(len(fn)+1,len(fn)+3))
		rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
		vt = [feature_AH(rns, feats, atms, mps, idx, cnrs, inv, das, fn) for idx in range(len(mps))]
	elif sum(tps[2:4])>0:
		fns = fns + list(range(len(fn)+1,len(fn)+6))
		rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
		if tps[2]==1:
			vt = [feature_AA(rns, feats, atms, mps, idx, cnrs, inv, das, fn) for idx in range(len(mps))]
		else:
			vt = [feature_PP(rns,feats, mmffp, mps, idx, cnrs, gap, inv, das, fn) for idx in range(len(mpd))]
	elif tps[4]==1:
		fns = fns + list(range(len(fn)+1,len(fn)+10))
		rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
		vt = [feature_AP(rns, feats, atms, mmffp, gap,mps, idx, cnrs, inv, das, fn) for idx in range(len(mpd))]			
	else:
		rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
		vt = [pairLr_dist(f, rns, cnrs, inv, das[i], fn) for i, f in enumerate(feats)]
	vdf = pd.DataFrame(vt)
	vdf.index = mn
	vdf['target'] = sdf.split('_')[-2]
	return vdf



def sdf2scf(info):
	cids, sdf, tps = info
	mm = Chem.SDMolSupplier(sdf, removeHs=True) 
	mm = [m for m in mm if m is not None]
	mm = [m for m in mm if m.GetPropsAsDict()['PUBCHEM_COMPOUND_CID'] in cids]
	mn =  [m.GetPropsAsDict()['PUBCHEM_COMPOUND_CID'] for m in mm]
	scf = [Chem.MolToSmiles(GetScaffoldForMol(m)) for m in mm]
	scfs = pd.DataFrame(scf)
	scfs.index = mn
	return scfs



def sdf2smi(info):
	cids, sdf, tps = info
	mm = Chem.SDMolSupplier(sdf, removeHs=True) 
	mm = [m for m in mm if m is not None]
	mm = [m for m in mm if m.GetPropsAsDict()['PUBCHEM_COMPOUND_CID'] in cids]
	mn =  [m.GetPropsAsDict()['PUBCHEM_COMPOUND_CID'] for m in mm]
	scf = [Chem.MolToSmiles(m) for m in mm]
	scfs = pd.DataFrame(scf)
	scfs.index = mn
	return scfs



vdf2 = mclapply(infos, sdf2smi, ncpu=100)



cidss = [[] for i in range(6358)]
for cid in cid_list:
	if cid < 158950000:
		cidss[int(cid // 25000)].append(cid)

sdfs = ['{}/{}_{}.sdf'.format(str_dir,str(i * 25000 +1).zfill(8), str((i + 1) * 25000).zfill(8)) for i,cd in enumerate(cidss)
 if len(cd)!=0]
cds = [cd for cd in cidss if len(cd)!=0]

infos = [[unique(dude[(dude.target==subdirs[i]) & (dude.label==1)].id), sdfs[i], tps] for i in range(len(subdirs))]

sdfs_acts= [DIR_AP+'/DB/DUDE/all/'+sub+'/actives_final.sdf' for sub in subdirs]
sdfs_decoys= [DIR_AP+'/DB/DUDE/all/'+sub+'/decoys_final.sdf' for sub in subdirs]

# dude=pd.read_csv(dtip +'/DUDE/dud_dataset.csv')

# dude['id'] = ['ZIN'+cd if len(cd)==9 else cd for cd in list(dude.id)]

# infos = [[unique(dude[(dude.target==subdirs[i]) & (dude.label==0)].id), sdfs[i], tps] for i in range(len(subdirs))]
infos = [[sdfs[i], tps] for i in range(len(subdirs))]



start_time = time.time()
ltis = mclapply(infos,sdf2fp, ncpu= 50) # 100 / 10 855s   / 700 / 20 --> 57m | 3466s
print("---{}s seconds---".format(time.time()-start_time))

pd.concat(ltis).to_pickle(dtip+'/pdbval/gpdb_cids_fp.pkl')

pd.concat(ltis).to_pickle(dtip+'/pdbval/sim_FP.pkl')



start_time = time.time()
ltis = mclapply(infos,AP_from_dud, ncpu= 102) # 100 / 10 855s   / 700 / 20 --> 57m | 3466s
print("---{}s seconds---".format(time.time()-start_time))

lltss = pd.concat(ltis)
lltss.to_pickle(dtip+'/DUDE/AP0_active_1.pkl')


lltss = pd.concat(ltis)
lltss.to_pickle(dtip+'/DUDE/AP0_decoy_1.pkl')


def sdf2scf(info):
	cids, sdf, tps = info
	mm = Chem.SDMolSupplier(sdf, removeHs=False) 
	mm = [m for m in mm if m is not None]
	mm = [m for m in mm if m.GetPropsAsDict()['PUBCHEM_COMPOUND_CID'] in cids]
	mn = [m.GetPropsAsDict()['PUBCHEM_COMPOUND_CID'] for m in mm]
	scf = [Chem.MolToSmiles(GetScaffoldForMol(m)) for m in mm]
	scfs = pd.DataFrame(scf)
	scfs.index = mn
	return scfs



def sdf2smi(infos):
	cids, sdf, tps = infos
	mm = Chem.SDMolSupplier(sdf, removeHs=True) 
	mm = [m for m in mm if m is not None]
	mm = [m for m in mm if m.GetPropsAsDict()['PUBCHEM_COMPOUND_CID'] in cids]
	mn = [m.GetPropsAsDict()['PUBCHEM_COMPOUND_CID'] for m in mm]
	smi = [Chem.MolToSmiles(m) for m in mm]
	df = pd.DataFrame({'cid':mn, 'smiles':smi})
	return df


ltis = mclapply(infos,sdf2smi, ncpu= 150) # 100 / 10 855s   / 700 / 20 --> 57m | 3466s

factory = ChemicalFeatures.BuildFeatureFactory('{}/etc/BaseFeatures_h.fdef'.format(DIR_AP))
fn = factory.GetFeatureFamilies()	
def process_apm(m):
	try:
		fns = [i for i, f in enumerate(fn) if f in use_feats]
		feats = factory.GetFeaturesForMol(m)
		mps = [m.GetConformer().GetPositions()]
		atms = [m.GetAtoms()]
		feats = [filter_feats_renew(feats)]
		fns = fns + [len(fn)]
		das = [[f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')] for feat in feats]
		das = [unique([x for x in das if das.count(x) > 1])]
		fns = fns + list(range(len(fn)+1,len(fn)+3))
		rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
		idx = 0
		vt = feature_AH(rns, feats, atms, mps, idx, cnrs, inv, das, fn)
	except:
		vt = [0]*len(rns)*len(cnrs)
	return tuple(vt)


start_time = time.time()
aa = mclapply(mm, process_apm, ncpu = 90)
print("---{}s seconds---".format(time.time()-start_time))

gap = 0.6
def AP_from_file_new(info):
	try:
		pps, tps = info
		factory = ChemicalFeatures.BuildFeatureFactory('{}/etc/BaseFeatures_h.fdef'.format(DIR_AP))
		fn = factory.GetFeatureFamilies()	
		fns = [i for i, f in enumerate(fn) if f in use_feats]
		mm = Chem.SDMolSupplier('{}/DB/refined-set/{}/{}_ligand.sdf'.format(DIR_AP,pps,pps), removeHs=True)
		if mm is None:
			mm = Chem.rdmolfiles.MolFromMol2File('{}/DB/refined-set/{}/{}_ligand.mol2'.format(DIR_AP,pps,pps), removeHs=True)
		# mm = Chem.SDMolSupplier(sdf, removeHs=True) 
		mm = [m for m in mm if m is not None]
		feats = [factory.GetFeaturesForMol(m) for m in mm]
		if sum(tps[1:5])>0:
			mps = [m.GetConformer().GetPositions() for m in mm]	
			if tps[3]!=1:
				atms = [m.GetAtoms() for m in mm]
			if sum(tps[3:5])>0:
				mpd = [m.GetPropsAsDict() for m in mm]
				mmffp = [m['PUBCHEM_MMFF94_PARTIAL_CHARGES'] if 'PUBCHEM_MMFF94_PARTIAL_CHARGES' in m.keys() else '' for m in mpd]
				mmffp = [m.split('\n') if (m!=0) & (m!='') else '' for m in mmffp]		
		# featur n selection
		if tps[0]==1:
			feats = [filter_feats_renew(feat) for feat in feats]
			fns = fns + [len(fn)]
			das = [[f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')] for feat in feats]
			das = [unique([x for x in da if da.count(x) > 1]) for da in das]		
		if tps[1]==1:
			fns = fns + list(range(len(fn)+1,len(fn)+3))
			rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
			vt = [feature_AH(rns, feats, atms, mps, idx, cnrs, inv, das, fn) for idx in range(len(mps))]
		elif sum(tps[2:4])>0:
			fns = fns + list(range(len(fn)+1,len(fn)+6))
			rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
			if tps[2]==1:
				vt = [feature_AA(rns, feats, atms, mps, idx, cnrs, inv, das, fn) for idx in range(len(mps))]
			else:
				vt = [feature_PP(rns,feats, mmffp, mps, idx, cnrs, gap, inv, das, fn) for idx in range(len(mpd))]
		elif tps[4]==1:
			fns = fns + list(range(len(fn)+1,len(fn)+10))
			rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
			vt = [feature_AP(rns, feats, atms, mmffp, gap,mps, idx, cnrs, inv, das, fn) for idx in range(len(mpd))]			
		else:
			rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
			vt = [pairLr_dist(f, rns, cnrs, inv, das[i], fn) for i, f in enumerate(feats)]
		vdf = pd.DataFrame(vt)
		vdf.index = [pps]
	except:vdf =pd.DataFrame()
	return vdf



gap = 0.6
def AP_from_file(info):
	try:
		pps, tps = info
		factory = ChemicalFeatures.BuildFeatureFactory('{}/etc/BaseFeatures_h.fdef'.format(DIR_AP))
		fn = factory.GetFeatureFamilies()	
		fns = [i for i, f in enumerate(fn) if f in use_feats]
		mm = Chem.SDMolSupplier('{}/DB/refined-set/{}/{}_ligand.sdf'.format(DIR_AP,pps,pps), removeHs=True)
		if mm is None:
			mm = Chem.rdmolfiles.MolFromMol2File('{}/DB/refined-set/{}/{}_ligand.mol2'.format(DIR_AP,pps,pps), removeHs=True)
		# mm = Chem.SDMolSupplier(sdf, removeHs=True) 
		mm = [m for m in mm if m is not None]
		feats = [factory.GetFeaturesForMol(m) for m in mm]
		if sum(tps[1:5])>0:
			mps = [m.GetConformer().GetPositions() for m in mm]	
			if tps[3]!=1:
				atms = [m.GetAtoms() for m in mm]
			if sum(tps[3:5])>0:
				mpd = [m.GetPropsAsDict() for m in mm]
				mmffp = [m['PUBCHEM_MMFF94_PARTIAL_CHARGES'] if 'PUBCHEM_MMFF94_PARTIAL_CHARGES' in m.keys() else '' for m in mpd]
				mmffp = [m.split('\n') if (m!=0) & (m!='') else '' for m in mmffp]		
		# featur n selection
		if tps[0]==1:
			feats = [filter_feats_renew(feat) for feat in feats]
			fns = fns + [len(fn)]
			das = [[f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')] for feat in feats]
			das = [unique([x for x in da if da.count(x) > 1]) for da in das]		
		if tps[1]==1:
			fns = fns + list(range(len(fn)+1,len(fn)+3))
			rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
			vt = [feature_AH(rns, feats, atms, mps, idx, cnrs, inv, das, fn) for idx in range(len(mps))]
		elif sum(tps[2:4])>0:
			fns = fns + list(range(len(fn)+1,len(fn)+6))
			rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
			if tps[2]==1:
				vt = [feature_AA(rns, feats, atms, mps, idx, cnrs, inv, das, fn) for idx in range(len(mps))]
			else:
				vt = [feature_PP(rns,feats, mmffp, mps, idx, cnrs, gap, inv, das, fn) for idx in range(len(mpd))]
		elif tps[4]==1:
			fns = fns + list(range(len(fn)+1,len(fn)+10))
			rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
			vt = [feature_AP(rns, feats, atms, mmffp, gap,mps, idx, cnrs, inv, das, fn) for idx in range(len(mpd))]			
		else:
			rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
			vt = [pairLr_dist(f, rns, cnrs, inv, das[i], fn) for i, f in enumerate(feats)]
		vdf = pd.DataFrame(vt)
		vdf.index = [pps]
	except:vdf =pd.DataFrame()
	return vdf


ltsi = mclapply([[pps, tps] for pps in cmpss], AP_from_file, ncpu = 100)

gap = 0.6
def AP_from_file(info):
	try:
		pps, tps = info
		factory = ChemicalFeatures.BuildFeatureFactory('{}/etc/BaseFeatures_h.fdef'.format(DIR_AP))
		fn = factory.GetFeatureFamilies()	
		fns = [i for i, f in enumerate(fn) if f in use_feats]
		mm = Chem.SDMolSupplier(pps, removeHs=True)
		mm = [m for m in mm if m is not None]
		feats = [factory.GetFeaturesForMol(m) for m in mm]
		if sum(tps[1:5])>0:
			mps = [m.GetConformer().GetPositions() for m in mm]	
			if tps[3]!=1:
				atms = [m.GetAtoms() for m in mm]
			if sum(tps[3:5])>0:
				mpd = [m.GetPropsAsDict() for m in mm]
				mmffp = [m['PUBCHEM_MMFF94_PARTIAL_CHARGES'] if 'PUBCHEM_MMFF94_PARTIAL_CHARGES' in m.keys() else '' for m in mpd]
				mmffp = [m.split('\n') if (m!=0) & (m!='') else '' for m in mmffp]		
		# featur n selection
		if tps[0]==1:
			feats = [filter_feats_renew(feat) for feat in feats]
			fns = fns + [len(fn)]
			das = [[f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')] for feat in feats]
			das = [unique([x for x in da if da.count(x) > 1]) for da in das]		
		if tps[1]==1:
			fns = fns + list(range(len(fn)+1,len(fn)+3))
			rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
			vt = [feature_AH(rns, feats, atms, mps, idx, cnrs, inv, das, fn) for idx in range(len(mps))]
		elif sum(tps[2:4])>0:
			fns = fns + list(range(len(fn)+1,len(fn)+6))
			rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
			if tps[2]==1:
				vt = [feature_AA(rns, feats, atms, mps, idx, cnrs, inv, das, fn) for idx in range(len(mps))]
			else:
				vt = [feature_PP(rns,feats, mmffp, mps, idx, cnrs, gap, inv, das, fn) for idx in range(len(mpd))]
		elif tps[4]==1:
			fns = fns + list(range(len(fn)+1,len(fn)+10))
			rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
			vt = [feature_AP(rns, feats, atms, mmffp, gap,mps, idx, cnrs, inv, das, fn) for idx in range(len(mpd))]			
		else:
			rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
			vt = [pairLr_dist(f, rns, cnrs, inv, das[i], fn) for i, f in enumerate(feats)]
		vdf = pd.DataFrame(vt)
		vdf.index = [pps.split('/')[-1].split('.')[0]]
	except:vdf =pd.DataFrame()
	return vdf


gap = 0.6
def AP_from_file2(info):
	try:
		pps, tps = info
		factory = ChemicalFeatures.BuildFeatureFactory('{}/etc/BaseFeatures_h.fdef'.format(DIR_AP))
		fn = factory.GetFeatureFamilies()	
		fns = [i for i, f in enumerate(fn) if f in use_feats]
		# mm = Chem.SDMolSupplier('{}/DB/refined-set/{}/{}_ligand.sdf'.format(DIR_AP,pps,pps), removeHs=True)
		# if mm is None:
		mm = Chem.rdmolfiles.MolFromMol2File('{}/DB/refined-set/{}/{}_ligand.mol2'.format(DIR_AP,pps,pps), removeHs=True)
		# mm = Chem.SDMolSupplier(sdf, removeHs=True) 
		# mm = [m for m in mm if m is not None]
		feats = [factory.GetFeaturesForMol(mm)]
		if sum(tps[1:5])>0:
			mps = [mm.GetConformer().GetPositions()]	
			if tps[3]!=1:
				atms = [mm.GetAtoms()]
			if sum(tps[3:5])>0:
				mpd = [mm.GetPropsAsDict()]
				mmffp = [m['PUBCHEM_MMFF94_PARTIAL_CHARGES'] if 'PUBCHEM_MMFF94_PARTIAL_CHARGES' in m.keys() else '' for m in mpd]
				mmffp = [m.split('\n') if (m!=0) & (m!='') else '' for m in mmffp]		
		# featur n selection
		if tps[0]==1:
			feats = [filter_feats_renew2(feat) for feat in feats]
			fns = fns + [len(fn)]
			das = [[f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')] for feat in feats]
			das = [unique([x for x in da if da.count(x) > 1]) for da in das]	
			feats = [filter_feats_both(feat) for feat in feats]	
		if tps[1]==1:
			fns = fns + list(range(len(fn)+1,len(fn)+3))
			rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
			vt = [feature_AH(rns, feats, atms, mps, idx, cnrs, inv, das, fn) for idx in range(len(mps))]
		elif sum(tps[2:4])>0:
			fns = fns + list(range(len(fn)+1,len(fn)+6))
			rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
			if tps[2]==1:
				vt = [feature_AA(rns, feats, atms, mps, idx, cnrs, inv, das, fn) for idx in range(len(mps))]
			else:
				vt = [feature_PP(rns,feats, mmffp, mps, idx, cnrs, gap, inv, das, fn) for idx in range(len(mpd))]
		elif tps[4]==1:
			fns = fns + list(range(len(fn)+1,len(fn)+10))
			rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
			vt = [feature_AP(rns, feats, atms, mmffp, gap,mps, idx, cnrs, inv, das, fn) for idx in range(len(mpd))]			
		else:
			rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
			vt = [pairLr_dist(f, rns, cnrs, inv, das[i], fn) for i, f in enumerate(feats)]
		vdf = pd.DataFrame(vt)
		vdf.index = [pps]
	except:vdf =pd.DataFrame()
	return vdf


gap = 0.6
def AP_from_file(info):
	try:
		sdf, tps = info
		factory = ChemicalFeatures.BuildFeatureFactory('{}/etc/BaseFeatures_h.fdef'.format(DIR_AP))
		fn = factory.GetFeatureFamilies()	
		fns = [i for i, f in enumerate(fn) if f in use_feats]
		# mm = Chem.SDMolSupplier('{}/DB/refined-set/{}/{}_ligand.sdf'.format(DIR_AP,pps,pps), removeHs=True)
		# if mm is None:
		# 	mm = Chem.rdmolfiles.MolFromMol2File('{}/DB/refined-set/{}/{}_ligand.mol2'.format(DIR_AP,pps,pps), removeHs=True)
		mm = Chem.SDMolSupplier(sdf, removeHs=True) 
		mm = [m for m in mm if m is not None]
		feats = [factory.GetFeaturesForMol(m) for m in mm]
		if sum(tps[1:5])>0:
			mps = [m.GetConformer().GetPositions() for m in mm]	
			if tps[3]!=1:
				atms = [m.GetAtoms() for m in mm]
			if sum(tps[3:5])>0:
				mpd = [m.GetPropsAsDict() for m in mm]
				mmffp = [m['PUBCHEM_MMFF94_PARTIAL_CHARGES'] if 'PUBCHEM_MMFF94_PARTIAL_CHARGES' in m.keys() else '' for m in mpd]
				mmffp = [m.split('\n') if (m!=0) & (m!='') else '' for m in mmffp]		
		# featur n selection
		if tps[0]==1:
			feats = [filter_feats_renew(feat) for feat in feats]
			fns = fns + [len(fn)]
			das = [[f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')] for feat in feats]
			das = [unique([x for x in da if da.count(x) > 1]) for da in das]		
		if tps[1]==1:
			fns = fns + list(range(len(fn)+1,len(fn)+3))
			rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
			vt = [feature_AH(rns, feats, atms, mps, idx, cnrs, inv, das, fn) for idx in range(len(mps))]
		elif sum(tps[2:4])>0:
			fns = fns + list(range(len(fn)+1,len(fn)+6))
			rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
			if tps[2]==1:
				vt = [feature_AA(rns, feats, atms, mps, idx, cnrs, inv, das, fn) for idx in range(len(mps))]
			else:
				vt = [feature_PP(rns,feats, mmffp, mps, idx, cnrs, gap, inv, das, fn) for idx in range(len(mpd))]
		elif tps[4]==1:
			fns = fns + list(range(len(fn)+1,len(fn)+10))
			rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
			vt = [feature_AP(rns, feats, atms, mmffp, gap,mps, idx, cnrs, inv, das, fn) for idx in range(len(mpd))]			
		else:
			rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
			vt = [pairLr_dist(f, rns, cnrs, inv, das[i], fn) for i, f in enumerate(feats)]
		vdf = pd.DataFrame(vt)
		# vdf.index = [pps]
		vdf.index = [sdf.split('/')[-1].split('_')[0]]
	except:vdf =pd.DataFrame()
	return vdf

gap = 0.6
def AP_from_file(info):
	try:
		sdf, tps = info
		factory = ChemicalFeatures.BuildFeatureFactory('{}/etc/BaseFeatures_h.fdef'.format(DIR_AP))
		fn = factory.GetFeatureFamilies()	
		fns = [i for i, f in enumerate(fn) if f in use_feats]
		# mm = Chem.SDMolSupplier('{}/DB/refined-set/{}/{}_ligand.sdf'.format(DIR_AP,pps,pps), removeHs=True)
		# if mm is None:
		# 	mm = Chem.rdmolfiles.MolFromMol2File('{}/DB/refined-set/{}/{}_ligand.mol2'.format(DIR_AP,pps,pps), removeHs=True)
		mm = Chem.SDMolSupplier(sdf, removeHs=True) 
		mm = [m for m in mm if m is not None]
		feats = [factory.GetFeaturesForMol(m) for m in mm]
		if sum(tps[1:5])>0:
			mps = [m.GetConformer().GetPositions() for m in mm]	
			if tps[3]!=1:
				atms = [m.GetAtoms() for m in mm]
			if sum(tps[3:5])>0:
				mpd = [m.GetPropsAsDict() for m in mm]
				mmffp = [m['PUBCHEM_MMFF94_PARTIAL_CHARGES'] if 'PUBCHEM_MMFF94_PARTIAL_CHARGES' in m.keys() else '' for m in mpd]
				mmffp = [m.split('\n') if (m!=0) & (m!='') else '' for m in mmffp]		
		# featur n selection
		if tps[0]==1:
			feats = [filter_feats_renew2(feat) for feat in feats]
			fns = fns + [len(fn)]
			das = [[f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')] for feat in feats]
			das = [unique([x for x in da if da.count(x) > 1]) for da in das]	
			feats = [filter_feats_both(feat) for feat in feats]
		if tps[1]==1:
			fns = fns + list(range(len(fn)+1,len(fn)+3))
			rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
			vt = [feature_AH(rns, feats, atms, mps, idx, cnrs, inv, das, fn) for idx in range(len(mps))]
		elif sum(tps[2:4])>0:
			fns = fns + list(range(len(fn)+1,len(fn)+6))
			rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
			if tps[2]==1:
				vt = [feature_AA(rns, feats, atms, mps, idx, cnrs, inv, das, fn) for idx in range(len(mps))]
			else:
				vt = [feature_PP(rns,feats, mmffp, mps, idx, cnrs, gap, inv, das, fn) for idx in range(len(mpd))]
		elif tps[4]==1:
			fns = fns + list(range(len(fn)+1,len(fn)+10))
			rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
			vt = [feature_AP(rns, feats, atms, mmffp, gap,mps, idx, cnrs, inv, das, fn) for idx in range(len(mpd))]			
		else:
			rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
			vt = [pairLr_dist(f, rns, cnrs, inv, das[i], fn) for i, f in enumerate(feats)]
		vdf = pd.DataFrame(vt)
		# vdf.index = [pps]
		vdf.index = [sdf.split('/')[-1].split('_')[0]]
	except:vdf =pd.DataFrame()
	return vdf

import glob
idd =  glob.glob(dtip+'/pdbval/sdf/*model*')



start_time = time.time()
ltis2 = mclapply([[md,[1,1,0,0,0]] for md in ppp],AP_from_file2, ncpu= 100) # 100 / 10 855s   / 700 / 20 --> 57m | 3466s
print("---{}s seconds---".format(time.time()-start_time))

start_time = time.time()
ltis = mclapply([[md,[1,1,0,0,0]] for md in idd],AP_from_file, ncpu= 100) # 100 / 10 855s   / 700 / 20 --> 57m | 3466s
print("---{}s seconds---".format(time.time()-start_time))


ppdf2 = pd.concat(ltis2)
ppdf2 = ppdf2[ppdf2.T.sum()!=0]

ppdf.to_pickle('/spstorage/USERS/gina/Project/AP/DB/InterpretableDTIP/pdbval/AP1_renew_pdb.pkl')
ppdf2.to_pickle('/spstorage/USERS/gina/Project/AP/DB/InterpretableDTIP/pdbval/AP1_renew2_both_pdb.pkl')


iidf = pd.concat(ltis)
iidf = iidf[iidf.T.sum()!=0]
iidf.to_pickle('/spstorage/USERS/gina/Project/AP/DB/InterpretableDTIP/pdbval/AP1_renew_ideal.pkl')
iidf.to_pickle('/spstorage/USERS/gina/Project/AP/DB/InterpretableDTIP/pdbval/AP1_renew2_no_both_ideal.pkl')
iidf.to_pickle('/spstorage/USERS/gina/Project/AP/DB/InterpretableDTIP/pdbval/AP1_renew2_both_ideal.pkl')



mmdf = pd.concat(ltis)
mmdf = mmdf[mmdf.T.sum()!=0]
mmdf.to_pickle('/spstorage/USERS/gina/Project/AP/DB/InterpretableDTIP/pdbval/AP1_renew2_both_model.pkl')


scfs6 = mclapply(infos, sdf2scf, ncpu = 50)
dd = pd.concat(scfs6)

cid_info = pd.read_csv(chemm+'/cid_info.csv')
cid_info = pd.read_csv(dtip+'/cid_list.csv')
subd = dd[dd.iloc[:,0]=='']
smis = [cid_info[cid_info.cid==cd]['isosmiles'].item() for cd in list(subd.index)]
smis = [list(set(allsets[allsets.cid==cd]['smiles']))[0] for cd in list(subd.index)]

subd[0] = smis
dd = pd.concat([dd[dd.iloc[:,0]!=''], subd])
c2s = pd.merge(dd, cid_info.loc[:,['cid','mw','isosmiles','rotbonds']], left_index = True, right_on = 'cid')
c2s = pd.merge(dd, allsets,  left_index = True, right_on = 'cid')

c2s.columns = ['scaffold','cid', 'pdb_id', 'label', 'tag', 'type', 'smiles', 'sequence']
c2s.columns = ['scaffold','cid','mw','isosmiles','rotbonds']

c2s.to_csv(chemm + '/cid2scf.csv')
c2s.to_csv(dtip + '/data/cid2scf.csv')
c2s.to_csv(dtip+'/human/cid2scf.csv')
c2s = pd.read_csv(chemm + '/cid2scf.csv', index_col =0)
c2s = pd.read_csv(dtip + '/human/cid2scf.csv', index_col =0)

allsets = pd.read_csv(chemm+'/dfs.csv', index_col = 0)

dfs2scf = pd.merge(allsets, c2s)
subscfs = []
subscf = pd.DataFrame()
while len(subscf)<len(dfs2scf)*0.8:
	subscfs.append(random.choice(list(dfs2scf.scaffold)))
	subscf = dfs2scf[dfs2scf.scaffold.isin(subscfs)]



test_subscfs = []
subscf = pd.DataFrame()
while len(subscf)<len(dfs2scf)*0.1:
	test_subscfs.append(random.choice(list(set(dfs2scf.scaffold) - set(subscfs))))
	subscf = dfs2scf[dfs2scf.scaffold.isin(test_subscfs)]

test_subscfs = test_subscfs[0:-5000]
dfs2scf['tag'] = ['train' if scf in subscfs else 'test' if scf in test_subscfs else 'dev' for scf in list(dfs2scf.scaffold)]
dfs2scf.to_csv(chemm+'/dataset_3D_scaffold.csv')
dfs2scf = dfs2scf.sample(frac= 1 ,random_state = 7)
dfs2scf.to_csv(dtip+'/human/human_dataset_3D_scaffold.csv')

gap = 0.6
def AP_from_dud(info):
	sdf, tps = info
	factory = ChemicalFeatures.BuildFeatureFactory('{}/etc/BaseFeatures_h.fdef'.format(DIR_AP))
	fn = factory.GetFeatureFamilies()	
	fns = [i for i, f in enumerate(fn) if f in use_feats]
	mm = Chem.SDMolSupplier(sdf, removeHs=True) 
	mm = [m for m in mm if m is not None]
	# mm = [m for m in mm if m.GetProp('_Name') in cids]
	mn =  [m.GetProp('_Name') for m in mm]
	feats = [factory.GetFeaturesForMol(m) for m in mm]
	if sum(tps[1:5])>0:
		mps = [m.GetConformer().GetPositions() for m in mm]	
		if tps[3]!=1:
			atms = [m.GetAtoms() for m in mm]
		if sum(tps[3:5])>0:
			mpd = [m.GetPropsAsDict() for m in mm]
			mmffp = [m['PUBCHEM_MMFF94_PARTIAL_CHARGES'] if 'PUBCHEM_MMFF94_PARTIAL_CHARGES' in m.keys() else '' for m in mpd]
			mmffp = [m.split('\n') if (m!=0) & (m!='') else '' for m in mmffp]		
	# featur n selection
	if tps[0]==1:
		feats = [filter_feats_renew(feat) for feat in feats]
		fns = fns + [len(fn)]
		das = [[f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')] for feat in feats]
		das = [unique([x for x in da if da.count(x) > 1]) for da in das]		
	if tps[1]==1:
		fns = fns + list(range(len(fn)+1,len(fn)+3))
		rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
		vt = [feature_AH(rns, feats, atms, mps, idx, cnrs, inv, das, fn) for idx in range(len(mps))]
	elif sum(tps[2:4])>0:
		fns = fns + list(range(len(fn)+1,len(fn)+6))
		rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
		if tps[2]==1:
			vt = [feature_AA(rns, feats, atms, mps, idx, cnrs, inv, das, fn) for idx in range(len(mps))]
		else:
			vt = [feature_PP(rns,feats, mmffp, mps, idx, cnrs, gap, inv, das, fn) for idx in range(len(mpd))]
	elif tps[4]==1:
		fns = fns + list(range(len(fn)+1,len(fn)+10))
		rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
		vt = [feature_AP(rns, feats, atms, mmffp, gap,mps, idx, cnrs, inv, das, fn) for idx in range(len(mpd))]			
	else:
		rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
		vt = [pairLr_dist(f, rns, cnrs, inv, das[i], fn) for i, f in enumerate(feats)]
	vdf = pd.DataFrame(vt)
	vdf.index = mn
	vdf['target'] = sdf.split('/')[-2]
	return vdf


gap = 0.6
def AP_from_dud(info):
	sdf, tps = info
	factory = ChemicalFeatures.BuildFeatureFactory('{}/etc/BaseFeatures.fdef'.format(DIR_AP))
	fn = factory.GetFeatureFamilies()	
	fns = [i for i, f in enumerate(fn) if f in use_feats]
	mm = Chem.SDMolSupplier(sdf, removeHs=True) 
	mm = [m for m in mm if m is not None]
	# mm = [m for m in mm if m.GetProp('_Name') in cids]
	mn =  [m.GetProp('_Name') for m in mm]
	feats = [factory.GetFeaturesForMol(m) for m in mm]
	das = [[f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')] for feat in feats]
	das = [unique([x for x in da if da.count(x) > 1]) for da in das]
	rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
	vt = [pairLr_dist(f, rns, cnrs, inv, das[i], fn) for i, f in enumerate(feats)]
	vdf = pd.DataFrame(vt)
	vdf.index = mn
	vdf['target'] = sdf.split('/')[-2]
	return vdf




gap = 0.6
def AP_from_cid(info):
	cids, sdf, tps = info
	factory = ChemicalFeatures.BuildFeatureFactory('{}/etc/BaseFeatures.fdef'.format(DIR_AP))
	fn = factory.GetFeatureFamilies()	
	fns = [i for i, f in enumerate(fn) if f in use_feats]
	mm = Chem.SDMolSupplier(sdf, removeHs=True) 
	mm = [m for m in mm if m is not None]
	mm = [m for m in mm if m.GetPropsAsDict()['PUBCHEM_COMPOUND_CID'] in cids]
	mn =  [m.GetPropsAsDict()['PUBCHEM_COMPOUND_CID'] for m in mm]
	feats = [factory.GetFeaturesForMol(m) for m in mm]
	das = [[f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')] for feat in feats]
	das = [unique([x for x in da if da.count(x) > 1]) for da in das]
	rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
	vt = [pairLr_dist(f, rns, cnrs, inv, das[i], fn) for i, f in enumerate(feats)]
	vdf = pd.DataFrame(vt)
	vdf.index = mn
	return vdf



def sdf2fp(info):
	cids, sdf, tps = info
	mm = Chem.SDMolSupplier(sdf, removeHs=False) 
	mm = [m for m in mm if m is not None]
	mm = [m for m in mm if m.GetPropsAsDict()['PUBCHEM_COMPOUND_CID'] in cids]
	fp = [list(Chem.rdMolDescriptors.GetMorganFingerprintAsBitVect(m, radius)) for m in mm]
	fps = pd.DataFrame(fp)
	fps.index = cids
	return fps


def sdf2fp(sdfs):
	mm = Chem.SDMolSupplier(sdfs, removeHs=False) 
	mn = [m.GetPropsAsDict()['PUBCHEM_COMPOUND_CID'] for m in mm if m is not None]
	fp = [list(Chem.rdMolDescriptors.GetMorganFingerprintAsBitVect(m, radius)) for m in mm if m is not None]
	fps = pd.DataFrame(fp)
	fps.index = mn
	return fps

gap = 0.6
def FP_from_cid(info):
	cids, sdf, tps = info
	mm = Chem.SDMolSupplier(sdf, removeHs=True) 
	mm = [m for m in mm if m is not None]
	mm = [m for m in mm if m.GetPropsAsDict()['PUBCHEM_COMPOUND_CID'] in cids]
	mn =  [m.GetPropsAsDict()['PUBCHEM_COMPOUND_CID'] for m in mm]
	fp = [list(Chem.rdMolDescriptors.GetMorganFingerprintAsBitVect(m, radius)) for m in mm if m is not None]
	fps = pd.DataFrame(fp)
	fps.index = mn
	return fps


cid_list = []
cidss = [[] for i in range(6358)]
for cid in cid_list:
	if cid < 158950000:
		cidss[int(cid // 25000)].append(cid)


str_dir = '/spstorage/DB/PUBCHEM/3D_1conf'

sdfs = ['{}/{}_{}.sdf'.format(str_dir,str(i * 25000 +1).zfill(8), str((i + 1) * 25000).zfill(8)) for i,cd in enumerate(cidss)
 if len(cd)!=0]
cds = [cd for cd in cidss if len(cd)!=0]

infos = [[cds[i], sdfs[i]] for i in range(len(cds))]


radius = 2 
def FP_from_cid(info):
	cids, sdf = info
	mm = Chem.SDMolSupplier(sdf, removeHs=True) 
	mm = [m for m in mm if m is not None]
	mm = [m for m in mm if m.GetPropsAsDict()['PUBCHEM_COMPOUND_CID'] in cids]
	mn =  [m.GetPropsAsDict()['PUBCHEM_COMPOUND_CID'] for m in mm]
	fp = [list(Chem.rdMolDescriptors.GetMorganFingerprintAsBitVect(m, radius)) for m in mm]
	fps = pd.DataFrame(fp)
	fps.index = mn
	return fps



start_time = time.time()
fps = mclapply(infos,FP_from_cid, ncpu= 100) # 100 / 10 855s   / 700 / 20 --> 57m | 3466s
print("---{}s seconds---".format(time.time()-start_time))

fff = pd.concat(fps)
fff.to_pickle(DIR_AP+'/dti_cadfs80_fp.pkl')
fff.to_pickle(dtip+'/FP.pkl')

lltt2.to_pickle(DIR_AP+'/dili_AP2.pkl')




	feats = [factory.GetFeaturesForMol(m) for m in mm]
	das = [[f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')] for feat in feats]
	das = [unique([x for x in da if da.count(x) > 1]) for da in das]
	if sum(tps[1:5])>0:
		mps = [m.GetConformer().GetPositions() for m in mm]	
		if tps[3]!=1:
			atms = [m.GetAtoms() for m in mm]
		if sum(tps[3:5])>0:
			mpd = [m.GetPropsAsDict() for m in mm]
			mmffp = [m['PUBCHEM_MMFF94_PARTIAL_CHARGES'] if 'PUBCHEM_MMFF94_PARTIAL_CHARGES' in m.keys() else '' for m in mpd]
			mmffp = [m.split('\n') if (m!=0) & (m!='') else '' for m in mmffp]		
	# featur n selection
	if tps[0]==1:
		fns = fns + [len(fn)]
	if tps[1]==1:
		fns = fns + list(range(len(fn)+1,len(fn)+3))
		rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
		vt = [feature_AH(rns, feats, atms, mps, idx, cnrs, inv, das, fn) for idx in range(len(mps))]
	elif sum(tps[2:4])>0:
		fns = fns + list(range(len(fn)+1,len(fn)+6))
		rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
		if tps[2]==1:
			vt = [feature_AA(rns, feats, atms, mps, idx, cnrs, inv, das, fn) for idx in range(len(mps))]
		else:
			vt = [feature_PP(rns,feats, mmffp, mps, idx, cnrs, gap, inv, das, fn) for idx in range(len(mpd))]
	elif tps[4]==1:
		fns = fns + list(range(len(fn)+1,len(fn)+10))
		rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
		vt = [feature_AP(rns, feats, atms, mmffp, gap,mps, idx, cnrs, inv, das, fn) for idx in range(len(mpd))]			
	else:
		vt = [pairLr_dist(f, rns, cnrs, inv, das[i], fn) for i, f in enumerate(feats)]
	vdf = pd.DataFrame(vt)
	vdf.index = mn
	return vdf





##### from ZINC DB

# zinc_list = ~~
# /spstorage/DB/ZINC/drug-like-clean-substance-3D/HC/AARL/HCAARL.xaa.sdf

use_feats = ['Donor', 'Acceptor', 'NegIonizable', 'PosIonizable', 'Halogen', 'Aromatic', 'Hydrophobe']
# D+A type use / AH / AA / PP / AP
tps = [0,0,1,0,1] 
tps = [1,1,0,0,0]

zinc_dir = '/spstorage/DB/ZINC/drug-like-clean-substance-3D'
zinc_list = list(zcdf.sample(frac=1)[0:100].zinc_id)
zinc_list = list(zinc1.zinc_id)

szcdf = szcdf[(szcdf.mw>180) & (szcdf.mw<1000)].sample(frac=1, random_state = 7)[0:10**6]
zinc_list = list(szcdf.zinc_id)
zcdf=pd.read_pickle('/spstorage/DB/ZINC/zinc_tranche_3D.pkl')
# 1 

start_time = time.time()
vvv = mclapply(infos[0:700],AP_from_cid, ncpu= 20) # 100 / 10 855s   / 700 / 20 --> 40m?
print("---{}s seconds---".format(time.time()-start_time))

start_time = time.time()
tranches = list(set(zcdf[zcdf.zinc_id.isin(zinc_list)].tranche_name_3D))
print("---{}s seconds---".format(time.time()-start_time))
start_time = time.time()
zidss = [list(zcdf[(zcdf.tranche_name_3D==tc)&(zcdf.zinc_id.isin(zinc_list))].zinc_id) for tc in tranches]
print("---{}s seconds---".format(time.time()-start_time))

# 2 --> 665
start_time = time.time()
tranches = list(set(szcdf.tranche_name_3D))
zidss = [[] for i in range(len(tranches))]
print("---{}s seconds---".format(time.time()-start_time))
start_time = time.time()
for zid in zinc_list:
	zidss[tranches.index(list(szcdf[szcdf.zinc_id==zid].tranche_name_3D)[0])].append(zid)



print("---{}s seconds---".format(time.time()-start_time))
start_time = time.time()
tranches = [tc for i,tc in enumerate(tranches) if len(zidss[i])!=0]
print("---{}s seconds---".format(time.time()-start_time))
ttranches = [tc if 'sdf' not in tc else tc.split('.')[0] for tc in tranches]


for cid in cid_list:
	if cid < 158950000:
		cidss[cid // 25000].append(cid)

# sdfs = ['{}/{}/{}/{}'.format(zinc_dir,tc[0:2],tc[2:6],tc) if 'sdf' in tc else '{}/{}/{}/{}.xaa.sdf'.format(zinc_dir,tc[0:2],tc[2:],tc) for tc in tranches]

sdfs = ['{}/{}/{}/{}.sdf'.format(zinc_dir,tc[0:2],tc[2:6],tc) for tc in tranches]

zds = [cd for cd in zidss if len(cd)!=0]

infos = [[['ZINC'+str(int(zz)).zfill(12) for zz in zds[i]], sdfs[i], tps] for i in range(len(zds))]

zinc_list = [5.606512e+06,7.261974e+08,1.197769e+09]
set(zcdf[zcdf.zinc_id.isin(zinc_list)].tranche_name_3D)

('/spstorage/DB/ZINC/drug-like-clean-substance-3D/FA/ABRL/FAABRL.xaa.sdf')
sdf = ('/spstorage/DB/ZINC/drug-like-clean-substance-3D/BC/ADRM/BCADRM.xaa.sdf')

gap = 0.6
def AP_from_zid(info):
	cids, sdf, tps = info
	factory = ChemicalFeatures.BuildFeatureFactory('{}/etc/BaseFeatures_h.fdef'.format(DIR_AP))
	fn = factory.GetFeatureFamilies()	
	fns = [i for i, f in enumerate(fn) if f in use_feats]
	mm = Chem.SDMolSupplier(sdf, removeHs=True) 
	mm = [m for m in mm if m is not None]
	mm = [m for m in mm if m.GetProp('_Name') in cids]
	mn =  [m.GetProp('_Name') for m in mm]
	feats = [factory.GetFeaturesForMol(m) for m in mm]
	das = [[f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')] for feat in feats]
	das = [unique([x for x in da if da.count(x) > 1]) for da in das]
	if tps[0]==1:
		fns = fns + [len(fn)]
	if sum(tps[1:3])>0:
		mps = [m.GetConformer().GetPositions() for m in mm]	
		atms = [m.GetAtoms() for m in mm]	
		# featur n selection
		if tps[1]==1:
			fns = fns + list(range(len(fn)+1,len(fn)+3))
			rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
			vt = [feature_AH(rns, feats, atms, mps, idx, cnrs, inv, das, fn) for idx in range(len(mps))]
		else:
			fns = fns + list(range(len(fn)+1,len(fn)+6))
			rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
			vt = [feature_AA(rns, feats, atms, mps, idx, cnrs, inv, das, fn) for idx in range(len(mps))]
	else:
		vt = [pairLr_dist(f, rns, cnrs, inv, das[i], fn) for i, f in enumerate(feats)]
	vdf = pd.DataFrame(vt)
	vdf.index = mn
	return vdf

len(infos)
with open(main+'/pdb_zinc_infos.pkl','wb') as f:
	pickle.dump(infos,f)

with open(main+'/pdb_zinc_infos.pkl','rb') as f:
	infos = pickle.load(f)


start_time = time.time()
zzz = mclapply(iinfos,AP_from_zid, ncpu= 100) # 100 / 10 855s   / 700 / 20 --> 40m?
print("---{}s seconds---".format(time.time()-start_time))



gap = 0.6
def FP_from_zid(info):
	cids, sdf, tps = info
	mm = Chem.SDMolSupplier(sdf, removeHs=True) 
	mm = [m for m in mm if m is not None]
	mm = [m for m in mm if m.GetProp('_Name') in cids]
	mn =  [m.GetProp('_Name') for m in mm]
	fp = [list(Chem.rdMolDescriptors.GetMorganFingerprintAsBitVect(m, radius)) for m in mm]
	fps = pd.DataFrame(fp)
	fps.index = mn
	return fps


iinfos = [[[a for a in info[0] if a in zdd],info[1],[1,1,0,0,0]] for info in infos if len([a for a in info[0] if a in zdd])!=0]

gap = 0.6
def FP_from_zid(info):
	cids, sdf, tps = info
	mm = Chem.SDMolSupplier(sdf, removeHs=True) 
	mm = [m for m in mm if m is not None]
	mm = [m for m in mm if m.GetProp('_Name') in cids]
	mn =  [m.GetProp('_Name') for m in mm]
	fp = [list(Chem.rdMolDescriptors.GetMorganFingerprintAsBitVect(m, radius)) for m in mm]
	fps = pd.DataFrame(fp)
	fps.index = mn
	scfm = [GetScaffoldForMol(m) for m in mm]
	scf = [list(Chem.rdMolDescriptors.GetMorganFingerprintAsBitVect(cff, radius)) for cff in scfm]
	scfs = pd.DataFrame(scf)
	scfs.index = mn
	return tuple([fps, scfs])

start_time = time.time()
fps = mclapply(iinfos,FP_from_zid, ncpu= 100) # 100 / 10 855s   / 700 / 20 --> 57m | 3466s
print("---{}s seconds---".format(time.time()-start_time))

fpss = pd.concat(fps)



fpss.to_pickle('/spstorage/USERS/gina/Project/AP/zzff.pkl')