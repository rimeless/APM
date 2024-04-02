
main = DIR_AP+'/final2'

pdb2uni = pd.read_csv(DIR_AP+'/DB/gpdb_ligand/pdb2uniprot.tab',sep='\t') #  4383

pdb2etid = select_in('BCDTv2_TARGET.EntrezUniprotKB2etid', '*', 'UniProtKB', unique(pdb2uni.Entry))
# uni2etid = select_in('BCDTv2_TARGET.EntrezUniprotKB2etid', '*', 'UniProtKB', unique(pdb2uni.Entry))


pdb2uni.columns= ['PDB','Entry','Entry name', 'Status', 'Protein names', 'Gene names', 'Organism','Length']

df=pdb2uni.join(pdb2uni.pop('PDB').str.split(',',expand=True))
df = df.set_index(list(df.columns)[0:7]).stack().reset_index(level=list(range(7)),name='PDB')

df2etid = pd.merge(df, pdb2etid, left_on='Entry', right_on='UniProtKB')

df2etid = df2etid.drop_duplicates()

pdb_dict = dict(zip(list(df2etid.PDB),list(df2etid.etid)))




gidx = pd.read_csv(DIR_AP+'/DB/pdb_ligand/INDEX_general_PL_data.2019')
gidx = pd.read_csv(DIR_AP+'/DB/v2020-other-PL/index/INDEX_general_PL_data.2020')



rgidx= pd.DataFrame([[st for st in et.split(' ') if (st!='')&(st!='//')] for et in list(gidx.iloc[:,0])])
rgidx = rgidx.iloc[:,list(range(7))]
rgidx.columns = gidx.columns
rgidx['ligand name']=[et.split('(')[1].split(')')[0] for et in list(rgidx['ligand name'])]
# rgidx['Kd/Ki']=[float(et.split('=')[1][:-2])*10**(3*(1-scl.index(et.split('=')[1][-2:]))) for et in list(rgidx['Kd/Ki'])]
rgidx['-logKd/Ki']=[float(et) for et in list(rgidx['-logKd/Ki'])]
rgidx.columns = ['PDB', 'resolution', 'release year', '-logKd/Ki', 'Kd/Ki', 'reference', 'cmp']


# with open(main+'/all_cmp.txt', "wt") as f:f.write(','.join([cd for cd in unique(rgidx.cmp) if len(cd)<=3]))


df2etid[df2etid.PDB.isin(rgidx[4285:].PDB)]

pdb2cmp = pd.merge(df2etid.loc[:,['PDB','etid']], rgidx.loc[:,['PDB','cmp','-logKd/Ki']]).drop_duplicates()
pdb2cmp = pdb2cmp.sort_values('-logKd/Ki')
apdb2cmp = pdb2cmp[pdb2cmp['-logKd/Ki']>5]
ipdb2cmp = pdb2cmp[pdb2cmp['-logKd/Ki']<=4.70]

acts_pdb = list(rgidx[rgidx['-logKd/Ki']>=5].PDB)
inacts_pdb = list(rgidx[rgidx['-logKd/Ki']<=4.70].PDB)

mers = [cd for cd in unique(pdb2cmp.cmp) if 'mer' in cd]
apdb2cmp[~apdb2cmp.cmp.isin(mers)].loc[:,['PDB','cmp']].drop_duplicates()
ipdb2cmp[~ipdb2cmp.cmp.isin(mers)].loc[:,['PDB','cmp']].drop_duplicates()
mapdb2cmp = apdb2cmp[~apdb2cmp.cmp.isin(mers)]

mipdb2cmp = ipdb2cmp[~ipdb2cmp.cmp.isin(mers)]



augipdb2cmp = pd.DataFrame()
for cd in unique(pdb2cmp.cmp):
	iet = ipdb2cmp[ipdb2cmp.cmp==cd]
	iet = iet[~iet.etid.isin(apdb2cmp[apdb2cmp.cmp==cd].etid)]
	newet= df2etid[(df2etid.etid.isin(iet[0:3].etid)) & (~df2etid.PDB.isin(iet.PDB))]
	if len(apdb2cmp[apdb2cmp.cmp==cd]) > len(newet):
		ett = list(newet[newet.etid==newet.etid].PDB) 
	else:
		ett = [list(newet[newet.etid==et].sample(3,random_state=7).PDB) if len(newet[newet.etid==et])>3 else list(newet[newet.etid==et].PDB) for et in unique(newet.etid)]
		ett = list(itertools.chain(*ett))
	augipdb2cmp = pd.concat([augipdb2cmp,pd.DataFrame({'PDB':ett, 'cmp':cd})])


augipdb2cmp = pd.concat([ipdb2cmp.loc[:,['PDB','cmp']], augipdb2cmp]).drop_duplicates()

augipdb2cmp[~augipdb2cmp.cmp.isin(mers)].loc[:,['PDB','cmp']].drop_duplicates()

gsets= gset[gset.cid.isin(app.index)]

augipdb2cid = pd.DataFrame()
cids = []
for et in unique(pdb2cmp.etid):
	iet = pdb2cmp[pdb2cmp.etid==et]
	iet = iet.sort_values('-logKd/Ki')
	ett = all_int[(all_int.etid==et)&(all_int.standard_value>20000)]
	ett = ett.sort_values('standard_value', ascending=False)
	ett = ett.loc[:,['Entrez','cid']].drop_duplicates()
	cds = list(ett.cid)
	cds = [cd for cd in cds if cd not in cids]
	hlen = len(mapdb2cmp[mapdb2cmp.etid==et]) - len(mipdb2cmp[mipdb2cmp.etid==et])
	hlen = int(hlen*1.5)
	if hlen>0:
		cds = cds[0:hlen]
		etp = list(list(iet.PDB) * hlen)[0:len(cds)]
		augipdb2cid = pd.concat([augipdb2cid,pd.DataFrame({'PDB':etp, 'etid':et, 'cid':cds})])
		cids.extend(cds)


fv = [str(int(cd)) for cd in list(pd.concat([augipdb2cid[~augipdb2cid.cid.isin(ff)],p5[~p5.cid.isin(apms.index)]]).cid) if not np.isnan(cd)]
fv = fv + [str(int(cd)) for cd in unique(kc[~kc.cid.isin(apms.index)].cid) if cd not in ff] + ['2051','5113385','3797', '5339183']
# with open(main+'/aug_cids2.txt', "wt") as f:f.write('\n'.join([str(int(cd)) for cd in list(augipdb2cid[~augipdb2cid.cid.isin(ff)].cid) if not np.isnan(cd)]))
with open(main+'/aug_cids2.txt', "wt") as f:f.write('\n'.join(fv))

with open(main+'/aug_cids.txt', "r") as f:ff= f.readlines()
ff=[int(fc) for fc in ff]

mapdb2cmp['label'] = 1
mipdb2cmp['label'] = 0
augipdb2cid['label'] = 0

set5 = pd.merge(pd.concat([mapdb2cmp, mipdb2cmp]),cmp2cid, on='cmp')

sets5 = pd.concat([set5, augipdb2cid]).loc[:,['PDB','cmp','cid','label']].drop_duplicates()
augipdb2cid

dfsets = pd.merge(df2etid, sets5)
dfsets = pd.merge(dfsets, dsts)
dfsets[dfsets.UniProtKB.isin(dsts[dsts.sims>0.40].UniProtKB)]

list(dfsets[dfsets.UniProtKB.isin(dsts[dsts.sims>0.22].UniProtKB)].label).count(1)
list(dfsets[dfsets.UniProtKB.isin(dsts[dsts.sims>0.22].UniProtKB)].label).count(0)
list(dfsets[~dfsets.UniProtKB.isin(dsts[dsts.sims>0.07].UniProtKB)].label).count(1)
list(dfsets[~dfsets.UniProtKB.isin(dsts[dsts.sims>0.07].UniProtKB)].label).count(0)

dfsets['tag']='test'
dfsets.loc[dfsets.sims>0.09,'tag']='val'
dfsets.loc[dfsets.sims>0.26,'tag']='train'

dfsets[dfsets.tag=='train']
dfsets[dfsets.tag=='val']

dfsets[dfsets.tag=='test']


tr = list(dfsets[dfsets.tag!='test'].cid)
dfsets[(dfsets.tag=='test') & (dfsets.cid.isin(tr))]

sets5 = dfsets.loc[:,['PDB','cmp','cid','tag','label']]
sets5 = sets5[~np.isnan(sets5.cid)]
sets5 = sets5.drop_duplicates()

list(sets5[sets5.tag=='train'].label).count(1)
list(sets5[sets5.tag=='val'].label).count(1)

list(sets5[sets5.tag=='test'].label).count(1)
sets5 = sets5[~((sets5.tag=='test') & (sets5.cid.isin(tr)))]



sets5.to_csv(DIR_AP+'/final2/ptncid52.csv')

# filtered
cid_info = pd.read_csv(main+'/cid_info.csv')
cid_info = pd.concat([cid_info, pd.read_csv(main+'/aug_cids_info.csv')])


cid_info = cid_info.loc[:,['cid', 'mw', 'rotbonds', 'hbonddonor','hbondacc','cidcdate', 'heavycnt', 'xlogp']]

subs = cid_info[cid_info.xlogp<=5]
subs = subs[subs.mw<=500]
subs = subs[subs.hbonddonor<=5]
subs = subs[subs.hbondacc<=10]
subs = subs[subs.mw>=180]
subs = subs[subs.rotbonds<=10]

fsets5 = sets5[sets5.cid.isin(subs.cid)]
fsets5[(fsets5.tag=='train') & (fsets5.label==1)]
fsets5[(fsets5.tag=='test') & (fsets5.label==1)]
fsets5[(fsets5.tag=='val') & (fsets5.label==1)]
fsets5.to_csv(DIR_AP+'/final2/fptncid5.csv')

## rigidity
sets5[sets5.cid.isin(cid_info[(cid_info.rotbonds>=5)&(cid_info.rotbonds<10)].cid)]
sets5[sets5.cid.isin(cid_info[(cid_info.rotbonds<5)].cid)]

list(sets5[sets5.cid.isin(cid_info[(cid_info.rotbonds<5)].cid)].label).count(0)
list(sets5[sets5.cid.isin(cid_info[(cid_info.rotbonds>=5)&(cid_info.rotbonds<10)].cid)].label).count(0)

list(sets5[(sets5.cid.isin(cid_info[(cid_info.rotbonds<5)].cid))&(sets5.tag=='train')].label).count(0)
list(sets5[sets5.cid.isin(cid_info[(cid_info.rotbonds>=5)&(cid_info.rotbonds<10)].cid)&(sets5.tag=='train')].label).count(0)

sets5[sets5.cid.isin(cid_info[(cid_info.rotbonds<5)].cid)].to_csv(DIR_AP+'/final2/ptncid5_rigid.csv')
sets5[sets5.cid.isin(cid_info[(cid_info.rotbonds>=5)&(cid_info.rotbonds<10)].cid)].to_csv(DIR_AP+'/final2/ptncid5_flex.csv')

augipdb2cmp = pd.concat([mipdb2cmp.loc[:,['PDB','cmp']], augipdb2cmp]).drop_duplicates()



## rigidity
# sets3 = pd.read_csv(DIR_AP+'/final2/ptncid3.csv')
# fsets3r = pd.read_csv(DIR_AP+'/final2/fptncid3_rigid.csv')
# sets3[sets3.cid.isin(cid_info[(cid_info.rotbonds>=5)&(cid_info.rotbonds<10)].cid)]
# sets3[sets3.cid.isin(cid_info[(cid_info.rotbonds<5)].cid)]

# list(sets3[sets3.cid.isin(cid_info[(cid_info.rotbonds<5)].cid)].label).count(0)
# list(sets3[sets3.cid.isin(cid_info[(cid_info.rotbonds>=5)&(cid_info.rotbonds<10)].cid)].label).count(0)

# list(sets3[(sets3.cid.isin(cid_info[(cid_info.rotbonds<5)].cid))&(sets3.tag=='train')].label).count(0)
# list(sets3[sets3.cid.isin(cid_info[(cid_info.rotbonds>=5)&(cid_info.rotbonds<10)].cid)&(sets3.tag=='train')].label).count(0)

# sets3[sets3.cid.isin(cid_info[(cid_info.rotbonds<5)].cid)].to_csv(DIR_AP+'/final2/ptncid5_rigid.csv')
# sets3[sets3.cid.isin(cid_info[(cid_info.rotbonds>=5)&(cid_info.rotbonds<10)].cid)].to_csv(DIR_AP+'/final2/ptncid5_flex.csv')

### 



pd.concat([rgidx[rgidx['-logKd/Ki']>5].loc[:,['PDB','cmp']], inacts_aug])

cmp2cid = pd.read_csv(main+'/cmp2cid.txt', sep='\t',header=None)
cmp2cid.columns= ['cmp','cid']

rgidx =pd.merge(rgidx,cmp2cid)
rgidx.index = list(rgidx.PDB)
all_cid = list(non_acts[non_acts.etid.isin(gset.etid)].cid)+list(gset.cid)
pdb_cid = list(set(rgidx.cid)-set(all_cid)-set(cid2ecid.cid))

ccc= [str(int(cd)) for cd in unique(pdb_cid) if not np.isnan(cd)]
ccc='\n'.join(ccc)


c2 =DIR_AP+'/DB/refined-set'

pds=[pps.split('/')[-1] for pps in glob.glob(c2+'/*')]

c2 = DIR_AP+'/DB/v2019-other-PL'
pds = pds + [pps.split('/')[-1] for pps in glob.glob(c2+'/*')]


a = [[8,20,40,80],[0,1],[0,1]]
oplist = list(itertools.product(*a))

# factory = ChemicalFeatures.BuildFeatureFactory('/spstorage/USERS/gina/Project/AP/etc/BaseFeatures.fdef')
# fn1= factory.GetFeatureFamilies()
rn6=[str(x)+str(y) for x,y in combinations_with_replacement([0,1,2,3,5,6],2)]
rn7 = [str(x)+str(y) for x,y in combinations_with_replacement([0,1,2,3,5,6,9],2)]
rnh = [str(x)+str(y) for x,y in combinations_with_replacement([0,1,2,3,4,6,7],2)]
rn_list = [rn6,rn7,rnh]
f_list=['.fdef','_h.fdef','_h.fdef']


for facts in oplist:
	rns = rn_list[facts[2]]
	factory = ChemicalFeatures.BuildFeatureFactory('/spstorage/USERS/gina/Project/AP/etc/BaseFeatures'+f_list[facts[2]])
	fn = factory.GetFeatureFamilies()
	bns = 30/facts[0]
	cns = [bns*g for g in list(range(facts[0]))]
	st = 30/(inv**facts[0]) # 16 bin
	cnrs = [st*(inv**g) for g in list(range(facts[0]))]
	bin_type = facts[1]
	def bulks_pdb(pps):
		try:
			mm = Chem.rdmolfiles.MolFromPDBFile(DIR_AP+'/DB/v2019-other-PL/'+pps+'/'+pps+'_pocket.pdb', removeHs=True)
			feats = factory.GetFeaturesForMol(mm)
			if bin_type==0:
				vt = pairLr_dist(feats, rns, cnrs, inv, fn)
			else:
				vt = pairLs_dist(feats, rns, cns, bns, fn)
		except:
			vt = [0]*len(rns)*len(cns)
		return tuple(vt)
	aps = mclapply(pds, bulks_pdb, ncpu=120)
	app = pd.DataFrame(aps)
	app.index= pds
	app.to_pickle('{}/DB/pdb_ligand/pocket/pocket_{}_{}_{}_30.pkl'.format(DIR_AP, facts[0], facts[1], facts[2]))
	del app


sym2feat = {'C':10,'O':11,'S':12}
def feature_AA(feats, atms, mps, col, inv, das, fn):
	# basic charac. extra atom to atom charac. 10 C / 11 O / 12 S / 13 etc
	row = [str(x)+str(y) for x,y in combinations_with_replacement([0,1,2,3,4,6,7,9,10,11,12,13],2)]
	ft_idx = unique(list(itertools.chain(*[f.GetAtomIds() for f in feats])))
	uu_idx = [i for i in range(len(mps)) if i not in ft_idx]
	use_pos = np.array([ft.GetPos() for ft in feats])
	use_fam = [fn.index(ft.GetFamily()) if not ft.GetAtomIds() in das else 9 for ft in feats]
	if (len(uu_idx) != 0) & (len(use_pos) !=0):
		use_idx = np.array(itemgetter(*uu_idx)(mps))
		if len(uu_idx)>1:
			use_atms = [sym2feat[at.GetSymbol()] if at.GetSymbol() in ['C','O','S'] else 13 for at in itemgetter(*uu_idx)(atms)]
			poses = np.concatenate((use_pos,use_idx))
		else:
			atsym = itemgetter(*uu_idx)(atms).GetSymbol()
			poses = np.concatenate((use_pos,np.array([use_idx])))
			if atsym in ['C','O','S']:
				use_atms = [sym2feat[atsym]]
			else:
				use_atms = [13]
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




for facts in [(10,0,2),(10,1,2),(12,0,2),(12,1,2)]:
	rns = rn_list[facts[2]]
	factory = ChemicalFeatures.BuildFeatureFactory('/spstorage/USERS/gina/Project/AP/etc/BaseFeatures'+f_list[facts[2]])
	fn = factory.GetFeatureFamilies()
	bns = 20/facts[0]
	cns = [bns*g for g in list(range(facts[0]))]
	st = 20/(inv**facts[0]) # 16 bin
	cnrs = [st*(inv**g) for g in list(range(facts[0]))]
	filter_type = facts[1]
	def bulk_cids(sdf):
		mm = Chem.SDMolSupplier(sdf, removeHs=True) 
		mpd = [m.GetPropsAsDict() for m in mm if m is not None]
		atms = [m.GetAtoms() for m in mm if m is not None]
		mn = [m['PUBCHEM_COMPOUND_CID'] for m in mpd]
		mps = [m.GetConformer().GetPositions() for m in mm if m is not None]
		feats = [factory.GetFeaturesForMol(m) for m in mm if m is not None]
		das = [[]]*len(feats)
		if filter_type==0:
			feats = [filter_feats_hyp(feat) for feat in feats]
		else:
			feats = [filter_feats_renew(feat) for feat in feats]
		if facts[2]==2:
			das = [[f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')] for feat in feats]
			das = [unique([x for x in da if da.count(x) > 1]) for da in das]
			feats = [filter_feats_both(feat) for feat in feats]
		vt = feature_AA(feats, atms, mps, cnrs, inv, das, fn)
		# vt = [pairLr_dist(f, rns, cnrs, inv, das[i], fn) for i, f in enumerate(feats)]
		vdf = pd.DataFrame(vt)
		vdf.index = mn
		# vdf['target'] = sdf.split('_')[-1].split('.')[0]
		return vdf
	aps = mclapply(sdfs, bulk_cids, ncpu=120)
	app = pd.concat(aps)
	# app = pd.concat([prep,app])
	app.to_pickle('{}/final2/dataset/mtxs/aa_apms_{}_{}_{}.pkl'.format(DIR_AP, facts[0], facts[1], facts[2]))
	del app




for facts in [(10,1,2,4),(12,1,2,4),(10,1,2,5),(12,1,2,5)]:
	rns = row = [str(x)+str(y) for x,y in combinations_with_replacement([0,1,2,3,4,6,7,9,10,11,12,13],2)]
	factory = ChemicalFeatures.BuildFeatureFactory('/spstorage/USERS/gina/Project/AP/etc/BaseFeatures'+f_list[facts[2]])
	fn = factory.GetFeatureFamilies()
	bns = 30/facts[0]
	cns = [bns*g for g in list(range(facts[0]))]
	st = 30/(inv**facts[0]) # 16 bin
	cnrs = [st*(inv**g) for g in list(range(facts[0]))]
	def bulks_pdb(pps):
		try:
			ppm = Chem.rdmolfiles.MolFromPDBFile(DIR_AP+'/DB/refined-set/'+pps+'/'+pps+'_pocket.pdb', removeHs=True)
			llm = Chem.SDMolSupplier(DIR_AP+'/DB/refined-set/'+pps+'/'+pps+'_ligand.sdf', removeHs=True)[0]
			if llm is None:
				llm = Chem.rdmolfiles.MolFromMol2File(DIR_AP+'/DB/refined-set/'+pps+'/'+pps+'_ligand.mol2', removeHs=True)
			llf = factory.GetFeaturesForMol(llm)
			lluse_pos = np.array([ft.GetPos() for ft in llf])
			ppf = factory.GetFeaturesForMol(ppm)
			ppuse_pos = np.array([ft.GetPos() for ft in ppf])
			atms = ppm.GetAtoms()
			mps = ppm.GetConformer().GetPositions()
			atms = [a for i,a in enumerate(atms) if min([dst([mps[i],lgd_po]) for lgd_po in lluse_pos])<facts[3]]
			mps = [m for i,m in enumerate(mps) if min([dst([m,lgd_po]) for lgd_po in lluse_pos])<facts[3]]
			feats = [f for i,f in enumerate(ppf) if min([dst([ppuse_pos[i],lgd_po]) for lgd_po in lluse_pos])<facts[3]]
			feats = filter_feats_renew(feats)
			if facts[2]==2:
				das = [f.GetAtomIds() for f in feats if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')]
				das = unique([x for x in das if das.count(x) > 1])
				feats = filter_feats_both(feats)
			vt = feature_AA(feats, atms, mps, cnrs, inv, das, fn)
		except:
			vt = [0]*len(rns)*len(cns)
		return tuple(vt)
	aps = mclapply(pds, bulks_pdb, ncpu=50)
	app = pd.DataFrame(aps)
	app.index= pds
	app.to_pickle('{}/final2/dataset/pdb/binding/pockets/aa_pockets_{}_{}_{}_{}.pkl'.format(DIR_AP, facts[0], facts[1], facts[2], facts[3]))
	del app


for facts in [(10,1,2,4),(12,1,2,4),(10,1,2,5),(12,1,2,5)]:
	rns = row = [str(x)+str(y) for x,y in combinations_with_replacement([0,1,2,3,4,6,7,9,10,11,12,13],2)]
	factory = ChemicalFeatures.BuildFeatureFactory('/spstorage/USERS/gina/Project/AP/etc/BaseFeatures'+f_list[facts[2]])
	fn = factory.GetFeatureFamilies()
	bns = 30/facts[0]
	cns = [bns*g for g in list(range(facts[0]))]
	st = 30/(inv**facts[0]) # 16 bin
	cnrs = [st*(inv**g) for g in list(range(facts[0]))]
	def bulks_pdb(pps):
		try:
			ppm = Chem.rdmolfiles.MolFromPDBFile(DIR_AP+'/DB/v2019-other-PL/'+pps+'/'+pps+'_pocket.pdb', removeHs=True)
			llm = Chem.SDMolSupplier(DIR_AP+'/DB/v2019-other-PL/'+pps+'/'+pps+'_ligand.sdf', removeHs=True)[0]
			if llm is None:
				llm = Chem.rdmolfiles.MolFromMol2File(DIR_AP+'/DB/v2019-other-PL/'+pps+'/'+pps+'_ligand.mol2', removeHs=True)
			llf = factory.GetFeaturesForMol(llm)
			lluse_pos = np.array([ft.GetPos() for ft in llf])
			ppf = factory.GetFeaturesForMol(ppm)
			ppuse_pos = np.array([ft.GetPos() for ft in ppf])
			atms = ppm.GetAtoms()
			mps = ppm.GetConformer().GetPositions()
			atms = [a for i,a in enumerate(atms) if min([dst([mps[i],lgd_po]) for lgd_po in lluse_pos])<facts[3]]
			mps = [m for i,m in enumerate(mps) if min([dst([m,lgd_po]) for lgd_po in lluse_pos])<facts[3]]
			feats = [f for i,f in enumerate(ppf) if min([dst([ppuse_pos[i],lgd_po]) for lgd_po in lluse_pos])<facts[3]]
			feats = filter_feats_renew(feats)
			if facts[2]==2:
				das = [f.GetAtomIds() for f in feats if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')]
				das = unique([x for x in das if das.count(x) > 1])
				feats = filter_feats_both(feats)
			vt = feature_AA(feats, atms, mps, cnrs, inv, das, fn)
		except:
			vt = [0]*len(rns)*len(cns)
		return tuple(vt)
	aps = mclapply(pds, bulks_pdb, ncpu=25)
	app = pd.DataFrame(aps)
	app.index= pds
	app.to_pickle('{}/final2/dataset/pdb/binding/pockets/aa_pockets_{}_{}_{}_{}_2.pkl'.format(DIR_AP, facts[0], facts[1], facts[2], facts[3]))
	del app

app = pd.read_pickle('{}/final2/dataset/pdb/binding/pockets/aa_pockets_{}_{}_{}_{}.pkl'.format(DIR_AP, facts[0], facts[1], facts[2], facts[3]))
	
for facts in [(10,1,1,4),(12,1,1,4),(10,1,2,4),(12,1,2,4),(10,1,2,5),(12,1,2,5)]:
	rns = rn_list[facts[2]]
	factory = ChemicalFeatures.BuildFeatureFactory('/spstorage/USERS/gina/Project/AP/etc/BaseFeatures'+f_list[facts[2]])
	fn = factory.GetFeatureFamilies()
	bns = 30/facts[0]
	cns = [bns*g for g in list(range(facts[0]))]
	st = 30/(inv**facts[0]) # 16 bin
	cnrs = [st*(inv**g) for g in list(range(facts[0]))]
	def bulks_pdb(pps):
		try:
			ppm = Chem.rdmolfiles.MolFromPDBFile(DIR_AP+'/DB/refined-set/'+pps+'/'+pps+'_pocket.pdb', removeHs=True)
			llm = Chem.SDMolSupplier(DIR_AP+'/DB/refined-set/'+pps+'/'+pps+'_ligand.sdf', removeHs=True)[0]
			if llm is None:
				llm = Chem.rdmolfiles.MolFromMol2File(DIR_AP+'/DB/refined-set/'+pps+'/'+pps+'_ligand.mol2', removeHs=True)
			llf = factory.GetFeaturesForMol(llm)
			lluse_pos = np.array([ft.GetPos() for ft in llf])
			ppf = factory.GetFeaturesForMol(ppm)
			ppuse_pos = np.array([ft.GetPos() for ft in ppf])
			feats = [f for i,f in enumerate(ppf) if min([dst([ppuse_pos[i],lgd_po]) for lgd_po in lluse_pos])<facts[3]]
			feats = filter_feats_renew(feats)
			if facts[2]==1:
				das = [f.GetAtomIds() for f in feats if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')]
				das = unique([x for x in das if das.count(x) > 1])
				feats = filter_feats_both(feats)
			vt = pairLr_dist(f, rns, cnrs, inv, das, fn)
		except:
			vt = [0]*len(rns)*len(cns)
		return tuple(vt)
	aps = mclapply(pds, bulks_pdb, ncpu=25)
	app = pd.DataFrame(aps)
	app.index= pds
	app.to_pickle('{}/final2/dataset/pdb/binding/pockets/pockets_{}_{}_{}_{}.pkl'.format(DIR_AP, facts[0], facts[1], facts[2], facts[3]))
	del app


for facts in [(10,1,4),(10,1,5),(12,1,4),(12,1,5)]:
	rns = rn7
	factory = ChemicalFeatures.BuildFeatureFactory('/spstorage/USERS/gina/Project/AP/etc/BaseFeatures.fdef')
	fn = factory.GetFeatureFamilies()
	bns = 30/facts[0]
	cns = [bns*g for g in list(range(facts[0]))]
	st = 30/(inv**facts[0]) # 16 bin
	cnrs = [st*(inv**g) for g in list(range(facts[0]))]
	def bulks_pdb(pps):
		try:
			if pps not in pds1:
				ppm = Chem.rdmolfiles.MolFromPDBFile(DIR_AP+'/DB/v2019-other-PL/'+pps+'/'+pps+'_pocket.pdb', removeHs=True)
			else:
				ppm = Chem.rdmolfiles.MolFromPDBFile(DIR_AP+'/DB/refined-set/'+pps+'/'+pps+'_pocket.pdb', removeHs=True)
			llm = Chem.SDMolSupplier(DIR_AP+'/DB/v2019-other-PL/'+pps+'/'+pps+'_ligand.sdf', removeHs=True)[0]
			if llm is None:
				llm = Chem.rdmolfiles.MolFromMol2File(DIR_AP+'/DB/v2019-other-PL/'+pps+'/'+pps+'_ligand.mol2', removeHs=True)
			llf = factory.GetFeaturesForMol(llm)
			lluse_pos = np.array([ft.GetPos() for ft in llf])
			ppf = factory.GetFeaturesForMol(ppm)
			ppuse_pos = np.array([ft.GetPos() for ft in ppf])
			feats = [f for i,f in enumerate(ppf) if min([dst([ppuse_pos[i],lgd_po]) for lgd_po in lluse_pos])<facts[2]]
			feats = filter_feats_renew(feats)
			das = [f.GetAtomIds() for f in feats if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')]
			das = unique([x for x in das if das.count(x) > 1])
			vt = pairLr_dist(f, rns, cnrs, inv, das, fn)
		except:
			vt = [0]*len(rns)*len(cns)
		return tuple(vt)
	aps = mclapply(pds, bulks_pdb, ncpu=50)
	app = pd.DataFrame(aps)
	app.index= pds
	app.to_pickle('{}/final2/dataset/pdb/binding/pockets/naive_pockets_{}_{}_{}.pkl'.format(DIR_AP, facts[0], facts[1], facts[2]))
	del app

ppm = Chem.rdmolfiles.MolFromMol2File(DIR_AP+'/kalasanty/predicted_pockets/'+poc+'/pocket0.mol2', removeHs=True)

c2 = DIR_AP+'/kalasanty/predicted_pockets'
# pds = glob.glob(DIR_AP+'/kalasanty/predicted_pockets/*')
pds = [pps.split('/')[-1] for pps in glob.glob(c2+'/*')]

pdb
for facts in [(10,1),(12,1)]:
	rns = rn7
	factory = ChemicalFeatures.BuildFeatureFactory('/spstorage/USERS/gina/Project/AP/etc/BaseFeatures'+f_list[0])
	fn = factory.GetFeatureFamilies()
	bns = 30/facts[0]
	cns = [bns*g for g in list(range(facts[0]))]
	st = 30/(inv**facts[0]) # 16 bin
	cnrs = [st*(inv**g) for g in list(range(facts[0]))]
	def bulks_pdb(pps):
		try:
			ppm = Chem.rdmolfiles.MolFromMol2File(DIR_AP+'/kalasanty/predicted_pockets/'+pps+'/pocket0.mol2', removeHs=True)
			ppf = factory.GetFeaturesForMol(ppm)
			ppa = ppm.GetAtoms()
			feats = filter_feats_renew(ppf)
			featss = []
			for f in feats:
				if len(f.GetAtomIds())==1:
					pa =  ppa[f.GetAtomIds()[0]]
					if pa.GetPropsAsDict()['_TriposAtomName'] not in ['N','O','CA']:
						if pa.GetPropsAsDict()['_TriposAtomName'][0]=='C':
							if len([aa for aa in pa.GetNeighbors() if aa.GetSymbol() in ['O','N']])==0:
								featss.append(f)
						else:
							featss.append(f)
				else:
					featss.append(f)
			das = [f.GetAtomIds() for f in featss if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')]
			das = unique([x for x in das if das.count(x) > 1])
			featss = filter_feats_both(featss)
			vt = pairLr_dist(featss, rns, cnrs, inv, das, fn)
		except:
			vt = [0]*len(rns)*len(cns)
		return tuple(vt)
	aps = mclapply(pds, bulks_pdb, ncpu=50)
	app = pd.DataFrame(aps)
	app.index= pds
	app.to_pickle('{}/final2/dataset/pdb/binding/pockets/predidcted_pockets_{}_{}.pkl'.format(DIR_AP, facts[0], facts[1]))
	del app



c2 =DIR_AP+'/DB/refined-set'

pds1=[pps.split('/')[-1] for pps in glob.glob(c2+'/*')]

c2 = DIR_AP+'/DB/v2020-other-PL'
pds = pds1 + [pps.split('/')[-1] for pps in glob.glob(c2+'/*')]



for facts in [(10,1),(12,1)]:
	rns = [str(x)+str(y) for x,y in combinations_with_replacement([0,1,2,3,5,6,9],2)]
	factory = ChemicalFeatures.BuildFeatureFactory('/spstorage/USERS/gina/Project/AP/etc/BaseFeatures'+f_list[0])
	fn = factory.GetFeatureFamilies()
	bns = 30/facts[0]
	cns = [bns*g for g in list(range(facts[0]))]
	st = 30/(inv**facts[0]) # 16 bin
	cnrs = [st*(inv**g) for g in list(range(facts[0]))]
	def bulks_pdb(pps):
		try:
			if pps not in pds1:
				ppm = Chem.rdmolfiles.MolFromPDBFile(DIR_AP+'/DB/v2019-other-PL/'+pps+'/'+pps+'_pocket.pdb', removeHs=True)
			else:
				ppm = Chem.rdmolfiles.MolFromPDBFile(DIR_AP+'/DB/refined-set/'+pps+'/'+pps+'_pocket.pdb', removeHs=True)
			ppf = factory.GetFeaturesForMol(ppm)
			ppa = ppm.GetAtoms()
			feats = filter_feats_renew(ppf)
			featss = []
			for f in feats:
				if len(f.GetAtomIds())==1:
					pa =  ppa[f.GetAtomIds()[0]]
					pi = pa.GetPDBResidueInfo()
					if pi.GetName().split(' ')[1] not in ['N','O','CA']:
						if pa.GetSymbol()=='C':
							if len([aa for aa in pa.GetNeighbors() if aa.GetSymbol() in ['O','N']])==0:
								featss.append(f)
						else:
							featss.append(f)
				else:
					featss.append(f)
			das = [f.GetAtomIds() for f in featss if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')]
			das = unique([x for x in das if das.count(x) > 1])
			featss = filter_feats_both(featss)
			vt = pairLr_dist(featss, rns, cnrs, inv, das, fn)
		except:
			vt = [0]*len(rns)*len(cns)
		return tuple(vt)
	aps = []
	for i in range(11):
		aps.extend(mclapply(pds[i*1610:(i+1)*1610], bulks_pdb, ncpu=75))
	#aps = mclapply(pds, bulks_pdb, ncpu=75)
	app = pd.DataFrame(aps)
	app.index= pds
	app.to_pickle('{}/final2/dataset/pdb/binding/pockets/pdb_pockets_{}_{}.pkl'.format(DIR_AP, facts[0], facts[1]))
	del app

pock = pd.read_pickle('{}/final2/dataset/pdb/binding/pockets/pdb_pockets_{}_{}.pkl'.format(DIR_AP, facts[0], facts[1]))

pdbs = glob.glob(DIR_AP+'/kalasanty/alphafold_pocket/*.pdb')
for facts in [(10,1),(12,1)]:
	rns = [str(x)+str(y) for x,y in combinations_with_replacement([0,1,2,3,5,6,9],2)]
	factory = ChemicalFeatures.BuildFeatureFactory('/spstorage/USERS/gina/Project/AP/etc/BaseFeatures'+f_list[0])
	fn = factory.GetFeatureFamilies()
	bns = 30/facts[0]
	cns = [bns*g for g in list(range(facts[0]))]
	st = 30/(inv**facts[0]) # 16 bin
	cnrs = [st*(inv**g) for g in list(range(facts[0]))]
	def bulks_pdb(pps):
		try:
			ppm = Chem.rdmolfiles.MolFromPDBFile(pps, removeHs=True)
			ppf = factory.GetFeaturesForMol(ppm)
			ppa = ppm.GetAtoms()
			feats = filter_feats_renew(ppf)
			featss = []
			for f in feats:
				if len(f.GetAtomIds())==1:
					pa =  ppa[f.GetAtomIds()[0]]
					pi = pa.GetPDBResidueInfo()
					if pi.GetName().split(' ')[1] not in ['N','O','CA']:
						if pa.GetSymbol()=='C':
							if len([aa for aa in pa.GetNeighbors() if aa.GetSymbol() in ['O','N']])==0:
								featss.append(f)
						else:
							featss.append(f)
				else:
					featss.append(f)
			das = [f.GetAtomIds() for f in featss if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')]
			das = unique([x for x in das if das.count(x) > 1])
			featss = filter_feats_both(featss)
			vt = pairLr_dist(featss, rns, cnrs, inv, das, fn)
		except:
			vt = [0]*len(rns)*len(cns)
		return tuple(vt)
	aps = mclapply(pdbs, bulks_pdb, ncpu=5)
	#aps = mclapply(pds, bulks_pdb, ncpu=75)
	app = pd.DataFrame(aps)
	app.index= [pds.split('/')[-1].split('.')[0] for pds in pdbs]
	app.to_pickle('{}/final2/dataset/pdb/binding/pockets/alpha2_pockets_{}_{}.pkl'.format(DIR_AP, facts[0], facts[1]))
	del app

app = pd.read_pickle('{}/final2/dataset/pdb/binding/pockets/pdb_pockets_{}_{}.pkl'.format(DIR_AP, facts[0], facts[1]))
	
augipdb2cid = pd.merge(augipdb2cmp, cmp2cid).drop_duplicates()
apdb2cid = pd.merge(apdb2cmp.loc[:,['PDB','cmp']], cmp2cid).drop_duplicates()
apdb2cid.to_csv(DIR_AP+'/final2/apdb2cid.csv')


### dataset ver.


# 1 random 
x = np.array(['train', 'val', 'test'])
augipdb2cid1 = augipdb2cid.sample(frac=1,random_state=7)
augipdb2cid1['tag'] = np.repeat(x, [int(len(augipdb2cid)*0.7),int(len(augipdb2cid)*0.15),len(augipdb2cid)-(int(len(augipdb2cid)*0.7)+int(len(augipdb2cid)*0.15))], axis=0)
augipdb2cid1 = augipdb2cid1.sample(frac=1,random_state=7)
augipdb2cid1['label'] = 0

apdb2cid1 = apdb2cid.sample(frac=1,random_state=7)
apdb2cid1['tag'] = np.repeat(x, [int(len(apdb2cid1)*0.7),int(len(apdb2cid1)*0.15),len(apdb2cid1)-(int(len(apdb2cid1)*0.7)+int(len(apdb2cid1)*0.15))], axis=0)
apdb2cid1 = apdb2cid1.sample(frac=1,random_state=7)
apdb2cid1['label'] = 1

sets11 = pd.concat([augipdb2cid1,apdb2cid1])
sets11 = sets11.sample(frac=1,random_state=7)

sets1.to_csv(DIR_AP+'/final2/ptncid1.csv')

app = pd.read_pickle('{}/final2/dataset/pdb/binding/pockets/pockets_{}_{}_{}_{}.pkl'.format(DIR_AP, facts[0], facts[1], facts[2], facts[3]))

# 2 even ratio for aug & target

augipdb2cid = pd.read_csv(DIR_AP+'/final2/augipdb2cid.csv', index_col=0)
apdb2cid = pd.read_csv(DIR_AP+'/final2/apdb2cid.csv',  index_col=0)
sets1 = pd.read_csv(DIR_AP+'/final2/ptncid1.csv',  index_col=0)


augipdb2cid2 = pd.merge(df2etid, apdb2cid)
apdb2cid2 = pd.merge(df2etid, apdb2cid)
for et in unique(augipdb2cid2.etid):
	augipdb2cid2[augipdb2cid2.etid==et]

# 3 non seen target
dsts = pd.read_csv(DIR_AP+'/final2/tree_result')
with open(main+'/tree_result', "rt") as f:
	dsts= f.readlines()

dsts = [d for d in dsts if '|' in d]
dsts = [(d.split('|')[1], float(d.split(':')[1].split('\n')[0])) for d in dsts]
dsts= pd.DataFrame(dsts)
dsts.columns=['UniProtKB','sims']
dsts=dsts.sort_values('sims', ascending=False)
dfsets = pd.merge(df2etid, sets11)
dfsets = pd.merge(dfsets, dsts)
dfsets[dfsets.UniProtKB.isin(dsts[dsts.sims>0.40].UniProtKB)]

list(dfsets[dfsets.UniProtKB.isin(dsts[dsts.sims>0.22].UniProtKB)].label).count(1)
list(dfsets[dfsets.UniProtKB.isin(dsts[dsts.sims>0.22].UniProtKB)].label).count(0)
list(dfsets[~dfsets.UniProtKB.isin(dsts[dsts.sims>0.07].UniProtKB)].label).count(1)
list(dfsets[~dfsets.UniProtKB.isin(dsts[dsts.sims>0.07].UniProtKB)].label).count(0)

dfsets.tag='test'
dfsets.loc[dfsets.sims>0.07,'tag']='val'
dfsets.loc[dfsets.sims>0.22,'tag']='train'

dfsets[dfsets.tag=='train']
dfsets[dfsets.tag=='val']

dfsets[dfsets.tag=='test']

sets3 = dfsets.loc[:,['PDB','cmp','cid','tag','label']]
sets3 = sets3.dropna()
sets3 = sets3.drop_duplicates()
sets3.to_csv(DIR_AP+'/final2/ptncid3.csv')


# filtered
cid_info = pd.read_csv(main+'/cid_info.csv')

cid_info = cid_info.loc[:,['cid', 'mw', 'rotbonds', 'hbonddonor','hbondacc','cidcdate', 'heavycnt', 'xlogp']]

subs = cid_info[cid_info.xlogp<=5]
subs = subs[subs.mw<=500]
subs = subs[subs.hbonddonor<=5]
subs = subs[subs.hbondacc<=10]
subs = subs[subs.mw>=180]
subs = subs[subs.rotbonds<=10]

fsets3 = sets3[sets3.cid.isin(subs.cid)]
fsets3[(fsets3.tag=='train') & (fsets3.label==1)]
fsets3[(fsets3.tag=='test') & (fsets3.label==1)]
fsets3[(fsets3.tag=='val') & (fsets3.label==1)]
fsets3.to_csv(DIR_AP+'/final2/fptncid3.csv')

## kinase

kk=pd.read_csv(main+'/kinase_screen.csv')
kc = pd.read_csv(main+'/kinase.csv')
pt = re.compile('[0-9]+-[0-9]{2}-[0-9]')

kkt = kk.T
kkt[kkt.iloc[:,1].isin(set(kk.iloc[1,:])-set(kc[(kc.cas!='') & (kc.cas.isin(kk.iloc[1,:]))].cas))]

for cs in list(kc.cmpdsynonym):
	re.search(pt, cs)[0]


kc['cas'] = [re.search(pt, cs)[0] if not re.search(pt, cs)==None else '' for cs in list(kc.cmpdsynonym)]
kc[(kc.cas!='') & (kc.cas.isin(kk.iloc[1,:]))]


with open(main+'/kinase_cids.txt', "wt") as f:f.write('\n'.join([str(int(cd)) for cd in list(kc.cid)]))

ksdf = Chem.SDMolSupplier(main+'/kinase.sdf', removeHs=False) 
kcds = [m.GetPropsAsDict()['PUBCHEM_COMPOUND_CID'] for m in ksdf]
kc[kc.cid.isin(kcds)]
# 2051 / 5113385 / 3797(176406) / 5339183 / 141246은 3D x / 249762-74-1은 삭제 

kns = list(itertools.chain(*[kn.split('/') for kn in list(kk.iloc[2:,0])]))
knsh = [kn + '_HUMAN' for kn in kns]
len(set(df2etid[df2etid['Entry name'].isin(knsh)]['Entry name'])) # 86 

pdbgns = list(itertools.chain(*[gn.split(' ') for gn in unique(df2etid['Gene names']) if not isinstance(gn, float)]))
kk[kk.iloc[:,0].isin(pdbgns)] # 90
df2etid[df2etid['Entry name'].isin(knsh)] # 86

list(kk[~kk.iloc[:,0].isin(pdbgns)].iloc[:,0])
kkns = list(itertools.chain(*[kn.split('/') for kn in list(kk[~kk.iloc[:,0].isin(pdbgns)].iloc[2:,0])]))
kknsh = [kn + '_HUMAN' for kn in kkns]
len(set(df2etid[df2etid['Entry name'].isin(kknsh)]['Entry name'])) # 21

df2etid[df2etid['Entry name'].isin(knsh)]


nknsh = [kn for kn in knsh if kn not in set(df2etid[df2etid['Entry name'].isin(knsh)]['Entry name'])]

# five conformer
for cd in list(sets3.cid):


# all conformer


for et in unique(dfsets.etid):
	augipdb2cid3[augipdb2cid3.etid==et]

cc=' '.join(unique(dfsets.UniProtKB))

with open(main+'/dfsets_unis.txt', "wt") as f:f.write(cc)

# 4 regression



### pocket from kalansky

llm = Chem.rdmolfiles.MolFromMol2File(DIR_AP+'/kalasanty/predicted_pockets/'+poc+'/pocket0.mol2', removeHs=True)


lls = llm.GetConformer().GetPositions()
# lls[8] = [22.723,   4.736,  34.281]
llf = factory.GetFeaturesForMol(llm)
lluse_pos = np.array([ft.GetPos() for ft in llf])
lluse_fam = [fn.index(ft.GetFamily()) for ft in llf]
llb = llm.GetBonds()
lla = llm.GetAtoms()



def filter_feats_renew(feat):
	arh = [f.GetAtomIds() for f in feat if (f.GetFamily()=='Aromatic')]
	hydp = [f.GetAtomIds() for f in feat if (f.GetFamily()=='Hydrophobe')]
	ions = [f.GetAtomIds() for f in feat if (f.GetFamily()=='NegIonizable')|(f.GetFamily()=='PosIonizable')]
	ions = list(itertools.chain(*ions))
	feat = [f for f in feat if not (f.GetFamily()=='Hydrophobe')&(f.GetAtomIds() in arh)]
	feat = [f for f in feat if not (f.GetFamily()=='Donor')&(len(intersect(set(f.GetAtomIds()),ions))!=0)]
	feat = [f for f in feat if not (f.GetFamily()=='Acceptor')&(len(intersect(set(f.GetAtomIds()),ions))!=0)]
	return feat


llff = filter_feats_renew(llf)


# 2 version. 1 = neutral poral not filtered / 2 = filtered
# 3 bound property cut 4 hydrophobe cut
llfff = []
for f in llff:
	if len(f.GetAtomIds())==1:
		la =  lla[f.GetAtomIds()[0]]
		if la.GetPropsAsDict()['_TriposAtomName'] not in ['N','O','CA']:
			if la.GetPropsAsDict()['_TriposAtomName'][0]=='C':
				if len([aa for aa in la.GetNeighbors() if aa.GetSymbol() in ['O','N']])==0:
					llfff.append(f)
			else:
				llfff.append(f)
	else:
		llfff.append(f)




[f for f in llff if lla[f.GetAtomIds()]]

close_feats = [f for i,f in enumerate(ppf) if min([dst([ppuse_pos[i],lgd_po]) for lgd_po in lls])<4.5]
# close_feats = [f for i,f in enumerate(ppf) if min([dst([ppuse_pos[i],lgd_po]) for lgd_po in lluse_pos])<4]
close_atoms = [f.GetAtomIds() for i,f in enumerate(ppf) if min([dst([ppuse_pos[i],lgd_po]) for lgd_po in lls])<4.5]
close_atoms = list(itertools.chain(*close_atoms))

pps, facts, rns, cnrs, fn = ssets
ppm = Chem.MolFromPDBFile(f'{dtip}/covid/{pps}_protein.pdb', removeHs=True)
llm = Chem.MolFromPDBFile(f'{dtip}/covid/{pps}_ligand.pdb', removeHs=True)
mm = Chem.SDMolSupplier(sdfs, removeHs=True) 
ppf = factory.GetFeaturesForMol(ppm)
ppuse_pos = np.array([ft.GetPos() for ft in ppf])
subps = ori[ori.iloc[:,1]==pps]
vts = []
for ilen in range(len(subps)):
	ldst, lx,ly,lz,x_min, x_max,y_min, y_max, z_min, z_max = subps.iloc[ilen, 2:12]
	feat = []
	for i, f in enumerate(ppf):
		if x_min < ppuse_pos[i][0] < x_max and y_min < ppuse_pos[i][1] < y_max and z_min < ppuse_pos[i][2] < z_max:
			feat.append(f)
	das = [[]]



from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure(figsize=(10, 5))
ax = fig.add_subplot(111, projection='3d') # Axe3D object
# ax.scatter(1.5277, -10.4733, 1.6311, c = 'k', s= 20, alpha=0.5)#, cmap=plt.cm.Greens)

xposs = [llm.GetConformer().GetAtomPosition(i).x for i in range(llm.GetNumAtoms())]
yposs = [llm.GetConformer().GetAtomPosition(i).y for i in range(llm.GetNumAtoms())]
zposs = [llm.GetConformer().GetAtomPosition(i).z for i in range(llm.GetNumAtoms())]


# xposs = [llm.GetConformer().GetAtomPosition(i).x if llm.GetConformer().GetAtomPosition(i).x!=723 else 22.723 for i in range(llm.GetNumAtoms())]
# yposs = [llm.GetConformer().GetAtomPosition(i).y if llm.GetConformer().GetAtomPosition(i).y!=736 else 4.736 for i in range(llm.GetNumAtoms())]
# zposs = [llm.GetConformer().GetAtomPosition(i).z if llm.GetConformer().GetAtomPosition(i).z!=281 else 34.281 for i in range(llm.GetNumAtoms())]

xatoms = [a.GetSymbol() for a in llm.GetAtoms()]

catoms = ['k' if xa=='C' else 'b' if xa=='N' else 'r' if xa=='O' else 'c' for xa in xatoms]

ax.scatter(xposs, yposs, zposs, c = catoms, s= 20, alpha=0.5)#, cmap=plt.cm.Greens)

lstdic = {1:'solid',2:'dashed',1.5:'dashdot'}
for bc in llb:
	xb = xposs[bc.GetBeginAtomIdx()], xposs[bc.GetEndAtomIdx()] 
	yb = yposs[bc.GetBeginAtomIdx()], yposs[bc.GetEndAtomIdx()] 
	zb = zposs[bc.GetBeginAtomIdx()], zposs[bc.GetEndAtomIdx()] 
	lst = lstdic[bc.GetBondTypeAsDouble()]
	ax.plot(xb, yb, zb, alpha=(bc.GetBondTypeAsDouble())/2, color = 'k', linestyle= lst)# marker='o')

for f in llfff:
	fam = f.GetFamily()
	if (fam!='ZnBinder')&(fam!='LumpedHydrophobe'):
		ax.text(f.GetPos().x,f.GetPos().y,f.GetPos().z, fam, (1,1,1))


for f in llf:
	fam = f.GetFamily()
	ax.text(f.GetPos().x,f.GetPos().y,f.GetPos().z, fam, (1,1,1))

# for i,mf in enumerate(use_mmff):
# 	ax.text(use_idx[i][0]+0.5,use_idx[i][1]+0.5,use_idx[i][2]+0.5, mf, (1,1,1))

# xposs = [ft.GetPos().x for ft in close_feats]
# yposs = [ft.GetPos().y for ft in close_feats]
# zposs = [ft.GetPos().z for ft in close_feats]

ppa = [a.GetIdx() for a in ppm.GetAtoms()]

xposs = [ppm.GetConformer().GetAtomPosition(i).x for i in ppa if i in close_atoms]
yposs = [ppm.GetConformer().GetAtomPosition(i).y for i in ppa if i in close_atoms]
zposs = [ppm.GetConformer().GetAtomPosition(i).z for i in ppa if i in close_atoms]

xatoms = [a.GetSymbol() for a in ppm.GetAtoms() if a.GetIdx() in close_atoms]

# xposs = [ppm.GetConformer().GetAtomPosition(i).x for i in ppa]
# yposs = [ppm.GetConformer().GetAtomPosition(i).y for i in ppa]
# zposs = [ppm.GetConformer().GetAtomPosition(i).z for i in ppa]

# xatoms = [a.GetSymbol() for a in ppm.GetAtoms()]

ppb = ppm.GetBonds()

catoms = ['k' if xa=='C' else 'b' if xa=='N' else 'r' if xa=='O' else 'c' for xa in xatoms]

ax.scatter(xposs, yposs, zposs, c = catoms, s= 20, alpha=0.5)#, cmap=plt.cm.Greens)

xposs = [ppm.GetConformer().GetAtomPosition(i).x for i in range(ppm.GetNumAtoms())]
yposs = [ppm.GetConformer().GetAtomPosition(i).y for i in range(ppm.GetNumAtoms())]
zposs = [ppm.GetConformer().GetAtomPosition(i).z for i in range(ppm.GetNumAtoms())]


lstdic = {1:'solid',2:'dashed',1.5:'dashdot'}
for bc in ppb:
	if (bc.GetBeginAtomIdx() in close_atoms)|(bc.GetEndAtomIdx() in close_atoms):
		xb = xposs[bc.GetBeginAtomIdx()], xposs[bc.GetEndAtomIdx()] 
		yb = yposs[bc.GetBeginAtomIdx()], yposs[bc.GetEndAtomIdx()] 
		zb = zposs[bc.GetBeginAtomIdx()], zposs[bc.GetEndAtomIdx()] 
		if bc.GetBondTypeAsDouble()!=0:
			lst = lstdic[bc.GetBondTypeAsDouble()]
			ax.plot(xb, yb, zb, alpha=(bc.GetBondTypeAsDouble())/2, color = 'limegreen', linestyle= lst)# marker='o')

for f in close_feats:
	fam = f.GetFamily()
	if fam not in ['ZnBinder','LumpedHydrophobe']:
		ax.text(f.GetPos().x,f.GetPos().y,f.GetPos().z, fam, (1,1,1), color='red')

plt.show()

###
poc='3s1y_S1Y_398'
llm = Chem.SDMolSupplier(DIR_AP+'/DB/refined-set/'+poc+'/'+poc+'_ligand.sdf', removeHs=True)[0]
llm = Chem.rdmolfiles.MolFromMol2File(DIR_AP+'/DB/refined-set/'+pps+'/'+pps+'_ligand.mol2', removeHs=True)
llm = Chem.rdmolfiles.MolFromPDBFile(DIR_AP+'/DB/betalactam/pdb/'+poc+'.pdb', removeHs=True)
llm = Chem.rdmolfiles.MolFromPDBFile(DIR_AP+'/DB/v2019-other-PL/'+poc+'/'+poc+'2_ligand.pdb', removeHs=True)
llm = Chem.rdmolfiles.MolFromMol2File(DIR_AP+'/DB/v2019-other-PL/'+poc+'/'+poc+'_ligand.mol2', removeHs=True)
llm = Chem.rdmolfiles.MolFromMol2File(DIR_AP+'/kalasanty/predicted_pockets/'+poc+'/pocket0.mol2', removeHs=True)


lls = llm.GetConformer().GetPositions()
# lls[8] = [22.723,   4.736,  34.281]
llf = factory.GetFeaturesForMol(llm)
lluse_pos = np.array([ft.GetPos() for ft in llf])
lluse_fam = [fn.index(ft.GetFamily()) for ft in llf]
llb = llm.GetBonds()


ppm = Chem.rdmolfiles.MolFromMol2File(DIR_AP+'/kalasanty/alphafold_pocket/'+'P25094_1'+'/pocket0.mol2', removeHs=True)

ppm = Chem.rdmolfiles.MolFromPDBFile(DIR_AP+'/DB/refined-set/'+poc+'/'+poc+'_pocket.pdb', removeHs=True)
ppm = Chem.rdmolfiles.MolFromPDBFile(DIR_AP+'/DB/v2019-other-PL/'+poc+'/'+poc+'_pocket.pdb', removeHs=True)
ppm = Chem.rdmolfiles.MolFromPDBFile(DIR_AP+'/DB/betalactam/pdb/'+poc[0:4]+'_protein.pdb', removeHs=True)
pps = ppm.GetConformer().GetPositions()
ppf = factory.GetFeaturesForMol(ppm)
ppuse_pos = np.array([ft.GetPos() for ft in ppf])
ppuse_fam = [fn.index(ft.GetFamily()) for ft in ppf]
close_feats = [f for i,f in enumerate(ppf) if min([dst([ppuse_pos[i],lgd_po]) for lgd_po in lls])<4.5]
# close_feats = [f for i,f in enumerate(ppf) if min([dst([ppuse_pos[i],lgd_po]) for lgd_po in lluse_pos])<4]
close_atoms = [f.GetAtomIds() for i,f in enumerate(ppf) if min([dst([ppuse_pos[i],lgd_po]) for lgd_po in lls])<4.5]
close_atoms = list(itertools.chain(*close_atoms))


fig = plt.figure(figsize=(10, 5))
ax = fig.add_subplot(111, projection='3d') # Axe3D object

xposs = [ppm.GetConformer().GetAtomPosition(i).x for i in range(ppm.GetNumAtoms())]
yposs = [ppm.GetConformer().GetAtomPosition(i).y for i in range(ppm.GetNumAtoms())]
zposs = [ppm.GetConformer().GetAtomPosition(i).z for i in range(ppm.GetNumAtoms())]


catoms = ['k' if xa=='C' else 'b' if xa=='N' else 'r' if xa=='O' else 'c' for xa in xatoms]

ax.scatter(xposs, yposs, zposs, c = catoms, s= 20, alpha=0.5)#, cmap=plt.cm.Greens)


lstdic = {1:'solid',2:'dashed',1.5:'dashdot'}
for bc in ppb:
	xb = xposs[bc.GetBeginAtomIdx()], xposs[bc.GetEndAtomIdx()] 
	yb = yposs[bc.GetBeginAtomIdx()], yposs[bc.GetEndAtomIdx()] 
	zb = zposs[bc.GetBeginAtomIdx()], zposs[bc.GetEndAtomIdx()] 
	if bc.GetBondTypeAsDouble()!=0:
		lst = lstdic[bc.GetBondTypeAsDouble()]
		ax.plot(xb, yb, zb, alpha=(bc.GetBondTypeAsDouble())/2, color = 'limegreen', linestyle= lst)# marker='o')

for f in feats:
	fam = f.GetFamily()
	if fam not in ['ZnBinder','LumpedHydrophobe']:
		ax.text(f.GetPos().x,f.GetPos().y,f.GetPos().z, fam, (1,1,1), color='black')

plt.show()

### 정리



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
	hydp = [f.GetAtomIds() for f in feat if (f.GetFamily()=='Hydrophobe')]
	ions = [f.GetAtomIds() for f in feat if (f.GetFamily()=='NegIonizable')|(f.GetFamily()=='PosIonizable')]
	ions = list(itertools.chain(*ions))
	feat = [f for f in feat if not (f.GetFamily()=='Hydrophobe')&(f.GetAtomIds() in arh)]
	feat = [f for f in feat if not (f.GetFamily()=='Donor')&(len(intersect(set(f.GetAtomIds()),ions))!=0)]
	feat = [f for f in feat if not (f.GetFamily()=='Acceptor')&(len(intersect(set(f.GetAtomIds()),ions))!=0)]
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




def feature_PP(row,feats, mmffp, mps, idx, col, gap, inv, das, fn):
	use_pos = np.array([ft.GetPos() for ft in feats[idx]])
	use_fam = [fn.index(ft.GetFamily()) if not ft.GetAtomIds() in das[idx] else 9 for ft in feats[idx]]
	mmffps = [m for m in mmffp[idx][1::] if len(mps[idx])>=int(m.split(' ')[0])]
	use_idx = np.array([mps[idx][int(m.split(' ')[0])-1] for m in mmffps if not bool(mps[idx][int(m.split(' ')[0])-1] in use_pos)])# not in f_idx]# if int(m.split(' ')[0])<mhac[idx]]
	use_mmff = [int(float(m.split(' ')[1])//gap + 12) for m in mmffps if not bool(mps[idx][int(m.split(' ')[0])-1] in use_pos)]# if int(m.split(' ')[0])<mhac[idx]]
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
	use_fam = [fn.index(ft.GetFamily()) if not ft.GetAtomIds() in das[idx] else 9 for ft in feats[idx]]
	if (len(uu_idx) != 0) & (len(use_pos) !=0):
		use_idx = np.array(itemgetter(*uu_idx)(mps[idx]))
		if len(uu_idx)>1:
			use_atms = [sym2feat[at.GetSymbol()] if at.GetSymbol() in ['C','O','N','S'] else 14 for at in itemgetter(*uu_idx)(atms[idx])]
			poses = np.concatenate((use_pos,use_idx))
		else:
			atsym = itemgetter(*uu_idx)(atms[idx]).GetSymbol()
			poses = np.concatenate((use_pos,np.array([use_idx])))
			if atsym in ['C','O','N','S']:
				use_atms = [sym2feat[atsym]]
			else:
				use_atms = [14]
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


sym2feat2 = {'C':7,'O':10,'N':10, 'S':10}
def feature_AH(row, feats, atms, mps, idx, col, inv, das, fn):
	# basic charac. extra atom to atom charac. 7 C hydrophobic / 10 O  N  S hydrophilic / 11 etc
	ft_idx = unique(list(itertools.chain(*[f.GetAtomIds() for f in feats[idx]])))
	uu_idx = [i for i in range(len(mps[idx])) if i not in ft_idx]
	use_pos = np.array([ft.GetPos() for ft in feats[idx]])
	use_fam = [fn.index(ft.GetFamily()) if not ft.GetAtomIds() in das[idx] else 9 for ft in feats[idx]]
	if (len(uu_idx) != 0) & (len(use_pos) !=0):
		use_idx = np.array(itemgetter(*uu_idx)(mps[idx]))
		if len(uu_idx)>1:
			use_atms = [sym2feat2[at.GetSymbol()] if at.GetSymbol() in ['C','O','N','S'] else 11 for at in itemgetter(*uu_idx)(atms[idx])]
			poses = np.concatenate((use_pos,use_idx))
		else:
			atsym = itemgetter(*uu_idx)(atms[idx]).GetSymbol()
			poses = np.concatenate((use_pos,np.array([use_idx])))
			if atsym in ['C','O','N','S']:
				use_atms = [sym2feat2[atsym]]
			else:
				use_atms = [11]
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



def feature_AP(row, feats, atms, mmffp, gap, mps, idx, col, inv, das, fn):
	# basic charac. extra atom to atom charac. 10 C / 11 O / 12 N / 13 S / 14 etc
	ft_idx = unique(list(itertools.chain(*[f.GetAtomIds() for f in feats[idx]])))
	uu_idx = [i for i in range(len(mps[idx])) if i not in ft_idx]
	use_pos = np.array([ft.GetPos() for ft in feats[idx]])
	use_fam = [fn.index(ft.GetFamily()) if not ft.GetAtomIds() in das[idx] else 9 for ft in feats[idx]]
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
	use_mmff = [int(float(m.split(' ')[1])//gap + 16) for m in mmffps if not bool(mps[idx][int(m.split(' ')[0])-1] in poses)]# if int(m.split(' ')[0])<mhac[idx]]
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


c2 =DIR_AP+'/DB/refined-set'

pds1=[pps.split('/')[-1] for pps in glob.glob(c2+'/*')]

c2 = DIR_AP+'/DB/v2020-other-PL'
pds = pds1 + [pps.split('/')[-1] for pps in glob.glob(c2+'/*')]

dstbin = 10
f_list=['','','','_all'] +['_h']*8
st = 30/(inv**dstbin) # 16 bin
cnrs = [st*(inv**g) for g in range(dstbin)]
fnss = [[0,1,2,3,5,6],[0,1,2,3,5,6],[0,1,2,3,5,6,9],[0,1,2,3,4,5,7,8,10],[0,1,2,3,4,6,7,9],[0,1,2,3,4,6,7,9,10,11,12,13,14],[0,1,2,3,4,6,7,9,10,11,12,13,14],[0,1,2,3,4,6,7,9,10,11,12,13,14,15,16,17,18],[],[0,1,2,3,4,6,7,9,10,11,12,13,14],[],[0,1,2,3,4,6,7,9,10,11]]

vidirs = glob.glob(DIR_AP+'/final3/dude/vina_docked/*')

for dd in [4,5]:
	for ii in [11]:
		facts = [ii,dd]
		factory = ChemicalFeatures.BuildFeatureFactory('{}/etc/BaseFeatures{}.fdef'.format(DIR_AP,f_list[facts[0]]))
		fn = factory.GetFeatureFamilies()
		fns = fnss[facts[0]]
		rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
		def bulks_pdb(ssets):
			try:
				pps, facts, rns, cnrs, fn = ssets
				if pps in pds1:didx = 'refined-set'
				else:didx = 'v2020-other-PL'
				ppm = Chem.rdmolfiles.MolFromPDBFile('{}/DB/{}/{}/{}_pocket.pdb'.format(DIR_AP,didx,pps,pps), removeHs=True)
				llm = Chem.SDMolSupplier('{}/DB/{}/{}/{}_ligand.sdf'.format(DIR_AP,didx,pps,pps), removeHs=True)[0]
				if llm is None:
					llm = Chem.rdmolfiles.MolFromMol2File('{}/DB/{}/{}/{}_ligand.mol2'.format(DIR_AP,didx,pps,pps), removeHs=True)
				llf = factory.GetFeaturesForMol(llm)
				lluse_pos = np.array([ft.GetPos() for ft in llf])
				ppf = factory.GetFeaturesForMol(ppm)
				ppuse_pos = np.array([ft.GetPos() for ft in ppf])
				feat = [f for i,f in enumerate(ppf) if min([dst([ppuse_pos[i],lgd_po]) for lgd_po in lluse_pos])<facts[1]]
				das = [[]]
				if facts[0] > 0:
					feat = filter_feats_renew(feat)
					if facts[0] > 1:
						das = [f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')]
						das = unique([x for x in das if das.count(x) > 1])
						if facts[0] > 4:
							mps = ppm.GetConformer().GetPositions()							
							atms = [[am for i,am in enumerate(ppm.GetAtoms()) if min([dst([mps[i],lgd_po]) for lgd_po in lluse_pos])<facts[1]]]
							mps = [pm for pm in ppm.GetConformer().GetPositions() if min([dst([pm,lgd_po]) for lgd_po in lluse_pos])<facts[1]]
							mps = [np.stack(mps)]
							if facts[0] ==6:
								vt = feature_AA(rns, [feat], atms, mps, 0, cnrs, inv, [das], fn)
							else:
								vt = feature_AH(rns, [feat], atms, mps, 0, cnrs, inv, [das], fn)
						else:
							vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
					else:
						vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
				else:
					vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
			except:
				vt = [0]*len(rns)*len(cnrs)
			return tuple(vt)
		aps = mclapply([[pps, facts, rns, cnrs, fn] for pps in pds], bulks_pdb, ncpu=100)
		app = pd.DataFrame(aps)
		app.index= pds
		app.to_pickle('{}/final3/pdbbind/dataset/pockets/pockets_{}_{}.pkl'.format(DIR_AP, facts[0], facts[1]))
		del app


for ii in [0,6,11]:
	facts = [ii]
	factory = ChemicalFeatures.BuildFeatureFactory('{}/etc/BaseFeatures{}.fdef'.format(DIR_AP,f_list[facts[0]]))
	fn = factory.GetFeatureFamilies()
	fns = fnss[facts[0]]
	rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
	def bulks_pdb(ssets):
		try:
			pps, facts, rns, cnrs, fn = ssets
			ppm = Chem.MolFromPDBFile('{}/PDB/{}_{}.pdb'.format(dtip,pps[0:4].upper(),pps[4].upper()), removeHs=True)
			ppf = factory.GetFeaturesForMol(ppm)
			ppuse_pos = np.array([ft.GetPos() for ft in ppf])
			subps = ori[ori.iloc[:,1]==pps]
			vts = []
			for ilen in range(len(subps)):
				x_min, x_max,y_min, y_max, z_min, z_max= subps.iloc[ilen, 2:8]
				feat = []
				for i, f in enumerate(ppf):
					if x_min < ppuse_pos[i][0] < x_max and y_min < ppuse_pos[i][1] < y_max and z_min < ppuse_pos[i][2] < z_max:
						feat.append(f)
				das = [[]]
				if len(feat)>5:
					if facts[0] > 0:
						feat = filter_feats_renew(feat)
						if facts[0] > 1:
							das = [f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')]
							das = unique([x for x in das if das.count(x) > 1])
							if facts[0] > 4:
								mps = ppm.GetConformer().GetPositions()							
								atms = [[am for i,am in enumerate(ppm.GetAtoms()) if x_min < mps[i][0] < x_max and y_min < mps[i][1] < y_max and z_min < mps[i][2] < z_max]]
								mps = [pm for pm in ppm.GetConformer().GetPositions() if x_min < pm[0] < x_max and y_min < pm[1] < y_max and z_min < pm[2] < z_max]
								mps = [np.stack(mps)]
								if facts[0] ==6:
									vt = feature_AA(rns, [feat], atms, mps, 0, cnrs, inv, [das], fn)
								else:
									vt = feature_AH(rns, [feat], atms, mps, 0, cnrs, inv, [das], fn)
							else:
								vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
						else:
							vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
					else:
						vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
					vts.append(vt)
			dfs = pd.DataFrame(vts).drop_duplicates()
			dfs['pdb'] = pps
		except:
			dfs = pd.DataFrame([[0]*len(rns)*len(cnrs)]).drop_duplicates()
			dfs['pdb'] = pps
		return dfs
	aps = mclapply([[pps, facts, rns, cnrs, fn] for pps in set(ori.iloc[:,1])], bulks_pdb, ncpu=104)
	app = pd.concat(aps)
	app.index= list(app.pdb)
	app = app.iloc[:,0:-1]
	app.to_pickle('{}/ori_pockets_{}.pkl'.format(dtip, facts[0]))
	del app

dd=4
for ii in [11]:
	facts = [ii,dd]
	factory = ChemicalFeatures.BuildFeatureFactory('{}/etc/BaseFeatures{}.fdef'.format(DIR_AP,f_list[facts[0]]))
	fn = factory.GetFeatureFamilies()
	fns = fnss[facts[0]]
	rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
	def bulks_pdb(ssets):
		try:
			pps, facts, rns, cnrs, fn = ssets
			if pps in pds1:didx = 'refined-set'
			else:didx = 'v2020-other-PL'
			ppm = Chem.rdmolfiles.MolFromPDBFile('{}/DB/{}/{}/{}_pocket.pdb'.format(DIR_AP,didx,pps,pps), removeHs=True)
			llm = Chem.SDMolSupplier('{}/DB/{}/{}/{}_ligand.sdf'.format(DIR_AP,didx,pps,pps), removeHs=True)[0]
			if llm is None:
				llm = Chem.rdmolfiles.MolFromMol2File('{}/DB/{}/{}/{}_ligand.mol2'.format(DIR_AP,didx,pps,pps), removeHs=True)
			llf = factory.GetFeaturesForMol(llm)
			lluse_pos = np.array([ft.GetPos() for ft in llf])
			ppa = ppm.GetAtoms()
			ppf = factory.GetFeaturesForMol(ppm)
			ppuse_pos = np.array([ft.GetPos() for ft in ppf])
			feat = [f for i,f in enumerate(ppf) if min([dst([ppuse_pos[i],lgd_po]) for lgd_po in lluse_pos])<facts[1]]
			das = [[]]
			feat = filter_feats_renew(feat)
			das = [f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')]
			das = unique([x for x in das if das.count(x) > 1])
			mps = ppm.GetConformer().GetPositions()							
			atms = [[am for i,am in enumerate(ppm.GetAtoms()) if min([dst([mps[i],lgd_po]) for lgd_po in lluse_pos])<facts[1]]]
			mps = [pm for pm in ppm.GetConformer().GetPositions() if min([dst([pm,lgd_po]) for lgd_po in lluse_pos])<facts[1]]
			mps = [np.stack(mps)]
			vt = feature_AH(rns, [feat], atms, mps, 0, cnrs, inv, [das], fn)
		except:
			vt = [0]*len(rns)*len(cnrs)
		return tuple(vt)
	aps = mclapply([[pps, facts, rns, cnrs, fn] for pps in pds], bulks_pdb, ncpu=50)
	app = pd.DataFrame(aps)
	app.index= pds
	app.to_pickle('{}/final3/pdbbind/dataset/pockets/pdbval_pockets_20_4_2.pkl'.format(DIR_AP))
	del app

for ii in [11]:
	facts = [ii]
	factory = ChemicalFeatures.BuildFeatureFactory('{}/etc/BaseFeatures{}.fdef'.format(DIR_AP,f_list[facts[0]]))
	fn = factory.GetFeatureFamilies()
	fns = fnss[facts[0]]
	rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
	def bulks_pdb(ssets):
		try:
			pps, facts, rns, cnrs, fn = ssets
			if pps in pds1:didx = 'refined-set'
			else:didx = 'v2020-other-PL'
			ppm = Chem.rdmolfiles.MolFromPDBFile('{}/DB/{}/{}/{}_pocket.pdb'.format(DIR_AP,didx,pps,pps), removeHs=True)
			ppf = factory.GetFeaturesForMol(ppm)
			ppa = ppm.GetAtoms()
			feats = [f for i,f in enumerate(ppf)]
			feat = []
			for f in feats:
				if len(f.GetAtomIds())==1:
					pa =  ppa[f.GetAtomIds()[0]]
					pi = pa.GetPDBResidueInfo()
					if pi.GetName().split(' ')[1] not in ['N','O','CA']:
						if pa.GetSymbol()=='C':
							if len([aa for aa in pa.GetNeighbors() if aa.GetSymbol() in ['O','N']])==0:
								feat.append(f)
						else:
							feat.append(f)
				else:
					feat.append(f)
feat = [f for i,f in enumerate(ppf) if min([dst([ppuse_pos[i],lgd_po]) for lgd_po in lluse_pos])<facts[1]]
			das = []
			feat = filter_feats_renew(feat)
			das = [f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')]
			das = unique([x for x in das if das.count(x) > 1])
			mps = [ppm.GetConformer().GetPositions()]						
			atms = [ppm.GetAtoms()]
			mps = [np.stack(mps)]
			vt = feature_AH(rns, [feat], atms, mps, 0, cnrs, inv, [das], fn)
		except:
			vt = [0]*len(rns)*len(cnrs)
		return tuple(vt)
	aps = mclapply([[pps, facts, rns, cnrs, fn] for pps in pds], bulks_pdb, ncpu=100)
	app = pd.DataFrame(aps)
	app.index= pds
	app.to_pickle('{}/final3/pdbbind/dataset/pockets/pdbval.pkl'.format(DIR_AP, facts[0]))
	del app



for ii in [6,11]:
	facts = [ii]
	factory = ChemicalFeatures.BuildFeatureFactory('{}/etc/BaseFeatures{}.fdef'.format(DIR_AP,f_list[facts[0]]))
	fn = factory.GetFeatureFamilies()
	fns = fnss[facts[0]]
	rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
	def bulks_pdb(ssets):
		try:
			pps, facts, rns, cnrs, fn = ssets
			if pps in pds1:didx = 'refined-set'
			else:didx = 'v2020-other-PL'
			ppm = Chem.rdmolfiles.MolFromPDBFile('{}/DB/{}/{}/{}_pocket.pdb'.format(DIR_AP,didx,pps,pps), removeHs=True)
			ppf = factory.GetFeaturesForMol(ppm)
			ppa = ppm.GetAtoms()
			feats = [f for i,f in enumerate(ppf)]
			feat = []
			for f in feats:
				if len(f.GetAtomIds())==1:
					pa =  ppa[f.GetAtomIds()[0]]
					pi = pa.GetPDBResidueInfo()
					if pi.GetName().split(' ')[1] not in ['N','O','CA']:
						if pa.GetSymbol()=='C':
							if len([aa for aa in pa.GetNeighbors() if aa.GetSymbol() in ['O','N']])==0:
								feat.append(f)
						else:
							feat.append(f)
				else:
					feat.append(f)
			das = []
			if facts[0] > 0:
				feat = filter_feats_renew(feat)
				if facts[0] > 1:
					das = [f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')]
					das = unique([x for x in das if das.count(x) > 1])
					if facts[0] > 4:
						mps = [ppm.GetConformer().GetPositions()]						
						atms = [ppm.GetAtoms()]
						mps = [np.stack(mps)]
						if facts[0] ==6:
							vt = feature_AA(rns, [feat], atms, mps, 0, cnrs, inv, [das], fn)
						else:
							vt = feature_AH(rns, [feat], atms, mps, 0, cnrs, inv, [das], fn)
					else:
						vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
				else:
					vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
			else:
				vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
		except:
			vt = [0]*len(rns)*len(cnrs)
		return tuple(vt)
	aps = mclapply([[pps, facts, rns, cnrs, fn] for pps in pds], bulks_pdb, ncpu=100)
	app = pd.DataFrame(aps)
	app.index= pds
	app.to_pickle('{}/final3/pdbbind/dataset/pockets/pdb_pockets_{}.pkl'.format(DIR_AP, facts[0]))
	del app





## alphafold pocket
for ii in [0,1,2,3,4]:
	facts = [ii,dd]
	factory = ChemicalFeatures.BuildFeatureFactory('{}/etc/BaseFeatures{}.fdef'.format(DIR_AP,f_list[facts[0]]))
	fn = factory.GetFeatureFamilies()
	fns = fnss[facts[0]]
	rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
	def bulks_pdb(ssets):
		try:
			pps, facts, rns, cnrs, fn = ssets
			ppm = Chem.rdmolfiles.MolFromPDBFile(pps, removeHs=True)
			ppf = factory.GetFeaturesForMol(ppm)
			ppa = ppm.GetAtoms()
			feats = [f for i,f in enumerate(ppf)]
			feat = []
			for f in feats:
				if len(f.GetAtomIds())==1:
					pa =  ppa[f.GetAtomIds()[0]]
					pi = pa.GetPDBResidueInfo()
					if pi.GetName().split(' ')[1] not in ['N','O','CA']:
						if pa.GetSymbol()=='C':
							if len([aa for aa in pa.GetNeighbors() if aa.GetSymbol() in ['O','N']])==0:
								feat.append(f)
						else:
							feat.append(f)
				else:
					feat.append(f)
			das = []
			if facts[0] > 0:
				feat = filter_feats_renew(feat)
				if facts[0] > 1:
					das = [f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')]
					das = unique([x for x in das if das.count(x) > 1])
					if facts[0] > 4:
						mps = [ppm.GetConformer().GetPositions()]						
						atms = [ppm.GetAtoms()]
						mps = [np.stack(mps)]
						vt = feature_AA(rns, [feat], atms, mps, 0, cnrs, inv, [das], fn)
					else:
						vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
				else:
					vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
			else:
				vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
		except:
			vt = [0]*len(rns)*len(cnrs)
		return tuple(vt)
	aps = mclapply([[pps, facts, rns, cnrs, fn] for pps in pdbs], bulks_pdb, ncpu=len(pdbs))
	app = pd.DataFrame(aps)
	app.index= [pdd.split('/')[-1].split('_')[0] for pdd in pdbs]
	app.to_pickle('{}/final3/pdbbind/dataset/pockets/alpha_pockets_{}.pkl'.format(DIR_AP, facts[0]))
	del app



for facts in [(10,1),(12,1)]:
	rns = [str(x)+str(y) for x,y in combinations_with_replacement([0,1,2,3,5,6,9],2)]
	factory = ChemicalFeatures.BuildFeatureFactory('/spstorage/USERS/gina/Project/AP/etc/BaseFeatures'+f_list[0])
	fn = factory.GetFeatureFamilies()
	bns = 30/facts[0]
	cns = [bns*g for g in list(range(facts[0]))]
	st = 30/(inv**facts[0]) # 16 bin
	cnrs = [st*(inv**g) for g in list(range(facts[0]))]
	def bulks_pdb(pps):
		try:
			if pps not in pds1:
				ppm = Chem.rdmolfiles.MolFromPDBFile(DIR_AP+'/DB/v2019-other-PL/'+pps+'/'+pps+'_pocket.pdb', removeHs=True)
			else:
				ppm = Chem.rdmolfiles.MolFromPDBFile(DIR_AP+'/DB/refined-set/'+pps+'/'+pps+'_pocket.pdb', removeHs=True)
			ppf = factory.GetFeaturesForMol(ppm)
			ppa = ppm.GetAtoms()
			feats = filter_feats_renew(ppf)
			featss = []
			for f in feats:
				if len(f.GetAtomIds())==1:
					pa =  ppa[f.GetAtomIds()[0]]
					pi = pa.GetPDBResidueInfo()
					if pi.GetName().split(' ')[1] not in ['N','O','CA']:
						if pa.GetSymbol()=='C':
							if len([aa for aa in pa.GetNeighbors() if aa.GetSymbol() in ['O','N']])==0:
								featss.append(f)
						else:
							featss.append(f)
				else:
					featss.append(f)
			das = [f.GetAtomIds() for f in featss if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')]
			das = unique([x for x in das if das.count(x) > 1])
			featss = filter_feats_both(featss)
			vt = pairLr_dist(featss, rns, cnrs, inv, das, fn)
		except:
			vt = [0]*len(rns)*len(cns)
		return tuple(vt)
	aps = []
	for i in range(11):
		aps.extend(mclapply(pds[i*1610:(i+1)*1610], bulks_pdb, ncpu=75))
	#aps = mclapply(pds, bulks_pdb, ncpu=75)
	app = pd.DataFrame(aps)
	app.index= pds
	app.to_pickle('{}/final2/dataset/pdb/binding/pockets/pdb_pockets_{}_{}.pkl'.format(DIR_AP, facts[0], facts[1]))
	del app



subdirs=list(set(subdirs)-set(['thb']))

dstbin = 10
f_list=['','','','_all'] +['_h']*4
st = 30/(inv**dstbin) # 16 bin
cnrs = [st*(inv**g) for g in range(dstbin)]
fnss = [[0,1,2,3,5,6],[0,1,2,3,5,6],[0,1,2,3,5,6,9],[0,1,2,3,4,5,7,8,10],[0,1,2,3,4,6,7,9],[0,1,2,3,4,6,7,9,10,11,12,13,14],[0,1,2,3,4,6,7,9,10,11,12,13,14],[0,1,2,3,4,6,7,9,10,11,12,13,14,15,16,17,18]]


for dd in [4,5]:
	for ii in [6,1,2,3,4]:
		facts = [ii,dd]
		factory = ChemicalFeatures.BuildFeatureFactory('{}/etc/BaseFeatures{}.fdef'.format(DIR_AP,f_list[facts[0]]))
		fn = factory.GetFeatureFamilies()
		fns = fnss[facts[0]]
		rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
		def bulks_pdb(ssets):
			try:
				pps, facts, rns, cnrs, fn = ssets
				ppm = Chem.rdmolfiles.MolFromPDBFile('{}/DB/DUDE/all/{}/receptor.pdb'.format(DIR_AP,pps), removeHs=True)
				if ppm is None:
					ppm = Chem.rdmolfiles.MolFromPDBFile('{}/DB/DUDE/all/{}/receptor1.pdb'.format(DIR_AP,pps), removeHs=True)
					if ppm is None:
						ppm = Chem.rdmolfiles.MolFromPDBFile('{}/DB/DUDE/all/{}/receptor2.pdb'.format(DIR_AP,pps), removeHs=True)
				llm = Chem.rdmolfiles.MolFromMol2File('{}/DB/DUDE/all/{}/crystal_ligand.mol2'.format(DIR_AP,pps), removeHs=True)
				llf = factory.GetFeaturesForMol(llm)
				lluse_pos = np.array([ft.GetPos() for ft in llf])
				ppf = factory.GetFeaturesForMol(ppm)
				ppuse_pos = np.array([ft.GetPos() for ft in ppf])
				feat = [f for i,f in enumerate(ppf) if min([dst([ppuse_pos[i],lgd_po]) for lgd_po in lluse_pos])<facts[1]]
				das = []
				if facts[0] > 0:
					feat = filter_feats_renew(feat)
					if facts[0] > 1:
						das = [f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')]
						das = unique([x for x in das if das.count(x) > 1])
						if facts[0] > 4:	
							mps = ppm.GetConformer().GetPositions()							
							atms = [[am for i,am in enumerate(ppm.GetAtoms()) if min([dst([mps[i],lgd_po]) for lgd_po in lluse_pos])<facts[1]]]
							mps = [pm for pm in ppm.GetConformer().GetPositions() if min([dst([pm,lgd_po]) for lgd_po in lluse_pos])<facts[1]]
							mps = [np.stack(mps)]
							vt = feature_AA(rns, [feat], atms, mps, 0, cnrs, inv, [das], fn)
						else:
							vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
					else:
						vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
				else:
					vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
			except:
				vt = [0]*len(rns)*len(cnrs)
			return tuple(vt)
		aps = mclapply([[pps, facts, rns, cnrs, fn] for pps in subdirs], bulks_pdb, ncpu=101)
		app = pd.DataFrame(aps)
		app.index= subdirs
		app.to_pickle('{}/final3/dude/dataset/pockets/pockets_{}_{}.pkl'.format(DIR_AP, facts[0], facts[1]))
		del app


for dd in [4]:
	for ii in [11]:
		facts = [ii,dd]
		factory = ChemicalFeatures.BuildFeatureFactory('{}/etc/BaseFeatures{}.fdef'.format(DIR_AP,f_list[facts[0]]))
		fn = factory.GetFeatureFamilies()
		fns = fnss[facts[0]]
		rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
		def bulks_pdb(ssets):
			pps, facts, rns, cnrs, fn = ssets
			ppm = Chem.rdmolfiles.MolFromPDBFile('{}/DB/DUDE/all/{}/receptor.pdb'.format(DIR_AP,pps), removeHs=True)
			if ppm is None:
				ppm = Chem.rdmolfiles.MolFromPDBFile('{}/DB/DUDE/all/{}/receptor1.pdb'.format(DIR_AP,pps), removeHs=True)
				if ppm is None:
					ppm = Chem.rdmolfiles.MolFromPDBFile('{}/DB/DUDE/all/{}/receptor2.pdb'.format(DIR_AP,pps), removeHs=True)
			fls = glob.glob('{}/final3/dude/vina_docked/{}/*'.format(DIR_AP,pps))
			vtss = pd.DataFrame()
			for adidx in ['actives_final','decoys_final_0','decoys_final_10']:
				if '{}/final3/dude/vina_docked/{}/{}_docked_vina.sdf'.format(DIR_AP, pps, adidx) not in fls:
					continue
				llms = Chem.SDMolSupplier('{}/final3/dude/vina_docked/{}/{}_docked_vina.sdf'.format(DIR_AP, pps, adidx), removeHs=True)
				lluse_pos = [llm.GetConformer().GetPositions() for llm in llms if llm is not None]
				lmn = [llm.GetProp('_Name') for llm in llms if llm is not None]
				ppf = factory.GetFeaturesForMol(ppm)
				ppuse_pos = np.array([ft.GetPos() for ft in ppf])
				feats = [[f for i,f in enumerate(ppf) if min([dst([ppuse_pos[i],lgd_po]) for lgd_po in lluse_po])<facts[1]] for lluse_po in lluse_pos]
				das = []
				if facts[0] > 0:
					if facts[0] ==6:
						feats = [filter_feats_renew(feat) for feat in feats]
					das = [[f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')] for feat in feats]
					das = [unique([x for x in da if da.count(x) > 1]) for da in das]
					mps = ppm.GetConformer().GetPositions()							
					atms = [[am for i,am in enumerate(ppm.GetAtoms()) if min([dst([mps[i],lgd_po]) for lgd_po in lluse_po])<facts[1]] for lluse_po in lluse_pos]
					mps = [[pm for pm in ppm.GetConformer().GetPositions() if min([dst([pm,lgd_po]) for lgd_po in lluse_po])<facts[1]] for lluse_po in lluse_pos]
					mps = [np.stack(mp) for mp in mps]
					if facts[0]<11:
						vt = [feature_AA(rns, feats, atms, mps, idx, cnrs, inv, das, fn) for idx, feat in enumerate(feats)]
					else:
						vt = [feature_AH(rns, feats, atms, mps, idx, cnrs, inv, das, fn) for idx, feat in enumerate(feats)]
				else:
					vt = [pairLr_dist(feat, rns, cnrs, inv, das, fn) for feat in feats]
				vts = pd.DataFrame(vt)
				vts.index = lmn
				vts['target'] = pps
				if 'actives' in adidx:
					vts['label'] =1
				else:
					vts['label'] = 0
				vtss = pd.concat([vtss, vts])
			return vtss
		aps = mclapply([[pps, facts, rns, cnrs, fn] for pps in subdirs], bulks_pdb, ncpu=101)
		app = pd.concat(aps)
		app.to_pickle('{}/final3/dude/dataset/pockets/vina_pockets_{}.pkl'.format(DIR_AP, facts[0]))
		del app




for dd in [4,5]:
	for ii in [9]:
		facts = [ii,dd]
		factory = ChemicalFeatures.BuildFeatureFactory('{}/etc/BaseFeatures{}.fdef'.format(DIR_AP,f_list[facts[0]-3]))
		fn = factory.GetFeatureFamilies()
		fns = fnss[facts[0]-3]
		rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
		def bulks_pdb(ssets):
			try:
				pps, facts, rns, cnrs, fn = ssets
				ppm = Chem.rdmolfiles.MolFromPDBFile('{}/DB/DUDE/all/{}/receptor.pdb'.format(DIR_AP,pps), removeHs=True)
				if ppm is None:
					ppm = Chem.rdmolfiles.MolFromPDBFile('{}/DB/DUDE/all/{}/receptor1.pdb'.format(DIR_AP,pps), removeHs=True)
					if ppm is None:
						ppm = Chem.rdmolfiles.MolFromPDBFile('{}/DB/DUDE/all/{}/receptor2.pdb'.format(DIR_AP,pps), removeHs=True)
				llm = Chem.rdmolfiles.MolFromMol2File('{}/DB/DUDE/all/{}/crystal_ligand.mol2'.format(DIR_AP,pps), removeHs=True)
				llf = factory.GetFeaturesForMol(llm)
				lluse_pos = np.array([ft.GetPos() for ft in llf])
				ppf = factory.GetFeaturesForMol(ppm)
				ppuse_pos = np.array([ft.GetPos() for ft in ppf])
				feat = [f for i,f in enumerate(ppf) if min([dst([ppuse_pos[i],lgd_po]) for lgd_po in lluse_pos])<facts[1]]
				das = []
				das = [f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')]
				das = unique([x for x in das if das.count(x) > 1])
				mps = ppm.GetConformer().GetPositions()							
				atms = [[am for i,am in enumerate(ppm.GetAtoms()) if min([dst([mps[i],lgd_po]) for lgd_po in lluse_pos])<facts[1]]]
				mps = [pm for pm in ppm.GetConformer().GetPositions() if min([dst([pm,lgd_po]) for lgd_po in lluse_pos])<facts[1]]
				mps = [np.stack(mps)]
				vt = feature_AA(rns, [feat], atms, mps, 0, cnrs, inv, [das], fn)
			except:
				vt = [0]*len(rns)*len(cnrs)
			return tuple(vt)
		aps = mclapply([[pps, facts, rns, cnrs, fn] for pps in subdirs], bulks_pdb, ncpu=51)
		app = pd.DataFrame(aps)
		app.index= subdirs
		app.to_pickle('{}/final3/dude/dataset/pockets/pockets_{}_{}.pkl'.format(DIR_AP, facts[0], facts[1]))
		del app


for ii in [1,2,3,4]:
	facts = [ii]
	factory = ChemicalFeatures.BuildFeatureFactory('{}/etc/BaseFeatures{}.fdef'.format(DIR_AP,f_list[facts[0]]))
	fn = factory.GetFeatureFamilies()
	fns = fnss[facts[0]]
	rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
	def bulks_pdb(ssets):
		try:
			pps, facts, rns, cnrs, fn = ssets
			ppm = Chem.rdmolfiles.MolFromPDBFile('{}/DB/DUDE/all/{}/receptor.pdb'.format(DIR_AP,pps), removeHs=True)
			if ppm is None:
				ppm = Chem.rdmolfiles.MolFromPDBFile('{}/DB/DUDE/all/{}/receptor1.pdb'.format(DIR_AP,pps), removeHs=True)
				if ppm is None:
					ppm = Chem.rdmolfiles.MolFromPDBFile('{}/DB/DUDE/all/{}/receptor2.pdb'.format(DIR_AP,pps), removeHs=True)
			ppf = factory.GetFeaturesForMol(ppm)
			ppa = ppm.GetAtoms()
			feats = [f for i,f in enumerate(ppf)]
			feat = []
			for f in feats:
				if len(f.GetAtomIds())==1:
					pa =  ppa[f.GetAtomIds()[0]]
					pi = pa.GetPDBResidueInfo()
					if pi.GetName().split(' ')[1] not in ['N','O','CA']:
						if pa.GetSymbol()=='C':
							if len([aa for aa in pa.GetNeighbors() if aa.GetSymbol() in ['O','N']])==0:
								feat.append(f)
						else:
							feat.append(f)
				else:
					feat.append(f)
			das = []
			if facts[0] > 0:
				feat = filter_feats_renew(feat)
				if facts[0] > 1:
					das = [f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')]
					das = unique([x for x in das if das.count(x) > 1])
					if facts[0] > 4:
						mps = [ppm.GetConformer().GetPositions()]						
						atms = [ppm.GetAtoms()]
						mps = [np.stack(mps)]
						vt = feature_AA(rns, [feat], atms, mps, 0, cnrs, inv, [das], fn)
					else:
						vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
				else:
					vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
			else:
				vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
		except:
			vt = [0]*len(rns)*len(cnrs)
		return tuple(vt)
	aps = mclapply([[pps, facts, rns, cnrs, fn] for pps in subdirs], bulks_pdb, ncpu=40)
	app = pd.DataFrame(aps)
	app.index= subdirs
	app.to_pickle('{}/final3/dude/dataset/pockets/pdb_pockets_{}.pkl'.format(DIR_AP, facts[0]))
	del app


dtip = '/spstorage/USERS/gina/Project/AP/DB/InterpretableDTIP'

ori = pd.read_csv(dtip+'/data/pdb_ptn_loc_gpdb.csv', index_col = 0)
ori = pd.read_csv(dtip+'/data/pdb_ptn_loc_betalactam.csv', index_col = 0)

ori = pd.read_csv(dtip+'/data/covid_loc10.csv', index_col = 0)
ori = ori[ori.iloc[:,0]=='sub']
allsets = pd.read_csv(dtip+'/data/data_dataset_3D.csv', index_col = 0)

ori = pd.read_csv(dtip+'/DUDE/dud_pocket_loc15.csv', index_col = 0)
# for chem
ori = pd.read_csv(dtip+'/data/pocket_loc20.csv', index_col = 0)

ori = pd.read_csv(dtip+'/data/ori_pocket_loc.csv', index_col = 0)
pp = pd.read_csv(dtip+'/data/pocekt_loc2.csv', index_col = 0)
pp = pd.read_csv(dtip+'/data/pocket_loc20_real.csv', index_col = 0)

pp = pd.read_csv(dtip+'/data/pdb_ptn_loc_refined_20.csv', index_col =0)

pp = pd.read_csv(dtip+'/data/pdb_ptn_loc_refined.csv', index_col =0)
ori = pp[pp.iloc[:,0]=='sub']


ori = pp[pp.iloc[:,0]=='original']
sub = pp[pp.iloc[:,0]=='sub']
pp = pd.read_csv(dtip+'/data/pocket_loc10.csv', index_col = 0)

ori = pd.read_csv(dtip+'/data/sub_pocket_loc2.csv', index_col = 0)


oori = pd.read_csv(dtip+'/human/human_pocket_loc.csv', index_col = 0)
oori = pd.concat([oori, pd.read_csv(dtip+'/human/human_pocket_loc2.csv', index_col = 0)])

				pps, facts, rns, cnrs, fn = ssets
				if pps in pds1:didx = 'refined-set'
				else:didx = 'v2020-other-PL'
				ppm = Chem.rdmolfiles.MolFromPDBFile('{}/DB/{}/{}/{}_pocket.pdb'.format(DIR_AP,didx,pps,pps), removeHs=True)
				llm = Chem.SDMolSupplier('{}/DB/{}/{}/{}_ligand.sdf'.format(DIR_AP,didx,pps,pps), removeHs=True)[0]


for ii in [0,6,11]:
	facts = [ii]
	factory = ChemicalFeatures.BuildFeatureFactory('{}/etc/BaseFeatures{}.fdef'.format(DIR_AP,f_list[facts[0]]))
	fn = factory.GetFeatureFamilies()
	fns = fnss[facts[0]]
	rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
	def bulks_pdb(ssets):
		try:
			pps, facts, rns, cnrs, fn = ssets
			ppm = Chem.MolFromPDBFile('{}/PDB/{}_{}.pdb'.format(dtip,pps[0:4].upper(),pps[4].upper()), removeHs=True)
			ppf = factory.GetFeaturesForMol(ppm)
			ppuse_pos = np.array([ft.GetPos() for ft in ppf])
			subps = ori[ori.iloc[:,1]==pps]
			vts = []
			for ilen in range(len(subps)):
				x_min, x_max,y_min, y_max, z_min, z_max= subps.iloc[ilen, 2:8]
				feat = []
				for i, f in enumerate(ppf):
					if x_min < ppuse_pos[i][0] < x_max and y_min < ppuse_pos[i][1] < y_max and z_min < ppuse_pos[i][2] < z_max:
						feat.append(f)
				das = [[]]
				if len(feat)>5:
					if facts[0] > 0:
						feat = filter_feats_renew(feat)
						if facts[0] > 1:
							das = [f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')]
							das = unique([x for x in das if das.count(x) > 1])
							if facts[0] > 4:
								mps = ppm.GetConformer().GetPositions()							
								atms = [[am for i,am in enumerate(ppm.GetAtoms()) if x_min < mps[i][0] < x_max and y_min < mps[i][1] < y_max and z_min < mps[i][2] < z_max]]
								mps = [pm for pm in ppm.GetConformer().GetPositions() if x_min < pm[0] < x_max and y_min < pm[1] < y_max and z_min < pm[2] < z_max]
								mps = [np.stack(mps)]
								if facts[0] ==6:
									vt = feature_AA(rns, [feat], atms, mps, 0, cnrs, inv, [das], fn)
								else:
									vt = feature_AH(rns, [feat], atms, mps, 0, cnrs, inv, [das], fn)
							else:
								vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
						else:
							vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
					else:
						vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
					vts.append(vt)
			dfs = pd.DataFrame(vts).drop_duplicates()
			dfs['pdb'] = pps
		except:
			dfs = pd.DataFrame([[0]*len(rns)*len(cnrs)]).drop_duplicates()
			dfs['pdb'] = pps
		return dfs
	aps = mclapply([[pps, facts, rns, cnrs, fn] for pps in set(ori.iloc[:,1])], bulks_pdb, ncpu=104)
	app = pd.concat(aps)
	app.index= list(app.pdb)
	app = app.iloc[:,0:-1]
	app.to_pickle('{}/ori_pockets_{}.pkl'.format(dtip, facts[0]))
	del app



for ii in [6,11,0]:
	facts = [ii]
	factory = ChemicalFeatures.BuildFeatureFactory('{}/etc/BaseFeatures{}.fdef'.format(DIR_AP,f_list[facts[0]]))
	fn = factory.GetFeatureFamilies()
	fns = fnss[facts[0]]
	rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
	def bulks_pdb(ssets):
		try:
			pps, facts, rns, cnrs, fn = ssets
			ppm = Chem.MolFromPDBFile('{}/human/{}.pdb'.format(dtip,pps), removeHs=True)
			ppf = factory.GetFeaturesForMol(ppm)
			ppuse_pos = np.array([ft.GetPos() for ft in ppf])
			subps = ori[ori.iloc[:,1]==pps]
			vts = []
			for ilen in range(len(subps)):
				x_min, x_max,y_min, y_max, z_min, z_max= subps.iloc[ilen, 2:8]
				feat = []
				for i, f in enumerate(ppf):
					if x_min < ppuse_pos[i][0] < x_max and y_min < ppuse_pos[i][1] < y_max and z_min < ppuse_pos[i][2] < z_max:
						feat.append(f)
				das = [[]]
				if len(feat)>5:
					if facts[0] > 0:
						feat = filter_feats_renew(feat)
						if facts[0] > 1:
							das = [f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')]
							das = unique([x for x in das if das.count(x) > 1])
							if facts[0] > 4:
								mps = ppm.GetConformer().GetPositions()							
								atms = [[am for i,am in enumerate(ppm.GetAtoms()) if x_min < mps[i][0] < x_max and y_min < mps[i][1] < y_max and z_min < mps[i][2] < z_max]]
								mps = [pm for pm in ppm.GetConformer().GetPositions() if x_min < pm[0] < x_max and y_min < pm[1] < y_max and z_min < pm[2] < z_max]
								mps = [np.stack(mps)]
								if facts[0] ==6:
									vt = feature_AA(rns, [feat], atms, mps, 0, cnrs, inv, [das], fn)
								else:
									vt = feature_AH(rns, [feat], atms, mps, 0, cnrs, inv, [das], fn)
							else:
								vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
						else:
							vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
					else:
						vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
					vts.append(vt)
			dfs = pd.DataFrame(vts).drop_duplicates()
			dfs['pdb'] = pps
		except:
			dfs = pd.DataFrame([[0]*len(rns)*len(cnrs)]).drop_duplicates()
			dfs['pdb'] = pps
		return dfs
	aps = mclapply([[pps, facts, rns, cnrs, fn] for pps in ppss], bulks_pdb, ncpu=100)
	app = pd.concat(aps)
	app.index= list(app.pdb)
	app = app.iloc[:,0:-1]
	app.to_pickle('{}/human/sub_pockets_{}_20_1.pkl'.format(dtip, facts[0]))
	del app



for ii in [11]:
	facts = [ii]
	factory = ChemicalFeatures.BuildFeatureFactory('{}/etc/BaseFeatures{}.fdef'.format(DIR_AP,f_list[facts[0]]))
	fn = factory.GetFeatureFamilies()
	fns = fnss[facts[0]]
	rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
	def bulks_pdb(ssets):
		try:
			pps, facts, rns, cnrs, fn = ssets
			ppm = Chem.rdmolfiles.MolFromPDBFile('{}/DB/v2020-other-PL/{}/{}_protein.pdb'.format(DIR_AP,pps,pps), removeHs=True)
			ppf = factory.GetFeaturesForMol(ppm)
			ppuse_pos = np.array([ft.GetPos() for ft in ppf])
			subps = ori[ori.iloc[:,1]==pps]
			vts = []
			for ilen in range(len(subps)):
				lgd, x_min, x_max,y_min, y_max, z_min, z_max= subps.iloc[ilen, 2:9]
				feat = []
				for i, f in enumerate(ppf):
					if x_min < ppuse_pos[i][0] < x_max and y_min < ppuse_pos[i][1] < y_max and z_min < ppuse_pos[i][2] < z_max:
						feat.append(f)
				das = [[]]
				if len(feat)>5:
					if facts[0] > 0:
						feat = filter_feats_renew(feat)
						if facts[0] > 1:
							das = [f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')]
							das = unique([x for x in das if das.count(x) > 1])
							if facts[0] > 4:
								mps = ppm.GetConformer().GetPositions()							
								atms = [[am for i,am in enumerate(ppm.GetAtoms()) if x_min < mps[i][0] < x_max and y_min < mps[i][1] < y_max and z_min < mps[i][2] < z_max]]
								mps = [pm for pm in ppm.GetConformer().GetPositions() if x_min < pm[0] < x_max and y_min < pm[1] < y_max and z_min < pm[2] < z_max]
								mps = [np.stack(mps)]
								if facts[0] ==6:
									vt = feature_AA(rns, [feat], atms, mps, 0, cnrs, inv, [das], fn)
								else:
									vt = feature_AH(rns, [feat], atms, mps, 0, cnrs, inv, [das], fn)
							else:
								vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
						else:
							vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
					else:
						vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
					vts.append(vt+[lgd])
			dfs = pd.DataFrame(vts).drop_duplicates()
			dfs['pdb'] = pps
		except:
			dfs = pd.DataFrame([[0]*(len(rns)*len(cnrs)+1)]).drop_duplicates()
			dfs['pdb'] = pps
		return dfs
	aps = mclapply([[pps, facts, rns, cnrs, fn] for pps in set(ori.iloc[:,1])], bulks_pdb, ncpu=100)
	app = pd.concat(aps)
	app.index= list(app.pdb)
	app = app.iloc[:,0:-1]
	app.to_pickle('{}/pdbval/sub_pockets_{}_20_3.pkl'.format(dtip, facts[0]))
	del app



for ii in [11]:
	facts = [ii]
	factory = ChemicalFeatures.BuildFeatureFactory('{}/etc/BaseFeatures{}.fdef'.format(DIR_AP,f_list[facts[0]]))
	fn = factory.GetFeatureFamilies()
	fns = fnss[facts[0]]
	rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
	def bulks_pdb(ssets):
		try:
			pps, facts, rns, cnrs, fn = ssets
			ppm = Chem.rdmolfiles.MolFromPDBFile('{}/DB/betalactam/pdb/{}_protein.pdb'.format(DIR_AP,pps,pps), removeHs=True)
			ppf = factory.GetFeaturesForMol(ppm)
			ppuse_pos = np.array([ft.GetPos() for ft in ppf])
			subps = ori[ori.iloc[:,1]==pps]
			vts = []
			for ilen in range(len(subps)):
				lgd_name, lgd, x_min, x_max,y_min, y_max, z_min, z_max= subps.iloc[ilen, 2:10]
				feat = []
				for i, f in enumerate(ppf):
					if x_min < ppuse_pos[i][0] < x_max and y_min < ppuse_pos[i][1] < y_max and z_min < ppuse_pos[i][2] < z_max:
						feat.append(f)
				das = [[]]
				if len(feat)>5:
					if facts[0] > 0:
						feat = filter_feats_renew2(feat)
						if facts[0] > 1:
							das = [f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')]
							das = unique([x for x in das if das.count(x) > 1])
							feat = filter_feats_both(feat)
							if facts[0] > 4:
								mps = ppm.GetConformer().GetPositions()							
								atms = [[am for i,am in enumerate(ppm.GetAtoms()) if x_min < mps[i][0] < x_max and y_min < mps[i][1] < y_max and z_min < mps[i][2] < z_max]]
								mps = [pm for pm in ppm.GetConformer().GetPositions() if x_min < pm[0] < x_max and y_min < pm[1] < y_max and z_min < pm[2] < z_max]
								mps = [np.stack(mps)]
								if facts[0] ==6:
									vt = feature_AA(rns, [feat], atms, mps, 0, cnrs, inv, [das], fn)
								else:
									vt = feature_AH(rns, [feat], atms, mps, 0, cnrs, inv, [das], fn)
							else:
								vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
						else:
							vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
					else:
						vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
					vts.append(vt+[lgd_name, lgd])
			dfs = pd.DataFrame(vts).drop_duplicates()
			dfs['pdb'] = pps
		except:
			dfs = pd.DataFrame([[0]*(len(rns)*len(cnrs)+2)]).drop_duplicates()
			dfs['pdb'] = pps
		return dfs
	aps = mclapply([[pps, facts, rns, cnrs, fn] for pps in set(ori.iloc[:,1])], bulks_pdb, ncpu=180)
	app = pd.concat(aps)
	app.index= list(app.pdb)
	app = app.iloc[:,0:-1]
	app.to_pickle('{}/pdbval/sub_pockets_{}_20_renew2.pkl'.format(dtip, facts[0]))
	del app

				llm = Chem.SDMolSupplier('{}/DB/{}/{}/{}_ligand.sdf'.format(DIR_AP,didx,pps,pps), removeHs=True)[0]
				if llm is None:
					llm = Chem.rdmolfiles.MolFromMol2File('{}/DB/{}/{}/{}_ligand.mol2'.format(DIR_AP,didx,pps,pps), removeHs=True)
				llf = factory.GetFeaturesForMol(llm)
				lluse_pos = np.array([ft.GetPos() for ft in llf])
				ppf = factory.GetFeaturesForMol(ppm)
				ppuse_pos = np.array([ft.GetPos() for ft in ppf])
				feat = [f for i,f in enumerate(ppf) if min([dst([ppuse_pos[i],lgd_po]) for lgd_po in lluse_pos])<facts[1]]


## llm 
for ii in [11]:
	facts = [ii,4]
	factory = ChemicalFeatures.BuildFeatureFactory('{}/etc/BaseFeatures{}.fdef'.format(DIR_AP,f_list[facts[0]]))
	fn = factory.GetFeatureFamilies()
	fns = fnss[facts[0]]
	rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
	def bulks_pdb(ssets):
		try:
			pps, facts, rns, cnrs, fn = ssets
			ppm = Chem.rdmolfiles.MolFromPDBFile('{}/DB/refined-set/{}/{}_protein.pdb'.format(DIR_AP,pps,pps), removeHs=True)
			llm = Chem.SDMolSupplier('{}/DB/refined-set/{}/{}_ligand.sdf'.format(DIR_AP,pps,pps), removeHs=True)[0]
			if llm is None:
				llm = Chem.rdmolfiles.MolFromMol2File('{}/DB/refined-set/{}/{}_ligand.mol2'.format(DIR_AP,pps,pps), removeHs=True)
			llf = factory.GetFeaturesForMol(llm)
			lluse_pos = np.array([ft.GetPos() for ft in llf])
			ppf = factory.GetFeaturesForMol(ppm)
			ppuse_pos = np.array([ft.GetPos() for ft in ppf])
			feat = [f for i,f in enumerate(ppf) if min([dst([ppuse_pos[i],lgd_po]) for lgd_po in lluse_pos])<facts[1]]
			das = [[]]
			if len(feat)>5:
				if facts[0] > 0:
					feat = filter_feats_renew(feat)
					if facts[0] > 1:
						das = [f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')]
						das = unique([x for x in das if das.count(x) > 1])
						if facts[0] > 4:
							mps = ppm.GetConformer().GetPositions()							
							atms = [[am for i,am in enumerate(ppm.GetAtoms()) if min([dst([mps[i],lgd_po]) for lgd_po in lluse_pos])<facts[1]]]
							mps = [pm for pm in ppm.GetConformer().GetPositions() if min([dst([pm,lgd_po]) for lgd_po in lluse_pos])<facts[1]]
							mps = [np.stack(mps)]
							if facts[0] ==6:
								vt = feature_AA(rns, [feat], atms, mps, 0, cnrs, inv, [das], fn)
							else:
								vt = feature_AH(rns, [feat], atms, mps, 0, cnrs, inv, [das], fn)
						else:
							vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
					else:
						vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
				else:
					vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
			else:
				vt = [0]*len(rns)*len(cnrs)
		except:
			vt = [0]*len(rns)*len(cnrs)
		return tuple(vt)
	aps = mclapply([[pps, facts, rns, cnrs, fn] for pps in pds], bulks_pdb, ncpu=110)
	app = pd.DataFrame(aps)
	app.index= pds
	app.to_pickle('{}/pdbval/pdb_sub_pockets_{}_20_4.pkl'.format(dtip, facts[0]))
	del app

## llm 
for ii in [11]:
	facts = [ii,5]
	factory = ChemicalFeatures.BuildFeatureFactory('{}/etc/BaseFeatures{}.fdef'.format(DIR_AP,f_list[facts[0]]))
	fn = factory.GetFeatureFamilies()
	fns = fnss[facts[0]]
	rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
	def bulks_pdb(ssets):
		try:
			pps, facts, rns, cnrs, fn = ssets
			ppm = Chem.rdmolfiles.MolFromPDBFile('{}/DB/refined-set/{}/{}_protein.pdb'.format(DIR_AP,pps,pps), removeHs=True)
			llm = Chem.SDMolSupplier('{}/DB/refined-set/{}/{}_ligand.sdf'.format(DIR_AP,pps,pps), removeHs=True)[0]
			if llm is None:
				llm = Chem.rdmolfiles.MolFromMol2File('{}/DB/refined-set/{}/{}_ligand.mol2'.format(DIR_AP,pps,pps), removeHs=True)
			llf = factory.GetFeaturesForMol(llm)
			lluse_pos = llm.GetConformer().GetPositions()
			ppf = factory.GetFeaturesForMol(ppm)
			ppuse_pos = np.array([ft.GetPos() for ft in ppf])
			feat = [f for i,f in enumerate(ppf) if min([dst([ppuse_pos[i],lgd_po]) for lgd_po in lluse_pos])<facts[1]]
			das = [[]]
			if len(feat)>5:
				if facts[0] > 0:
					feat = filter_feats_renew(feat)
					if facts[0] > 1:
						das = [f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')]
						das = unique([x for x in das if das.count(x) > 1])
						if facts[0] > 4:
							mps = ppm.GetConformer().GetPositions()							
							atms = [[am for i,am in enumerate(ppm.GetAtoms()) if min([dst([mps[i],lgd_po]) for lgd_po in lluse_pos])<facts[1]]]
							mps = [pm for pm in ppm.GetConformer().GetPositions() if min([dst([pm,lgd_po]) for lgd_po in lluse_pos])<facts[1]]
							mps = [np.stack(mps)]
							if facts[0] ==6:
								vt = feature_AA(rns, [feat], atms, mps, 0, cnrs, inv, [das], fn)
							else:
								vt = feature_AH(rns, [feat], atms, mps, 0, cnrs, inv, [das], fn)
						else:
							vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
					else:
						vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
				else:
					vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
		except:
			vt = [0]*len(rns)*len(cnrs)
		return tuple(vt)
	aps = mclapply([[pps, facts, rns, cnrs, fn] for pps in pds], bulks_pdb, ncpu=80)
	app = pd.DataFrame(aps)
	app.index= pds
	app.to_pickle('{}/pdbval/pdb_atom_sub_pockets_{}_20.pkl'.format(dtip, facts[0]))
	del app

for ii in [11]:
	facts = [ii]
	factory = ChemicalFeatures.BuildFeatureFactory('{}/etc/BaseFeatures{}.fdef'.format(DIR_AP,f_list[facts[0]]))
	fn = factory.GetFeatureFamilies()
	fns = fnss[facts[0]]
	rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
	def bulks_pdb(ssets):
		try:
			pps, facts, rns, cnrs, fn = ssets
			ppm = Chem.rdmolfiles.MolFromPDBFile('{}/DB/v2020-other-PL/{}/{}_protein.pdb'.format(DIR_AP,pps,pps), removeHs=True)
			ppf = factory.GetFeaturesForMol(ppm)
			ppuse_pos = np.array([ft.GetPos() for ft in ppf])
			subps = pori[pori.iloc[:,1]==pps]
			vts = []
			for ilen in range(len(subps)):
				lgd, x_min, x_max,y_min, y_max, z_min, z_max= subps.iloc[ilen, 2:9]
				feat = []
				for i, f in enumerate(ppf):
					if x_min < ppuse_pos[i][0] < x_max and y_min < ppuse_pos[i][1] < y_max and z_min < ppuse_pos[i][2] < z_max:
						feat.append(f)
				das = [[]]
				if len(feat)>5:
					if facts[0] > 0:
						feat = filter_feats_renew2(feat)
						if facts[0] > 1:
							das = [f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')]
							das = unique([x for x in das if das.count(x) > 1])
							# feat = filter_feats_both(feat)
							if facts[0] > 4:
								mps = ppm.GetConformer().GetPositions()							
								atms = [[am for i,am in enumerate(ppm.GetAtoms()) if x_min < mps[i][0] < x_max and y_min < mps[i][1] < y_max and z_min < mps[i][2] < z_max]]
								mps = [pm for pm in ppm.GetConformer().GetPositions() if x_min < pm[0] < x_max and y_min < pm[1] < y_max and z_min < pm[2] < z_max]
								mps = [np.stack(mps)]
								if facts[0] ==6:
									vt = feature_AA(rns, [feat], atms, mps, 0, cnrs, inv, [das], fn)
								else:
									vt = feature_AH(rns, [feat], atms, mps, 0, cnrs, inv, [das], fn)
							else:
								vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
						else:
							vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
					else:
						vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
					vts.append(vt+[lgd])
			dfs = pd.DataFrame(vts).drop_duplicates()
			dfs['pdb'] = pps
		except:
			dfs = pd.DataFrame([[0]*(len(rns)*len(cnrs)+1)]).drop_duplicates()
			dfs['pdb'] = pps
		return dfs
	aps = mclapply([[pps, facts, rns, cnrs, fn] for pps in set(pori.iloc[:,1])], bulks_pdb, ncpu=2)
	app = pd.concat(aps)
	app.index= list(app.pdb)
	app.index =[cd.upper() for cd in list(app.index)]
	app = app.iloc[:,0:-1]
	app.to_pickle('{}/sub_pockets_{}_20_renew2_no_both2.pkl'.format(chemm, facts[0]))
	del app

app3.to_pickle('{}/sub_pockets_{}_20_renew2_no_both.pkl'.format(chemm, facts[0]))


for ii in [11]:
	facts = [ii]
	factory = ChemicalFeatures.BuildFeatureFactory('{}/etc/BaseFeatures{}.fdef'.format(DIR_AP,f_list[facts[0]]))
	fn = factory.GetFeatureFamilies()
	fns = fnss[facts[0]]
	rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
	def bulks_pdb(ssets):
		try:
			pps, facts, rns, cnrs, fn = ssets
			ppm = Chem.rdmolfiles.MolFromPDBFile('{}/DB/refined-set/{}/{}_protein.pdb'.format(DIR_AP,pps,pps), removeHs=True)
			ppf = factory.GetFeaturesForMol(ppm)
			ppuse_pos = np.array([ft.GetPos() for ft in ppf])
			subps = ori[ori.iloc[:,1]==pps]
			vts = []
			for ilen in range(len(subps)):
				lgd, x_min, x_max,y_min, y_max, z_min, z_max= subps.iloc[ilen, 2:9]
				feat = []
				for i, f in enumerate(ppf):
					if x_min < ppuse_pos[i][0] < x_max and y_min < ppuse_pos[i][1] < y_max and z_min < ppuse_pos[i][2] < z_max:
						feat.append(f)
				das = [[]]
				if len(feat)>5:
					if facts[0] > 0:
						feat = filter_feats_renew2(feat)
						if facts[0] > 1:
							das = [f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')]
							das = unique([x for x in das if das.count(x) > 1])
							# feat = filter_feats_both(feat)
							if facts[0] > 4:
								mps = ppm.GetConformer().GetPositions()							
								atms = [[am for i,am in enumerate(ppm.GetAtoms()) if x_min < mps[i][0] < x_max and y_min < mps[i][1] < y_max and z_min < mps[i][2] < z_max]]
								mps = [pm for pm in ppm.GetConformer().GetPositions() if x_min < pm[0] < x_max and y_min < pm[1] < y_max and z_min < pm[2] < z_max]
								mps = [np.stack(mps)]
								if facts[0] ==6:
									vt = feature_AA(rns, [feat], atms, mps, 0, cnrs, inv, [das], fn)
								else:
									vt = feature_AH(rns, [feat], atms, mps, 0, cnrs, inv, [das], fn)
							else:
								vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
						else:
							vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
					else:
						vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
					vts.append(vt + [lgd])
			dfs = pd.DataFrame(vts).drop_duplicates()
			dfs['pdb'] = pps
		except:
			dfs = pd.DataFrame([[0]*(len(rns)*len(cnrs)+1)]).drop_duplicates()
			dfs['pdb'] = pps
		return dfs
	aps = mclapply([[pps, facts, rns, cnrs, fn] for pps in set(ori.iloc[:,1])], bulks_pdb, ncpu=100)
	app = pd.concat(aps)
	app.index= list(app.pdb)
	app = app.iloc[:,0:-1]
	app.to_pickle('{}/pdbval/sub_pockets_{}_20_renew2_no_both.pkl'.format(dtip, facts[0]))
	del app

dtip='/spstorage/USERS/gina/Project/AP/DB/InterpretableDTIP'
rrridx = pd.read_pickle(dtip+'/rridx_13K.pkl')
cridx = pd.read_pickle(dtip+'/cridx_220.pkl')

for ii in [11]:
	facts = [ii]
	factory = ChemicalFeatures.BuildFeatureFactory('{}/etc/BaseFeatures{}.fdef'.format(DIR_AP,f_list[facts[0]]))
	fn = factory.GetFeatureFamilies()
	fns = fnss[facts[0]]
	rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
	def bulks_pdb(ssets):
		try:
			pps, facts, rns, cnrs, fn = ssets
			ppm = Chem.rdmolfiles.MolFromPDBFile('{}/DB/v2020-other-PL/{}/{}_protein.pdb'.format(DIR_AP,pps,pps), removeHs=True)
			ppf = factory.GetFeaturesForMol(ppm)
			ppuse_pos = np.array([ft.GetPos() for ft in ppf])
			subps = ori[ori.iloc[:,1]==pps]
			vts = []
			for ilen in range(len(subps)):
				lgd, x_min, x_max,y_min, y_max, z_min, z_max= subps.iloc[ilen, 2:9]
				feat = []
				for i, f in enumerate(ppf):
					if x_min < ppuse_pos[i][0] < x_max and y_min < ppuse_pos[i][1] < y_max and z_min < ppuse_pos[i][2] < z_max:
						feat.append(f)
				das = [[]]
				if len(feat)>5:
					if facts[0] > 0:
						feat = filter_feats_renew2(feat)
						if facts[0] > 1:
							das = [f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')]
							das = unique([x for x in das if das.count(x) > 1])
							# feat = filter_feats_both(feat)
							if facts[0] > 4:
								mps = ppm.GetConformer().GetPositions()							
								atms = [[am for i,am in enumerate(ppm.GetAtoms()) if x_min < mps[i][0] < x_max and y_min < mps[i][1] < y_max and z_min < mps[i][2] < z_max]]
								mps = [pm for pm in ppm.GetConformer().GetPositions() if x_min < pm[0] < x_max and y_min < pm[1] < y_max and z_min < pm[2] < z_max]
								mps = [np.stack(mps)]
								if facts[0] ==6:
									vt = feature_AA(rns, [feat], atms, mps, 0, cnrs, inv, [das], fn)
								else:
									vt = feature_AH(rns, [feat], atms, mps, 0, cnrs, inv, [das], fn)
							else:
								vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
						else:
							vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
					else:
						vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
					vts.append(vt + [lgd])
			dfs = pd.DataFrame(vts).drop_duplicates()
			dfs['pdb'] = pps
		except:
			dfs = pd.DataFrame([[0]*(len(rns)*len(cnrs)+1)]).drop_duplicates()
			dfs['pdb'] = pps
		return dfs
	aps = mclapply([[pps, facts, rns, cnrs, fn] for pps in ppdd], bulks_pdb, ncpu=65)
	app = pd.concat(aps)
	app.index= list(app.pdb)
	app = app.iloc[:,0:-1]
	app.to_pickle('{}/pdbval/sub_pockets_{}_20_renew2_no_both_v2020.pkl'.format(dtip, facts[0]))
	del app



for ii in [11]:
	facts = [ii]
	factory = ChemicalFeatures.BuildFeatureFactory('{}/etc/BaseFeatures{}.fdef'.format(DIR_AP,f_list[facts[0]]))
	fn = factory.GetFeatureFamilies()
	fns = fnss[facts[0]]
	rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
	def bulks_pdb(ssets):
		try:
			pps, facts, rns, cnrs, fn = ssets
			ppm = Chem.rdmolfiles.MolFromPDBFile('{}/DB/refined-set/{}/{}_protein.pdb'.format(DIR_AP,pps,pps), removeHs=True)
			ppf = factory.GetFeaturesForMol(ppm)
			ppuse_pos = np.array([ft.GetPos() for ft in ppf])
			subps = ori[ori.iloc[:,1]==pps]
			vts = []
			for ilen in range(len(subps)):
				lgd, x_min, x_max,y_min, y_max, z_min, z_max= subps.iloc[ilen, 2:9]
				feat = []
				for i, f in enumerate(ppf):
					if x_min < ppuse_pos[i][0] < x_max and y_min < ppuse_pos[i][1] < y_max and z_min < ppuse_pos[i][2] < z_max:
						feat.append(f)
				das = [[]]
				if len(feat)>5:
					if facts[0] > 0:
						feat = filter_feats_renew2(feat)
						if facts[0] > 1:
							das = [f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')]
							das = unique([x for x in das if das.count(x) > 1])
							if facts[0] > 4:
								mps = ppm.GetConformer().GetPositions()							
								atms = [[am for i,am in enumerate(ppm.GetAtoms()) if x_min < mps[i][0] < x_max and y_min < mps[i][1] < y_max and z_min < mps[i][2] < z_max]]
								mps = [pm for pm in ppm.GetConformer().GetPositions() if x_min < pm[0] < x_max and y_min < pm[1] < y_max and z_min < pm[2] < z_max]
								mps = [np.stack(mps)]
								if facts[0] ==6:
									vt = feature_AA(rns, [feat], atms, mps, 0, cnrs, inv, [das], fn)
								else:
									vt = feature_AH(rns, [feat], atms, mps, 0, cnrs, inv, [das], fn)
							else:
								vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
						else:
							vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
					else:
						vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
					vts.append(vt + [lgd])
			dfs = pd.DataFrame(vts).drop_duplicates()
			dfs['pdb'] = pps
		except:
			dfs = pd.DataFrame([[0]*(len(rns)*len(cnrs)+1)]).drop_duplicates()
			dfs['pdb'] = pps
		return dfs
	aps = mclapply([[pps, facts, rns, cnrs, fn] for pps in set(ori.iloc[:,1])], bulks_pdb, ncpu=120)
	app = pd.concat(aps)
	app.index= list(app.pdb)
	app = app.iloc[:,0:-1]
	app.to_pickle('{}/pdbval/sub_pockets_{}_20_renew2_lgd.pkl'.format(dtip, facts[0]))
	del app


for ii in [11]:
	facts = [ii]
	factory = ChemicalFeatures.BuildFeatureFactory('{}/etc/BaseFeatures{}.fdef'.format(DIR_AP,f_list[facts[0]]))
	fn = factory.GetFeatureFamilies()
	fns = fnss[facts[0]]
	rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
	def bulks_pdb(ssets):
		try:
			pps, facts, rns, cnrs, fn = ssets
			ppm = Chem.MolFromPDBFile('{}/PDB/{}_{}.pdb'.format(dtip,pps[0:4].upper(),pps[4].upper()), removeHs=True)
			ppf = factory.GetFeaturesForMol(ppm)
			ppa = ppm.GetAtoms()
			ppuse_pos = np.array([ft.GetPos() for ft in ppf])
			subps = pp[pp.iloc[:,1]==pps]
			vts = []
			for ilen in range(len(subps)):
				x_min, x_max,y_min, y_max, z_min, z_max= subps.iloc[ilen, 2:8]
				feats = []
				for i, f in enumerate(ppf):
					if x_min < ppuse_pos[i][0] < x_max and y_min < ppuse_pos[i][1] < y_max and z_min < ppuse_pos[i][2] < z_max:
						feats.append(f)
				feat = []
				for f in feats:
					if len(f.GetAtomIds())==1:
						pa =  ppa[f.GetAtomIds()[0]]
						pi = pa.GetPDBResidueInfo()
						if pi.GetName().split(' ')[1] not in ['N','O','CA']:
							if pa.GetSymbol()=='C':
								if len([aa for aa in pa.GetNeighbors() if aa.GetSymbol() in ['O','N']])==0:
									feat.append(f)
							else:
								feat.append(f)
					else:
						feat.append(f)
				das = [[]]
				if len(feat)>5:
					if facts[0] > 0:
						feat = filter_feats_renew2(feat)
						if facts[0] > 1:
							das = [f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')]
							das = unique([x for x in das if das.count(x) > 1])
							# feat = filter_feats_both(feat)
							if facts[0] > 4:
								mps = ppm.GetConformer().GetPositions()							
								atms = [[am for i,am in enumerate(ppm.GetAtoms()) if x_min < mps[i][0] < x_max and y_min < mps[i][1] < y_max and z_min < mps[i][2] < z_max]]
								mps = [pm for pm in ppm.GetConformer().GetPositions() if x_min < pm[0] < x_max and y_min < pm[1] < y_max and z_min < pm[2] < z_max]
								mps = [np.stack(mps)]
								if facts[0] ==6:
									vt = feature_AA(rns, [feat], atms, mps, 0, cnrs, inv, [das], fn)
								else:
									vt = feature_AH(rns, [feat], atms, mps, 0, cnrs, inv, [das], fn)
							else:
								vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
						else:
							vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
					else:
						vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
					vts.append(vt)
			dfs = pd.DataFrame(vts).drop_duplicates()
			dfs['pdb'] = pps
		except:
			dfs = pd.DataFrame([[0]*len(rns)*len(cnrs)]).drop_duplicates()
			dfs['pdb'] = pps
		return dfs
	aps = mclapply([[pps, facts, rns, cnrs, fn] for pps in set(allsets.pdb_id)], bulks_pdb, ncpu=100)
	app = pd.concat(aps)
	app.index= list(app.pdb)
	app = app.iloc[:,0:-1]
	app.to_pickle('{}/sub_pockets_{}_20_side_renew2_z_no_both_20.pkl'.format(dtip, facts[0]))
	del app


for ii in [11]:
	facts = [ii]
	factory = ChemicalFeatures.BuildFeatureFactory('{}/etc/BaseFeatures{}.fdef'.format(DIR_AP,f_list[facts[0]]))
	fn = factory.GetFeatureFamilies()
	fns = fnss[facts[0]]
	rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
	def bulks_pdb(ssets):
		try:
			pps, facts, rns, cnrs, fn = ssets
			ppm = Chem.MolFromPDBFile('{}/PDB/{}_{}.pdb'.format(dtip,pps[0:4].upper(),pps[4].upper()), removeHs=True)
			ppf = factory.GetFeaturesForMol(ppm)
			ppa = ppm.GetAtoms()
			ppuse_pos = np.array([ft.GetPos() for ft in ppf])
			subps = ori[ori.iloc[:,1]==pps]
			vts = []
			for ilen in range(len(subps)):
				x_min, x_max,y_min, y_max, z_min, z_max= subps.iloc[ilen, 2:8]
				feats = []
				for i, f in enumerate(ppf):
					if x_min < ppuse_pos[i][0] < x_max and y_min < ppuse_pos[i][1] < y_max and z_min < ppuse_pos[i][2] < z_max:
						feats.append(f)
				feat = []
				for f in feats:
					if len(f.GetAtomIds())!=1:
						feat.append(f)
					# if len(f.GetAtomIds())==1:
					# 	pa =  ppa[f.GetAtomIds()[0]]
					# 	pi = pa.GetPDBResidueInfo()
					# 	if pi.GetName().split(' ')[1] not in ['N','O','CA']:
					# 		if pa.GetSymbol()=='C':
					# 			if len([aa for aa in pa.GetNeighbors() if aa.GetSymbol() in ['O','N']])==0:
					# 				feat.append(f)
					# 		else:
					# 			feat.append(f)
					# else:
						# feat.append(f)
				das = [[]]
				if len(feat)>5:
					if facts[0] > 0:
						feat = filter_feats_renew2(feat)
						if facts[0] > 1:
							das = [f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')]
							das = unique([x for x in das if das.count(x) > 1])
							if facts[0] > 4:
								mps = ppm.GetConformer().GetPositions()							
								atms = [[am for i,am in enumerate(ppm.GetAtoms()) if x_min < mps[i][0] < x_max and y_min < mps[i][1] < y_max and z_min < mps[i][2] < z_max]]
								mps = [pm for pm in ppm.GetConformer().GetPositions() if x_min < pm[0] < x_max and y_min < pm[1] < y_max and z_min < pm[2] < z_max]
								mps = [np.stack(mps)]
								if facts[0] ==6:
									vt = feature_AA(rns, [feat], atms, mps, 0, cnrs, inv, [das], fn)
								else:
									vt = feature_AH(rns, [feat], atms, mps, 0, cnrs, inv, [das], fn)
							else:
								vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
						else:
							vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
					else:
						vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
					vts.append(vt)
			dfs = pd.DataFrame(vts).drop_duplicates()
			dfs['pdb'] = pps
		except:
			dfs = pd.DataFrame([[0]*len(rns)*len(cnrs)]).drop_duplicates()
			dfs['pdb'] = pps
		return dfs
	aps = mclapply([[pps, facts, rns, cnrs, fn] for pps in set(allsets.pdb_id)], bulks_pdb, ncpu=100)
	app = pd.concat(aps)
	app.index= list(app.pdb)
	app = app.iloc[:,0:-1]
	app.to_pickle('{}/sub_pockets_{}_20_side2_renew2.pkl'.format(dtip, facts[0]))
	del app



for ii in [11]:
	facts = [ii]
	factory = ChemicalFeatures.BuildFeatureFactory('{}/etc/BaseFeatures{}.fdef'.format(DIR_AP,f_list[facts[0]]))
	fn = factory.GetFeatureFamilies()
	fns = fnss[facts[0]]
	rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
	def bulks_pdb(ssets):
		try:
			pps, facts, rns, cnrs, fn = ssets
			ppm = Chem.MolFromPDBFile('{}/PDB/{}_{}.pdb'.format(dtip,pps[0:4].upper(),pps[4].upper()), removeHs=True)
			ppf = factory.GetFeaturesForMol(ppm)
			ppuse_pos = np.array([ft.GetPos() for ft in ppf])
			subps = ori[ori.iloc[:,1]==pps]
			vts = []
			for ilen in range(len(subps)):
				x_min, x_max,y_min, y_max, z_min, z_max= subps.iloc[ilen, 2:8]
				feat = []
				for i, f in enumerate(ppf):
					if x_min < ppuse_pos[i][0] < x_max and y_min < ppuse_pos[i][1] < y_max and z_min < ppuse_pos[i][2] < z_max:
						feat.append(f)
				das = [[]]
				if len(feat)>5:
					if facts[0] > 0:
						feat = filter_feats_renew2(feat)
						if facts[0] > 1:
							das = [f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')]
							das = unique([x for x in das if das.count(x) > 1])
							feat = filter_feats_both(feat)
							if facts[0] > 4:
								mps = ppm.GetConformer().GetPositions()							
								atms = [[am for i,am in enumerate(ppm.GetAtoms()) if x_min < mps[i][0] < x_max and y_min < mps[i][1] < y_max and z_min < mps[i][2] < z_max]]
								mps = [pm for pm in ppm.GetConformer().GetPositions() if x_min < pm[0] < x_max and y_min < pm[1] < y_max and z_min < pm[2] < z_max]
								mps = [np.stack(mps)]
								if facts[0] ==6:
									vt = feature_AA(rns, [feat], atms, mps, 0, cnrs, inv, [das], fn)
								else:
									vt = feature_AH(rns, [feat], atms, mps, 0, cnrs, inv, [das], fn)
							else:
								vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
						else:
							vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
					else:
						vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
					vts.append(vt)
			dfs = pd.DataFrame(vts).drop_duplicates()
			dfs['pdb'] = pps
		except:
			dfs = pd.DataFrame([[0]*len(rns)*len(cnrs)]).drop_duplicates()
			dfs['pdb'] = pps
		return dfs
	aps = mclapply([[pps, facts, rns, cnrs, fn] for pps in set(allsets.pdb_id)], bulks_pdb, ncpu=100)
	app = pd.concat(aps)
	app.index= list(app.pdb)
	app = app.iloc[:,0:-1]
	app.to_pickle('{}/sub_pockets_{}_20_renew2_filter.pkl'.format(dtip, facts[0]))
	del app

for ii in [11]:
	facts = [ii]
	factory = ChemicalFeatures.BuildFeatureFactory('{}/etc/BaseFeatures{}.fdef'.format(DIR_AP,f_list[facts[0]]))
	fn = factory.GetFeatureFamilies()
	fns = fnss[facts[0]]
	rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
	def bulks_pdb(ssets):
		try:
			pps, facts, rns, cnrs, fn = ssets
			ppm = Chem.MolFromPDBFile('{}/PDB/{}_{}.pdb'.format(dtip,pps[0:4].upper(),pps[4].upper()), removeHs=True)
			ppf = factory.GetFeaturesForMol(ppm)
			ppuse_pos = np.array([ft.GetPos() for ft in ppf])
			subps = ori[ori.iloc[:,1]==pps]
			vts = []
			for ilen in range(len(subps)):
				x_min, x_max,y_min, y_max, z_min, z_max= subps.iloc[ilen, 2:8]
				feat = []
				for i, f in enumerate(ppf):
					if x_min < ppuse_pos[i][0] < x_max and y_min < ppuse_pos[i][1] < y_max and z_min < ppuse_pos[i][2] < z_max:
						feat.append(f)
				das = [[]]
				if len(feat)>5:
					if facts[0] > 0:
						feat = filter_feats_renew_posarm(feat)
						if facts[0] > 1:
							das = [f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')]
							das = unique([x for x in das if das.count(x) > 1])
							posarm = [f.GetAtomIds() for f in feat if (f.GetFamily()=='Aromatic')|(f.GetFamily()=='PosIonizable')]
							posarm = [a for b in posarm for a in b]
							posarm = unique([x for x in posarm if posarm.count(x) > 1])
							feat = filter_feats_dup(feat)
							if facts[0] > 4:
								mps = ppm.GetConformer().GetPositions()							
								atms = [[am for i,am in enumerate(ppm.GetAtoms()) if x_min < mps[i][0] < x_max and y_min < mps[i][1] < y_max and z_min < mps[i][2] < z_max]]
								mps = [pm for pm in ppm.GetConformer().GetPositions() if x_min < pm[0] < x_max and y_min < pm[1] < y_max and z_min < pm[2] < z_max]
								mps = [np.stack(mps)]
								if facts[0] ==6:
									vt = feature_AA(rns, [feat], atms, mps, 0, cnrs, inv, [das], fn)
								else:
									vt = feature_AHP(rns, [feat], atms, mps, 0, cnrs, inv, [das], [posarm],fn)
							else:
								vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
						else:
							vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
					else:
						vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
					vts.append(vt)
			dfs = pd.DataFrame(vts).drop_duplicates()
			dfs['pdb'] = pps
		except:
			dfs = pd.DataFrame([[0]*len(rns)*len(cnrs)]).drop_duplicates()
			dfs['pdb'] = pps
		return dfs
	aps = mclapply([[pps, facts, rns, cnrs, fn] for pps in set(allsets.pdb_id)], bulks_pdb, ncpu=150)
	app = pd.concat(aps)
	app.index= list(app.pdb)
	app = app.iloc[:,0:-1]
	app.to_pickle('{}/sub_pockets_{}_20_renewposarm.pkl'.format(dtip, facts[0]))
	del app


for ii in [11]:
	facts = [ii]
	factory = ChemicalFeatures.BuildFeatureFactory('{}/etc/BaseFeatures{}.fdef'.format(DIR_AP,f_list[facts[0]]))
	fn = factory.GetFeatureFamilies()
	fns = fnss[facts[0]]
	rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
	def bulks_pdb(ssets):
		try:
			pps, facts, rns, cnrs, fn = ssets
			ppm = Chem.MolFromPDBFile('{}/PDB/{}_{}.pdb'.format(dtip,pps[0:4].upper(),pps[4].upper()), removeHs=True)
			ppf = factory.GetFeaturesForMol(ppm)
			ppuse_pos = np.array([ft.GetPos() for ft in ppf])
			subps = ori[ori.iloc[:,1]==pps]
			vts = []
			for ilen in range(len(subps)):
				x_min, x_max,y_min, y_max, z_min, z_max= subps.iloc[ilen, 2:8]
				feat = []
				for i, f in enumerate(ppf):
					if x_min < ppuse_pos[i][0] < x_max and y_min < ppuse_pos[i][1] < y_max and z_min < ppuse_pos[i][2] < z_max:
						feat.append(f)
				das = [[]]
				if len(feat)>5:
					if facts[0] > 0:
						feat = filter_feats_renew3(feat)
						if facts[0] > 1:
							das = [f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')]
							das = unique([x for x in das if das.count(x) > 1])
							feat = filter_feats_both(feat)
							if facts[0] > 4:
								mps = ppm.GetConformer().GetPositions()							
								atms = [[am for i,am in enumerate(ppm.GetAtoms()) if x_min < mps[i][0] < x_max and y_min < mps[i][1] < y_max and z_min < mps[i][2] < z_max]]
								mps = [pm for pm in ppm.GetConformer().GetPositions() if x_min < pm[0] < x_max and y_min < pm[1] < y_max and z_min < pm[2] < z_max]
								mps = [np.stack(mps)]
								if facts[0] ==6:
									vt = feature_AA(rns, [feat], atms, mps, 0, cnrs, inv, [das], fn)
								else:
									vt = feature_AH(rns, [feat], atms, mps, 0, cnrs, inv, [das], fn)
							else:
								vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
						else:
							vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
					else:
						vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
					vts.append(vt)
			dfs = pd.DataFrame(vts).drop_duplicates()
			dfs['pdb'] = pps
		except:
			dfs = pd.DataFrame([[0]*len(rns)*len(cnrs)]).drop_duplicates()
			dfs['pdb'] = pps
		return dfs
	aps = mclapply([[pps, facts, rns, cnrs, fn] for pps in set(allsets.pdb_id)], bulks_pdb, ncpu=120)
	app = pd.concat(aps)
	app.index= list(app.pdb)
	app = app.iloc[:,0:-1]
	app.to_pickle('{}/sub_pockets_{}_10_renew3.pkl'.format(dtip, facts[0]))
	del app


for ii in [11]:
	facts = [ii]
	factory = ChemicalFeatures.BuildFeatureFactory('{}/etc/BaseFeatures{}.fdef'.format(DIR_AP,f_list[facts[0]]))
	fn = factory.GetFeatureFamilies()
	fns = fnss[facts[0]]
	rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
	def bulks_pdb(ssets):
		try:
			pps, facts, rns, cnrs, fn = ssets
			ppm = Chem.MolFromPDBFile('{}/PDB/{}_{}.pdb'.format(dtip,pps[0:4].upper(),pps[4].upper()), removeHs=True)
			ppf = factory.GetFeaturesForMol(ppm)
			ppuse_pos = np.array([ft.GetPos() for ft in ppf])
			subps = pp[pp.iloc[:,1]==pps]
			vts = []
			for ilen in range(len(subps)):
				x_min, x_max,y_min, y_max, z_min, z_max= subps.iloc[ilen, 2:8]
				feat = []
				for i, f in enumerate(ppf):
					if x_min < ppuse_pos[i][0] < x_max and y_min < ppuse_pos[i][1] < y_max and z_min < ppuse_pos[i][2] < z_max:
						feat.append(f)
				das = [[]]
				if len(feat)>5:
					if facts[0] > 0:
						feat = filter_feats_renew2(feat)
						if facts[0] > 1:
							das = [f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')]
							das = unique([x for x in das if das.count(x) > 1])
							# feat = filter_feats_both(feat)
							if facts[0] > 4:
								mps = ppm.GetConformer().GetPositions()							
								atms = [[am for i,am in enumerate(ppm.GetAtoms()) if x_min < mps[i][0] < x_max and y_min < mps[i][1] < y_max and z_min < mps[i][2] < z_max]]
								mps = [pm for pm in ppm.GetConformer().GetPositions() if x_min < pm[0] < x_max and y_min < pm[1] < y_max and z_min < pm[2] < z_max]
								mps = [np.stack(mps)]
								if facts[0] ==6:
									vt = feature_AA(rns, [feat], atms, mps, 0, cnrs, inv, [das], fn)
								else:
									vt = feature_AH(rns, [feat], atms, mps, 0, cnrs, inv, [das], fn)
							else:
								vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
						else:
							vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
					else:
						vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
					vts.append(vt)
			dfs = pd.DataFrame(vts).drop_duplicates()
			dfs['pdb'] = pps
		except:
			dfs = pd.DataFrame([[0]*len(rns)*len(cnrs)]).drop_duplicates()
			dfs['pdb'] = pps
		return dfs
	aps = mclapply([[pps, facts, rns, cnrs, fn] for pps in set(allsets.pdb_id)], bulks_pdb, ncpu=150)
	app = pd.concat(aps)
	app.index= list(app.pdb)
	app = app.iloc[:,0:-1]
	app.to_pickle('{}/sub_pockets_{}_20_renew2_z_no_both_20.pkl'.format(dtip, facts[0]))
	del app


for ii in [11]:
	facts = [ii]
	factory = ChemicalFeatures.BuildFeatureFactory('{}/etc/BaseFeatures{}.fdef'.format(DIR_AP,f_list[facts[0]]))
	fn = factory.GetFeatureFamilies()
	fns = fnss[facts[0]]
	rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
	def bulks_pdb(ssets):
		try:
			pps, facts, rns, cnrs, fn = ssets
			ppm = Chem.MolFromPDBFile('{}/PDB/{}_{}.pdb'.format(dtip,pps[0:4].upper(),pps[5].upper()), removeHs=True)
			ppf = factory.GetFeaturesForMol(ppm)
			ppuse_pos = np.array([ft.GetPos() for ft in ppf])
			subps = ori[ori.iloc[:,1]==pps]
			vts = []
			for ilen in range(len(subps)):
				x_min, x_max,y_min, y_max, z_min, z_max= subps.iloc[ilen, 2:8]
				feat = []
				for i, f in enumerate(ppf):
					if x_min < ppuse_pos[i][0] < x_max and y_min < ppuse_pos[i][1] < y_max and z_min < ppuse_pos[i][2] < z_max:
						feat.append(f)
				das = [[]]
				if len(feat)>5:
					if facts[0] > 0:
						feat = filter_feats_renew2(feat)
						if facts[0] > 1:
							das = [f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')]
							das = unique([x for x in das if das.count(x) > 1])
							# feat = filter_feats_both(feat)
							if facts[0] > 4:
								mps = ppm.GetConformer().GetPositions()							
								atms = [[am for i,am in enumerate(ppm.GetAtoms()) if x_min < mps[i][0] < x_max and y_min < mps[i][1] < y_max and z_min < mps[i][2] < z_max]]
								mps = [pm for pm in ppm.GetConformer().GetPositions() if x_min < pm[0] < x_max and y_min < pm[1] < y_max and z_min < pm[2] < z_max]
								mps = [np.stack(mps)]
								if facts[0] ==6:
									vt = feature_AA(rns, [feat], atms, mps, 0, cnrs, inv, [das], fn)
								else:
									vt = feature_AH(rns, [feat], atms, mps, 0, cnrs, inv, [das], fn)
							else:
								vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
						else:
							vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
					else:
						vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
					vts.append(vt)
			dfs = pd.DataFrame(vts).drop_duplicates()
			dfs['pdb'] = pps
		except:
			dfs = pd.DataFrame([[0]*len(rns)*len(cnrs)]).drop_duplicates()
			dfs['pdb'] = pps
		return dfs
	aps = mclapply([[pps, facts, rns, cnrs, fn] for pps in set(ori.iloc[:,1])], bulks_pdb, ncpu=10)
	app = pd.concat(aps)
	app.index= list(app.pdb)
	app = app.iloc[:,0:-1]
	app.to_pickle('{}/sub_pockets_{}_20_renew2_no_both_20_left.pkl'.format(dtip, facts[0]))
	del app



chemm = '/spstorage/USERS/gina/Project/AP/DB/ChEMBL'
ori = pd.read_csv('/spstorage/USERS/gina/Project/AP/DB/ChEMBL/PDB/pdb_pockets.csv',index_col=0)

chemm = '/spstorage/USERS/gina/Project/AP/DB/ChEMBL'
ori = pd.read_csv('/spstorage/USERS/gina/Project/AP/DB/ChEMBL/PDB/pdb_pockets_10.csv',index_col=0)
ori['PDB'] = [aa.split('_')[0] for aa in list(ori.iloc[:,1])]
ori = ori[ori.PDB.isin(list(uu.To)+list(cc.To))]
pori = pd.read_csv(dtip+'/data/pdb_ptn_loc_gpdb.csv', index_col = 0)
pori['PDB'] = [aa.upper() for aa in list(pori.iloc[:,1])]
pori = pori[pori.PDB.isin(uu.To)]
pori = pori[pori.iloc[:,0]=='sub']


for ii in [11]:
	facts = [ii]
	factory = ChemicalFeatures.BuildFeatureFactory('{}/etc/BaseFeatures{}.fdef'.format(DIR_AP,f_list[facts[0]]))
	fn = factory.GetFeatureFamilies()
	fns = fnss[facts[0]]
	rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
	def bulks_pdb(ssets):
		try:
			pps, facts, rns, cnrs, fn = ssets
			ppm = Chem.MolFromPDBFile('{}/PDB/{}_{}.pdb'.format(chemm,pps[0:4].upper(),pps[5].upper()), removeHs=True)
			ppf = factory.GetFeaturesForMol(ppm)
			ppuse_pos = np.array([ft.GetPos() for ft in ppf])
			subps = ori[ori.iloc[:,1]==pps]
			vts = []
			for ilen in range(len(subps)):
				x_min, x_max,y_min, y_max, z_min, z_max, pck = subps.iloc[ilen, 2:9]
				feat = []
				for i, f in enumerate(ppf):
					if x_min < ppuse_pos[i][0] < x_max and y_min < ppuse_pos[i][1] < y_max and z_min < ppuse_pos[i][2] < z_max:
						feat.append(f)
				das = [[]]
				if len(feat)>5:
					if facts[0] > 0:
						feat = filter_feats_renew2(feat)
						if facts[0] > 1:
							das = [f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')]
							das = unique([x for x in das if das.count(x) > 1])
							# feat = filter_feats_both(feat)
							if facts[0] > 4:
								mps = ppm.GetConformer().GetPositions()							
								atms = [[am for i,am in enumerate(ppm.GetAtoms()) if x_min < mps[i][0] < x_max and y_min < mps[i][1] < y_max and z_min < mps[i][2] < z_max]]
								mps = [pm for pm in ppm.GetConformer().GetPositions() if x_min < pm[0] < x_max and y_min < pm[1] < y_max and z_min < pm[2] < z_max]
								mps = [np.stack(mps)]
								if facts[0] ==6:
									vt = feature_AA(rns, [feat], atms, mps, 0, cnrs, inv, [das], fn)
								else:
									vt = feature_AH(rns, [feat], atms, mps, 0, cnrs, inv, [das], fn)
							else:
								vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
						else:
							vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
					else:
						vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
					vts.append(vt+[pck])
			dfs = pd.DataFrame(vts).drop_duplicates()
			dfs['pdb'] = pps
		except:
			dfs = pd.DataFrame([[0]*len(rns)*len(cnrs)]).drop_duplicates()
			dfs['pdb'] = pps
		return dfs
	aps = mclapply([[pps, facts, rns, cnrs, fn] for pps in set(ori.iloc[:,1])], bulks_pdb, ncpu=62)
	app = pd.concat(aps)
	app.index= list(app.pdb)
	app = app.iloc[:,0:-1]
	app.to_pickle('{}/sub_pockets_{}_20_renew2_no_both_20.pkl'.format(chemm, facts[0]))
	del app

for ii in [11]:
	facts = [ii]
	factory = ChemicalFeatures.BuildFeatureFactory('{}/etc/BaseFeatures{}.fdef'.format(DIR_AP,f_list[facts[0]]))
	fn = factory.GetFeatureFamilies()
	fns = fnss[facts[0]]
	rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
	def bulks_pdb(ssets):
		try:
			pps, facts, rns, cnrs, fn = ssets
			ppm = Chem.MolFromPDBFile('{}/PDB/{}_{}.pdb'.format(chemm,pps[0:4].upper(),pps[5].upper()), removeHs=True)
			ppf = factory.GetFeaturesForMol(ppm)
			ppuse_pos = np.array([ft.GetPos() for ft in ppf])
			subps = ori[ori.iloc[:,1]==pps]
			vts = []
			for ilen in range(len(subps)):
				x_min, x_max,y_min, y_max, z_min, z_max, pck = subps.iloc[ilen, 2:9]
				feat = []
				for i, f in enumerate(ppf):
					if x_min < ppuse_pos[i][0] < x_max and y_min < ppuse_pos[i][1] < y_max and z_min < ppuse_pos[i][2] < z_max:
						feat.append(f)
				das = [[]]
				if len(feat)>5:
					if facts[0] > 0:
						feat = filter_feats_renew2(feat)
						if facts[0] > 1:
							das = [f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')]
							das = unique([x for x in das if das.count(x) > 1])
							# feat = filter_feats_both(feat)
							if facts[0] > 4:
								mps = ppm.GetConformer().GetPositions()							
								atms = [[am for i,am in enumerate(ppm.GetAtoms()) if x_min < mps[i][0] < x_max and y_min < mps[i][1] < y_max and z_min < mps[i][2] < z_max]]
								mps = [pm for pm in ppm.GetConformer().GetPositions() if x_min < pm[0] < x_max and y_min < pm[1] < y_max and z_min < pm[2] < z_max]
								mps = [np.stack(mps)]
								if facts[0] ==6:
									vt = feature_AA(rns, [feat], atms, mps, 0, cnrs, inv, [das], fn)
								else:
									vt = feature_AH(rns, [feat], atms, mps, 0, cnrs, inv, [das], fn)
							else:
								vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
						else:
							vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
					else:
						vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
					vts.append(vt+[pck])
			dfs = pd.DataFrame(vts).drop_duplicates()
			dfs['pdb'] = pps
		except:
			dfs = pd.DataFrame([[0]*len(rns)*len(cnrs)]).drop_duplicates()
			dfs['pdb'] = pps
		return dfs
	aps = mclapply([[pps, facts, rns, cnrs, fn] for pps in set(ori.iloc[:,1])], bulks_pdb, ncpu=150)
	app = pd.concat(aps)
	app.index= list(app.pdb)
	app = app.iloc[:,0:-1]
	app.to_pickle('{}/sub_pockets_{}_20_renew2.pkl'.format(chemm, facts[0]))
	del app




for ii in [0,6,11]:
	facts = [ii]
	factory = ChemicalFeatures.BuildFeatureFactory('{}/etc/BaseFeatures{}.fdef'.format(DIR_AP,f_list[facts[0]]))
	fn = factory.GetFeatureFamilies()
	fns = fnss[facts[0]]
	rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
	def bulks_pdb(ssets):
		try:
			pps, facts, rns, cnrs, fn = ssets
			ppm = Chem.MolFromPDBFile('{}/PDB/{}_{}.pdb'.format(chemm,pps[0:4].upper(),pps[5].upper()), removeHs=True)
			ppf = factory.GetFeaturesForMol(ppm)
			ppuse_pos = np.array([ft.GetPos() for ft in ppf])
			subps = ori[ori.iloc[:,1]==pps]
			vts = []
			for ilen in range(len(subps)):
				x_min, x_max,y_min, y_max, z_min, z_max, pck = subps.iloc[ilen, 2:9]
				feat = []
				for i, f in enumerate(ppf):
					if x_min < ppuse_pos[i][0] < x_max and y_min < ppuse_pos[i][1] < y_max and z_min < ppuse_pos[i][2] < z_max:
						feat.append(f)
				das = [[]]
				if len(feat)>5:
					if facts[0] > 0:
						feat = filter_feats_renew(feat)
						if facts[0] > 1:
							das = [f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')]
							das = unique([x for x in das if das.count(x) > 1])
							if facts[0] > 4:
								mps = ppm.GetConformer().GetPositions()							
								atms = [[am for i,am in enumerate(ppm.GetAtoms()) if x_min < mps[i][0] < x_max and y_min < mps[i][1] < y_max and z_min < mps[i][2] < z_max]]
								mps = [pm for pm in ppm.GetConformer().GetPositions() if x_min < pm[0] < x_max and y_min < pm[1] < y_max and z_min < pm[2] < z_max]
								mps = [np.stack(mps)]
								if facts[0] ==6:
									vt = feature_AA(rns, [feat], atms, mps, 0, cnrs, inv, [das], fn)
								else:
									vt = feature_AH(rns, [feat], atms, mps, 0, cnrs, inv, [das], fn)
							else:
								vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
						else:
							vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
					else:
						vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
					vts.append(vt+[pck])
			dfs = pd.DataFrame(vts).drop_duplicates()
			dfs['pdb'] = pps
		except:
			dfs = pd.DataFrame([[0]*len(rns)*len(cnrs)]).drop_duplicates()
			dfs['pdb'] = pps
		return dfs
	aps = mclapply([[pps, facts, rns, cnrs, fn] for pps in set(ori.iloc[:,1])], bulks_pdb, ncpu=100)
	app = pd.concat(aps)
	app.index= list(app.pdb)
	app = app.iloc[:,0:-1]
	app.to_pickle('{}/sub_pockets_{}_20.pkl'.format(chemm, facts[0]))
	del app


for ii in [11]:
	facts = [ii]
	factory = ChemicalFeatures.BuildFeatureFactory('{}/etc/BaseFeatures{}.fdef'.format(DIR_AP,f_list[facts[0]]))
	fn = factory.GetFeatureFamilies()
	fns = fnss[facts[0]]
	rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
	def bulks_pdb(ssets):
		try:
			pps, facts, rns, cnrs, fn = ssets
			ppm = Chem.MolFromPDBFile(f'{dtip}/covid/{pps}_protein.pdb', removeHs=True)
			ppf = factory.GetFeaturesForMol(ppm)
			ppuse_pos = np.array([ft.GetPos() for ft in ppf])
			subps = ori[ori.iloc[:,1]==pps]
			vts = []
			for ilen in range(len(subps)):
				ldst, lx,ly,lz,x_min, x_max,y_min, y_max, z_min, z_max = subps.iloc[ilen, 2:12]
				feat = []
				for i, f in enumerate(ppf):
					if x_min < ppuse_pos[i][0] < x_max and y_min < ppuse_pos[i][1] < y_max and z_min < ppuse_pos[i][2] < z_max:
						feat.append(f)
				das = [[]]
				if len(feat)>5:
					if facts[0] > 0:
						feat = filter_feats_renew(feat)
						if facts[0] > 1:
							das = [f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')]
							das = unique([x for x in das if das.count(x) > 1])
							if facts[0] > 4:
								mps = ppm.GetConformer().GetPositions()							
								atms = [[am for i,am in enumerate(ppm.GetAtoms()) if x_min < mps[i][0] < x_max and y_min < mps[i][1] < y_max and z_min < mps[i][2] < z_max]]
								mps = [pm for pm in ppm.GetConformer().GetPositions() if x_min < pm[0] < x_max and y_min < pm[1] < y_max and z_min < pm[2] < z_max]
								mps = [np.stack(mps)]
								if facts[0] ==6:
									vt = feature_AA(rns, [feat], atms, mps, 0, cnrs, inv, [das], fn)
								else:
									vt = feature_AH(rns, [feat], atms, mps, 0, cnrs, inv, [das], fn)
							else:
								vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
						else:
							vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
					else:
						vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
					vts.append(vt+[ldst])
			dfs = pd.DataFrame(vts).drop_duplicates()
			dfs['pdb'] = pps
		except:
			dfs = pd.DataFrame([[0]*len(rns)*len(cnrs)]+['-']).drop_duplicates()
			dfs['pdb'] = pps
		return dfs
	aps = mclapply([[pps, facts, rns, cnrs, fn] for pps in c2], bulks_pdb, ncpu=len(c2))
	app = pd.concat(aps)
	app.index= list(app.pdb)
	# app = app.iloc[:,0:-1]
	app.to_pickle('{}/covid/sub_pockets_{}_20.pkl'.format(dtip, facts[0]))
	del app

for ii in [11]:
	facts = [ii]
	factory = ChemicalFeatures.BuildFeatureFactory('{}/etc/BaseFeatures{}.fdef'.format(DIR_AP,f_list[facts[0]]))
	fn = factory.GetFeatureFamilies()
	fns = fnss[facts[0]]
	rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
	def bulks_pdb(ssets):
		try:
			pps, facts, rns, cnrs, fn = ssets
			ppm = Chem.MolFromPDBFile(f'{dtip}/covid/{pps}_protein.pdb', removeHs=True)
			ppf = factory.GetFeaturesForMol(ppm)
			ppuse_pos = np.array([ft.GetPos() for ft in ppf])
			subps = ori[ori.iloc[:,1]==pps]
			vts = []
			for ilen in range(len(subps)):
				ldst, lx,ly,lz,x_min, x_max,y_min, y_max, z_min, z_max = subps.iloc[ilen, 2:12]
				feat = []
				for i, f in enumerate(ppf):
					if x_min < ppuse_pos[i][0] < x_max and y_min < ppuse_pos[i][1] < y_max and z_min < ppuse_pos[i][2] < z_max:
						feat.append(f)
				das = [[]]
				if len(feat)>5:
					if facts[0] > 0:
						feat = filter_feats_renew2(feat)
						if facts[0] > 1:
							das = [f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')]
							das = unique([x for x in das if das.count(x) > 1])
							feat = filter_feats_both(feat)
							if facts[0] > 4:
								mps = ppm.GetConformer().GetPositions()							
								atms = [[am for i,am in enumerate(ppm.GetAtoms()) if x_min < mps[i][0] < x_max and y_min < mps[i][1] < y_max and z_min < mps[i][2] < z_max]]
								mps = [pm for pm in ppm.GetConformer().GetPositions() if x_min < pm[0] < x_max and y_min < pm[1] < y_max and z_min < pm[2] < z_max]
								mps = [np.stack(mps)]
								if facts[0] ==6:
									vt = feature_AA(rns, [feat], atms, mps, 0, cnrs, inv, [das], fn)
								else:
									vt = feature_AH(rns, [feat], atms, mps, 0, cnrs, inv, [das], fn)
							else:
								vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
						else:
							vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
					else:
						vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
					vts.append(vt+[ldst])
			dfs = pd.DataFrame(vts).drop_duplicates()
			dfs['pdb'] = pps
		except:
			dfs = pd.DataFrame([[0]*len(rns)*len(cnrs)]+['-']).drop_duplicates()
			dfs['pdb'] = pps
		return dfs
	aps = mclapply([[pps, facts, rns, cnrs, fn] for pps in list(set(ori.iloc[:,1]))], bulks_pdb, ncpu=31)
	app = pd.concat(aps)
	app.index= list(app.pdb)
	# app = app.iloc[:,0:-1]
	app.to_pickle('{}/covid/sub_pockets_{}_20_renew2_both.pkl'.format(dtip, facts[0]))
	del app


for ii in [11]:
	facts = [ii]
	factory = ChemicalFeatures.BuildFeatureFactory('{}/etc/BaseFeatures{}.fdef'.format(DIR_AP,f_list[facts[0]]))
	fn = factory.GetFeatureFamilies()
	fns = fnss[facts[0]]
	rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
	def bulks_pdb(ssets):
		try:
			pps, facts, rns, cnrs, fn = ssets
			ppm = Chem.MolFromPDBFile(f'{dtip}/covid/{pps}_protein.pdb', removeHs=True)
			ppf = factory.GetFeaturesForMol(ppm)
			ppa = ppm.GetAtoms()
			ppuse_pos = np.array([ft.GetPos() for ft in ppf])
			subps = ori[ori.iloc[:,1]==pps]
			vts = []
			for ilen in range(len(subps)):
				ldst, lx,ly,lz,x_min, x_max,y_min, y_max, z_min, z_max = subps.iloc[ilen, 2:12]
				feats = []
				for i, f in enumerate(ppf):
					if x_min < ppuse_pos[i][0] < x_max and y_min < ppuse_pos[i][1] < y_max and z_min < ppuse_pos[i][2] < z_max:
						feats.append(f)
				feat = []
				for f in feats:
					if len(f.GetAtomIds())!=1:
						feat.append(f)
					if len(f.GetAtomIds())==1:
						pa =  ppa[f.GetAtomIds()[0]]
						pi = pa.GetPDBResidueInfo()
						if pi.GetName().split(' ')[1] not in ['N','O','CA']:
							if pa.GetSymbol()=='C':
								if len([aa for aa in pa.GetNeighbors() if aa.GetSymbol() in ['O','N']])==0:
									feat.append(f)
							else:
								feat.append(f)
					else:
						feat.append(f)
				das = [[]]
				if len(feat)>5:
					if facts[0] > 0:
						feat = filter_feats_renew(feat)
						if facts[0] > 1:
							das = [f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')]
							das = unique([x for x in das if das.count(x) > 1])
							if facts[0] > 4:
								mps = ppm.GetConformer().GetPositions()							
								atms = [[am for i,am in enumerate(ppm.GetAtoms()) if x_min < mps[i][0] < x_max and y_min < mps[i][1] < y_max and z_min < mps[i][2] < z_max]]
								mps = [pm for pm in ppm.GetConformer().GetPositions() if x_min < pm[0] < x_max and y_min < pm[1] < y_max and z_min < pm[2] < z_max]
								mps = [np.stack(mps)]
								if facts[0] ==6:
									vt = feature_AA(rns, [feat], atms, mps, 0, cnrs, inv, [das], fn)
								else:
									vt = feature_AH(rns, [feat], atms, mps, 0, cnrs, inv, [das], fn)
							else:
								vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
						else:
							vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
					else:
						vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
					vts.append(vt+[ldst])
			dfs = pd.DataFrame(vts).drop_duplicates()
			dfs['pdb'] = pps
		except:
			dfs = pd.DataFrame([[0]*len(rns)*len(cnrs)]+['-']).drop_duplicates()
			dfs['pdb'] = pps
		return dfs
	aps = mclapply([[pps, facts, rns, cnrs, fn] for pps in c2], bulks_pdb, ncpu=len(c2))
	app = pd.concat(aps)
	app.index= list(app.pdb)
	# app = app.iloc[:,0:-1]
	app.to_pickle('{}/covid/sub_pockets_{}_20_side.pkl'.format(dtip, facts[0]))
	del app

app = pd.read_pickle('{}/sub_pockets_{}_20_side_new.pkl'.format(dtip, facts[0]))


for ii in [0,6,11]:
	facts = [ii]
	factory = ChemicalFeatures.BuildFeatureFactory('{}/etc/BaseFeatures{}.fdef'.format(DIR_AP,f_list[facts[0]]))
	fn = factory.GetFeatureFamilies()
	fns = fnss[facts[0]]
	rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
	def bulks_pdb(ssets):
		try:
			pps, facts, rns, cnrs, fn = ssets
			ppm = Chem.MolFromPDBFile('{}/DUDE/{}.pdb'.format(dtip,pps), removeHs=True)
			ppf = factory.GetFeaturesForMol(ppm)
			ppuse_pos = np.array([ft.GetPos() for ft in ppf])
			subps = ori[ori.iloc[:,1]==pps]
			vts = []
			for ilen in range(len(subps)):
				x_min, x_max,y_min, y_max, z_min, z_max= subps.iloc[ilen, 2:8]
				feat = []
				for i, f in enumerate(ppf):
					if x_min < ppuse_pos[i][0] < x_max and y_min < ppuse_pos[i][1] < y_max and z_min < ppuse_pos[i][2] < z_max:
						feat.append(f)
				das = [[]]
				if len(feat)>5:
					if facts[0] > 0:
						feat = filter_feats_renew(feat)
						if facts[0] > 1:
							das = [f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')]
							das = unique([x for x in das if das.count(x) > 1])
							if facts[0] > 4:
								mps = ppm.GetConformer().GetPositions()							
								atms = [[am for i,am in enumerate(ppm.GetAtoms()) if x_min < mps[i][0] < x_max and y_min < mps[i][1] < y_max and z_min < mps[i][2] < z_max]]
								mps = [pm for pm in ppm.GetConformer().GetPositions() if x_min < pm[0] < x_max and y_min < pm[1] < y_max and z_min < pm[2] < z_max]
								mps = [np.stack(mps)]
								if facts[0] ==6:
									vt = feature_AA(rns, [feat], atms, mps, 0, cnrs, inv, [das], fn)
								else:
									vt = feature_AH(rns, [feat], atms, mps, 0, cnrs, inv, [das], fn)
							else:
								vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
						else:
							vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
					else:
						vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
					vts.append(vt)
			dfs = pd.DataFrame(vts).drop_duplicates()
			dfs['pdb'] = pps
		except:
			dfs = pd.DataFrame([[0]*len(rns)*len(cnrs)]).drop_duplicates()
			dfs['pdb'] = pps
		return dfs
	aps = mclapply([[pps, facts, rns, cnrs, fn] for pps in set(ori.iloc[:,1])], bulks_pdb, ncpu=102)
	app = pd.concat(aps)
	app.index= list(app.pdb)
	app = app.iloc[:,0:-1]
	app.to_pickle('{}/DUDE/sub_pockets_{}_20_15space.pkl'.format(dtip, facts[0]))
	del app

for ii in [11]:
	facts = [ii]
	factory = ChemicalFeatures.BuildFeatureFactory('{}/etc/BaseFeatures{}.fdef'.format(DIR_AP,f_list[facts[0]]))
	fn = factory.GetFeatureFamilies()
	fns = fnss[facts[0]]
	rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
	def bulks_pdb(ssets):
		try:
			pps, facts, rns, cnrs, fn = ssets
			ppm = Chem.MolFromPDBFile('{}/DUDE/{}.pdb'.format(dtip,pps), removeHs=True)
			ppf = factory.GetFeaturesForMol(ppm)
			ppuse_pos = np.array([ft.GetPos() for ft in ppf])
			subps = ori[ori.iloc[:,1]==pps]
			vts = []
			for ilen in range(len(subps)):
				x_min, x_max,y_min, y_max, z_min, z_max= subps.iloc[ilen, 2:8]
				feats = []
				for i, f in enumerate(ppf):
					if x_min < ppuse_pos[i][0] < x_max and y_min < ppuse_pos[i][1] < y_max and z_min < ppuse_pos[i][2] < z_max:
						feats.append(f)
				feat = []
				for f in feats:
					if len(f.GetAtomIds())==1:
						pa =  ppa[f.GetAtomIds()[0]]
						pi = pa.GetPDBResidueInfo()
						if pi.GetName().split(' ')[1] not in ['N','O','CA']:
							if pa.GetSymbol()=='C':
								if len([aa for aa in pa.GetNeighbors() if aa.GetSymbol() in ['O','N']])==0:
									feat.append(f)
							else:
								feat.append(f)
					else:
						feat.append(f)
				das = [[]]
				if len(feat)>5:
					if facts[0] > 0:
						feat = filter_feats_renew2(feat)
						if facts[0] > 1:
							das = [f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')]
							das = unique([x for x in das if das.count(x) > 1])
							if facts[0] > 4:
								mps = ppm.GetConformer().GetPositions()							
								atms = [[am for i,am in enumerate(ppm.GetAtoms()) if x_min < mps[i][0] < x_max and y_min < mps[i][1] < y_max and z_min < mps[i][2] < z_max]]
								mps = [pm for pm in ppm.GetConformer().GetPositions() if x_min < pm[0] < x_max and y_min < pm[1] < y_max and z_min < pm[2] < z_max]
								mps = [np.stack(mps)]
								if facts[0] ==6:
									vt = feature_AA(rns, [feat], atms, mps, 0, cnrs, inv, [das], fn)
								else:
									vt = feature_AH(rns, [feat], atms, mps, 0, cnrs, inv, [das], fn)
							else:
								vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
						else:
							vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
					else:
						vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
					vts.append(vt)
			dfs = pd.DataFrame(vts).drop_duplicates()
			dfs['pdb'] = pps
		except:
			dfs = pd.DataFrame([[0]*len(rns)*len(cnrs)]).drop_duplicates()
			dfs['pdb'] = pps
		return dfs
	aps = mclapply([[pps, facts, rns, cnrs, fn] for pps in set(ori.iloc[:,1])], bulks_pdb, ncpu=100)
	app = pd.concat(aps)
	app.index= list(app.pdb)
	app = app.iloc[:,0:-1]
	app.to_pickle('{}/sub_pockets_{}_20_side_renew2.pkl'.format(dtip, facts[0]))
	del app



for ii in [0,6,11]:
	facts = [ii]
	factory = ChemicalFeatures.BuildFeatureFactory('{}/etc/BaseFeatures{}.fdef'.format(DIR_AP,f_list[facts[0]]))
	fn = factory.GetFeatureFamilies()
	fns = fnss[facts[0]]
	rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
	def bulks_pdb(ssets):
		try:
			pps, facts, rns, cnrs, fn = ssets
			ppm = Chem.MolFromPDBFile('{}/DUDE/{}.pdb'.format(dtip,pps), removeHs=True)
			ppf = factory.GetFeaturesForMol(ppm)
			ppuse_pos = np.array([ft.GetPos() for ft in ppf])
			subps = ori[ori.iloc[:,1]==pps]
			vts = []
			for ilen in range(len(subps)):
				x_min, x_max,y_min, y_max, z_min, z_max= subps.iloc[ilen, 2:8]
				feats = []
				for i, f in enumerate(ppf):
					if x_min < ppuse_pos[i][0] < x_max and y_min < ppuse_pos[i][1] < y_max and z_min < ppuse_pos[i][2] < z_max:
						feats.append(f)
				feat = []
				for f in feats:
					if len(f.GetAtomIds())==1:
						pa =  ppa[f.GetAtomIds()[0]]
						pi = pa.GetPDBResidueInfo()
						if pi.GetName().split(' ')[1] not in ['N','O','CA']:
							if pa.GetSymbol()=='C':
								if len([aa for aa in pa.GetNeighbors() if aa.GetSymbol() in ['O','N']])==0:
									feat.append(f)
							else:
								feat.append(f)
					else:
						feat.append(f)
				das = [[]]
				if len(feat)>5:
					if facts[0] > 0:
						feat = filter_feats_renew(feat)
						if facts[0] > 1:
							das = [f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')]
							das = unique([x for x in das if das.count(x) > 1])
							if facts[0] > 4:
								mps = ppm.GetConformer().GetPositions()							
								atms = [[am for i,am in enumerate(ppm.GetAtoms()) if x_min < mps[i][0] < x_max and y_min < mps[i][1] < y_max and z_min < mps[i][2] < z_max]]
								mps = [pm for pm in ppm.GetConformer().GetPositions() if x_min < pm[0] < x_max and y_min < pm[1] < y_max and z_min < pm[2] < z_max]
								mps = [np.stack(mps)]
								if facts[0] ==6:
									vt = feature_AA(rns, [feat], atms, mps, 0, cnrs, inv, [das], fn)
								else:
									vt = feature_AH(rns, [feat], atms, mps, 0, cnrs, inv, [das], fn)
							else:
								vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
						else:
							vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
					else:
						vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
					vts.append(vt)
			dfs = pd.DataFrame(vts).drop_duplicates()
			dfs['pdb'] = pps
		except:
			dfs = pd.DataFrame([[0]*len(rns)*len(cnrs)]).drop_duplicates()
			dfs['pdb'] = pps
		return dfs
	aps = mclapply([[pps, facts, rns, cnrs, fn] for pps in set(ori.iloc[:,1])], bulks_pdb, ncpu=102)
	app = pd.concat(aps)
	app.index= list(app.pdb)
	app = app.iloc[:,0:-1]
	app.to_pickle('{}/DUDE/sub_pockets_{}_20_side.pkl'.format(dtip, facts[0]))
	del app



#with open("train_new_dude_balanced_all2_active.pkl", "rb") as fp:
#    active = pickle.load(fp)

#with open("train_new_dude_balanced_all2_decoy.pkl", "rb") as fp:
#    inactive = pickle.load(fp)

with open("/spstorage/USERS/gina/Project/AP/DB/DUDE/train_new_dude_balanced_all2_active.pkl", "rb") as fp:
   active = pickle.load(fp)

with open("/spstorage/USERS/gina/Project/AP/DB/DUDE/train_new_dude_balanced_all2_decoy.pkl", "rb") as fp:
   inactive = pickle.load(fp)



with open("/spstorage/USERS/gina/Project/AP/DB/DUDE/test_new_dude_all2_active_none_pdb.pkl", "rb") as fp:
   active = pickle.load(fp)

with open("/spstorage/USERS/gina/Project/AP/DB/DUDE/test_new_dude_all2_decoy_none_pdb.pkl", "rb") as fp:
   inactive = pickle.load(fp)

ds_test = active + inactive
random.shuffle(ds_test)


X_test = [i[0] for i in ds_test]
y_test = [i[1][0] for i in ds_test]





random_selector = np.random.randint(len(inactive) - len(active))
a = int(len(inactive) / (len(active)))
ds = active + inactive

random.shuffle(ds)

X = [i[0] for i in ds]
y = [i[1][0] for i in ds]

model = DTITAG()


def fwd_pass(X, y, train=False):
    if train:
        model.zero_grad()
    out = []
    for item in X:
        x = [0, 0]
        x[0] = item[0].to(device)
        x[1] = item[1].to(device)
        out.append(model(x))
        del x
    out = th.stack(out, 0).view(-1, 1).to(device)
    y = th.Tensor(y).view(-1, 1).to(device)
    loss = criterion(out, y)
    matches = [th.round(i) == th.round(j) for i, j in zip(out, y)]
    acc = matches.count(True) / len(matches)
    if train:
        loss.backward()
        optimizer.step()
    return acc, loss, out


def train(net):
    EPOCHS = 30
    BATCH_SIZE = 80
    with open("model_dude_fold.log", "a") as f:
        for epoch in range(EPOCHS):
            losses = []
            accs = []
            with tqdm(range(0, len(X), BATCH_SIZE)) as tepoch:
                for i in tepoch:
                    tepoch.set_description(f"Epoch {epoch + 1}")
                    try:
                        batch_X = X[i: i+BATCH_SIZE]
                        batch_y = y[i: i+BATCH_SIZE]
                    except:
                        gc.collect()
                        continue
                    acc, loss, _ = fwd_pass(batch_X, batch_y, train=True)
                    losses.append(loss.item())
                    accs.append(acc)
                    acc_mean = np.array(accs).mean()
                    loss_mean = np.array(losses).mean()
                    tepoch.set_postfix(loss=loss_mean, accuracy=100. * acc_mean)
                    if i % 100000 == 0:
                        test_func(model, y_test, X_test)
                        # print(f'Average Loss: {val_loss}')
                        # print(f'Average Accuracy: {val_acc}')
                        f.write(
                            f"{MODEL_NAME},{round(time.time(), 3)},{round(float(acc), 2)},{round(float(loss), 4)}\n")
                scheduler.step()
            print(f'Average Loss: {np.array(losses).mean()}')
            print(f'Average Accuracy: {np.array(accs).mean()}')
            #dt = time.strftime("%Y_%m_%d-%H_%M_%S")
            #fn = "without_batching" + str(dt) + str("-") + \
            #     str(epoch) + "_checkpoint.pt"
            #info_dict = {
            #    'epoch': epoch,
            #    'net_state': model.state_dict(),
            #    'optimizer_state': optimizer.state_dict()
            #}
            #th.save(info_dict, fn)



pd.DataFrame(np.array(list(vvt.max())[0:550]).reshape([55,10]))

ft=itertools.combinations(fts, 2)
max([cal_dst(f) for f in ft])

for ii in [0,6,11]:
	facts = [ii]
	factory = ChemicalFeatures.BuildFeatureFactory('{}/etc/BaseFeatures{}.fdef'.format(DIR_AP,f_list[facts[0]]))
	fn = factory.GetFeatureFamilies()
	fns = fnss[facts[0]]
	rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
	def bulks_pdb(ssets):
		try:
			pps, facts, rns, cnrs, fn = ssets
			ppm = Chem.MolFromPDBFile('{}/PDB/{}_{}.pdb'.format(dtip,pps[0:4].upper(),pps[4].upper()), removeHs=True)
			ppf = factory.GetFeaturesForMol(ppm)
			ppuse_pos = np.array([ft.GetPos() for ft in ppf])
			subps = ori[ori.iloc[:,1]==pps]
			vts = []
			for ilen in range(len(subps)):
				x_min, x_max,y_min, y_max, z_min, z_max= subps.iloc[ilen, 2:8]
				feats = []
				for i, f in enumerate(ppf):
					if x_min < ppuse_pos[i][0] < x_max and y_min < ppuse_pos[i][1] < y_max and z_min < ppuse_pos[i][2] < z_max:
						feats.append(f)
				feat = []
				for f in feats:
					if len(f.GetAtomIds())==1:
						pa =  ppa[f.GetAtomIds()[0]]
						pi = pa.GetPDBResidueInfo()
						if pi.GetName().split(' ')[1] not in ['N','O','CA']:
							if pa.GetSymbol()=='C':
								if len([aa for aa in pa.GetNeighbors() if aa.GetSymbol() in ['O','N']])==0:
									feat.append(f)
							else:
								feat.append(f)
					else:
						feat.append(f)
				das = [[]]
				if len(feat)>5:
					if facts[0] > 0:
						feat = filter_feats_renew(feat)
						if facts[0] > 1:
							das = [f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')]
							das = unique([x for x in das if das.count(x) > 1])
							if facts[0] > 4:
								mps = ppm.GetConformer().GetPositions()							
								atms = [[am for i,am in enumerate(ppm.GetAtoms()) if x_min < mps[i][0] < x_max and y_min < mps[i][1] < y_max and z_min < mps[i][2] < z_max]]
								mps = [pm for pm in ppm.GetConformer().GetPositions() if x_min < pm[0] < x_max and y_min < pm[1] < y_max and z_min < pm[2] < z_max]
								mps = [np.stack(mps)]
								if facts[0] ==6:
									vt = feature_AA(rns, [feat], atms, mps, 0, cnrs, inv, [das], fn)
								else:
									vt = feature_AH(rns, [feat], atms, mps, 0, cnrs, inv, [das], fn)
							else:
								vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
						else:
							vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
					else:
						vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
					vts.append(vt)
			dfs = pd.DataFrame(vts).drop_duplicates()
			dfs['pdb'] = pps
		except:
			dfs = pd.DataFrame([[0]*len(rns)*len(cnrs)]).drop_duplicates()
			dfs['pdb'] = pps
		return dfs
	aps = mclapply([[pps, facts, rns, cnrs, fn] for pps in set(ori.iloc[:,1])], bulks_pdb, ncpu=104)
	app = pd.concat(aps)
	app.index= list(app.pdb)
	app = app.iloc[:,0:-1]
	app.to_pickle('{}/sub_pockets_{}_20_side.pkl'.format(dtip, facts[0]))
	del app


for ii in [0,6,11]:
	facts = [ii]
	factory = ChemicalFeatures.BuildFeatureFactory('{}/etc/BaseFeatures{}.fdef'.format(DIR_AP,f_list[facts[0]]))
	fn = factory.GetFeatureFamilies()
	fns = fnss[facts[0]]
	rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
	def bulks_pdb(ssets):
		try:
			pps, facts, rns, cnrs, fn = ssets
			ppm = Chem.MolFromPDBFile('{}/PDB/{}_{}.pdb'.format(dtip,pps[0:4].upper(),pps[4].upper()), removeHs=True)
			ppf = factory.GetFeaturesForMol(ppm)
			ppuse_pos = np.array([ft.GetPos() for ft in ppf])
			subps = pp[pp.iloc[:,1]==pps]
			vts = []
			for ilen in range(len(subps)):
				x_min, x_max,y_min, y_max, z_min, z_max= subps.iloc[ilen, 2:8]
				feats = []
				for i, f in enumerate(ppf):
					if x_min < ppuse_pos[i][0] < x_max and y_min < ppuse_pos[i][1] < y_max and z_min < ppuse_pos[i][2] < z_max:
						feats.append(f)
				feat = []
				for f in feats:
					if len(f.GetAtomIds())==1:
						pa =  ppa[f.GetAtomIds()[0]]
						pi = pa.GetPDBResidueInfo()
						if pi.GetName().split(' ')[1] not in ['N','O','CA']:
							if pa.GetSymbol()=='C':
								if len([aa for aa in pa.GetNeighbors() if aa.GetSymbol() in ['O','N']])==0:
									feat.append(f)
							else:
								feat.append(f)
					else:
						feat.append(f)
				das = [[]]
				if len(feat)>5:
					if facts[0] > 0:
						feat = filter_feats_renew(feat)
						if facts[0] > 1:
							das = [f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')]
							das = unique([x for x in das if das.count(x) > 1])
							if facts[0] > 4:
								mps = ppm.GetConformer().GetPositions()							
								atms = [[am for i,am in enumerate(ppm.GetAtoms()) if x_min < mps[i][0] < x_max and y_min < mps[i][1] < y_max and z_min < mps[i][2] < z_max]]
								mps = [pm for pm in ppm.GetConformer().GetPositions() if x_min < pm[0] < x_max and y_min < pm[1] < y_max and z_min < pm[2] < z_max]
								mps = [np.stack(mps)]
								if facts[0] ==6:
									vt = feature_AA(rns, [feat], atms, mps, 0, cnrs, inv, [das], fn)
								else:
									vt = feature_AH(rns, [feat], atms, mps, 0, cnrs, inv, [das], fn)
							else:
								vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
						else:
							vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
					else:
						vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
					vts.append(vt)
			dfs = pd.DataFrame(vts).drop_duplicates()
			dfs['pdb'] = pps
		except:
			dfs = pd.DataFrame([[0]*len(rns)*len(cnrs)]).drop_duplicates()
			dfs['pdb'] = pps
		return dfs
	aps = mclapply([[pps, facts, rns, cnrs, fn] for pps in set(allsets.pdb_id)], bulks_pdb, ncpu=100)
	app = pd.concat(aps)
	app.index= list(app.pdb)
	app = app.iloc[:,0:-1]
	app.to_pickle('{}/sub_pockets_{}_20_side_new.pkl'.format(dtip, facts[0]))
	del app



for ii in [11,6,0]:
	print(ii)
	facts = [ii]
	factory = ChemicalFeatures.BuildFeatureFactory('{}/etc/BaseFeatures{}.fdef'.format(DIR_AP,f_list[facts[0]]))
	fn = factory.GetFeatureFamilies()
	fns = fnss[facts[0]]
	rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
	def bulks_pdb(ssets):
		try:
			pps, facts, rns, cnrs, fn = ssets
			ppm = Chem.MolFromPDBFile('{}/human/{}.pdb'.format(dtip,pps), removeHs=True)
			ppf = factory.GetFeaturesForMol(ppm)
			ppuse_pos = np.array([ft.GetPos() for ft in ppf])
			subps = ori[ori.iloc[:,1]==pps]
			vts = []
			for ilen in range(len(subps)):
				x_min, x_max,y_min, y_max, z_min, z_max= subps.iloc[ilen, 2:8]
				feats = []
				for i, f in enumerate(ppf):
					if x_min < ppuse_pos[i][0] < x_max and y_min < ppuse_pos[i][1] < y_max and z_min < ppuse_pos[i][2] < z_max:
						feats.append(f)
				feat = []
				for f in feats:
					if len(f.GetAtomIds())==1:
						ppa = ppm.GetAtoms()
						pa =  ppa[f.GetAtomIds()[0]]
						pi = pa.GetPDBResidueInfo()
						if len(pi.GetName().split(' '))==1:
							feat.append(f)
						else:
							if pi.GetName().split(' ')[1] not in ['N','O','CA']:
								if pa.GetSymbol()=='C':
									if len([aa for aa in pa.GetNeighbors() if aa.GetSymbol() in ['O','N']])==0:
										feat.append(f)
								else:
									feat.append(f)
					else:
						feat.append(f)
				das = [[]]
				if len(feat)>5:
					if facts[0] > 0:
						feat = filter_feats_renew(feat)
						if facts[0] > 1:
							das = [f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')]
							das = unique([x for x in das if das.count(x) > 1])
							if facts[0] > 4:
								mps = ppm.GetConformer().GetPositions()							
								atms = [[am for i,am in enumerate(ppm.GetAtoms()) if x_min < mps[i][0] < x_max and y_min < mps[i][1] < y_max and z_min < mps[i][2] < z_max]]
								mps = [pm for pm in ppm.GetConformer().GetPositions() if x_min < pm[0] < x_max and y_min < pm[1] < y_max and z_min < pm[2] < z_max]
								mps = [np.stack(mps)]
								if facts[0] ==6:
									vt = feature_AA(rns, [feat], atms, mps, 0, cnrs, inv, [das], fn)
								else:
									vt = feature_AH(rns, [feat], atms, mps, 0, cnrs, inv, [das], fn)
							else:
								vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
						else:
							vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
					else:
						vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
					vts.append(vt)
			dfs = pd.DataFrame(vts).drop_duplicates()
			dfs['pdb'] = pps
		except:
			dfs = pd.DataFrame([[0]*len(rns)*len(cnrs)]).drop_duplicates()
			dfs['pdb'] = pps
		return dfs
	aps = mclapply([[pps, facts, rns, cnrs, fn] for pps in set(ori.iloc[:,1])], bulks_pdb, ncpu=150)
	app = pd.concat(aps)
	app.index= list(app.pdb)
	app = app.iloc[:,0:-1]
	app.to_pickle('{}/human/sub_pockets_{}_20_3.pkl'.format(dtip, facts[0]))
	del app



app = pd.read_pickle('{}/human/sub_pockets_{}_20_2.pkl'.format(dtip, '0'))
app = pd.concat([pd.read_pickle('{}/human/sub_pockets_{}_20.pkl'.format(dtip, '0')),app])
ppock = pock[(pock.T.sum()>10) & (pock.T.sum()<100)]


dtip = '/spstorage/USERS/gina/Project/AP/DB/InterpretableDTIP'


ori = pd.read_csv(dtip+'/human/human_pocket_loc.csv', index_col = 0)

for ii in [0,6,11]:
	facts = [ii]
	factory = ChemicalFeatures.BuildFeatureFactory('{}/etc/BaseFeatures{}.fdef'.format(DIR_AP,f_list[facts[0]]))
	fn = factory.GetFeatureFamilies()
	fns = fnss[facts[0]]
	rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
	def bulks_pdb(ssets):
		try:
			pps, facts, rns, cnrs, fn = ssets
			ppm = Chem.MolFromPDBFile('{}/human/{}.pdb'.format(dtip,pps), removeHs=True)
			ppf = factory.GetFeaturesForMol(ppm)
			ppuse_pos = np.array([ft.GetPos() for ft in ppf])
			subps = ori[ori.iloc[:,1]==pps]
			vts = []
			for ilen in range(len(subps)):
				x_min, x_max,y_min, y_max, z_min, z_max= subps.iloc[ilen, 2:8]
				feats = []
				for i, f in enumerate(ppf):
					if x_min < ppuse_pos[i][0] < x_max and y_min < ppuse_pos[i][1] < y_max and z_min < ppuse_pos[i][2] < z_max:
						feats.append(f)
				feat = []
				for f in feats:
					if len(f.GetAtomIds())==1:
						ppa = ppm.GetAtoms()
						pa =  ppa[f.GetAtomIds()[0]]
						pi = pa.GetPDBResidueInfo()
						if pi.GetName().split(' ')[1] not in ['N','O','CA']:
							if pa.GetSymbol()=='C':
								if len([aa for aa in pa.GetNeighbors() if aa.GetSymbol() in ['O','N']])==0:
									feat.append(f)
							else:
								feat.append(f)
					else:
						feat.append(f)
				das = [[]]
				if len(feat)>5:
					if facts[0] > 0:
						feat = filter_feats_renew(feat)
						if facts[0] > 1:
							das = [f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')]
							das = unique([x for x in das if das.count(x) > 1])
							if facts[0] > 4:
								mps = ppm.GetConformer().GetPositions()							
								atms = [[am for i,am in enumerate(ppm.GetAtoms()) if x_min < mps[i][0] < x_max and y_min < mps[i][1] < y_max and z_min < mps[i][2] < z_max]]
								mps = [pm for pm in ppm.GetConformer().GetPositions() if x_min < pm[0] < x_max and y_min < pm[1] < y_max and z_min < pm[2] < z_max]
								mps = [np.stack(mps)]
								if facts[0] ==6:
									vt = feature_AA(rns, [feat], atms, mps, 0, cnrs, inv, [das], fn)
								else:
									vt = feature_AH(rns, [feat], atms, mps, 0, cnrs, inv, [das], fn)
							else:
								vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
						else:
							vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
					else:
						vt = pairLr_dist(feat, rns, cnrs, inv, das, fn)
					vts.append(vt)
			dfs = pd.DataFrame(vts).drop_duplicates()
			dfs['pdb'] = pps
		except:
			dfs = pd.DataFrame([[0]*len(rns)*len(cnrs)]).drop_duplicates()
			dfs['pdb'] = pps
		return dfs
	aps = mclapply([[pps, facts, rns, cnrs, fn] for pps in set(ori.iloc[:,1])], bulks_pdb, ncpu=102)
	app = pd.concat(aps)
	app.index= list(app.pdb)
	app = app.iloc[:,0:-1]
	app.to_pickle('{}/human/sub_pockets_{}_20.pkl'.format(dtip, facts[0]))
	del app
