


types= ['fp','f6_ratio16','f7_ratio16','f6_gap40','f6_ratio16_filtered','f6_ratio8']

tag_df= pd.read_csv(DIR_AP+'/target50_mc_chembl.csv',index_col=0)

confs = pd.read_pickle(DIR_AP+'/APMs/apm_'+tps+'_conf.pkl') for tps in types[1:]
ltss = [pd.read_pickle(DIR_AP+'/APMs/apm_'+tps+'.pkl') for tps in types]


fps = pd.read_pickle(DIR_AP+'/APMs/apm_'+'fp'+'.pkl')
ltss = [pd.read_pickle(DIR_AP+'/APMs/apm_'+tps+'.pkl') for tps in types]

lts = pd.read_pickle(DIR_AP+'/APMs/apm_'+'f6_ratio16'+'.pkl')


sub_int = all_int[all_int.etid.isin(aa)]


fps[fps.index.isin(sub_int.cid)]


etids = unique(sub_int.etid)
abc = []
for i in range(len(etids)):
	sp = sub_int[(sub_int.etid==etids[i]) & (sub_int.standard_value<=1000)]
	len(sp)
	len(fps[fps.index.isin(sp.cid)])
	if len(fps[fps.index.isin(sp.cid)])>50:abc.append(etids[i])


# pca 
from sklearn.decomposition import PCA
rcms = random.sample(cms, 100)

for et in abc:
	cid_set = [] 
	for ecut in [10000,1000]:
		sp = sub_int[(sub_int.etid==et) & (sub_int.standard_value<=ecut)]
		cid_set.append(unique(sp.cid))
        
        
        
	for eecut in [20000,50000]:
		sp = sub_int[(sub_int.etid==et) & (sub_int.standard_value>eecut)]
		cid_set.append(unique(sp.cid))
        
cmaps = [plt.cm.get_cmap('Greens'),plt.cm.get_cmap('Greens'),plt.cm.get_cmap('Reds'),plt.cm.get_cmap('Reds')]
plt.subplot(1,2,1)
all_cid = list(itertools.chain(*cid_set)) +rcms
al_fps = fps[fps.index.isin(all_cid)]
pca = PCA(n_components=2) # 주성분을 몇개로 할지 결정
printcipalComponents = pca.fit_transform(al_fps)
principalDf = pd.DataFrame(data=printcipalComponents, columns = ['PC1', 'PC2'])
principalDf.index = al_fps.index
rdf = principalDf[principalDf.index.isin(rcms)]
tdf = principalDf[~principalDf.index.isin(rcms)]	
plt.scatter(rdf.loc[:, 'PC1'], rdf.loc[:, 'PC2'], c = 'black')
for cids, cm in zip(cid_set,cmaps):
    subpca = tdf[tdf.index.isin(cids)]
    svs = [np.mean(sub_int[sub_int.cid==ec]['standard_value']) for ec in list(subpca.index)]
    plt.scatter(subpca.loc[:, 'PC1']
               , subpca.loc[:, 'PC2']
               , c = svs
               , cmap = cm
               , s = 20)	


cbb = list(principalDf[(principalDf.PC1>3)&(principalDf.PC2>2)].index)

plt.subplot(1,2,2)
all_cid = list(itertools.chain(*cid_set)) +rcms
al_fps = lts[lts.index.isin(all_cid)]
pca = PCA(n_components=2) # 주성분을 몇개로 할지 결정
printcipalComponents = pca.fit_transform(al_fps)
principalDf = pd.DataFrame(data=printcipalComponents, columns = ['PC1', 'PC2'])
principalDf.index = al_fps.index
rdf = principalDf[principalDf.index.isin(rcms)]
tdf = principalDf[~principalDf.index.isin(rcms)]	
plt.scatter(rdf.loc[:, 'PC1'], rdf.loc[:, 'PC2'], c = 'black')
for cids, cm in zip(cid_set,cmaps):
    subpca = tdf[tdf.index.isin(cids)]
    svs = [np.mean(sub_int[sub_int.cid==ec]['standard_value']) for ec in list(subpca.index)]
    plt.scatter(subpca.loc[:, 'PC1']
               , subpca.loc[:, 'PC2']
               , c = svs
               , cmap = cm
               , s = 20)	



subpca = tdf[tdf.index.isin(cbb)]
svs = [np.mean(sub_int[sub_int.cid==ec]['standard_value']) for ec in list(subpca.index)]
plt.scatter(subpca.loc[:, 'PC1']
           , subpca.loc[:, 'PC2']
           , c = svs
           , cmap = plt.cm.get_cmap('Blues')
           , s = 20)	


plt.show()

print(facts, pca.explained_variance_ratio_, round(sum(pca.explained_variance_ratio_),2))
# pca.components_
pca.components_[1].argmax()






subpcas = pd.DataFrame()
for beta, color in zip(betas,colors):
    subpca = principalDf[principalDf.index.isin(pock[[pidx.split('_')[0] in list(betalm[betalm.Entry.isin(beta.Entry)].PDB) for pidx in list(pock.index)]].index)]
    subpcas = pd.concat([subpcas, subpca])
    # ttx = list(rgidx.loc[list(subpca.index),:][mms[idx]])
    # dds = [et.split('_')[0] for et in list(df.loc[list(subpca.index),:]['PDB'])]
    # for x,y,tt in zip(subpca.loc[:, 'PC1'],subpca.loc[:, 'PC2'],dds):
    #     if tt in list(classdfs.PDB):
    #         tx = classdfs[classdfs.PDB==tt]['name'].item()
    #         ttxs.append(plt.text(x,y,tx))
    # for x,y,tt in zip(subpca.loc[:, 'PC1'],subpca.loc[:, 'PC2'],list(subpca.index)):
    #     plt.text(x,y,tt)  		
    mttx = [float(tt)/mmi[idx] for tt in ttx]
    plt.scatter(subpca.loc[:, 'PC1']
               , subpca.loc[:, 'PC2']
               , c = color
               # ,c = mttx
               # , cmap = cm
               , s = 20)


plt.xlabel(pcs[0])
plt.ylabel(pcs[1])
plt.legend(betat)
plt.show()

	plt.subplot(3,1,tidx+1)
	plt.plot(list(dfdf.qq), list(dfdf.aa))
	plt.title(subdirs[nns]+'; corr. :' + str(round(correlation,2))+ ' ' + ttid[tidx])

plt.show()


ccids =unique(sub_int.cid)
abb='\n'.join([str(int(cd)) for cd in ccids if not np.isnan(cd)])
with open(DIR_AP+'/cid_list_chem98.txt','wt') as f:
	f.write(abb)


pcids = []
for pk in pkls:
	pcids.extend(list(pk[pk.etid.isin(cc)].ecid))


pcids = unique(pcids)

ecid2inck = select_in('DTIMAP_COMPOUND.inchikey2ecid', '*', 'ecid',unique(pcids))

incks=pd.concat([ecid2inck,test_ecid2inck])
ccc= ','.join(list(incks.inchikey))
with open(DIR_AP+'/pkls_inchk.txt','wt') as f:f.write(ccc)

inchik2cid = pd.read_csv(DIR_AP+'/etc2/inchik2cid.txt',header=None, sep='\t')

ncid = list(set(inchi2cid.iloc[:,1])-set(ccids))
abb='\n'.join([str(int(cd)) for cd in ncid if not np.isnan(cd)])
with open(DIR_AP+'/ncids_dti98.txt','wt') as f:f.write(ccc)

for i in range(31):
	abb='\n'.join([str(int(cd)) for cd in ncid[13000+i*2000:13000+(i+1)*2000] if not np.isnan(cd)])
	with open(DIR_AP+'/ncids_dti98_'+str(i)+'.txt','wt') as f:f.write(abb)

with open(DIR_AP+'/final/dti_sdf/ncids_dti98.txt','r') as f:abb= f.readlines()

with open(DIR_AP+'/ncids_dti98.txt','wt') as f:f.write(abb)

rfacts = ['PUBCHEM_COMPOUND_CID','PUBCHEM_EFFECTIVE_ROTOR_COUNT','PUBCHEM_MMFF94_ENERGY','PUBCHEM_CONFORMER_RMSD','PUBCHEM_HEAVY_ATOM_COUNT']
def rigidity_check(sdfs):
	mm = Chem.SDMolSupplier(sdfs, removeHs=False) 
	mn = [[m.GetPropsAsDict()[rfs] for rfs in rfacts] for m in mm if len(set(rfacts)-set(m.GetPropsAsDict().keys()))==0]
	vdf = pd.DataFrame(mn)
	vdf.columns = ['cid','er','mmff','rmsd', 'hac']
	return vdf	


vdf = mclapply(sdfs, rigidity_check, ncpu=len(sdfs))
app = pd.concat(vdf)



app.to_pickle('{}/final/dataset/rigids.pkl'.format(DIR_AP))


c2 = DIR_AP+'/final2/dataset/sdfs/'
sdfs = glob.glob(c2+'*3.sdf')

split_number= 5000 # (number of molecules in each file--delete this line in bracket)
number_of_sdfs = split_number
i=0
j=0
f2=open('/spstorage/USERS/gina/Project/AP/final2/dataset/sdfs/splits/cid_'+str(j)+'.sdf','w')
for sdf in sdfs:
	for line in open(sdf):
		_ = f2.write(line)
		if line[:4] == "$$$$":
			i+=1
		if i > number_of_sdfs:
			number_of_sdfs += split_number 
			f2.close()
			j+=1
			f2=open('/spstorage/USERS/gina/Project/AP/final2/dataset/sdfs/splits/cid_'+str(j)+'.sdf','w')
			print(j)

f2.close()




c2 = DIR_AP+'/final3/chembl/dataset/'
sdfs = glob.glob(c2+'*extend_1.sdf')

# sdfs = glob.glob(c2+'sdfs/*2[0-9].sdf') +glob.glob(c2+'sdfs/*3[0-9].sdf')


split_number= 5000 # (number of molecules in each file--delete this line in bracket)
number_of_sdfs = split_number
i=0
j=0
f2=open(c2+'sdfs/extend_1_'+str(j)+'.sdf','w')
for sdf in sdfs:
	for line in open(sdf):
		_ = f2.write(line)
		if line[:4] == "$$$$":
			i+=1
		if i > number_of_sdfs:
			number_of_sdfs += split_number 
			f2.close()
			j+=1
			f2=open(c2+'sdfs/extend_1_'+str(j)+'.sdf','w')
			print(j)

f2.close()


c2 = DIR_AP+'/final3/dude/dataset/'
split_number= 5000 # (number of molecules in each file--delete this line in bracket)
number_of_sdfs = split_number
i=0
j=0
f2=open(c2+'decoy_sdfs/decoys_'+str(j)+'.sdf','w')
for sdf in sdfs:
	for line in open(sdf):
		_ = f2.write(line)
		if line[:4] == "$$$$":
			i+=1
		if i > number_of_sdfs:
			number_of_sdfs += split_number 
			f2.close()
			j+=1
			f2=open(c2+'decoy_sdfs/decoys_'+str(j)+'.sdf','w')
			print(j)

f2.close()



c2 = DIR_AP+'/final/dti_sdf'
sdfs = glob.glob(c2+'/*')


split_number= 5000 # (number of molecules in each file--delete this line in bracket)
number_of_sdfs = split_number
i=0
j=0
f2=open('/spstorage/USERS/gina/Project/AP/final/dataset/dti_cid'+str(j)+'.sdf','w')
for sdf in sdfs:
	for line in open(sdf):
		_ = f2.write(line)
		if line[:4] == "$$$$":
			i+=1
		if i > number_of_sdfs:
			number_of_sdfs += split_number 
			f2.close()
			j+=1
			f2=open('/spstorage/USERS/gina/Project/AP/final/dataset/dti_cid'+str(j)+'.sdf','w')
			print(j)

f2.close()

c2 = DIR_AP+'/final2/dataset/sdfs/splits'
sdfs = glob.glob(c2+'/add*')

c2 = DIR_AP+ '/final3/pdbbind/dataset/sdfs/'
# c2 = c2+'sdfs/aa_cids_10_'
sdfs = glob.glob(c2+'*')



def sdf2scf(sdf):
	mm = Chem.SDMolSupplier(sdf, removeHs=False) 
	mn = [m.GetPropsAsDict()['PUBCHEM_COMPOUND_CID'] for m in mm if m is not None]
	scf = [Chem.MolToSmiles(GetScaffoldForMol(m)) for m in mm]
	scfs = pd.DataFrame(scf)
	scfs.index = mn
	return scfs

vdf = mclapply(sdfs, sdf2scf, ncpu=48)


app = pd.concat(vdf)
app = app[~app.index.duplicated()]
app.columns = ['scf']

app.to_pickle('{}/final3/chembl/dataset/mtxs/dti_scfs.pkl'.format(DIR_AP))
app.to_pickle('{}/final3/chembl/dataset/mtxs/scfs.pkl'.format(DIR_AP))
scfs = pd.read_pickle('{}/final3/chembl/dataset/mtxs/scfs.pkl'.format(DIR_AP))
dscfs = pd.read_pickle('{}/final3/chembl/dataset/mtxs/scfs.pkl'.format(DIR_AP))

radius = 2

def sdf2fp(pps):
	fls = glob.glob('{}/final3/dude/vina_docked/{}/*'.format(DIR_AP,pps))
	vtss = pd.DataFrame()
	for adidx in ['actives_final','decoys_final_0','decoys_final_10']:
		if '{}/final3/dude/vina_docked/{}/{}_docked_vina.sdf'.format(DIR_AP, pps, adidx) not in fls:
			continue
		mm = Chem.SDMolSupplier('{}/final3/dude/vina_docked/{}/{}_docked_vina.sdf'.format(DIR_AP, pps, adidx), removeHs=True)
		mn = [m.GetProp('_Name') for m in mm if m is not None]
		fp = [list(Chem.rdMolDescriptors.GetMorganFingerprintAsBitVect(m, radius)) for m in mm if m is not None]
		fps = pd.DataFrame(fp)
		fps.index = mn
		fps['target'] = pps
		if 'actives' in adidx:
			fps['label'] =1
		else:
			fps['label'] = 0				
		vtss = pd.concat([vtss, fps])
	return vtss	


vdf = mclapply(sdfs, sdf2fp, ncpu=33)

app = pd.concat(vdf)

app[app.sum()==0]

app = pd.concat(vdf)
app = app[~app.index.duplicated()]

app1 = pd.read_pickle('{}/final3/pdbbind/dataset/mtxs/fps.pkl'.format(DIR_AP))

app.to_pickle('{}/final3/chembl/dataset/mtxs/extend_fps.pkl'.format(DIR_AP))
			
def sdf2fp(sdfs):
	mm = Chem.SDMolSupplier(sdfs, removeHs=False) 
	mn = [m.GetPropsAsDict()['PUBCHEM_COMPOUND_CID'] for m in mm if m is not None]
	fp = [list(Chem.rdMolDescriptors.GetMorganFingerprintAsBitVect(m, radius)) for m in mm if m is not None]
	fps = pd.DataFrame(fp)
	fps.index = mn
	return fps

def sdf2fp(sdfs):
	mm = Chem.SDMolSupplier(sdfs, removeHs=False) 
	mn = [m.GetPropsAsDict()['source_label'] for m in mm if m is not None]
	fp = [list(Chem.rdMolDescriptors.GetMorganFingerprintAsBitVect(m, radius)) for m in mm if m is not None]
	fps = pd.DataFrame(fp)
	fps.index = mn
	return fps


c2 =DIR_AP+'/DB/refined-set'
pds1=[pps.split('/')[-1] for pps in glob.glob(c2+'/*')]
c2 = DIR_AP+'/DB/v2020-other-PL'
pds = pds1 + [pps.split('/')[-1] for pps in glob.glob(c2+'/*')]


			
def sdf2fp(pps):
	try:
		if pps in pds1:didx = 'refined-set'
		else:didx = 'v2020-other-PL'
		m = Chem.SDMolSupplier('{}/DB/{}/{}/{}_ligand.sdf'.format(DIR_AP,didx,pps,pps), removeHs=True)[0]
		if m is None:
			m = Chem.rdmolfiles.MolFromMol2File('{}/DB/{}/{}/{}_ligand.mol2'.format(DIR_AP,didx,pps,pps), removeHs=True)
		fps = list(Chem.rdMolDescriptors.GetMorganFingerprintAsBitVect(m, radius))	
	except:
		fps = [0]*2048
	return tuple(fps)


vdf = mclapply(sdfs, sdf2fp, ncpu=20)


app = pd.DataFrame(vdf)
app.index = pds

app.sum()

app = pd.concat(vdf)
app = app[~app.index.duplicated()]

app1 = pd.read_pickle('{}/final3/pdbbind/dataset/mtxs/fps.pkl'.format(DIR_AP))
app = pd.concat([app1,app])
app.to_pickle('{}/final3/pdbbind/dataset/mtxs/fps.pkl'.format(DIR_AP))


app.to_pickle('{}/final3/chembl/dataset/mtxs/dti_fps.pkl'.format(DIR_AP))

app = app[app.T.sum()!=0]
app.to_pickle('{}/final3/pdbbind/dataset/mtxs/pdb_fps.pkl'.format(DIR_AP))


vdf = mclapply(sdfs, sdf2fp, ncpu=len(sdfs))
app = pd.concat(vdf)



app.to_pickle('{}/final/dataset/add3_fps.pkl'.format(DIR_AP))



### APM generate --1 

def filter_feats_both(feat):
	dns = [f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')]
	feat = [f for f in feat if not (f.GetFamily()=='Acceptor')&(f.GetAtomIds() in dns)]
	return feat

def filter_feats_hyp(feat):
	arh = [f.GetAtomIds() for f in feat if (f.GetFamily()=='Aromatic')]
	feat = [f for f in feat if not (f.GetFamily()=='Hydrophobe')&(f.GetAtomIds() in arh)]
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

def bulk_cids_ratio(sdfs):
	mm = Chem.SDMolSupplier(sdfs, removeHs=True) 
	mn = [m.GetPropsAsDict()['PUBCHEM_COMPOUND_CID'] for m in mm]
	das = [f.GetAtomIds() for f in feats if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')]
	das = unique([x for x in das if das.count(x) > 1])
	feats = [filter_feats_hyp_acc(m) for m in mm]
	vt = [pairLr_dist(f, rn6, cn_r, inv, das) for f in feats]
	vdf = pd.DataFrame(vt)
	vdf.index = mn
	return vdf


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


def pairLs_dist(fts, row, col, bn, das, fn):
	ft = itertools.combinations(fts, 2)
	lt = [0]*len(row)*len(col)
	for f in ft:
		dst=cal_dst(f)	
		if dst!=0:
			cv= int(dst//bn)
			if (cv>0) & (cv<len(col)-1):n1=col[cv-1];n2=col[cv+1]
			else:
				if cv<0:n1=0;n2=col[0]
				else:n1=col[-1];n2=20			
			b=(stats.norm.cdf(n2,dst,0.5)-stats.norm.cdf(n1,dst,0.5))
			b1=stats.norm.cdf(n1,dst,0.5)
			b2=1-(stats.norm.cdf(n2,dst,0.5))
			aType = fn.index(f[0].GetFamily())
			bType = fn.index(f[1].GetFamily())
			if f[0].GetAtomIds() in das:aType=9
			if f[1].GetAtomIds() in das:bType=9
			Tpair = str(min(aType,bType))+str(max(aType,bType))
			if (dst!=0) & (Tpair in row) & (dst<20):
				tidx = row.index(Tpair)
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

'CIT, 5NI, 5EN, OTA, AXN, ADP, 6ZX, GC2, 164, 8FM, G04, ANV, 0HO, 5O5, J27, L14, 5RO, F6W, KMN, ET, 27O, OMD, 2J5, N5N, F1Q, VAB, 12G, GA3, DCL, G88, 72D, 2-mer, I2E, V49, AEZ, 9PB, 5EF, P34, 1XT, AVJ, PRO, WHA, MET, A8Q, VD9, D07, LNR, A2P, LZ1, 0QY, 3EF, B55, F52, SAM, HRG, BTN, AVX, TDR, AKG, 6RB, 4H2, 23L, MFU, IO2, UL7, 1B3, 5ZE, UX9, V51, G07, 5ON, ARG, 0HN, 3NY, 98, 6FG, A8K, 3P3, K19, AV9, 9LN, TR7, 7FF, BDP, 36Z, RZH, G08, BS5, GRI, 0FT, NXG, V26, 0J9, 7QY, 2C9, A98, SO8, JG1, FYM, 8HD, QIC, 1WN, OAN, 3C7, TVZ, 96Y, 1JR, A7K, DE5, K82, G64, WTZ, 777, G39, HSX, 0JV, 4ME, FKE, A65, 0G0, 520, GDP, TTN, LOC, 39R, 1BW, ABQ, KB1, 6XC, F0Y, I29, SRO, M86, 9UN, SW1, 5-mer, 4L7, M73, 3UD, KTS, 7-mer, HCE, HWY, 017, 6KK, CHT, 9GQ, D9B, K0I, 5FX, 88W, 27N, GTP, 8JQ, SSI, P51, 44O, 5XS, 1K3, F1E, 149, 5CX, R47, 4CC, SAH, IFM, 2VD, BEA, C90, 4UY, 2ZM, 5P8, MBN, 14F, 5O6, AB1, 72H, 88R, 33R, RGJ, 3NG, 1NP, 2L2, 6Z9, 56O, N7P, 6DP, FB2, P74, B7I, SKM, K9S, 4FH, ANP, 2VF, PL8, 1D6, 3-mer, 2KT, 3QA, 10O, 9UW, PHE, 091, A4Q, IPH, 3XN, 4ZE, 949, 5OJ, 4-mer, 79F, 4UX, ZB6, O4N, 23K, 6JY, 7MH, ROC, 1U5, G61, 5MI, 16G, 9UT, YJA, 60P, BEW, 833, 37N, 44N, C5P, 91B, 5WD, GSY, SFB, G8Z, E7E, 7R7, TI8, 7K0, M4Z, 2L1, FKQ, G89, PHB, PYR, 8ZE, MN0, MPV, 6MS, 9GN, 9GW, JPZ, K13, R40, Y0V, MK5, N44, E1B, 913, A8H, K1T, 44L, THM, 6O5, ZT2, 6K8, 5RL, NOK, 5H1, J4K, GR7, 7FH, 44K, 3F0, 18U, MRW, ICO, 0LX, 2ZV, SU3, GR8, TPV, ZT4, ADK, B50, S3G, U6M, DSN, QFH, FT8, ZH2, J0S, 0OC, ASP, 1DL, 478, FKK, CWJ, VHE, 24W, YKI, VXO, SER, 06P, TDI, 8V5, PC, 1U1, 58T, MYC, 0NL, ZMR, A3S, 6VD, A60, AMG, DTP, DQT, GLN, MEW, 1OS, AKH, 54F, BR0, A2G, R4B, CTT, B04, BNZ, 1TS, 15T, 6KT, G53, 8FP, 6FJ, A9K, N8P, YPW, 5MH, 5GP, IMP, MCF, U5P, 56Y, G79, ABV, DTB, ES3, 76X, JB7, A61, MOF, 94W, 6P4, 031, RPI, 2AP, RAM, ZOF, 8ZT, 5WZ, 70A, PGA, 51J, 89J, 1XM, 513, 0KJ, P84, 6MH, ZSV, HIC, I31, 4SZ, 4QJ, O82, 03V, L4T, UFV, 5B7, BBY, SDS, TSN, LOG, 0W1, K14, VJJ, BIG, AZZ, 9ZE, G10, X6P, X0P, AMP, F0W, HRD, FZM, 6FR, 147, BFM, MRZ, BP7, CYT, 23B, 594, GR5, P86, 6CE, GLU, VXQ, EDG, MGT, 5B5, 8W9, A3T, 0XR, SHH, GAL-MHD, 44Q, 5OO, HQT, SAL, UO1, 5GO, IPT, 3EB, 7GR, 1XN, 9-mer, 6-mer, BZI, EO1, AU8, 1DK, 8V8, 7ZE, L36, E88'
use_mmff = [int(float(m.split(' ')[1])//gap + 10) for m in mmffp[idx][1::] if bool(mps[idx][int(m.split(' ')[0])-1] not in use_pos)]# if int(m.split(' ')[0])<mhac[idx]]
	

def MMFFLr_dist_gap_feats(feats, mmffp, mps, idx, col, gap, inv, das):
	row = [str(x)+str(y) for x,y in combinations_with_replacement([0,1,2,3,5,6,7,8,9,10,11,12,13],2)]
	use_pos = np.array([ft.GetPos() for ft in feats[idx]])
	use_fam = [fn.index(ft.GetFamily()) for ft in feats[idx]]
	mmffps = [m for m in mmffp[idx][1::] if len(mps[idx])>=int(m.split(' ')[0])]
	use_idx = np.array([mps[idx][int(m.split(' ')[0])-1] for m in mmffps])# if bool(mps[idx][int(m.split(' ')[0])-1] in use_pos)])# not in f_idx]# if int(m.split(' ')[0])<mhac[idx]]
	use_mmff = [int(float(m.split(' ')[1])//gap + 11) for m in mmffp[idx][1::]]# if bool(mps[idx][int(m.split(' ')[0])-1] in use_pos)]# if int(m.split(' ')[0])<mhac[idx]]
	if len(use_idx) != 0:poses = np.concatenate((use_pos,use_idx))
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


gap=0.6
def feature_PP(feats, mmffp, mps, idx, col, gap, inv, das, fn):
	row = [str(x)+str(y) for x,y in combinations_with_replacement([0,1,2,3,5,6,9,10,11,12,13,14],2)]
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

sym2feat = {'C':10,'O':11,'S':12}
def feature_AA(feats, atms, mps, idx, col, inv, das, fn):
	# basic charac. extra atom to atom charac. 10 C / 11 O / 12 S / 13 etc
	row = [str(x)+str(y) for x,y in combinations_with_replacement([0,1,2,3,4,6,7,9,10,11,12,13],2)]
	ft_idx = unique(list(itertools.chain(*[f.GetAtomIds() for f in feats[idx]])))
	uu_idx = [i for i in range(len(mps[idx])) if i not in ft_idx]
	use_pos = np.array([ft.GetPos() for ft in feats[idx]])
	use_fam = [fn.index(ft.GetFamily()) if not ft.GetAtomIds() in das[idx] else 9 for ft in feats[idx]]
	if (len(uu_idx) != 0) & (len(use_pos) !=0):
		use_idx = np.array(itemgetter(*uu_idx)(mps[idx]))
		if len(uu_idx)>1:
			use_atms = [sym2feat[at.GetSymbol()] if at.GetSymbol() in ['C','O','S'] else 13 for at in itemgetter(*uu_idx)(atms[idx])]
			poses = np.concatenate((use_pos,use_idx))
		else:
			atsym = itemgetter(*uu_idx)(atms[idx]).GetSymbol()
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


sym2feat = {'C':10,'O':11,'S':12}
def feature_AP(feats, atms, mmffp, gap, mps, idx, col, inv, das, fn):
	# basic charac. extra atom to atom charac. 10 C / 11 O / 12 S / 13 etc
	row = [str(x)+str(y) for x,y in combinations_with_replacement([0,1,2,3,4,6,7,9,10,11,12,13,14,15,16,17,18],2)]
	ft_idx = unique(list(itertools.chain(*[f.GetAtomIds() for f in feats[idx]])))
	uu_idx = [i for i in range(len(mps[idx])) if i not in ft_idx]
	use_pos = np.array([ft.GetPos() for ft in feats[idx]])
	use_fam = [fn.index(ft.GetFamily()) if not ft.GetAtomIds() in das[idx] else 9 for ft in feats[idx]]
	if (len(uu_idx) != 0) & (len(use_pos) !=0):
		use_idx = np.array(itemgetter(*uu_idx)(mps[idx]))
		if len(uu_idx)>1:
			use_atms = [sym2feat[at.GetSymbol()] if at.GetSymbol() in ['C','O','S'] else 13 for at in itemgetter(*uu_idx)(atms[idx])]
			poses = np.concatenate((use_pos,use_idx))
		else:
			atsym = itemgetter(*uu_idx)(atms[idx]).GetSymbol()
			poses = np.concatenate((use_pos,np.array([use_idx])))
			if atsym in ['C','O','S']:
				use_atms = [sym2feat[atsym]]
			else:
				use_atms = [13]
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


sfds= pdbs[0]
def bulk_cids_m_f(sdfs):
	mm = Chem.SDMolSupplier(sdfs, removeHs=True) 
	mpd = [m.GetPropsAsDict() for m in mm]
	mmffp = [m['PUBCHEM_MMFF94_PARTIAL_CHARGES'].split('\n') if 'PUBCHEM_MMFF94_PARTIAL_CHARGES' in m.keys() else '' for m in mpd]
	mn = [m['PUBCHEM_COMPOUND_CID'] for m in mpd]
	mps = [m.GetConformer().GetPositions() for m in mm]
	feats = [factory.GetFeaturesForMol(m) for m in mm]
	vt = [MMFFLr_dist_gap_feats(feats, mmffp, mps, idx, cn_r, gap,inv) for idx in range(len(mpd))]
	vdf = pd.DataFrame(vt)
	vdf.index = mn
	return vdf


def bulk_cids(sdf):
	mm = Chem.SDMolSupplier(sdf, removeHs=True) 
	mpd = [m.GetPropsAsDict() for m in mm if m is not None]
	mmffp = [m['PUBCHEM_MMFF94_PARTIAL_CHARGES'] if 'PUBCHEM_MMFF94_PARTIAL_CHARGES' in m.keys() else '' for m in mpd]
	mmffp = [m.split('\n') if (m!=0) & (m!='') else '' for m in mmffp]
	mn = [m['PUBCHEM_COMPOUND_CID'] for m in mpd]
	mps = [m.GetConformer().GetPositions() for m in mm if m is not None]
	feats = [factory.GetFeaturesForMol(m) for m in mm if m is not None]
	das = [[]]*len(feats)
	if filter_type==0:
		feats = [filter_feats_hyp(feat) for feat in feats]
	else:
		feats = [filter_feats_renew(feat) for feat in feats]
	if facts[2]==1:
		das = [[f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')] for feat in feats]
		das = [unique([x for x in da if da.count(x) > 1]) for da in das]
		feats = [filter_feats_both(feat) for feat in feats]
	vt = [feature_PP(feats, mmffp, mps, idx, cnrs, gap, inv, das, fn) for idx in range(len(mpd))]
	vt = [pairLr_dist(f, rns, cnrs, inv, das[i], fn) for i, f in enumerate(feats)]
	vdf = pd.DataFrame(vt)
	vdf.index = mn
	# vdf['target'] = sdf.split('_')[-1].split('.')[0]
	return vdf


# generate by oplist options
a = [[8,10,12],[0,1],[0,1,2]]
oplist = list(itertools.product(*a))
					
# factory = ChemicalFeatures.BuildFeatureFactory('/spstorage/USERS/gina/Project/AP/etc/BaseFeatures.fdef')
# fn1= factory.GetFeatureFamilies()
rn6=[str(x)+str(y) for x,y in combinations_with_replacement([0,1,2,3,5,6],2)]
rn7 = [str(x)+str(y) for x,y in combinations_with_replacement([0,1,2,3,5,6,9],2)]
# factory2 = ChemicalFeatures.BuildFeatureFactory('/spstorage/USERS/gina/Project/AP/etc/BaseFeatures_h.fdef')
# fn2= factory2.GetFeatureFamilies()
rnh = [str(x)+str(y) for x,y in combinations_with_replacement([0,1,2,3,4,6,7],2)]
rn_list = [rn6,rn7,rnh]
# factory_list = [factory, factory, factory2]
# fn_list = [fn1,fn1,fn2]
f_list=['.fdef','.fdef','_h.fdef']
# [coridx[i][cidx] for i, cidx in enumerate([corp[0][corp[1].index(max(corp[1]))] for corp in coroplist])]
# [min(corp[1]) for corp in coroplist]


# c2 = DIR_AP+'/final/dataset/'
# sdfs = glob.glob(c2+'*cid*.sdf')



for facts in [(10,0,1),(10,1,1),(12,0,1),(12,1,1)]:
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
		mmffp = [m['PUBCHEM_MMFF94_PARTIAL_CHARGES'] if 'PUBCHEM_MMFF94_PARTIAL_CHARGES' in m.keys() else '' for m in mpd]
		mmffp = [m.split('\n') if (m!=0) & (m!='') else '' for m in mmffp]
		mn = [m['PUBCHEM_COMPOUND_CID'] for m in mpd]
		mps = [m.GetConformer().GetPositions() for m in mm if m is not None]
		feats = [factory.GetFeaturesForMol(m) for m in mm if m is not None]
		das = [[]]*len(feats)
		if filter_type==0:
			feats = [filter_feats_hyp(feat) for feat in feats]
		else:
			feats = [filter_feats_renew(feat) for feat in feats]
		if facts[2]==1:
			das = [[f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')] for feat in feats]
			das = [unique([x for x in da if da.count(x) > 1]) for da in das]
			feats = [filter_feats_both(feat) for feat in feats]
		vt = [feature_PP(feats, mmffp, mps, idx, cnrs, gap, inv, das, fn) for idx in range(len(mpd))]
		# vt = [pairLr_dist(f, rns, cnrs, inv, das[i], fn) for i, f in enumerate(feats)]
		vdf = pd.DataFrame(vt)
		vdf.index = mn
		# vdf['target'] = sdf.split('_')[-1].split('.')[0]
		return vdf
	aps = mclapply(sdfs, bulk_cids, ncpu=len(sdfs))
	app = pd.concat(aps)
	# app = pd.concat([prep,app])
	app.to_pickle('{}/final2/dataset/mtxs/pp_apms_{}_{}_{}.pkl'.format(DIR_AP, facts[0], facts[1], facts[2]))
	del app



for facts in [(10,1,2),(12,1,2)]:
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
		vt = [feature_AA(feats, atms, mps, idx, cnrs, inv, das, fn) for idx in range(len(mpd))]
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


sdfs =  glob.glob(c2+'splits/add2*.sdf')

for facts in [(10,1,2),(12,1,2)]:
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
		vt = [feature_AA(feats, atms, mps, idx, cnrs, inv, das, fn) for idx in range(len(mpd))]
		# vt = [pairLr_dist(f, rns, cnrs, inv, das[i], fn) for i, f in enumerate(feats)]
		vdf = pd.DataFrame(vt)
		vdf.index = mn
		# vdf['target'] = sdf.split('_')[-1].split('.')[0]
		return vdf
	aps = mclapply(sdfs, bulk_cids, ncpu=len(sdfs))
	app = pd.concat(aps)
	# app = pd.concat([prep,app])
	app.to_pickle('{}/final2/dataset/mtxs/add2_aa_apms_{}_{}_{}.pkl'.format(DIR_AP, facts[0], facts[1], facts[2]))
	del app


types = ['fps',(10,1,2,4), (10,1,2,5),(12,1,2,4)]

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
		mmffp = [m['PUBCHEM_MMFF94_PARTIAL_CHARGES'] if 'PUBCHEM_MMFF94_PARTIAL_CHARGES' in m.keys() else '' for m in mpd]
		mmffp = [m.split('\n') if (m!=0) & (m!='') else '' for m in mmffp]
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
		vt = [feature_AP(feats, atms, mmffp, gap,mps, idx, cnrs, inv, das, fn) for idx in range(len(mpd))]
		# vt = [pairLr_dist(f, rns, cnrs, inv, das[i], fn) for i, f in enumerate(feats)]
		vdf = pd.DataFrame(vt)
		vdf.index = mn
		# vdf['target'] = sdf.split('_')[-1].split('.')[0]
		return vdf
	aps = mclapply(sdfs, bulk_cids, ncpu=50)
	app = pd.concat(aps)
	# app = pd.concat([prep,app])
	app.to_pickle('{}/final2/dataset/mtxs/ap_apms_{}_{}_{}.pkl'.format(DIR_AP, facts[0], facts[1], facts[2]))
	del app

sdfs =  glob.glob(c2+'splits/add3*.sdf')
for facts in [(10,1,2),(12,1,2)]:
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
		mmffp = [m['PUBCHEM_MMFF94_PARTIAL_CHARGES'] if 'PUBCHEM_MMFF94_PARTIAL_CHARGES' in m.keys() else '' for m in mpd]
		mmffp = [m.split('\n') if (m!=0) & (m!='') else '' for m in mmffp]
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
		vt = [feature_AP(feats, atms, mmffp, gap,mps, idx, cnrs, inv, das, fn) for idx in range(len(mpd))]
		# vt = [pairLr_dist(f, rns, cnrs, inv, das[i], fn) for i, f in enumerate(feats)]
		vdf = pd.DataFrame(vt)
		vdf.index = mn
		# vdf['target'] = sdf.split('_')[-1].split('.')[0]
		return vdf
	aps = mclapply(sdfs, bulk_cids, ncpu=len(sdfs))
	app = pd.concat(aps)
	# app = pd.concat([prep,app])
	app.to_pickle('{}/final2/dataset/mtxs/add3_ap_apms_{}_{}_{}.pkl'.format(DIR_AP, facts[0], facts[1], facts[2]))
	del app



# for facts in oplist[3:]:

for facts in [(10,1,0),(10,1,1)]:
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
		mn = [m.GetPropsAsDict()['PUBCHEM_COMPOUND_CID'] for m in mm if m is not None]
		feats = [factory.GetFeaturesForMol(m) for m in mm if m is not None]
		das = [[]]*len(feats)
		if filter_type==0:
			feats = [filter_feats_hyp(feat) for feat in feats]
		else:
			feats = [filter_feats_renew(feat) for feat in feats]
		if facts[2]==1:
			das = [[f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')] for feat in feats]
			das = [unique([x for x in da if da.count(x) > 1]) for da in das]
			feats = [filter_feats_both(feat) for feat in feats]
		vt = [pairLr_dist(f, rns, cnrs, inv, das[i], fn) for i, f in enumerate(feats)]
		vdf = pd.DataFrame(vt)
		vdf.index = mn
		# vdf['target'] = sdf.split('_')[-1].split('.')[0]
		return vdf
	aps = mclapply(sdfs, bulk_cids, ncpu=50)
	app = pd.concat(aps)
	# app = pd.concat([prep,app])
	app.to_pickle('{}/final2/dataset/mtxs/apms_{}_{}_{}.pkl'.format(DIR_AP, facts[0], facts[1], facts[2]))
	del app

from rdkit.Chem.Scaffolds.MurckoScaffold import GetScaffoldForMol

c2 = DIR_AP+'/final/dataset/'
sdfs = glob.glob(c2+'*cid*.sdf')



def Get_scf(sdf):
	mm = Chem.SDMolSupplier(sdf, removeHs=False) 
	scf = [Chem.MolToSmiles(GetScaffoldForMol(m)) for m in mm if m is not None]
	mn = [m.GetPropsAsDict()['PUBCHEM_COMPOUND_CID'] for m in mm if m is not None]
	scfs = pd.DataFrame(scf)
	scfs.index= mn
	return scfs

aps = mclapply(sdfs, Get_scf, ncpu=50)
app = pd.concat(aps)
# app = pd.concat([prep,app])
app.to_pickle('{}/final2/dataset/scf.pkl'.format(DIR_AP))


###### 

scfs = pd.read_pickle('{}/final/dataset/scf.pkl'.format(DIR_AP))

with open(DIR_AP+'/final/dataset/chembl_split.pkl', 'rb') as f:foldss=pickle.load(f)



## 정리


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

## 기본 

import glob
c2 = DIR_AP+'/DB/DUDE/'
sdfs = glob.glob(c2+'/*.sdf')
sdfs = sorted(sdfs)
subdirs = [x for x in os.walk(c2)][2][1]


c2 = DIR_AP+'/final3/chembl/dataset/'
sdfs = glob.glob(c2+'sdfs/extend_1_*.sdf')
sdfs = glob.glob(DIR_AP+'/DB/MUV/muv_[0-9].sdf')+glob.glob(DIR_AP+'/DB/MUV/muv_[0-9][0-9].sdf')


sdfs = ['/spstorage/USERS/gina/Project/AP/final3/pdbbind/dataset/sdfs/cids_10_{}.sdf'.format(str(nn)) for nn in range(32,38)]




dstbin = 10
f_list=['','','','_all'] +['_h']*4
st = 20/(inv**dstbin) # 16 bin
cnrs = [st*(inv**g) for g in range(dstbin)]
fnss = [[0,1,2,3,5,6],[0,1,2,3,5,6],[0,1,2,3,5,6,9],[0,1,2,3,4,5,7,8,10],[0,1,2,3,4,6,7,9],[0,1,2,3,4,6,7,9,10,11,12,13,14],[0,1,2,3,4,6,7,9,10,11,12,13,14],[0,1,2,3,4,6,7,9,10,11,12,13,14,15,16,17,18]]
for fact in range(8):
	factory = ChemicalFeatures.BuildFeatureFactory('{}/etc/BaseFeatures{}.fdef'.format(DIR_AP,f_list[fact]))
	fn = factory.GetFeatureFamilies()
	fns = fnss[fact]
	rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
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
					mps = [m.GetConformer().GetPositions() for m in mm if m is not None]		
					if fact > 5:
						atms = [m.GetAtoms() for m in mm if m is not None]
						if fact==7:
							mpd = [m.GetPropsAsDict() for m in mm if m is not None]
							mmffp = [m['PUBCHEM_MMFF94_PARTIAL_CHARGES'] if 'PUBCHEM_MMFF94_PARTIAL_CHARGES' in m.keys() else '' for m in mpd]
							mmffp = [m.split('\n') if (m!=0) & (m!='') else '' for m in mmffp]
							vt = [feature_AP(rns, feats, atms, mmffp, gap,mps, idx, cnrs, inv, das, fn) for idx in range(len(mpd))]			
						else:
							vt = [feature_AA(rns, feats, atms, mps, idx, cnrs, inv, das, fn) for idx in range(len(mps))]
					else:
						mpd = [m.GetPropsAsDict() for m in mm if m is not None]
						mmffp = [m['PUBCHEM_MMFF94_PARTIAL_CHARGES'] if 'PUBCHEM_MMFF94_PARTIAL_CHARGES' in m.keys() else '' for m in mpd]
						mmffp = [m.split('\n') if (m!=0) & (m!='') else '' for m in mmffp]		
						vt = [feature_PP(rns,feats, mmffp, mps, idx, cnrs, gap, inv, das, fn) for idx in range(len(mpd))]
				else:
					vt = [pairLr_dist(f, rns, cnrs, inv, das[i], fn) for i, f in enumerate(feats)]
			else:
				vt = [pairLr_dist(f, rns, cnrs, inv, das[i], fn) for i, f in enumerate(feats)]
		else:
			vt = [pairLr_dist(f, rns, cnrs, inv, das[i], fn) for i, f in enumerate(feats)]
		vdf = pd.DataFrame(vt)
		vdf.index = mn
		return vdf
	aps = mclapply([[sdf, fact, rns, cnrs, fn] for sdf in sdfs], bulk_cids, ncpu=len(sdfs))
	app = pd.concat(aps)
	# app2 = pd.read_pickle('{}/final3/pdbbind/dataset/mtxs/apms_{}.pkl'.format(DIR_AP, fact))
	# app = pd.concat([app2,app])
	app.to_pickle('{}/final3/chembl/dataset/mtxs/dti_apms_{}.pkl'.format(DIR_AP, fact))
	del app




dstbin = 10
f_list=['','','','_all'] +['_h']*8
st = 20/(inv**dstbin) # 16 bin
cnrs = [st*(inv**g) for g in range(dstbin)]
fnss = [[0,1,2,3,5,6],[0,1,2,3,5,6],[0,1,2,3,5,6,9],[0,1,2,3,4,5,7,8,10],[0,1,2,3,4,6,7,9],[0,1,2,3,4,6,7,9,10,11,12,13,14],[0,1,2,3,4,6,7,9,10,11,12,13,14],[0,1,2,3,4,6,7,9,10,11,12,13,14,15,16,17,18],[0,1,2,3,4,6,7,9,10,11,12,13,14],[0,1,2,3,4,6,7,9,10,11,12,13,14],[0,1,2,3,4,6,7,9,10,11,12,13,14,15,16,17,18],[0,1,2,3,4,6,7,9,10,11]]
for fact in [0,6,8,9,10,11]:
	factory = ChemicalFeatures.BuildFeatureFactory('{}/etc/BaseFeatures{}.fdef'.format(DIR_AP,f_list[fact]))
	fn = factory.GetFeatureFamilies()
	fns = fnss[fact]
	rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
	def bulk_cids(ssets):
		sdf, fact, rns, cnrs, fn = ssets
		mm = Chem.SDMolSupplier(sdf, removeHs=True) 
		mn = [m.GetPropsAsDict()['PUBCHEM_COMPOUND_CID'] for m in mm if m is not None]
		feats = [factory.GetFeaturesForMol(m) for m in mm if m is not None]
		das = [[]]*len(feats)
		if fact == 0:
			vt = [pairLr_dist(f, rns, cnrs, inv, das[i], fn) for i, f in enumerate(feats)]
		else:
			if fact ==6:
				feats = [filter_feats_renew(feat) for feat in feats]
			das = [[f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')] for feat in feats]
			das = [unique([x for x in da if da.count(x) > 1]) for da in das]
			mps = [m.GetConformer().GetPositions() for m in mm if m is not None]		
			if fact != 8:
				atms = [m.GetAtoms() for m in mm if m is not None]
				if fact == 9:
					vt = [feature_AA(rns, feats, atms, mps, idx, cnrs, inv, das, fn) for idx in range(len(mps))]
				elif fact==10:
					mpd = [m.GetPropsAsDict() for m in mm if m is not None]
					mmffp = [m['PUBCHEM_MMFF94_PARTIAL_CHARGES'] if 'PUBCHEM_MMFF94_PARTIAL_CHARGES' in m.keys() else '' for m in mpd]
					mmffp = [m.split('\n') if (m!=0) & (m!='') else '' for m in mmffp]
					vt = [feature_AP(rns, feats, atms, mmffp, gap,mps, idx, cnrs, inv, das, fn) for idx in range(len(mpd))]			
				else:
					vt = [feature_AH(rns, feats, atms, mps, idx, cnrs, inv, das, fn) for idx in range(len(mps))]
			else:
				mpd = [m.GetPropsAsDict() for m in mm if m is not None]
				mmffp = [m['PUBCHEM_MMFF94_PARTIAL_CHARGES'] if 'PUBCHEM_MMFF94_PARTIAL_CHARGES' in m.keys() else '' for m in mpd]
				mmffp = [m.split('\n') if (m!=0) & (m!='') else '' for m in mmffp]		
				vt = [feature_PP(rns,feats, mmffp, mps, idx, cnrs, gap, inv, das, fn) for idx in range(len(mpd))]
		vdf = pd.DataFrame(vt)
		vdf.index = mn
		return vdf
	aps = mclapply([[sdf, fact, rns, cnrs, fn] for sdf in sdfs], bulk_cids, ncpu=33)
	app = pd.concat(aps)
	app.to_pickle('{}/final3/chembl/dataset/mtxs/extend_apms_{}.pkl'.format(DIR_AP, fact))
	del app

dstbin = 10
f_list=['','','','_all'] +['_h']*4
st = 20/(inv**dstbin) # 16 bin
cnrs = [st*(inv**g) for g in range(dstbin)]
fnss = [[0,1,2,3,5,6],[0,1,2,3,5,6],[0,1,2,3,5,6,9],[0,1,2,3,4,5,7,8,10],[0,1,2,3,4,6,7,9],[0,1,2,3,4,6,7,9,10,11,12,13,14],[0,1,2,3,4,6,7,9,10,11,12,13,14],[0,1,2,3,4,6,7,9,10,11,12,13,14,15,16,17,18]]
for fact in range(8,9):
	factory = ChemicalFeatures.BuildFeatureFactory('{}/etc/BaseFeatures{}.fdef'.format(DIR_AP,f_list[fact-3]))
	fn = factory.GetFeatureFamilies()
	fns = fnss[fact-3]
	rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
	def bulk_cids(ssets):
		sdf, fact, rns, cnrs, fn = ssets
		mm = Chem.SDMolSupplier(sdf, removeHs=True) 
		mn = [m.GetPropsAsDict()['PUBCHEM_COMPOUND_CID'] for m in mm if m is not None]
		feats = [factory.GetFeaturesForMol(m) for m in mm if m is not None]
		das = [[]]*len(feats)
		das = [[f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')] for feat in feats]
		das = [unique([x for x in da if da.count(x) > 1]) for da in das]
		mps = [m.GetConformer().GetPositions() for m in mm if m is not None]		
		if fact > 8:
			atms = [m.GetAtoms() for m in mm if m is not None]
			if fact==10:
				mpd = [m.GetPropsAsDict() for m in mm if m is not None]
				mmffp = [m['PUBCHEM_MMFF94_PARTIAL_CHARGES'] if 'PUBCHEM_MMFF94_PARTIAL_CHARGES' in m.keys() else '' for m in mpd]
				mmffp = [m.split('\n') if (m!=0) & (m!='') else '' for m in mmffp]
				vt = [feature_AP(rns, feats, atms, mmffp, gap,mps, idx, cnrs, inv, das, fn) for idx in range(len(mpd))]			
			else:
				vt = [feature_AA(rns, feats, atms, mps, idx, cnrs, inv, das, fn) for idx in range(len(mps))]
		else:
			mpd = [m.GetPropsAsDict() for m in mm if m is not None]
			mmffp = [m['PUBCHEM_MMFF94_PARTIAL_CHARGES'] if 'PUBCHEM_MMFF94_PARTIAL_CHARGES' in m.keys() else '' for m in mpd]
			mmffp = [m.split('\n') if (m!=0) & (m!='') else '' for m in mmffp]		
			vt = [feature_PP(rns,feats, mmffp, mps, idx, cnrs, gap, inv, das, fn) for idx in range(len(mpd))]
		vdf = pd.DataFrame(vt)
		vdf.index = mn
		return vdf
	aps = mclapply([[sdf, fact, rns, cnrs, fn] for sdf in sdfs], bulk_cids, ncpu=46)
	app = pd.concat(aps)
	# app2 = pd.read_pickle('{}/final3/pdbbind/dataset/mtxs/apms_{}.pkl'.format(DIR_AP, fact))
	# app = pd.concat([app2,app])
	app.to_pickle('{}/final3/chembl/dataset/mtxs/apms_{}.pkl'.format(DIR_AP, fact))
	del app



aps = mclapply([[sdf, fact, rns, cnrs, fn] for sdf in sdfs[7:]], bulk_cids, ncpu=len(sdfs[7:]))
app = pd.concat([app, pd.concat(aps)])
app.to_pickle('{}/final3/chembl/dataset/mtxs/apms_{}.pkl'.format(DIR_AP, fact))

app.to_pickle('{}/final3/pdbbind/dataset/mtxs/apms_{}.pkl'.format(DIR_AP, fact))
app.to_pickle('{}/final3/dude/dataset/mtxs/apms_{}.pkl'.format(DIR_AP, fact))

sdfs= [DIR_AP+'/DB/DUDE/all/'+sub+'/actives_final.sdf' for sub in subdirs]
sdfs= [DIR_AP+'/DB/DUDE/all/'+sub+'/decoys_final.sdf' for sub in subdirs]

dstbin = 10
f_list=['','','','_all'] +['_h']*4
st = 20/(inv**dstbin) # 16 bin
cnrs = [st*(inv**g) for g in range(dstbin)]
fnss = [[0,1,2,3,5,6],[0,1,2,3,5,6],[0,1,2,3,5,6,9],[0,1,2,3,4,5,7,8,10],[0,1,2,3,4,6,7,9],[0,1,2,3,4,6,7,9,10,11,12,13,14],[0,1,2,3,4,6,7,9,10,11,12,13,14],[0,1,2,3,4,6,7,9,10,11,12,13,14,15,16,17,18]]
for fact in [0,1,2,3,4,6]:
	factory = ChemicalFeatures.BuildFeatureFactory('{}/etc/BaseFeatures{}.fdef'.format(DIR_AP,f_list[fact]))
	fn = factory.GetFeatureFamilies()
	fns = fnss[fact]
	rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
	def bulk_cids(ssets):
		sdf, fact, rns, cnrs, fn = ssets
		mm = Chem.SDMolSupplier(sdf, removeHs=True) 
		mn = [m.GetProp('_Name') for m in mm]
		feats = [factory.GetFeaturesForMol(m) for m in mm if m is not None]
		das = [[]]*len(feats)
		if fact > 0:
			feats = [filter_feats_renew(feat) for feat in feats]
			if fact > 1:
				das = [[f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')] for feat in feats]
				das = [unique([x for x in da if da.count(x) > 1]) for da in das]
				if fact > 4:	
					mps = [m.GetConformer().GetPositions() for m in mm if m is not None]		
					if fact > 5:
						atms = [m.GetAtoms() for m in mm if m is not None]
						if fact==7:
							mpd = [m.GetPropsAsDict() for m in mm if m is not None]
							mmffp = [m['PUBCHEM_MMFF94_PARTIAL_CHARGES'] if 'PUBCHEM_MMFF94_PARTIAL_CHARGES' in m.keys() else '' for m in mpd]
							mmffp = [m.split('\n') if (m!=0) & (m!='') else '' for m in mmffp]
							vt = [feature_AP(rns, feats, atms, mmffp, gap,mps, idx, cnrs, inv, das, fn) for idx in range(len(mpd))]			
						else:
							vt = [feature_AA(rns, feats, atms, mps, idx, cnrs, inv, das, fn) for idx in range(len(mps))]
					else:
						mpd = [m.GetPropsAsDict() for m in mm if m is not None]
						mmffp = [m['PUBCHEM_MMFF94_PARTIAL_CHARGES'] if 'PUBCHEM_MMFF94_PARTIAL_CHARGES' in m.keys() else '' for m in mpd]
						mmffp = [m.split('\n') if (m!=0) & (m!='') else '' for m in mmffp]		
						vt = [feature_PP(rns,feats, mmffp, mps, idx, cnrs, gap, inv, das, fn) for idx in range(len(mpd))]
				else:
					vt = [pairLr_dist(f, rns, cnrs, inv, das[i], fn) for i, f in enumerate(feats)]
			else:
				vt = [pairLr_dist(f, rns, cnrs, inv, das[i], fn) for i, f in enumerate(feats)]
		else:
			vt = [pairLr_dist(f, rns, cnrs, inv, das[i], fn) for i, f in enumerate(feats)]
		vdf = pd.DataFrame(vt)
		vdf.index = mn
		return vdf
	aps = mclapply([[sdf, fact, rns, cnrs, fn] for sdf in sdfs], bulk_cids, ncpu=len(sdfs))
	app = pd.concat(aps)
	app.to_pickle('{}/final3/dude/dataset/mtxs/dud_apms_{}.pkl'.format(DIR_AP, fact))
	del app


dstbin = 10
f_list=['','','','_all'] +['_h']*8
st = 20/(inv**dstbin) # 16 bin
cnrs = [st*(inv**g) for g in range(dstbin)]
fnss = [[0,1,2,3,5,6],[0,1,2,3,5,6],[0,1,2,3,5,6,9],[0,1,2,3,4,5,7,8,10],[0,1,2,3,4,6,7,9],[0,1,2,3,4,6,7,9,10,11,12,13,14],[0,1,2,3,4,6,7,9,10,11,12,13,14],[0,1,2,3,4,6,7,9,10,11,12,13,14,15,16,17,18],[],[0,1,2,3,4,6,7,9,10,11,12,13,14],[],[0,1,2,3,4,6,7,9,10,11]]
for fact in [6,9,11]:
	factory = ChemicalFeatures.BuildFeatureFactory('{}/etc/BaseFeatures{}.fdef'.format(DIR_AP,f_list[fact]))
	fn = factory.GetFeatureFamilies()
	fns = fnss[fact]
	rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
	def bulk_cids(ssets):
		pps, fact, rns, cnrs, fn = ssets
		fls = glob.glob('{}/final3/dude/vina_docked/{}/*'.format(DIR_AP,pps))
		vtss = pd.DataFrame()
		for adidx in ['actives_final','decoys_final_0','decoys_final_10']:
			if '{}/final3/dude/vina_docked/{}/{}_docked_vina.sdf'.format(DIR_AP, pps, adidx) not in fls:
				continue
			mm = Chem.SDMolSupplier('{}/final3/dude/vina_docked/{}/{}_docked_vina.sdf'.format(DIR_AP, pps, adidx), removeHs=True)
			mn = [m.GetProp('_Name') for m in mm if m is not None]
			feats = [factory.GetFeaturesForMol(m) for m in mm if m is not None]
			das = [[]]*len(feats)
			if fact > 0:
				if fact ==6:
					feats = [filter_feats_renew(feat) for feat in feats]
				das = [[f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')] for feat in feats]
				das = [unique([x for x in da if da.count(x) > 1]) for da in das]
				mps = [m.GetConformer().GetPositions() for m in mm if m is not None]		
				atms = [m.GetAtoms() for m in mm if m is not None]
				if fact ==9:
					vt = [feature_AA(rns, feats, atms, mps, idx, cnrs, inv, das, fn) for idx in range(len(mps))]
				else:
					vt = [feature_AH(rns, feats, atms, mps, idx, cnrs, inv, das, fn) for idx in range(len(mps))]
			else:
				vt = [pairLr_dist(f, rns, cnrs, inv, das[i], fn) for i, f in enumerate(feats)]
			vdf = pd.DataFrame(vt)
			vdf.index = mn
			vdf['target'] = pps
			if 'actives' in adidx:
				vdf['label'] =1
			else:
				vdf['label'] = 0		
			vtss = pd.concat([vtss, vdf])	
		return vtss
	aps = mclapply([[pps, fact, rns, cnrs, fn] for pps in subdirs], bulk_cids, ncpu=51)
	app = pd.concat(aps)
	app.to_pickle('{}/final3/dude/dataset/mtxs/vina_apms_{}.pkl'.format(DIR_AP, fact))
	del app


for fact in [9]:
	factory = ChemicalFeatures.BuildFeatureFactory('{}/etc/BaseFeatures{}.fdef'.format(DIR_AP,f_list[fact-3]))
	fn = factory.GetFeatureFamilies()
	fns = fnss[fact-3]
	rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
	def bulk_cids(ssets):
		sdf, fact, rns, cnrs, fn = ssets
		mm = Chem.SDMolSupplier(sdf, removeHs=True) 
		mn = [m.GetProp('_Name') for m in mm]
		feats = [factory.GetFeaturesForMol(m) for m in mm if m is not None]
		das = [[]]*len(feats)
		das = [[f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')] for feat in feats]
		das = [unique([x for x in da if da.count(x) > 1]) for da in das]
		mps = [m.GetConformer().GetPositions() for m in mm if m is not None]		
		atms = [m.GetAtoms() for m in mm if m is not None]
		vt = [feature_AA(rns, feats, atms, mps, idx, cnrs, inv, das, fn) for idx in range(len(mps))]
		vdf = pd.DataFrame(vt)
		vdf.index = mn
		return vdf
	aps = mclapply([[sdf, fact, rns, cnrs, fn] for sdf in sdfs], bulk_cids, ncpu=51)
	app = pd.concat(aps)
	app.to_pickle('{}/final3/dude/dataset/mtxs/dud_apms_{}.pkl'.format(DIR_AP, fact))
	del app

c2 = glob.glob(DIR_AP+'/final3/dude/dataset/decoy_sdfs/*')


dstbin = 10
f_list=['','','','_all'] +['_h']*4
st = 20/(inv**dstbin) # 16 bin
cnrs = [st*(inv**g) for g in range(dstbin)]
fnss = [[0,1,2,3,5,6],[0,1,2,3,5,6],[0,1,2,3,5,6,9],[0,1,2,3,4,5,7,8,10],[0,1,2,3,4,6,7,9],[0,1,2,3,4,6,7,9,10,11,12,13,14],[0,1,2,3,4,6,7,9,10,11,12,13,14],[0,1,2,3,4,6,7,9,10,11,12,13,14,15,16,17,18]]
for fact in [0,1,2,3,4,6]:
	factory = ChemicalFeatures.BuildFeatureFactory('{}/etc/BaseFeatures{}.fdef'.format(DIR_AP,f_list[fact]))
	fn = factory.GetFeatureFamilies()
	fns = fnss[fact]
	rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
	def bulk_cids(ssets):
		sdf, fact, rns, cnrs, fn = ssets
		mm = Chem.SDMolSupplier(sdf, removeHs=True) 
		mm = [m for m in mm if m is not None]
		if len(mm) >5000:
			random.seed(a=7)
			mm = random.sample(mm,5000)
		mn = [m.GetProp('_Name') for m in mm if m is not None]
		feats = [factory.GetFeaturesForMol(m) for m in mm if m is not None]
		das = [[]]*len(feats)
		if fact > 0:
			feats = [filter_feats_renew(feat) for feat in feats]
			if fact > 1:
				das = [[f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')] for feat in feats]
				das = [unique([x for x in da if da.count(x) > 1]) for da in das]
				if fact > 4:	
					mps = [m.GetConformer().GetPositions() for m in mm if m is not None]		
					if fact > 5:
						atms = [m.GetAtoms() for m in mm if m is not None]
						if fact==7:
							mpd = [m.GetPropsAsDict() for m in mm if m is not None]
							mmffp = [m['PUBCHEM_MMFF94_PARTIAL_CHARGES'] if 'PUBCHEM_MMFF94_PARTIAL_CHARGES' in m.keys() else '' for m in mpd]
							mmffp = [m.split('\n') if (m!=0) & (m!='') else '' for m in mmffp]
							vt = [feature_AP(rns, feats, atms, mmffp, gap,mps, idx, cnrs, inv, das, fn) for idx in range(len(mpd))]			
						else:
							vt = [feature_AA(rns, feats, atms, mps, idx, cnrs, inv, das, fn) for idx in range(len(mps))]
					else:
						mpd = [m.GetPropsAsDict() for m in mm if m is not None]
						mmffp = [m['PUBCHEM_MMFF94_PARTIAL_CHARGES'] if 'PUBCHEM_MMFF94_PARTIAL_CHARGES' in m.keys() else '' for m in mpd]
						mmffp = [m.split('\n') if (m!=0) & (m!='') else '' for m in mmffp]		
						vt = [feature_PP(rns,feats, mmffp, mps, idx, cnrs, gap, inv, das, fn) for idx in range(len(mpd))]
				else:
					vt = [pairLr_dist(f, rns, cnrs, inv, das[i], fn) for i, f in enumerate(feats)]
			else:
				vt = [pairLr_dist(f, rns, cnrs, inv, das[i], fn) for i, f in enumerate(feats)]
		else:
			vt = [pairLr_dist(f, rns, cnrs, inv, das[i], fn) for i, f in enumerate(feats)]
		vdf = pd.DataFrame(vt)
		vdf.index = mn
		return vdf
	aps = mclapply([[sdf, fact, rns, cnrs, fn] for sdf in sdfs], bulk_cids, ncpu=len(sdfs))
	app = pd.concat(aps)
	app.to_pickle('{}/final3/dude/dataset/mtxs/dud_decoy_apms_{}.pkl'.format(DIR_AP, fact))
	del app


for fact in [9]:
	factory = ChemicalFeatures.BuildFeatureFactory('{}/etc/BaseFeatures{}.fdef'.format(DIR_AP,f_list[fact-3]))
	fn = factory.GetFeatureFamilies()
	fns = fnss[fact-3]
	rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
	def bulk_cids(ssets):
		sdf, fact, rns, cnrs, fn = ssets
		mm = Chem.SDMolSupplier(sdf, removeHs=True) 
		mm = [m for m in mm if m is not None]
		if len(mm) >5000:
			random.seed(a=7)
			mm = random.sample(mm,5000)
		mn = [m.GetProp('_Name') for m in mm if m is not None]
		feats = [factory.GetFeaturesForMol(m) for m in mm if m is not None]
		das = [[]]*len(feats)
		das = [[f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')] for feat in feats]
		das = [unique([x for x in da if da.count(x) > 1]) for da in das]
		mps = [m.GetConformer().GetPositions() for m in mm if m is not None]		
		atms = [m.GetAtoms() for m in mm if m is not None]
		vt = [feature_AA(rns, feats, atms, mps, idx, cnrs, inv, das, fn) for idx in range(len(mps))]
		vdf = pd.DataFrame(vt)
		vdf.index = mn
		return vdf
	aps = mclapply([[sdf, fact, rns, cnrs, fn] for sdf in sdfs], bulk_cids, ncpu=51)
	app = pd.concat(aps)
	app.to_pickle('{}/final3/dude/dataset/mtxs/dud_decoy_apms_{}.pkl'.format(DIR_AP, fact))
	del app


sdfs= [DIR_AP+'/DB/DUDE/all/'+sub+'/decoys_final.sdf' for sub in subdirs]
def sdf2fp_dd(sdf):
	mm = Chem.SDMolSupplier(sdf, removeHs=True) 
	mm = [m for m in mm if m is not None]
	if len(mm) >5000:
		random.seed(a=7)
		mm = random.sample(mm,5000)
	mn = [m.GetProp('_Name') for m in mm if m is not None]
	fp = [list(Chem.rdMolDescriptors.GetMorganFingerprintAsBitVect(m, radius)) for m in mm]
	fps = pd.DataFrame(fp)
	fps.index = mn
	fps['target'] = sdf.split('/')[-2]
	return fps

aps = mclapply(sdfs, sdf2fp_dd, ncpu=51)
app = pd.concat(aps)
app.to_pickle('{}/final3/dude/dataset/mtxs/dud_decoy_fps.pkl'.format(DIR_AP))



c2 =DIR_AP+'/DB/refined-set'
pds1=[pps.split('/')[-1] for pps in glob.glob(c2+'/*')]
c2 = DIR_AP+'/DB/v2020-other-PL'
pds = pds1 + [pps.split('/')[-1] for pps in glob.glob(c2+'/*')]

dstbin = 10
f_list=['','','','_all'] +['_h']*4
st = 20/(inv**dstbin) # 16 bin
cnrs = [st*(inv**g) for g in range(dstbin)]
fnss = [[0,1,2,3,5,6],[0,1,2,3,5,6],[0,1,2,3,5,6,9],[0,1,2,3,4,5,7,8,10],[0,1,2,3,4,6,7,9],[0,1,2,3,4,6,7,9,10,11,12,13,14],[0,1,2,3,4,6,7,9,10,11,12,13,14],[0,1,2,3,4,6,7,9,10,11,12,13,14,15,16,17,18]]
for fact in [0,1,2,3,4,6]:
	factory = ChemicalFeatures.BuildFeatureFactory('{}/etc/BaseFeatures{}.fdef'.format(DIR_AP,f_list[fact]))
	fn = factory.GetFeatureFamilies()
	fns = fnss[fact]
	rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
	def bulk_cids(ssets):
		try:
			pps, fact, rns, cnrs, fn = ssets
			if pps in pds1:didx = 'refined-set'
			else:didx = 'v2020-other-PL'
			m = Chem.SDMolSupplier('{}/DB/{}/{}/{}_ligand.sdf'.format(DIR_AP,didx,pps,pps), removeHs=True)[0]
			if m is None:
				m = Chem.rdmolfiles.MolFromMol2File('{}/DB/{}/{}/{}_ligand.mol2'.format(DIR_AP,didx,pps,pps), removeHs=True)
			feat = [factory.GetFeaturesForMol(m)]
			das = [[]]
			if fact > 0:
				feat = [filter_feats_renew(feat[idx])]
				if fact > 1:
					das = [f.GetAtomIds() for f in feat[idx] if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')]
					das = [unique([x for x in das if das.count(x) > 1])]
					if fact > 4:	
						mps = [m.GetConformer().GetPositions()]
						atms = [m.GetAtoms()]
						vt = feature_AA(rns, feat, atms, mps, idx, cnrs, inv, das, fn)
					else:
						vt = pairLr_dist(feat[idx], rns, cnrs, inv, das[idx], fn)
				else:
					vt = pairLr_dist(feat[idx], rns, cnrs, inv, das[idx], fn)
			else:
				vt = pairLr_dist(feat[idx], rns, cnrs, inv, das[idx], fn)
		except:
			vt = [0]*len(rns)*len(cnrs)
		return tuple(vt)
	aps = mclapply([[pps, fact, rns, cnrs, fn] for pps in pds], bulk_cids, ncpu=40)
	app = pd.DataFrame(aps)
	app.index= pds
	app.to_pickle('{}/final3/pdbbind/dataset/mtxs/pdb_apms_{}.pkl'.format(DIR_AP, fact))
	del app




##### from puchem DB

# cid_list = ~~


cid_list = set(cid_info[0:10].cid)

start_time = time.time()
mm = Chem.SDMolSupplier(sdf[1], removeHs=True) 
print("---{}s seconds---".format(time.time()-start_time))
start_time = time.time()
mm = [m for m in mm if m is not None]
print("---{}s seconds---".format(time.time()-start_time))
start_time = time.time()
mm = [m for m in mm if m.GetPropsAsDict()['PUBCHEM_COMPOUND_CID'] in cids]
print("---{}s seconds---".format(time.time()-start_time))

start_time = time.time()
mn = [m.GetPropsAsDict()['PUBCHEM_COMPOUND_CID'] for m in mm]
print("---{}s seconds---".format(time.time()-start_time))

start_time = time.time()
mn = [m.GetPropsAsDict()['PUBCHEM_COMPOUND_CID'] for m in mm]
print("---{}s seconds---".format(time.time()-start_time))
start_time = time.time()
feats = [factory.GetFeaturesForMol(m) for m in mm]
print("---{}s seconds---".format(time.time()-start_time))

for fact in range(8,9):
	factory = ChemicalFeatures.BuildFeatureFactory('{}/etc/BaseFeatures{}.fdef'.format(DIR_AP,f_list[fact-3]))
	fn = factory.GetFeatureFamilies()
	fns = fnss[fact-3]
	rns = [str(x)+str(y) for x,y in combinations_with_replacement(fns,2)]
	def bulk_cids(ssets):
		sdf, fact, rns, cnrs, fn = ssets
		mm = Chem.SDMolSupplier(sdf, removeHs=True) 
		mn = [m.GetPropsAsDict()['PUBCHEM_COMPOUND_CID'] for m in mm if m is not None]
		feats = [factory.GetFeaturesForMol(m) for m in mm if m is not None]
		das = [[]]*len(feats)
		das = [[f.GetAtomIds() for f in feat if (f.GetFamily()=='Donor')|(f.GetFamily()=='Acceptor')] for feat in feats]
		das = [unique([x for x in da if da.count(x) > 1]) for da in das]
		mps = [m.GetConformer().GetPositions() for m in mm if m is not None]		
		if fact > 8:
			atms = [m.GetAtoms() for m in mm if m is not None]
			if fact==10:
				mpd = [m.GetPropsAsDict() for m in mm if m is not None]
				mmffp = [m['PUBCHEM_MMFF94_PARTIAL_CHARGES'] if 'PUBCHEM_MMFF94_PARTIAL_CHARGES' in m.keys() else '' for m in mpd]
				mmffp = [m.split('\n') if (m!=0) & (m!='') else '' for m in mmffp]
				vt = [feature_AP(rns, feats, atms, mmffp, gap,mps, idx, cnrs, inv, das, fn) for idx in range(len(mpd))]			
			else:
				vt = [feature_AA(rns, feats, atms, mps, idx, cnrs, inv, das, fn) for idx in range(len(mps))]
		else:
			mpd = [m.GetPropsAsDict() for m in mm if m is not None]
			mmffp = [m['PUBCHEM_MMFF94_PARTIAL_CHARGES'] if 'PUBCHEM_MMFF94_PARTIAL_CHARGES' in m.keys() else '' for m in mpd]
			mmffp = [m.split('\n') if (m!=0) & (m!='') else '' for m in mmffp]		
			vt = [feature_PP(rns,feats, mmffp, mps, idx, cnrs, gap, inv, das, fn) for idx in range(len(mpd))]
		vdf = pd.DataFrame(vt)
		vdf.index = mn
		return vdf
	aps = mclapply([[sdf, fact, rns, cnrs, fn] for sdf in sdfs], bulk_cids, ncpu=46)
	app = pd.concat(aps)
	# app2 = pd.read_pickle('{}/final3/pdbbind/dataset/mtxs/apms_{}.pkl'.format(DIR_AP, fact))
	# app = pd.concat([app2,app])
	app.to_pickle('{}/final3/chembl/dataset/mtxs/apms_{}.pkl'.format(DIR_AP, fact))
	del app



start_time = time.time()
mm = [mm for m in mm if m is not None]
mn = [m.GetPropsAsDict()['PUBCHEM_COMPOUND_CID'] for m in mm if m is not None]
feats = [factory.GetFeaturesForMol(m) for m in mm if m is not None]	
print("---{}s seconds---".format(time.time()-start_time))

# last cid idx --> 158925001
use_feats = ['Donor', 'Acceptor', 'NegIonizable', 'PosIonizable', 'Aromatic', 'Hydrophobe']

use_feats = ['Donor', 'Acceptor', 'NegIonizable', 'PosIonizable', 'Halogen', 'Aromatic', 'Hydrophobe']
use_feats = ['Donor', 'Acceptor', 'NegIonizable', 'PosIonizable', 'Halogen', 'Aromatic', 'Hydrophobe', 'LumpedHydrophobe']
# D+A type use / AH / AA / PP / AP
tps = [1,0,1,0,0] # 2
tps = [1,1,0,0,0] # 1
# tps = [1,0,0,0,0] # 0
tps = [0,0,0,0,0] # 0

str_dir = '/spstorage/DB/PUBCHEM/3D_1conf'
str_dir_files = glob.glob(f'{str_dir}/*')
allsets = pd.read_csv(dtip+'/data/data_dataset_3D.csv', index_col = 0)

ttgss = pd.read_csv(DIR_AP+'/ttgs.csv', index_col=0)
spkdfc = pd.read_csv(DIR_AP+'/spkdfc.csv', index_col=0)
spkdfc = pd.read_csv('/spstorage/USERS/gina/Project/AP/DB/InterpretableDTIP/pdbval/cmp2cid2.txt',sep='\t', header =None)

cid_list = set(allsets.cid)
cid_list = set(spkdfc.iloc[:,1])
cid_list = set(datadf.cid)

cid_list = set(spkdfc.cid)
cid_list = set(ff.index)
cid_list = [331300,19767005,441384,469588]
cid_list = [int(cd) for cd in livercid]
cid_list = set(uun.cid)

random.seed(7)
cid_list = random.sample(range(158950000),10**6)
cid_list = cid_list[10**5:10**5+5*10**4]
# cid_list = [cd for cd in cid_list if cd in str_dir_files]
del fpss

del aps 
del zzz

del fff 
import gc

gc.collect()

cidss = [[] for i in range(6358)]
for cid in cid_list:
	if cid < 158950000:
		cidss[int(cid // 25000)].append(cid)


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

start_time = time.time()
fiss = mclapply(infos,FP_from_cid, ncpu= 75) # 100 / 10 855s   / 700 / 20 --> 57m | 3466s
print("---{}s seconds---".format(time.time()-start_time))


# cid_list=[208898.0, 6918155.0, 16220172.0, 3598.0, 5284360.0, 57379345.0, 441362.0, 5702159.0, 44247568.0, 5284373.0, 24748573.0, 2724385.0, 2724387.0, 406563.0, 441383.0, 10223146.0, 213039.0, 54726191.0, 49803313.0, 9556529.0, 3121.0, 680502.0, 92727.0, 5360696.0, 11978813.0, 44093.0, 5309446.0, 6857793.0, 12761153.0, 447043.0, 3652.0, 2049.0, 121926.0, 3246155.0, 16726095.0, 197712.0, 9889366.0, 5280343.0, 66576990.0, 5282398.0, 6240.0, 118598754.0, 5282407.0, 9852519.0, 108144.0, 129138801.0, 468595.0, 30323.0, 2165.0, 11570805.0, 3191.0, 444025.0, 10047612.0, 46926973.0, 56959.0, 3715.0, 3005572.0, 3005573.0, 64646.0, 4746.0, 72331.0, 71821.0, 53389.0, 9817231.0, 154257.0, 676352.0, 12947.0, 4756.0, 166548.0, 9829523.0, 151193.0, 166558.0, 2719.0, 253602.0, 37542.0, 3033767.0, 71549093.0, 126961335.0, 6918328.0, 71172.0, 9910986.0, 5311180.0, 22049997.0, 123596.0, 123601.0, 41684.0, 46220502.0, 4829.0, 148195.0, 4581100.0, 237.0, 275182.0, 17134.0, 2800.0, 3823.0, 135449332.0, 2805.0, 656630.0, 656641.0, 5284613.0, 86271238.0, 446727.0, 53469448.0, 5284616.0, 71496458.0, 16219921.0, 2723601.0, 67356.0, 3038494.0, 3085092.0, 73051434.0, 1069873.0, 9887537.0, 135413553.0, 11681588.0, 22324.0, 49806644.0, 10182969.0, 49843517.0, 9864510.0, 4413.0, 2880.0, 130881.0, 3396.0, 214347.0, 5452.0, 11513676.0, 60749.0, 122201421.0, 132010322.0, 131411.0, 3926.0, 135398745.0, 24775005.0, 33630.0, 11103.0, 10474335.0, 5475.0, 285033.0, 83818.0, 11282283.0, 26987.0, 467825.0, 11167602.0, 3955.0, 49830258.0, 3957.0, 73078.0, 441207.0, 6364534.0, 25151352.0, 121427831.0, 104827.0, 4477.0, 6014.0, 135539077.0, 11955716.0, 6536.0, 16779.0, 5351307.0, 46907787.0, 169870.0, 6918543.0, 25126798.0, 4239764.0, 71587.0, 57389999.0, 10113978.0, 68539.0, 95168.0, 6084.0, 5362119.0, 58298316.0, 135564749.0, 4046.0, 21127119.0, 121304016.0, 45375953.0, 5280723.0, 2733525.0, 10206.0, 51166.0, 1712095.0, 4066.0, 5478883.0, 16195554.0, 3559.0, 2536.0, 16362.0, 146047983.0, 53232.0, 146047984.0, 4594.0, 146047985.0, 6413301.0, 6911989.0, 4091.0]
(71172.0,4631),(3005572.0,3005573),(5309446.0,2782),(64646.0,2165),(5284360.0,2536)

lltss2 = pd.read_pickle(dtip+'/pdbval/AP1_beta.pkl')
lltss2 = pd.concat(ltiss)
lltss2.to_pickle(dtip+'/data/AP1_no_both_left.pkl')

lltss2.to_pickle(dtip+'/covid/AP1_no_both.pkl')
lltss2.to_pickle(dtip+'/covid/AP1_no_both_assay_and_paper.pkl')

lltss = pd.concat([lltss,lltss2],axis=1)
lltss.to_pickle(dtip+'/covid/AP1_both.pkl')

lltss.to_pickle(dtip+'/data/AP1_z_no_both.pkl')
lltss.to_pickle(dtip+'/pdbval/AP1_both.pkl')
lltss.to_pickle(dtip+'/data/AP1_no_both_random2.pkl')


lltss.to_pickle(dtip+'/data/AP1_covid_z_no_both.pkl')

lltss.to_pickle(dtip+'/data/AP1_renewposarm_all.pkl')
lltss.to_pickle(dtip+'/data/AP1_renew2_highes.pkl')


lltss = pd.concat([lltss,lltss2])
ltt.to_pickle(dtip+'/pdbval/AP1_renew22.pkl')

lltss.to_pickle(dtip+'/pdbval/AP1_renew2.pkl')
lltss.to_pickle(dtip+'/pdbval/AP1_beta2.pkl')


lltss.to_pickle('/spstorage/USERS/gina/Project/AP/DB/ChEMBL/AP_2.pkl')
lltss.to_pickle('/spstorage/USERS/gina/Project/AP/DB/ChEMBL/AP_1_renew2.pkl')

lltss.to_pickle('/spstorage/USERS/gina/Project/AP/DB/InterpretableDTIP/pdbval/AP1_renew.pkl')
lltss.to_pickle('/spstorage/USERS/gina/Project/AP/DB/InterpretableDTIP/pdbval/AP1_renew_pdb.pkl')


lltss = pd.concat(ltis)
lltss.to_pickle(dtip+'/human/AP0.pkl')


lltss.to_pickle(dtip+'/human/AP2_renew.pkl')
aa = pd.read_pickle(dtip+'/data/AP2_renew.pkl')

# llts.to_pickle(DIR_AP+'/dti_cadfs80_AP1.pkl')


start_time = time.time()
ltissd = mclapply(infos,AP_from_dud, ncpu= 102) # 100 / 10 855s   / 700 / 20 --> 57m | 3466s
print("---{}s seconds---".format(time.time()-start_time))


start_time = time.time()
ltissd = mclapply(infos,FP_from_dud, ncpu= 102) # 100 / 10 855s   / 700 / 20 --> 57m | 3466s
print("---{}s seconds---".format(time.time()-start_time))

start_time = time.time()
ltisc = mclapply(infos,FP_from_cid, ncpu= 100) # 100 / 10 855s   / 700 / 20 --> 57m | 3466s
print("---{}s seconds---".format(time.time()-start_time))


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


def scfFP_from_cid(info):
    cids, sdf, tps = info
    mm = Chem.SDMolSupplier(sdf, removeHs=True) 
    mm = [m for m in mm if m is not None]
    mm = [m for m in mm if m.GetPropsAsDict()['PUBCHEM_COMPOUND_CID'] in cids]
    mn =  [m.GetPropsAsDict()['PUBCHEM_COMPOUND_CID'] for m in mm]
    scfm = [GetScaffoldForMol(m) for m in mm]
    scf = [list(Chem.rdMolDescriptors.GetMorganFingerprintAsBitVect(cff, radius)) for cff in scfm]
    scfs = pd.DataFrame(scf)
    scfs.index = mn
    return scfs


def FP_from_dud(info):
	cids, sdf, tps = info
	mm = Chem.SDMolSupplier(sdf, removeHs=True) 
	mm = [m for m in mm if m is not None]
	mm = [m for m in mm if m.GetProp('_Name') in cids]
	mn =  [m.GetProp('_Name') for m in mm]
	fp = [list(Chem.rdMolDescriptors.GetMorganFingerprintAsBitVect(m, radius)) for m in mm if m is not None]
	fps = pd.DataFrame(fp)
	fps.index = mn
	return fps



import glob
c2 = DIR_AP+'/DB/DUDE/'
sdfs = glob.glob(c2+'/*.sdf')
sdfs = sorted(sdfs)
subdirs = [x for x in os.walk(c2)][2][1]

d
glob.glob()

# 57 cid --> 48s
# 454 cid --> 84s
# 

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