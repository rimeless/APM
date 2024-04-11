from apms import *
from utility import *

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--cmp_apm_file',type=str,
                        default='data/apm_compound.csv')
    parser.add_argument('--pck_apm_file',type=str,
                        default='data/apm_protein.csv')
    parser.add_argument('--checkpoint_file',type=str,
                        default='model_weight/ap_model.pth')
    parser.add_argument('--result_dir',type=str,
                        default='result')
    args = parser.parse_args()



model = CNN_model(55, 7, 0.3, 30)
model.to(device)
model.load_state_dict(torch.load(args.checkpoint_file))

cmp_apm = pd.read_csv(args.cmp_apm_file, index_col =0)
pck_apm = pd.read_csv(args.pck_apm_file, index_col =0)

apm = np.array(cmp_apm) 
apm = torch.tensor(apm, dtype = torch.float).view(len(apm),1,55,10)
ppock = torch.tensor(np.array(pck_apm, dtype = float).reshape([-1,55,10]), dtype = torch.float)
X = [[apm[i], ppock] for i in range(len(apm))]

probs = []
for idx, X_test in enumerate(X):
    x = [0, 0]
    x[0] = X_test[0].to(device)
    x[1] = X_test[1].to(device)
    res = model(x)
    prob = res[0].item()
    probs.append([idx, prob])



dti_table = pd.DataFrame(probs)
dti_table.columns = ['compound','probability']
dti_table.compound = list(cmp_apm.index)
dti_table.to_csv(f'{args.result_dir}/result.csv')

