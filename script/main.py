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


# Initialize the convolutional neural network (CNN) model
model = CNN_model(55, 7, 0.3, 30)  # Model with input size 55, 7 filters, dropout rate 0.3, and kernel size 30
model.to(device)
model.load_state_dict(torch.load(args.checkpoint_file)) # Load pre-trained weights into the model

# Load compound and protein APM 
cmp_apm = pd.read_csv(args.cmp_apm_file, index_col =0)
pck_apm = pd.read_csv(args.pck_apm_file, index_col =0)

# Convert compound APM data into a 4D tensor (for CNN input)
apm = np.array(cmp_apm) 
apm = torch.tensor(apm, dtype = torch.float).view(len(apm),1,55,10)

# Convert protein pocket APM data into a tensor
pocket = torch.tensor(np.array(pck_apm, dtype = float).reshape([-1,55,10]), dtype = torch.float)
X = [[apm[i], pocket] for i in range(len(apm))]

# Make predictions using the trained model
probs = []
for idx, X_test in enumerate(X):
    x = [0, 0]
    x[0] = X_test[0].to(device)
    x[1] = X_test[1].to(device)
    res = model(x)
    prob = res[0].item()
    probs.append([idx, prob])

# Create a DataFrame to store prediction results
dti_table = pd.DataFrame(probs)
dti_table.columns = ['compound','probability']
dti_table.compound = list(cmp_apm.index)
dti_table.to_csv(f'{args.result_dir}/result.csv')

