from apms import *
from utility import *

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--cmp_apm_file',type=str,
                        default='data/apm_compound.csv')
    parser.add_argument('--pck_apm_file',type=str,
                        default='data/apm_protein.csv')
    parser.add_argument('--label_file',type=str,
                        default='data/interaction.csv')
    parser.add_argument('--checkpoint_dir',type=str,
                        default='model_weight')
    parser.add_argument('--log_dir',type=str,
                        default='result')
    args = parser.parse_args()

# Load datasets
dtis = pd.read_csv(args.label_file, header=None)
cmp_apm = pd.read_csv(args.cmp_apm_file, index_col =0)
pck_apm = pd.read_csv(args.pck_apm_file, index_col =0)

# Prepare training, validation, and test datasets
X_train,y_train = setting_dataset(cmp_apm, pck_apm, dtis, 'train',30, 55)
X_dev,y_dev = setting_dataset(cmp_apm, pck_apm, dtis, 'dev',30, 55)
X_test,y_test = setting_dataset(cmp_apm, pck_apm, dtis, 'test',30, 55)

# Initialize the CNN model
torch.manual_seed(7)
model = CNN_model(55, 7, 0.3, 30) # Model with specific parameters (input size, filters, dropout, kernel size)
model.to(device)

# Define optimizer, loss function, and learning rate scheduler
optimizer = torch.optim.Adam(model.parameters(), lr=0.001)
criterion = torch.nn.BCELoss()
scheduler = torch.optim.lr_scheduler.StepLR(optimizer, step_size=10, gamma=0.9)
model.zero_grad()


# Training loop with early stopping
with open(f'{args.log_dir}/log', 'a') as f:
    best_model_wts = copy.deepcopy(model.state_dict())
    best_val_rocauc = 0.5
    early = 0
    epochs_early_stop = 10
    for epoch in range(20):
        losses = []
        accs = []
        # Mini-batch training
        with tqdm(range(0, len(X_train), 80)) as tepoch:
            for i in tepoch:
                tepoch.set_description(f"Epoch {epoch + 1}")
                try:
                    batch_X = X_train[i: i+80]
                    batch_y = y_train[i: i+80]
                except:
                    gc.collect()
                    continue
                # Forward pass and compute loss
                acc, loss, _ = fwd_pass(batch_X, batch_y, model, optimizer, criterion, train=True)
                losses.append(loss.item())
                accs.append(acc)

                # Update progress bar with metrics
                acc_mean = np.array(accs).mean()
                loss_mean = np.array(losses).mean()
                tepoch.set_postfix(loss=loss_mean, accuracy=100. * acc_mean)
            scheduler.step()
        early +=1
        val_rocauc = test_func(model, y_dev, X_dev)
        test_rocauc = test_func(model, y_test, X_test)
        f.write(
            f"{epoch},{round(float(acc), 5)},{round(float(val_rocauc), 5)},{round(float(test_rocauc), 5)}\n")       
        # Early stopping logic
        if early >= epochs_early_stop:
            model.load_state_dict(best_model_wts)
            torch.save(model.state_dict(), f"{arg.checkpoint_dir}/train_model.pth")
            break
        if float(val_rocauc) > best_val_rocauc:
            best_val_rocauc = float(val_rocauc)
            best_model_wts = copy.deepcopy(model.state_dict())
            early=0
        # Save model weights after each epoch
        model.load_state_dict(best_model_wts)
        torch.save(model.state_dict(), f"{arg.checkpoint_dir}/train_model.pth")


