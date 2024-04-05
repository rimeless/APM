from utility import *

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--cmp_apm_file',type=str,
                        default='data/apm_compound.csv')
    parser.add_argument('--pck_apm_file',type=str,
                        default='data/apm_protein.csv')
    parser.add_argument('--checkpoint_dir',type=str,
                        default='model_weight')
    parser.add_argument('--log_dir',type=str,
                        default='result')
    args = parser.parse_args()









model = CNN_model(55, 7, 0.3, 30)
model.to(device)
model.load_state_dict(torch.load(args.checkpoint_file, map_location=torch.device('cpu')))

cmp_apm = pd.read_csv(cmp_apm_file)
pck_apm = pd.read_csv(pck_apm_file)

apm = np.array(cmp_apm) 
apm = torch.tensor(apm, dtype = torch.float).view(len(vdf),1,rn,10)
ppock = torch.tensor(np.array(pck_apm, dtype = float).reshape([-1,rn,10]), dtype = torch.float)
X = [[apm[i], ppock] for i in range(len(vdf))]


EPOCHS = 20
for BATCH_SIZE in [512,256,128]:
    for pns in [7]:
        for dropout in [0.1,0.3,0.5]:
            for lrt in [1e-4, 1e-3]:
                torch.manual_seed(123)
                model = CNN_model(rn, pns, dropout, pockn)
                model.to(device)
                optimizer = torch.optim.Adam(model.parameters(), lr=lrt)
                criterion = torch.nn.BCELoss()
                # criterion = nn.CrossEntropyLoss()
                scheduler = torch.optim.lr_scheduler.StepLR(optimizer, step_size=10, gamma=0.9)
                model.zero_grad()
                MODEL_NAME = 'renew3_ori_ap{}_lr{}_dr{}_{}pock_{}pns_{}_{}batch'.format(str(idx),str(lrt),str(dropout),str(pockn),str(pns),str(norm),str(BATCH_SIZE))
                with open(f'{arg.log_dir}/log', 'a') as f:
                    best_model_wts = copy.deepcopy(model.state_dict())
                    best_val_rocauc = 0.5
                    early = 0
                    epochs_early_stop = 10
                    for epoch in range(EPOCHS):
                        losses = []
                        accs = []
                        with tqdm(range(0, len(X_train), BATCH_SIZE)) as tepoch:
                            for i in tepoch:
                                tepoch.set_description(f"Epoch {epoch + 1}")
                                try:
                                    batch_X = X_train[i: i+BATCH_SIZE]
                                    batch_y = y_train[i: i+BATCH_SIZE]
                                except:
                                    gc.collect()
                                    continue
                                acc, loss, _ = fwd_pass(batch_X, batch_y, model, train=True)
                                losses.append(loss.item())
                                accs.append(acc)
                                acc_mean = np.array(accs).mean()
                                loss_mean = np.array(losses).mean()
                                tepoch.set_postfix(loss=loss_mean, accuracy=100. * acc_mean)
                            scheduler.step()
                        early +=1
                        val_rocauc = test_func(model, y_dev, X_dev)
                        test_rocauc = test_func(model, y_test, X_test)
                        test_seen_rocauc = test_func(model, y_seen_test, X_seen_test)
                        test_unseen_rocauc = test_func(model, y_unseen_test, X_unseen_test)
                        f.write(
                            f"{MODEL_NAME},{epoch},{rn},{dropout},{pockn},{pns},{round(float(acc), 5)},{round(float(val_rocauc), 5)},{round(float(test_rocauc), 5)},{round(float(test_seen_rocauc), 5)},{round(float(test_unseen_rocauc), 5)}\n")
                        print(str(val_rocauc))
                        print(str(test_rocauc))
                        print(str(test_seen_rocauc))
                        print(str(test_unseen_rocauc))            
                        if early >= epochs_early_stop:
                            model.load_state_dict(best_model_wts)
                            torch.save(model.state_dict(),
                                "{}/data/{}.pth".format(dtip,MODEL_NAME))
                            break
                        if float(val_rocauc) > best_val_rocauc:
                            best_val_rocauc = float(val_rocauc)
                            best_model_wts = copy.deepcopy(model.state_dict())
                            early=0




probdf = pd.DataFrame(probs)
probdf.columns = ['compound','probability']
dti_table.to_csv(f'{result_file}/result.csv')

