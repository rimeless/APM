from apms import *

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_type',type=str,
                        default='compound')
    parser.add_argument('--input_file',type=str,
                        default='data/compounds.sdf')
    parser.add_argument('--out_path',type=str,
                        default='data')
    parser.add_argument('--use_feats',type=str,
                        default='Donor,Acceptor,NegIonizable,PosIonizable,Halogen,Aromatic,Hydrophobe')
    parser.add_argument('--distbin',type=int,
                        default=10)
    parser.add_argument('--pocketN',type=int,
                        default=29)
    args = parser.parse_args()


apmdf = generate_APM(args.input_type, args.input_file, args.distbin, args.use_feats, args.pocketN)
apmdf.to_csv(f'{args.out_path}/apm_{args.input_type}.csv')
