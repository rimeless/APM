from apms import *

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_type',type=str,
                        default='compound')
    parser.add_argument('--input_path',type=str,
                        default='data/compounds.sdf')
    parser.add_argument('--out_path',type=str,
                        default='result/res.pkl')
    parser.add_argument('--use_feats',type=str,
                        default='Donor,Acceptor,NegIonizable,PosIonizable,Halogen,Aromatic,Hydrophobe')
    parser.add_argument('--distbin',type=str,
                        default=10)
    parser.add_argument('--pocketN',type=int,
                        default=29)
    args = parser.parse_args()


apmdf = generate_APM(args.input_type, args.input_path, args.dstbin, args.use_feats)
apmdf.to_pickle(args.out_path)
