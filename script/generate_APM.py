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


# Generate Atom-Pair Matrix (APM) DataFrame
apmdf = generate_APM(
    args.input_type,  # The type of input (e.g., 'compound' or 'protein')
    args.input_file,  # Path to the input file (e.g., compound or protein data)
    args.distbin,     # Number of distance bins used in APM calculation
    args.use_feats,   # Features to be used in APM generation (e.g., 'donor,acceptor')
    args.pocketN      # Number of protein pockets to consider (for protein APM)
)

# Save the generated APM to a CSV file
apmdf.to_csv(f'{args.out_path}/apm_{args.input_type}.csv')


# Explanation:
# 1. `generate_APM` is a function that processes input data to create an Atom-Pair Matrix (APM).
#    - This matrix encodes spatial relationships and feature interactions (e.g., donor-acceptor pairs) 
#      for compounds or proteins.
#    - Parameters like `distbin` determine the granularity of distance bins, while `use_feats` specifies
#      the chemical features to include.
# 2. The resulting APM DataFrame (`apmdf`) is saved to a CSV file in the specified output directory (`args.out_path`).
# 3. The file is named `apm_<input_type>.csv`, where `<input_type>` distinguishes between compound and protein APMs.
