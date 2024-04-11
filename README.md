**Atom-Pair Map**
=============

![figure1](https://github.com/rimeless/APM/assets/48581374/f0dcb2a3-6785-4988-8fea-b4c9b88b56d1)

A improved molecule representation method that captures both the physicochemical properties and the spatial arrangements of atoms

**Installation**
-------------
```
git clone https://github.com/rimeless/APM.git
cd APM
conda env create -f env.yaml
conda activate apm
```

**Genarate APM**
-------------
To generate an Atom Pair Map (APM), input your 3D structure files (SDF or PDB) for the compound and protein, selecting the atom pair types and the number of distance bins. Configure these parameters to tailor the APM generation to your specific research needs.

```
python script/generate_APM.py --input_type compound --input_file data/compound.sdf --out_path data --distbin 10
```

**Training with new data**
-------------
Train APNet using your data, with compound and PDB ID as dataframe indices.
```
python script/train.py --cmp_apm_file data/apm_compound.csv --pck_apm_path data/apm_protein.csv --labels data/interaction.csv --checkpoint_path model_weight
```

**Predicting with saved_checkpoint**
-------------
Obtain prediction results as a CSV file detailing interactions between combination of given compounds and proteins.
```
python script/main.py --cmp_apm_file data/apm_compound.csv --pck_apm_file data/apm_protein.csv --checkpoint_file model_weight/ap_model.pth --result_dir result/result.csv
```
