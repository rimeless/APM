**Atom-Pair Map**
=============

![figure1](https://github.com/rimeless/APM/assets/48581374/f0dcb2a3-6785-4988-8fea-b4c9b88b56d1)

A improved molecule representation method that captures both the physicochemical properties and the spatial arrangements of atoms

**Installation**
-------------
```
git clone https://github.com/rimeless/APM.git
cd APM
conda env create -f env.yml
conda activate apm
```

**Genarate APM**
-------------
```
python generate_APM.py --input_type compound --file_type sdf
```

**Training**
-------------
```
python train_model.py --cmp_path data/compound --pck_path data/pocket --labels data/interaction.csv --save_path model_weight
```

**Predicting**
-------------
```
python main.py --cmp_path data/compound --pck_path data/pocket --model_name model_weight/ap_model.pth --result_path result
```
