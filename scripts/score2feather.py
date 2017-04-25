import argparse
parser = argparse.ArgumentParser(prog='score2feather.py', description='''
	Given model prediction (all labels) extract labels of interest and
    save them as feather file.
	''')
parser.add_argument('--allele1', help='''
	Allele1 prediction
	''')
parser.add_argument('--allele2', help='''
	Allele2 prediction
	''')
parser.add_argument('--label_name', nargs='+', help='''
    The names of each label
    ''')
parser.add_argument('--label_num', nargs='+', type=int, help='''
    The column number of each label (1-based)
    ''')
parser.add_argument('--id')
parser.add_argument('--rsid')
parser.add_argument('--out')
args = parser.parse_args()

import h5py
import numpy as np
import pandas as pd
import sys
if 'scripts' not in sys.path:
    sys.path.insert(0, 'scripts')
import my_python
import feather
ids = np.loadtxt(args.id)
rsid = np.loadtxt(args.rsid)
a1 = my_python.getData(args.allele1, 'y_pred')
a2 = my_python.getData(args.allele2, 'y_pred')
total_table = pd.DataFrame()
for i in range(len(args.label_num)):
    value1 = a1[:,args.label_num[i] - 1]
    value2 = a2[:,args.label_num[i] - 1]
    sub_table = pd.DataFrame({
				'Varient.ID' : ids.astype(int),
				'rsID' : rsid.astype(str)
				'Annotation' : args.label_name[i],
				'Allele1' : value1,
                'Allele2' : value2,
				})
    total_table = pd.concat([total_table, sub_table])
feather.write_dataframe(total_table, args.out)
