import argparse
parser = argparse.ArgumentParser(prog='merge_rsid_formatting.py', description='''
	Given two SNP lists and merge them according to rsID
	''')
parser.add_argument('--input_list', help='''
    SNP list input filtered by extracted SNPs from database
	''')
parser.add_argument('--input_db', help='''
    SNP list extracted from database
	''')
parser.add_argument('--col_list', type=int, help='''
    Column index of rsID in input_list
	''')
parser.add_argument('--col_db', type=int, help='''
    Column index of rsID in input_db
	''')
parser.add_argument('--output')
args = parser.parse_args()

import h5py
import numpy as np
import pandas as pd
import sys
import gzip

list_dic = {}
counter = 0
with gzip.open(args.input_db,'rb') as f:
    for i in f:
        i = i.decode().split('\t')
        rsID = i[args.col_db - 1]
        list_dic[rsID] = counter
        counter += 1
reorder = []
with gzip.open(args.input_list,'rb') as f:
    for i in f:
        i = i.decode().split('\t')
        rsID = i[args.col_list - 1]
        reorder.append(list_dic[rsID])
table = pd.read_table(args.input_db, header = None, compression='gzip')
table = table.reindex(reorder).reset_index(drop=True)
table2 = pd.read_table(args.input_list, header=None, compression='gzip')
# merged = pd.concat([table2, table], axis=1)
merged = pd.concat([table2, table], axis=1, ignore_index=True)
merged.to_csv(args.output, sep = '\t', compression='gzip', header=False, index=False)
