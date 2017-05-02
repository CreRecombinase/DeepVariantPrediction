import argparse
parser = argparse.ArgumentParser(prog='bed2feather.py', description='''
	Given a collection of intersected bed files, summarize the result into
    a single feather file
	''')
parser.add_argument('--inputs', nargs='+', help='''
	A collection of input bed files
	''')
parser.add_argument('--names', nargs='+', help='''
	The corresponding names of the input files
	''')
parser.add_argument('--output')
args = parser.parse_args()

import h5py
import numpy as np
import pandas as pd
import sys
import feather

def get_id(instr):
    e = instr.split(':')
    return(e[0])

total_table = pd.DataFrame()
for inbed in zip(args.inputs, args.names):
    data = pd.read_csv(inbed[0], compression='gzip', delimiter='\t', header=None)
    ids = data[3].apply(get_id)
    rs = data[4]
    sub_table = pd.DataFrame({
				'Varient.ID' : ids.astype(int),
				'rsID' : rs,
				'Annotation' : inbed[1],
				})
    total_table = pd.concat([total_table, sub_table])
feather.write_dataframe(total_table, args.output)
