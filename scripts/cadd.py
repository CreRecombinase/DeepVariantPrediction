import argparse
parser = argparse.ArgumentParser(prog='cadd.py', description='''
	Given a bed file that contains variant information, output CADD
    score for each variant summarized in a feather file
	''')
parser.add_argument('--cadd_path', help='''
	Path to CADD score database
	''')
parser.add_argument('--input_bed')
parser.add_argument('--output')
args = parser.parse_args()

import h5py
import numpy as np
import pandas as pd
import feather
import sys
if 'scripts' not in sys.path:
    sys.path.insert(0, 'scripts')
import my_python

ids = []
rss = []
phreds = []

with open(args.input_bed, 'r') as f:
    for snp in f:
        snp = snp.strip()
        snp = snp.split('\t')
        chrm = snp[0]
        pos1 = snp[2]
        info = snp[3].split(':')
        ref = info[1]
        alt = info[2]
        idx = info[0]
        rs = snp[4]
        query = '{chr}:{pos}-{pos}'.format(chr=re.sub('chr', '', chrm), pos=pos1)
        cmd = 'tabix {database} {query}'.format(database=args.cadd_path, query=query)
        result = my_python.mySubprocess(cmd, False)

        if result.decode("utf-8") == '':
            pass
        else:
            for i in result.decode("utf-8").split('\n'):
                detail = i.split('\t')
                if detail[2] != ref:
                    my_python.eprint('Wrong ref at {chr}:{pos} {rs}'.format(chr=chrm, pos=pos1, rs=rs))
                    break
                if detail[3] != alt:
                    continue
                else:
                    phred = float(detail[-1])
                    ids.append(idx)
                    rss.append(rs)
                    phreds.append(phred)

total_table = pd.DataFrame({
			'Varient.ID' : ids.astype(int),
			'rsID' : rss,
			'CADD.Phred' : phreds,
			})
feather.write_dataframe(total_table, args.output)
