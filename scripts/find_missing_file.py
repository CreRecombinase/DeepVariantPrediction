import argparse
parser = argparse.ArgumentParser(prog='find_missing_file.py', description='''
    Given the filename and the naming convention, find the missing file by ID
''')
parser.add_argument('--prefix', help='''
    Prefix of filename (exclude dirname)
''')
parser.add_argument('--suffix')
parser.add_argument('--dir', help = '''
    The directory to explore
''')
parser.add_argument('--target_id_range_from', type=int)
parser.add_argument('--target_id_range_to', type=int)
args = parser.parse_args()

import glob
import ntpath
import re
import os

names = glob.glob(args.dir + os.sep + '*')
exists = []
# print(names)
for n in names:
    # print(n)    
    n = ntpath.basename(n)
    n = re.sub(args.prefix, '', n)
    n = re.sub(args.suffix, '', n)
    try:
        if n == '':
            continue
        idx = int(n)
        # print(idx)
    except ValueError:
        continue
    exists.append(idx)

targets = [ i for i in range(args.target_id_range_from, args.target_id_range_to + 1) ]
missing = []
for n in targets:
    if n not in exists:
        missing.append(str(n))
print(','.join(missing))
