import sys
import pandas as pd
import numpy as np
if 'scripts' not in sys.path:
    sys.path.insert(0, 'scripts')
import my_python
# extract sequence from genome
cmd = '''bedtools getfasta -fi {genome} -bed {input} -fo {out} \
-name -tab'''.format(genome=snakemake.params.assembly, input=snakemake.input.a1, out=snakemake.output.o3)
my_python.myOsSystem(cmd, False)

allele1s = []
allele2s = []
refs = []
ids = []
with open(snakemake.output.o3, 'r') as o:
    for i in o:
        i = i.strip()
        info = i.split('\t')
        allele = info[0].split(':')
        allele1s.append(allele[1])
        allele2s.append(allele[2])
        refs.append(info[1].upper())
        ids.append(allele[0])
data = np.vstack((ids, allele1s, allele2s, refs)).T
dtf = pd.DataFrame(data=data, columns=['ID', 'Allele1', 'Allele2', 'Ref'])
match_allele1 = (dtf.Allele1 == dtf.Ref).sum()
match_allele2 = (dtf.Allele2 == dtf.Ref).sum()
total = dtf.shape[0]
dtf.to_csv(path_or_buf=snakemake.output.o2, sep='\t', index=False)
from snakemake.utils import report

report("""
Variant Formatting Report
=========================

The allele1, allele2, reference allele are summarized here table_
The formatting step is done by {snakemake.params.formatting}

* total variants = {total}
* variants match allele1 = {match_allele1}
* variants match allele2 = {match_allele2}

""", snakemake.output.o1, table=snakemake.output.o2)
