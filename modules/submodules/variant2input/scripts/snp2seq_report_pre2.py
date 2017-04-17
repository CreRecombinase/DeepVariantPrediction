import pandas as pd
import numpy as np
# extract sequence from genome

allele1s = []
allele2s = []
refs = []
ids = []
with open(snakemake.input.a, 'r') as o:
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
dtf.to_csv(path_or_buf=snakemake.output.o, sep='\t', index=False)
