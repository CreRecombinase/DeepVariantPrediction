import pandas as pd
import numpy as np
# extract sequence from genome

allele1s = []
allele2s = []
refs = []
ids = []
allele1s2 = []
allele2s2 = []
refs2 = []
ids2 = []
with open(snakemake.input.a, 'r') as o:
    for i in o:
        i = i.strip()
        info = i.split('\t')
        allele = info[0].split(':')
        allele1s.append(allele[1])
        allele2s.append(allele[2])
        refs.append(info[1].upper())
        ids.append(allele[0])
        if info[1].upper() != allele[1] and info[1].upper() != allele[2]:
            allele1s2.append(allele[1])
            allele2s2.append(allele[2])
            refs2.append(info[1].upper())
            ids2.append(allele[0])
data = np.vstack((ids, allele1s, allele2s, refs)).T
data2 = np.vstack((ids2, allele1s2, allele2s2, refs2)).T
dtf = pd.DataFrame(data=data, columns=['ID', 'Allele1', 'Allele2', 'Ref'])
dtf.to_csv(path_or_buf=snakemake.output.o, sep='\t', index=False)
dtf = pd.DataFrame(data=data2, columns=['ID', 'Allele1', 'Allele2', 'Ref'])
dtf.to_csv(path_or_buf=snakemake.output.o2, sep='\t', index=False)
