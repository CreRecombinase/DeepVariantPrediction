import pandas as pd
import numpy as np

behavior = []
with open(snakemake.input.fasta) as o:
    for i in o:
        i = i.strip()
        info = i.split('\t')
        allele = info[0].split(':')
        if info[1].upper() != allele[1] and info[1].upper() != allele[2]:
            behavior.append(info[1].upper())
            continue
        if info[1].upper() == allele[2]:
            behavior.append(1)
            continue
        behavior.append(0)

good = open(snakemake.output.o, 'w')
bad = open(snakemake.output.o2, 'w')
counter = 0
with open(snakemake.input.a, 'r') as o:
    for i in o:
        i = i.strip()
        beh = behavior[counter]
        counter += 1

        if type(beh) is str:
            bad.write(i + '\t' + beh + '\n')
        elif beh == 0:
            good.write(i + '\n')
        elif beh == 1:
            info = i.split('\t')
            allele = info[3].split(':')
            allele = ':'.join([allele[0], allele[2], allele[1], allele[3]])
            i = '\t'.join(info[:3] + [allele] + info[4:])
            good.write(i + '\n')
good.close()
bad.close()
