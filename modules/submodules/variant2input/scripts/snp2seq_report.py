from snakemake.utils import report
import pandas as pd

# derived from http://stackoverflow.com/questions/845058/how-to-get-line-count-cheaply-in-python
def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    try:
        return i + 1
    except NameError:
        return 0
# end

outof_size = file_len(snakemake.input.size)
notmatch = file_len(snakemake.input.reorder)
passed = file_len(snakemake.input.passed)
total =  file_len(snakemake.input.total)
attention = 'nothing'
if total != outof_size + notmatch + passed:
    attention = 'total != passed + not_matched + out_of_genome'

report("""
Variant Formatting Report
=========================

The out of range variants are summarized here table1_
The unmatched allele1, allele2, reference allele are summarized here table2_
The formatting step is done by `{snakemake.params.formatting}`

* total variants = {total}
* variants that passed QC = {passed}
* variants that are out of genome range = {outof_size}
* variants that cannot match reference allele = {notmatch}

ATTENTION: {attention}

""", snakemake.output.o, table1 = snakemake.input.size, table2=snakemake.input.reorder)
