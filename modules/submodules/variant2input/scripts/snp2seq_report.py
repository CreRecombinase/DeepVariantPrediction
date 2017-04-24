from snakemake.utils import report
import pandas as pd

dtf = pd.read_csv(snakemake.input.a, sep = '\t')
match_allele1 = (dtf.Allele1 == dtf.Ref).sum()
match_allele2 = (dtf.Allele2 == dtf.Ref).sum()
total = dtf.shape[0]
notmatch = total - match_allele1 - match_allele2

report("""
Variant Formatting Report
=========================

The unmatched allele1, allele2, reference allele are summarized here table_
The formatting step is done by `{snakemake.params.formatting}`

* total variants = {total}
* variants match allele1 = {match_allele1}
* variants match allele2 = {match_allele2}
* number of unmatched variants = {notmatch}

""", snakemake.output.o, table=snakemake.input.a2)
