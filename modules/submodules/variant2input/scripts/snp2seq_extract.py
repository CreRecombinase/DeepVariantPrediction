import my_python
# extract sequence from genome
cmd = '''bedtools getfasta -fi {genome} -bed {input} -fo {out} \
-name -tab'''.format(genome=snakemake.params.assembly, input=snakemake.input.a1, out=snakemake.output.out)
my_python.myOsSystem(cmd, False)
