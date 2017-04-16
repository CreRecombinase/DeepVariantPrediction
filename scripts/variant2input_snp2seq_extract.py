import my_python
# extract sequence from genome
cmd = '''bedtools getfasta -fi {genome} -bed {in} -fo {out} \
-name -tab'''.format(genome=params.assembly, in=input, out=output)
my_python.myOsSystem(cmd, False)
