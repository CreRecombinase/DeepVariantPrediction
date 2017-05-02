import sys
if 'scripts' not in sys.path:
    sys.path.insert(0, 'scripts')
import my_python
if snakemake.params.window % 2 == 0:
	even_flag = 0
	start = snakemake.params.window / 2 - 1
	end = snakemake.params.window / 2
else:
	start = (snakemake.params.window - 1) / 2
	end = start
	even_flag = 1
midpos = start
cmd = '''cat {input} | awk -F"\\t" '{{ printf $1"\\t"$2-{start}"\\t"$3+{end}; \
 {{ for (i=4;i<=NF;i++) printf "\\t%s", $i }}; print ""}}\' > {out}'''.format(input=snakemake.input.a1, out=snakemake.output.out, start=start, end=end)
my_python.myOsSystem(cmd, False)
