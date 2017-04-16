import my_python
if params.window % 2 == 0:
	even_flag = 0
	start = params.window / 2 - 1
	end = params.window / 2
else:
	start = (params.window - 1) / 2
	end = -start
	even_flag = 1
midpos = start
cmd = '''cat {in} | awk -F"\\t" '{{ printf $1"\\t"$2-{start}"\\t"$3+{end}; \
 {{ for (i=4;i<=NF;i++) printf "\\t%s", $i }}; print ""}}\' > {out}'''.format(in=input, out=output, start=start, end=end)
my_python.myOsSystem(cmd, False)
