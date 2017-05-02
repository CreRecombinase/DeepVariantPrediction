def read_size(size_file):
    dic = {}
    with open(size_file, 'r') as o:
        for i in o:
            i = i.strip()
            i = i.split('\t')
            dic[i[0]] = int(i[1])
    return dic

if snakemake.params.window % 2 == 0:
	end_expand = snakemake.params.window / 2
else:
	end_expand = (snakemake.params.window - 1) / 2

size_dic = read_size(snakemake.params.size)
passed = open(snakemake.output.o1, 'w')
notpassed = open(snakemake.output.o2, 'w')

with open(snakemake.input.a1, 'r') as o:
    for i in o:
        i = i.strip()
        info = i.split('\t')
        chrm = info[0]
        end = int(info[2])
        if size_dic[chrm] < end:
            notpassed.write(i + '\tInvalid Site\n')
            continue
        if size_dic[chrm] < end + end_expand:
            notpassed.write(i + '\tWindow Is Too big\n')
            continue
        passed.write(i + '\n')
passed.close()
notpassed.close()
