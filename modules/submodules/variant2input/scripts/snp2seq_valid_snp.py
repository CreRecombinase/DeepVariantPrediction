passed = open(snakemake.output.o1, 'w')
notpassed = open(snakemake.output.o2, 'w')

chr_dic = [ 'chr' + str(i) for i in range(1, 24) ]
chr_dic.append('chrX', 'chrY')

with open(snakemake.input.a1, 'r') as o:
    for i in o:
        i = i.strip()
        info = i.split('\t')
        chrm = info[0]
        start = int(info[1])
        end = int(info[2])
        if chrm not in chr_dic:
            notpassed.write(i + '\tNot in Chromosome list\n')
            continue
        if int(end) - int(start) != 1:
            notpassed.write(i + '\tEnd - Start != 1\n')
            continue
        passed.write(i + '\n')
passed.close()
notpassed.close()
