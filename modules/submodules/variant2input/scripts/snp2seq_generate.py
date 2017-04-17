# generate allele1.fa and allele2.fa

allele1_o = open(snakemake.output.a1, 'w')
allele2_o = open(snakemake.output.a2, 'w')
extracted_fa = snakemake.input.a1

# some intermediate dependencies
if snakemake.params.window % 2 == 0:
	even_flag = 0
	start = snakemake.params.window / 2 - 1
	end = snakemake.params.window / 2
else:
	start = (snakemake.params.window - 1) / 2
	end = -start
	even_flag = 1
midpos = int(start)
# end

with open(extracted_fa, 'r') as infile:
	for line in infile:
		line = line.strip().upper()
		info = line.split('\t')
		name_chunk = info[0].split(':')
		# print(name_chunk)
		strand = name_chunk[3]
		allele1 = name_chunk[1]
		allele2 = name_chunk[2]

		mid = info[1][midpos]
		allel1_seq = info[1][:midpos] + allele1 + info[1][midpos + 1:]
		allel2_seq = info[1][:midpos] + allele2 + info[1][midpos + 1:]

		## write allele1
		allele1_o.write('\t'.join([info[0], allel1_seq]) + '\n')
		## write allele2
		allele2_o.write('\t'.join([info[0], allel2_seq]) + '\n')

allele1_o.close()
allele2_o.close()
