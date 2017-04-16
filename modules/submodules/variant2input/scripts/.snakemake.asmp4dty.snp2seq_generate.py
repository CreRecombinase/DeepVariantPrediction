
######## Snakemake header ########
import sys; sys.path.insert(0, "/project2/xinhe/yanyul/softwares/Anaconda2/envs/deepvarpred/lib/python3.5/site-packages"); import pickle; snakemake = pickle.loads(b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\t\x00\x00\x00resourcesq\x03csnakemake.io\nResources\nq\x04)\x81q\x05(K\x01K\x01e}q\x06(X\x06\x00\x00\x00_nodesq\x07K\x01X\x06\x00\x00\x00_namesq\x08}q\t(X\x06\x00\x00\x00_coresq\nK\x00N\x86q\x0bh\x07K\x01N\x86q\x0cuh\nK\x01ubX\t\x00\x00\x00wildcardsq\rcsnakemake.io\nWildcards\nq\x0e)\x81q\x0fX\x0c\x00\x00\x00pranav_test1q\x10a}q\x11(h\x08}q\x12X\x04\x00\x00\x00dataq\x13K\x00N\x86q\x14sX\x04\x00\x00\x00dataq\x15h\x10ubX\x06\x00\x00\x00configq\x16}q\x17(X\x05\x00\x00\x00labelq\x18}q\x19(X\x06\x00\x00\x00group1q\x1a}q\x1b(X\x04\x00\x00\x00E081q\x1cK\xc0X\x04\x00\x00\x00E082q\x1dK\xc1uX\x06\x00\x00\x00group2q\x1e}q\x1f(X\x04\x00\x00\x00E129q K\xc2X\x06\x00\x00\x00Noonanq!K\xc3uuX\x05\x00\x00\x00modelq"}q#(X\x07\x00\x00\x00workdirq$X\x10\x00\x00\x00path_to/work_dirq%X\t\x00\x00\x00snakemakeq&X\x11\x00\x00\x00path_to/snakemakeq\'X\x04\x00\x00\x00nameq(X\x1d\x00\x00\x00path_to_model/model_name.hdf5q)uX\x0f\x00\x00\x00genome_assemblyq*X-\x00\x00\x00/project2/xinhe/yanyul/databases/hg19/hg19.faq+X\x0b\x00\x00\x00window_sizeq,M\xe8\x03X\x0b\x00\x00\x00performanceq-}q.(X\x05\x00\x00\x00mode2q/}q0(X\x06\x00\x00\x00scriptq1X\x1e\x00\x00\x00path_to/performance_script1.pyq2X\x06\x00\x00\x00paramsq3}q4(X\x08\x00\x00\x00thresoldq5X\x04\x00\x00\x001e-5q6X\n\x00\x00\x00annotationq7X\x16\x00\x00\x00path_to/annotation.bedq8uuX\x05\x00\x00\x00mode1q9}q:(X\x06\x00\x00\x00scriptq;X\x1e\x00\x00\x00path_to/performance_script1.pyq<X\x06\x00\x00\x00paramsq=X\x15\x00\x00\x00some other input hereq>uuX\x04\x00\x00\x00dataq?}q@(X\x0c\x00\x00\x00pranav_test2qA}qB(X\x04\x00\x00\x00nameqCX\x15\x00\x00\x00test/pranav_test2.txtqDX\x06\x00\x00\x00methodqEX\x1c\x00\x00\x00_formatting_pranav.snakemakeqFuX\x0c\x00\x00\x00pranav_test1qG}qH(X\x04\x00\x00\x00nameqIX\x15\x00\x00\x00test/pranav_test1.txtqJX\x06\x00\x00\x00methodqKX\x1c\x00\x00\x00_formatting_pranav.snakemakeqLuuuX\x07\x00\x00\x00threadsqMK\x01X\x04\x00\x00\x00ruleqNX\x18\x00\x00\x00generate_allele1_allele2qOX\x03\x00\x00\x00logqPcsnakemake.io\nLog\nqQ)\x81qR}qSh\x08}qTsbX\x06\x00\x00\x00outputqUcsnakemake.io\nOutputFiles\nqV)\x81qW(X\x1c\x00\x00\x00data/pranav_test1_allele1.faqXX\x1c\x00\x00\x00data/pranav_test1_allele2.faqYe}qZ(X\x02\x00\x00\x00a2q[hYh\x08}q\\(X\x02\x00\x00\x00a1q]K\x00N\x86q^h[K\x01N\x86q_uh]hXubX\x06\x00\x00\x00paramsq`csnakemake.io\nParams\nqa)\x81qbM\xe8\x03a}qc(h\x08}qdX\x06\x00\x00\x00windowqeK\x00N\x86qfsheM\xe8\x03ubX\x05\x00\x00\x00inputqgcsnakemake.io\nInputFiles\nqh)\x81qiX\'\x00\x00\x00data/pranav_test1.formatted.expended.faqja}qk(h\x08}qlh]K\x00N\x86qmsh]hjubub.')
######## Original script #########
# generate allele1.fa and allele2.fa

allele1_o = open(snakemake.output.a1, 'w')
allele2_o = open(snakemake.output.a2, 'w')
extracted_fa = snakemake.input.a1

# some intermediate dependencies
reverse_dic = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
if snakemake.params.window % 2 == 0:
	even_flag = 0
	start = snakemake.params.window / 2 - 1
	end = snakemake.params.window / 2
else:
	start = (snakemake.params.window - 1) / 2
	end = -start
	even_flag = 1
midpos = start
print(midpos)
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
		if strand == '-':
			mid = info[1][midpos + even_flag * 2]
			allel1_seq = info[1][:midpos + even_flag * 2] + allele1 + info[1][midpos + even_flag * 2 + 1:]
			allel2_seq = info[1][:midpos + even_flag * 2] + allele2 + info[1][midpos + even_flag * 2 + 1:]
			# print(mid)
		else:
			mid = info[1][midpos]
			allel1_seq = info[1][:midpos] + allele1 + info[1][midpos + 1:]
			allel2_seq = info[1][:midpos] + allele2 + info[1][midpos + 1:]

		## generate checking info for allele1
		if mid == allele1:
			check1 = '0'
		elif mid == reverse_dic[allele1]:
			check1 = '1'
		else:
			check1 = '2'

		## generate checking info for allele2
		if mid == allele2:
			check2 = '0'
		elif mid == reverse_dic[allele2]:
			check2 = '1'
		else:
			check2 = '2'

		## write allele1
		allele1_o.write('\t'.join([info[0], allel1_seq, check1]) + '\n')
		## write allele2
		allele2_o.write('\t'.join([info[0], allel2_seq, check2]) + '\n')
