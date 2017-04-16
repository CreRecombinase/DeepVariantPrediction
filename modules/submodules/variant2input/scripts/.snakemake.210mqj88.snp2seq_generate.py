
######## Snakemake header ########
import sys; sys.path.insert(0, "/project2/xinhe/yanyul/softwares/Anaconda2/envs/deepvarpred/lib/python3.5/site-packages"); import pickle; snakemake = pickle.loads(b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\x04\x00\x00\x00ruleq\x03X\x18\x00\x00\x00generate_allele1_allele2q\x04X\x05\x00\x00\x00inputq\x05csnakemake.io\nInputFiles\nq\x06)\x81q\x07X\'\x00\x00\x00data/pranav_test1.formatted.expended.faq\x08a}q\tX\x06\x00\x00\x00_namesq\n}q\x0bsbX\x06\x00\x00\x00paramsq\x0ccsnakemake.io\nParams\nq\r)\x81q\x0eM\xe8\x03a}q\x0f(X\x06\x00\x00\x00windowq\x10M\xe8\x03h\n}q\x11h\x10K\x00N\x86q\x12subX\t\x00\x00\x00wildcardsq\x13csnakemake.io\nWildcards\nq\x14)\x81q\x15X\x0c\x00\x00\x00pranav_test1q\x16a}q\x17(X\x04\x00\x00\x00dataq\x18h\x16h\n}q\x19X\x04\x00\x00\x00dataq\x1aK\x00N\x86q\x1bsubX\t\x00\x00\x00resourcesq\x1ccsnakemake.io\nResources\nq\x1d)\x81q\x1e(K\x01K\x01e}q\x1f(X\x06\x00\x00\x00_nodesq K\x01X\x06\x00\x00\x00_coresq!K\x01h\n}q"(h K\x00N\x86q#h!K\x01N\x86q$uubX\x07\x00\x00\x00threadsq%K\x01X\x06\x00\x00\x00configq&}q\'(X\x0f\x00\x00\x00genome_assemblyq(X-\x00\x00\x00/project2/xinhe/yanyul/databases/hg19/hg19.faq)X\x05\x00\x00\x00labelq*}q+(X\x06\x00\x00\x00group1q,}q-(X\x04\x00\x00\x00E081q.K\xc0X\x04\x00\x00\x00E082q/K\xc1uX\x06\x00\x00\x00group2q0}q1(X\x06\x00\x00\x00Noonanq2K\xc3X\x04\x00\x00\x00E129q3K\xc2uuX\x04\x00\x00\x00dataq4}q5(X\x0c\x00\x00\x00pranav_test1q6}q7(X\x06\x00\x00\x00methodq8X\x1c\x00\x00\x00_formatting_pranav.snakemakeq9X\x04\x00\x00\x00nameq:X\x15\x00\x00\x00test/pranav_test1.txtq;uX\x0c\x00\x00\x00pranav_test2q<}q=(X\x06\x00\x00\x00methodq>X\x1c\x00\x00\x00_formatting_pranav.snakemakeq?X\x04\x00\x00\x00nameq@X\x15\x00\x00\x00test/pranav_test2.txtqAuuX\x0b\x00\x00\x00performanceqB}qC(X\x05\x00\x00\x00mode1qD}qE(X\x06\x00\x00\x00scriptqFX\x1e\x00\x00\x00path_to/performance_script1.pyqGX\x06\x00\x00\x00paramsqHX\x15\x00\x00\x00some other input hereqIuX\x05\x00\x00\x00mode2qJ}qK(X\x06\x00\x00\x00scriptqLX\x1e\x00\x00\x00path_to/performance_script1.pyqMX\x06\x00\x00\x00paramsqN}qO(X\x08\x00\x00\x00thresoldqPX\x04\x00\x00\x001e-5qQX\n\x00\x00\x00annotationqRX\x16\x00\x00\x00path_to/annotation.bedqSuuuX\x0b\x00\x00\x00window_sizeqTM\xe8\x03X\x05\x00\x00\x00modelqU}qV(X\t\x00\x00\x00snakemakeqWX\x11\x00\x00\x00path_to/snakemakeqXX\x07\x00\x00\x00workdirqYX\x10\x00\x00\x00path_to/work_dirqZX\x04\x00\x00\x00nameq[X\x1d\x00\x00\x00path_to_model/model_name.hdf5q\\uuX\x03\x00\x00\x00logq]csnakemake.io\nLog\nq^)\x81q_}q`h\n}qasbX\x06\x00\x00\x00outputqbcsnakemake.io\nOutputFiles\nqc)\x81qd(X\x1c\x00\x00\x00data/pranav_test1_allele1.faqeX\x1c\x00\x00\x00data/pranav_test1_allele2.faqfe}qg(X\x02\x00\x00\x00a1qhheX\x02\x00\x00\x00a2qihfh\n}qj(hhK\x00N\x86qkhiK\x01N\x86qluubub.')
######## Original script #########
# generate allele1.fa and allele2.fa

allele1_o = open(snakemake.output.a1, 'w')
allele2_o = open(snakemake.output.a2, 'w')
extracted_fa = snakemake.input

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
