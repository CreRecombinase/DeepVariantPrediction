
######## Snakemake header ########
import sys; sys.path.insert(0, "/project2/xinhe/yanyul/softwares/Anaconda2/envs/deepvarpred/lib/python3.5/site-packages"); import pickle; snakemake = pickle.loads(b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\x04\x00\x00\x00ruleq\x03X\x18\x00\x00\x00generate_allele1_allele2q\x04X\t\x00\x00\x00wildcardsq\x05csnakemake.io\nWildcards\nq\x06)\x81q\x07X\x0c\x00\x00\x00pranav_test1q\x08a}q\t(X\x04\x00\x00\x00dataq\nh\x08X\x06\x00\x00\x00_namesq\x0b}q\x0cX\x04\x00\x00\x00dataq\rK\x00N\x86q\x0esubX\x06\x00\x00\x00paramsq\x0fcsnakemake.io\nParams\nq\x10)\x81q\x11M\xe8\x03a}q\x12(X\x06\x00\x00\x00windowq\x13M\xe8\x03h\x0b}q\x14h\x13K\x00N\x86q\x15subX\x06\x00\x00\x00outputq\x16csnakemake.io\nOutputFiles\nq\x17)\x81q\x18(X\x1c\x00\x00\x00data/pranav_test1_allele1.faq\x19X\x1c\x00\x00\x00data/pranav_test1_allele2.faq\x1ae}q\x1b(X\x02\x00\x00\x00a1q\x1ch\x19X\x02\x00\x00\x00a2q\x1dh\x1ah\x0b}q\x1e(h\x1cK\x00N\x86q\x1fh\x1dK\x01N\x86q uubX\t\x00\x00\x00resourcesq!csnakemake.io\nResources\nq")\x81q#(K\x01K\x01e}q$(X\x06\x00\x00\x00_nodesq%K\x01X\x06\x00\x00\x00_coresq&K\x01h\x0b}q\'(h%K\x00N\x86q(h&K\x01N\x86q)uubX\x03\x00\x00\x00logq*csnakemake.io\nLog\nq+)\x81q,}q-h\x0b}q.sbX\x07\x00\x00\x00threadsq/K\x01X\x05\x00\x00\x00inputq0csnakemake.io\nInputFiles\nq1)\x81q2X\'\x00\x00\x00data/pranav_test1.formatted.expended.faq3a}q4h\x0b}q5sbX\x06\x00\x00\x00configq6}q7(X\x05\x00\x00\x00modelq8}q9(X\t\x00\x00\x00snakemakeq:X\x11\x00\x00\x00path_to/snakemakeq;X\x07\x00\x00\x00workdirq<X\x10\x00\x00\x00path_to/work_dirq=X\x04\x00\x00\x00nameq>X\x1d\x00\x00\x00path_to_model/model_name.hdf5q?uX\x0f\x00\x00\x00genome_assemblyq@X-\x00\x00\x00/project2/xinhe/yanyul/databases/hg19/hg19.faqAX\x0b\x00\x00\x00window_sizeqBM\xe8\x03X\x04\x00\x00\x00dataqC}qD(X\x0c\x00\x00\x00pranav_test2qE}qF(X\x06\x00\x00\x00methodqGX\x1c\x00\x00\x00_formatting_pranav.snakemakeqHX\x04\x00\x00\x00nameqIX\x15\x00\x00\x00test/pranav_test2.txtqJuX\x0c\x00\x00\x00pranav_test1qK}qL(X\x06\x00\x00\x00methodqMX\x1c\x00\x00\x00_formatting_pranav.snakemakeqNX\x04\x00\x00\x00nameqOX\x15\x00\x00\x00test/pranav_test1.txtqPuuX\x0b\x00\x00\x00performanceqQ}qR(X\x05\x00\x00\x00mode1qS}qT(X\x06\x00\x00\x00scriptqUX\x1e\x00\x00\x00path_to/performance_script1.pyqVX\x06\x00\x00\x00paramsqWX\x15\x00\x00\x00some other input hereqXuX\x05\x00\x00\x00mode2qY}qZ(X\x06\x00\x00\x00scriptq[X\x1e\x00\x00\x00path_to/performance_script1.pyq\\X\x06\x00\x00\x00paramsq]}q^(X\n\x00\x00\x00annotationq_X\x16\x00\x00\x00path_to/annotation.bedq`X\x08\x00\x00\x00thresoldqaX\x04\x00\x00\x001e-5qbuuuX\x05\x00\x00\x00labelqc}qd(X\x06\x00\x00\x00group1qe}qf(X\x04\x00\x00\x00E082qgK\xc1X\x04\x00\x00\x00E081qhK\xc0uX\x06\x00\x00\x00group2qi}qj(X\x04\x00\x00\x00E129qkK\xc2X\x06\x00\x00\x00NoonanqlK\xc3uuuub.')
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

with open(extracted_fa[0], 'r') as infile:
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
