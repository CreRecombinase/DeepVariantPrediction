
######## Snakemake header ########
import sys; sys.path.insert(0, "/project2/xinhe/yanyul/softwares/Anaconda2/envs/deepvarpred/lib/python3.5/site-packages"); import pickle; snakemake = pickle.loads(b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\t\x00\x00\x00wildcardsq\x03csnakemake.io\nWildcards\nq\x04)\x81q\x05X\x0c\x00\x00\x00pranav_test1q\x06a}q\x07(X\x04\x00\x00\x00dataq\x08h\x06X\x06\x00\x00\x00_namesq\t}q\nX\x04\x00\x00\x00dataq\x0bK\x00N\x86q\x0csubX\x03\x00\x00\x00logq\rcsnakemake.io\nLog\nq\x0e)\x81q\x0f}q\x10h\t}q\x11sbX\x06\x00\x00\x00outputq\x12csnakemake.io\nOutputFiles\nq\x13)\x81q\x14(X\x1c\x00\x00\x00data/pranav_test1_allele1.faq\x15X\x1c\x00\x00\x00data/pranav_test1_allele2.faq\x16e}q\x17(h\t}q\x18(X\x02\x00\x00\x00a1q\x19K\x00N\x86q\x1aX\x02\x00\x00\x00a2q\x1bK\x01N\x86q\x1cuh\x19h\x15h\x1bh\x16ubX\t\x00\x00\x00resourcesq\x1dcsnakemake.io\nResources\nq\x1e)\x81q\x1f(K\x01K\x01e}q (X\x06\x00\x00\x00_nodesq!K\x01X\x06\x00\x00\x00_coresq"K\x01h\t}q#(h"K\x00N\x86q$h!K\x01N\x86q%uubX\x07\x00\x00\x00threadsq&K\x01X\x06\x00\x00\x00paramsq\'csnakemake.io\nParams\nq()\x81q)M\xe8\x03a}q*(X\x06\x00\x00\x00windowq+M\xe8\x03h\t}q,h+K\x00N\x86q-subX\x05\x00\x00\x00inputq.csnakemake.io\nInputFiles\nq/)\x81q0X\'\x00\x00\x00data/pranav_test1.formatted.expended.faq1a}q2(h\x19h1h\t}q3h\x19K\x00N\x86q4subX\x04\x00\x00\x00ruleq5X\x18\x00\x00\x00generate_allele1_allele2q6X\x06\x00\x00\x00configq7}q8(X\x05\x00\x00\x00modelq9}q:(X\x04\x00\x00\x00nameq;X\x1d\x00\x00\x00path_to_model/model_name.hdf5q<X\x07\x00\x00\x00workdirq=X\x10\x00\x00\x00path_to/work_dirq>X\t\x00\x00\x00snakemakeq?X\x11\x00\x00\x00path_to/snakemakeq@uX\x05\x00\x00\x00labelqA}qB(X\x06\x00\x00\x00group2qC}qD(X\x06\x00\x00\x00NoonanqEK\xc3X\x04\x00\x00\x00E129qFK\xc2uX\x06\x00\x00\x00group1qG}qH(X\x04\x00\x00\x00E082qIK\xc1X\x04\x00\x00\x00E081qJK\xc0uuX\x04\x00\x00\x00dataqK}qL(X\x0c\x00\x00\x00pranav_test2qM}qN(X\x04\x00\x00\x00nameqOX\x15\x00\x00\x00test/pranav_test2.txtqPX\x06\x00\x00\x00methodqQX\x1c\x00\x00\x00_formatting_pranav.snakemakeqRuX\x0c\x00\x00\x00pranav_test1qS}qT(X\x04\x00\x00\x00nameqUX\x15\x00\x00\x00test/pranav_test1.txtqVX\x06\x00\x00\x00methodqWX\x1c\x00\x00\x00_formatting_pranav.snakemakeqXuuX\x0b\x00\x00\x00performanceqY}qZ(X\x05\x00\x00\x00mode2q[}q\\(X\x06\x00\x00\x00paramsq]}q^(X\n\x00\x00\x00annotationq_X\x16\x00\x00\x00path_to/annotation.bedq`X\x08\x00\x00\x00thresoldqaX\x04\x00\x00\x001e-5qbuX\x06\x00\x00\x00scriptqcX\x1e\x00\x00\x00path_to/performance_script1.pyqduX\x05\x00\x00\x00mode1qe}qf(X\x06\x00\x00\x00paramsqgX\x15\x00\x00\x00some other input hereqhX\x06\x00\x00\x00scriptqiX\x1e\x00\x00\x00path_to/performance_script1.pyqjuuX\x0f\x00\x00\x00genome_assemblyqkX-\x00\x00\x00/project2/xinhe/yanyul/databases/hg19/hg19.faqlX\x0b\x00\x00\x00window_sizeqmM\xe8\x03uub.')
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
