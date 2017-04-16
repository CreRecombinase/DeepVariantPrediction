
######## Snakemake header ########
import sys; sys.path.insert(0, "/project2/xinhe/yanyul/softwares/Anaconda2/envs/deepvarpred/lib/python3.5/site-packages"); import pickle; snakemake = pickle.loads(b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\x06\x00\x00\x00paramsq\x03csnakemake.io\nParams\nq\x04)\x81q\x05M\xe8\x03a}q\x06(X\x06\x00\x00\x00_namesq\x07}q\x08X\x06\x00\x00\x00windowq\tK\x00N\x86q\nsh\tM\xe8\x03ubX\x04\x00\x00\x00ruleq\x0bX\x18\x00\x00\x00generate_allele1_allele2q\x0cX\x06\x00\x00\x00outputq\rcsnakemake.io\nOutputFiles\nq\x0e)\x81q\x0f(X\x1c\x00\x00\x00data/pranav_test1_allele2.faq\x10X\x1c\x00\x00\x00data/pranav_test1_allele1.faq\x11e}q\x12(X\x02\x00\x00\x00a2q\x13h\x10h\x07}q\x14(h\x13K\x00N\x86q\x15X\x02\x00\x00\x00a1q\x16K\x01N\x86q\x17uh\x16h\x11ubX\t\x00\x00\x00wildcardsq\x18csnakemake.io\nWildcards\nq\x19)\x81q\x1aX\x0c\x00\x00\x00pranav_test1q\x1ba}q\x1c(h\x07}q\x1dX\x04\x00\x00\x00dataq\x1eK\x00N\x86q\x1fsX\x04\x00\x00\x00dataq h\x1bubX\x05\x00\x00\x00inputq!csnakemake.io\nInputFiles\nq")\x81q#X\'\x00\x00\x00data/pranav_test1.formatted.expended.faq$a}q%h\x07}q&sbX\x03\x00\x00\x00logq\'csnakemake.io\nLog\nq()\x81q)}q*h\x07}q+sbX\x06\x00\x00\x00configq,}q-(X\x0b\x00\x00\x00performanceq.}q/(X\x05\x00\x00\x00mode2q0}q1(X\x06\x00\x00\x00scriptq2X\x1e\x00\x00\x00path_to/performance_script1.pyq3X\x06\x00\x00\x00paramsq4}q5(X\n\x00\x00\x00annotationq6X\x16\x00\x00\x00path_to/annotation.bedq7X\x08\x00\x00\x00thresoldq8X\x04\x00\x00\x001e-5q9uuX\x05\x00\x00\x00mode1q:}q;(X\x06\x00\x00\x00scriptq<X\x1e\x00\x00\x00path_to/performance_script1.pyq=X\x06\x00\x00\x00paramsq>X\x15\x00\x00\x00some other input hereq?uuX\x05\x00\x00\x00modelq@}qA(X\x07\x00\x00\x00workdirqBX\x10\x00\x00\x00path_to/work_dirqCX\x04\x00\x00\x00nameqDX\x1d\x00\x00\x00path_to_model/model_name.hdf5qEX\t\x00\x00\x00snakemakeqFX\x11\x00\x00\x00path_to/snakemakeqGuX\x0f\x00\x00\x00genome_assemblyqHX-\x00\x00\x00/project2/xinhe/yanyul/databases/hg19/hg19.faqIX\x05\x00\x00\x00labelqJ}qK(X\x06\x00\x00\x00group1qL}qM(X\x04\x00\x00\x00E082qNK\xc1X\x04\x00\x00\x00E081qOK\xc0uX\x06\x00\x00\x00group2qP}qQ(X\x06\x00\x00\x00NoonanqRK\xc3X\x04\x00\x00\x00E129qSK\xc2uuX\x0b\x00\x00\x00window_sizeqTM\xe8\x03X\x04\x00\x00\x00dataqU}qV(X\x0c\x00\x00\x00pranav_test2qW}qX(X\x06\x00\x00\x00methodqYX\x1c\x00\x00\x00_formatting_pranav.snakemakeqZX\x04\x00\x00\x00nameq[X\x15\x00\x00\x00test/pranav_test2.txtq\\uX\x0c\x00\x00\x00pranav_test1q]}q^(X\x06\x00\x00\x00methodq_X\x1c\x00\x00\x00_formatting_pranav.snakemakeq`X\x04\x00\x00\x00nameqaX\x15\x00\x00\x00test/pranav_test1.txtqbuuuX\t\x00\x00\x00resourcesqccsnakemake.io\nResources\nqd)\x81qe(K\x01K\x01e}qf(h\x07}qg(X\x06\x00\x00\x00_coresqhK\x00N\x86qiX\x06\x00\x00\x00_nodesqjK\x01N\x86qkuhjK\x01hhK\x01ubX\x07\x00\x00\x00threadsqlK\x01ub.')
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

with open(extracted_fa) as infile:
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
