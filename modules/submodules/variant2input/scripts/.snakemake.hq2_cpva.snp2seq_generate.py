
######## Snakemake header ########
import sys; sys.path.insert(0, "/project2/xinhe/yanyul/softwares/Anaconda2/envs/deepvarpred/lib/python3.5/site-packages"); import pickle; snakemake = pickle.loads(b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\x06\x00\x00\x00paramsq\x03csnakemake.io\nParams\nq\x04)\x81q\x05M\xe8\x03a}q\x06(X\x06\x00\x00\x00windowq\x07M\xe8\x03X\x06\x00\x00\x00_namesq\x08}q\th\x07K\x00N\x86q\nsubX\x05\x00\x00\x00inputq\x0bcsnakemake.io\nInputFiles\nq\x0c)\x81q\rX\'\x00\x00\x00data/pranav_test1.formatted.expended.faq\x0ea}q\x0fh\x08}q\x10sbX\x04\x00\x00\x00ruleq\x11X\x18\x00\x00\x00generate_allele1_allele2q\x12X\x07\x00\x00\x00threadsq\x13K\x01X\x06\x00\x00\x00configq\x14}q\x15(X\x04\x00\x00\x00dataq\x16}q\x17(X\x0c\x00\x00\x00pranav_test2q\x18}q\x19(X\x04\x00\x00\x00nameq\x1aX\x15\x00\x00\x00test/pranav_test2.txtq\x1bX\x06\x00\x00\x00methodq\x1cX\x1c\x00\x00\x00_formatting_pranav.snakemakeq\x1duX\x0c\x00\x00\x00pranav_test1q\x1e}q\x1f(X\x04\x00\x00\x00nameq X\x15\x00\x00\x00test/pranav_test1.txtq!X\x06\x00\x00\x00methodq"X\x1c\x00\x00\x00_formatting_pranav.snakemakeq#uuX\x0b\x00\x00\x00window_sizeq$M\xe8\x03X\x0f\x00\x00\x00genome_assemblyq%X-\x00\x00\x00/project2/xinhe/yanyul/databases/hg19/hg19.faq&X\x05\x00\x00\x00labelq\'}q((X\x06\x00\x00\x00group2q)}q*(X\x04\x00\x00\x00E129q+K\xc2X\x06\x00\x00\x00Noonanq,K\xc3uX\x06\x00\x00\x00group1q-}q.(X\x04\x00\x00\x00E082q/K\xc1X\x04\x00\x00\x00E081q0K\xc0uuX\x05\x00\x00\x00modelq1}q2(X\x04\x00\x00\x00nameq3X\x1d\x00\x00\x00path_to_model/model_name.hdf5q4X\t\x00\x00\x00snakemakeq5X\x11\x00\x00\x00path_to/snakemakeq6X\x07\x00\x00\x00workdirq7X\x10\x00\x00\x00path_to/work_dirq8uX\x0b\x00\x00\x00performanceq9}q:(X\x05\x00\x00\x00mode2q;}q<(X\x06\x00\x00\x00paramsq=}q>(X\x08\x00\x00\x00thresoldq?X\x04\x00\x00\x001e-5q@X\n\x00\x00\x00annotationqAX\x16\x00\x00\x00path_to/annotation.bedqBuX\x06\x00\x00\x00scriptqCX\x1e\x00\x00\x00path_to/performance_script1.pyqDuX\x05\x00\x00\x00mode1qE}qF(X\x06\x00\x00\x00paramsqGX\x15\x00\x00\x00some other input hereqHX\x06\x00\x00\x00scriptqIX\x1e\x00\x00\x00path_to/performance_script1.pyqJuuuX\t\x00\x00\x00wildcardsqKcsnakemake.io\nWildcards\nqL)\x81qMX\x0c\x00\x00\x00pranav_test1qNa}qO(h\x08}qPX\x04\x00\x00\x00dataqQK\x00N\x86qRsX\x04\x00\x00\x00dataqShNubX\x03\x00\x00\x00logqTcsnakemake.io\nLog\nqU)\x81qV}qWh\x08}qXsbX\x06\x00\x00\x00outputqYcsnakemake.io\nOutputFiles\nqZ)\x81q[(X\x1c\x00\x00\x00data/pranav_test1_allele1.faq\\X\x1c\x00\x00\x00data/pranav_test1_allele2.faq]e}q^(h\x08}q_(X\x02\x00\x00\x00a1q`K\x00N\x86qaX\x02\x00\x00\x00a2qbK\x01N\x86qcuh`h\\hbh]ubX\t\x00\x00\x00resourcesqdcsnakemake.io\nResources\nqe)\x81qf(K\x01K\x01e}qg(h\x08}qh(X\x06\x00\x00\x00_nodesqiK\x00N\x86qjX\x06\x00\x00\x00_coresqkK\x01N\x86qluhkK\x01hiK\x01ubub.')
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
import os 
os.system('pwd')
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
