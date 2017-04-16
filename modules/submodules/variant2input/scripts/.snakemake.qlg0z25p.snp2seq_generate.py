
######## Snakemake header ########
import sys; sys.path.insert(0, "/project2/xinhe/yanyul/softwares/Anaconda2/envs/deepvarpred/lib/python3.5/site-packages"); import pickle; snakemake = pickle.loads(b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\x04\x00\x00\x00ruleq\x03X\x18\x00\x00\x00generate_allele1_allele2q\x04X\x06\x00\x00\x00outputq\x05csnakemake.io\nOutputFiles\nq\x06)\x81q\x07(X\x1c\x00\x00\x00data/pranav_test1_allele1.faq\x08X\x1c\x00\x00\x00data/pranav_test1_allele2.faq\te}q\n(X\x06\x00\x00\x00_namesq\x0b}q\x0c(X\x02\x00\x00\x00a1q\rK\x00N\x86q\x0eX\x02\x00\x00\x00a2q\x0fK\x01N\x86q\x10uh\x0fh\th\rh\x08ubX\x03\x00\x00\x00logq\x11csnakemake.io\nLog\nq\x12)\x81q\x13}q\x14h\x0b}q\x15sbX\x07\x00\x00\x00threadsq\x16K\x01X\t\x00\x00\x00resourcesq\x17csnakemake.io\nResources\nq\x18)\x81q\x19(K\x01K\x01e}q\x1a(h\x0b}q\x1b(X\x06\x00\x00\x00_nodesq\x1cK\x01N\x86q\x1dX\x06\x00\x00\x00_coresq\x1eK\x00N\x86q\x1fuh\x1eK\x01h\x1cK\x01ubX\x06\x00\x00\x00paramsq csnakemake.io\nParams\nq!)\x81q"M\xe8\x03a}q#(h\x0b}q$X\x06\x00\x00\x00windowq%K\x00N\x86q&sh%M\xe8\x03ubX\x05\x00\x00\x00inputq\'csnakemake.io\nInputFiles\nq()\x81q)X\'\x00\x00\x00data/pranav_test1.formatted.expended.faq*a}q+h\x0b}q,sbX\t\x00\x00\x00wildcardsq-csnakemake.io\nWildcards\nq.)\x81q/X\x0c\x00\x00\x00pranav_test1q0a}q1(h\x0b}q2X\x04\x00\x00\x00dataq3K\x00N\x86q4sX\x04\x00\x00\x00dataq5h0ubX\x06\x00\x00\x00configq6}q7(X\x0b\x00\x00\x00performanceq8}q9(X\x05\x00\x00\x00mode1q:}q;(X\x06\x00\x00\x00scriptq<X\x1e\x00\x00\x00path_to/performance_script1.pyq=X\x06\x00\x00\x00paramsq>X\x15\x00\x00\x00some other input hereq?uX\x05\x00\x00\x00mode2q@}qA(X\x06\x00\x00\x00scriptqBX\x1e\x00\x00\x00path_to/performance_script1.pyqCX\x06\x00\x00\x00paramsqD}qE(X\n\x00\x00\x00annotationqFX\x16\x00\x00\x00path_to/annotation.bedqGX\x08\x00\x00\x00thresoldqHX\x04\x00\x00\x001e-5qIuuuX\x04\x00\x00\x00dataqJ}qK(X\x0c\x00\x00\x00pranav_test1qL}qM(X\x06\x00\x00\x00methodqNX\x1c\x00\x00\x00_formatting_pranav.snakemakeqOX\x04\x00\x00\x00nameqPX\x15\x00\x00\x00test/pranav_test1.txtqQuX\x0c\x00\x00\x00pranav_test2qR}qS(X\x06\x00\x00\x00methodqTX\x1c\x00\x00\x00_formatting_pranav.snakemakeqUX\x04\x00\x00\x00nameqVX\x15\x00\x00\x00test/pranav_test2.txtqWuuX\x0f\x00\x00\x00genome_assemblyqXX-\x00\x00\x00/project2/xinhe/yanyul/databases/hg19/hg19.faqYX\x05\x00\x00\x00labelqZ}q[(X\x06\x00\x00\x00group1q\\}q](X\x04\x00\x00\x00E082q^K\xc1X\x04\x00\x00\x00E081q_K\xc0uX\x06\x00\x00\x00group2q`}qa(X\x04\x00\x00\x00E129qbK\xc2X\x06\x00\x00\x00NoonanqcK\xc3uuX\x05\x00\x00\x00modelqd}qe(X\x07\x00\x00\x00workdirqfX\x10\x00\x00\x00path_to/work_dirqgX\t\x00\x00\x00snakemakeqhX\x11\x00\x00\x00path_to/snakemakeqiX\x04\x00\x00\x00nameqjX\x1d\x00\x00\x00path_to_model/model_name.hdf5qkuX\x0b\x00\x00\x00window_sizeqlM\xe8\x03uub.')
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
os.system('pwd > out')
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
