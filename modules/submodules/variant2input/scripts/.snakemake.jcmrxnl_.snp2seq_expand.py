
######## Snakemake header ########
import sys; sys.path.insert(0, "/project2/xinhe/yanyul/softwares/Anaconda2/envs/deepvarpred/lib/python3.5/site-packages"); import pickle; snakemake = pickle.loads(b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\x03\x00\x00\x00logq\x03csnakemake.io\nLog\nq\x04)\x81q\x05}q\x06X\x06\x00\x00\x00_namesq\x07}q\x08sbX\t\x00\x00\x00resourcesq\tcsnakemake.io\nResources\nq\n)\x81q\x0b(K\x01K\x01e}q\x0c(X\x06\x00\x00\x00_coresq\rK\x01X\x06\x00\x00\x00_nodesq\x0eK\x01h\x07}q\x0f(h\rK\x00N\x86q\x10h\x0eK\x01N\x86q\x11uubX\x05\x00\x00\x00inputq\x12csnakemake.io\nInputFiles\nq\x13)\x81q\x14X\x1b\x00\x00\x00data/pranav_test1.formattedq\x15a}q\x16h\x07}q\x17sbX\x06\x00\x00\x00paramsq\x18csnakemake.io\nParams\nq\x19)\x81q\x1aM\xe8\x03a}q\x1b(X\x06\x00\x00\x00windowq\x1cM\xe8\x03h\x07}q\x1dh\x1cK\x00N\x86q\x1esubX\x07\x00\x00\x00threadsq\x1fK\x01X\x04\x00\x00\x00ruleq X\n\x00\x00\x00expand_bedq!X\x06\x00\x00\x00outputq"csnakemake.io\nOutputFiles\nq#)\x81q$X$\x00\x00\x00data/pranav_test1.formatted.expendedq%a}q&h\x07}q\'sbX\t\x00\x00\x00wildcardsq(csnakemake.io\nWildcards\nq))\x81q*X\x0c\x00\x00\x00pranav_test1q+a}q,(X\x04\x00\x00\x00dataq-h+h\x07}q.X\x04\x00\x00\x00dataq/K\x00N\x86q0subX\x06\x00\x00\x00configq1}q2(X\x05\x00\x00\x00labelq3}q4(X\x06\x00\x00\x00group1q5}q6(X\x04\x00\x00\x00E081q7K\xc0X\x04\x00\x00\x00E082q8K\xc1uX\x06\x00\x00\x00group2q9}q:(X\x06\x00\x00\x00Noonanq;K\xc3X\x04\x00\x00\x00E129q<K\xc2uuX\x0f\x00\x00\x00genome_assemblyq=X-\x00\x00\x00/project2/xinhe/yanyul/databases/hg19/hg19.faq>X\x05\x00\x00\x00modelq?}q@(X\x07\x00\x00\x00workdirqAX\x10\x00\x00\x00path_to/work_dirqBX\x04\x00\x00\x00nameqCX\x1d\x00\x00\x00path_to_model/model_name.hdf5qDX\t\x00\x00\x00snakemakeqEX\x11\x00\x00\x00path_to/snakemakeqFuX\x04\x00\x00\x00dataqG}qH(X\x0c\x00\x00\x00pranav_test1qI}qJ(X\x06\x00\x00\x00methodqKX\x1c\x00\x00\x00_formatting_pranav.snakemakeqLX\x04\x00\x00\x00nameqMX\x15\x00\x00\x00test/pranav_test1.txtqNuX\x0c\x00\x00\x00pranav_test2qO}qP(X\x06\x00\x00\x00methodqQX\x1c\x00\x00\x00_formatting_pranav.snakemakeqRX\x04\x00\x00\x00nameqSX\x15\x00\x00\x00test/pranav_test2.txtqTuuX\x0b\x00\x00\x00performanceqU}qV(X\x05\x00\x00\x00mode1qW}qX(X\x06\x00\x00\x00scriptqYX\x1e\x00\x00\x00path_to/performance_script1.pyqZX\x06\x00\x00\x00paramsq[X\x15\x00\x00\x00some other input hereq\\uX\x05\x00\x00\x00mode2q]}q^(X\x06\x00\x00\x00scriptq_X\x1e\x00\x00\x00path_to/performance_script1.pyq`X\x06\x00\x00\x00paramsqa}qb(X\n\x00\x00\x00annotationqcX\x16\x00\x00\x00path_to/annotation.bedqdX\x08\x00\x00\x00thresoldqeX\x04\x00\x00\x001e-5qfuuuX\x0b\x00\x00\x00window_sizeqgM\xe8\x03uub.')
######## Original script #########
import my_python
if params.window % 2 == 0:
	even_flag = 0
	start = params.window / 2 - 1
	end = params.window / 2
else:
	start = (params.window - 1) / 2
	end = -start
	even_flag = 1
midpos = start
cmd = '''cat {input} | awk -F"\\t" '{{ printf $1"\\t"$2-{start}"\\t"$3+{end}; \
 {{ for (i=4;i<=NF;i++) printf "\\t%s", $i }}; print ""}}\' > {out}'''.format(input=input, out=output, start=start, end=end)
my_python.myOsSystem(cmd, False)
