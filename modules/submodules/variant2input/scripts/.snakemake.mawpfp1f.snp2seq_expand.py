
######## Snakemake header ########
import sys; sys.path.insert(0, "/project2/xinhe/yanyul/softwares/Anaconda2/envs/deepvarpred/lib/python3.5/site-packages"); import pickle; snakemake = pickle.loads(b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\t\x00\x00\x00resourcesq\x03csnakemake.io\nResources\nq\x04)\x81q\x05(K\x01K\x01e}q\x06(X\x06\x00\x00\x00_namesq\x07}q\x08(X\x06\x00\x00\x00_coresq\tK\x00N\x86q\nX\x06\x00\x00\x00_nodesq\x0bK\x01N\x86q\x0cuh\tK\x01h\x0bK\x01ubX\x04\x00\x00\x00ruleq\rX\n\x00\x00\x00expand_bedq\x0eX\x05\x00\x00\x00inputq\x0fcsnakemake.io\nInputFiles\nq\x10)\x81q\x11X\x1b\x00\x00\x00data/pranav_test1.formattedq\x12a}q\x13h\x07}q\x14sbX\t\x00\x00\x00wildcardsq\x15csnakemake.io\nWildcards\nq\x16)\x81q\x17X\x0c\x00\x00\x00pranav_test1q\x18a}q\x19(h\x07}q\x1aX\x04\x00\x00\x00dataq\x1bK\x00N\x86q\x1csX\x04\x00\x00\x00dataq\x1dh\x18ubX\x06\x00\x00\x00paramsq\x1ecsnakemake.io\nParams\nq\x1f)\x81q M\xe8\x03a}q!(h\x07}q"X\x06\x00\x00\x00windowq#K\x00N\x86q$sh#M\xe8\x03ubX\x07\x00\x00\x00threadsq%K\x01X\x03\x00\x00\x00logq&csnakemake.io\nLog\nq\')\x81q(}q)h\x07}q*sbX\x06\x00\x00\x00configq+}q,(X\x0f\x00\x00\x00genome_assemblyq-X-\x00\x00\x00/project2/xinhe/yanyul/databases/hg19/hg19.faq.X\x04\x00\x00\x00dataq/}q0(X\x0c\x00\x00\x00pranav_test2q1}q2(X\x06\x00\x00\x00methodq3X\x1c\x00\x00\x00_formatting_pranav.snakemakeq4X\x04\x00\x00\x00nameq5X\x15\x00\x00\x00test/pranav_test2.txtq6uX\x0c\x00\x00\x00pranav_test1q7}q8(X\x06\x00\x00\x00methodq9X\x1c\x00\x00\x00_formatting_pranav.snakemakeq:X\x04\x00\x00\x00nameq;X\x15\x00\x00\x00test/pranav_test1.txtq<uuX\x0b\x00\x00\x00window_sizeq=M\xe8\x03X\x0b\x00\x00\x00performanceq>}q?(X\x05\x00\x00\x00mode2q@}qA(X\x06\x00\x00\x00scriptqBX\x1e\x00\x00\x00path_to/performance_script1.pyqCX\x06\x00\x00\x00paramsqD}qE(X\x08\x00\x00\x00thresoldqFX\x04\x00\x00\x001e-5qGX\n\x00\x00\x00annotationqHX\x16\x00\x00\x00path_to/annotation.bedqIuuX\x05\x00\x00\x00mode1qJ}qK(X\x06\x00\x00\x00scriptqLX\x1e\x00\x00\x00path_to/performance_script1.pyqMX\x06\x00\x00\x00paramsqNX\x15\x00\x00\x00some other input hereqOuuX\x05\x00\x00\x00modelqP}qQ(X\x07\x00\x00\x00workdirqRX\x10\x00\x00\x00path_to/work_dirqSX\t\x00\x00\x00snakemakeqTX\x11\x00\x00\x00path_to/snakemakeqUX\x04\x00\x00\x00nameqVX\x1d\x00\x00\x00path_to_model/model_name.hdf5qWuX\x05\x00\x00\x00labelqX}qY(X\x06\x00\x00\x00group2qZ}q[(X\x06\x00\x00\x00Noonanq\\K\xc3X\x04\x00\x00\x00E129q]K\xc2uX\x06\x00\x00\x00group1q^}q_(X\x04\x00\x00\x00E082q`K\xc1X\x04\x00\x00\x00E081qaK\xc0uuuX\x06\x00\x00\x00outputqbcsnakemake.io\nOutputFiles\nqc)\x81qdX$\x00\x00\x00data/pranav_test1.formatted.expendedqea}qfh\x07}qgsbub.')
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
cmd = '''cat {in} | awk -F"\\t" '{{ printf $1"\\t"$2-{start}"\\t"$3+{end}; \
 {{ for (i=4;i<=NF;i++) printf "\\t%s", $i }}; print ""}}\' > {out}'''.format(in=input, out=output, start=start, end=end)
my_python.myOsSystem(cmd, False)
