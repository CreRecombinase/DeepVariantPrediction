
######## Snakemake header ########
import sys; sys.path.insert(0, "/project2/xinhe/yanyul/softwares/Anaconda2/envs/deepvarpred/lib/python3.5/site-packages"); import pickle; snakemake = pickle.loads(b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\x06\x00\x00\x00configq\x03}q\x04(X\x05\x00\x00\x00labelq\x05}q\x06(X\x06\x00\x00\x00group1q\x07}q\x08(X\x04\x00\x00\x00E082q\tK\xc1X\x04\x00\x00\x00E081q\nK\xc0uX\x06\x00\x00\x00group2q\x0b}q\x0c(X\x06\x00\x00\x00Noonanq\rK\xc3X\x04\x00\x00\x00E129q\x0eK\xc2uuX\x04\x00\x00\x00dataq\x0f}q\x10(X\x0c\x00\x00\x00pranav_test2q\x11}q\x12(X\x04\x00\x00\x00nameq\x13X\x15\x00\x00\x00test/pranav_test2.txtq\x14X\x06\x00\x00\x00methodq\x15X\x1c\x00\x00\x00_formatting_pranav.snakemakeq\x16uX\x0c\x00\x00\x00pranav_test1q\x17}q\x18(X\x04\x00\x00\x00nameq\x19X\x15\x00\x00\x00test/pranav_test1.txtq\x1aX\x06\x00\x00\x00methodq\x1bX\x1c\x00\x00\x00_formatting_pranav.snakemakeq\x1cuuX\x0f\x00\x00\x00genome_assemblyq\x1dX-\x00\x00\x00/project2/xinhe/yanyul/databases/hg19/hg19.faq\x1eX\x0b\x00\x00\x00window_sizeq\x1fM\xe8\x03X\x05\x00\x00\x00modelq }q!(X\x04\x00\x00\x00nameq"X\x1d\x00\x00\x00path_to_model/model_name.hdf5q#X\t\x00\x00\x00snakemakeq$X\x11\x00\x00\x00path_to/snakemakeq%X\x07\x00\x00\x00workdirq&X\x10\x00\x00\x00path_to/work_dirq\'uX\x0b\x00\x00\x00performanceq(}q)(X\x05\x00\x00\x00mode2q*}q+(X\x06\x00\x00\x00scriptq,X\x1e\x00\x00\x00path_to/performance_script1.pyq-X\x06\x00\x00\x00paramsq.}q/(X\n\x00\x00\x00annotationq0X\x16\x00\x00\x00path_to/annotation.bedq1X\x08\x00\x00\x00thresoldq2X\x04\x00\x00\x001e-5q3uuX\x05\x00\x00\x00mode1q4}q5(X\x06\x00\x00\x00scriptq6X\x1e\x00\x00\x00path_to/performance_script1.pyq7X\x06\x00\x00\x00paramsq8X\x15\x00\x00\x00some other input hereq9uuuX\x04\x00\x00\x00ruleq:X\n\x00\x00\x00expand_bedq;X\t\x00\x00\x00resourcesq<csnakemake.io\nResources\nq=)\x81q>(K\x01K\x01e}q?(X\x06\x00\x00\x00_coresq@K\x01X\x06\x00\x00\x00_nodesqAK\x01X\x06\x00\x00\x00_namesqB}qC(hAK\x00N\x86qDh@K\x01N\x86qEuubX\t\x00\x00\x00wildcardsqFcsnakemake.io\nWildcards\nqG)\x81qHX\x0c\x00\x00\x00pranav_test1qIa}qJ(X\x04\x00\x00\x00dataqKhIhB}qLX\x04\x00\x00\x00dataqMK\x00N\x86qNsubX\x07\x00\x00\x00threadsqOK\x01X\x05\x00\x00\x00inputqPcsnakemake.io\nInputFiles\nqQ)\x81qRX\x1b\x00\x00\x00data/pranav_test1.formattedqSa}qThB}qUsbX\x03\x00\x00\x00logqVcsnakemake.io\nLog\nqW)\x81qX}qYhB}qZsbX\x06\x00\x00\x00outputq[csnakemake.io\nOutputFiles\nq\\)\x81q]X$\x00\x00\x00data/pranav_test1.formatted.expendedq^a}q_hB}q`sbX\x06\x00\x00\x00paramsqacsnakemake.io\nParams\nqb)\x81qcM\xe8\x03a}qd(X\x06\x00\x00\x00windowqeM\xe8\x03hB}qfheK\x00N\x86qgsubub.')
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
