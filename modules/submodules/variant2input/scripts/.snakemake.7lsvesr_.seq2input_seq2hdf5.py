
######## Snakemake header ########
import sys; sys.path.insert(0, "/project2/xinhe/yanyul/softwares/Anaconda2/envs/deepvarpred/lib/python3.5/site-packages"); import pickle; snakemake = pickle.loads(b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\x04\x00\x00\x00ruleq\x03X\x08\x00\x00\x00seq2hdf5q\x04X\x06\x00\x00\x00paramsq\x05csnakemake.io\nParams\nq\x06)\x81q\x07M\xe8\x03a}q\x08(X\x06\x00\x00\x00_namesq\t}q\nX\x06\x00\x00\x00windowq\x0bK\x00N\x86q\x0csh\x0bM\xe8\x03ubX\x05\x00\x00\x00inputq\rcsnakemake.io\nInputFiles\nq\x0e)\x81q\x0f(X\x1c\x00\x00\x00data/pranav_test1_allele2.faq\x10X\x1c\x00\x00\x00data/pranav_test1_allele1.faq\x11e}q\x12(h\t}q\x13(X\x02\x00\x00\x00a2q\x14K\x00N\x86q\x15X\x02\x00\x00\x00a1q\x16K\x01N\x86q\x17uh\x14h\x10h\x16h\x11ubX\x03\x00\x00\x00logq\x18csnakemake.io\nLog\nq\x19)\x81q\x1a}q\x1bh\t}q\x1csbX\t\x00\x00\x00wildcardsq\x1dcsnakemake.io\nWildcards\nq\x1e)\x81q\x1fX\x0c\x00\x00\x00pranav_test1q a}q!(h\t}q"X\x04\x00\x00\x00dataq#K\x00N\x86q$sX\x04\x00\x00\x00dataq%h ubX\t\x00\x00\x00resourcesq&csnakemake.io\nResources\nq\')\x81q((K\x01K\x01e}q)(h\t}q*(X\x06\x00\x00\x00_nodesq+K\x00N\x86q,X\x06\x00\x00\x00_coresq-K\x01N\x86q.uh+K\x01h-K\x01ubX\x06\x00\x00\x00configq/}q0(X\x05\x00\x00\x00modelq1}q2(X\t\x00\x00\x00snakemakeq3X\x11\x00\x00\x00path_to/snakemakeq4X\x04\x00\x00\x00nameq5X\x1d\x00\x00\x00path_to_model/model_name.hdf5q6X\x07\x00\x00\x00workdirq7X\x10\x00\x00\x00path_to/work_dirq8uX\x04\x00\x00\x00dataq9}q:(X\x0c\x00\x00\x00pranav_test1q;}q<(X\x04\x00\x00\x00nameq=X\x15\x00\x00\x00test/pranav_test1.txtq>X\x06\x00\x00\x00methodq?X\x1c\x00\x00\x00_formatting_pranav.snakemakeq@uX\x0c\x00\x00\x00pranav_test2qA}qB(X\x04\x00\x00\x00nameqCX\x15\x00\x00\x00test/pranav_test2.txtqDX\x06\x00\x00\x00methodqEX\x1c\x00\x00\x00_formatting_pranav.snakemakeqFuuX\x0b\x00\x00\x00window_sizeqGM\xe8\x03X\x0f\x00\x00\x00genome_assemblyqHX-\x00\x00\x00/project2/xinhe/yanyul/databases/hg19/hg19.faqIX\x0b\x00\x00\x00performanceqJ}qK(X\x05\x00\x00\x00mode2qL}qM(X\x06\x00\x00\x00scriptqNX\x1e\x00\x00\x00path_to/performance_script1.pyqOX\x06\x00\x00\x00paramsqP}qQ(X\n\x00\x00\x00annotationqRX\x16\x00\x00\x00path_to/annotation.bedqSX\x08\x00\x00\x00thresoldqTX\x04\x00\x00\x001e-5qUuuX\x05\x00\x00\x00mode1qV}qW(X\x06\x00\x00\x00scriptqXX\x1e\x00\x00\x00path_to/performance_script1.pyqYX\x06\x00\x00\x00paramsqZX\x15\x00\x00\x00some other input hereq[uuX\x05\x00\x00\x00labelq\\}q](X\x06\x00\x00\x00group1q^}q_(X\x04\x00\x00\x00E081q`K\xc0X\x04\x00\x00\x00E082qaK\xc1uX\x06\x00\x00\x00group2qb}qc(X\x06\x00\x00\x00NoonanqdK\xc3X\x04\x00\x00\x00E129qeK\xc2uuuX\x06\x00\x00\x00outputqfcsnakemake.io\nOutputFiles\nqg)\x81qh(X\x1f\x00\x00\x00input/pranav_test1_allele2.hdf5qiX\x1f\x00\x00\x00input/pranav_test1_allele1.hdf5qje}qk(h\x14hih\t}ql(h\x14K\x00N\x86qmh\x16K\x01N\x86qnuh\x16hjubX\x07\x00\x00\x00threadsqoK\x01ub.')
######## Original script #########
import h5py
import numpy as np

def seq2hdf5(fasta, out, window, huge_array, encode):
    cmd = '''cat {file} | wc -l'''.format(file=fasta)
    length = my_python.mySubprocess(cmd, False)
    huge_array = np.zeros((int(length), window, 4), np.bool_)
    counter = 0
    with open(fasta) as infile:
    	for line in infile:
    		seq = line.strip().split('\t')[1].upper()
    		digit_seq = np.zeros((window, 4))
    		for i in range(len(seq)):
    			if seq[i] not in encode:
    				continue
    			else:
    				digit_seq[i, encode[seq[i]]] = 1
    		huge_array[counter] = digit_seq
    		counter += 1

    huge_array_flip = huge_array[:,::-1,::-1];
    huge_array = np.concatenate([huge_array, huge_array_flip],axis=0)
    huge_array = huge_array.astype(np.uint8)


    f = h5py.File(out, 'w')
    f.create_dataset('x', data=huge_array)
    f.close()

encode = {'A': 0, 'T': 3, 'G': 1, 'C': 2}
seq2hdf5(snakemake.input.a1, snakemake.output.a1, snakemake.params.window, huge_array, encode)
seq2hdf5(snakemake.input.a2, snakemake.output.a2, snakemake.params.window, huge_array, encode)
