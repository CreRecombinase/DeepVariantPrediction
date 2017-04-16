
######## Snakemake header ########
import sys; sys.path.insert(0, "/project2/xinhe/yanyul/softwares/Anaconda2/envs/deepvarpred/lib/python3.5/site-packages"); import pickle; snakemake = pickle.loads(b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\x04\x00\x00\x00ruleq\x03X\x08\x00\x00\x00seq2hdf5q\x04X\x07\x00\x00\x00threadsq\x05K\x01X\t\x00\x00\x00resourcesq\x06csnakemake.io\nResources\nq\x07)\x81q\x08(K\x01K\x01e}q\t(X\x06\x00\x00\x00_nodesq\nK\x01X\x06\x00\x00\x00_coresq\x0bK\x01X\x06\x00\x00\x00_namesq\x0c}q\r(h\nK\x01N\x86q\x0eh\x0bK\x00N\x86q\x0fuubX\x06\x00\x00\x00paramsq\x10csnakemake.io\nParams\nq\x11)\x81q\x12M\xe8\x03a}q\x13(X\x06\x00\x00\x00windowq\x14M\xe8\x03h\x0c}q\x15h\x14K\x00N\x86q\x16subX\x06\x00\x00\x00configq\x17}q\x18(X\x04\x00\x00\x00dataq\x19}q\x1a(X\x0c\x00\x00\x00pranav_test1q\x1b}q\x1c(X\x04\x00\x00\x00nameq\x1dX\x15\x00\x00\x00test/pranav_test1.txtq\x1eX\x06\x00\x00\x00methodq\x1fX\x1c\x00\x00\x00_formatting_pranav.snakemakeq uX\x0c\x00\x00\x00pranav_test2q!}q"(X\x04\x00\x00\x00nameq#X\x15\x00\x00\x00test/pranav_test2.txtq$X\x06\x00\x00\x00methodq%X\x1c\x00\x00\x00_formatting_pranav.snakemakeq&uuX\x05\x00\x00\x00labelq\'}q((X\x06\x00\x00\x00group1q)}q*(X\x04\x00\x00\x00E081q+K\xc0X\x04\x00\x00\x00E082q,K\xc1uX\x06\x00\x00\x00group2q-}q.(X\x04\x00\x00\x00E129q/K\xc2X\x06\x00\x00\x00Noonanq0K\xc3uuX\x05\x00\x00\x00modelq1}q2(X\x07\x00\x00\x00workdirq3X\x10\x00\x00\x00path_to/work_dirq4X\x04\x00\x00\x00nameq5X\x1d\x00\x00\x00path_to_model/model_name.hdf5q6X\t\x00\x00\x00snakemakeq7X\x11\x00\x00\x00path_to/snakemakeq8uX\x0f\x00\x00\x00genome_assemblyq9X-\x00\x00\x00/project2/xinhe/yanyul/databases/hg19/hg19.faq:X\x0b\x00\x00\x00window_sizeq;M\xe8\x03X\x0b\x00\x00\x00performanceq<}q=(X\x05\x00\x00\x00mode1q>}q?(X\x06\x00\x00\x00paramsq@X\x15\x00\x00\x00some other input hereqAX\x06\x00\x00\x00scriptqBX\x1e\x00\x00\x00path_to/performance_script1.pyqCuX\x05\x00\x00\x00mode2qD}qE(X\x06\x00\x00\x00paramsqF}qG(X\x08\x00\x00\x00thresoldqHX\x04\x00\x00\x001e-5qIX\n\x00\x00\x00annotationqJX\x16\x00\x00\x00path_to/annotation.bedqKuX\x06\x00\x00\x00scriptqLX\x1e\x00\x00\x00path_to/performance_script1.pyqMuuuX\x03\x00\x00\x00logqNcsnakemake.io\nLog\nqO)\x81qP}qQh\x0c}qRsbX\t\x00\x00\x00wildcardsqScsnakemake.io\nWildcards\nqT)\x81qUX\x0c\x00\x00\x00pranav_test1qVa}qW(X\x04\x00\x00\x00dataqXhVh\x0c}qYX\x04\x00\x00\x00dataqZK\x00N\x86q[subX\x06\x00\x00\x00outputq\\csnakemake.io\nOutputFiles\nq])\x81q^(X\x1f\x00\x00\x00input/pranav_test1_allele1.hdf5q_X\x1f\x00\x00\x00input/pranav_test1_allele2.hdf5q`e}qa(h\x0c}qb(X\x02\x00\x00\x00a2qcK\x01N\x86qdX\x02\x00\x00\x00a1qeK\x00N\x86qfuhch`heh_ubX\x05\x00\x00\x00inputqgcsnakemake.io\nInputFiles\nqh)\x81qi(X\x1c\x00\x00\x00data/pranav_test1_allele1.faqjX\x1c\x00\x00\x00data/pranav_test1_allele2.faqke}ql(hchkhehjh\x0c}qm(heK\x00N\x86qnhcK\x01N\x86qouubub.')
######## Original script #########
import h5py
import numpy as np

def seq2hdf5(fasta, out, window, encode):
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
seq2hdf5(snakemake.input.a1, snakemake.output.a1, snakemake.params.window, encode)
seq2hdf5(snakemake.input.a2, snakemake.output.a2, snakemake.params.window, encode)
