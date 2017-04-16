
######## Snakemake header ########
import sys; sys.path.insert(0, "/project2/xinhe/yanyul/softwares/Anaconda2/envs/deepvarpred/lib/python3.5/site-packages"); import pickle; snakemake = pickle.loads(b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\x06\x00\x00\x00paramsq\x03csnakemake.io\nParams\nq\x04)\x81q\x05M\xe8\x03a}q\x06(X\x06\x00\x00\x00windowq\x07M\xe8\x03X\x06\x00\x00\x00_namesq\x08}q\th\x07K\x00N\x86q\nsubX\x05\x00\x00\x00inputq\x0bcsnakemake.io\nInputFiles\nq\x0c)\x81q\r(X\x1c\x00\x00\x00data/pranav_test1_allele2.faq\x0eX\x1c\x00\x00\x00data/pranav_test1_allele1.faq\x0fe}q\x10(X\x02\x00\x00\x00a2q\x11h\x0eX\x02\x00\x00\x00a1q\x12h\x0fh\x08}q\x13(h\x11K\x00N\x86q\x14h\x12K\x01N\x86q\x15uubX\x06\x00\x00\x00configq\x16}q\x17(X\x0f\x00\x00\x00genome_assemblyq\x18X-\x00\x00\x00/project2/xinhe/yanyul/databases/hg19/hg19.faq\x19X\x04\x00\x00\x00dataq\x1a}q\x1b(X\x0c\x00\x00\x00pranav_test2q\x1c}q\x1d(X\x06\x00\x00\x00methodq\x1eX\x1c\x00\x00\x00_formatting_pranav.snakemakeq\x1fX\x04\x00\x00\x00nameq X\x15\x00\x00\x00test/pranav_test2.txtq!uX\x0c\x00\x00\x00pranav_test1q"}q#(X\x06\x00\x00\x00methodq$X\x1c\x00\x00\x00_formatting_pranav.snakemakeq%X\x04\x00\x00\x00nameq&X\x15\x00\x00\x00test/pranav_test1.txtq\'uuX\x05\x00\x00\x00modelq(}q)(X\t\x00\x00\x00snakemakeq*X\x11\x00\x00\x00path_to/snakemakeq+X\x04\x00\x00\x00nameq,X\x1d\x00\x00\x00path_to_model/model_name.hdf5q-X\x07\x00\x00\x00workdirq.X\x10\x00\x00\x00path_to/work_dirq/uX\x05\x00\x00\x00labelq0}q1(X\x06\x00\x00\x00group2q2}q3(X\x04\x00\x00\x00E129q4K\xc2X\x06\x00\x00\x00Noonanq5K\xc3uX\x06\x00\x00\x00group1q6}q7(X\x04\x00\x00\x00E081q8K\xc0X\x04\x00\x00\x00E082q9K\xc1uuX\x0b\x00\x00\x00performanceq:}q;(X\x05\x00\x00\x00mode1q<}q=(X\x06\x00\x00\x00paramsq>X\x15\x00\x00\x00some other input hereq?X\x06\x00\x00\x00scriptq@X\x1e\x00\x00\x00path_to/performance_script1.pyqAuX\x05\x00\x00\x00mode2qB}qC(X\x06\x00\x00\x00paramsqD}qE(X\n\x00\x00\x00annotationqFX\x16\x00\x00\x00path_to/annotation.bedqGX\x08\x00\x00\x00thresoldqHX\x04\x00\x00\x001e-5qIuX\x06\x00\x00\x00scriptqJX\x1e\x00\x00\x00path_to/performance_script1.pyqKuuX\x0b\x00\x00\x00window_sizeqLM\xe8\x03uX\t\x00\x00\x00wildcardsqMcsnakemake.io\nWildcards\nqN)\x81qOX\x0c\x00\x00\x00pranav_test1qPa}qQ(X\x04\x00\x00\x00dataqRhPh\x08}qSX\x04\x00\x00\x00dataqTK\x00N\x86qUsubX\x07\x00\x00\x00threadsqVK\x01X\x04\x00\x00\x00ruleqWX\x08\x00\x00\x00seq2hdf5qXX\x06\x00\x00\x00outputqYcsnakemake.io\nOutputFiles\nqZ)\x81q[(X\x1f\x00\x00\x00input/pranav_test1_allele2.hdf5q\\X\x1f\x00\x00\x00input/pranav_test1_allele1.hdf5q]e}q^(h\x11h\\h\x12h]h\x08}q_(h\x11K\x00N\x86q`h\x12K\x01N\x86qauubX\x03\x00\x00\x00logqbcsnakemake.io\nLog\nqc)\x81qd}qeh\x08}qfsbX\t\x00\x00\x00resourcesqgcsnakemake.io\nResources\nqh)\x81qi(K\x01K\x01e}qj(X\x06\x00\x00\x00_nodesqkK\x01X\x06\x00\x00\x00_coresqlK\x01h\x08}qm(hlK\x00N\x86qnhkK\x01N\x86qouubub.')
######## Original script #########
import h5py
import numpy as np

def seq2hdf5(fasta, out, window, huge_array, encode):
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

huge_array = np.zeros((int(length), snakemake.params.window, 4), np.bool_)
encode = {'A': 0, 'T': 3, 'G': 1, 'C': 2}
seq2hdf5(snakemake.input.a1, snakemake.output.a1, snakemake.params.window, huge_array, encode)
seq2hdf5(snakemake.input.a2, snakemake.output.a2, snakemake.params.window, huge_array, encode)
