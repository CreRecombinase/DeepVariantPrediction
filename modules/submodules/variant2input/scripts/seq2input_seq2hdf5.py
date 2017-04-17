import h5py
import numpy as np
import sys
if 'scripts' not in sys.path:
    sys.path.insert(0, 'scripts')
import my_python

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
