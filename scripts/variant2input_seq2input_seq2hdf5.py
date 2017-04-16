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

huge_array = np.zeros((int(length), args.window, 4), np.bool_)
encode = {'A': 0, 'T': 3, 'G': 1, 'C': 2}
seq2hdf5(input.a1, output.a1, params.window, huge_array, encode)
seq2hdf5(input.a2, output.a2, params.window, huge_array, encode)
