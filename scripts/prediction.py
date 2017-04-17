import argparse
parser = argparse.ArgumentParser(prog='prediction.py', description='''
	Given HDF5 input (dataset named x) and model, make the prediction
	''')
parser.add_argument('--model', help='''
	HDF5 which can be loaded by keras.models.load_model and it takes
    raw data as input
	''')
parser.add_argument('--data', help='''
	HDF5 file where input data is saved in x dataset and the dimension
    is (nsamples, window_size, 4)
	''')
parser.add_argument('--out')
args = parser.parse_args()

from keras.models import load_model
import h5py
import numpy as np
import sys
if 'scripts' not in sys.path:
    sys.path.insert(0, 'scripts')
import my_python

print('Loading model')
model = load_model(args.model)

print('Loading test data')
testmat = h5py.File(args.data,'r')
x = testmat['x'][()]
testmat.close()

outfile = args.out

print('Predicting on test sequences')
y = model.predict(x, verbose=2)
ny = int(y.shape[0] / 2)
y = (y[:ny] + y[ny:]) / 2
out = h5py.File(outfile, 'w')
out.create_dataset('y_pred',data=y)
out.close()
