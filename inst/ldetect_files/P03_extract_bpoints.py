#!/usr/bin/env python3.4

import matplotlib.pyplot as pt
import csv
import numpy as np
import pickle
import itertools

import commanderline.commander_line as cl
import ldetect.baselib.flat_file_consts as cnst
import ldetect.baselib.flat_file as flat

def chr_bpoints_to_bed(name, dataset_path, subset, input_pickle_fname):
	'''
	subset is one of ['fourier', 'fourier_ls', 'uniform', 'uniform_ls']
	'''
	
	# input_config = cnst.const['orig_data_'+dataset]
	input_config = cnst.return_conf(dataset_path)

	partitions = flat.read_partitions(name, input_config)

	with open(input_pickle_fname, 'rb') as f_in:
		loaded = pickle.load(f_in)

		# print(loaded)

		loci = loaded[subset]['loci']

		first = partitions[0][0]
		last = partitions[len(partitions)-1][1]

		# print(loci)

		print('chr','\t','start','\t','stop')

		print(name,'\t',first,'\t',loci[0])

		for i in range(0,len(loci)-1):
			print(name,'\t',loci[i],'\t',loci[i+1])

		print(name,'\t',loci[len(loci)-1],'\t',last+1)


if __name__ == '__main__':
    cl.commander_line((chr_bpoints_to_bed, ), print_argv_to_output=False, print_done=False)