#!/usr/bin/env python3

'''
Created on Aug 7, 2014

@author: tberisa
'''

import sys
import getopt
import decimal
import pickle
import datetime

import ldetect.baselib.flat_file_consts as cnst
import ldetect.baselib.flat_file as flat
import ldetect.baselib.read_data as rd
import ldetect.baselib.filters as filt

import ldetect.pipeline_elements.E03_matrix_to_vector as matrix_to_vector

import commanderline.commander_line as cl

def pipeline_lean(dataset_path, name, out_fname, begin=-1, end=-1, img='no', orient='diag', red='sum', dataset_name='NONAME'):
    '''
    pipeline_lean(dataset_path, name, begin=-1, end=-1, img='no', orient='diag', red='sum')
    '''
    
    # analysis = matrix_to_vector.MatrixAnalysis(name, cnst.const[dataset], begin, end)
    analysis = matrix_to_vector.MatrixAnalysis(name, cnst.return_conf(dataset_path), begin, end)

    print(analysis.snp_first)
    print(analysis.snp_last)

    t = datetime.datetime.now() 
    t_formatted = t.strftime('%Y_%m_%d_%H_%M_%S')

    # out_fname = 'vector-'+dataset_name+'-'+name+'-'+str(analysis.snp_first)+'-'+str(analysis.snp_last)+'-'+orient+'-'+red+'-img_'+img+'-'+t_formatted
    # out_fname += '.txt.gz'
    flat.print_log_msg('out_fname: '+out_fname)

    if(img=='yes'):
        generate_img = True
    elif(img=='no'):
        generate_img = False
    else:
        raise Exception('Error: Unknown argument: '+img)

    if(orient=='vert'):
        analysis.calc_vert(not generate_img) 
    elif(orient=='diag'):
        analysis.calc_diag_lean(out_fname, cnst.const['out_delim'], not generate_img) 
    else:
        raise Exception('Error: Unknown argument: '+orient)

    if(red=='avg'):
        avg = True
        raise Exception('Average used, but its output is not always consistent - especially for diag!')
    elif(red=='sum'):
        avg = False
    else:
        raise Exception('Error: Unknown argument: '+red)

    # Output is done step-by-step
    # analysis.write_output_to_file(out_fname+'.txt.gz', cnst.const['out_delim'], avg)

    if generate_img:
        analysis.generate_img(out_fname+cnst.const['img_out_ext'])

    flat.print_log_msg('Done')

def print_opt_arg_error():
    print('For help use --help')
    return(2)

# if __name__ == '__main__':
#     main()

if __name__ == '__main__':
    cl.commander_line((pipeline_lean, )) 
