#!/usr/bin/env python3

'''
Created on Aug 7, 2014

@author: tberisa
'''

import sys
import getopt
import decimal
import pickle
import numpy as np
import datetime
import math
import json

import ldetect.baselib.flat_file_consts as cnst
import ldetect.baselib.flat_file as flat
import ldetect.baselib.read_data as rd
import ldetect.baselib.filters as filt
import ldetect.baselib.binary_search as binsrch

import ldetect.pipeline_elements.E05_find_minima as find_minima
import ldetect.pipeline_elements.E07_metric as metric
import ldetect.pipeline_elements.E08_local_search as local_search

import commanderline.commander_line as cl

def pipeline(input_fname, chr_name, dataset_path, n_snps_bw_bpoints, out_fname, begin=-1, end=-1, trackback_delta=200, trackback_step=20, init_search_location=1000):
    config=cnst.return_conf(dataset_path)
    # begin, end = flat.first_last(chr_name, cnst.const[dataset], begin, end)
    begin, end = flat.first_last(chr_name, config, begin, end)
    # READ DATA
    flat.print_log_msg('* Reading data')
    init_array, init_array_x = rd.read_data_raw(input_fname) 

    # Clip the input data to the required range and convert to numpy array
    begin_ind = binsrch.find_ge_ind(init_array_x, begin) # = init_array_x.index(begin)
    end_ind = binsrch.find_le_ind(init_array_x, end) # = init_array_x.index(end)

    np_init_array = np.array(init_array[begin_ind:(end_ind+1)])
    np_init_array_x = np.array(init_array_x[begin_ind:(end_ind+1)])

    # DETERMINE NUMBER OF BREAKPOINTS
    n_bpoints = int(math.ceil( len(np_init_array_x) / n_snps_bw_bpoints - 1 ))
    flat.print_log_msg('* Number of breakpoints: '+repr(n_bpoints))

    # SEARCH FOR FILTER WIDTH
    flat.print_log_msg('* Starting search...')
    found_width = find_minima.custom_binary_search_with_trackback(np_init_array, filt.apply_filter_get_minima, n_bpoints, trackback_delta=trackback_delta, trackback_step=trackback_step, init_search_location=init_search_location)
    flat.print_log_msg('* Found_width: ' + repr(found_width))
    
    # GET MINIMA LOCATIONS
    flat.print_log_msg('* Applying filter and getting minima locations...')
    g = filt.apply_filter(np_init_array, found_width)
    breakpoint_loci = filt.get_minima_loc(g, np_init_array_x)
    
    # METRIC
    flat.print_log_msg('* Calculating metric for non-uniform breakpoints (minima of filtered data)...')
        
    # metric_out = apply_metric(chr_name, begin, end, cnst.const[dataset], breakpoint_loci)
    metric_out = apply_metric(chr_name, begin, end, config, breakpoint_loci)
    flat.print_log_msg('Global metric:')
    print_metric(metric_out)
    
    # METRIC FOR UNIFORM BREAKPOINTS
    flat.print_log_msg('* Calculating metric for uniform breakpoints...')
    # step = int((end-begin)/(len(breakpoint_loci)+1))
    # breakpoint_loci_uniform = [l for l in range(begin+step, end-step+1, step)] 
    step = int(len(init_array_x)/(len(breakpoint_loci)+1))
    breakpoint_loci_uniform = [init_array_x[i] for i in range(step, len(init_array_x)-step+1, step)]

    # metric_out_uniform = apply_metric(chr_name, begin, end, cnst.const[dataset], breakpoint_loci_uniform)
    metric_out_uniform = apply_metric(chr_name, begin, end, config, breakpoint_loci_uniform)
    flat.print_log_msg('Global metric:')
    print_metric(metric_out_uniform)
    
    # LOCAL SEARCH ON FOURIER - missing N runs
    flat.print_log_msg('* Running local search for fourier...')

    # breakpoint_loci_local_search = run_local_search_complete(chr_name, breakpoint_loci, begin, end, cnst.const[dataset], metric_out)
    breakpoint_loci_local_search = run_local_search_complete(chr_name, breakpoint_loci, begin, end, config, metric_out)
    
    # RUN METRIC AGAIN W/ NEW BREAKPOINTS FROM FOURIER LOCAL SEARCH
    flat.print_log_msg('* Calculating metric for new fourier breakpoints...')
    
    # metric_out_local_search = apply_metric(chr_name, begin, end, cnst.const[dataset], breakpoint_loci_local_search['loci'])    
    metric_out_local_search = apply_metric(chr_name, begin, end, config, breakpoint_loci_local_search['loci'])
    flat.print_log_msg('Global metric:')
    print_metric(metric_out_local_search)
    
    # LOCAL SEARCH ON UNIFORM - missing N runs
    flat.print_log_msg('* Running local search for uniform breakpoints...')

    # breakpoint_loci_uniform_local_search = run_local_search_complete(chr_name, breakpoint_loci_uniform, begin, end, cnst.const[dataset], metric_out_uniform)
    breakpoint_loci_uniform_local_search = run_local_search_complete(chr_name, breakpoint_loci_uniform, begin, end, config, metric_out_uniform)
    
    # RUN METRIC AGAIN W/ NEW BREAKPOINTS FROM UNIFORM
    flat.print_log_msg('* Calculating metric for new uniform breakpoints...')
    
    # metric_out_uniform_local_search = apply_metric(chr_name, begin, end, cnst.const[dataset], breakpoint_loci_uniform_local_search['loci'])    
    metric_out_uniform_local_search = apply_metric(chr_name, begin, end, config, breakpoint_loci_uniform_local_search['loci'])
    flat.print_log_msg('Global metric:')
    print_metric(metric_out_uniform_local_search)
    
    # DUMP DATA INTO PICKLE SO IT CAN BE ANALYZED AND LOOKED AT WITHOUT RE-RUNNING EVERYTHING
    pickle_out = {}
    pickle_out['argv'] = sys.argv
    pickle_out['n_bpoints'] = n_bpoints
    pickle_out['found_width'] = found_width
    pickle_out['fourier'] = {}
    pickle_out['fourier']['loci'] = breakpoint_loci
    pickle_out['fourier']['metric'] = metric_out
    pickle_out['uniform'] = {}
    pickle_out['uniform']['loci'] = breakpoint_loci_uniform
    pickle_out['uniform']['metric'] = metric_out_uniform
    pickle_out['fourier_ls'] = breakpoint_loci_local_search # Yes, breakpoint_loci_local_search is already a dict with 'loci' and 'metrics' keys
    pickle_out['fourier_ls']['metric'] = metric_out_local_search
    pickle_out['uniform_ls'] = breakpoint_loci_uniform_local_search 
    pickle_out['uniform_ls']['metric'] = metric_out_uniform_local_search

    t = datetime.datetime.now()
    t_formatted = t.strftime('%Y_%m_%d_%H_%M_%S')

    # pickle_dump_fname = 'pickle-'+dataset+'-'+chr_name+'-'+str(n_bpoints)+'-'+str(begin)+'-'+str(end)+'-'+t_formatted+'.pickle'
    with open(out_fname, 'wb') as f_out:
        pickle.dump(pickle_out, f_out)
        
    flat.print_log_msg('Done')

def midpoint(a, b):
    if a>b:
        first = b 
        second = a
    else:
        first = a
        second = b

    return first + (second-first)/2 # Takes care of huge values (overflow)

def run_local_search_complete(chr_name, breakpoint_loci, begin, end, input_config, metric_out):
    breakpoint_loci_local_search = {}
    breakpoint_loci_local_search['loci'] = []
    breakpoint_loci_local_search['metrics'] = []
    
    
    total_sum = metric_out['sum']
    total_N = metric_out['N_zero']

    # Search between begin and first midpoint
    b_stop = int(midpoint(breakpoint_loci[0], breakpoint_loci[1])) #-1 # -1 so as to not overlap with next region! -> this is taken care of in local search
    
    new_breakpoint, new_metric = run_local_search_single(chr_name, breakpoint_loci, 0, begin, b_stop, total_sum, total_N, input_config, metric_out)
    breakpoint_loci_local_search['loci'].append(new_breakpoint)
    breakpoint_loci_local_search['metrics'].append(new_metric)
    
    for locus_index in range(1, len(breakpoint_loci)-1):
        b_start = int(midpoint(breakpoint_loci[locus_index-1], breakpoint_loci[locus_index]))
        b_stop = int(midpoint(breakpoint_loci[locus_index], breakpoint_loci[locus_index+1])) #-1 # -1 so as to not overlap with next region! -> this is taken care of in local search
        
        new_breakpoint, new_metric = run_local_search_single(chr_name, breakpoint_loci, locus_index, b_start, b_stop, total_sum, total_N, input_config, metric_out)
        breakpoint_loci_local_search['loci'].append(new_breakpoint)
        breakpoint_loci_local_search['metrics'].append(new_metric)
#         local_search_run = local_search.LocalSearch(chr_name, breakpoint_loci[locus_index-1], breakpoint_loci[locus_index+1], locus_index, breakpoint_loci, total_sum, total_N, input_config)       
#         
#         new_breakpoint, new_metric = local_search_run.search()
#         
#         print_breakpoint_comparison(new_breakpoint, new_metric, breakpoint_loci[locus_index], metric_out)
# #         print(new_breakpoint, new_metric['sum']/new_metric['N_zero'])
# #         print(breakpoint_loci[locus_index], total_sum/total_N)
#         
#         breakpoint_loci_local_search['loci'].append(new_breakpoint)
#         breakpoint_loci_local_search['metrics'].append(new_metric)

    # Search between last midpoint and end
    b_start = int(midpoint(breakpoint_loci[len(breakpoint_loci)-2], breakpoint_loci[len(breakpoint_loci)-1]))
    
    new_breakpoint, new_metric = run_local_search_single(chr_name, breakpoint_loci, len(breakpoint_loci)-1, b_start, end, total_sum, total_N, input_config, metric_out)
    breakpoint_loci_local_search['loci'].append(new_breakpoint)
    breakpoint_loci_local_search['metrics'].append(new_metric)
    
    flat.print_log_msg('New breakpoints:')
    print(breakpoint_loci_local_search)

    return breakpoint_loci_local_search 

def run_local_search_single(chr_name, breakpoint_loci, locus_index, start, stop, total_sum, total_N, input_config, metric_out):
    try:
        local_search_run = local_search.LocalSearch(chr_name, start, stop, locus_index, breakpoint_loci, total_sum, total_N, input_config)       
        
        new_breakpoint, new_metric = local_search_run.search()
        
        print_breakpoint_comparison(new_breakpoint, new_metric, breakpoint_loci[locus_index], metric_out)
        
        return new_breakpoint, new_metric
    except Exception as e:
        flat.print_log_msg('Error!')
        flat.print_log_msg('start: '+repr(start))
        flat.print_log_msg('stop: '+repr(stop))
        flat.print_log_msg('local_search.__dict__: '+repr(local_search.__dict__))
        flat.print_log_msg('Continuing...')
        return breakpoint_loci[locus_index], None


def print_breakpoint_comparison(breakpoint1, metric1, breakpoint2, metric2):
    flat.print_log_msg('Breakpoint 1: '+repr(breakpoint1))
    flat.print_log_msg('Metric 1:')
    print_metric(metric1)
    flat.print_log_msg('Breakpoint 2: '+repr(breakpoint2))
    flat.print_log_msg('Metric 2:')
    print_metric(metric2)

def apply_metric(chr_name, begin, end, input_config, loci): 
    metric_out = metric.Metric(chr_name, input_config, loci, begin, end)
    out = metric_out.calc_metric()
    
    return out

def print_metric(metric_out):
    flat.print_log_msg('Sum: '+repr(metric_out['sum']))
    flat.print_log_msg('N (w/ zero\'s): '+repr(metric_out['N_zero']))
    flat.print_log_msg('Metric: '+repr(metric_out['sum']/metric_out['N_zero']))

def print_opt_arg_error():
    print('For help use --help')
    return(2)

# if __name__ == '__main__':
#     main()

if __name__ == '__main__':
    cl.commander_line((pipeline))
