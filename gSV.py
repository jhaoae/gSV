import sys
import os
import pysam
from datetime import datetime
import argparse
import multiprocessing
import json
import numpy as np
from main_function import p_encode, graph_cut, detect, summarize

def parse_arguments(arguments = sys.argv[1:]):
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                    description="""GenSV {0} \n \nShort Usage: GenSV [parameters] -o <output path> -b <input bam path> -r <reference>""")

    required_params = parser.add_argument_group("Input/Output parameters")
    required_params.add_argument('-r', dest="ref_path",type = os.path.abspath, required=True, help='Path to reference')
    required_params.add_argument('-b', dest="bam_path",type = os.path.abspath, required=True, help='Path to bam file')
    required_params.add_argument('-o', dest="out_path",type = os.path.abspath, required=True, help='Path to output')

    optional_params = parser.add_argument_group("Parallel parameters")
    optional_params.add_argument('-t', dest="thread_num", type=int, default=1, help='Thread numbers (default: %(default)s)')
    optional_params.add_argument('--grachcut_thread', dest="graph_cut_thread_num", type=int, default=1, help='Thread numbers for Graph Cut (default: %(default)s)')

    encoding_params = parser.add_argument_group("Parameters related to encoding step")
    encoding_params.add_argument('--min_mapq', dest="min_mapping_quality", type=int, default=10, help='Minimum mapping quality required for encoding alignment results (default: %(default)s)')
    encoding_params.add_argument('--min_sig_len', dest="min_signature_length", type=int, default=10, help='Minimum length for signature recording during encoding step (default: %(default)s)')
    
    graphcut_params = parser.add_argument_group("Parameters related to detection (graphcut) step")
    graphcut_params.add_argument('--beta', dest="beta", type=int, default=5, help='Graphcut-related parameter beta (default: %(default)s)')
    graphcut_params.add_argument('--gamma', dest="gamma", type=int, default=1, help='Graphcut-related parameter gamma (default: %(default)s)')
    graphcut_params.add_argument('--amplify', dest="amplify", type=int, default=7, help='Graphcut-related parameters amplify (default: %(default)s)')
    graphcut_params.add_argument('--combine_distance', dest="combine_distance", type=int, default=500, help='Threshold for merging candidate regions detected by graphcut (default: %(default)s bps)')
    
    complex_params = parser.add_argument_group("Parameters related to determining whether a candidate region is COMPLEX or SIMPLE")
    complex_params.add_argument('--Complex', dest="complex_mode", type=bool, default=False, help='All candidate regions enter COMPLEX path (default: %(default)s)')
    complex_params.add_argument('--position_range', dest="pos_range", type=int, default=100, help='Extended range of candidate regions (default: %(default)s bps)')
    complex_params.add_argument('-L', dest="L", type=float, default=3, help='The signal length is a fraction (1/L) of the candidate region (default: %(default)s)')
    complex_params.add_argument('--support_signal', dest="support_signal", type=int, default=3, help='Minimum support number of signal in a candidate region (default: %(default)s)')
    

    cluster_params = parser.add_argument_group("Parameters related to clustering")
    cluster_params.add_argument('--cluster_similarity', dest="cluster_similarity", type=float, default=0.9, help='Similarity threshold for clustering two reads together (default: %(default)s)')

    
    realign_params = parser.add_argument_group("Parameters related to realigning")
    realign_params.add_argument('--mempath', dest="mempath", type = os.path.abspath, required=True, help='Path of CopMEM2')
    realign_params.add_argument('--memlen', dest="memlen", type=float, default=100, help='Minimal MEM length (default: %(default)s). This value must be >= 50')

    final_params = parser.add_argument_group("Parameters related to finalizing")
    final_params.add_argument('--min_SV_len', dest="min_SV_len", type=int, default=50, help='Minimal SV length (default: %(default)s bps)')

    
    options = parser.parse_args(arguments)

    return options


if __name__ == '__main__':
    options = parse_arguments()

    out_put = options.out_path
    if not os.path.exists(out_put):
        os.mkdir(out_put)

    ref_path = options.ref_path
    bam_path = options.bam_path
    out_path = options.out_path

    #multiprocessing
    chromosome = []
    bam_file = pysam.AlignmentFile(bam_path)
    for r in bam_file:
        if r.reference_name not in chromosome:
            chromosome.append(r.reference_name)

    currentDateAndTime = datetime.now()
    print("The current date and time is", currentDateAndTime) 
    print(chromosome) 
    ''' 
    #Encoding
    process_pool = multiprocessing.Pool(processes=options.thread_num)
    pool_rets = []
    
    for chrom in chromosome:
        pool_rets.append([process_pool.apply_async(p_encode, (chrom, bam_path, ref_path, out_path, options.min_mapping_quality, options.min_signature_length)), chrom])
    
    process_pool.close()
    process_pool.join()
    
    
    #Graph_cut
    process_pool = multiprocessing.Pool(processes=options.graph_cut_thread_num)
    pool_rets = []
    for chrom in chromosome:
        pool_rets.append([process_pool.apply_async(graph_cut,(chrom, out_path, options.combine_distance, options.beta, options.gamma, options.amplify)),chrom]) 
    process_pool.close()
    process_pool.join()
    ''' 
    
    #Detection
    process_pool = multiprocessing.Pool(processes=options.thread_num)
    pool_rets = []
    detect('20', bam_path, ref_path, out_path, options.complex_mode, options.pos_range, options.L, options.support_signal, options.cluster_similarity, options.mempath, options.memlen, options.min_SV_len)
    
    for chrom in chromosome:
        pool_rets.append([process_pool.apply_async(detect, (chrom, bam_path, ref_path, out_path, options.complex_mode, options.pos_range, options.L, options.support_signal, options.cluster_similarity, options.mempath, options.memlen, options.min_SV_len)), chrom])

    process_pool.close()
    process_pool.join()
      
    #Summarize
    summarize(ref_path, out_path, options.min_SV_len)
   
    
