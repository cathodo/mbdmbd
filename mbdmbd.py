#!/usr/bin/env python3

import sys
import os
import glob
from argparse import ArgumentParser
import math
import joblib
import time
# dataframes
import pandas as pd
import numpy as np
# bars
from tqdm import tqdm
# specialized datatypes
import h5py
import pysam
import ast

def get_args(argv):
    parser = ArgumentParser()
    optional_args = parser._action_groups.pop()
    required_args = parser.add_argument_group('required arguments')
    ### req
    required_args.add_argument('--bamdir', dest='bamdir', required=True, help='directory of bam files output from guppy. PARENT of pass/*.bam')
    required_args.add_argument('--fdir', dest='f5dir', required=True, help='directory of fast5 files input to guppy. PARENT of fast5_pass and fast5_fail directories')
    # optional
    optional_args.add_argument('--verbose', dest='verbose', action='store_true', help='call this to make program loud')
    optional_args.add_argument('--targetbase', dest='targetBase', default='A', help='you can change which base is in the middle of NN.NN kmer motifs that we compare, default A')
    optional_args.add_argument('--contextlen', dest='ctlen', default=1, type=int, help='total size of base context, 5 for 5mer, etc, default 1')
    ##########
    parser._action_groups.append(optional_args)
    parser.set_defaults()
    return parser.parse_args(argv[1:])

def extract_f5_rawseq(f, rn):
    return f[rn]['Raw']['Signal'][...]

def extract_bam_data(r, f5dir):
    read = str(r).split('\t')
    # for getting f5 raw sig
    read_name = "read_" + read[0]
    # move info and fields added to bam via guppy 6.3
    guppy_data = ast.literal_eval(ast.literal_eval(read[11].__repr__().replace('array', '')))
    # indices of signal array which correspond to read (number of 1's = read len)
    move_table = guppy_data[0][1][1]
    # length of each signal segment per base
    stride_len = move_table.pop(0)
    # number of signal ind to start at
    trimmed_samp = guppy_data[8][1]
    # filename for the f5 which corresponds to this (glob should never resolve > 1 file)
    f5_fn = glob.glob(os.path.join(f5dir, "fast5_*", guppy_data[6][1]))[0]
    f = h5py.File(f5_fn, 'r')
    full_sig = extract_f5_rawseq(f, read_name)
    # start and end for the fulllength seq
    sig_s = trimmed_samp
    sig_e = trimmed_samp + len(move_table)*stride_len
    sub_sig = full_sig[sig_s:sig_e]
    # reshaping into stride-len bins
    stride_sig = sub_sig.reshape((sub_sig.size//stride_len, stride_len))
    # get indices of move table
    ss_ind = np.where([b==1 for b in move_table])[0]
    sig = stride_sig[ss_ind]
    ## return seq, sig
    return read[9], sig

def make_kmer_table_from_bamfile(bamfile, f5dir, klen):
    #########################
    ####### iter bams #######
    #########################
    # decide kmer or single base
    negshift = math.floor(klen/2)
    posshift = math.floor(klen/2)+1
    # iter
    big_list = []
    ## this verbosity thing is to try and avoid the missing idx msg 
    ## solution according to https://github.com/pysam-developers/pysam/issues/939#issuecomment-669016051
    save = pysam.set_verbosity(0)
    bam = pysam.AlignmentFile(bamfile, 'rb', check_sq=False)
    pysam.set_verbosity(save)
    for r in bam.fetch(until_eof=True):
        ## just A for now 
        seq_all, sig_all = extract_bam_data(r, f5dir)
        rn = "read_" + str(r).split('\t')[0]
        seq_idx = [args.targetBase == b for b in seq_all]
        iw = np.where(seq_idx)[0]
        ## vector of individual kmer seq/sig pairs
        seq = [seq_all[i-negshift:i+posshift] for i in iw]
        sig = [sig_all[i-negshift:i+posshift] for i in iw]
        # remove truncated
        iw2 = np.where([len(s) == klen for s in seq])[0]
        # check
        if len(seq) == 0 | len(sig) == 0:
            print("ACHTUNG, ME VECTOR IS EMPTY")
            continue
        if len(seq) != len(sig):
            print("AVAST AND CURSES, ME VECTORS' LENGTHS ARE ASKEW")
            continue
        # make serieses
        big_list.append(pd.DataFrame([pd.Series(data={'signal':sig[i].flatten(),'kmer':seq[i],'read_name':rn,'seqloc':str(iw[i])}) for i in iw2]))
    ##############################
    kmer_table = pd.concat(big_list).reset_index()
    del big_list, kmer_table['index']
    return kmer_table


def build_kmer_table(bamdir, f5dir, klen):
    all_bam_files = glob.glob(os.path.join(bamdir, "pass/*.bam"))
    kmer_table_list = joblib.Parallel(n_jobs=50, verbose=10)(joblib.delayed(make_kmer_table_from_bamfile)(file_name,f5dir,klen) for file_name in all_bam_files)
    kmer_table = pd.concat(kmer_table_list).reset_index()
    return kmer_table

def main(argv=sys.argv):
    st = time.time()
    ### handle args
    args = get_args(argv)
    # set bars verbosity
    tqdm.pandas(disable = not args.verbose)
    ### load data
    kmer_table = build_kmer_table(args.bamdir, args.f5dir, args.ctlen)
    print("it's been {} days since you looked at me".format((time.time()-st)/60/60/24))
    # to get raw signal do this
    # X = np.array([x for x in kmer_table['signal']])

if __name__=="__main__":
    main(sys.argv)
