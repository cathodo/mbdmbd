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
    # for reporting sequence
    seq = read[9]
    # move info and fields added to bam via guppy 6.4
    guppy_data = ast.literal_eval(ast.literal_eval(read[11].__repr__().replace('array', '')))
    ######################################################
    ############# PARSE the added fields #################
    ######################################################
    for g in guppy_data:
        code = g[0]
        value = g[1]
        # see https://github.com/nanoporetech/bonito/blob/master/documentation/SAM.md#read-tags
        if code == "mv":
            # indices of signal array which correspond to read (number of 1's = read len)
            move_table = value[1]
            # length of each signal segment per base
            stride_len = move_table.pop(0)
        elif code == "ts":
            # number of signal ind to start at
            trimmed_samp = value
        elif code == "f5":
            # filename for the f5 which corresponds to this (glob should never resolve > 1 file)
            f5_fn = glob.glob(os.path.join(f5dir, "fast5_*", value))[0]
        elif code == "sv":
            # make sure it's the scaling we expect
            if value != "med_mad":
                print("Warning: unexpected scaling version, data scaling may not make sense")
        elif code == "sm":
            # scaling midpoint
            scale_mid = value
        elif code == "sd":
            scale_dis = value
        elif code == "st":
            start_time = value
        else:
            continue
    #####################################################
    ########### reshape data into sig + seq #############
    #####################################################
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
    # output
    return seq, sig, scale_mid, scale_dis, start_time

def make_kmer_table_from_bamfile(bamfile, f5dir, targetBase, klen, verbose):
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
    if not verbose:
        save = pysam.set_verbosity(0)
    bam = pysam.AlignmentFile(bamfile, 'rb', check_sq=False)
    if not verbose:
        pysam.set_verbosity(save)
    for r in bam.fetch(until_eof=True):
        ## get bam data
        seq_all, sig_all, sm, sd, st = extract_bam_data(r, f5dir)
        ## handle bam data
        rn = "read_" + str(r).split('\t')[0]
        seq_idx = [targetBase == b for b in seq_all]
        iw = np.where(seq_idx)[0]
        ## vector of individual kmer seq/sig pairs
        seq = [seq_all[i-negshift:i+posshift] for i in iw]
        sig = [sig_all[i-negshift:i+posshift] for i in iw]
        # remove truncated
        iw2 = np.where([len(s) == klen for s in seq])[0]
        # check
        if len(seq) == 0 | len(sig) == 0:
            if verbose:
                print("ACHTUNG, ME VECTOR IS EMPTY")
            continue
        if len(seq) != len(sig):
            if verbose:
                print("AVAST AND CURSES, ME VECTORS' LENGTHS ARE ASKEW")
            continue
        # make serieses
        big_list.append(pd.DataFrame([pd.Series(data={
            'signal':sig[i].flatten(),
            'kmer':seq[i],
            'read_name':rn,
            'seqloc':str(iw[i]),
            'scale_mid':sm,
            'scale_dis':sd,
            'start_time':st
            }) for i in iw2]))
    ##############################
    kmer_table = pd.concat(big_list).reset_index()
    del big_list, kmer_table['index']
    return kmer_table

def build_kmer_table(args):
    all_bam_files = glob.glob(os.path.join(args.bamdir, "pass/*.bam"))
    kmer_table_list = joblib.Parallel(n_jobs=50, verbose=10)(
        joblib.delayed(make_kmer_table_from_bamfile)
            (file_name,
            args.f5dir,
            args.targetBase,
            args.ctlen,
            args.verbose) for file_name in all_bam_files
    )
    kmer_table = pd.concat(kmer_table_list).reset_index()
    return kmer_table

def get_kmer_table(argv=sys.argv):
    st = time.time()
    ### handle args
    args = get_args(argv)
    # set bars verbosity
    tqdm.pandas(disable = not args.verbose)
    ### load data
    kmer_table = build_kmer_table(args)
    print("it's been {} days since you looked at me".format((time.time()-st)/60/60/24))
    # NOTE to get raw signal do this:
    # X = np.array([x for x in kmer_table['signal']])
    return kmer_table

if __name__=="__main__":
    kmer_table = get_kmer_table(sys.argv)
