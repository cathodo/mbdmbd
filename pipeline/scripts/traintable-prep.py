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
import pyarrow

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
    optional_args.add_argument('--contextlen', dest='ctlen', default=5, type=int, help='total size of base context, 5 for 5mer, 1 for single base, default 5')
    optional_args.add_argument('--outfile', dest='outfile', help='specify this to name and make output file')
    optional_args.add_argument('--readname_file', dest='readname_file', help='specify this for readnames file data object')
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
    ###### see https://github.com/nanoporetech/bonito/blob/master/documentation/SAM.md#read-tags
    #####################
    ### guppy indices ###
    #####################
    # structure is list of tuples, guppy_data[N] = (key, value)
    # since we want value, we do guppy_data[N][1] mostly
    # move table needs another [1] for some reason
    # list of ind:
    # 0 mv
    # 1 qs
    # 2 mx
    # 3 ch
    # 4 rn
    # 5 st
    # 6 f5
    # 7 ns
    # 8 ts
    # 9 sm
    # 10 sd
    # 11 sv
    # 12 RG
    ######### CORE EXTRACTION #########
    # indices of signal array which correspond to read (number of 1's = read len)
    move_table = guppy_data[0][1][1]
    # length of each signal segment per base
    stride_len = move_table.pop(0)
    # number of signal ind to start at
    trimmed_samp = guppy_data[8][1]
    # filename for the f5 which corresponds to this (glob should never resolve > 1 file)
    f5_fn = glob.glob(os.path.join(f5dir, "fast5_*", guppy_data[6][1]))[0]
    ####### optional exraction ########
    start_time = guppy_data[5][1]
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
    #####################################################
    # output, return optional as needed
    return seq, sig, start_time

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
        seq_all, sig_all, st = extract_bam_data(r, f5dir)
        ## handle bam data
        rn = "read_" + r.qname
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
            'start_time':st
            }) for i in iw2]))
    ##############################
    kmer_table = pd.concat(big_list).reset_index()
    del big_list, kmer_table['index']
    return kmer_table

def check_kmer(x, seq_all):
    k = x['kmer'] 
    s = int(x['move_start'])
    e = int(x['move_end'])
    k2 = seq_all[s:e]
    if k == k2:
        return True
    else:
        return False

def get_signal(x, sig_all):
    s = int(x['move_start'])
    e = int(x['move_end'])
    return sig_all[s:e]

def build_kmer_table_from_coordinates(bamfile, f5dir, readinfo, verbose):
    #########################
    ####### iter bams #######
    #########################
    unique_reads = list(readinfo['read_name'].drop_duplicates())
    # coord data preset kmer size
    negshift = 2
    posshift = 3
    # iter
    big_list = []
    ## this verbosity thing is to try and avoid the missing idx msg 
    ## solution according to https://github.com/pysam-developers/pysam/issues/939#issuecomment-669016051
    if not verbose:
        save = pysam.set_verbosity(0)
    bam = pysam.AlignmentFile(bamfile, 'rb', check_sq=False)
    if not verbose:
        pysam.set_verbosity(save)
    # check if any reads match
    any_match = False
    for r in bam.fetch(until_eof=True):
        if r.qname in unique_reads:
            any_match = True
            ## get bam data
            seq_all, sig_all, st = extract_bam_data(r, f5dir)
            ## handle bam data
            rn = "read_" + r.qname
            ##### check kmer match
            readmatch = readinfo[readinfo['read_name']==r.qname]
            matchidx = readmatch.apply(lambda x: check_kmer(x, seq_all), axis=1)
            readhits = readmatch[matchidx].reset_index()
            if len(readhits) == 0:
                # continue skips this iter if none of the kmers match
                continue
            readhits.insert(len(readhits.columns), "signal", readhits.apply(lambda x: get_signal(x, sig_all), axis=1))
            # make serieses
            big_list.append(pd.DataFrame([pd.Series(data={
                'signal':readhits.loc[i]['signal'].flatten(),
                'kmer':readhits.loc[i]['kmer'],
                'read_name':rn,
                'seqloc':str(readhits.loc[i]['move_start']),
                'start_time':st
                }) for i in range(len(readhits))]))
    ##############################
    # big list is empty in bamfile misses
    if any_match:
        kmer_table = pd.concat(big_list).reset_index()
        del big_list, kmer_table['index']
        return kmer_table

def get_kmer_table_coord(readinfo, argv=sys.argv):
    st = time.time()
    ### handle args
    args = get_args(argv)
    # set bars verbosity
    tqdm.pandas(disable = not args.verbose)
    ### load data
    all_bam_files = glob.glob(os.path.join(args.bamdir, "pass/*.bam"))
    kmer_table_list = joblib.Parallel(n_jobs=50, verbose=10)(
        joblib.delayed(build_kmer_table_from_coordinates)
            (file_name,
            args.f5dir,
            readinfo,
            args.verbose) for file_name in all_bam_files
    )
    kmer_table = pd.concat(kmer_table_list).reset_index()
    print("it's been {} days since you looked at me".format((time.time()-st)/60/60/24))
    # NOTE to get raw signal do this:
    # X = np.array([x for x in kmer_table['signal']])
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
    if args.readname_file:
        readnames = pd.read_csv(args.readname_file)
        kmer_table = get_kmer_table_coord(readnames, argv)
    else:
        kmer_table = build_kmer_table(args)
    print("it's been {} days since you looked at me".format((time.time()-st)/60/60/24))
    # NOTE to get raw signal do this:
    # X = np.array([x for x in kmer_table['signal']])
    if args.outfile:
        kmer_table.to_parquet(args.outfile, engine='pyarrow')
    else:
        return kmer_table

if __name__=="__main__":
    kmer_table = get_kmer_table()
