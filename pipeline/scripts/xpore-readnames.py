#!/usr/bin/env python3

import sys
import os
from argparse import ArgumentParser
import time
import re
import glob
# big data
import pandas as pd
# to disable copy warning cuz i'm overwriting df in the bam function
pd.options.mode.chained_assignment = None
import numpy as np
from dask import dataframe as dd
import dask
# bioinfo
import pysam
# bars
from tqdm import tqdm

def get_args(argv):
    parser = ArgumentParser()
    optional_args = parser._action_groups.pop()
    required_args = parser.add_argument_group('required arguments')
    ### required
    required_args.add_argument('--xpore_table', dest='xporeFn', required=True, help='output from xpore: the diffmod.table file')
    required_args.add_argument('--sampleA_dir', dest='sampleAdir', required=True, help='output from xpore pipeline: the folder for the individual sample A')
    required_args.add_argument('--sampleB_dir', dest='sampleBdir', required=True, help='output from xpore pipeline: the folder for the individual sample B')
    ### optional
    optional_args.add_argument('--verbose', dest='verbose', action='store_true', help='use to give more detailed output msgs')
    optional_args.add_argument('--outfile', dest='outfile', help='filename to write xpore readnames to')
    ################
    parser._action_groups.append(optional_args)
    parser.set_defaults()
    return parser.parse_args(argv[1:])

def get_xpipe_outpaths(dir):
    npevel_path = os.path.join(dir,"np","eventalign.txt")
    npsumm_path = os.path.join(dir,"np","summary.txt")
    bam_path = glob.glob(os.path.join(dir,"bam","*.bam"))[0]
    return npevel_path, npsumm_path, bam_path

def eventalign_cols(path):
    col_list = ["contig", "position", "reference_kmer", "read_index"]
    col_typedict = {'contig': str, 'position': int, 'reference_kmer': str, 'read_index': int}
    eventalign = dd.read_csv(path, delimiter='\t', usecols = col_list, dtype=col_typedict)
    return eventalign

def summary_cols(path):
    col_list = ["read_index", "read_name", "fast5_path"]
    col_typedict = {'read_index': int, 'read_name': str, 'fast5_path': str}
    summary = dd.read_csv(path, delimiter='\t', usecols=col_list, dtype=col_typedict)
    return summary

def eventalign_filter(eventy, xpore):
    filter1 = eventy['contig'].isin(xpore['id'])
    filter2 = eventy['position'].isin(xpore['position'])
    filter3 = eventy['reference_kmer'].isin(xpore['kmer'])
    return eventy[filter1 & filter3]

def summary_filter(summ, ridx):
    return summ[summ['read_index'].isin(ridx)]

def column_renamer(rn):
    rn.columns = rn.columns.str.replace('contig', 'id')
    rn.columns = rn.columns.str.replace('reference_kmer', 'kmer')
    return rn

def attach_bam_data(rn, bamfile):
    unique_reads = list(rn['read_name'].drop_duplicates())
    ## this verbosity thing is to try and avoid the missing idx msg 
    ## solution according to https://github.com/pysam-developers/pysam/issues/939#issuecomment-669016051
    save = pysam.set_verbosity(0)
    bam = pysam.AlignmentFile(bamfile, 'rb', check_sq=False)
    pysam.set_verbosity(save)
    for r in bam.fetch(until_eof=True):
        if r.qname in unique_reads:
            # rows to touch
            ridx = np.where(rn['read_name'] == r.qname)[0]
            # get first cig index, where the mapped region starts in the read seq
            cigsidx = r.cigar[0][1]
            # cig start + position of row - reference_start
            kmer_starts = cigsidx + rn.loc[ridx,'position'] - r.reference_start
            # send out
            rn.loc[ridx,'move_start'] = kmer_starts
            rn.loc[ridx,'move_end'] = kmer_starts+5
    return rn

def get_xpore_readnames(argv=sys.argv):
    st = time.time()
    args = get_args(argv)
    # set bars verbosity
    tqdm.pandas(disable = not args.verbose)
    # filenames
    xpore_path = args.xporeFn
    evalA, summA, bamA = get_xpipe_outpaths(os.path.dirname(args.sampleAdir))
    evalB, summB, bamB = get_xpipe_outpaths(os.path.dirname(args.sampleBdir))
    # filtered xpore table
    xpore = pd.read_csv(args.xporeFn)
    # filtered nanopolish eventaligns
    evA = eventalign_filter(eventalign_cols(evalA), xpore).compute()
    evB = eventalign_filter(eventalign_cols(evalB), xpore).compute()
    # get readnames from summary
    smA = summary_filter(summary_cols(summA), evA['read_index']).compute()
    smB = summary_filter(summary_cols(summB), evB['read_index']).compute()
    # merge 1
    rnA = evA.merge(smA, on='read_index', how='outer').drop_duplicates().reset_index()
    rnB = evB.merge(smB, on='read_index', how='outer').drop_duplicates().reset_index()
    del rnA['index'], rnB['index']
    # get mapping pos from bamfile
    rnA[['move_start','move_end']] = np.nan
    rnB[['move_start','move_end']] = np.nan
    rnA = attach_bam_data(rnA, bamA)
    rnB = attach_bam_data(rnB, bamB)
    # rename cols so they match xpore format
    rnA = column_renamer(rnA)
    rnB = column_renamer(rnB)
    print("that took so long, I've been waiting forever. {} whole seconds".format(time.time()-st))
    outDF = pd.concat([rnA, rnB])
    if args.outfile:
        outDF.to_csv(args.outfile, index=False)
    else:
        return outDF

if __name__=="__main__":
    get_xpore_readnames()
    
