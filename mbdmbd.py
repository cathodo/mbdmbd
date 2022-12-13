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
from dask import dataframe as dd
from dask import compute
# bars
from tqdm import tqdm
from dask.diagnostics import ProgressBar
# stats
from scipy import stats
from sklearn.mixture import GaussianMixture
# kmeans
from kneed import KneeLocator
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler
# specialized datatypes
import h5py
import pysam
import ast
# plots
from matplotlib import pyplot as plt
import seaborn as sns

def get_args(argv):
    parser = ArgumentParser()
    optional_args = parser._action_groups.pop()
    required_args = parser.add_argument_group('required arguments')
    ### req
    required_args.add_argument('--dirA', dest='f5dirA', required=True, help='directory of fast5 files from sample A')
    required_args.add_argument('--dirB', dest='f5dirB', required=True, help='directory of fast5 files from sample B')
    required_args.add_argument('--output', dest='outDir', required=True, help='where the output files are to be stored')
    # optional
    optional_args.add_argument('--verbose', dest='verbose', action='store_true', help='call this to make program loud')
    optional_args.add_argument('--targetbase', dest='targetBase', default='A', help='you can change which base is in the middle of NN.NN kmer motifs that we compare, default A')
    optional_args.add_argument('--test', dest='test', action='store_true', help='run on only one fast5 as a test')
    ##########
    parser._action_groups.append(optional_args)
    parser.set_defaults()
    return parser.parse_args(argv[1:])

def extract_f5_data(f, rn):
    sequence = str(f[rn]['Analyses']['Basecall_1D_001']['BaseCalled_template']['Fastq'][...]).split("\\n")[1]
    trace = f[rn]['Analyses']['Basecall_1D_001']['BaseCalled_template']['Trace'][...]
    move = f[rn]['Analyses']['Basecall_1D_001']['BaseCalled_template']['Move'][...]
    # return trace where move == 1, so we only get basecalls present in the seq
    return sequence.replace("U","T"), trace[move==1]

def extract_f5_rawseq(f, rn):
    return f[rn]['Raw']['Signal'][...]

def extract_bam_data(r, filedir):
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
    f5_fn = glob.glob(os.path.join(filedir, "workspace/fast5_*", guppy_data[6][1]))[0]
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

def make_kmer_table_from_bamfile(bamfile, filedir, klen):
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
        seq_all, sig_all = extract_bam_data(r, filedir)
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
            print("ACHTUN, ME VECTOR IS EMPTY")
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


def k_means_cluster_plot(X, plotPath, numClust):
    # scale data (mean to 0 stdev to 1)
    scaler = StandardScaler()
    scaled_features = scaler.fit_transform(X)
    # choosing best clusters
    kmeans_kwargs = {
        'init': "random",
        'n_init': 10,
        'max_iter': 300,
        'random_state': 42
    }
    # sse vals for each k
    sse = []
    for k in range(1,numClust):
        kmeans = KMeans(n_clusters=k, **kmeans_kwargs)
        kmeans.fit(scaled_features)
        sse.append(kmeans.inertia_)
    # plot
    plt.style.use("fivethirtyeight")
    plt.plot(range(1, numClust), sse)
    plt.xticks(range(1, numClust))
    plt.xlabel("Number of Clusters")
    plt.ylabel("SSE")
    plt.savefig(plotPath)
    return sse

def build_kmer_table(filedir, bamstr, klen):
    all_bam_files = glob.glob(os.path.join(filedir, bamstr))
    kmer_table_list = joblib.Parallel(n_jobs=50, verbose=10)(joblib.delayed(make_kmer_table_from_bamfile)(file_name,filedir,klen) for file_name in all_bam_files)
    kmer_table = pd.concat(kmer_table_list).reset_index()
    return kmer_table

#def main(argv=sys.argv):
if True:
    st = time.time()
    argv=sys.argv
    ###
    args = get_args(argv)
    # set bars verbosity
    tqdm.pandas(disable = not args.verbose)
    ### output storage
    if not os.path.exists(args.outDir):
        os.mkdir(args.outDir)
    ### load data
    # decide bam
    bamstr = "pass/*_11*.bam" if args.test else "pass/*.bam"
    kmer_table = build_kmer_table(args.f5dirA, bamstr, 1)
    print("it's been {} days since you looked at me".format((time.time()-st)/60/60/24))
    X = np.array([x for x in kmer_table['signal']])
    exit()
if False:
    ########################################
    ############### analysis ###############
    ########################################
    # https://realpython.com/k-means-clustering-python/
    # Nondeterministic machine learning algorithms like k-means are difficult to reproduce. 
    # The random_state parameter is set to an integer value so you can follow the data presented in the tutorial. 
    # In practice, it’s best to leave random_state as the default value, None.
    #kmer_subset = kmer_table[kmer_table['kmer']=='CGACG']
    #
    k_attempt = 21
    sse = k_means_cluster_plot(X, os.path.join(args.outDir, "test_big_clust.png"), k_attempt)
    # find knee programmatically
    kneedle = KneeLocator(range(len(sse)), sse, S=1.0, curve='convex', direction='decreasing')
    print(args.targetBase, kneedle.knee)
    #
if True:
    ### gmm ###
    gmm = GaussianMixture(3, covariance_type='full', random_state=0).fit(X)
    #kmer_subset['predicted_clust'] = gmm.predict(X)
    #
    ### plot for sanity checking ###
    plotdata = pd.DataFrame(data={'signal_mean': np.mean(X, axis=1), 'cluster': gmm.predict(X)})
    plot = sns.violinplot(data = plotdata, x='cluster', y='signal_mean')
    fig = plot.get_figure()
    fig.savefig(os.path.join(args.outDir, "test_big_A.png"))

#if __name__=="__main__":
#    main(sys.argv)
