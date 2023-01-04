#!/usr/bin/env python

import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--in', dest='path', default=None)
parser.add_argument('--out', dest='out', default=None)
args = parser.parse_args()

######

col_list = ["id", "position", "kmer", "diff_mod_rate_A_vs_B", "mod_rate_A-rep1", "mod_rate_B-rep1", "pval_A_vs_B"]
xpore = pd.read_csv(args.path, delimiter=',', usecols=col_list)
###############
### filters ###
###############
#kmer
filter1 = xpore['kmer'].str.contains(r'..A..')
# exteme modrates
mr_thresh_lower = 0.001
mr_thresh_upper = 0.99
filter2 = (xpore['mod_rate_A-rep1'] < mr_thresh_upper) & (xpore['mod_rate_A-rep1'] > mr_thresh_lower)
filter3 = (xpore['mod_rate_B-rep1'] < mr_thresh_upper) & (xpore['mod_rate_B-rep1'] > mr_thresh_lower)
# middling diffmod
dm_thresh = 0.6
filter4 = (xpore['diff_mod_rate_A_vs_B'] > dm_thresh) & (xpore['diff_mod_rate_A_vs_B'] < -dm_thresh)
# pval
pval_thresh = 0.01
filter5 = (xpore['pval_A_vs_B'] < pval_thresh)
### apply filters ###
#NOTE: choose which ones you want here. 1 & 5 for B2 at moment
xporeFilter = (filter1 & filter5)
xpore = xpore[xporeFilter].reset_index()
xpore['index'] = xpore.index.values

xpore.to_csv(args.out)
