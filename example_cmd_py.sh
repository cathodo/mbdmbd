#!/usr/bin/bash

bamdir="/storage/common/rebasecall/p48_anti-adar_guppy6.4"
fdir="/data/nanopore/data/p48-51_20210210_Adar-Abeta-Aza_R1_SQK_RNA002_HT22/P48_Anti-adar/20210210_2130_1-A1-D1_PAE86248_66c1fc9e/"

python -i mbdmbd.py --bamdir $bamdir --fdir $fdir
