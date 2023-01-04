#!/bin/bash

sampleAin="/data/p48-51_20210210_Adar-Abeta-Aza_R1_SQK_RNA002_HT22/P49_Rev-Aza/20210210_2130_1-E1-H1_PAE86504_0221a67e/"
sampleBin="/data/p48-51_20210210_Adar-Abeta-Aza_R1_SQK_RNA002_HT22/P51_Rev-DMSO/20210210_2130_2-A3-D3_PAE87936_546f1bab/"
sampleAout="/data/p49_Rev-Aza_guppy6.4/"
sampleBout="/data/p51_Rev-DMSO_guppy6.4/"

if guppy_basecaller --version | grep -q "6.4"; then
  mkdir $sampleAout $sampleBout
  guppy_basecaller --flowcell FLO-PRO002 --kit SQK-RNA002 -i $sampleAin -r -s $sampleAout -x 'cuda:0,1,2,3' --trim_adapters --bam_out --moves_out
  guppy_basecaller --flowcell FLO-PRO002 --kit SQK-RNA002 -i $sampleBin -r -s $sampleBout -x 'cuda:0,1,2,3' --trim_adapters --bam_out --moves_out
else
  echo "guppy basecaller not version 6.4, you need that for fast5 file formatting"
fi
