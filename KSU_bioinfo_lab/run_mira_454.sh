#!/bin/sh
########## Runs Mira 3.4 for the 454 data alone ###########
cd $HOME/bb_pyro/
/homes/bioinfo/bioinfo_software/mira/bin/mira --project=bb_pyro --job=denovo,est,accurate,454 COMMON_SETTINGS -GE:not=4 454_SETTINGS -LR:lsd=yes:ft=fastq:mxti=no