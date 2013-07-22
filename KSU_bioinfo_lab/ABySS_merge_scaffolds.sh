#!/bin/sh
export PATH=$(find /homes/bjsco/abyss-1.3.4 -type d | tr '
' ':' | sed 's/:$//'):${PATH}
#  ABySS_merge_scaffolds.sh
#
#
#  Created by jennifer shelton on 7/22/13.
#
## Script to merge ASSEMBLIES of illumina LJD libraries http://ngs-expert.com/tag/long-jumping-distance-libraries/ using parameters listed on http://manpages.ubuntu.com/manpages/raring/en/man1/abyss-pe.1.html ##
## Script runs on beocat

kmer=63
mkdir Kmer_merge
cd Kmer_merge
/homes/bjsco/abyss-1.3.4/bin/abyss-pe name=Kmer_merge scaffolds k=${kmer} n=10 N=10 np=32 se="/path/to/scaffolds/Kmer_opt*-scaffolds.fa"

## qsub -l h_rt=300:00:00,mem=12G -pe single 32  ABySS_merge_scaffolds.sh