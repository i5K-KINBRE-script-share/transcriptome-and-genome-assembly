#!/bin/bash
set -o verbose
PATH=/homes/sheltonj/abjc/velvet_1.2.08:/homes/sheltonj/abjc/oases_0.2.08:${PATH}
export PATH
cd /homes/sheltonj/bb_velvet
velveth mergedAssembly 27 -long velv_bb_*/transcripts.fa
velvetg mergedAssembly -read_trkg yes -conserveLong yes
oases mergedAssembly -merge yes
### qsub -l h_rt=100:00:00,mem=100G /homes/sheltonj/bb_scripts/oases_merge.sh ###