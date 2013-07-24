#!/bin/bash

export PATH=$PATH:/homes/bjsco/bin

/homes/bjsco/trinityrnaseq_r2013-02-25/Trinity.pl --seqType fq --JM 400G --left /homes/bjsco/iochromas/Project_SmithS_01_SOL/blue_combined/all_blue_R1.fq  --right  /homes/bjsco/iochromas/Project_SmithS_01_SOL/blue_combined/all_blue_R2.fq  --single  /homes/bjsco/iochromas/read/reads.combined.fna --single read/reads.s1r1.fna --output  /homes/bjsco/iochromas/trinity-Icyan_ALL --CPU $NSLOTS --inchworm_cpu $NSLOTS --bflyCPU $NSLOTS
