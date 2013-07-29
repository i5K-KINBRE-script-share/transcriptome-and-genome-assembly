#!/bin/bash

#  run_velvet_oases_k-mer_assemblies.sh
#  This script will write single k-mer transcriptome assembly scripts and a list of qsubs that can be used to call them on beocat. Comment out line 10 if your sequences are already shuffled. After read cleaning I load my broken pairs (now single end reads) as "-short" and pairs as "-shortPaired". The script calculates memory request for each assembly as "Ram required for velvetg = -109635 + 18977*ReadSize + 86326*GenomeSize + 233353*NumReads - 51092*K" which is taken from http://listserver.ebi.ac.uk/pipermail/velvet-users/2009-June/000359.html. You should change the ReadSize, GenomeSize, and NumReads to reflect your data.
#
#  Created by jennifer shelton on 12/18/12.
#
#!/bin/bash
######### shuffle sequences (if your pairs are unbroken but in two fastq files) ##########
perl /homes/sheltonj/abjc/velvet_1.2.08/contrib/shuffleSequences_fasta/shuffleSequences_fastq.pl HB_IB_R1_paired.fastq HB_IB_R2_paired.fastq HB_IB_shuff_pairs.fastq
######### assemble ##########
for i in {25..61}
do
out=$(( $i % 2 ))
if [ $out -eq 1 ]
then
echo "#!/bin/bash" > velv${i}HI.sh
echo "set -o verbose" >> velv${i}HI.sh
echo "PATH=/homes/sheltonj/abjc/velvet_1.2.08:/homes/sheltonj/abjc/oases_0.2.08:${PATH}" >> velv${i}HI.sh
echo "export PATH" >> velv${i}HI.sh
echo "cd bb_velvet" >> velv${i}HI.sh
echo "velveth velv_bb_${i} ${i} -fastq -short HB_IB_single.fastq -shortPaired -interleaved -fastq HB_IB_shuff_pairs.fastq " >> velv${i}HI.sh
echo "velvetg velv_bb_${i} -read_trkg yes -ins_length 240" >> velv${i}HI.sh
echo "oases velv_bb_${i}" >> velv${i}HI.sh

######### estimates memory requirements and write qsubs for beocat ###
kmem=$((-109635 + 18977*100 + 86326*177 + 233353*476 - 51092*${i}))
mem=$((${kmem}/1000000))
echo "qsub -ckpt dmtcp -l h_rt=100:00:00,mem=${mem}G velv${i}HI.sh"
fi
done