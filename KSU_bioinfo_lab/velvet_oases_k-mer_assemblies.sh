#!/bin/bash

#  velvet_script_writer.sh
#
#
#  Created by jennifer shelton on 12/18/12.
#
#!/bin/bash
######### shuffle sequences ##########
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


kmem=$((-109635 + 18977*100 + 86326*177 + 233353*476 - 51092*${i}))
mem=$((${kmem}/1000000))
echo "qsub -ckpt dmtcp -l h_rt=100:00:00,mem=${mem}G velv${i}HI.sh"
fi
done