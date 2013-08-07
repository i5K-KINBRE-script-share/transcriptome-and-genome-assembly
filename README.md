transcriptome-and-genome-assembly
=================================

Scripts in this repository have been used to assemble de novo transcriptomes and genomes 

KSU_bioinfo_lab folder
----------------------

Count_fastas.pl - see assembly_quality_stats_for_multiple_assemblies.pl

assembly_quality_stats_for_multiple_assemblies.pl - This script runs a slightly modified version of Joseph Fass' Count_fasta.pl (original available at http://wiki.bioinformatics.ucdavis.edu/index.php/Count_fasta.pl ) on a fasta file from each assembly. It then creates comma separated file called assembly_metrics.csv listing the N25,N50,N75, cumulative contig length, and number of contigs for each assembly (also download Count_fastas.pl and change $path_to_Count_fastas on line 13 of assembly_quality_stats_for_multiple_assemblies.pl).

pre_post_cleaning_metrics.pl - This script summarizes average read lengths etc before and after cleaning for multiple single end or paired end files. Output log files from prinseq (v prinseq-lite-0.20.3) are taken as the input and cleaning metrics are written to the file pre_post_clean_reads.csv.

run_Trinity.sh - This script calls Trinity to assemble with a single k-mer. After read cleaning broken pairs can be loaded as (now single end reads) as "--single".

run_abyss.pl - This script runs an installation of ABySS-1.3.4 on beocat. Step 1 is the assembly of multiple single k-mer assemblies. Step 2 is the merging of these single k-mer assemblies by running ABySS on reads and the unitigs from the multi-kmer assemblies. Step 3 writes a script to calculate N25,N50,N75, cumulative scaffold length, number of scaffolds for all assemblies.

run_broken_pair_removal.sh - Runs the broken pair removal script from the velvet package on multiple files. This script can convert illumina headers so that they end in "/1" or "/2" if your reads are not illumina cassava (post v1.8) you should comment out lines 28 and 29. This script uses a script called stats.sh to monitor and report resource usage on a SGE cluster like Beocat at KSU. But this section (lines 36-42) also can be commented out.

run_mira_454.sh - Runs Mira 3.4 for 454 data alone. If reads have been cleaned than the trace .xml file extracted from the .sff file will be missing so the parameters in this script are set with mxti=no to reflect the lack of xml trace files.

run_velvet_oases_k-mer_assemblies.sh - This script will write single k-mer transcriptome assembly scripts and a list of qsubs that can be used to call them on beocat. Comment out line 10 if your sequences are already shuffled. After read cleaning I load my broken pairs (now single end reads) as "-short" and pairs as "-shortPaired". The script calculates memory request for each assembly as "Ram required for velvetg = -109635 + 18977*ReadSize + 86326*GenomeSize + 233353*NumReads - 51092*K" which is taken from http://listserver.ebi.ac.uk/pipermail/velvet-users/2009-June/000359.html. You should change the ReadSize, GenomeSize, and NumReads to reflect your data.

run_oases_merge.sh - This script calls velvet and oases to merge many single k-mer assemblies (see velvet_oases_k-mer_assemblies.sh).

run_prinseq.sh - paired end read cleaning with prinseq. This script can convert illumina headers so that they end in "/1" or "/2" if your reads are not illumina cassava (post v1.8) you should comment out lines 28 and 29. This script uses a script called stats.sh to monitor and report resource usage on a SGE cluster like Beocat at KSU. But this section (lines 36-42) also can be commented out.  

stats.sh - see prinseq.sh
 
 





