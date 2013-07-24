transcriptome-and-genome-assembly
=================================

Scripts in this repository have been used to assemble de novo transcriptomes and genomes 

KSU_bioinfo_lab folder
----------------------
prinseq.sh - paired end read cleaning with prinseq. If your reads are not illumina cassava (post v1.8) you should comment out lines 28 and 29. This script uses a script called stats.sh to monitor and report resource usage on a SGE cluster like Beocat at KSU. But this section (lines 36-42) also can be commented out.  

stats.sh - see prinseq.sh
 
 
mira_454.sh - Runs Mira 3.4 for 454 data alone. If reads have been cleaned than the trace .xml file extracted from the .sff file will be missing so the parameters in this script are set with mxti=no to reflect the lack of xml trace files.

velvet_oases_k-mer_assemblies.sh - This script will write single k-mer transcriptome assembly scripts and a list of qsubs that can be used to call them on beocat. Comment out line 10 if your sequences are already shuffled. After read cleaning I load my broken pairs (now single end reads) as "-short" and pairs as "-shortPaired". The script calculates memory request for each assembly as "Ram required for velvetg = -109635 + 18977*ReadSize + 86326*GenomeSize + 233353*NumReads - 51092*K" which is taken from http://listserver.ebi.ac.uk/pipermail/velvet-users/2009-June/000359.html. You should change the ReadSize, GenomeSize, and NumReads to reflect your data.

oases_merge.sh - This script calls velvet and oases to merge many single k-mer assemblies (see velvet_oases_k-mer_assemblies.sh).

run_Trinity.sh - This script calls Trinity to assemble with a single k-mer. After read cleaning broken pairs can be loaded as (now single end reads) as "--single".

master_ABySS.pl - this scripts call ABySS_script_writer.pl to write single k-mer genome assembly scripts then executes them (on an SGE cluster). Comment out line 15 and uncomment line 16 to run assemblies every odd k-mer between k=51 and k=71 (the range often seen to have the highest N50).

ABySS_script_writer.pl - see master_ABySS.pl

ABySS_merge_scaffolds.sh - This script merges ASSEMBLIES of illumina LJD libraries ( http://ngs-expert.com/tag/long-jumping-distance-libraries/ ) using parameters listed on in the abyss-pe manual (http://manpages.ubuntu.com/manpages/raring/en/man1/abyss-pe.1.html). Run this after single kmer assemblies have completed (see master_ABySS.pl). 

assembly_quality_stats_for_multiple_assemblies.pl - This script runs a slightly modified version of Joseph Fass' Count_fasta.pl (original available at http://wiki.bioinformatics.ucdavis.edu/index.php/Count_fasta.pl ) on a fasta file from each assembly. It then creates comma separated file called assembly_metrics.csv listing the N25,N50,N75, cumulative contig length, and number of contigs for each assembly (also download Count_fastas.pl and change $path_to_Count_fastas on line 13 of assembly_quality_stats_for_multiple_assemblies.pl).

Count_fastas.pl - see assembly_quality_stats_for_multiple_assemblies.pl



