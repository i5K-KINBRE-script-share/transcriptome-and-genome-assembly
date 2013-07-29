#!/usr/bin/perl

use warnings;
use strict;
use 5.010;
#  run_abyss.pl
####### Modified from a script by Bradley Olson http://www.k-state.edu/olsonlab/projects.html ######
#  USAGE: perl run_abyss.pl ## alter variables (project name, base directory, output directory etc.) within your copy of the script and run. 
# This script runs an installation of ABySS-1.3.4 on beocat. Step 1 is the assembly of multiple single k-mer assemblies. Step 2 is the merging of these single k-mer assemblies by running ABySS on reads and the unitigs from the multi-kmer assemblies. Step 3 writes a script to calculate N25,N50,N75, cumulative scaffold length, number of scaffolds for all assemblies.
#
#  Created by jennifer shelton on 6/29/13.
##################  define variables #################################################

#STEP 1: EDIT THIS TO POINT TO YOUR READS FOR THE ORIGINAL SINGLE K ASSEMBLIES, USE FULL PATH TO THE READS. 
#STEP 2: EDIT THIS TO POINT TO YOUR READS AND THE UNITIGS (FROM YOUR SINGLE K ASSEMBLIES) USE FULL PATH TO THE READS AND UNITIGS. THIS MERGE STRATEGY IS BASED ON THE FOLLOWING POST BY SHAUN JACKMAN "Another approach that I’ve used before and works reasonably well is to run the individual k assemblies to the unitigs stage, then reassemble the reads plus the unitigs from the first assemblies at a larger value of k." https://groups.google.com/d/msg/abyss-users/RXIbiucgmPs/VotHOWKcjhMJ
 
#Base directory where reads are located, must end in a / charachter
my $base_dir='/homes/bjsco/abyss_test/';

#file names in the above directory, must be paired end
my $read_1='test_R1.fastq';
my $read_2='test_R2.fastq';
my $unpaired_reads='test_R1.fastq test_R2.fastq test_R1.fastq'; # SPACE SEPARATED LIST OF YOU SINGLE END READS (READS MUST ALSO BE IN ${base_dir})
my $name='test'; # YOUR PROJECT NAME

my $n=64; #Number of processors to use
my $mem_per_core=8; #Gigs of memory to use


my @kmers = (21, 31, 41, 51, 53, 55, 57, 59); #Kmer values for step1a (these kmers require a higher value for ${n}. ${n}=64 worked for me.)
#my @kmers = (61, 63, 65, 67, 69, 71, 81, 91); #Kmer values for step1b (these kmers require a lower value for ${n}. ${n}=32 worked for me.)
#my @kmers = (61, 71, 81, 91); #Kmer values for step2 (comment out line 31 and uncomment line 32 to run step 2) (these kmers require a lower value for ${n}. ${n}=32 worked for me.)

#Write shell scripts and qsub them
foreach my $k (@kmers)
{
  my $outdir="${base_dir}${name}-${k}";
  my $script="run_abyss_k${k}.sh";
  open(SHELL, ">$script");
  say SHELL '#!/bin/sh';
  say SHELL "export PATH=\$(find /homes/bjsco/abyss-1.3.4 -type d | tr '\n' ':' | sed 's/:\$//'):\${PATH}\n";
  say SHELL "/homes/bjsco/local/bin/abyss-pe name=${name}-${k} k=${k} np=\$NSLOTS lib=\'pe1 pe2\' pe1=\'${base_dir}${read_1}\' pe2=\'${base_dir}${read_2}\' se=\'${unpaired_reads}\' -C ${outdir}";
  close(SHELL);
  `chmod +x ./${script}`;
  `mkdir ${outdir}`;
  `qsub -l h_rt=48:00:00,mem=${mem_per_core}G -pe single ${n} ./${script}`;
}
#STEP 3: RUN ON ASSEMBLIES TO GENERATE ASSEMBLY STATS ON YOUR SCAFFOLDS. 
my $assembly_quality_stats_script="abyss_assemblies_quality_stats.sh";
open(ASSEMBLY_STATS_SHELL, ">>$assembly_quality_stats_script");
say ASSEMBLY_STATS_SHELL '#!/bin/sh';
say ASSEMBLY_STATS_SHELL "perl /homes/bioinfo/bioinfo_software/github_scripts/transcriptome-and-genome-assembly/KSU_bioinfo_lab/assembly_quality_stats_for_multiple_assemblies.pl ${base_dir}${name}-*/${name}-*-scaffolds.fa";
close(ASSEMBLY_STATS_SHELL);