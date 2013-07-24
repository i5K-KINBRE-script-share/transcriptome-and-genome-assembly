#!/bin/perl

#  assembly_quality_stats_for_multiple_assemblies.pl
#  
#
#  Created by jennifer shelton on 7/24/13.
# This script runs a slightly modified version of Joseph Fass' Count_fasta.pl (original available at http://wiki.bioinformatics.ucdavis.edu/index.php/Count_fasta.pl ) on a fasta file from each assembly. It then creates comma separated file called assembly_metrics.csv listing the N25,N50,N75, cumulative contig length, and number of contigs for each assembly. To use also download Count_fastas.pl https://github.com/i5K-KINBRE-script-share/transcriptome-and-genome-assembly/blob/master/KSU_bioinfo_lab/Count_fastas.pl and change $path_to_Count_fastas on line 13 of assembly_quality_stats_for_multiple_assemblies.pl.
# usage: assembly_quality_stats_for_multiple_assemblies.pl [fasta_file or files]

use strict;
use warnings;
my $out;
my $path_to_Count_fastas="/homes/bioinfo/bioinfo_software/github_scripts/transcriptome-and-genome-assembly/KSU_bioinfo_lab";
my $outfile="assembly_metrics.csv";
open (METRICS, ">$outfile");
print METRICS "Filename,Cumulative length of contigs(bp),Number of contigs,N25(bp),N50(bp),N75(bp),Number of contigs longer than N25,Number of contigs longer than N50,Number of contigs longer than N75\n";
foreach my $f (@ARGV)
{
		my $out=`perl ${path_to_Count_fastas}/Count_fastas.pl ${f}`;
		print METRICS "${f},${out}\n";

}
close (METRICS);

