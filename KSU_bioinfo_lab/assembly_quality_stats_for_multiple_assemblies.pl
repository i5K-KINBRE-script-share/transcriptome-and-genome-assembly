#!/bin/perl

#  assembly_quality_stats_for_multiple_assemblies.pl
#  
#
#  Created by jennifer shelton on 7/24/13.
# This script runs a slightly modified version of Joseph Fass' Count_fasta.pl (original available at http://wiki.bioinformatics.ucdavis.edu/index.php/Count_fasta.pl ) on a fasta file from each assembly. It then creates comma separated file called assembly_metrics.txt listing the N25,N50,N75, cumulative contig length, and number of contigs for each assembly.
# usage: assembly_quality_stats_for_multiple_assemblies.pl [fasta_file or files]

use strict;
use warnings;
my $out;
my $path_to_Count_fastas="/Users/jennifershelton/Desktop/Perl_course_texts";
my $outfile="assembly_metrics.txt";
open (METRICS, ">$outfile");
print METRICS "Filename,Cumulative length of contigs(bp),Number of contigs,N25(bp),N50(bp),N75(bp),Number of contigs longer than N25,Number of contigs longer than N50,Number of contigs longer than N75\n";
foreach my $f (@ARGV)
{
		my $out=`perl ${path_to_Count_fastas}/Count_fastas.pl $f`;
		print METRICS "${f},${out}\n";

}
close (METRICS);

