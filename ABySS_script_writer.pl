#!/usr/local/bin/perl
use strict;
use warnings;

# E R Hanschen
# This script writes a ABySS shell script

my $pe_file_path = $ARGV[0];
my $single_file_path = $ARGV[1];
my $name = $ARGV[2];
my $kmer = $ARGV[3];
my $contig_pairs = $ARGV[4];
my $scaffold_pairs = $ARGV[5];
my $num_threads = $ARGV[6];

if ( !defined $pe_file_path || !defined $single_file_path || !defined $name || !defined $kmer ) {
	die "Necessary input: perl $0 pe_file_path single_file_path name kmer\n";
};

if ( !defined $contig_pairs ) { $contig_pairs = 10; };
if ( !defined $scaffold_pairs ) { $scaffold_pairs = 10; };
if ( !defined $num_threads ) { $num_threads = 64; };\

open (SCRIPT, ">$name.sh") || die "Couldn't open $name: $!\n";

print SCRIPT "#!/bin/sh\n";

print SCRIPT "/homes/bjsco/local/bin/abyss-pe name=$name k=$kmer n=$contig_pairs N=$scaffold_pairs np=\$NSLOTS lib='pe1 pe2' pe1='$pe_file_path/Iloxense5_rep1_ACAGTG_L005_R1_001.fastq.gz  $pe_file_path/Iloxense5_rep1_ACAGTG_L005_R2_001.fastq.gz' pe2='$pe_file_path/Iloxense5_rep2_GCCAAT_L005_R1_001.fastq.gz $pe_file_path/Iloxense5_rep2_GCCAAT_L005_R2_001.fastq.gz' se='$single_file_path/reads.s1r2.fna $single_file_path/reads.s1r1.fna'\n";

close SCRIPT;
