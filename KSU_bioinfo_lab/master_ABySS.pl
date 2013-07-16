#!/usr/bin/perl
use strict;
use warnings;

# E R Hanschen
# This script is the master runner for making ABySS shell scripts

my $pe_file_path = "/homes/bjsco/iochromas/Project_SmithS_01_SOL/Sample_Iloxense_combined";
my $single_file_path = "/homes/bjsco/iochromas/read";
my $short_name = "Kmer_opt";
my $contig_pairs = "";
my $scaffold_pairs = "";
my $num_threads = "80";

my @number_array = (21, 31, 41, 51, 61, 71, 81, 91);
#my @number_array = (53, 55, 57, 59, 63, 65, 67, 69);
foreach my $i (@number_array) {
	my $long_name = "$short_name$i";
	&WritingScript($long_name, $i);
};

sub WritingScript {

my $name = $_[0];
my $j = $_[1];

`perl ABySS_script_writer.pl $pe_file_path $single_file_path $name $j $contig_pairs $scaffold_pairs $num_threads`;
`qsub -l h_rt=72:00:00,mem=12G -pe single $num_threads $name.sh`;

};
