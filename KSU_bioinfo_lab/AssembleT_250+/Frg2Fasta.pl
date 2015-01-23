#!/usr/bin/perl
###############################################################################
#
#	USAGE: perl Frg2Fasta.pl <FRG_file>
#   DESCRIPTION: Script converts Celera assembly FRG format files to FASTA
#
#  Created by jennifer shelton
#
###############################################################################
use strict;
use warnings;
use File::Basename; # enable manipulating of the full path

my $frag_file = $ARGV[0];
open (my $frag, "<", $frag_file) or die "Can't open $frag_file!\n";
my (${filename}, ${directories}, ${suffix}) = fileparse($ARGV[0],'\..*'); # break appart filenames
my $fasta_file = "${directories}/${filename}.fa";
open (my $fasta, ">", $fasta_file) or die "Can't open $fasta_file!\n";
my $seq_count =1;
while (<$frag>)
{
    if (/^seq:/)
    {
        my $seq .= <$frag>;
        print $fasta ">Superead_${seq_count}\n";
        print $fasta "$seq";
        ++$seq_count;
    }
}
