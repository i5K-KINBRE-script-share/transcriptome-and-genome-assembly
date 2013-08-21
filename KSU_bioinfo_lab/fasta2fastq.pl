#!/usr/bin/env perl

##### Script converts FAST and QAUL files to FASTQ format. It is a slight modification of a script posted to SeqAnwsers http://seqanswers.com/forums/showthread.php?t=2775 
##### Usage: fasta2fastq [fasta.file] [qual.file] 

use strict;
use warnings;

my $offset = 33; # I think this was 33 for sanger FASTQ, change this if required!
my $count = 0;
my %seqs;
my $fastq = '/homes/bioinfo/Tcas/sanger_reads/tribolium_castaneum.001_R2.fastq';
die ("Usage: fasta2fastq [fasta.file] [qual.file]") unless  (scalar @ARGV) == 2;

open FASTA, $ARGV[0] or die "cannot open fasta: $!\n";
open QUAL, $ARGV[1] or die "cannot open qual: $!\n";
open FASTQ, ">$fastq";

$/ = ">";
while (<FASTA>)
{
	unless ($count == 0)
	{
		chomp;
		my ($fdef, @seqLines) = split /\n/;
		my $seq = join '', @seqLines;
		$seqs{$fdef} = $seq;
	}
	$count++;
}
close FASTA;
$count = 0;
my $qual0 =33;
while (<QUAL>)
{
	unless ($count == 0)
	{
		chomp;
		my ($qdef, @qualLines) = split /\n/;
		my $qualString = join ' ', @qualLines;
		$qualString =~ s/\s+/ /g;
		my @quals = split / /, $qualString;
		print FASTQ "@","$qdef\n";
		print FASTQ "$seqs{$qdef}\n";
		print FASTQ "+\n";
		foreach my $qual (@quals)
		{
            unless ($qual eq ''){print FASTQ chr($qual + $qual0)};
		}
		print FASTQ "\n";
	}
	$count++;
}

close QUAL;