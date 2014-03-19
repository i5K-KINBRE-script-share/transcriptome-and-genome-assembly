#!/bin/perl
###############################################################################
#   
#	USAGE: perl transcriptome_assembly_pipeline.pl [options]
#
#  Created by jennifer shelton
#
###############################################################################
use strict;
use warnings;
# use File::Basename; # enable manipulating of the full path
# use List::Util qw(max);
# use List::Util qw(sum);
# use Bio::SeqIO;
# use Bio::Seq;
# use Bio::DB::Fasta;
use Getopt::Long;
use Pod::Usage;
###############################################################################
##############         Print informative message             ##################
###############################################################################
print "########################################################################\n";
print "#  transcriptome_assembly_pipeline.pl                                  #\n";
print "#                                                                      #\n";
print "#  Created by Jennifer Shelton 03/19/14                                #\n";
print "#  github.com/i5K-KINBRE-script-share/genome-annotation-and-comparison #\n";
print "#  perl transcriptome_assembly_pipeline.pl -help # for usage/options   #\n";
print "#  perl transcriptome_assembly_pipeline.pl -man # for more details     #\n";
print "########################################################################\n";
###############################################################################
##############                get arguments                  ##################
###############################################################################
my ();

my $man = 0;
my $help = 0;
GetOptions (
			  'help|?' => \$help, 
			  'man' => \$man,
			  'r|r_cmap:s' => \$r_cmap,    
              's_algn|sa:f' => \$second_min_per_aligned  
              )  
or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;
# my $dirname = dirname(__FILE__);
##################################################################################
##############                        run                       ##################
##################################################################################

print "done\n";
##################################################################################
##############                  Documentation                   ##################
##################################################################################
## style adapted from http://www.perlmonks.org/?node_id=489861 
__END__

=head1 NAME

transcriptome_assembly_pipeline.pl - a package of scripts that ...

=head1 USAGE

perl transcriptome_assembly_pipeline.pl [options]

 Documentation options:
   -help    brief help message
   -man	    full documentation
 Required options:
   -r	     reference CMAP
 Filtering options:
   --s_algn	 second minimum % of possible alignment   
   
=head1 OPTIONS

=over 8

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the more detailed manual page with output details and examples and exits.

=item B<-r, --r_cmap>

The reference CMAP produced by IrysView when you create an XMAP. It can be found in the "Imports" folder within a workspace.

=item B<--f_algn, --sa>

The minimum percent of the full potential length of the alignment allowed for the second round of filtering. This should be higher than the setting for the first round of filtering.

=back

=head1 DESCRIPTION

B<OUTPUT DETAILS:>

This appears when the manual is viewed!!!!

B<Test with sample datasets:>


perl transcriptome_assembly_pipeline.pl -r sample_data/sample.r.cmap --s_algn .9

=cut