#!/usr/bin/perl
###############################################################################
#   
#	USAGE: perl AssembleG.pl [options]
#
#  Created by jennifer shelton
#
###############################################################################
use strict;
use warnings;
use File::Basename; # enable manipulating of the full path
use Cwd;
use lib '/homes/bioinfo/bioinfo_software/perl_modules/File-Slurp-9999.19/lib';
use File::Slurp;
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
print "#  AssembleG.pl Version 1.3                                            #\n";
print "#                                                                      #\n";
print "#  Created by Jennifer Shelton 03/21/14                                #\n";
print "# github.com/i5K-KINBRE-script-share/transcriptome-and-genome-assembly #\n";
print "#  perl AssembleG.pl -help # for usage/options                         #\n";
print "#  perl AssembleG.pl -man # for more details                           #\n";
print "########################################################################\n";
###############################################################################
##############                get arguments                  ##################
###############################################################################
my ($r_list,$clean_read_file1,$clean_read_file2,$clean_read_singletons,$lib_name,$new_nodes,$text_out);
my $project_name = "my_project";
my $convert_header = 0;
my $shortest_k = 21; # must be odd
my $longest_k = 91; # must be odd
my $increment_k = 10; # must be even
my $merge_k = 61; # must be odd
my $min_read_length = 93;
my $count = 0; # count reads in files
my $nodes=32; #Number of processors to use
my $mem_per_core=8; #Gigs of memory to use
my $lib_count = 1;
my $libx_code = ''; #list of libraries and fastq files
my $lib_code = "lib=\'"; #list of libraries
my $man = 0;
my $help = 0;
GetOptions (
			  'help|?' => \$help,
			  'man' => \$man,
              'r|r_list:s' => \$r_list,
              'p|project_name:s' => \$project_name,
              'c|convert_header' => \$convert_header,
              's|shortest_k:i' => \$shortest_k,
              'l|longest_k:i' => \$longest_k,
              'i|increment_k:i' => \$increment_k,
              'm|merge_k:i' => \$merge_k,
              'n|min_read_length:i' => \$min_read_length,
              'nodes:i' => \$nodes,
              'mem_per_core:i' => \$mem_per_core
              )
or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;
sub quote { qq!"$_[0]"! } ## interpolate slurped text
my $dirname = dirname(__FILE__); # github directories (all github directories must be in the same directory) no trailing slash
my $home = getcwd; # working directory (this is where output files will be printed)
#HOME = /homes/bioinfo/test_git
mkdir "${home}/${project_name}_scripts";
mkdir "${home}/${project_name}_qsubs";
mkdir "${home}/${project_name}_prinseq";
###############################################################################
############## Create array of the sample names and read files    #############
###############################################################################
my @reads;
open (READ_LIST, '<', $r_list) or die "Can't open $r_list!\n";
while (<READ_LIST>)
{
    chomp;
    push @reads , [split];
}
###############################################################################
##############                Open QSUB scripts                   #############
###############################################################################
open (QSUBS_CLEAN, '>', "${home}/${project_name}_qsubs/${project_name}_qsubs_clean.sh") or die "Can't open ${home}/${project_name}_qsubs/${project_name}_qsubs_clean.sh!\n";
print QSUBS_CLEAN "#!/bin/bash\n";
###############################################################################
##############     Write scripts for each sample             ##################
###############################################################################
for my $samples (@reads)
{
    $clean_read_file1 = ''; # initialize R1 and R2 for each library
    $clean_read_file2 = '';
    my @r1 = split(',',$samples->[1]); # get list of forward reads
    my @r2 = split(',',$samples->[2]); # get list of reverse reads
    #######################################################################
    ######### Check that read fastq files each have a matching pair #######
    #######################################################################
    if (scalar(@r1) != scalar(@r2))
    {
        print "Error the number of forward and reverse read files does not match!\n"; # each forward file must have a corresponding reverse file
        exit;
    }
    #######################################################################
    #########   Check that library type was correctly specified   #########
    #######################################################################
    unless (($samples->[0] eq "pe") || ($samples->[0] eq "mp"))
    {
        print "Error \"pe\" or \"mp\" are the only library types that this script accepts. \"$samples->[0]\" was used instead in $r_list. See description of the \"-r\" parameter by typing \"perl AssembleG.pl -man\"\n";
        exit;
    }
    #######################################################################
    ############ Convert headers of illumina paired-end data ##############
    #######################################################################
    for my $file (0..$#r1)
    {
        #######################################################################
        ###### Split read filenames into useful parts for renaming    #########
        # and avoiding relative paths (some software disliked relative paths) #
        #######################################################################
        my (${filename}, ${directories}, ${suffix}) = fileparse($r1[$file],'\..*'); # break appart filenames
        my (${filename2}, ${directories2}, ${suffix2}) = fileparse($r2[$file],'\..*'); # break appart filenames
        open (SCRIPT, '>', "${home}/${project_name}_scripts/${filename}_clean.sh") or die "Can't open ${home}/${project_name}_scripts/${filename}_clean.sh!\n"; # create a shell script for each read-pair set
        if ($convert_header)
        {
            print SCRIPT "#!/bin/bash\n";
            print SCRIPT "#######################################################################\n############ Convert headers of illumina paired-end data ##############\n#######################################################################\n";
            print SCRIPT "cat $r1[$file] | awk \'{if (NR % 4 == 1) {split(\$1, arr, \":\"); printf \"%s_%s:%s:%s:%s:%s#0/%s\\n\", arr[1], arr[3], arr[4], arr[5], arr[6], arr[7], substr(\$2, 1, 1), \$0} else if (NR % 4 == 3){print \"+\"} else {print \$0} }\' > ${home}/${filename}_header.fastq\n"; # Convert headers for R1
            $r1[$file] = "${home}/${filename}_header.fastq"; # redefine R1
            print SCRIPT "cat $r2[$file] | awk \'{if (NR % 4 == 1) {split(\$1, arr, \":\"); printf \"%s_%s:%s:%s:%s:%s#0/%s\\n\", arr[1], arr[3], arr[4], arr[5], arr[6], arr[7], substr(\$2, 1, 1), \$0} else if (NR % 4 == 3){print \"+\"} else {print \$0} }\' > ${home}/${filename2}_header.fastq\n"; # Convert headers for R2
            $r2[$file] = "${home}/${filename2}_header.fastq"; #redefine R2            
        }
        #######################################################################
        ######### Clean reads for low quality without de-duplicating ##########
        #######################################################################
        print QSUBS_CLEAN "qsub -l h_rt=24:00:00,mem=10G ${home}/${project_name}_scripts/${filename}_clean.sh\n";
        my $text_out = read_file("${dirname}/Prinseq_template.txt"); ## read shell template with slurp
        print SCRIPT eval quote($text_out);
        print SCRIPT "\n";
        #######################################################################
        #########      List reads files to concatinante later        ##########
        #######################################################################
        if ($clean_read_file1)
        {
            $clean_read_file1 = "$clean_read_file1"." ${home}/${filename}_good_1.fastq";
            $clean_read_file2 = "$clean_read_file2"." ${home}/${filename}_good_2.fastq";
            $clean_read_singletons = "$clean_read_singletons". " ${home}/${filename}_good_1_singletons.fastq ${home}/${filename}_good_2_singletons.fastq";
        }
        else
        {
            $clean_read_file1 = " ${home}/${filename}_good_1.fastq";
            $clean_read_file2 = " ${home}/${filename}_good_2.fastq";
            $clean_read_singletons = " ${home}/${filename}_good_1_singletons.fastq ${home}/${filename}_good_2_singletons.fastq";
        }
    }
    close (SCRIPT);
    #######################################################################
    #########            Concantinate clean reads                ##########
    #######################################################################
    open (SCRIPT, '>>', "${home}/${project_name}_scripts/cat_reads.sh") or die "Can't open ${home}/${project_name}_scripts/cat_reads.sh!\n"; # create a shell script
    print SCRIPT "#!/bin/bash\n";
    $lib_name= "$samples->[0]${lib_count}";
    print SCRIPT "cat$clean_read_file1 > ${home}/${project_name}_${lib_name}_good_1.fastq # concatenate fasta\n";
    print SCRIPT "cat$clean_read_file2 > ${home}/${project_name}_${lib_name}_good_2.fastq # concatenate fasta\n";
    #######################################################################
    #########       Write lists of reads for each library        ##########
    #######################################################################
    ++$lib_count;
    $libx_code = "$libx_code"."${lib_name}=\'${home}/${project_name}_${lib_name}_good_1.fastq ${home}/${project_name}_${lib_name}_good_2.fastq\' "; # append to list of pe or mp libraries
    $lib_code = "$lib_code"."${lib_name} ";
}
my $se_lib_code = "se=\'$clean_read_singletons\'"; # write list of se libraries
$lib_code = "$lib_code"."\' "; # close list of library names
#######################################################################
#########         Assemble single k-mer assemblies           ##########
#######################################################################
open (QSUBS_SINGLEK, '>', "${home}/${project_name}_qsubs/${project_name}_qsubs_singlek.sh") or die "Can't open ${home}/${project_name}_qsubs/${project_name}_qsubs_map.sh!\n";
print QSUBS_SINGLEK "#!/bin/bash\n";
#######################################################################
#########       Write scripts for single k-mer assemblies    ##########
#######################################################################
for ( my $k = $shortest_k; $k <= $longest_k; $k += $increment_k )
{
    if ( $k > ($shortest_k + (($longest_k - $shortest_k)/2)))
    {
        $new_nodes = $nodes/2; # Shorter kmers require a higher value for ${nodes}. ${nodes}=64 worked for a 200Mb genome when k = 21 to 59. ${nodes}=32 worked for the same genome when k = 61 to 91.Therefore, we divide $nodes by 2 for the longest kmer values.
    }
    else
    {
        $new_nodes = $nodes;
    }
    #######################################################################
    #########       Adjust number of nodes for longer kmers      ##########
    #######################################################################
    close (SCRIPT);
    open (SCRIPT, '>', "${home}/${project_name}_scripts/${project_name}_${k}_assemble.sh") or die "Can't open ${home}/${project_name}_scripts/${project_name}_${k}_assemble.sh!\n"; # create a shell script for each read-pair set
    print SCRIPT "#!/bin/bash\n";
    print SCRIPT "#######################################################################\n";
    print SCRIPT "#########         Assemble single k-mer assemblies  k=$k     ##########\n";
    print SCRIPT "#######################################################################\n";
    print SCRIPT "set -o verbose\n";
    $text_out = read_file("${dirname}/Abyss_singlek_template.txt"); ## read shell template with slurp
    print SCRIPT eval quote($text_out);
    print SCRIPT "\n";
    print QSUBS_SINGLEK "qsub -l h_rt=48:00:00,mem=${mem_per_core}G -pe single ${new_nodes} ${home}/${project_name}_scripts/${project_name}_${k}_assemble.sh\n";
}
#######################################################################
#########   Assemble merged k-mer assemblies  k=${merge_k}   ##########
#######################################################################
$se_lib_code = "se=\'$clean_read_singletons ${home}/${project_name}_*/${project_name}-*-unitigs.fa\'"; # write list of se libraries
open (QSUBS_MERGE, '>', "${home}/${project_name}_qsubs/${project_name}_qsubs_merge.sh") or die "Can't open ${home}/${project_name}_qsubs/${project_name}_qsubs_merge.sh!\n";
print QSUBS_MERGE "#!/bin/bash\n";
close (SCRIPT);
open (SCRIPT, '>', "${home}/${project_name}_scripts/${project_name}_merge_${merge_k}_assemble.sh") or die "Can't open ${home}/${project_name}_scripts/${project_name}_merge_${merge_k}_assemble.sh!\n"; # create a shell script for each read-pair set
print SCRIPT "#!/bin/bash\n";
print SCRIPT "#######################################################################\n";
print SCRIPT "#########     Assemble merged k-mer assemblies  k=${merge_k}  #########\n";
print SCRIPT "#######################################################################\n";
print SCRIPT "set -o verbose\n";
$text_out = read_file("${dirname}/Abyss_mergek_template.txt"); ## read shell template with slurp
print SCRIPT eval quote($text_out);
print SCRIPT "\n";
print QSUBS_MERGE "qsub -l h_rt=48:00:00,mem=${mem_per_core}G -pe single ${new_nodes} ${home}/${project_name}_scripts/${project_name}_merge_${merge_k}_assemble.sh\n";
close (SCRIPT);
#######################################################################
#########    QC assemblies and summarize cleaning steps      ##########
#######################################################################
open (QC, '>', "${home}/${project_name}_scripts/${project_name}_qc_assemblies.sh") or die "Can't open ${home}/${project_name}_scripts/${project_name}_qc_assemblies.sh!\n";
print QC "#!/bin/bash\n";
print QC "set -o verbose\n";
print QC "cd ${home}\n";
print QC "#######################################################################\n#########    QC assemblies and summarize cleaning steps      ##########\n#######################################################################\n";
print QC "perl ~/read-cleaning-format-conversion/KSU_bioinfo_lab/pre_post_cleaning_metrics.pl ${home}/${project_name}_prinseq/*_paired.log\n";
print QC "perl ~/genome-annotation-and-comparison/KSU_bioinfo_lab/assembly_quality_stats_for_multiple_assemblies.pl ${home}/${project_name}_*/${project_name}-*-scaffolds.fa ${home}/${project_name}_merge_${merge_k}/${project_name}-merge-*-scaffolds.fa\n";
open (QSUBS_QC, '>', "${home}/${project_name}_qsubs/${project_name}_qsubs_qc.sh") or die "Can't open ${home}/${project_name}_qsubs/${project_name}_qsubs_qc.sh!\n";
print QSUBS_QC "#!/bin/bash\n";
print QSUBS_QC "qsub -l h_rt=300:00:00,mem=2G ${home}/${project_name}_scripts/${project_name}_qc_assemblies.sh\n";
print "done\n";
##################################################################################
##############                  Documentation                   ##################
##################################################################################
## style adapted from http://www.perlmonks.org/?node_id=489861 
__END__

=head1 SYNOPSIS

AssembleG.pl - The script writes scripts and qsubs to assemble illumina paired end reads into a de novo genome. The script 1) converts illumina headers if the "-c" parameter is used, 2) cleans and deduplicates raw reads using Prinseq http://prinseq.sourceforge.net/manual.html, 3) reads are the assembled multiple times with a range of values of k, 4) finally assembly metrics are generated for all assemblies and read length and number are summarized before and after cleaning.

=head1 USAGE

perl AssembleG.pl [options]

 Documentation options:
   -help    brief help message
   -man	    full documentation
 Required options:
   -r	     filename for file with tab separated list of fastq forward and reverse read files (one line per library insert size)
 Recommended options:
   -p	     project name (no spaces)(default = my_project)
   -s	     shortest kmer (default = 21)
   -l	     longest kmer (default = 91)
   -i	     kmer increments (default = 10)
   -m	     merge kmer (default = 61)
 Filtering options:
   -n	     minimum read length (default = 93)
 Fastq format options:
   -c	     convert fastq headers
   
=head1 OPTIONS

=over 8

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the more detailed manual page with output details and examples and exits.

=item B<-r, --r_list>
 
The filename of the user provided list of read files. The file should be tab separated. The first column should be the library type (either "pe" for paired end or "mp" for mate paired). Paired end or "pe" libraries are used to assemble sequence. Your library should be labeled "mp" if it is a long distance mate-pair library. "mp" libraries are not used for assembly, only for scaffolding. See https://github.com/bcgsc/abyss#assembling-multiple-libraries.

Note: If you have multiple insert lengths read files from each should be listed on separate lines so that ABySS can estimate insert length separately.
 
The second column is the first read file, then the third column is the second read file. Example:
 
 pe sample_data/sample_300bp_1_R1.fastq   sample_data/sample_300bp_1_R2.fastq
 pe sample_data/sample_500bp_1_R1.fastq   sample_data/sample_500bp_1_R2.fastq
 mp sample_data/sample_3_R1.fastq   sample_data/sample_3_R2.fastq
 
If a library has multiple fastq files for R1 and R2 separate these with commas. Example:
 
 pe sample_1a_R1.fastq,sample_1b_R1.fastq,sample_1c_R1.fastq   sample_1a_R2.fastq,sample_1b_R2.fastq,sample_1c_R2.fastq
 pe sample_2a_R1.fastq,sample_2b_R1.fastq,sample_2c_R1.fastq   sample_2a_R2.fastq,sample_2b_R2.fastq,sample_2c_R2.fastq
 mp sample_3a_R1.fastq,sample_3b_R1.fastq,sample_3c_R1.fastq   sample_3a_R2.fastq,sample_3b_R2.fastq,sample_3c_R2.fastq
 
=item B<-s, --shortest_k>
 
The minimum kmer length for single kmer assemblies. Default minimum kmer is 21bp.
 
=item B<-l, --longtest_k>
 
The maximum kmer length for single kmer assemblies. Default maximum kmer is 91bp.This value must be even.
 
=item B<-i, --increments_k>
 
The length by which the value of k increases for the next single kmer assembly. Default kmer is 10bp. This value must be even.
 
=item B<-m, --merge_k>
 
The kmer length used when merging single kmer assemblies. Default merge kmer is 61bp. This value must be odd.
 
=item B<-n, --min_read_length>
 
The minimum read length. Reads shorter than this after cleaning will be discarded. Default minimum length is 93bp.

=item B<-c, --convert_header>
 
If the illumina headers do not end in /1 or /2 use this parameter to indicat that headers need to be converted. Check your headers by typing "head [fasta filename]" and read more about illumina headers at http://en.wikipedia.org/wiki/Fastq#Illumina_sequence_identifiers.
 
 
=back

=head1 DESCRIPTION

B<OUTPUT DETAILS:>

 See: https://github.com/i5K-KINBRE-script-share/transcriptome-and-genome-assembly/blob/master/KSU_bioinfo_lab/AssembleG/AssembleG_LAB.md

B<Test with sample datasets:>
 
# Find a more detailed instructions at https://github.com/i5K-KINBRE-script-share/transcriptome-and-genome-assembly/blob/master/KSU_bioinfo_lab/AssembleG/AssembleG_LAB.md
 
# Clone the Git repositories
 
git clone https://github.com/i5K-KINBRE-script-share/transcriptome-and-genome-assembly
git clone https://github.com/i5K-KINBRE-script-share/genome-annotation-and-comparison

# Make a working directory.
 
mkdir de_novo_genome
cd de_novo_genome

# Create symbolic links to subsampled raw RNA reads from the human breast cancer cell lines.
 
ln -s /homes/bioinfo/pipeline_datasets/AssembleG/* ~/de_novo_genome/
 
# Write assembly scripts

 perl ~/transcriptome-and-genome-assembly/KSU_bioinfo_lab/AssembleG/AssembleG.pl -r ~/de_novo_genome/S_aureus_reads.txt -p S_aureus -n 35 --nodes 8 --mem_per_core 3 -s 21 -l 45 -i 2 -m 31
 
# Clean raw reads. When these jobs are complete go to next step. Test completion by typing "status" in a Beocat session.
 
bash ~/de_novo_genome/S_aureus_qsubs/S_aureus_qsubs_clean.sh
 
# Concatenate cleaned reads
 
bash ~/de_novo_genome/S_aureus_scripts/cat_reads.sh
 
# Assemble single kmer transcriptomes. When these jobs are complete go to next step. Test completion by typing "status" in a Beocat session.
 
bash ~/de_novo_genome/S_aureus_qsubs/S_aureus_qsubs_singlek.sh
 
# Merge single kmer transcriptomes. When these jobs are complete go to next step. Test completion by typing "status" in a Beocat session.
 
bash ~/de_novo_genome/S_aureus_qsubs/S_aureus_qsubs_merge.sh
 
# This step will generate assembly metrics and summarize the cleaning step results.
 
bash ~/de_novo_genome/S_aureus_qsubs/S_aureus_qsubs_qc.sh

=cut
