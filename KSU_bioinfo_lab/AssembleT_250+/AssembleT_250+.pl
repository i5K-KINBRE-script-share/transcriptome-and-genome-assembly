#!/usr/bin/perl
###############################################################################
#   
#	USAGE: perl AssembleT_250+.pl [options]
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
use List::Util qw(max);
use List::Util qw(sum);
# use Bio::SeqIO;
# use Bio::Seq;
# use Bio::DB::Fasta;
use Term::ANSIColor;
use Getopt::Long;
use Pod::Usage;
###############################################################################
##############         Print informative message             ##################
###############################################################################
print "########################################################################\n";
print "#  AssembleT_250+.pl Version 1.0                                       #\n";
print colored ("#                 !!!!!WARNING UNDER DEVELOPMENT!!!!                   #", 'bold white on_blue'), "\n";
print "#                                                                      #\n";
print "#  Created by Jennifer Shelton 01/20/15                                #\n";
print "# github.com/i5K-KINBRE-script-share/transcriptome-and-genome-assembly #\n";
print "#  perl AssembleT_250+.pl -help # for usage/options                    #\n";
print "#  perl AssembleT_250+.pl -man # for more details                      #\n";
print "########################################################################\n";
###############################################################################
##############                get arguments                  ##################
###############################################################################
my ($r_list,$clean_read_file1,$clean_read_file2,$clean_read_singletons,$text_out);
my $project_name = "my_project";
my $convert_header = 0;
my $min_read_length = 90;
my $count = 0; # count reads in files
my $max_threads=16; #default max threads
my $man = 0;
my $help = 0;
GetOptions (
			  'help|?' => \$help,
			  'man' => \$man,
              'r|r_list:s' => \$r_list,
              'p|project_name:s' => \$project_name,
              'c|convert_header' => \$convert_header,
              'n|min_read_length:i' => \$min_read_length,
              'j|max_threads:i' => \$max_threads
              )
or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;
my $dirname = dirname(__FILE__); # github directories (all github directories must be in the same directory) no trailing slash
sub quote { qq!"$_[0]"! } ## interpolate slurped text
my $home = getcwd; # working directory (this is where output files will be printed)
#HOME = /homes/bioinfo/test_git
mkdir "${home}/${project_name}_scripts";
mkdir "${home}/${project_name}_qsubs";
mkdir "${home}/${project_name}_prinseq";
###############################################################################
##############                Subroutines                    ##################
###############################################################################
sub mean { return @_ ? sum(@_) / @_ : 0 } # calculate mean of array

sub std_dev #calculate standard deviation of the reverse read lengths
{
    my ($sub_mean, @sub_lengths) = @_;
    my $sub_std_dev_sum = 0;
    $sub_std_dev_sum += ($_ - $sub_mean) ** 2 for @sub_lengths;
    return @sub_lengths ? sqrt($sub_std_dev_sum / @sub_lengths) : 0;
}

###############################################################################
############## Create array of the sample names and read files    #############
###############################################################################
my @reads;
open (READ_LIST, '<', $r_list) or die "Can't open $r_list!\n";
while (<READ_LIST>)
{
    chomp;
    if (/\S/)
    {
        push @reads , [split (/\t/)]; #Unless the line is blank
    }
}
#######################################################################
#########          Prep Masurca assembly scripts             ##########
#######################################################################
open (SCRIPT_MAS, '>', "${home}/${project_name}_scripts/${project_name}_Masurca_Mira.sh") or die "Can't open ${home}/${project_name}_scripts/${project_name}_Masurca_Mira.sh!\n"; # create a shell script for super read assembly with masurca
print SCRIPT_MAS "#!/bin/bash\n";
open (SCRIPT_MAS_CONFIG, '>', "${home}/${project_name}_scripts/${project_name}_make_superead_config.txt") or die "Can't open ${home}/${project_name}_scripts/${project_name}_make_superead_config.txt!\n"; # create a shell script for super read assembly with masurca
print SCRIPT_MAS_CONFIG "DATA\n";
###############################################################################
##############     Write scripts for each sample             ##################
###############################################################################
my $singleton_count = 1;
for my $samples (@reads)
{
    my  @seq_lengths; ## INITIALIZE ARRAY FOR LENGTH AND STD CALCULATION FOR READS
    my $sample_code = $samples->[0]; # get two character sample code
    if ($sample_code =~ /^S/)
    {
        print "Warning, your library prefix may not start with a \"S\".\n";
        exit;
    }
    my @r1 = split(',',$samples->[1]); # get list of forward reads
    my @r2 = split(',',$samples->[2]); # get list of reverse reads
    if (scalar(@r1) != scalar(@r2))
    {
        print "Error the number of forward and reverse read files does not match!\n"; # each forward file must have a corresponding reverse file
        exit;
    }
    #######################################################################
    ############ Convert headers of illumina paired-end data ##############
    #######################################################################
    open (QSUBS_CLEAN, '>', "${home}/${project_name}_qsubs/${project_name}_qsubs_clean.sh") or die "Can't open ${home}/${project_name}_qsubs/${project_name}_qsubs_clean.sh!\n";
    print QSUBS_CLEAN '#!/bin/bash';
    print QSUBS_CLEAN "\n";
    for my $file (0..$#r1)
    {
        #######################################################################
        ######                 Check forward file paths               #########
        #######################################################################
        if (($r1[$file] =~ /^\s/) || ($r1[$file] =~ /\s$/))
        {
            $r1[$file] =~ s/ //g; ## removed white space because bioperl doesn't allow it in filenames
        }
        elsif ($r1[$file] =~ /\s/)
        {
            print "Error filename $r1[$file] includes spaces. Remove spaces from your filenames, update your readlist $r_list accordingly and re-run command\n";
            exit;
        }
        #######################################################################
        ######                 Check reverse file paths               #########
        #######################################################################
        if (($r2[$file] =~ /^\s/) || ($r2[$file] =~ /\s$/))
        {
            $r2[$file] =~ s/ //g; ## removed white space because bioperl doesn't allow it in filenames
        }
        elsif ($r2[$file] =~ /\s/)
        {
            print "Error filename $r2[$file] includes spaces. Remove spaces from your filenames, update your readlist $r_list accordingly and re-run command\n";
            exit;
        }
        #######################################################################
        ###### Split read filenames into usefull parts for renaming   #########
        # and avoiding relative paths (some software disliked relative paths) #
        #######################################################################
        my (${filename}, ${directories}, ${suffix}) = fileparse($r1[$file],'\..*'); # break appart filenames
        my (${filename2}, ${directories2}, ${suffix2}) = fileparse($r2[$file],'\..*'); # break appart filenames
        open (SCRIPT, '>', "${home}/${project_name}_scripts/${filename}_clean.sh") or die "Can't open ${home}/${project_name}_scripts/${filename}_clean.sh!\n"; # create a shell script for each read-pair set
        print SCRIPT "#!/bin/bash\n";
        if ($convert_header)
        {
            print SCRIPT "#######################################################################\n############ Convert headers of illumina paired-end data ##############\n#######################################################################\n";
            print SCRIPT "cat $r1[$file] | awk \'{if (NR % 4 == 1) {split(\$1, arr, \":\"); printf \"%s_%s:%s:%s:%s:%s#0/%s\\n\", arr[1], arr[3], arr[4], arr[5], arr[6], arr[7], substr(\$2, 1, 1), \$0} else if (NR % 4 == 3){print \"+\"} else {print \$0} }\' > ${home}/${filename}_header.fastq\n"; # Convert headers for R1
            $r1[$file] = "${home}/${filename}_header.fastq"; # redefine R1
            print SCRIPT "cat $r2[$file] | awk \'{if (NR % 4 == 1) {split(\$1, arr, \":\"); printf \"%s_%s:%s:%s:%s:%s#0/%s\\n\", arr[1], arr[3], arr[4], arr[5], arr[6], arr[7], substr(\$2, 1, 1), \$0} else if (NR % 4 == 3){print \"+\"} else {print \$0} }\' > ${home}/${filename2}_header.fastq\n"; # Convert headers for R2
            $r2[$file] = "${home}/${filename2}_header.fastq"; #redefine R2
        }
        #######################################################################
        ######     Get read size of R2 library to pass to Masurca     #########
        #######################################################################
        open (my $test_read,'<',"${home}/${filename}${suffix2}") or die "can't open ${home}/${filename}${suffix2}. You must use absolute paths in the read list file \"-r\" or cd to the directory with you reads before you call this script!\n";
        ## ADD LENGTH AND STD CALCULATION FOR READS
        my $first_line = 1;
        while( <$test_read> )
        {
            if ($first_line)
            {
                $first_line = 0; #skip first line
                next;
            }
            chomp;
            my @seq = split ('');
            push (@seq_lengths,scalar(@seq)); # grab length of sequence line and store in a list
            foreach (1..3) # skip the rest of the lines
            {
                my $next_lines .= <$test_read>;
            }
        }
#        1 while( <TEST_READ> );
#        my $count = ($. + $count); ##Get number of lines in a file
        #######################################################################
        ######### Clean reads for low quality without de-duplicating ##########
        #######################################################################
        print QSUBS_CLEAN "qsub -l h_rt=24:00:00,mem=10G ${home}/${project_name}_scripts/${filename}_clean.sh\n";
        $text_out = read_file("${dirname}/Prinseq_template.txt"); ## read shell template with slurp
        print SCRIPT eval quote($text_out);
        print SCRIPT "\n";
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
#    my  @seq_lengths; ## INITIALIZE ARRAY FOR LENGTH AND STD CALCULATION FOR READS
#    my $sample_code = $samples->[0]; # get two character sample code
    print SCRIPT_MAS "cat$clean_read_file1 > ${home}/${project_name}_${sample_code}_good_1.fastq # concatenate fastq\n";
    print SCRIPT_MAS "cat$clean_read_file2 > ${home}/${project_name}_${sample_code}_good_2.fastq # concatenate fastq\n";
    print SCRIPT_MAS "cat$clean_read_singletons > ${home}/${project_name}_${sample_code}_good_singletons.fastq # concatenate single fastq\n";
    #######################################################################
    #########            Write Masurca Config file               ##########
    #######################################################################
    my $read_mean = &mean(@seq_lengths); #find mean read length for a file
    my $read_standard_deviation = &std_dev($read_mean,@seq_lengths); #find standard deviation of read length for a file
    print SCRIPT_MAS_CONFIG "PE= $sample_code $read_mean $read_standard_deviation ${home}/${project_name}_${sample_code}_good_1.fastq ${home}/${project_name}_${sample_code}_good_2.fastq\n"; # add to Masurca Config file
    print SCRIPT_MAS_CONFIG "PE= S${singleton_count} $read_mean $read_standard_deviation ${home}/${project_name}_${sample_code}_good_singletons.fastq\n"; # add to Masurca Config file
    ++$singleton_count;
    if ($singleton_count >= 10) ## Warn about auto generated library prefix
    {
        print "Warning, your scripts have more than 9 singleton libraries. You will need to make the library prefix unique and only two letters manually in the script ${home}/${project_name}_scripts/${project_name}_make_superead.sh before running. The two letter prefix for a library must be unique for Masurca (see manual ftp://ftp.genome.umd.edu/pub/MaSuRCA/MaSuRCA_QuickStartGuide.pdf)."
    }
    $clean_read_file1 = '';
    $clean_read_file2 =  '';
    $clean_read_singletons = '';
}
#######################################################################
## Write parameter section of masurca_config_parameter_template.txt  ##
#######################################################################
$text_out = read_file("${dirname}/masurca_config_parameter_template.txt"); ## read shell template with slurp
print SCRIPT_MAS_CONFIG eval quote($text_out);
print SCRIPT_MAS_CONFIG "\n";

#######################################################################
########               Write Masurca and Mira script          #########
#######################################################################
$text_out = read_file("${dirname}/Masurca_Mira_template.txt"); ## read shell template with slurp
print SCRIPT_MAS eval quote($text_out);
print SCRIPT_MAS "\n";
open (QSUBS_MAS_MIRA, '>', "${home}/${project_name}_qsubs/${project_name}_qsubs_Masurca_Mira.sh") or die "Can't open ${home}/${project_name}_qsubs/${project_name}_qsubs_Masurca_Mira.sh!\n";
print QSUBS_MAS_MIRA "#!/bin/bash\n";
print QSUBS_MAS_MIRA "qsub -l h_rt=100:00:00,mem=4G -pe single ${max_threads} ${home}/${project_name}_scripts/${project_name}_Masurca_Mira.sh\n";


print "Done\n";
##################################################################################
##############                  Documentation                   ##################
##################################################################################
## style adapted from http://www.perlmonks.org/?node_id=489861 
__END__

=head1 SYNOPSIS

AssembleT_250+.pl - The script writes scripts and qsubs to assemble Illumina paired end reads that are 250 bp or longer into a de novo transcriptome. The script 1) converts illumina headers if the "-c" parameter is used, 2) cleans and deduplicates raw reads using Prinseq http://prinseq.sourceforge.net/manual.html, 3) reads are the assembled into supereads wuth Masurca, 5) these supereads are converted to FASTA file format and assembled with Mira.

=head1 USAGE

perl AssembleT_250+.pl [options]

 Documentation options:
   -help    brief help message
   -man	    full documentation
 Required options:
   -r        full path for file with tab separated list of two-letter Masurca library prefix codes, fastq forward and fastq reverse read files
 Recommended options:
   -p	     project name (no spaces)(default = my_project)
 Filtering options:
   -n	     minimum read length (default = 90)
 Fastq format options:
   -c	     convert fastq headers
 Additional options:
   -j	     max number of threads for superead assembly
=head1 OPTIONS

=over 8

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the more detailed manual page with output details and examples and exits.

=item B<-r, --r_list>
 
This is the the full path (path and filename) of the user provided list of read files. The file should be tab separated with the two-letter Masurca library prefix codes, then the first read file, then the second read file. Example:
 
D1  sample_data/sample_1_R1.fastq   sample_data/sample_1_R2.fastq
 
If a sample has multiple fastq files for R1 and R2 separate these with commas. Example:
 
D2  sample_data/sample_1a_R1.fastq,sample_data/sample_1b_R1.fastq,sample_data/sample_1c_R1.fastq   sample_data/sample_1a_R2.fastq,sample_data/sample_1b_R2.fastq,sample_data/sample_1c_R2.fastq
 
=item B<-j, --max_threads>
 
The max number of threads for superead assembly on Beocat (default=16).
 
=item B<-n, --min_read_length>
 
The minimum read length. Reads shorter than this after cleaning will be discarded. Default minimum length is 90bp.

=item B<-c, --convert_header>
 
If the illumina headers do not end in /1 or /2 use this parameter to indicat that headers need to be converted. Check your headers by typing "head [fasta full path]" and read more about illumina headers at http://en.wikipedia.org/wiki/Fastq#Illumina_sequence_identifiers.
 
 
=back

=head1 DESCRIPTION

B<OUTPUT DETAILS:>

see: https://github.com/i5K-KINBRE-script-share/transcriptome-and-genome-assembly/blob/master/KSU_bioinfo_lab/AssembleT/AssembleT_LAB.md

B<Test with sample datasets:>
 
# Find a more detailed instructions at https://github.com/i5K-KINBRE-script-share/transcriptome-and-genome-assembly/blob/master/KSU_bioinfo_lab/AssembleT/AssembleT_LAB.md
 
# Clone the Git repositories
 
git clone https://github.com/i5K-KINBRE-script-share/transcriptome-and-genome-assembly

# Make a working directory.
 
mkdir de_novo_transcriptome
cd de_novo_transcriptome

# Create symbolic links to subsampled raw RNA reads from the human breast cancer cell lines.
 
ln -s /homes/bioinfo/pipeline_datasets/AssembleT/* ~/de_novo_transcriptome/
 
# Write assembly scripts

perl ~/transcriptome-and-genome-assembly/KSU_bioinfo_lab/AssembleT/AssembleT_250+.pl -r cell_line_reads_assembly.txt -p cell_line
 
# Clean raw reads. When these jobs are complete go to next step. Test completion by typing "status" in a Beocat session.
 
bash ~/de_novo_transcriptome/cell_line_qsubs/cell_line_qsubs_clean.sh
 
# Concatenate cleaned reads and Assemble transcriptomes. Test completion by typing "status" in a Beocat session.
 
 bash ~/de_novo_transcriptome/cell_line_qsubs/testing_qsubs_Masurca_Mira.sh
 

=cut
