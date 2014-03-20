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
use File::Basename; # enable manipulating of the full path
use Cwd;
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
my ($r_list,$clean_read_file1,$clean_read_file2,$clean_read_singletons);
my $project_name = "my_project";
my $convert_header = 0;
my $shortest_k = 25,
my $longest_k = 65,
my $increment_k = 2,
my $merge_k = 39;
my $ins_length = 250;
my $count = 0;

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
              'x|ins_length:i' => \$ins_length
              )
or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;
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
##############     Write scripts for each sample             ##################
###############################################################################
for my $samples (@reads)
{
    my @r1 = split(',',$samples->[0]); # get list of forward reads
    my @r2 = split(',',$samples->[1]); # get list of reverse reads
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
    open (QSUBS_MAP, '>', "${home}/${project_name}_qsubs/${project_name}_qsubs_map.sh") or die "Can't open ${home}/${project_name}_qsubs/${project_name}_qsubs_map.sh!\n";
    print QSUBS_MAP '#!/bin/bash';
    print QSUBS_MAP "\n";
    for my $file (0..$#r1)
    {
        open (TEST_READ,'<',"$r1[$file]") or die "can't open $r1[$file]!\n";
        1 while( <TEST_READ> );
        my $count = ($. + $count);
        
        my (${filename}, ${directories}, ${suffix}) = fileparse($r1[$file],'\..*'); # break appart filenames
        my (${filename2}, ${directories2}, ${suffix2}) = fileparse($r2[$file],'\..*'); # break appart filenames
        open (SCRIPT, '>', "${home}/${project_name}_scripts/${filename}_clean.sh") or die "Can't open ${home}/${project_name}_scripts/${filename}_clean.sh!\n"; # create a shell script for each read-pair set
        print SCRIPT '#!/bin/bash';
        print SCRIPT "\n";
        if ($convert_header)
        {
            print SCRIPT "#######################################################################\n############ Convert headers of illumina paired-end data ##############\n#######################################################################\n";
            print SCRIPT "cat $r1[$file] | awk \'{if (NR % 4 == 1) {split(\$1, arr, \":\"); printf \"%s_%s:%s:%s:%s:%s#0/%s\\n\", arr[1], arr[3], arr[4], arr[5], arr[6], arr[7], substr(\$2, 1, 1), \$0} else if (NR % 4 == 3){print \"+\"} else {print \$0} }\' > ${home}/${filename}_header.fastq\n";
            $r1[$file] = "${home}/${filename}_header.fastq";
            print SCRIPT "cat $r2[$file] | awk \'{if (NR % 4 == 1) {split(\$1, arr, \":\"); printf \"%s_%s:%s:%s:%s:%s#0/%s\\n\", arr[1], arr[3], arr[4], arr[5], arr[6], arr[7], substr(\$2, 1, 1), \$0} else if (NR % 4 == 3){print \"+\"} else {print \$0} }\' > ${home}/${filename2}_header.fastq\n";
            $r2[$file] = "${home}/${filename2}_header.fastq";
            
        }
        #######################################################################
        ######### Clean reads for low quality without de-duplicating ##########
        #######################################################################
        print SCRIPT "#######################################################################\n######### Clean reads for low quality without de-duplicating ##########\n#######################################################################\n";
        print QSUBS_CLEAN "qsub -l h_rt=48:00:00,mem=40G ${home}/${project_name}_scripts/${filename}_clean.sh\n";
        print SCRIPT "perl /homes/sheltonj/abjc/prinseq-lite-0.20.3/prinseq-lite.pl -verbose -fastq $r1[$file] -fastq2 $r2[$file] -min_len 90 -min_qual_mean 25 -trim_qual_type mean -trim_qual_rule lt -trim_qual_window 2 -trim_qual_step 1 -trim_qual_left 20 -trim_qual_right 20 -ns_max_p 1 -trim_ns_left 5 -trim_ns_right 5 -lc_method entropy -lc_threshold 70 -out_format 3 -no_qual_header -log ${home}/${project_name}_prinseq/${filename}_paired.log\ -graph_data ${home}/${project_name}_prinseq/${filename}_raw.gd -out_good ${home}/${filename}_good -out_bad ${home}/${filename}_bad\n";
        print SCRIPT "perl /homes/sheltonj/abjc/prinseq-lite-0.20.3/prinseq-lite.pl -verbose -fastq ${home}/${filename}_good_1.fastq -fastq2 ${home}/${filename}_good_2.fastq -out_good null -graph_data ${home}/${project_name}_prinseq/${filename}_cleaned.gd -out_bad null\n";
        print SCRIPT "perl /homes/sheltonj/abjc/prinseq-lite-0.20.3/prinseq-lite.pl -verbose -fastq ${home}/${filename}_good_1_singletons.fastq -out_good null -graph_data ${home}/${project_name}_prinseq/${filename}_cleaned_1_singletons.gd -out_bad null\n";
        print SCRIPT "perl /homes/sheltonj/abjc/prinseq-lite-0.20.3/prinseq-lite.pl -verbose -fastq ${home}/${filename}_good_2_singletons.fastq -out_good null -graph_data ${home}/${project_name}_prinseq/${filename}_cleaned_2_singletons.gd -out_bad null\n";
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
    open (SCRIPT, '>', "${home}/${project_name}_scripts/cat_reads.sh") or die "Can't open ${home}/${project_name}_scripts/cat_reads.sh!\n"; # create a shell script
    print SCRIPT "#!/bin/bash\n";
    print SCRIPT "cat$clean_read_file1 > ${home}/${project_name}_good_1.fastq # concatenate fasta\n";
    print SCRIPT "cat$clean_read_file2 > ${home}/${project_name}_good_2.fastq # concatenate fasta\n";
    print SCRIPT "cat$clean_read_singletons > ${home}/${project_name}_good_singletons.fastq # concatenate single fasta\n";
    
    #######################################################################
    #########         Assemble single k-mer assemblies           ##########
    #######################################################################
    open (QSUBS_SINGLEK, '>', "${home}/${project_name}_qsubs/${project_name}_qsubs_singlek.sh") or die "Can't open ${home}/${project_name}_qsubs/${project_name}_qsubs_map.sh!\n";
    print QSUBS_SINGLEK "#!/bin/bash\n";


    for ( my $k = $shortest_k; $k <= $longest_k; $k += $increment_k )
    {

        close (SCRIPT);
        open (SCRIPT, '>', "${home}/${project_name}_scripts/${project_name}_$k_assemble.sh") or die "Can't open ${home}/${project_name}_scripts/${project_name}_$k_assemble.sh!\n"; # create a shell script for each read-pair set
        print SCRIPT "#!/bin/bash\n";
        print SCRIPT "#######################################################################\n#########         Assemble single k-mer assemblies  k=$k     ##########\n#######################################################################\n";
        print SCRIPT "set -o verbose\n";
        print SCRIPT "PATH=/homes/sheltonj/abjc/velvet_1.2.08:/homes/sheltonj/abjc/oases_0.2.08:\${PATH}\n";
        print SCRIPT "export PATH\n";
        ######### shuffle sequences (if your pairs are unbroken but in two fastq files) ##########
        print SCRIPT "perl /homes/sheltonj/abjc/velvet_1.2.08/contrib/shuffleSequences_fasta/shuffleSequences_fastq.pl ${home}/${project_name}_good_1.fastq ${home}/${project_name}_good_2.fastq ${home}/${project_name}_good_shuff_pairs.fastq\n";
        print SCRIPT "cd ${home}\n";
        print SCRIPT "velveth ${project_name}_${k} ${k} -fastq -short ${home}/${project_name}_good_singletons.fastq -shortPaired -interleaved -fastq ${home}/${project_name}_good_shuff_pairs.fastq\n";
        print SCRIPT "velvetg ${project_name}_${k} -read_trkg yes -ins_length ${ins_length}\n";
        print SCRIPT "oases ${project_name}_${k}\n";
        ######### estimates memory requirements and write qsubs for beocat ###
        my $mem=30;
        my $kmem=(-109635 + 18977*100 + 86326*177 + 233353*$count*3 - 51092*${k});
        $mem=(${kmem}/1000000);
        print QSUBS_SINGLEK "qsub -l h_rt=100:00:00,mem=${mem}G ${home}/${project_name}_scripts/${project_name}_${k}_assemble.sh\n";
#        Ram required for velvetg = -109635 + 18977*ReadSize + 86326*GenomeSize + 233353*NumReads - 51092*K

    }
    #######################################################################
    #########   Assemble merged k-mer assemblies  k=${merge_k}   ##########
    #######################################################################
    close (SCRIPT);
    open (SCRIPT, '>', "${home}/${project_name}_scripts/${project_name}_merge_${merge_k}_assemble.sh") or die "Can't open ${home}/${project_name}_scripts/${project_name}_merge_${merge_k}_assemble.sh!\n"; # create a shell script for each read-pair set
    print SCRIPT "#!/bin/bash\n";
    print SCRIPT "#######################################################################\n#########         Assemble merged k-mer assemblies  k=${merge_k}     ##########\n#######################################################################\n";
    print SCRIPT "set -o verbose\n";
    print SCRIPT "PATH=/homes/sheltonj/abjc/velvet_1.2.08:/homes/sheltonj/abjc/oases_0.2.08:\${PATH}\n";
    print SCRIPT "export PATH\n";
    print SCRIPT "cd ${home}\n";
    print SCRIPT "velveth mergedAssembly ${merge_k} -long ${project_name}_*/transcripts.fa\n";
    print SCRIPT "velvetg mergedAssembly -read_trkg yes -conserveLong yes\n";
    print SCRIPT "oases mergedAssembly -merge yes\n";
    close (SCRIPT);
}

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

=item B<-r, --r_list>
 
The filename of the user provided list of read files. The file should be tab separated with the first read file, then the second read file. Example:
sample_data/sample_1_R1.fastq   sample_data/sample_1_R2.fastq
 
If a sample has multiple fastq files for R1 and R2 separate these with commas. Example:
sample_data/sample_1a_R1.fastq,sample_data/sample_1b_R1.fastq,sample_data/sample_1c_R1.fastq   sample_data/sample_1a_R2.fastq,sample_data/sample_1b_R2.fastq,sample_data/sample_1c_R2.fastq
 
 
 
=back

=head1 DESCRIPTION

B<OUTPUT DETAILS:>

This appears when the manual is viewed!!!!

B<Test with sample datasets:>


perl transcriptome_assembly_pipeline.pl -r sample_data/sample.r.cmap --s_algn .9

=cut