#!/bin/perl
###############################################################################
#   
#	USAGE: perl AssembleT.pl [options]
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
print "#  AssembleT.pl                                  #\n";
print "#                                                                      #\n";
print "#  Created by Jennifer Shelton 03/19/14                                #\n";
print "# github.com/i5K-KINBRE-script-share/transcriptome-and-genome-assembly #\n";
print "#  perl AssembleG.pl -help # for usage/options                         #\n";
print "#  perl AssembleG.pl -man # for more details                           #\n";
print "########################################################################\n";
###############################################################################
##############                get arguments                  ##################
###############################################################################
my ($r_list,$clean_read_file1,$clean_read_file2,$clean_read_singletons);
my $project_name = "my_project";
my $convert_header = 0;
my $shortest_k = 25; # must by odd
my $longest_k = 65; # must by odd
my $increment_k = 2; # must by even
my $merge_k = 39; # must by odd
my $min_read_length = 90;
my $count = 0; # count reads in files
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
              'n|min_read_length:i' => \$min_read_length
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
    for my $file (0..$#r1)
    {
        #######################################################################
        ###### Split read filenames into usefull parts for renaming   #########
        # and avoiding relative paths (some software disliked relative paths) #
        #######################################################################
        my (${filename}, ${directories}, ${suffix}) = fileparse($r1[$file],'\..*'); # break appart filenames
        my (${filename2}, ${directories2}, ${suffix2}) = fileparse($r2[$file],'\..*'); # break appart filenames
        open (SCRIPT, '>', "${home}/${project_name}_scripts/${filename}_clean.sh") or die "Can't open ${home}/${project_name}_scripts/${filename}_clean.sh!\n"; # create a shell script for each read-pair set
        print SCRIPT '#!/bin/bash';
        print SCRIPT "\n";
        if ($convert_header)
        {
            print SCRIPT "#######################################################################\n############ Convert headers of illumina paired-end data ##############\n#######################################################################\n";
            print SCRIPT "cat $r1[$file] | awk \'{if (NR % 4 == 1) {split(\$1, arr, \":\"); printf \"%s_%s:%s:%s:%s:%s#0/%s\\n\", arr[1], arr[3], arr[4], arr[5], arr[6], arr[7], substr(\$2, 1, 1), \$0} else if (NR % 4 == 3){print \"+\"} else {print \$0} }\' > ${home}/${filename}_header.fastq\n"; # Convert headers for R1
            $r1[$file] = "${home}/${filename}_header.fastq"; # redefine R1
            print SCRIPT "cat $r2[$file] | awk \'{if (NR % 4 == 1) {split(\$1, arr, \":\"); printf \"%s_%s:%s:%s:%s:%s#0/%s\\n\", arr[1], arr[3], arr[4], arr[5], arr[6], arr[7], substr(\$2, 1, 1), \$0} else if (NR % 4 == 3){print \"+\"} else {print \$0} }\' > ${home}/${filename2}_header.fastq\n"; # Convert headers for R2
            $r2[$file] = "${home}/${filename2}_header.fastq"; #redefine R2
            
        }
        #######################################################################
        ###### Estimate size of R1 library to estimate the mem needed #########
        #######################################################################
        open (TEST_READ,'<',"${home}/${filename}${suffix2}") or die "can't open ${home}/${filename}${suffix2}. You must use absolute paths in the read list file \"-r\" or cd to the directory with you reads before you call this script!\n";
        1 while( <TEST_READ> );
        my $count = ($. + $count);
        #######################################################################
        ######### Clean reads for low quality without de-duplicating ##########
        #######################################################################
        print SCRIPT "#######################################################################\n######### Clean reads for low quality without de-duplicating ##########\n#######################################################################\n";
        print QSUBS_CLEAN "qsub -l h_rt=24:00:00,mem=10G ${home}/${project_name}_scripts/${filename}_clean.sh\n";
        print SCRIPT "perl /homes/sheltonj/abjc/prinseq-lite-0.20.3/prinseq-lite.pl -verbose -fastq $r1[$file] -fastq2 $r2[$file] -min_len ${min_read_length} -min_qual_mean 25 -trim_qual_type mean -trim_qual_rule lt -trim_qual_window 2 -trim_qual_step 1 -derep 1 -trim_qual_left 20 -trim_qual_right 20 -ns_max_p 1 -trim_ns_left 5 -trim_ns_right 5 -lc_method entropy -lc_threshold 70 -out_format 3 -no_qual_header -log ${home}/${project_name}_prinseq/${filename}_paired.log -graph_data ${home}/${project_name}_prinseq/${filename}_raw.gd -out_good ${home}/${filename}_good -out_bad ${home}/${filename}_bad\n"; # run prinseq to filter low quality reads
        print SCRIPT "perl /homes/sheltonj/abjc/prinseq-lite-0.20.3/prinseq-lite.pl -verbose -fastq ${home}/${filename}_good_1.fastq -fastq2 ${home}/${filename}_good_2.fastq -out_good null -graph_data ${home}/${project_name}_prinseq/${filename}_cleaned.gd -out_bad null\n"; # cleaned metrics on paired reads
        print SCRIPT "perl /homes/sheltonj/abjc/prinseq-lite-0.20.3/prinseq-lite.pl -verbose -fastq ${home}/${filename}_good_1_singletons.fastq -out_good null -graph_data ${home}/${project_name}_prinseq/${filename}_cleaned_1_singletons.gd -out_bad null\n"; # cleaned metrics on singletons (reads where only one mate passed the qc)
        print SCRIPT "perl /homes/sheltonj/abjc/prinseq-lite-0.20.3/prinseq-lite.pl -verbose -fastq ${home}/${filename}_good_2_singletons.fastq -out_good null -graph_data ${home}/${project_name}_prinseq/${filename}_cleaned_2_singletons.gd -out_bad null\n"; # cleaned metrics on singletons (reads where only one mate passed the qc)

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
    ######### shuffle sequences (if your pairs are unbroken but in two fastq files) ##########
    print SCRIPT "perl /homes/sheltonj/abjc/velvet_1.2.08/contrib/shuffleSequences_fasta/shuffleSequences_fastq.pl ${home}/${project_name}_good_1.fastq ${home}/${project_name}_good_2.fastq ${home}/${project_name}_good_shuff_pairs.fastq\n";
    
    #######################################################################
    #########         Assemble single k-mer assemblies           ##########
    #######################################################################
    open (QSUBS_SINGLEK, '>', "${home}/${project_name}_qsubs/${project_name}_qsubs_singlek.sh") or die "Can't open ${home}/${project_name}_qsubs/${project_name}_qsubs_map.sh!\n";
    print QSUBS_SINGLEK "#!/bin/bash\n";


    for ( my $k = $shortest_k; $k <= $longest_k; $k += $increment_k )
    {

        close (SCRIPT);
        open (SCRIPT, '>', "${home}/${project_name}_scripts/${project_name}_${k}_assemble.sh") or die "Can't open ${home}/${project_name}_scripts/${project_name}_${k}_assemble.sh!\n"; # create a shell script for each read-pair set
        print SCRIPT "#!/bin/bash\n";
        print SCRIPT "#######################################################################\n#########         Assemble single k-mer assemblies  k=$k     ##########\n#######################################################################\n";
        print SCRIPT "set -o verbose\n";
        print SCRIPT "PATH=/homes/sheltonj/abjc/velvet_1.2.08:/homes/sheltonj/abjc/oases_0.2.08:\${PATH}\n";
        print SCRIPT "export PATH\n";
        print SCRIPT "cd ${home}\n";
        print SCRIPT "velveth ${project_name}_${k} ${k} -fastq -short ${home}/${project_name}_good_singletons.fastq -shortPaired -interleaved -fastq ${home}/${project_name}_good_shuff_pairs.fastq\n";
        print SCRIPT "velvetg ${project_name}_${k} -read_trkg yes\n";
        print SCRIPT "oases ${project_name}_${k}\n";
        ######### estimates memory requirements and write qsubs for beocat ###
        my $mem=30;
        my $kmem=(-109635 + 18977*100 + 86326*400 + 233353*$count*2 - 51092*${k});
        $mem=(${kmem}/1000000);
        print QSUBS_SINGLEK "qsub -l h_rt=100:00:00,mem=${mem}G ${home}/${project_name}_scripts/${project_name}_${k}_assemble.sh\n";
#        Ram required for velvetg = -109635 + 18977*ReadSize + 86326*GenomeSize + 233353*NumReads - 51092*K

    }
    #######################################################################
    #########   Assemble merged k-mer assemblies  k=${merge_k}   ##########
    #######################################################################
    open (QSUBS_MERGE, '>', "${home}/${project_name}_qsubs/${project_name}_qsubs_merge.sh") or die "Can't open ${home}/${project_name}_qsubs/${project_name}_qsubs_merge.sh!\n";
    print QSUBS_MERGE "#!/bin/bash\n";
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
    my $mem=30;
    my $kmem=(-109635 + 18977*100 + 86326*400 + 233353*$count*2 - 51092*${merge_k});
    $mem=(${kmem}/1000000);
    print QSUBS_MERGE "qsub -l h_rt=100:00:00,mem=${mem}G ${home}/${project_name}_scripts/${project_name}_merge_${merge_k}_assemble.sh\n";
    #######################################################################
    #########           Cluster merged assembly with CDH         ##########
    #######################################################################
    open (CLUSTER_QC, '>', "${home}/${project_name}_scripts/${project_name}_cluster_and_qc_assemblies.sh") or die "Can't open ${home}/${project_name}_scripts/${project_name}_cluster_and_qc_assemblies.sh!\n";
    print CLUSTER_QC "#!/bin/bash\n";
    print CLUSTER_QC "set -o verbose\n";
    print CLUSTER_QC "cd ${home}\n";
#    print CLUSTER_QC '#######  non-redudant from http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0056217#s4 The transcripts from three individual assemblies were clustered (CD-HIT v4.5.4 http://www.bioinformatics.org/cd-hit/) [56] in order to generate a comprehensive reference. Sequence identity threshold and alignment coverage (for the shorter sequence) were both set as 80% to generate clusters. Such clustered transcripts were defined as reference transcripts in this work.###########';
    print CLUSTER_QC "\n/homes/sheltonj/abjc/cd-hit-v4.6.1/cd-hit-est -i ${home}/mergedAssembly/transcripts.fa -o ${home}/mergedAssembly/CDH_clustermergedAssembly_${project_name}_${merge_k}.fa -c .80 -aS .80s\n";

    #######################################################################
    #########    QC assemblies and summarize cleaning steps      ##########
    #######################################################################
    print CLUSTER_QC "#######################################################################\n#########    QC assemblies and summarize cleaning steps      ##########\n#######################################################################\n";
    print CLUSTER_QC "perl ~/read-cleaning-format-conversion/KSU_bioinfo_lab/pre_post_cleaning_metrics.pl ${home}/${project_name}_prinseq/*_paired.log\n";
    print CLUSTER_QC "perl ~/genome-annotation-and-comparison/KSU_bioinfo_lab/assembly_quality_stats_for_multiple_assemblies.pl ${project_name}_*/transcripts.fa mergedAssembly/transcripts.fa mergedAssembly/CDH_clustermergedAssembly_${project_name}_${merge_k}.fa\n";
    open (QSUBS_CLUSTER_QC, '>', "${home}/${project_name}_qsubs/${project_name}_qsubs_cluster_and_qc.sh") or die "Can't open ${home}/${project_name}_qsubs/${project_name}_qsubs_cluster_and_qc.sh!\n";
    print QSUBS_CLUSTER_QC "#!/bin/bash\n";
    print QSUBS_CLUSTER_QC "qsub -l h_rt=300:00:00,mem=2G ${home}/${project_name}_scripts/${project_name}_cluster_and_qc_assemblies.sh\n";
}

print "done\n";
##################################################################################
##############                  Documentation                   ##################
##################################################################################
## style adapted from http://www.perlmonks.org/?node_id=489861 
__END__

=head1 SYNOPSIS

AssembleT.pl - The script writes scripts and qsubs to assemble illumina paired end reads into a de novo transcriptome. The script 1) converts illumina headers if the "-c" parameter is used, 2) cleans and deduplicates raw reads using Prinseq http://prinseq.sourceforge.net/manual.html, 3) index the reference genome for mapping, 4) reads are the assembled multiple times with a range of values of k, 5) these assemblies are merged with Oases using a merge kmer value, 6) then the merged assembly is clusted with CDHit to take the longest of similar putative transcripts, 7) finally assembly metrics are generated for all assemblies and read length and number are summarized before and after cleaning.

=head1 USAGE

perl AssembleT.pl [options]

 Documentation options:
   -help    brief help message
   -man	    full documentation
 Recommended options:
   -p	     project name (no spaces)(default = my_project)
   -s	     shortest kmer (default = 25)
   -l	     longest kmer (default = 65)
   -i	     kmer increments (default = 2)
   -m	     merge kmer (default = 39)
 Filtering options:
   -n	     minimum read length (default = 90)
 Fastq format options:
   -c	     convert fastq headers
   
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
 
=item B<-s, --shortest_k>
 
The minimum kmer length for single kmer assemblies. Default minimum kmer is 25bp. This value must by odd.
 
=item B<-l, --longtest_k>
 
The maximum kmer length for single kmer assemblies. Default maximum kmer is 65bp. This value must by odd.
 
=item B<-i, --increments_k>
 
The length by which the value of k increases for the next single kmer assembly. Default kmer is 2bp. This value must by even.
 
=item B<-m, --merge_k>
 
The kmer length used when merging single kmer assemblies. Default merge kmer is 39bp. This value must by odd.
 
=item B<-n, --min_read_length>
 
The minimum read length. Reads shorter than this after cleaning will be discarded. Default minimum length is 90bp.

=item B<-c, --convert_header>
 
If the illumina headers do not end in /1 or /2 use this parameter to indicat that headers need to be converted. Check your headers by typing "head [fasta filename]" and read more about illumina headers at http://en.wikipedia.org/wiki/Fastq#Illumina_sequence_identifiers.
 
 
=back

=head1 DESCRIPTION

B<OUTPUT DETAILS:>

see: https://github.com/i5K-KINBRE-script-share/transcriptome-and-genome-assembly/blob/master/KSU_bioinfo_lab/AssembleT/AssembleT_LAB.md

B<Test with sample datasets:>
 
# Find a more detailed instructions at https://github.com/i5K-KINBRE-script-share/transcriptome-and-genome-assembly/blob/master/KSU_bioinfo_lab/AssembleT/AssembleT_LAB.md
 
# Clone the Git repositories
 
git clone https://github.com/i5K-KINBRE-script-share/transcriptome-and-genome-assembly
git clone https://github.com/i5K-KINBRE-script-share/read-cleaning-format-conversion
git clone https://github.com/i5K-KINBRE-script-share/genome-annotation-and-comparison

# Make a working directory.
 
mkdir de_novo_transcriptome
cd de_novo_transcriptome

# Create symbolic links to subsampled raw RNA reads from the human breast cancer cell lines.
 
ln -s /homes/bioinfo/RNA-Seq_sample/* ~/de_novo_transcriptome/
 
# Write assembly scripts

perl ~/transcriptome-and-genome-assembly/KSU_bioinfo_lab/AssembleT/AssembleT.pl -r cell_line_reads_assembly.txt -p cell_line -s 25 -l 39 -i 2 -n 35 -m 33
 
# Clean raw reads. When these jobs are complete go to next step. Test completion by typing "status" in a Beocat session.
 
bash ~/de_novo_transcriptome/cell_line_qsubs/cell_line_qsubs_clean.sh
 
# Concatenate cleaned reads and shuffle sequences for Oases
 
bash ~/de_novo_transcriptome/cell_line_scripts/cat_reads.sh
 
# Assemble single kmer transcriptomes. When these jobs are complete go to next step. Test completion by typing "status" in a Beocat session.
 
 bash ~/de_novo_transcriptome/cell_line_qsubs/cell_line_qsubs_singlek.sh
 
# Merge single kmer transcriptomes. When these jobs are complete go to next step. Test completion by typing "status" in a Beocat session.
 
 bash ~/de_novo_transcriptome/cell_line_qsubs/cell_line_qsubs_merge.sh
 
# Cluster merged assembly with CDH. Putative transcripts that share 80% identity over 80% of the length are clustered and the longest transcript is printed in the clustered fasta file. This step will also generate assembly metrics and summarize the cleaning step results.
 
 bash ~/de_novo_transcriptome/cell_line_qsubs/cell_line_qsubs_cluster_and_qc.sh

=cut