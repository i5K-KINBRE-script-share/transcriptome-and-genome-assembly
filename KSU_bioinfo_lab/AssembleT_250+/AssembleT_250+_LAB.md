![alttext](https://raw.githubusercontent.com/i5K-KINBRE-script-share/transcriptome-and-genome-assembly/master/images/ngs_pipelines_on_beocat.png)
##Transcriptome assembly with Masurca and Mira

All of the scripts you will need to complete this lab as well as the sample dataset will be copied to your Beocat directory as you follow the instructions below. You should type or paste the text in the beige code block into your terminal as you follow along with the instructions below. If you are not used to commandline, practice with real data is one of the best ways to learn.

If you would like a quick primer on basic linux commands try these 10 minute lessons from Software Carpentry http://software-carpentry.org/v4/shell/index.html. To learn to start using Beocat and begin using the terminal go to https://github.com/i5K-KINBRE-script-share/FAQ/blob/master/UsingBeocat.md. Learn how to download files from Beocat at https://github.com/i5K-KINBRE-script-share/FAQ/blob/master/BeocatEditingTransferingFiles.md.

We will be using the script "AssembleT.pl" to organize our working directory and write scripts to clean our reads, assemble our de novo transcriptome for human breast cancer cell lines, and summarize our assembly metrics.

To begin this lab your should read about the software we will be using. Prinseq will be used to clean raw reads. Priseq cleaning is highly customizable. You can see a detailed parameter list by typing `perl /homes/bioinfo/bioinfo_software/prinseq-lite.pl -h` or by visiting their manual at http://prinseq.sourceforge.net/manual.html. You can read a detailed list of parameter options for Masurca superead assembler by typing `/homes/bioinfo/bioinfo_software/MaSuRCA/bin/masurca -help` or visit the Masurca manual at ftp://ftp.genome.umd.edu/pub/MaSuRCA/MaSuRCA_QuickStartGuide.pdf. 

To find out more about the parameters for "AssembleT_250+.pl" run "perl ~/transcriptome-and-genome-assembly/KSU_bioinfo_lab/AssembleT_250+/AssembleT_250+.pl -man" or visit its manual at https://github.com/i5K-KINBRE-script-share/transcriptome-and-genome-assembly/tree/master/KSU_bioinfo_lab/AssembleT_250+/AssembleT_250+_MANUAL.md.

###Step 1: Clone the Git repositories 

        git clone https://github.com/i5K-KINBRE-script-share/transcriptome-and-genome-assembly
        
###Step 2: Create project directory and add your input data to it

Make a working directory.

        mkdir de_novo_transcriptome
        cd de_novo_transcriptome

Create symbolic links to raw RNA reads from the human breast cancer cell lines. Creating a symbolic link rather than copying avoids wasting disk space and protects your raw data from being altered. We are using Illumina paired end reads from four breast cancer cell lines. This data was used in a biological visualization competition that illumina held a few years ago http://blog.expressionplot.com/2011/03/15/idea-challenge-2011-illuminas-data-excellence-award/.

        ln -s /homes/bioinfo/pipeline_datasets/AssembleT/* ~/de_novo_transcriptome/
        
Step 3: Write assembly scripts

Check to see if your fastq headers end in "/1" or "/2" (if they do not you must add the parameter "-c" when you run "AssembleT.pl"

        head ~/de_novo_transcriptome/*_1.fastq
        
Your output will look similar to the output below for the sample data. Because these reads end in "/1" or "/2" we will not add "-c" when we call "AssembleT.pl".

        ==> /homes/bioinfo/de_novo_transcriptome/BT20_paired-end_RNA-seq_subsampled_1.fastq <==
        @HWUSI-EAS1794_0001_FC61KOJ:4:30:19389:13787#0/1
        GCGGCCCGGCCCCGGCCCCCTGCTCGTTGGCTGTGGCAGGGCCGCCGTGG
        +
        HHHHHHHGEHHHHDHDHHHHBDGBBC@CAC?8C><AAAACD>DDB?####
        @HWUSI-EAS1794_0001_FC61KOJ:4:57:10821:2162#0/1
        CAGATATCGAAGATGAAGACTTAAAGTTAGAGCTGCGACGACTACGAGAT
        +
        IIHIIIIIIHHIHIIIIIEIIIIIIIIIIIHHII@IHIIIIHHEIIHIID
        @HWUSI-EAS1794_0001_FC61KOJ:4:75:5014:13576#0/1
        CTCAGCCACCAGCAGCGGCACCCCCATCTGCAGTTGGCTCTTCTGCTGCT
        
Call "AssembleT.pl". Our reads are only 50 bp long so we are setting our minimum read length to 35 bp. Generally you want to keep this length ~10 bp shorter than our read length. We would also raise the longest kmer value "-l" to ~61 if our reads were 100bp.

        perl ~/transcriptome-and-genome-assembly/KSU_bioinfo_lab/AssembleT/AssembleT.pl -r cell_line_reads_assembly.txt -p cell_line

###Step 4: Run prinseq and the assembly scripts

Clean raw reads. When these jobs are complete go to next step. Test completion by typing `status` in a Beocat session. Download the `.gd` files in the `~/de_novo_transcriptome/cell_line__prinseq` directory and upload them to http://edwards.sdsu.edu/cgi-bin/prinseq/prinseq.cgi?report=1 to evaluate read quality pre and post cleaning.

        bash ~/de_novo_transcriptome/cell_line_qsubs/cell_line_qsubs_clean.sh

Concatenate cleaned reads and Assemble transcriptomes. Test completion by
       typing "status" in a Beocat session.

        bash ~/de_novo_transcriptome/cell_line_qsubs/testing_qsubs_Masurca_Mira.sh
        
        

        
        
        


