To begin this lab your should read about the software we will be using. Prinseq will be used to clean raw reads. Priseq cleaning is highly customizable. You can see a detailed parameter list by typing "perl /homes/sheltonj/abjc/prinseq-lite-0.20.3/prinseq-lite.pl -h" or by visiting their manual at http://prinseq.sourceforge.net/manual.html. You can read a detailed list of parameter options for Oases transcriptome assembler by typing "/homes/sheltonj/abjc/oases_0.2.08/oases" or visit the Oases manual at https://www.ebi.ac.uk/~zerbino/oases/. 

We will be using the script "transcriptome_assembly_pipeline.pl" to organize our working directory and write scripts to clean our reads, assemble our de novo transcriptome for human breast cancer cell lines, and summarize our assembly metrics.

To find out more about the parameters for "transcriptome_assembly_pipeline.pl" run "perl ~/transcriptome-and-genome-assembly/KSU_bioinfo_lab/transcriptome_assembly_pipeline/transcriptome_assembly_pipeline.pl -man" or visit its manual at https://github.com/i5K-KINBRE-script-share/transcriptome-and-genome-assembly/tree/master/KSU_bioinfo_lab/transcriptome_assembly_pipeline/transcriptome_assembly_README.md.

###Step 1: Clone the Git repository 

        git clone https://github.com/i5K-KINBRE-script-share/transcriptome-and-genome-assembly
        
###Step 2: Create project directory and add your input data to it

Make a working directory.

        mkdir de_novo_transcriptome
        cd de_novo_transcriptome

Create symbolic links to raw RNA reads from the human breast cancer cell lines. We are using Illumina paired end reads from four breast cancer cell lines. This data was used in a biological visualization competition that illumina held a few years ago http://blog.expressionplot.com/2011/03/15/idea-challenge-2011-illuminas-data-excellence-award/.

        ln -s /homes/bioinfo/RNA-Seq_sample/* ~/de_novo_transcriptome/
        
Step 3: Write assembly scripts

Check to see if your fastq headers end in "/1" or "/2" (if they do not you must add the parameter "-c" when you run "transcriptome_assembly_pipeline.pl"

        head ~/de_novo_transcriptome/*_1.fastq
        
Your output will look similar to the output below for the sample data. Because these reads end in "/1" or "/2" we will not add "-c" when we call "transcriptome_assembly_pipeline.pl".

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

Step 4: Run prinseq and the assembly scripts

        perl ~/transcriptome-and-genome-assembly/KSU_bioinfo_lab/transcriptome_assembly_pipeline/transcriptome_assembly_pipeline.pl -r cell_line_reads_assembly.txt -p cell_line -s 25 -l 41 -i 2 


## remember to change prinseq commands!!!!


transcriptome_assembly_pipeline.pl
