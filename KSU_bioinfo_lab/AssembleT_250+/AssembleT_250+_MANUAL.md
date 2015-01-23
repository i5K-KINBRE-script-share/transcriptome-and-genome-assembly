SYNOPSIS

AssembleT_250+.pl - The script writes scripts and qsubs to assemble
    Illumina paired end reads that are 250 bp or longer into a de novo
    transcriptome. The script 1) converts illumina headers if the "-c"
    parameter is used, 2) cleans and deduplicates raw reads using Prinseq
    http://prinseq.sourceforge.net/manual.html, 3) reads are the assembled
    into supereads wuth Masurca, 5) these supereads are converted to FASTA
    file format and assembled with Mira.
       

UPDATES

####AssembleT_250+.pl Version 1.0 01/20/15


USAGE

       perl AssembleT_250+.pl [options]

        Documentation options:
          -help    brief help message
          -man     full documentation
        Required options:
          -r        full path for file with tab separated list of fastq forward and reverse read files 
        Recommended options:
          -p        project name (no spaces)(default = my_project)
        Filtering options:
          -n        minimum read length (default = 90)
        Fastq format options:
          -c        convert fastq headers
        Additional options:
       	  -j        max number of threads for superead assembly

OPTIONS

       -help   Print a brief help message and exits.

       -man    Prints the more detailed manual page with output details and
               examples and exits.
               
       -r, --r_list
               This is the the full path (path and filename) of the user
            	provided list of read files. The file should be tab separated
            	with the two-letter Masurca library prefix codes, then the first
            	read file, then the second read file. Example:

                D1	sample_data/sample_1_R1.fastq   sample_data/sample_1_R2.fastq

               If a sample has multiple fastq files for R1 and R2 separate
               these with commas. Example:

                D2	sample_data/sample_1a_R1.fastq,sample_data/sample_1b_R1.fastq,sample_data/sample_1c_R1.fastq   sample_data/sample_1a_R2.fastq,sample_data/sample_1b_R2.fastq,sample_data/sample_1c_R2.fastq

       -n, --min_read_length
               The minimum read length. Reads shorter than this after cleaning
               will be discarded. Default minimum length is 90bp.
               
       -c, --convert_header
               If the illumina headers do not end in /1 or /2 use this
               parameter to indicat that headers need to be converted. Check
               your headers by typing "head [fasta filename]" and read more
               about illumina headers at
               http://en.wikipedia.org/wiki/Fastq#Illumina_sequence_identifiers.
       -j, --max_threads
            The max number of threads for superead assembly on Beocat
            (default=16).

DESCRIPTION


OUTPUT DETAILS:

See:

https://github.com/i5K-KINBRE-script-share/transcriptome-and-genome-assembly/blob/master/KSU_bioinfo_lab/AssembleT_250+/AssembleT_250+_LAB.md


Test with sample datasets:

See:

https://github.com/i5K-KINBRE-script-share/transcriptome-and-genome-assembly/blob/master/KSU_bioinfo_lab/AssembleT_250+/AssembleT_250+_LAB.md
