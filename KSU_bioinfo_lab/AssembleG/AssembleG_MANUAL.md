SYNOPSIS
       AssembleG.pl - The script writes scripts and qsubs to assemble illumina
       paired end reads into a de novo genome. The script 1) converts illumina
       headers if the "-c" parameter is used, 2) cleans and deduplicates raw
       reads using Prinseq http://prinseq.sourceforge.net/manual.html, 3)
       reads are the assembled multiple times with a range of values of k, 4)
       finally assembly metrics are generated for all assemblies and read
       length and number are summarized before and after cleaning.

USAGE
       perl AssembleG.pl [options]

        Documentation options:
          -help    brief help message
          -man     full documentation
        Recommended options:
          -p        project name (no spaces)(default = my_project)
          -s        shortest kmer (default = 21)
          -l        longest kmer (default = 91)
          -i        kmer increments (default = 10)
          -m        merge kmer (default = 61)
        Filtering options:
          -n        minimum read length (default = 93)
        Fastq format options:
          -c        convert fastq headers

OPTIONS
       -help   Print a brief help message and exits.

       -man    Prints the more detailed manual page with output details and
               examples and exits.

       -r, --r_list
               The filename of the user provided list of read files. The file
               should be tab separated. The first column should be the library
               type (either "pe" for paired end or "mp" for mate paired).
               Paired end or "pe" libraries are used to assemble sequence.
               Your library should be labeled "mp" if it is a long distance
               mate-pair library. "mp" libraries are not used for assembly,
               only for scaffolding. See
               https://github.com/bcgsc/abyss#assembling-multiple-libraries.

               Note: If you have multple insert lengths read files from each
               should be listed on separate lines so that Abyss can estimate
               insert length separately.

               The second column is the first read file, then the third column
               is the second read file. Example:

                pe sample_data/sample_300bp_1_R1.fastq   sample_data/sample_300bp_1_R2.fastq
                pe sample_data/sample_500bp_1_R1.fastq   sample_data/sample_500bp_1_R2.fastq
                mp sample_data/sample_3_R1.fastq   sample_data/sample_3_R2.fastq

               If a library has multiple fastq files for R1 and R2 separate
               these with commas. Example:

                pe sample_1a_R1.fastq,sample_1b_R1.fastq,sample_1c_R1.fastq   sample_1a_R2.fastq,sample_1b_R2.fastq,sample_1c_R2.fastq
                pe sample_2a_R1.fastq,sample_2b_R1.fastq,sample_2c_R1.fastq   sample_2a_R2.fastq,sample_2b_R2.fastq,sample_2c_R2.fastq
                mp sample_3a_R1.fastq,sample_3b_R1.fastq,sample_3c_R1.fastq   sample_3a_R2.fastq,sample_3b_R2.fastq,sample_3c_R2.fastq

       -s, --shortest_k
               The minimum kmer length for single kmer assemblies. Default
               minimum kmer is 21bp.

       -l, --longtest_k
               The maximum kmer length for single kmer assemblies. Default
               maximum kmer is 91bp.This value must be even.

       -i, --increments_k
               The length by which the value of k increases for the next
               single kmer assembly. Default kmer is 10bp. This value must be
               even.

       -m, --merge_k
               The kmer length used when merging single kmer assemblies.
               Default merge kmer is 61bp. This value must be odd.

       -n, --min_read_length
               The minimum read length. Reads shorter than this after cleaning
               will be discarded. Default minimum length is 93bp.

       -c, --convert_header
               If the illumina headers do not end in /1 or /2 use this
               parameter to indicat that headers need to be converted. Check
               your headers by typing "head [fasta filename]" and read more
               about illumina headers at
               http://en.wikipedia.org/wiki/Fastq#Illumina_sequence_identifiers.

DESCRIPTION
       OUTPUT DETAILS:

        See: https://github.com/i5K-KINBRE-script-share/transcriptome-and-genome-assembly/blob/master/KSU_bioinfo_lab/AssembleG/AssembleG_LAB.md

       Test with sample datasets:

       # Find a detailed instructions at
       https://github.com/i5K-KINBRE-script-share/transcriptome-and-genome-assembly/blob/master/KSU_bioinfo_lab/AssembleG/AssembleG_LAB.md