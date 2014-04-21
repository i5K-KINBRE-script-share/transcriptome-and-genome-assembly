SYNOPSIS

AssembleT.pl - The script writes scripts and qsubs to assemble illumina
       paired end reads into a de novo transcriptome. The script 1) converts
       illumina headers if the "-c" parameter is used, 2) cleans and
       deduplicates raw reads using Prinseq
       http://prinseq.sourceforge.net/manual.html, 3) index the reference
       genome for mapping, 4) reads are the assembled multiple times with a
       range of values of k, 5) these assemblies are merged with Oases using a
       merge kmer value, 6) then the merged assembly is clusted with CDHit to
       take the longest of similar putative transcripts, 7) finally assembly
       metrics are generated for all assemblies and read length and number are
       summarized before and after cleaning.

USAGE

       perl AssembleT.pl [options]

        Documentation options:
          -help    brief help message
          -man     full documentation
        Required options:
          -r        filename for file with tab separated list of fastq forward and reverse read files 
        Recommended options:
          -p        project name (no spaces)(default = my_project)
          -s        shortest kmer (default = 25)
          -l        longest kmer (default = 65)
          -i        kmer increments (default = 2)
          -m        merge kmer (default = 39)
        Filtering options:
          -n        minimum read length (default = 90)
        Fastq format options:
          -c        convert fastq headers

OPTIONS

       -help   Print a brief help message and exits.

       -man    Prints the more detailed manual page with output details and
               examples and exits.
               
       -r, --r_list
               The filename of the user provided list of read files. The file
               should be tab separated with the first read file, then the
               second read file. Example:

                sample_data/sample_1_R1.fastq   sample_data/sample_1_R2.fastq

               If a sample has multiple fastq files for R1 and R2 separate
               these with commas. Example:

                sample_data/sample_1a_R1.fastq,sample_data/sample_1b_R1.fastq,sample_data/sample_1c_R1.fastq   sample_data/sample_1a_R2.fastq,sample_data/sample_1b_R2.fastq,sample_data/sample_1c_R2.fastq

       -s, --shortest_k
               The minimum kmer length for single kmer assemblies. Default
               minimum kmer is 25bp. This value must by odd.

       -l, --longtest_k
               The maximum kmer length for single kmer assemblies. Default
               maximum kmer is 65bp. This value must by odd.

       -i, --increments_k
               The length by which the value of k increases for the next
               single kmer assembly. Default kmer is 2bp. This value must by
               even.

       -m, --merge_k
               The kmer length used when merging single kmer assemblies.
               Default merge kmer is 39bp. This value must by odd.

       -n, --min_read_length
               The minimum read length. Reads shorter than this after cleaning
               will be discarded. Default minimum length is 90bp.
               
       -c, --convert_header
               If the illumina headers do not end in /1 or /2 use this
               parameter to indicat that headers need to be converted. Check
               your headers by typing "head [fasta filename]" and read more
               about illumina headers at
               http://en.wikipedia.org/wiki/Fastq#Illumina_sequence_identifiers.

DESCRIPTION


OUTPUT DETAILS:

See:

https://github.com/i5K-KINBRE-script-share/transcriptome-and-genome-assembly/blob/master/KSU_bioinfo_lab/AssembleT/AssembleT_LAB.md


Test with sample datasets:

See:

https://github.com/i5K-KINBRE-script-share/transcriptome-and-genome-assembly/blob/master/KSU_bioinfo_lab/AssembleT/AssembleT_LAB.md
