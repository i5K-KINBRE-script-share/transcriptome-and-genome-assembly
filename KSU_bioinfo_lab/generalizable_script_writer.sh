#!/bin/sh

#  USAGE: generalizable_script_writer.sh [base_filename_for_reads]
#  This script takes the full path for the input data (e.g. reads) and makes output directories and scripts with rational names. For example, the current "call program" section assumes one input file ends in "_single.fastq" and one in "_paired.fastq" and creates directories based on the name of the data and the parameter. The script can also find a file extension in the argument you pass to it, this is stored as the ${ext} variable. Change lines 36 to 41 to call your program.
#
#  Example USAGE: bash generalizable_script_writer.sh /test/good_Bittersweet_seed_aril 
#
#  Created by jennifer shelton on 4/12/13.
#

##################  define variables #################################################
job_description="oases"                         # name for all jobs
mem=120
for fullpath in "$@"
do
filename="${fullpath##*/}"                      # Strip longest match of */ from start
dir="${fullpath:0:${#fullpath} - ${#filename}}" # Substring from 0 thru pos of filename
base="${filename%.[^.]*}"                       # Strip shortest match of . plus at least one non-dot char from end
ext="${filename:${#base} + 1}"                  # Substring from len of base thru end
if [[ -z "$base" && -n "$ext" ]]; then          # If we have an extension and no base, it's really the base
base=".$ext"
ext=""
fi

echo -e "$fullpath:\n\tdir  = \"$dir\"\n\tbase = \"$base\"\n\text  = \"$ext\"\n\tfilename  = \"$filename\""

######################  add code to shell scripts to make directories with readable names (e.g. balsam_27 for the balsam read files and k=27)  #################################################
for i in {25..61}								# 25 is the lower limit for the parameter, 61 is the upper and the next four lines increment by odd numbers
do
out=$(( $i % 2 ))
if [ $out -eq 1 ]
then
echo "#!/bin/bash" > ${base}_${job_description}${i}.sh
echo "set -o verbose" >> ${base}_${job_description}${i}.sh
echo "cd ${dir}" >> ${base}_${job_description}${i}.sh
echo "mkdir ${base}_${i}" >> ${base}_${job_description}${i}.sh
echo "cd ${base}_${i}" >> ${base}_${job_description}${i}.sh

######################  call program  #################################################
echo "PATH=/homes/sheltonj/abjc/velvet_1.2.08:/homes/sheltonj/abjc/oases_0.2.08:/usr/local/hadoop/bin:/opt/sge/bin/lx24-amd64:${PATH}" >> ${base}_${job_description}${i}.sh
echo "export PATH" >> ${base}_${job_description}${i}.sh
echo "velveth oases${i}Bitter ${i} -fastq -short ${filename}_single.fastq -shortPaired -interleaved -fastq ${filename}_paired.fastq " >> ${base}_${job_description}${i}.sh
echo "velvetg oases${i}Bitter -read_trkg yes -ins_length 332" >> ${base}_${job_description}${i}.sh
echo "oases oases${i}Bitter" >> ${base}_${job_description}${i}.sh

######################  write qsubs all to one file (e.g. jobs_oases.txt) #################################################

echo "qsub -l h_rt=250:00:00,mem=${mem}G ${base}_${job_description}${i}.sh" >> jobs_${job_description}.txt
fi
done
done
