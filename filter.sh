#!/bin/bash -l
#PBS -l nodes=1:ppn=6,walltime=06:00:00

#	First take trimmed reads and make temprary file, then use temp to make file file.
module load cutadapt/1.2.1
input_dir=${fastq_dir}
output_dir=${input_dir}/filter_36 #make filter_36 directory

#generic
cutadapt --minimum-length 36 -o ${output_dir}/temp_${g}_R1_trimmed_filtered.fastq -p ${output_dir}/temp_${g}_R2_trimmed_filtered.fastq \
${input_dir}/${g}_R1.fastq ${input_dir}/${g}_R2.fastq

cutadapt --minimum-length 36 -o ${output_dir}/${g}_R2_trimmed_filtered.fastq -p ${output_dir}/${g}_R1_trimmed_filtered.fastq \
${output_dir}/temp_${g}_R2_trimmed_filtered.fastq ${output_dir}/temp_${g}_R1_trimmed_filtered.fastq

