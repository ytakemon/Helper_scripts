#!/bin/bash -l
#PBS-l nodes=1:ppn=3,walltime=00:30:00
module load fastqc/0.11.3
cd ${fastq_dir}

fastqc -o fastqc_check/ ${sample}_R1*.fastq
fastqc -o fastqc_check/ ${sample}_R2*.fastq
