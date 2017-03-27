##bam2bai.sh

#!/bin/bash -l
#PBS -l nodes=1:ppn=5,walltime=24:00:00

module load samtools/1.2 
cd ${tophat_dir}

samtools index ${sample}/accepted_hits.bam 
