#!/bin/bash           # the shell language when run outside of the job scheduler
#                     # lines starting with #$ is an instruction to the job scheduler
#$ -S /bin/bash       # the shell language when run via the job scheduler [IMPORTANT]
#$ -cwd               # job should run in the current working directory
#$ -j y               # STDERR and STDOUT should be joined
#$ -l mem_free=2G     # job requires up to 2 GiB of RAM per slot
#$ -l h_rt=24:00:00   # job requires up to 24 hours of runtime

module load CBI miniconda3/23.5.2-0-py311

conda activate kallisto

INDEX="/wynton/home/marshall/aralbright/deepa/scoe.idx"

FILES="/wynton/home/marshall/aralbright/deepa/raw/"

OUTPUTDIR="kallisto_output"

for first in $FILES*1_001.fastq.gz
do
  output=$(echo $first | cut -c44-46)
  echo "$OUTPUTDIR/$output"
  echo "$first"
  echo "${first/1_001.fastq.gz/2_001.fastq.gz}"

  kallisto quant -i $INDEX -o $OUTPUTDIR/$output -b30 -t 16 $first ${first/1_001.fastq.gz/2_001.fastq.gz}

done