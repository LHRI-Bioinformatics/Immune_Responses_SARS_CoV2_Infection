#!/bin/sh

#$ -N testarray
#$ -cwd
#$ -V
#$ -j y
#$ -o sge-outputNov11_log
#$ -m e
#$ -M lei.xu2@nih.gov
#$ -t 1-40
module load cellranger/7.0.0
var=(2_2_10_21)

base=${var[$SGE_TASK_ID - 1]}

echo $base >> base.txt

cellranger multi --id=$base --csv=configs_1/${base}_config.csv --jobmode=sge
