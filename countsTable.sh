#!/bin/bash

###PURPOSE: Make counts table from seacr peaks###

#SBATCH --partition          exacloud                # partition (queue)
#SBATCH --nodes              1                       # number of nodes
#SBATCH --ntasks             1                       # number of "tasks" to be allocated for the job
#SBATCH --ntasks-per-core    1                       # Max number of "tasks" per core.
#SBATCH --cpus-per-task      1                       # Set if you know a task requires multiple processors
##SBATCH --mem               8000                  # memory pool for each node
#SBATCH --time               0-24:00                 # time (D-HH:MM)
#SBATCH --output             counts_%A.out     # Standard output
#SBATCH --error              counts_%A.err     # Standard error

### Executable
BEDTOOLS=/home/groups/MaxsonLab/smithb/KLHOXB_TAG_09_19/Dense_ChromHMM/bedtools2/bin/bedtools

### SET I/O VARIABLES
PROJECT=/home/groups/MaxsonLab/smithb/KASUMI_TAG_12_19
MARK=H3K4me1
IN=$PROJECT/process/30_downsampled/seacr
IN2=$PROJECT/process/30_downsampled/bams
OUT=$PROJECT/process/30_downsampled/counts
mkdir -p $OUT

echo "Counts table:"
cmd="$BEDTOOLS multicov -bams $IN2/*$MARK\.ds.sorted.bam -bed $IN/$MARK\_bed_for_multicov.bed > $OUT/$MARK\_counts.txt"
echo $cmd
eval $cmd
