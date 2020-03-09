#!/bin/bash

###PURPOSE: Make counts table from seacr peaks###

#SBATCH --partition          exacloud                # partition (queue)
#SBATCH --nodes              1                       # number of nodes
#SBATCH --ntasks             1                       # number of "tasks" to be allocated for the job
#SBATCH --ntasks-per-core    1                       # Max number of "tasks" per core.
#SBATCH --cpus-per-task      1                       # Set if you know a task requires multiple processors
##SBATCH --mem               16000                  # memory pool for each node
#SBATCH --time               0-24:00                 # time (D-HH:MM)
#SBATCH --output             counts_%A.out     # Standard output
#SBATCH --error              counts_%A.err     # Standard error


PROJECT=/home/groups/MaxsonLab/smithb/KASUMI_TAG_12_19

########################################################

source $PROJECT/cutAnd_seacr/cutAndConfig.sh

IN=$PROJECT/process/cutAnd_seacr/seacr
IN2=$PROJECT/process/cutAnd_seacr/bams
OUT=$PROJECT/process/cutAnd_seacr/counts
mkdir -p $OUT

echo "Superset total:"
cmd="`cat $IN/$MARK\_merge.bed | awk '{$3=$3"\t""peak_"NR}1' OFS="\t" | tee $IN/$MARK\_bed_for_multicov.bed| wc -l`"
echo $cmd

echo "Counts table:"
cmd="$BEDTOOLS multicov -bams $IN2/*$MARK\.ds.sorted.bam -bed $IN/$MARK\_bed_for_multicov.bed > $OUT/$MARK\_counts.txt"
echo $cmd
eval $cmd

ls $IN2/*$MARK.ds.sorted.bam | grep -o '[^/]*$' | cut -d_ -f1 | tr "\n" "\t" | awk '{print "\t\t\t\t" $0}' | cat -  $OUT/$MARK\_counts.txt > $OUT/$MARK\_counts_labeled.txt

