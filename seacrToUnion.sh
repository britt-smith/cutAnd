#!/bin/bash

###PURPOSE: find SEACR peaks from fragment bedgraphs, make consensus peak file, result in counts table###

#SBATCH --partition          exacloud                # partition (queue)
#SBATCH --nodes              1                       # number of nodes
#SBATCH --ntasks             1                       # number of "tasks" to be allocated for the job
#SBATCH --ntasks-per-core    1                       # Max number of "tasks" per core.
#SBATCH --cpus-per-task      1                       # Set if you know a task requires multiple processors
#SBATCH --mem-per-cpu        4000                   # Memory required per allocated CPU (mutually exclusive with mem)
##SBATCH --mem               16000                  # memory pool for each node
#SBATCH --time               0-24:00                 # time (D-HH:MM)
#SBATCH --output             seacr_%A_%a.out     # Standard output
#SBATCH --error              seacr_%A_%a.err     # Standard error
#SBATCH --array              1-8                    # sets number of jobs in array


PROJECT=/home/groups/MaxsonLab/smithb/KASUMI_TAG_12_19

###############################################################

source $PROJECT/cutAnd_seacr/cutAndConfig.sh
IN=$PROJECT/process/cutAnd_seacr/beds
IN2=$PROJECT/process/cutAnd_seacr/bams
OUT=$PROJECT/process/cutAnd_seacr/seacr
TODO=$PROJECT/cuttag/todo/30_downsampleTodo.txt
mkdir -p $OUT

### Record slurm info
echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID

### Get file info
currINFO=`awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}' $TODO`

### Set variables
NAME=${currINFO%%.bam}
DATA=$NAME.ds.bedgraph

### Merge peaks
if [[ "$DATA" == *1_* ]]; then
	REF1=$NAME.relaxed.bed
	REF2=`echo $REF1 | sed 's/1_/2_/'`
	BASE=`echo $NAME | sed 's/1_/_/'`
fi

cmd="$BEDTOOLS intersect -a $OUT/$REF1 -b $OUT/$REF2 -wa | cut -f1-3 | sort | uniq > $OUT/$BASE\_Ref1_merge.bed"
echo "Merge -wa"
echo $cmd
eval $cmd

cmd="$BEDTOOLS intersect -a $OUT/$REF1 -b $OUT/$REF2 -wb | cut -f1-3 | sort | uniq > $OUT/$BASE\_Ref2_merge.bed"	
echo "Merge -wb"
echo $cmd
eval $cmd

### Make superset of peaks
echo "Replicate merge total:"
cmd="cat $OUT/*$MARK\_*_merge.bed | sort -k1,1 -k2,2n | $BEDTOOLS merge | tee $OUT/$MARK\_merge.bed | wc -l"
echo $cmd
eval $cmd

