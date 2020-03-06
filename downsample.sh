#!/bin/bash

### PURPOSE: Filter bam files then downsample to lowest number of reads (within reason, will most likely be ~3 million) 
###				and produce fragment bedgraphs for SEACR analysis ###

#SBATCH --partition          exacloud                # partition (queue)
#SBATCH --nodes              1                       # number of nodes
#SBATCH --ntasks             1                       # number of "tasks" to be allocated for the job
#SBATCH --ntasks-per-core    1                       # Max number of "tasks" per core.
#SBATCH --cpus-per-task      4                       # Set if you know a task requires multiple processors
#SBATCH --mem                16000                  # memory pool for each node
#SBATCH --time               0-24:00                 # time (D-HH:MM)
#SBATCH --output             downsample_%A_%a.out           # Standard output
#SBATCH --error              downsample_%A_%a.err           # Standard error
#SBATCH --array              1-8                     # sets number of jobs in array

### Executable
SAMTOOLS=/home/groups/MaxsonLab/software/miniconda3/bin/samtools
BEDTOOLS=/home/groups/MaxsonLab/smithb/KLHOXB_TAG_09_19/Dense_ChromHMM/bedtools2/bin/bedtools
SEACR=/home/groups/MaxsonLab/software/SEACR/SEACR_1.1.sh

### SET VARIABLES
#Set your project directory
PROJECT=/home/groups/MaxsonLab/smithb/KASUMI_TAG_12_19

#Choose mouse or human, place # infront of genome you aren't using
REF=/home/groups/MaxsonLab/software/ChromHMM/CHROMSIZES/hg38.txt
#REF=/home/groups/MaxsonLab/software/ChromHMM/CHROMSIZES/mm10.txt

#change read number (RN) to what you want to downsample to
RN=3000000

### Other arguments
NORM="norm"
THRESH="relaxed"

#These don't need to change
TODO=$PROJECT/cuttag/todo/30_downsampleH3K4me3Todo.txt
IN=$PROJECT/process/20_alignments
OUT1=$PROJECT/process/30_downsampled/bams
OUT2=$PROJECT/process/30_downsampled/beds
OUT3=$PROJECT/process/30_downsampled/seacr
mkdir -p $OUT1
mkdir -p $OUT2
mkdir -p $OUT3

### Record slurm info
echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID

### Get file info
currINFO=`awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}' $TODO`
NAME=${currINFO%%.bam}
echo "Name:"
echo $NAME
CTL=${NAME%%_*}_IgG.ds.bedgraph
DATA=$NAME.ds.bedgraph

#Sort bam files
cmd="$SAMTOOLS sort $IN/$currINFO -o $OUT1/$NAME\.sorted.bam"
echo $cmd
#eval $cmd


### Calculate the fraction of reads to downsample to a certain number (i.e. 3 million reads) (Eye Bioinformatician)
# ds in the output name stands for downsample, if fraction >1 it will print .99, adding 42 for the random seed
#frac=$(samtools idxstats $OUT1/$NAME\.sorted.bam | cut -f3 | awk -v DS="$RN" 'BEGIN {total=0} {total += $1} END {frac=DS/total; if (frac > 1) {print .99} else {print frac}}')
#scale=`echo "42+$frac" | bc`
echo "Scale:"
echo $scale
cmd="$SAMTOOLS view -bs $scale $OUT1/$NAME\.sorted.bam > $OUT1/$NAME\.ds.bam"
echo "Downsample"
echo $cmd
#eval $cmd


### Sort bam by locus
cmd="$SAMTOOLS sort $OUT1/$NAME\.ds.bam > $OUT1/$NAME\.ds.sorted.bam"
echo "Sort bam"
echo $cmd
#eval $cmd

### Index bam files
cmd="$SAMTOOLS index $OUT1/$NAME\.ds.sorted.bam"
echo "Index bam"
echo $cmd
#eval $cmd

### Sort bam by name
cmd="$SAMTOOLS sort -n $OUT1/$NAME\.ds.bam > $OUT1/$NAME\.ds.name.sorted.bam"
echo "Sort bam"
echo $cmd
#eval $cmd

### Bam to bed
cmd="$BEDTOOLS bamtobed -bedpe -i $OUT1/$NAME\.ds.name.sorted.bam > $OUT2/$NAME\.ds.bed"
echo "Bam to bed"
echo $cmd
#eval $cmd

### Commands
cleanBed="awk '\$1==\$4 && \$6-\$2 < 1000 {print \$0}' $OUT2/$NAME\.ds.bed > $OUT2/$NAME.ds.clean.bed"
getFrag="cut -f 1,2,6 $OUT2/$NAME.ds.clean.bed > $OUT2/$NAME.ds.fragments.bed"
sortFrag="sort -k1,1 -k2,2n -k3,3n $OUT2/$NAME.ds.fragments.bed > $OUT2/$NAME.ds.sortfragments.bed"
bedgraph="$BEDTOOLS genomecov -bg -i $OUT2/$NAME.ds.sortfragments.bed -g $REF > $OUT2/$NAME.ds.bedgraph"

### Run
echo "Clean bed"
echo $cleanBed
#eval $cleanBed

echo "Get fragments"
echo $getFrag
#eval $getFrag

echo "Sort fragments"
echo $sortFrag
#eval $sortFrag

echo "Convert to bedgraph"
echo $bedgraph
#eval $bedgraph

### Seacr Peaks
cmd="$SEACR $OUT2/$DATA $OUT2/$CTL $NORM $THRESH $OUT3/$NAME"
echo $cmd
eval $cmd

### Number of peaks
echo "Number of peaks:"
cmd="cat $OUT3/$NAME\.relaxed.bed | wc -l"
eval $cmd
