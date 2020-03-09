#!/bin/sh

### Set these THREE variables for your experiment
### Make sure the project path is correct in the three scripts
### Change array number

###TODO FILE = 30_downsampleTodo.txt
#This file should be a list of .bam files from alignment. Include IgGs. Put in cutAnd_seacr directory

###METADATA FILE = metadata.txt
#This file should have one column of samples and one column of IgG. 
#If you are using IgGs from a different experiment use cp to copy them into the directory with other .bam files

##1##Choose mouse or human, place # infront of genome you aren't using
REF=/home/groups/MaxsonLab/software/ChromHMM/CHROMSIZES/hg38.txt
#REF=/home/groups/MaxsonLab/software/ChromHMM/CHROMSIZES/mm10.txt

##2##Change read number (RN) to what you want to downsample to
RN=3000000

##3##Set mark
MARK=H3K27Ac

#Executables these don't need to change
SAMTOOLS=/home/groups/MaxsonLab/software/miniconda3/bin/samtools
BEDTOOLS=/home/groups/MaxsonLab/smithb/KLHOXB_TAG_09_19/Dense_ChromHMM/bedtools2/bin/bedtools
SEACR=/home/groups/MaxsonLab/software/SEACR/SEACR_1.1.sh
umask 007
