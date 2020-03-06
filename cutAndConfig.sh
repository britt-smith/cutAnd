#!/bin/sh

### Set variables for your experiment, you shouldn't need to change anything in the scripts except the array number ###

	#TODO=30_downsampleTodo.txt This file should be a list of .bam files from alignment. Include IgGs or put the path to IgG in section indicated below.

### downsampling.sh ###

#Executables these don't need to change
SAMTOOLS=/home/groups/MaxsonLab/software/miniconda3/bin/samtools
BEDTOOLS=/home/groups/MaxsonLab/smithb/KLHOXB_TAG_09_19/Dense_ChromHMM/bedtools2/bin/bedtools
SEACR=/home/groups/MaxsonLab/software/SEACR/SEACR_1.1.sh

#Set your project directory
PROJECT=/home/groups/MaxsonLab/smithb/KASUMI_TAG_12_19

#Choose mouse or human, place # infront of genome you aren't using
REF=/home/groups/MaxsonLab/software/ChromHMM/CHROMSIZES/hg38.txt
#REF=/home/groups/MaxsonLab/software/ChromHMM/CHROMSIZES/mm10.txt

#Change read number (RN) to what you want to downsample to
RN=3000000

#Path to IgG

### seacrToUnion.sh and countsTable.sh ###

#Set mark
MARK=H3K27Ac

