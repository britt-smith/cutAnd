# cutAnd_seacr
A pipeline to analyze CUT&amp;Tag/RUN with seacr peak calling.

**Exploratory downsampling**
Ted found that the Cell Signaling CUT&RUN kit suggests to use the lowest number of mapped E. coli reads to calculate the fraction of reads to downsample to

sample | mappedreads_sample | mappedreads_ecoli | fraction | downsample_reads
-------|--------------------|-------------------|----------|-----------------
AG2_H3K27Ac | 2439841 | 85854 | 1 | 2439841
ND1_H3K27Ac | 5110814 | 1380136 | 0.06 | 317927

However, the ratios were all over the place, so to us it makes more since to bring everything down to the same number of reads instead of having variable number of reads.

Using technique found on https://davemcg.github.io/post/easy-bam-downsampling/ to downsample to number of reads opposed to percentage

*Troubleshooting filtering*

I used 31_sbatchBam2Bed.sh and 32_sbatchBed2BG_woEcoli.sh to generate bedgraphs. I noticed that it did not look right at all after downsampling

Using A1_H3K27Ac as test

`samtools view -c example.bam`

file | reads 
-------|------
original bam | 6711828
sorted bam | 6711828
filtered bam | 6711828
downsampled bam | 2997574
sorted bam | 2997574
bedpe | 1498787
bedgraph | 411734

The **-n** is *very* important for bam to bed, but sorting before making the index requires sorting *without* the -n

srun samtools sort -n A1_H3K27Ac.ds.bam > A1_H3K27Ac.ds.sorted.bam

**SEACR peaks**

sample | seacr peaks rep1 | seacr peaks rep2 | ref 1 merge | ref 2 merge
-------|--------------------|-------------------|----------|-----------------
ND_H3K4me3 | 8491 | 8588 | 8107 | 8118
A_H3K4me3 | 9549 | 9982 | 9328 | 9338
G_H3K4me3 | 8346 | 9770 | 8277 | 8304
AG_H3K4me3 | 6189 | 9864 | 6152 | 6155
**TOTAL** | 9534
