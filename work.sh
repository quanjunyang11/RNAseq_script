#!/bin/bash
echo ==========start at : `date` ========== ;
set -ex ;
hisat2 -p 8 --dta -x Mus_musculus.genome.fa -1 CleanData/SampleName/SampleName.clean_1.fq.gz -2 CleanData/SampleName/SampleName.clean_2.fq.gz -S SampleName.sam --new-summary --summary-file SampleName.Mapping.Summary.txt ;
samtools sort -@ 8 -o SampleName.bam SampleName.sam ;
rm SampleName.sam ;
mkdir -p Stringtie_eb/sample_SampleName ;
/data/pipeline/RNAseq-pipe/software/stringtie SampleName.bam -e -B -p 8 -G Mus_musculus.gtf -o Stringtie_eb/sample_SampleName/SampleName.gtf ;

echo ==========end at : `date` ==========
