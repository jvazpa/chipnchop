#! /bin/bash

## Sample processing
## Developed by: González-Toro A., Tamargo-Azpilicueta J., Vazquez-Pacheco J.

#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -j yes


## READING PARAMS

SAMPLE=$1
i=$2
TOTNUM=$3
INSDIR=$4
RESDIR=$5
SAMDIR=$6

## QUALITY CONTROL
fastqc $SAMPLE

## MAPPING READS
bowtie2 -x ../../genome/index -U $SAMPLE -S sample_$i.sam

samtools sort -o sample_$i.bam sample_$i.sam
rm sample_$i.sam
samtools index sample_$i.bam


## TRANSCRIPT ASSEMBLY
stringtie -G ../../annotation/annotation.gtf -o sample_$i.gtf -l sample_$i sample_$i.bam


## Preparing merge list file for transcriptome merging
echo ${SAMDIR}/sample_$i.gtf >> ../../results/merge_list.txt


## BLACKBOARD -> This would be useful for parallelization
echo "sample_$i done!" >> ../../results/blackboard


echo ""
echo "Sample $i processing is done"
echo ""
