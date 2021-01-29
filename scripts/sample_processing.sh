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

## QUALITY CONTROL
fastqc $SAMPLE

## MAPPING READS
bowtie2 -x ../../genome/index -U $SAMPLE -S sample_$i.sam

samtools sort -o sample_$i.bam sample_$i.sam
rm sample_$i.sam
rm $SAMPLE
samtools index sample_$i.bam


## TRANSCRIPT ASSEMBLY
stringtie -G ../../annotation/annotation.gtf -o sample_$i.gtf -l sample_$i sample_$i.bam


## Preparing merge list file for transcriptome merging
echo ${SAMPLE_DIR}/sample_$i.gtf >> ../../results/merge_list.txt


## BLACKBOARD -> This would be useful for parallelization
echo "sample_$i done!" >> ../../results/blackboard


echo ""
echo "Sample $i processing is done"
echo ""

## Reading the blackboard (counting all lines)
NUM_PROC=$(wc -l ../../results/blackboard | awk '{ print $1 }')


if [ $NUM_PROC -eq $TOTNUM ]
then
	echo "All sample processing (both CHIP and INPUT) have succesfully finished. Transcriptome merging is running."
	echo ""
	
	bash $INSDIR/scripts/transcriptome_merging.sh ${SAMPLE}/../../results
	
	## (error-prone): qsub -o merge -N merge $INSDIR/chipnchop/transcriptome_merging.sh ${SAMPLE}/../../results $INSDIR 
fi
