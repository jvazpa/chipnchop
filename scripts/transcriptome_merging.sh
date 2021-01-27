#! /bin/bash

## Transcriptome merging
## Developed by: González-Toro A., Tamargo-Azpilicueta J., Vazquez-Pacheco J.

#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -j yes

## Reading parameters

RESDIR=$1

echo ""
echo "============================"
echo "|  Transcriptome merging   |"
echo "============================"
echo ""

echo "Transcriptome merging script is running. This is about to end. Just hang on!"


## Accessing results folder
cd $RESDIR

## Merging sample transcriptomes
stringtie --merge -G ../annotation/annotation.gtf -o stringtie_merged.gtf merge_list.txt

## Comparing our assembly with the reference
gffcompare -r ../annotation/annotation.gtf -G -o comparison stringtie_merged.gtf

## Picking percentage of the previous comparison. Cheking differences with reference annotation.

DIFF=$(grep "Novel loci:" comparison.stats | awk '{ print $5 }' | cut -d"%" -f 1 | cut -d"." -f 1)

if [ $DIFF -ge 5 ]
then
	echo "========================================================="
	echo "| WARNING: Reference transcriptome differs more than 5% |"
	echo "========================================================="
fi

echo ""
echo "Transcriptome merging have finished."

