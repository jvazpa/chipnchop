## Developed by: González-Toro A., Tamargo-Azpilicueta J., Vazquez-Pacheco J.

#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -j yes


## READING PARAMS

SAMPLE=$1
j=$2
NUMOFSAMP=$3
INSDIR=$4

## QUALITY CONTROL

fastqc $SAMPLE

## MAPPING READS

bowtie2 -x ../../genome/index -U $SAMPLE -S chip.sam

## BLACKBOARD


##############################################################






## Sample quality control and read mapping to reference genome
if [ -f <accession_sra>_2.fastq ]
then
   fastqc <accession_sra>_1.fastq
   fastqc <accession_sra>_2.fastq

   hisat2 --dta -x ../../genome/index -1 <accession_sra>_1.fastq -2 <accession_sra>_2.fastq -S sample<N>.sam
else
   fastqc <accession_sra>_1.fastq

   hisat2 --dta -x ../../genome/index -U <accession_sra>_1.fastq -S sample<N>.sam
fi

## Generting sorted bam file
samtools sort -o sample<N>.bam sample<N>.sam
rm sample<N>.sam
rm *.fastq
samtools index sample<N>.bam
bamCoverage -bs 10 --normalizeUsing CPM --bam sample<N>.bam -o sample<N>.bw


## Transcript assembly
stringtie -G ../../annotation/annotation.gtf -o sample<N>.gtf -l sample<N> sample<N>.bam

## Preparing merge list file for transcriptome merging
echo /home/<grupo>/<exp>/samples/sample<N>/sample<N>.gtf >> ../../results/merge_list.txt

## Gene Expression Quantification
stringtie -e -B -G ../../annotation/annotation.gtf -o sample<N>.gtf sample<N>.bam
