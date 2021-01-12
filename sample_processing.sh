## Developed by: González-Toro A., Tamargo-Azpilicueta J., Vazquez-Pacheco J.

## READING PARAMS

SAMPLE=$1

## QUALITY CONTROL

fastqc $SAMPLE

## MAPPING READS

bowtie2 -x ../../genome/index -U $SAMPLE -S chip.sam

## BLACKBOARD
