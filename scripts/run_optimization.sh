#!/bin/bash
#module load R/4.0.4
module load R/4.1.0
#module load R

inputfile=$1
errormatrix=$2
outputfile=$3
mymodel=$4
myn=$5
errorflag=$6
inductioncurve=$7

./DSBtimecourse_optimizer.R -m ${mymodel} -n $myn -E $errormatrix -T $inputfile -o $outputfile -l $errorflag -z ${inductioncurve}
