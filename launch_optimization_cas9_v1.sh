#!/bin/bash
module load R/4.0.4

inputfile=$1
errormatrix=$2
outputfile=$3
mymodel=$4
#target=$5
#induction=$6
#estimaterr=$7


#for target in CRTISO Psy1 CRTISO_49and50bp;do
#for induction in gRNA RNP;do
#bsub -q new.long -J ${target}${induction} -e ${target}${induction}.e -o ${target}${induction}.o -R rusage[mem=16000] 
#./optimize_model_backbone.R  -t ${target} -i $induction -m modelDSBs1i1_fullimpreciseDSB -e E_errorsfromunbroken -n 10000 -o ${target}_${induction}_optimization.tsv
./optimize_model_backbone.v1.R -m ${mymodel} -n 1000 -E $errormatrix -T $inputfile -o $outputfile
#done
#done
