#!/bin/bash
module load R/4.0.4

inputf=$1
NMAX=$2
#for target in CRTISO Psy1 CRTISO_49and50bp;do
#for induction in gRNA RNP;do
#bsub -q new.long -J ${target}${induction} -e ${target}${induction}.e -o ${target}${induction}.o -R rusage[mem=16000] 
#./optimize_model_backbone.R  -t ${target} -i $induction -m modelDSBs1i1_fullimpreciseDSB -e E_errorsfromunbroken -n 10000 -o ${target}_${induction}_optimization.tsv
outputf=$( basename $inputf .tsv )
./calculate_CI.R -i ${inputf} -o ${outputf}.CI -n ${NMAX}
#done
#done
