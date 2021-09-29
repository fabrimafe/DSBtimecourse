#!/bin/bash
module load R/4.0.4

mymodel=$1
target=$2
induction=$3
estimaterr=$4
#for target in CRTISO Psy1 CRTISO_49and50bp;do
#for induction in gRNA RNP;do
#bsub -q new.long -J ${target}${induction} -e ${target}${induction}.e -o ${target}${induction}.o -R rusage[mem=16000] 
#./optimize_model_backbone.R  -t ${target} -i $induction -m modelDSBs1i1_fullimpreciseDSB -e E_errorsfromunbroken -n 10000 -o ${target}_${induction}_optimization.tsv
./optimize_model_backbone.R -d $estimaterr -t ${target} -i $induction -m ${mymodel} -e E_errorsfromunbroken -n 1000 -o ${mymodel}_${target}_${induction}_ester${estimaterr}_optimization.tsv
#done
#done
