#!/bin/bash
module load R/4.0.4

inputf=$1 #likelihood file
data_file=$2
ERRORMATRIX=$3
NMAX=$4
MYMODEL=$5
errorflag=$6
#for target in CRTISO Psy1 CRTISO_49and50bp;do
#for induction in gRNA RNP;do
#bsub -q new.long -J ${target}${induction} -e ${target}${induction}.e -o ${target}${induction}.o -R rusage[mem=16000] 
#./optimize_model_backbone.R  -t ${target} -i $induction -m modelDSBs1i1_fullimpreciseDSB -e E_errorsfromunbroken -n 10000 -o ${target}_${induction}_optimization.tsv
outputf=$( basename $inputf .tsv )
outputf2=$( dirname $inputf )
#if [ ! -f ${outputf2}/${outputf}_n${NMAX}.CI ];then
./calculate_CI.v1.R -i ${inputf} -o ${outputf2}/${outputf}_n${NMAX}.CI -d $data_file -n ${NMAX} -E ${ERRORMATRIX} -m $MYMODEL -l $errorflag
#fi
#done
#done
