#!/bin/bash
module load R/4.0.4

inputf=$1 #likelihood file
data_file=$2
ERRORMATRIX=$3
model=$4
NMAX=$5
OUTPUT=$6
likfun=$7
INDUCTIONCURVE=$8
TIMERESOLUTION=$9
CALCULATEFLOW=${10}
outputf=$( basename $inputf .tsv )
outputf2=$( dirname $inputf )
#if [ ! -f ${outputf2}/${outputf}_plot.flow.tab ];then
./plot_bootstraps.R -i ${inputf} -d $data_file -n ${NMAX} -E ${ERRORMATRIX} -m ${model} -o ${OUTPUT} -l $likfun -z ${INDUCTIONCURVE} -r ${TIMERESOLUTION} -w $CALCULATEFLOW #{outputf2}/${outputf}
#fi
#done
#done
