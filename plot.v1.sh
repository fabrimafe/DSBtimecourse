#!/bin/bash
module load R/4.0.4

inputf=$1 #likelihood file
data_file=$2
ERRORMATRIX=$3
model=$4
NMAX=$5
OUTPUT=$6
outputf=$( basename $inputf .tsv )
outputf2=$( dirname $inputf )
#if [ ! -f ${outputf2}/${outputf}_plot.flow.tab ];then
./plot_CI.v1.R -i ${outputf2}/${outputf}.CI.RData -d $data_file -n ${NMAX} -E ${ERRORMATRIX} -m ${model} -o ${OUTPUT} #{outputf2}/${outputf}
#fi
#done
#done
