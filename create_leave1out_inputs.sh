#!/bin/bash
INPUTFILE=$1
OUTPATH=$2
NLINES=$( cat $INPUTFILE | wc -l )
for i in `seq 2 $NLINES`;do
NAMELINE=$(( $i - 1 ))
cat $INPUTFILE | awk -v NOTIN=${i} '{if (NR!=NOTIN){print}}' > ${OUTPATH}_jk${NAMELINE}.txt
done
