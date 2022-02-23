#!/bin/bash
#module load R/4.0.4

#new README created 20220106 for new induction curves

#current issue: 
runs with 2 params induction curve had 0 r11, even with delay. To understand why, plotted trajectories. Fit badly. Induction clearly different.
Problem might be that it cannot go up sharply, so no need for repair to keep cut high. Explore function to see if it can actually go up fast.
If real and no bug, perhaps for 2 I should try delay. Otherwise I stick to 3.
From plots, it looks real. Now
####---OPTIMIZATION---#############################
for MYDELAY in ''  _mydelay0.25;do #''
for MYINDUCTION in RNP; do # gRNA;do
for MYINDUCTION_CURVE in 2 3;do # 3 4; do # 2 4
for mymodel in modelDSBs1i1_mini modelDSBs1i1_3x4 modelDSBs1i1_nok12; do 
for MYTARGET in CRTISO.cleanedandnov2021 CRTISO.nov2021.49and50 Psy1 PhyB2.2.nov2021 CRTISO.nov2021; do 
MYTARGETM=$( echo $MYTARGET | sed 's/0h//' | sed 's/m$//' | sed 's/72h//' | sed 's/_mydelay0.25//' )
ERRORMATRIX=~/workspace/daniela/error_matrices/error_matrix4_${MYTARGETM}_errorsfromunbroken.tsv
TIMECOURSE=~/workspace/daniela/input_datasets/timecourse_${MYINDUCTION}_${MYTARGET}${MYDELAY}.txt
if [ $mymodel == "model5i1" ];then
ERRORMATRIX=~/workspace/daniela/error_matrices/error_matrix3_${MYTARGETM}_errorsfromunbroken.tsv
TIMECOURSE=~/workspace/daniela/input_datasets/timecourse_${MYINDUCTION}_${MYTARGET}${MYDELAY}_3states.txt
fi
OUTPUT=~/workspace/daniela/resultsv4/results_${mymodel}_${MYINDUCTION}_ind.c${MYINDUCTION_CURVE}_${MYTARGET}${MYDELAY}
likfun=0; if [ $mymodel == "modelDSBs1i1_3x4" ];then likfun=3; fi
echo $OUTPUTA; NAMESIMS=v4 #${MYTARGET};
MEMS=18000; NITER=1250;
for IITER in a b c d;do
if [ ! -f ${OUTPUT}.${IITER} ];then bsub -q new-long -J ${NAMESIMS} -e ~/mylogs/${NAMESIMS}-%J.e -o ~/mylogs/${NAMESIMS}-%J.o -R rusage[mem=${MEMS}] ./launch_optimization_cas9_v1b.sh $TIMECOURSE $ERRORMATRIX ${OUTPUT}.${IITER} ${mymodel} ${NITER} $likfun ${MYINDUCTION_CURVE}; fi
done;done;done;done;done;done

for MYDELAY in '';do # _mydelay0.25;do
for MYINDUCTION in RNP; do # gRNA;do
for MYINDUCTION_CURVE in 2 3;do # 3 4; do # 2 4
for mymodel in modelDSBs1i1_mini modelDSBs1i1_3x4 modelDSBs1i1_nok12; do 
for MYTARGET in CRTISO.cleanedandnov2021;do
for bootstraptype in stratifiedbootstraps;do #stationarybootstraps MEbootstraps; do #MEbootstraps stationarybootstraps
for ip in `seq 1 100`;do
MYTARGETM=$( echo $MYTARGET | sed 's/0h//' | sed 's/m$//' | sed 's/72h//' )
ERRORMATRIX=~/workspace/daniela/error_matrices/error_matrix4_${MYTARGETM}_errorsfromunbroken.tsv
TIMECOURSEPATH=~/workspace/daniela/input_datasets/${bootstraptype}/timecourse_${MYINDUCTION}_${MYTARGET}${MYDELAY}
TIMECOURSE=${TIMECOURSEPATH}/timecourse.bs${ip}.txt
OUTPUTA=~/workspace/daniela/resultsv4/${bootstraptype}/timecourse_${MYINDUCTION}_${MYTARGET}${MYDELAY}/${mymodel}
OUTPUT=~/workspace/daniela/resultsv4/${bootstraptype}/timecourse_${MYINDUCTION}_${MYTARGET}${MYDELAY}/${mymodel}/results_${mymodel}_${MYINDUCTION}_ind.c${MYINDUCTION_CURVE}_${MYTARGET}${MYDELAY}
mkdir -p $OUTPUT;
OUTPUTTOT=${OUTPUT}/results.tot
OUTPUTA=${OUTPUT}/results_${ip}.txt
likfun=0; if [ $mymodel == "modelDSBs1i1_3x4" ];then likfun=3; fi
echo $OUTPUTA; NAMESIMS=MEboot #${MYTARGET};
MEMS=18000; NITER=1250;
#if [ -f $OUTPUTA ];then rm $OUTPUTA;fi
#bsub -q new-long -J ${NAMESIMS} -e ~/mylogs/${NAMESIMS}-%J.e -o ~/mylogs/${NAMESIMS}-%J.o -R rusage[mem=${MEMS}] ./launch_optimization_cas9_v1b.sh $TIMECOURSE $ERRORMATRIX ${OUTPUTA} ${mymodel} ${NITER} $likfun ${MYINDUCTION_CURVE} 
if [ $ip -eq 1 ]; then rm ${OUTPUTTOT};fi
cat ${OUTPUTA} | awk -v var1=$ip -v OFS='\t' '{if (NR==1){print $0,"bootstrap"} else {print $0,var1}}' >> ${OUTPUTTOT}
done;done;done;done;done;done;done



#CALCULATE CI
#----likelihood based
for MYN in 20000 50000 500000;do # 500000;do # 150000 2000000
for MYDELAY in '' _mydelay0.25;do #''
for MYINDUCTION in RNP;do #gRNA
for MYINDUCTION_CURVE in 2 3; do #4; do # gRNA;do
for mymodel in modelDSBs1i1_mini modelDSBs1i1_3x4 modelDSBs1i1_nok12; do #modelDSBs1i1_realimprecise model5i1 modelDSBs1i1_realimpnor11;do
for MYTARGET in CRTISO.nov2021.49and50 Psy1 PhyB2.2.nov2021 CRTISO.nov2021;do 
MYTARGETM=$( echo $MYTARGET | sed 's/0h//' | sed 's/m$//' | sed 's/72h//' )
ERRORMATRIX=~/workspace/daniela/error_matrices/error_matrix4_${MYTARGETM}_errorsfromunbroken.tsv
TIMECOURSE=~/workspace/daniela/input_datasets/timecourse_${MYINDUCTION}_${MYTARGET}${MYDELAY}.txt
if [ $mymodel == "model5i1" ];then
ERRORMATRIX=~/workspace/daniela/error_matrices/error_matrix3_${MYTARGETM}_errorsfromunbroken.tsv
TIMECOURSE=~/workspace/daniela/input_datasets/timecourse_${MYINDUCTION}_${MYTARGET}${MYDELAY}_3states.txt
fi
OUTPUT=~/workspace/daniela/resultsv4/results_${mymodel}_${MYINDUCTION}_ind.c${MYINDUCTION_CURVE}_${MYTARGET}${MYDELAY}
OUTPUTALL=${OUTPUT}.all.tsv
likfun=0; if [ $mymodel == "modelDSBs1i1_3x4" ];then likfun=3; fi
if [ -f ${OUTPUT}.a ];then
	echo $MYN $TIMECOURSE
	cat ${OUTPUT}.a | head -1 > temp.tab
	for IITER in a b c d;do
		cat ${OUTPUT}.${IITER} | grep -v k11 >> temp.tab
	done
	cat temp.tab | head -5001 > $OUTPUTALL
	bsub -q new-long -J CI -e ~/mylogs2/confint.e%J -o ~/mylogs2/confint.o%J -R rusage[mem=18000] ./calculate_CI_v1.sh $OUTPUTALL $TIMECOURSE $ERRORMATRIX $MYN $mymodel $likfun $MYINDUCTION_CURVE
fi
done;done;done;done;done;done


#in plot use z!=1 not to normalize flow on induction curve (with new induction function)
TIMERESOLUTION=0.002
echo "aaa" | awk '{printf "max\tCIlow\tCIhigh\trate\ttarget\tinduction\tmodel\terror\ttimecourse\n"}' > ~/workspace/daniela/resultsv2/table.CI
for MYDELAY in '' _mydelay0.25;do #''
for MYINDUCTION in RNP;do #gRNA
for MYINDUCTION_CURVE in 2 3; do # gRNA;do
for mymodel in modelDSBs1i1_mini modelDSBs1i1_3x4 modelDSBs1i1_nok12;do 
for MYTARGET in Psy1 PhyB2.2.nov2021 CRTISO.nov2021 CRTISO.cleanedandnov2021 CRTISO.nov2021.49and50 ;do 
MYTARGETM=$( echo $MYTARGET | sed 's/0h//' | sed 's/m$//' | sed 's/72h//' )
ERRORMATRIX=~/workspace/daniela/error_matrices/error_matrix4_${MYTARGETM}_errorsfromunbroken.tsv
TIMECOURSE=~/workspace/daniela/input_datasets/timecourse_${MYINDUCTION}_${MYTARGET}${MYDELAY}.txt
likfun=0; if [ $mymodel == "modelDSBs1i1_3x4" ];then likfun=3; fi
if [ $mymodel == "model5i1" ];then
ERRORMATRIX=~/workspace/daniela/error_matrices/error_matrix3_${MYTARGETM}_errorsfromunbroken.tsv
TIMECOURSE=~/workspace/daniela/input_datasets/timecourse_${MYINDUCTION}_${MYTARGET}${MYDELAY}_3states.txt
fi
OUTPUT=~/workspace/daniela/resultsv4/results_${mymodel}_${MYINDUCTION}_ind.c${MYINDUCTION_CURVE}_${MYTARGET}${MYDELAY}
if [ -f ${OUTPUT}.a ];then
        #ls ~/workspace/daniela/resultsv4/results_${mymodel}_${MYINDUCTION}_${MYTARGET}${MYDELAY}.all_n*
        myfile="notexist"
        if [ -f ${OUTPUT}.all_n20000.CI.RData ];then myfile=${OUTPUT}.all_n20000.CI.RData;fi
        if [ -f ${OUTPUT}.all_n50000.CI.RData ];then myfile=${OUTPUT}.all_n50000.CI.RData;fi
        if [ -f ${OUTPUT}.all_n150000.CI.RData ];then myfile=${OUTPUT}.all_n150000.CI.RData;fi
        if [ -f ${OUTPUT}.all_n500000.CI.RData ];then myfile=${OUTPUT}.all_n500000.CI.RData;fi
        ERRORMATRIXM=$( basename $ERRORMATRIX .tsv )
        TIMECOURSEM=$( basename $TIMECOURSE .txt )
        if [ -f ${myfile} ];then
                #cat ${myfile}.CI | awk -v var1=$MYTARGET -v var2=$MYINDUCTION -v var3=$mymodel -v var4=$ERRORMATRIXM -v var5=$TIMECOURSEM '{print $0,var1,var2,var3,var4,var5}' | grep -v max >> ~/workspace/daniela/resultsv2/table.CI
                bsub -q new-short -J plotshort -e ~/mylogs/plot.e%J -o ~/mylogs/plot.o%J -R rusage[mem=16000] ./plot.v1.sh $myfile $TIMECOURSE $ERRORMATRIX $mymodel 1000 ~/workspace/daniela/resultsv4/plots/plot_${mymodel}_${MYINDUCTION}_ind.c${MYINDUCTION_CURVE}_${MYTARGET}${MYDELAY} $likfun ${MYINDUCTION_CURVE} ${TIMERESOLUTION}
        fi
fi
done;done;done;done;done


