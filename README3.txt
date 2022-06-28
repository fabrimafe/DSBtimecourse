#!/bin/bash
#module load R/4.0.4

#new README created 20220106 for new induction curves

############################################################################
#############################---OPTIMIZATION---#############################
############################################################################

#-------------------------------SINGLE RUNS---------------------------------

for MYDELAY in '';do # _mydelay0.25;do #''
for MYINDUCTION in RNP;do #RNP; do # gRNA;do #Psy1 PhyB2 CRTISO.cleaned for gRNA
for MYINDUCTION_CURVE in 3;do # 2
for mymodel in modelDSBs1i1_3x4nor11;do #modelDSBs1i1_realimprecise;do #modelDSBs1i1_3x4 modelDSBs1i1_nok12;do #modelDSBs1i1_mini
for MYTARGET in Psy1 CRTISO.nov2021.49and50 CRTISO.nov2021 PhyB2.2.nov2021 Psy1_allb PhyB2.2_allb;do #PhyB2 CRTISO.cleaned CRTISO_allb CRTISO.49and50bp_allb Psy1_allb PhyB2.2_allb Psy1 PhyB2.2.nov2021 PhyB2.2_R2_Feb2022 Psy1_R2_Feb2022 CRTISO.nov2021.49and50 CRTISO.nov2021;do  
MYTARGETM=$( echo $MYTARGET | sed 's/_allb//' | sed 's/0h//' | sed 's/m$//' | sed 's/72h//' | sed 's/_mydelay0.25//' | sed 's/_R2_Feb2022//' )
ERRORMATRIX=~/workspace/daniela/error_matrices/error_matrix4_${MYTARGETM}_errorsfromunbroken.tsv
TIMECOURSE=~/workspace/daniela/input_datasets/timecourse_${MYINDUCTION}_${MYTARGET}${MYDELAY}.txt
if [ $mymodel == "model5i1" ];then
ERRORMATRIX=~/workspace/daniela/error_matrices/error_matrix3_${MYTARGETM}_errorsfromunbroken.tsv
TIMECOURSE=~/workspace/daniela/input_datasets/timecourse_${MYINDUCTION}_${MYTARGET}${MYDELAY}_3states.txt
fi
OUTPUT=~/workspace/daniela/resultsv5/results_${mymodel}_${MYINDUCTION}_ind.c${MYINDUCTION_CURVE}_${MYTARGET}${MYDELAY}
likfun=0; #if [ $mymodel == "modelDSBs1i1_3x4" ];then likfun=3; fi
echo $OUTPUT; NAMESIMS=v4 #${MYTARGET};
MEMS=18000; NITER=1250;
for IITER in a b c d;do
if [ ! -f ${OUTPUT}.${IITER} ];then bsub -q new-long -J ${NAMESIMS} -e ~/mylogs/${NAMESIMS}-%J.e -o ~/mylogs/${NAMESIMS}-%J.o -R rusage[mem=${MEMS}] ./launch_optimization_cas9_v1b.sh $TIMECOURSE $ERRORMATRIX ${OUTPUT}.${IITER} ${mymodel} ${NITER} $likfun ${MYINDUCTION_CURVE}; fi
done;done;done;done;done;done

#-----------------------------------SINGLE RUNS ALL READS----------------------------------

for MYDELAY in '';do # _mydelay0.25;do #''
for MYINDUCTION in RNP; do # gRNA;do
for MYINDUCTION_CURVE in 2;do
for mymodel in modelDSBs1i1_3x4.bytarget modelDSBs1i1_mini.bytarget modelDSBs1i1_nok12.bytarget; do 
for MYTARGET in Psy1_all PhyB2.2_all CRTISO_all CRTISO.49and50bp_all;do 
MYTARGETM=$( echo $MYTARGET | sed 's/_all//' | sed 's/0h//' | sed 's/m$//' | sed 's/72h//' | sed 's/_mydelay0.25//' | sed 's/_R2_Feb2022//' )
ERRORMATRIX=~/workspace/daniela/error_matrices/error_matrix4_${MYTARGETM}_errorsfromunbroken.tsv
TIMECOURSE=~/workspace/daniela/input_datasets/timecourse_${MYINDUCTION}_${MYTARGET}${MYDELAY}.txt
if [ $mymodel == "model5i1" ];then
ERRORMATRIX=~/workspace/daniela/error_matrices/error_matrix3_${MYTARGETM}_errorsfromunbroken.tsv
TIMECOURSE=~/workspace/daniela/input_datasets/timecourse_${MYINDUCTION}_${MYTARGET}${MYDELAY}_3states.txt
fi
OUTPUT=~/workspace/daniela/resultsv5/results_${mymodel}_${MYINDUCTION}_ind.c${MYINDUCTION_CURVE}_${MYTARGET}${MYDELAY}
likfun=0; #if [ $mymodel == "modelDSBs1i1_3x4" || $mymodel == "modelDSBs1i1_3x4.bytarget" ];then likfun=3; fi
echo $OUTPUTA; NAMESIMS=v5 #${MYTARGET};
MEMS=18000; NITER=1250;
for IITER in a b c d;do
if [ ! -f ${OUTPUT}.${IITER} ];then bsub -q new-long -J ${NAMESIMS} -e ~/mylogs/${NAMESIMS}-%J.e -o ~/mylogs/${NAMESIMS}-%J.o -R rusage[mem=${MEMS}] ./launch_optimization_cas9_v1b.sh $TIMECOURSE $ERRORMATRIX ${OUTPUT}.${IITER} ${mymodel} ${NITER} $likfun ${MYINDUCTION_CURVE}; fi
done;done;done;done;done;done

#------------------------------------BOOTSTRAPS----------------------------------

#stratifiedbootstraps2 goes with all and modelDSBs1i1_mini.bytarget; stratifiedbootstraps with allb
for MYDELAY in '';do # _mydelay0.25;do
for MYINDUCTION in RNP; do # gRNA;do #Psy1 PhyB2 CRTISO.cleaned x gRNA
for MYINDUCTION_CURVE in 3;do # 2
#for mymodel in modelDSBs1i1_3x4.bytarget;do #modelDSBs1i1_nok12.bytarget;do #modelDSBs1i1_3x4.bytarget modelDSBs1i1_mini.bytarget;do 
#for MYTARGET in Psy1_all PhyB2.2_all CRTISO_all;do # CRTISO.49and50bp_all;do #Psy1_allb PhyB2.2_allb;do #CRTISO_allb Psy1_allb PhyB2.2_allb;do #CRTISO.cleanedandnov2021;do
#for bootstraptype in stratifiedbootstraps2;do 
for mymodel in modelDSBs1i1_3x4nor11;do #modelDSBs1i1_realimprecise;do #modelDSBs1i1_3x4 modelDSBs1i1_nok12 modelDSBs1i1_mini; do # modelDSBs1i1_3x4;do #modelDSBs1i1_mini modelDSBs1i1_nok12 modelDSBs1i1_mini modelDSBs1i1_3x4 modelDSBs1i1_nok12; do 
for MYTARGET in Psy1 CRTISO.nov2021.49and50 CRTISO.nov2021 PhyB2.2.nov2021 Psy1_allb PhyB2.2_allb;do #CRTISO.cleaned PhyB2.2.nov2021 PhyB2.2_R2_Feb2022 Psy1_R2_Feb2022 CRTISO.nov2021.49and50 CRTISO.nov2021 CRTISO_allb CRTISO.49and50bp_allb Psy1_allb PhyB2.2_allb;do
for bootstraptype in stratifiedbootstraps;do #stratifiedbootstraps;do #stationarybootstraps MEbootstraps; do #MEbootstraps stationarybootstraps
for ip in `seq 1 100`;do
MYTARGETM=$( echo $MYTARGET | sed 's/_allb//' |  sed 's/_all//' | sed 's/0h//' | sed 's/m$//' | sed 's/72h//' | sed 's/_R2_Feb2022//' )
ERRORMATRIX=~/workspace/daniela/error_matrices/error_matrix4_${MYTARGETM}_errorsfromunbroken.tsv
TIMECOURSEPATH=~/workspace/daniela/input_datasets/${bootstraptype}/timecourse_${MYINDUCTION}_${MYTARGET}${MYDELAY}
TIMECOURSE=${TIMECOURSEPATH}/timecourse.bs${ip}.txt
OUTPUTA=~/workspace/daniela/resultsv5/${bootstraptype}/timecourse_${MYINDUCTION}_${MYTARGET}${MYDELAY}/${mymodel}
OUTPUT=~/workspace/daniela/resultsv5/${bootstraptype}/timecourse_${MYINDUCTION}_${MYTARGET}${MYDELAY}/${mymodel}/results_${mymodel}_${MYINDUCTION}_ind.c${MYINDUCTION_CURVE}_${MYTARGET}${MYDELAY}
mkdir -p $OUTPUT;
OUTPUTTOT=${OUTPUT}/results.tot
OUTPUTA=${OUTPUT}/results_${ip}.txt
likfun=0; #if [ $mymodel == "modelDSBs1i1_3x4" ];then likfun=3; fi
echo $OUTPUTA; NAMESIMS=boot3 #MEboot #${MYTARGET};
MEMS=18000; NITER=1250;
if [ ! -f $OUTPUTA ];then 
rm $OUTPUTA
bsub -q new-long -J ${NAMESIMS} -e ~/mylogs/${NAMESIMS}-%J.e -o ~/mylogs/${NAMESIMS}-%J.o -R rusage[mem=${MEMS}] ./launch_optimization_cas9_v1b.sh $TIMECOURSE $ERRORMATRIX ${OUTPUTA} ${mymodel} ${NITER} $likfun ${MYINDUCTION_CURVE} 
fi
#if [ $ip -eq 1 ]; then rm ${OUTPUTTOT};fi
#cat ${OUTPUTA} | awk -v var1=$ip -v OFS='\t' '{if (NR==1){print $0,"bootstrap"} else {print $0,var1}}' >> ${OUTPUTTOT}
done;done;done;done;done;done;done


#########################################################################################
################################### CALCULATE CI ########################################
#########################################################################################

#-------------------------------------SINGLE RUNS----------------------------------------

for MYN in 50000 500000;do # 150000 2000000
for MYDELAY in '';do # _mydelay0.25;do #''
for MYINDUCTION in RNP;do #RNP
for MYINDUCTION_CURVE in 3; do #3 4; do # gRNA;do
#for mymodel in modelDSBs1i1_nok12.bytarget;do #modelDSBs1i1_mini.bytarget modelDSBs1i1_3x4.bytarget modelDSBs1i1_nok12.bytarget;do 
#for MYTARGET in CRTISO_all CRTISO.49and50bp_all Psy1_all PhyB2.2_all;do # CRTISO_all;do #PhyB2.2_R2_Feb2022 Psy1_R2_Feb2022 Psy1_allb PhyB2.2_allb CRITSO_allb CRTISO Psy1;do #CRTISO.nov2021.49and50 Psy1 PhyB2.2.nov2021 CRTISO.nov2021;do 
for mymodel in modelDSBs1i1_3x4nor11 modelDSBs1i1_realimprecise;do #modelDSBs1i1_nok12 modelDSBs1i1_mini modelDSBs1i1_3x4; do #modelDSBs1i1_realimprecise model5i1 modelDSBs1i1_realimpnor11;do
for MYTARGET in Psy1 CRTISO.nov2021.49and50 CRTISO.nov2021 PhyB2.2.nov2021 Psy1_allb PhyB2.2_allb;do #CRTISO.cleaned PhyB2.2.nov2021 PhyB2.2_R2_Feb2022 Psy1_R2_Feb2022 CRTISO.nov2021.49and50 CRTISO.nov2021 CRTISO_allb CRTISO.49and50bp_allb Psy1_allb PhyB2.2_allb;do
MYTARGETM=$( echo $MYTARGET | sed 's/_allb//' |  sed 's/_all//' | sed 's/_R2_Feb2022//' | sed 's/0h//' | sed 's/m$//' | sed 's/72h//' )
ERRORMATRIX=~/workspace/daniela/error_matrices/error_matrix4_${MYTARGETM}_errorsfromunbroken.tsv
TIMECOURSE=~/workspace/daniela/input_datasets/timecourse_${MYINDUCTION}_${MYTARGET}${MYDELAY}.txt
if [ $mymodel == "model5i1" ];then
ERRORMATRIX=~/workspace/daniela/error_matrices/error_matrix3_${MYTARGETM}_errorsfromunbroken.tsv
TIMECOURSE=~/workspace/daniela/input_datasets/timecourse_${MYINDUCTION}_${MYTARGET}${MYDELAY}_3states.txt
fi
OUTPUT=~/workspace/daniela/resultsv5/results_${mymodel}_${MYINDUCTION}_ind.c${MYINDUCTION_CURVE}_${MYTARGET}${MYDELAY}
OUTPUTALL=${OUTPUT}.all.tsv
likfun=0; if [ $mymodel == "modelDSBs1i1_3x4" ];then likfun=3; fi
if [ -f ${OUTPUT}.a ];then
	echo $MYN $TIMECOURSE
	cat ${OUTPUT}.a | head -1 > temp.tab
	for IITER in a b c d;do
		cat ${OUTPUT}.${IITER} | grep -v k11 >> temp.tab
	done
	cat temp.tab | head -5001 > $OUTPUTALL
	bsub -q new-long -J CIs -e ~/mylogs2/confint.e%J -o ~/mylogs2/confint.o%J -R rusage[mem=18000] ./calculate_CI_v1.sh $OUTPUTALL $TIMECOURSE $ERRORMATRIX $MYN $mymodel $likfun $MYINDUCTION_CURVE
fi
done;done;done;done;done;done

#------------------------------------BOOTSTRAPS----------------------------------

#stratifiedbootstraps2 goes with all and modelDSBs1i1_mini.bytarget; stratifiedbootstraps with allb
for MYN in 50000;do #20000 150000 2000000
for MYDELAY in '';do # _mydelay0.25;do
for MYINDUCTION in RNP;do # RNP;do
for MYINDUCTION_CURVE in 3;do # 3 4; do # 2 4
#for mymodel in modelDSBs1i1_nok12.bytarget;do #modelDSBs1i1_mini.bytarget modelDSBs1i1_3x4.bytarget;do 
#for MYTARGET in Psy1_all PhyB2.2_all CRTISO_all CRTISO.49and50bp_all;do
#for bootstraptype in stratifiedbootstraps2;do #stratifiedbootstraps2;do #stationarybootstraps MEbootstraps; do #MEbootstraps stationarybootstraps
for mymodel in modelDSBs1i1_3x4nor11;do #modelDSBs1i1_realimprecise;do #modelDSBs1i1_3x4 modelDSBs1i1_mini modelDSBs1i1_nok12;do # modelDSBs1i1_3x4;do #modelDSBs1i1_mini modelDSBs1i1_nok12;do 
for MYTARGET in Psy1 CRTISO.nov2021.49and50 CRTISO.nov2021 PhyB2.2.nov2021 Psy1_allb PhyB2.2_allb;do
for bootstraptype in stratifiedbootstraps;do #stratifiedbootstraps2;do #stationarybootstraps MEbootstraps; do #MEbootstraps stationarybootstraps
for ip in `seq 1 100`;do
MYTARGETM=$( echo $MYTARGET | sed 's/_allb//' |  sed 's/_all//' | sed 's/0h//' | sed 's/m$//' | sed 's/72h//' | sed 's/_R2_Feb2022//' )
ERRORMATRIX=~/workspace/daniela/error_matrices/error_matrix4_${MYTARGETM}_errorsfromunbroken.tsv
TIMECOURSEPATH=~/workspace/daniela/input_datasets/${bootstraptype}/timecourse_${MYINDUCTION}_${MYTARGET}${MYDELAY}
TIMECOURSE=${TIMECOURSEPATH}/timecourse.bs${ip}.txt
OUTPUTA=~/workspace/daniela/resultsv5/${bootstraptype}/timecourse_${MYINDUCTION}_${MYTARGET}${MYDELAY}/${mymodel}
OUTPUT=~/workspace/daniela/resultsv5/${bootstraptype}/timecourse_${MYINDUCTION}_${MYTARGET}${MYDELAY}/${mymodel}/results_${mymodel}_${MYINDUCTION}_ind.c${MYINDUCTION_CURVE}_${MYTARGET}${MYDELAY}
mkdir -p $OUTPUT;
OUTPUTTOT=${OUTPUT}/results.tot
OUTPUTA=${OUTPUT}/results_${ip}.txt
likfun=0; #if [ $mymodel == "modelDSBs1i1_3x4" ];then likfun=3; fi
echo $OUTPUTA; NAMESIMS=boot2 #MEboot #${MYTARGET};
MEMS=18000; NITER=1250;
cat ${OUTPUTA} | head -1 > ${OUTPUTA}.parsed.tab 
cat ${OUTPUTA} | grep -v k11 >> ${OUTPUTA}.parsed.tab
bsub -q new-medium -J CI -e ~/mylogs3/confint.e%J -o ~/mylogs3/confint.o%J -R rusage[mem=18000] ./calculate_CI_v1.sh ${OUTPUTA}.parsed.tab $TIMECOURSE $ERRORMATRIX $MYN $mymodel $likfun $MYINDUCTION_CURVE
done;done;done;done;done;done;done;done


#########################################################################################
################################### PLOTTING / FLOW #####################################
#########################################################################################

#-------------------------------------SINGLE RUNS----------------------------------------

#in plot use z!=1 not to normalize flow on induction curve (with new induction function)
TIMERESOLUTION=0.002
#echo "aaa" | awk '{printf "max\tCIlow\tCIhigh\trate\ttarget\tinduction\tmodel\terror\ttimecourse\n"}' > ~/workspace/daniela/resultsv2/table.CI
for MYDELAY in '';do # _mydelay0.25;do #''
for MYINDUCTION in RNP;do #RNP
for MYINDUCTION_CURVE in 3; do #3 gRNA;do
#for mymodel in modelDSBs1i1_mini.bytarget modelDSBs1i1_3x4.bytarget modelDSBs1i1_nok12.bytarget;do
#for MYTARGET in Psy1_all CRTISO_all CRTISO.49and50bp_all PhyB2.2_all;do #PhyB2.2_R2_Feb2022 Psy1_R2_Feb2022 Psy1_allb PhyB2.2_allb CRITSO_allb; do #Psy1 PhyB2.2.nov2021 CRTISO.nov2021 CRTISO.cleanedandnov2021 CRTISO.nov2021.49and50 ;do 
#modelDSBs1i1_mini modelDSBs1i1_3x4 modelDSBs1i1_nok12;do 
for mymodel in modelDSBs1i1_realimprecise;do #modelDSBs1i1_3x4nor11;do #modelDSBs1i1_realimprecise;do #modelDSBs1i1_3x4 modelDSBs1i1_nok12 modelDSBs1i1_mini;do
for MYTARGET in Psy1 CRTISO.nov2021.49and50 CRTISO.nov2021 PhyB2.2.nov2021;do #PhyB2 CRTISO.cleaned CRTISO_allb CRTISO.49and50bp_allb Psy1_allb PhyB2.2_allb Psy1 PhyB2.2.nov2021 PhyB2.2_R2_Feb2022 Psy1_R2_Feb2022 CRTISO.nov2021.49and50 CRTISO.nov2021;do 
MYTARGETM=$( echo $MYTARGET | sed 's/_allb//' |  sed 's/_all//' | sed 's/_R2_Feb2022//' | sed 's/0h//' | sed 's/m$//' | sed 's/72h//' )
ERRORMATRIX=~/workspace/daniela/error_matrices/error_matrix4_${MYTARGETM}_errorsfromunbroken.tsv
TIMECOURSE=~/workspace/daniela/input_datasets/timecourse_${MYINDUCTION}_${MYTARGET}${MYDELAY}.txt
likfun=0; if [ $mymodel == "modelDSBs1i1_3x4" ];then likfun=3; fi
if [ $mymodel == "model5i1" ];then
ERRORMATRIX=~/workspace/daniela/error_matrices/error_matrix3_${MYTARGETM}_errorsfromunbroken.tsv
TIMECOURSE=~/workspace/daniela/input_datasets/timecourse_${MYINDUCTION}_${MYTARGET}${MYDELAY}_3states.txt
fi
OUTPUT=~/workspace/daniela/resultsv5/results_${mymodel}_${MYINDUCTION}_ind.c${MYINDUCTION_CURVE}_${MYTARGET}${MYDELAY}
if [ -f ${OUTPUT}.a ];then
        #ls ~/workspace/daniela/resultsv4/results_${mymodel}_${MYINDUCTION}_${MYTARGET}${MYDELAY}.all_n*
        myfile="notexist"
        if [ -f ${OUTPUT}.all_n20000.CI.RData ];then myfile=${OUTPUT}.all_n20000.CI.RData;fi
        if [ -f ${OUTPUT}.all_n50000.CI.RData ];then myfile=${OUTPUT}.all_n50000.CI.RData;fi
        if [ -f ${OUTPUT}.all_n150000.CI.RData ];then myfile=${OUTPUT}.all_n150000.CI.RData;fi
#        if [ -f ${OUTPUT}.all_n500000.CI.RData ];then myfile=${OUTPUT}.all_n500000.CI.RData;fi
        ERRORMATRIXM=$( basename $ERRORMATRIX .tsv )
        TIMECOURSEM=$( basename $TIMECOURSE .txt )
        if [ -f ${myfile} ];then
		echo $myfile
                # cat ${myfile}.CI | awk -v var1=$MYTARGET -v var2=$MYINDUCTION -v var3=$mymodel -v var4=$ERRORMATRIXM -v var5=$TIMECOURSEM '{print $0,var1,var2,var3,var4,var5}' | grep -v max >> ~/workspace/daniela/resultsv5/table.CI
#               bsub -q new-short -J plotshort -e ~/mylogs/plot.e%J -o ~/mylogs/plot.o%J -R rusage[mem=16000] ./plot.v1.sh $myfile $TIMECOURSE $ERRORMATRIX $mymodel 1000 ~/workspace/daniela/resultsv4/plots/plot_${mymodel}_${MYINDUCTION}_ind.c${MYINDUCTION_CURVE}_${MYTARGET}${MYDELAY} $likfun ${MYINDUCTION_CURVE} ${TIMERESOLUTION}
	bsub -q new-short -J plotshort -e ~/mylogs/plot.e%J -o ~/mylogs/plot.o%J -R rusage[mem=16000] ./plot.v2.sh $myfile $TIMECOURSE $ERRORMATRIX $mymodel 1000 ~/workspace/daniela/resultsv5/plots/plot_${mymodel}_${MYINDUCTION}_ind.c${MYINDUCTION_CURVE}_${MYTARGET}${MYDELAY} $likfun ${MYINDUCTION_CURVE} ${TIMERESOLUTION} 1
        fi
fi
done;done;done;done;done

#-------------------------------------BOOTSTRAPS----------------------------------------
#/home/labs/alevy/fabrizio/workspace/daniela/resultsv4/stratifiedbootstraps/timecourse_RNP_PhyB2.2_allb/modelDSBs1i1_mini/results_modelDSBs1i1_mini_RNP_ind.c2_PhyB2.2_allb/



TIMERESOLUTION=0.002
#stratifiedbootstraps2 goes with all and modelDSBs1i1_mini.bytarget; stratifiedbootstraps with allb
for MYDELAY in '';do # _mydelay0.25;do
for MYINDUCTION in RNP; do # RNP;do
for MYINDUCTION_CURVE in 3;do # 3 4; do # 2 4
#TWO INDUCTION CURVES
#for mymodel in modelDSBs1i1_nok12.bytarget;do # modelDSBs1i1_mini.bytarget modelDSBs1i1_3x4.bytarget;do  
#for MYTARGET in CRTISO.49and50bp_all CRTISO_all Psy1_all PhyB2.2_all;do 
#for bootstraptype in stratifiedbootstraps2;do 
#SINGLE INDUCTION CURVE
#for mymodel in modelDSBs1i1_realimprecise modelDSBs1i1_3x4 modelDSBs1i1_mini modelDSBs1i1_nok12;do 
#for MYTARGET in PhyB2.2_R2_Feb2022 Psy1_R2_Feb2022 Psy1 PhyB2 CRTISO.cleaned Psy1 PhyB2.2.nov2021 PhyB2.2_R2_Feb2022 Psy1_R2_Feb2022 CRTISO.nov2021.49and50 CRTISO.nov2021 CRTISO.49and50bp_allb CRTISO_allb Psy1_allb PhyB2.2_allb;do 
for mymodel in modelDSBs1i1_realimprecise;do #modelDSBs1i1_3x4nor11;do #modelDSBs1i1_realimprecise;do #modelDSBs1i1_3x4 modelDSBs1i1_nok12 modelDSBs1i1_mini;do
for MYTARGET in Psy1 CRTISO.nov2021.49and50 CRTISO.nov2021 PhyB2.2.nov2021;do #PhyB2 CRTISO.cleaned CRTISO_allb CRTISO.49and50bp_allb Psy1_allb PhyB2.2_allb Psy1 PhyB2.2.nov2021 PhyB2.2_R2_Feb2022 Psy1_R2_Feb2022 CRTISO.nov2021.49and50 CRTISO.nov2021;do 
for bootstraptype in stratifiedbootstraps;do #stationarybootstraps MEbootstraps; do #MEbootstraps stationarybootstraps
for ip in `seq 1 100`;do
MYTARGETM=$( echo $MYTARGET | sed 's/_allb//' |  sed 's/_all//' | sed 's/0h//' | sed 's/m$//' | sed 's/72h//'  | sed 's/_R2_Feb2022//' )
ERRORMATRIX=~/workspace/daniela/error_matrices/error_matrix4_${MYTARGETM}_errorsfromunbroken.tsv
TIMECOURSEPATH=~/workspace/daniela/input_datasets/${bootstraptype}/timecourse_${MYINDUCTION}_${MYTARGET}${MYDELAY}
TIMECOURSE=${TIMECOURSEPATH}/timecourse.bs${ip}.txt
OUTPUTA=~/workspace/daniela/resultsv5/${bootstraptype}/timecourse_${MYINDUCTION}_${MYTARGET}${MYDELAY}/${mymodel}
OUTPUT=~/workspace/daniela/resultsv5/${bootstraptype}/timecourse_${MYINDUCTION}_${MYTARGET}${MYDELAY}/${mymodel}/results_${mymodel}_${MYINDUCTION}_ind.c${MYINDUCTION_CURVE}_${MYTARGET}${MYDELAY}
mkdir -p $OUTPUT;
OUTPUTTOT=${OUTPUT}/results.tot
OUTPUTA=${OUTPUT}/results_${ip}.txt
likfun=0; #if [ $mymodel == "modelDSBs1i1_3x4" ];then likfun=3; fi
#echo $OUTPUTA; NAMESIMS=boot2 #MEboot #${MYTARGET};
MEMS=18000; NITER=1250;
        #ls ~/workspace/daniela/resultsv4/results_${mymodel}_${MYINDUCTION}_${MYTARGET}${MYDELAY}.all_n*
        myfile="notexist"
        if [ -f ${OUTPUTA}.parsed.tab_n20000.CI.RData ];then myfile=${OUTPUTA}.parsed.tab_n20000.CI.RData;fi
        if [ -f ${OUTPUTA}.parsed.tab_n50000.CI.RData ];then myfile=${OUTPUTA}.parsed.tab_n50000.CI.RData;fi
        #if [ -f ${OUTPUTA}.parsed.tab_n150000.CI.RData ];then myfile=${OUTPUTA}.parsed.tab_n150000.CI.RData;fi
        #if [ -f ${OUTPUTA}.parsed.tab_n500000.CI.RData ];then myfile=${OUTPUTA}.parsed.tab_n500000.CI.RData;fi
        ERRORMATRIXM=$( basename $ERRORMATRIX .tsv )
        TIMECOURSEM=$( basename $TIMECOURSE .txt )
#        if [ -f ${myfile} ];then
                #cat ${myfile}.CI | awk -v var1=$MYTARGET -v var2=$MYINDUCTION -v var3=$mymodel -v var4=$ERRORMATRIXM -v var5=$TIMECOURSEM '{print $0,var1,var2,var3,var4,var5}' | grep -v max >> ~/workspace/daniela/resultsv5/table.CI
#WITH STRATIFIED
#bsub -q new-short -J plotshort -e ~/mylogs4/plot.e%J -o ~/mylogs4/plot.o%J -R rusage[mem=16000] ./plot.v1.sh $myfile $TIMECOURSE $ERRORMATRIX $mymodel 1000 ${OUTPUTA}.plot $likfun ${MYINDUCTION_CURVE} ${TIMERESOLUTION}
#WITH STRATIFIED2
bsub -q new-short -J plotshort -e ~/mylogs4/plot.e%J -o ~/mylogs4/plot.o%J -R rusage[mem=16000] ./plot.v2.sh $myfile $TIMECOURSE $ERRORMATRIX $mymodel 100 ${OUTPUTA}.plot $likfun ${MYINDUCTION_CURVE} ${TIMERESOLUTION} 1

#        fi
done
#####build flow tables
cat ${OUTPUT}/results_*.txt.plot_plot.flow.tab > ${OUTPUT}/results.tot_plot.flow.tab 
echo ${OUTPUT}/results.tot_plot.flow.tab
#####plots
TIMECOURSE=~/workspace/daniela/input_datasets/timecourse_${MYINDUCTION}_${MYTARGET}.txt
listCIFILE=/home/labs/alevy/fabrizio/workspace/daniela/resultsv5/${bootstraptype}/timecourse_${MYINDUCTION}_${MYTARGET}/${mymodel}/results_${mymodel}_${MYINDUCTION}_ind.c${MYINDUCTION_CURVE}_${MYTARGET}/listCIfiles.txt
myfiles=$( ls /home/labs/alevy/fabrizio/workspace/daniela/resultsv5/results_${mymodel}_${MYINDUCTION}_ind.c${MYINDUCTION_CURVE}_${MYTARGET}.all_n50000.CI /home/labs/alevy/fabrizio/workspace/daniela/resultsv5/${bootstraptype}/timecourse_${MYINDUCTION}_${MYTARGET}/${mymodel}/results_${mymodel}_${MYINDUCTION}_ind.c${MYINDUCTION_CURVE}_${MYTARGET}/results*.tab_n50000.CI )
for ifile in $myfiles;do
cat $ifile | awk -v OFS='\t' '{if (ar[$4]==1){$4=$4".1"}; ar[$4]=1; print $0}' > ${ifile}.t.CI 
echo $ifile
done
ls /home/labs/alevy/fabrizio/workspace/daniela/resultsv5/results_${mymodel}_${MYINDUCTION}_ind.c${MYINDUCTION_CURVE}_${MYTARGET}.all_n50000.CI.t.CI > $listCIFILE 
ls /home/labs/alevy/fabrizio/workspace/daniela/resultsv5/${bootstraptype}/timecourse_${MYINDUCTION}_${MYTARGET}/${mymodel}/results_${mymodel}_${MYINDUCTION}_ind.c${MYINDUCTION_CURVE}_${MYTARGET}/results*.tab_n50000.CI.t.CI >> $listCIFILE
bsub -q new-short -J plotshort -e ~/mylogs4/plot.e%J -o ~/mylogs4/plot.o%J -R rusage[mem=16000] ./plot.v2.sh $listCIFILE $TIMECOURSE $ERRORMATRIX $mymodel 1000 ${OUTPUTTOT}.plot $likfun ${MYINDUCTION_CURVE} ${TIMERESOLUTION} 1
done;done;done;done;done;done




