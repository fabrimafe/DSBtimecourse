#!/bin/bash
#module load R/4.0.4

#README for the analyses of DSB repair for Daniela's project

#####################rerun everything for final numbers########################

./prepare_data.R

cat ~/workspace/daniela/input_datasets/timecourse_RNP_Psy1.txt ~/workspace/daniela/input_datasets/timecourse_RNP_Psy1_R2_Feb2022.txt | awk 'BEGIN{count=0}{if ($1=="time"){count=count+1}; if (! (count>1 && $1=="time")){print $0}}' > ~/workspace/daniela/input_datasets/timecourse_RNP_Psy1_all.txt
cat ~/workspace/daniela/input_datasets/timecourse_RNP_PhyB2.2.nov2021.txt ~/workspace/daniela/input_datasets/timecourse_RNP_PhyB2.2_R2_Feb2022.txt | awk 'BEGIN{count=0}{if ($1=="time"){count=count+1}; if (! (count>1 && $1=="time")){print $0}}' > ~/workspace/daniela/input_datasets/timecourse_RNP_PhyB2.2_all.txt

### for tomato PhyB2 no controls yet, so only from time 0 
for IGENE in PhyB2.1 PhyB2.2 PhyB2.3; do 
cat ~/workspace/daniela/input_datasets/timecourse_RNP_${IGENE}.txt | awk -v OFS='\t' '{if ($1=="time"){printf "time\ty1\ty2\ty3\n"} else if ($1==0){print $1,$2,$4+$5,$3}}' > ~/workspace/daniela/input_datasets/control_${IGENE}_time0_3states.txt
#./calculate_error_matrix.R -i ~/workspace/daniela/input_datasets/control_${IGENE}_time0.txt -o ~/workspace/daniela/error_matrices/error_matrix4_${IGENE}_ -f 1
./calculate_error_matrix.R -i ~/workspace/daniela/input_datasets/control_${IGENE}_time0_3states.txt -o ~/workspace/daniela/error_matrices/error_matrix3_${IGENE}_ -f 1 -d 3
done

inputfolder=/home/labs/alevy/fabrizio/workspace/daniela/input_datasets
csvfolder=/home/labs/alevy/fabrizio/workspace/daniela/csv
cp ${csvfolder}/CRTISO.cleanedandnov2021.49and50_control_Types_MH_df.csv ${inputfolder}/control_CRTISO.cleanedandnov2021.49and50_Types_MH_df.csv
cp ${csvfolder}/CRTISO.cleanedandnov2021.49and50_control_Types_MH_df.csv ${inputfolder}/control_CRTISO.nov2021.49and50_Types_MH_df.csv

#need the cleaned control CRTISO 49 from daniela
for IGENE in CRTISO.cleanedandnov2021.49and50 CRTISO.nov2021.49and50;do #  CRTISO.cleanedandnov2021 PhyB2.2.nov2021;do  #CRTISOcleaned Psy1 CRTISO CRTISO.49and50bp; do 
./calculate_error_matrix.R -i ~/workspace/daniela/input_datasets/control_${IGENE}_Types_MH_df.csv -o ~/workspace/daniela/error_matrices/error_matrix4_${IGENE}_ -f 2
./calculate_error_matrix.R -i ~/workspace/daniela/input_datasets/control_${IGENE}_Types_MH_df.csv -o ~/workspace/daniela/error_matrices/error_matrix3_${IGENE}_ -f 2 -d 3
done

myfiles=$( ls ~/workspace/daniela/input_datasets/timecourse_* )
for i in $myfiles;do 
rootname=$( basename $i .txt )
cat $i | awk -v OFS='\t' '{if (NR==1){print $1,$2,$3,$4} else {print $1,$2,$4+$5,$3}}' > ~/workspace/daniela/input_datasets/${rootname}_3states.txt
done

#create inputs for jacknife
mkdir -p ~/workspace/daniela/input_datasets/leave1out_jk
myfiles=$( ls ~/workspace/daniela/input_datasets/timecourse_*txt )
for i in $myfiles;do
xname=$(basename $i .txt )
newpath=~/workspace/daniela/input_datasets/leave1out_jk
newpath=${newpath}/${xname}
rm -r $newpath
mkdir -p $newpath
./create_leave1out_inputs.sh $i ${newpath}/${xname}
done
#create inputs for bootstraps
module load R/4.1.0
mkdir -p ~/workspace/daniela/input_datasets/stationarybootstraps
mkdir -p ~/workspace/daniela/input_datasets/MEbootstraps
mkdir -p ~/workspace/daniela/input_datasets/stratifiedbootstraps
mkdir -p ~/workspace/daniela/input_datasets/stratifiedbootstraps2
myfiles=$( ls ~/workspace/daniela/input_datasets/timecourse_*Psy1x5.txt )
myfiles=$( ls ~/workspace/daniela/input_datasets/timecourse_RNP_CRTISO.cleanedandnov2021*.txt )
myfiles=$( ls ~/workspace/daniela/input_datasets/timecourse_*allb.txt )
#myfiles=$( ls ~/workspace/daniela/input_datasets/timecourse_*all.txt )
for i in $myfiles;do
echo $i
xname=$(basename $i .txt )
newpath=~/workspace/daniela/input_datasets/stratifiedbootstraps/${xname}
#newpath=~/workspace/daniela/input_datasets/stratifiedbootstraps2/${xname}
rm -r $newpath;mkdir -p $newpath
./timeseriesbootstraps.R -i $i -o ${newpath} -n 100 -m 2
#./timeseriesbootstraps.R -i $i -o ${newpath} -n 100 -m 3
#newpath=~/workspace/daniela/input_datasets/stationarybootstraps/${xname}
#rm -r $newpath;mkdir -p $newpath
#./timeseriesbootstraps.R -i $i -o ${newpath} -n 100 -m 0
#newpath=~/workspace/daniela/input_datasets/MEbootstraps/${xname}
#rm -r $newpath;mkdir -p $newpath
#./timeseriesbootstraps.R -i $i -o ${newpath} -n 100 -m 1 -r 10
done

#3states
#for MYDELAY in '' _mydelay0.25;do
#for MYINDUCTION in gRNA RNP;do
#for mymodel in model5i1;do #modelDSBs1i1_realimprecise modelDSBs1i1_realimpnor11;do
#for MYTARGET in Psy1 CRTISO CRTISO.49and50bp PhyB2.2 PhyB2.1m PhyB2.3m;do # CRTISO0h;do #CRTISO.49and50bp0h PhyB2.1m PhyB2.3m;do
#MYTARGETM=$( echo $MYTARGET | sed 's/0h//' | sed 's/m$//' | sed 's/72h//' )
#ERRORMATRIX=~/workspace/daniela/error_matrices/error_matrix3_${MYTARGETM}_errorsfromunbroken.tsv
#TIMECOURSE=~/workspace/daniela/input_datasets/timecourse_${MYINDUCTION}_${MYTARGET}${MYDELAY}_3states.txt
#OUTPUTA=~/workspace/daniela/resultsv2/results_${mymodel}_${MYINDUCTION}_${MYTARGET}${MYDELAY}.a
#OUTPUTB=~/workspace/daniela/resultsv2/results_${mymodel}_${MYINDUCTION}_${MYTARGET}${MYDELAY}.b
#OUTPUTC=~/workspace/daniela/resultsv2/results_${mymodel}_${MYINDUCTION}_${MYTARGET}${MYDELAY}.c
#OUTPUTD=~/workspace/daniela/resultsv2/results_${mymodel}_${MYINDUCTION}_${MYTARGET}${MYDELAY}.d
#echo $OUTPUTA
#if [ ! -f $OUTPUTA ];then bsub -q new-long -J m3 -e ~/mylogs/m3.e%J -o ~/mylogs/m3.o%J -R rusage[mem=16000] ./launch_optimization_cas9_v1.sh $TIMECOURSE $ERRORMATRIX $OUTPUTA ${mymodel} 1250 ; fi
#if [ ! -f $OUTPUTB ];then bsub -q new-long -J m3 -e ~/mylogs/m3.e%J -o ~/mylogs/m3.o%J -R rusage[mem=16000] ./launch_optimization_cas9_v1.sh $TIMECOURSE $ERRORMATRIX $OUTPUTB ${mymodel} 1250 ; fi
#if [ ! -f $OUTPUTC ];then bsub -q new-long -J m3 -e ~/mylogs/m3.e%J -o ~/mylogs/m3.o%J -R rusage[mem=16000] ./launch_optimization_cas9_v1.sh $TIMECOURSE $ERRORMATRIX $OUTPUTC ${mymodel} 1250 ; fi
#if [ ! -f $OUTPUTD ];then bsub -q new-long -J m3 -e ~/mylogs/m3.e%J -o ~/mylogs/m3.o%J -R rusage[mem=16000] ./launch_optimization_cas9_v1.sh $TIMECOURSE $ERRORMATRIX $OUTPUTD ${mymodel} 1250 ; fi
#done;done;done;done

#for MYDELAY in '' _mydelay0.25;do
#for MYINDUCTION in gRNA RNP;do
#for mymodel in modelDSBs1i1_realnor21;do #modelDSBs1i1_realimprecise modelDSBs1i1_realimpnor11;do
#for MYTARGET in Psy1 CRTISO CRTISO.49and50bp PhyB2.2 PhyB2.1m PhyB2.3m;do # CRTISO0h;do #CRTISO.49and50bp0h PhyB2.1m PhyB2.3m;do
#MYTARGETM=$( echo $MYTARGET | sed 's/0h//' | sed 's/m$//' | sed 's/72h//' )
#ERRORMATRIX=~/workspace/daniela/error_matrices/error_matrix4_${MYTARGETM}_errorsfromunbroken.tsv
#TIMECOURSE=~/workspace/daniela/input_datasets/timecourse_${MYINDUCTION}_${MYTARGET}${MYDELAY}.txt
#OUTPUTA=~/workspace/daniela/resultsv2/results_${mymodel}_${MYINDUCTION}_${MYTARGET}${MYDELAY}.a
#OUTPUTB=~/workspace/daniela/resultsv2/results_${mymodel}_${MYINDUCTION}_${MYTARGET}${MYDELAY}.b
#OUTPUTC=~/workspace/daniela/resultsv2/results_${mymodel}_${MYINDUCTION}_${MYTARGET}${MYDELAY}.c
#OUTPUTD=~/workspace/daniela/resultsv2/results_${mymodel}_${MYINDUCTION}_${MYTARGET}${MYDELAY}.d
#echo $OUTPUTA
#if [ ! -f $OUTPUTA ];then bsub -q new-long -J m21 -e ~/mylogs/m1.e%J -o ~/mylogs/m1.o%J -R rusage[mem=16000] ./launch_optimization_cas9_v1.sh $TIMECOURSE $ERRORMATRIX $OUTPUTA ${mymodel} 1250 ; fi
#if [ ! -f $OUTPUTB ];then bsub -q new-long -J m21 -e ~/mylogs/m1.e%J -o ~/mylogs/m1.o%J -R rusage[mem=16000] ./launch_optimization_cas9_v1.sh $TIMECOURSE $ERRORMATRIX $OUTPUTB ${mymodel} 1250 ; fi
#if [ ! -f $OUTPUTC ];then bsub -q new-long -J m21 -e ~/mylogs/m1.e%J -o ~/mylogs/m1.o%J -R rusage[mem=16000] ./launch_optimization_cas9_v1.sh $TIMECOURSE $ERRORMATRIX $OUTPUTC ${mymodel} 1250 ; fi
#if [ ! -f $OUTPUTD ];then bsub -q new-long -J m21 -e ~/mylogs/m1.e%J -o ~/mylogs/m1.o%J -R rusage[mem=16000] ./launch_optimization_cas9_v1.sh $TIMECOURSE $ERRORMATRIX $OUTPUTD ${mymodel} 1250 ; fi
#done;done;done;done

####---OPTIMIZATION---#############################
for MYDELAY in '';do # _mydelay0.25;do
for MYINDUCTION in RNP; do # gRNA;do
for mymodel in modelDSBs1i1_mini;do #modelDSBs1i1_3x4 modelDSBs1i1_nok12; do # modelDSBs1i1_realimprecise model5i1 modelDSBs1i1_realimpnor11;do # modelDSBs1i1_3x4 modelDSBs1i1_realnor21 modelDSBs1i1_realimprecise modelDSBs1i1_realimpnor11;do
for MYTARGET in CRTISO.nov2021.49and50 Psy1 PhyB2.2.nov2021 CRTISO.nov2021; do #CRTISO.nov2021 CRTISO.cleanedandnov2021 Psy1 PhyB2.2.nov2021;do #CRTISOcleaned; do #Psy1 CRTISO CRTISO.49and50bp PhyB2.2 PhyB2.1m PhyB2.3m;do # CRTISO0h;do #CRTISO.49and50bp0h PhyB2.1m PhyB2.3m;do
MYTARGETM=$( echo $MYTARGET | sed 's/0h//' | sed 's/m$//' | sed 's/72h//' )
ERRORMATRIX=~/workspace/daniela/error_matrices/error_matrix4_${MYTARGETM}_errorsfromunbroken.tsv
TIMECOURSE=~/workspace/daniela/input_datasets/timecourse_${MYINDUCTION}_${MYTARGET}${MYDELAY}.txt
if [ $mymodel == "model5i1" ];then
ERRORMATRIX=~/workspace/daniela/error_matrices/error_matrix3_${MYTARGETM}_errorsfromunbroken.tsv
TIMECOURSE=~/workspace/daniela/input_datasets/timecourse_${MYINDUCTION}_${MYTARGET}${MYDELAY}_3states.txt
fi
OUTPUTA=~/workspace/daniela/resultsv2/results_${mymodel}_${MYINDUCTION}_${MYTARGET}${MYDELAY}.a
OUTPUTB=~/workspace/daniela/resultsv2/results_${mymodel}_${MYINDUCTION}_${MYTARGET}${MYDELAY}.b
OUTPUTC=~/workspace/daniela/resultsv2/results_${mymodel}_${MYINDUCTION}_${MYTARGET}${MYDELAY}.c
OUTPUTD=~/workspace/daniela/resultsv2/results_${mymodel}_${MYINDUCTION}_${MYTARGET}${MYDELAY}.d
likfun=0; if [ $mymodel == "modelDSBs1i1_3x4" ];then likfun=3; fi
echo $OUTPUTA; NAMESIMS=nok12 #${MYTARGET}; 
MEMS=18000; NITER=1250;
if [ ! -f $OUTPUTA ];then bsub -q new-long -J ${NAMESIMS} -e ~/mylogs/${NAMESIMS}-%J.e -o ~/mylogs/${NAMESIMS}-%J.o -R rusage[mem=${MEMS}] ./launch_optimization_cas9_v1b.sh $TIMECOURSE $ERRORMATRIX $OUTPUTA ${mymodel} ${NITER} $likfun ; fi
if [ ! -f $OUTPUTB ];then bsub -q new-long -J ${NAMESIMS} -e ~/mylogs/${NAMESIMS}-%J.e -o ~/mylogs/${NAMESIMS}.-%J.o -R rusage[mem=${MEMS}] ./launch_optimization_cas9_v1b.sh $TIMECOURSE $ERRORMATRIX $OUTPUTB ${mymodel} ${NITER} $likfun ; fi
if [ ! -f $OUTPUTC ];then bsub -q new-long -J ${NAMESIMS} -e ~/mylogs/${NAMESIMS}-%J.e -o ~/mylogs/${NAMESIMS}-%J.o -R rusage[mem=${MEMS}] ./launch_optimization_cas9_v1b.sh $TIMECOURSE $ERRORMATRIX $OUTPUTC ${mymodel} ${NITER} $likfun ; fi
if [ ! -f $OUTPUTD ];then bsub -q new-long -J ${NAMESIMS} -e ~/mylogs/${NAMESIMS}-%J.e -o ~/mylogs/${NAMESIMS}-%J.o -R rusage[mem=${MEMS}] ./launch_optimization_cas9_v1b.sh $TIMECOURSE $ERRORMATRIX $OUTPUTD ${mymodel} ${NITER} $likfun ; fi
done;done;done;done

####---OPTIMIZATION jk ---###########################
for MYDELAY in '';do # _mydelay0.25;do
for MYINDUCTION in RNP; do # gRNA;do
for mymodel in modelDSBs1i1_3x4;do #modelDSBs1i1_mini modelDSBs1i1_3x4; do #modelDSBs1i1_nok12; do # modelDSBs1i1_realimprecise model5i1 modelDSBs1i1_realimpnor11;do # modelDSBs1i1_3x4 modelDSBs1i1_realnor21 modelDSBs1i1_realimprecise modelDSBs1i1_realimpnor11;
# f [ ! -f $OUTPUTA ];then bsub -q new-long -J ${NAMESIMS} -e ~/mylogs/${NAMESIMS}-%J.e -o ~/mylogs/${NAMESIMS}-%J.o -R rusage[mem=${MEMS}] ./launch_optimization_cas9_v1b.sh $TIMECOURSE $ERRORMATRIX $OUTPUTA ${mymodel} ${NITER} $likfun ; fi
for MYTARGET in Psy1x5;do #Psy1x2  Psy1.MockEarlyTimePoints;do # CRTISO.nov2021.49and50; do # CRTISO.nov2021 Psy1 PhyB2.2.nov2021;do #CRTISOcleaned; do #Psy1 CRTISO CRTISO.49and50bp PhyB2.2 PhyB2.1m PhyB2.3m;do # CRTISO0h;do #CRTISO.49and50bp0h PhyB2.1m PhyB2.3m;do
for bootstraptype in MEbootstraps stationarybootstraps;do
for ip in `seq 1 100`;do
MYTARGETM=$( echo $MYTARGET | sed 's/0h//' | sed 's/m$//' | sed 's/72h//' | sed 's/Psy1.*/Psy1/' )
ERRORMATRIX=~/workspace/daniela/error_matrices/error_matrix4_${MYTARGETM}_errorsfromunbroken.tsv
TIMECOURSEPATH=~/workspace/daniela/input_datasets/${bootstraptype}/timecourse_${MYINDUCTION}_${MYTARGET}${MYDELAY}
TIMECOURSE=${TIMECOURSEPATH}/timecourse.bs${ip}.txt 
#if [ $mymodel == "model5i1" ];then
#ERRORMATRIX=~/workspace/daniela/error_matrices/error_matrix3_${MYTARGETM}_errorsfromunbroken.tsv
#TIMECOURSE=~/workspace/daniela/input_datasets/timecourse_${MYINDUCTION}_${MYTARGET}${MYDELAY}_3states.txt
#fi
#mkdir -p ~/workspace/daniela/resultsv2/leave1out_jk/timecourse_${MYINDUCTION}_${MYTARGET}${MYDELAY}
#OUTPUTA=~/workspace/daniela/resultsv2/leave1out_jk/timecourse_${MYINDUCTION}_${MYTARGET}${MYDELAY}/results_jk${ijk}.a
#likfun=0; if [ $mymodel == "modelDSBs1i1_3x4" ];then likfun=3; fi
#echo $OUTPUTA; NAMESIMS=noJK12 #${MYTARGET}; 
#MEMS=18000; NITER=1250;
#if [ ! -f $OUTPUTA ];then bsub -q new-long -J ${NAMESIMS} -e ~/mylogs/${NAMESIMS}-%J.e -o ~/mylogs/${NAMESIMS}-%J.o -R rusage[mem=${MEMS}] ./launch_optimization_cas9_v1b.sh $TIMECOURSE $ERRORMATRIX $OUTPUTA ${mymodel} ${NITER} $likfun ; fi
OUTPUTA=~/workspace/daniela/resultsv3/${bootstraptype}/timecourse_${MYINDUCTION}_${MYTARGET}${MYDELAY}/${mymodel}
mkdir -p $OUTPUTA; OUTPUTA=${OUTPUTA}/results_${ip}.txt
likfun=0; if [ $mymodel == "modelDSBs1i1_3x4" ];then likfun=3; fi
echo $OUTPUTA; NAMESIMS=MEboot #${MYTARGET}; 
MEMS=18000; NITER=1250;
if [ -f $OUTPUTA ];then rm $OUTPUTA;fi
bsub -q new-long -J ${NAMESIMS} -e ~/mylogs/${NAMESIMS}-%J.e -o ~/mylogs/${NAMESIMS}-%J.o -R rusage[mem=${MEMS}] ./launch_optimization_cas9_v1b.sh $TIMECOURSE $ERRORMATRIX $OUTPUTA ${mymodel} ${NITER} $likfun
done;done;done;done;done;done


for MYDELAY in '';do # _mydelay0.25;do
for MYINDUCTION in RNP; do # gRNA;do
for mymodel in modelDSBs1i1_mini modelDSBs1i1_3x4; do 
for MYTARGET in Psy1 CRTISO.nov2021.49and50; do # CRTISO.nov2021 Psy1 PhyB2.2.nov2021;do #CRTISOcleaned; do #Psy1 CRTISO CRTISO.49and50bp PhyB2.2 PhyB2.1m PhyB2.3m;do # CRTISO0h;do #CRTISO.49and50bp0h PhyB2.1m PhyB2.3m;do
for bootstraptype in stationarybootstraps MEbootstraps; do #MEbootstraps stationarybootstraps
for ip in `seq 1 100`;do
MYTARGETM=$( echo $MYTARGET | sed 's/0h//' | sed 's/m$//' | sed 's/72h//' )
ERRORMATRIX=~/workspace/daniela/error_matrices/error_matrix4_${MYTARGETM}_errorsfromunbroken.tsv
TIMECOURSEPATH=~/workspace/daniela/input_datasets/${bootstraptype}/timecourse_${MYINDUCTION}_${MYTARGET}${MYDELAY}
TIMECOURSE=${TIMECOURSEPATH}/timecourse.bs${ip}.txt 
OUTPUTA=~/workspace/daniela/resultsv3/${bootstraptype}/timecourse_${MYINDUCTION}_${MYTARGET}${MYDELAY}/${mymodel}
mkdir -p $OUTPUTA; 
OUTPUTTOT=${OUTPUTA}/results.tot
OUTPUTA=${OUTPUTA}/results_${ip}.txt
likfun=0; if [ $mymodel == "modelDSBs1i1_3x4" ];then likfun=3; fi
#echo $OUTPUTA; NAMESIMS=MEboot #${MYTARGET}; 
#MEMS=18000; NITER=1250;
#if [ -f $OUTPUTA ];then rm $OUTPUTA;fi
#bsub -q new-long -J ${NAMESIMS} -e ~/mylogs/${NAMESIMS}-%J.e -o ~/mylogs/${NAMESIMS}-%J.o -R rusage[mem=${MEMS}] ./launch_optimization_cas9_v1b.sh $TIMECOURSE $ERRORMATRIX $OUTPUTA ${mymodel} ${NITER} $likfun
if [ $ip -eq 1 ]; then rm ${OUTPUTTOT};fi
cat ${OUTPUTA} | awk -v var1=$ip -v OFS='\t' '{if (NR==1){print $0,"bootstrap"} else {print $0,var1}}' >> ${OUTPUTTOT}
done;done;done;done;done;done


####---CALCULATE CI LIKELIHODD---#############################
for MYN in 20000 50000 500000;do # 150000 2000000
for MYDELAY in '';do #'' _mydelay0.25;do
for MYINDUCTION in RNP;do #gRNA
for mymodel in modelDSBs1i1_mini modelDSBs1i1_3x4 modelDSBs1i1_nok12; do #modelDSBs1i1_realimprecise model5i1 modelDSBs1i1_realimpnor11;do
for MYTARGET in CRTISO.nov2021.49and50 Psy1 PhyB2.2.nov2021 CRTISO.nov2021 CRTISO.cleanedandnov2021;do #CRTISOcleaned; do #Psy1 CRTISO CRTISO.49and50bp PhyB2.2 PhyB2.1m PhyB2.3m;do # CRTISO0h;do #CRTISO.49and50bp0h PhyB2.1m PhyB2.3m;do
#for mymodel in modelDSBs1i1_3x4 model5i1 modelDSBs1i1_realnor21 modelDSBs1i1_realimprecise modelDSBs1i1_realimpnor11;do
#for MYTARGET in Psy1 CRTISO CRTISO.49and50bp PhyB2.2 PhyB2.1m PhyB2.3m;do # CRTISO0h;do #CRTISO.49and50bp0h PhyB2.1m PhyB2.3m;do
MYTARGETM=$( echo $MYTARGET | sed 's/0h//' | sed 's/m$//' | sed 's/72h//' )
ERRORMATRIX=~/workspace/daniela/error_matrices/error_matrix4_${MYTARGETM}_errorsfromunbroken.tsv
TIMECOURSE=~/workspace/daniela/input_datasets/timecourse_${MYINDUCTION}_${MYTARGET}${MYDELAY}.txt
if [ $mymodel == "model5i1" ];then
ERRORMATRIX=~/workspace/daniela/error_matrices/error_matrix3_${MYTARGETM}_errorsfromunbroken.tsv
TIMECOURSE=~/workspace/daniela/input_datasets/timecourse_${MYINDUCTION}_${MYTARGET}${MYDELAY}_3states.txt
fi
OUTPUTA=~/workspace/daniela/resultsv2/results_${mymodel}_${MYINDUCTION}_${MYTARGET}${MYDELAY}.a
OUTPUTB=~/workspace/daniela/resultsv2/results_${mymodel}_${MYINDUCTION}_${MYTARGET}${MYDELAY}.b
OUTPUTC=~/workspace/daniela/resultsv2/results_${mymodel}_${MYINDUCTION}_${MYTARGET}${MYDELAY}.c
OUTPUTD=~/workspace/daniela/resultsv2/results_${mymodel}_${MYINDUCTION}_${MYTARGET}${MYDELAY}.d
OUTPUTALL=~/workspace/daniela/resultsv2/results_${mymodel}_${MYINDUCTION}_${MYTARGET}${MYDELAY}.all.tsv
likfun=0; if [ $mymodel == "modelDSBs1i1_3x4" ];then likfun=3; fi
if [ -f $OUTPUTA ];then
echo $MYN $TIMECOURSE
cat $OUTPUTA | head -1 > $OUTPUTALL 
cat $OUTPUTA $OUTPUTB $OUTPUTC $OUTPUTD | grep -v k11 | head -5000 >> $OUTPUTALL
bsub -q new-short -J CI -e ~/mylogs2/confint.e%J -o ~/mylogs2/confint.o%J -R rusage[mem=10000] ./calculate_CI_v1.sh $OUTPUTALL $TIMECOURSE $ERRORMATRIX $MYN $mymodel $likfun
fi
done;done;done;done;done

echo "aaa" | awk '{printf "max\tCIlow\tCIhigh\trate\ttarget\tinduction\tmodel\terror\ttimecourse\n"}' > ~/workspace/daniela/resultsv2/table.CI
for MYDELAY in '';do # _mydelay0.25;do
for MYINDUCTION in RNP;do #gRNA
for mymodel in modelDSBs1i1_mini modelDSBs1i1_3x4 modelDSBs1i1_nok12 modelDSBs1i1_realimprecise;do # modelDSBs1i1_realimprecise model5i1 modelDSBs1i1_realimpnor11;do
for MYTARGET in Psy1 PhyB2.2.nov2021 CRTISO.nov2021 CRTISO.cleanedandnov2021 CRTISO.nov2021.49and50 ;do # Psy1 PhyB2.2.nov2021 CRTISOcleaned; do #Psy1 CRTISO CRTISO.49and50bp PhyB2.2 PhyB2.1m PhyB2.3m;do # CRTISO0h;do #CRTISO.49and50bp0h PhyB2.1m PhyB2.3m;do
MYTARGETM=$( echo $MYTARGET | sed 's/0h//' | sed 's/m$//' | sed 's/72h//' )
ERRORMATRIX=~/workspace/daniela/error_matrices/error_matrix4_${MYTARGETM}_errorsfromunbroken.tsv
TIMECOURSE=~/workspace/daniela/input_datasets/timecourse_${MYINDUCTION}_${MYTARGET}${MYDELAY}.txt
likfun=0; if [ $mymodel == "modelDSBs1i1_3x4" ];then likfun=3; fi
if [ $mymodel == "model5i1" ];then
ERRORMATRIX=~/workspace/daniela/error_matrices/error_matrix3_${MYTARGETM}_errorsfromunbroken.tsv
TIMECOURSE=~/workspace/daniela/input_datasets/timecourse_${MYINDUCTION}_${MYTARGET}${MYDELAY}_3states.txt
fi
OUTPUTA=~/workspace/daniela/resultsv2/results_${mymodel}_${MYINDUCTION}_${MYTARGET}${MYDELAY}.a
if [ -f $OUTPUTA ];then
	trialfile1=~/workspace/daniela/resultsv2/results_${mymodel}_${MYINDUCTION}_${MYTARGET}${MYDELAY}.all_n500000
	trialfile2=~/workspace/daniela/resultsv2/results_${mymodel}_${MYINDUCTION}_${MYTARGET}${MYDELAY}.all_n150000
	trialfile3=~/workspace/daniela/resultsv2/results_${mymodel}_${MYINDUCTION}_${MYTARGET}${MYDELAY}.all_n50000
	trialfile4=~/workspace/daniela/resultsv2/results_${mymodel}_${MYINDUCTION}_${MYTARGET}${MYDELAY}.all_n20000
	ls ~/workspace/daniela/resultsv2/results_${mymodel}_${MYINDUCTION}_${MYTARGET}${MYDELAY}.all_n*
	myfile="notexist"
	if [ -f ${trialfile4}.CI.RData ];then myfile=$trialfile4;fi 
	if [ -f ${trialfile3}.CI.RData ];then myfile=$trialfile3;fi 
	if [ -f ${trialfile2}.CI.RData ];then myfile=$trialfile2;fi 
	if [ -f ${trialfile1}.CI.RData ];then myfile=$trialfile1;fi 
	ERRORMATRIXM=$( basename $ERRORMATRIX .tsv )
	TIMECOURSEM=$( basename $TIMECOURSE .txt )
	myfile=${myfile}
	if [ -f ${myfile}.CI.RData ];then
		#cat ${myfile}.CI | awk -v var1=$MYTARGET -v var2=$MYINDUCTION -v var3=$mymodel -v var4=$ERRORMATRIXM -v var5=$TIMECOURSEM '{print $0,var1,var2,var3,var4,var5}' | grep -v max >> ~/workspace/daniela/resultsv2/table.CI 
		bsub -q new-short -J plot -e ~/mylogs/plot.e%J -o ~/mylogs/plot.o%J -R rusage[mem=8000] ./plot.v1.sh $myfile $TIMECOURSE $ERRORMATRIX $mymodel 1000 ~/workspace/daniela/resultsv2/plots/plot_${mymodel}_${MYINDUCTION}_${MYTARGET}${MYDELAY} $likfun
	fi
fi
done;done;done;done

MYDELAY='' # _mydelay0.25;do
MYINDUCTION=RNP #gRNA
mymodel=modelDSBs1i1_3x4; # modelDSBs1i1_realimprecise model5i1 modelDSBs1i1_realimpnor11;do
MYTARGET=Psy1; #PhyB2.2.nov2021 CRTISO.nov2021 CRTISO.cleanedandnov2021;do #CRTISOcleaned; do #Psy1 CRTISO CRTISO.49and50bp PhyB2.2 PhyB2.1m PhyB2.3m;do # CRTISO0h;do #CRTISO.49and50bp0h PhyB2.1m PhyB2.3m;do
MYTARGETM=$( echo $MYTARGET | sed 's/0h//' | sed 's/m$//' | sed 's/72h//' )
ERRORMATRIX=~/workspace/daniela/error_matrices/error_matrix4_${MYTARGETM}_errorsfromunbroken.tsv
TIMECOURSE=~/workspace/daniela/input_datasets/timecourse_${MYINDUCTION}_${MYTARGET}${MYDELAY}.txt
if [ $mymodel == "model5i1" ];then
ERRORMATRIX=~/workspace/daniela/error_matrices/error_matrix3_${MYTARGETM}_errorsfromunbroken.tsv
TIMECOURSE=~/workspace/daniela/input_datasets/timecourse_${MYINDUCTION}_${MYTARGET}${MYDELAY}_3states.txt
fi
OUTPUTA=~/workspace/daniela/resultsv2/results_${mymodel}_${MYINDUCTION}_${MYTARGET}${MYDELAY}.a
	trialfile3=~/workspace/daniela/resultsv2/results_${mymodel}_${MYINDUCTION}_${MYTARGET}${MYDELAY}.all_n50000
	if [ -f ${trialfile3}.CI.RData ];then myfile=$trialfile3;fi 
	ERRORMATRIXM=$( basename $ERRORMATRIX .tsv )
	TIMECOURSEM=$( basename $TIMECOURSE .txt )
	myfile=${myfile}
./plot.v1.sh $myfile $TIMECOURSE $ERRORMATRIX $mymodel 100 ~/workspace/daniela/resultsv2/plots/plot_${mymodel}_${MYINDUCTION}_${MYTARGET}${MYDELAY}

for MYINDUCTION in gRNA RNP;do
for mymodel in modelDSBs1i1_realimprecise modelDSBs1i1_realimpnor11;do
for MYTARGET in Psy1 CRTISO CRTISO.49and50bp Psy172h CRTISO0h CRTISO.49and50bp0h PhyB2.1 PhyB2.2 PhyB2.3 PhyB2.1m PhyB2.3m;do
MYTARGETM=$( echo $MYTARGET | sed 's/0h//' | sed 's/m$//' | sed 's/72h$//' )
for MYDELAY in '' '_mydelay0.25';do # '';do
ERRORMATRIX=~/workspace/daniela/error_matrices/error_matrix4_${MYTARGETM}_E_errorsfromunbroken.txt
TIMECOURSE=~/workspace/daniela/input_datasets/timecourse_${MYINDUCTION}_${MYTARGET}${MYDELAY}.txt
OUTPUT=~/workspace/daniela/resultsv1/results_${mymodel}_${MYINDUCTION}_${MYTARGET}${MYDELAY}
#bsub -q new-long -J m -e ~/mylogs/m.e%J -o ~/mylogs/m.o%J -R rusage[mem=16000] ./launch_optimization_cas9_v1.sh $TIMECOURSE $ERRORMATRIX $OUTPUT ${mymodel}
INPUT=${OUTPUT}.parsed.tsv
bsub -q new-medium -J plot -e ~/mylogs/plot.e%J -o ~/mylogs/plot.o%J -R rusage[mem=24000] ./plot.v1.sh $INPUT $TIMECOURSE $ERRORMATRIX $mymodel 1000
done;done;done;done

###################################BY INDUCTION###########################################

cat ~/workspace/daniela/input_datasets/timecourse_RNP_Psy1.txt ~/workspace/daniela/input_datasets/timecourse_RNP_CRTISO.txt ~/workspace/daniela/input_datasets/timecourse_RNP_PhyB2.2.txt | head -1 > ~/workspace/daniela/input_datasets/timecourse_RNP_Psy1CRTISOPhyB2.2.txt
cat ~/workspace/daniela/input_datasets/timecourse_RNP_Psy1.txt ~/workspace/daniela/input_datasets/timecourse_RNP_CRTISO.txt ~/workspace/daniela/input_datasets/timecourse_RNP_PhyB2.2.txt | grep -v time >> ~/workspace/daniela/input_datasets/timecourse_RNP_Psy1CRTISOPhyB2.2.txt
cat ~/workspace/daniela/error_matrices/error_matrix4_CRTISO_errorsfromunbroken.tsv ~/workspace/daniela/error_matrices/error_matrix4_Psy1_errorsfromunbroken.tsv ~/workspace/daniela/error_matrices/error_matrix4_PhyB2.2_errorsfromunbroken.tsv > ~/workspace/daniela/error_matrices/error_matrix4_Psy1CRTISOPhyB2.2_errorsfromunbroken.tsv

#CI to be rerun because not a long of iterations with induction, since it is slow. opt for RNP still running as mrep
MYN=50000
MYDELAY=''
mymodel=modelDSBs1i1_realimprecise.inductionx3
MYINDUCTION=RNP #gRNA #RNP #gRNA #RNP
MYTARGET=Psy1CRTISOPhyB2.2
TIMECOURSE=~/workspace/daniela/input_datasets/timecourse_RNP_Psy1CRTISOPhyB2.2.txt
ERRORMATRIX=~/workspace/daniela/error_matrices/error_matrix4_Psy1CRTISOPhyB2.2_errorsfromunbroken.tsv
OUTPUTA=~/workspace/daniela/resultsv2/results_${mymodel}_${MYINDUCTION}_${MYTARGET}${MYDELAY}.a
OUTPUTB=~/workspace/daniela/resultsv2/results_${mymodel}_${MYINDUCTION}_${MYTARGET}${MYDELAY}.b
OUTPUTC=~/workspace/daniela/resultsv2/results_${mymodel}_${MYINDUCTION}_${MYTARGET}${MYDELAY}.c
OUTPUTD=~/workspace/daniela/resultsv2/results_${mymodel}_${MYINDUCTION}_${MYTARGET}${MYDELAY}.d
OUTPUTALL=~/workspace/daniela/resultsv2/results_${mymodel}_${MYINDUCTION}_${MYTARGET}${MYDELAY}.all.tsv
#bsub -q new-long -J mind -e ~/mylogs/m1.e%J -o ~/mylogs/m1.o%J -R rusage[mem=16000] ./launch_optimization_cas9_v1.sh $TIMECOURSE $ERRORMATRIX $OUTPUTA ${mymodel} 1250
#bsub -q new-long -J mind -e ~/mylogs/m1.e%J -o ~/mylogs/m1.o%J -R rusage[mem=16000] ./launch_optimization_cas9_v1.sh $TIMECOURSE $ERRORMATRIX $OUTPUTB ${mymodel} 1250
#bsub -q new-long -J mind -e ~/mylogs/m1.e%J -o ~/mylogs/m1.o%J -R rusage[mem=16000] ./launch_optimization_cas9_v1.sh $TIMECOURSE $ERRORMATRIX $OUTPUTC ${mymodel} 1250
#bsub -q new-long -J mind -e ~/mylogs/m1.e%J -o ~/mylogs/m1.o%J -R rusage[mem=16000] ./launch_optimization_cas9_v1.sh $TIMECOURSE $ERRORMATRIX $OUTPUTD ${mymodel} 1250
cat $OUTPUTA | head -1 > $OUTPUTALL 
cat $OUTPUTA $OUTPUTB $OUTPUTC $OUTPUTD | grep -v k11 | head -5000 >> $OUTPUTALL
bsub -q new-long -J CI.ind -e ~/mylogs/CI.ind.e%J -o ~/mylogs/CI.ind.o%J -R rusage[mem=48000] ./calculate_CI_v1.sh $OUTPUTALL $TIMECOURSE $ERRORMATRIX $MYN $mymodel 0



#check systematically what CI files are missing
myfiles=$( ls *.[a])
for i in $myfiles;do xfile=$( basename $i .a);ls ${xfile}*CI;done 2> errors.log








##########################################
######### RUN OPTIMIZATION ###############
##########################################

MYINDUCTION=RNP
for mymodel in modelDSBs1i1_realimprecise modelDSBs1i1_realimpnor11;do
for MYTARGET in Psy1 CRTISO CRTISO.49and50bp CRTISO0h CRTISO.49and50bp0h PhyB2.1 PhyB2.2 PhyB2.3 PhyB2.1m PhyB2.3m;do	 
MYTARGETM=$( echo $MYTARGET | sed 's/0h//' | sed 's/m$//' )
MYDELAY=_mydelay0.25
ERRORMATRIX=~/workspace/daniela/error_matrices/error_matrix4_${MYTARGETM}_E_errorsfromunbroken.txt
TIMECOURSE=~/workspace/daniela/input_datasets/timecourse_${MYINDUCTION}_${MYTARGET}${MYDELAY}.txt
OUTPUT=~/workspace/daniela/resultsv1/results_${mymodel}_${MYINDUCTION}_${MYTARGET}${MYDELAY}
bsub -q new-long -J m -e ~/mylogs/m.e%J -o ~/mylogs/m.o%J -R rusage[mem=16000] ./launch_optimization_cas9_v1.sh $TIMECOURSE $ERRORMATRIX $OUTPUT ${mymodel}
done;done

MYINDUCTION=gRNA
for mymodel in modelDSBs1i1_realimprecise modelDSBs1i1_realimpnor11;do
for MYTARGET in Psy1 CRTISO CRTISO.49and50bp CRTISO0h CRTISO.49and50bp0h;do	 
MYTARGETM=$( echo $MYTARGET | sed 's/0h//' | sed 's/m$//' )
MYDELAY=_mydelay0.25
ERRORMATRIX=~/workspace/daniela/error_matrices/error_matrix4_${MYTARGETM}_E_errorsfromunbroken.txt
TIMECOURSE=~/workspace/daniela/input_datasets/timecourse_${MYINDUCTION}_${MYTARGET}${MYDELAY}.txt
OUTPUT=~/workspace/daniela/resultsv1/results_${mymodel}_${MYINDUCTION}_${MYTARGET}${MYDELAY}
bsub -q new-long -J m -e ~/mylogs/m.e%J -o ~/mylogs/m.o%J -R rusage[mem=16000] ./launch_optimization_cas9_v1.sh $TIMECOURSE $ERRORMATRIX $OUTPUT ${mymodel}
done;done

MYDELAY=''
MYINDUCTION=RNP
for mymodel in modelDSBs1i1_realimprecise modelDSBs1i1_realimpnor11;do
for MYTARGET in Psy1 CRTISO CRTISO.49and50bp;do #CRTISO0h CRTISO.49and50bp0h PhyB2.1m PhyB2.3m;do	 
MYTARGETM=$( echo $MYTARGET | sed 's/0h//' | sed 's/m$//' )
ERRORMATRIX=~/workspace/daniela/error_matrices/error_matrix4_${MYTARGETM}_E_errorsfromunbroken.txt
TIMECOURSE=~/workspace/daniela/input_datasets/timecourse_${MYINDUCTION}_${MYTARGET}${MYDELAY}.txt
OUTPUT=~/workspace/daniela/resultsv1/results_${mymodel}_${MYINDUCTION}_${MYTARGET}${MYDELAY}
bsub -q new-long -J m -e ~/mylogs/m1.e%J -o ~/mylogs/m1.o%J -R rusage[mem=16000] ./launch_optimization_cas9_v1.sh $TIMECOURSE $ERRORMATRIX $OUTPUT ${mymodel}
done;done
MYINDUCTION=gRNA
for mymodel in modelDSBs1i1_realimprecise modelDSBs1i1_realimpnor11;do
for MYTARGET in Psy1 CRTISO CRTISO.49and50bp;do #CRTISO0h CRTISO.49and50bp0h;do	 
MYTARGETM=$( echo $MYTARGET | sed 's/0h//' | sed 's/m$//' )
ERRORMATRIX=~/workspace/daniela/error_matrices/error_matrix4_${MYTARGETM}_E_errorsfromunbroken.txt
TIMECOURSE=~/workspace/daniela/input_datasets/timecourse_${MYINDUCTION}_${MYTARGET}${MYDELAY}.txt
OUTPUT=~/workspace/daniela/resultsv1/results_${mymodel}_${MYINDUCTION}_${MYTARGET}${MYDELAY}
bsub -q new-long -J m -e ~/mylogs/m1.e%J -o ~/mylogs/m1.o%J -R rusage[mem=16000] ./launch_optimization_cas9_v1.sh $TIMECOURSE $ERRORMATRIX $OUTPUT ${mymodel}
done;done


MYDELAY=''
for MYINDUCTION in gRNA RNP;do
for mymodel in modelDSBs1i1_realimprecise modelDSBs1i1_realimpnor11;do
for MYTARGET in Psy172h;do #Psy1 CRTISO CRTISO.49and50bp PhyB2.2;do #CRTISO0h CRTISO.49and50bp0h PhyB2.1m PhyB2.3m;do
MYTARGETM=$( echo $MYTARGET | sed 's/0h//' | sed 's/m$//' | sed 's/72h//' )
ERRORMATRIX=~/workspace/daniela/error_matrices/error_matrix4_${MYTARGETM}_E_errorsfromunbroken.txt
TIMECOURSE=~/workspace/daniela/input_datasets/timecourse_${MYINDUCTION}_${MYTARGET}${MYDELAY}.txt
OUTPUTA=~/workspace/daniela/resultsv1/results_${mymodel}_${MYINDUCTION}_${MYTARGET}${MYDELAY}.a
OUTPUTB=~/workspace/daniela/resultsv1/results_${mymodel}_${MYINDUCTION}_${MYTARGET}${MYDELAY}.b
OUTPUTC=~/workspace/daniela/resultsv1/results_${mymodel}_${MYINDUCTION}_${MYTARGET}${MYDELAY}.c
bsub -q new-long -J mrep -e ~/mylogs/m1.e%J -o ~/mylogs/m1.o%J -R rusage[mem=16000] ./launch_optimization_cas9_v1.sh $TIMECOURSE $ERRORMATRIX $OUTPUTA ${mymodel}
bsub -q new-long -J mrep -e ~/mylogs/m1.e%J -o ~/mylogs/m1.o%J -R rusage[mem=16000] ./launch_optimization_cas9_v1.sh $TIMECOURSE $ERRORMATRIX $OUTPUTB ${mymodel}
bsub -q new-long -J mrep -e ~/mylogs/m1.e%J -o ~/mylogs/m1.o%J -R rusage[mem=16000] ./launch_optimization_cas9_v1.sh $TIMECOURSE $ERRORMATRIX $OUTPUTC ${mymodel}
done;done;done


#######################################
############## CALCULATE CI ###########
#######################################
for MYINDUCTION in gRNA RNP;do
for mymodel in modelDSBs1i1_realimprecise modelDSBs1i1_realimpnor11;do
for MYTARGET in Psy1 CRTISO CRTISO.49and50bp CRTISO0h CRTISO.49and50bp0h PhyB2.1 PhyB2.2 PhyB2.3 PhyB2.1m PhyB2.3m;do	 
MYTARGETM=$( echo $MYTARGET | sed 's/0h//' | sed 's/m$//' )
for MYDELAY in '' '_mydelay0.25' '';do
ERRORMATRIX=~/workspace/daniela/error_matrices/error_matrix4_${MYTARGETM}_E_errorsfromunbroken.txt
TIMECOURSE=~/workspace/daniela/input_datasets/timecourse_${MYINDUCTION}_${MYTARGET}${MYDELAY}.txt
OUTPUT=~/workspace/daniela/resultsv1/results_${mymodel}_${MYINDUCTION}_${MYTARGET}${MYDELAY}
#bsub -q new-long -J m -e ~/mylogs/m.e%J -o ~/mylogs/m.o%J -R rusage[mem=16000] ./launch_optimization_cas9_v1.sh $TIMECOURSE $ERRORMATRIX $OUTPUT ${mymodel}
INPUT=${OUTPUT}.parsed.tsv
if [ -f ${INPUT} ];then
	echo $INPUT
	cat ${OUTPUT} | sed 's/^\t\t/model\t/' | sed 's/^NA\t//' | awk 'BEGIN{count=0}{if ($1=="model"){count=count+1;if (count==1){print}} else print}' > ${INPUT}
	bsub -q new-medium -J CI.2 -e ~/mylogs/confint.e%J -o ~/mylogs/confint.o%J -R rusage[mem=24000] ./calculate_CI_v1.sh $INPUT $TIMECOURSE $ERRORMATRIX 50000
fi
done;done;done;done

#some errors with modified PhyB, can't understand why. So less iterations.
for MYINDUCTION in gRNA RNP;do
for mymodel in modelDSBs1i1_realimprecise modelDSBs1i1_realimpnor11;do
for MYTARGET in PhyB2.1m PhyB2.3m;do
MYTARGETM=$( echo $MYTARGET | sed 's/0h//' | sed 's/m$//' )
for MYDELAY in '_mydelay0.25' '';do
ERRORMATRIX=~/workspace/daniela/error_matrices/error_matrix4_${MYTARGETM}_E_errorsfromunbroken.txt
TIMECOURSE=~/workspace/daniela/input_datasets/timecourse_${MYINDUCTION}_${MYTARGET}${MYDELAY}.txt
OUTPUT=~/workspace/daniela/resultsv1/results_${mymodel}_${MYINDUCTION}_${MYTARGET}${MYDELAY}
#bsub -q new-long -J m -e ~/mylogs/m.e%J -o ~/mylogs/m.o%J -R rusage[mem=16000] ./launch_optimization_cas9_v1.sh $TIMECOURSE $ERRORMATRIX $OUTPUT ${mymodel}
INPUT=${OUTPUT}.parsed.tsv
if [ -f ${INPUT} ];then
        echo $INPUT
        cat ${OUTPUT} | sed 's/^\t\t/model\t/' | sed 's/^NA\t//' | awk 'BEGIN{count=0}{if ($1=="model"){count=count+1;if (count==1){print}} else print}' > ${INPUT}
        bsub -q new-medium -J CI -e ~/mylogs/confint.e%J -o ~/mylogs/confint.o%J -R rusage[mem=24000] ./calculate_CI_v1.sh $INPUT $TIMECOURSE $ERRORMATRIX 20000
fi
done;done;done;done

#for old calculations (results)
#for MYINDUCTION in gRNA RNP;do
#for mymodel in modelDSBs1i1_realimprecise modelDSBs1i1_realimpnor11;do
#for MYTARGET in Psy1 CRTISO CRTISO.49and50bp PhyB2.1 PhyB2.2 PhyB2.3;do	 
#MYTARGETM=$( echo $MYTARGET | sed 's/0h//' | sed 's/m$//' )
#ERRORMATRIX=~/workspace/daniela/error_matrices/error_matrix4_${MYTARGETM}_E_errorsfromunbroken.txt
#TIMECOURSE=~/workspace/daniela/input_datasets/timecourse_${MYINDUCTION}_${MYTARGET}${MYDELAY}.txt
#OUTPUT=~/workspace/daniela/${mymodel}_${MYTARGET}_${MYINDUCTION}_optimization.tsv
#bsub -q new-long -J m -e ~/mylogs/m.e%J -o ~/mylogs/m.o%J -R rusage[mem=16000] ./launch_optimization_cas9_v1.sh $TIMECOURSE $ERRORMATRIX $OUTPUT ${mymodel}
#INPUT=${OUTPUT}.parsed.tsv
#if [ -f ${OUTPUT} ];then
#	echo $INPUT
#cat ${OUTPUT} | sed 's/^\t\t\t\t/target\tinduction\tmodel\terror\t/g' | sed 's/^NA\t//' | awk 'BEGIN{count=0}{if ($1=="target"){count=count+1;if (count==1){print}} else print}' > ${INPUT}
#	bsub -q new-medium -J CI -e ~/mylogs/confint.e%J -o ~/mylogs/confint.o%J -R rusage[mem=24000] ./calculate_CI_v1.sh $INPUT $TIMECOURSE $ERRORMATRIX 500000
#fi
#done;done;done;

#for run in parallel

MYDELAY=''
for MYINDUCTION in gRNA RNP;do
for mymodel in modelDSBs1i1_realimprecise modelDSBs1i1_realimpnor11;do
for MYTARGET in Psy1 CRTISO CRTISO.49and50bp PhyB2.2 Psy72h;do #Psy172h CRTISO0h CRTISO.49and50bp0h PhyB2.1m PhyB2.3m;do
MYTARGETM=$( echo $MYTARGET | sed 's/0h//' | sed 's/m$//' )
ERRORMATRIX=~/workspace/daniela/error_matrices/error_matrix4_${MYTARGETM}_E_errorsfromunbroken.txt
TIMECOURSE=~/workspace/daniela/input_datasets/timecourse_${MYINDUCTION}_${MYTARGET}${MYDELAY}.txt
OUTPUT=~/workspace/daniela/resultsv1/results_${mymodel}_${MYINDUCTION}_${MYTARGET}${MYDELAY}
OUTPUTA=~/workspace/daniela/resultsv1/results_${mymodel}_${MYINDUCTION}_${MYTARGET}${MYDELAY}.a
OUTPUTB=~/workspace/daniela/resultsv1/results_${mymodel}_${MYINDUCTION}_${MYTARGET}${MYDELAY}.b
OUTPUTC=~/workspace/daniela/resultsv1/results_${mymodel}_${MYINDUCTION}_${MYTARGET}${MYDELAY}.c
INPUT=${OUTPUT}.parsed.tsv
cat $OUTPUTA $OUTPUTB $OUTPUTC | awk 'BEGIN{count=0}{if ($1=="k11"){count=count+1;if (count==1){print}} else print}' > $OUTPUT 
cat ${OUTPUT} | sed 's/^\t\t/model\t/' | sed 's/^NA\t//' | awk 'BEGIN{count=0}{if ($1=="model"){count=count+1;if (count==1){print}} else print}' > ${INPUT}
	bsub -q new-medium -J CI -e ~/mylogs/confint.e%J -o ~/mylogs/confint.o%J -R rusage[mem=24000] ./calculate_CI_v1.sh $INPUT $TIMECOURSE $ERRORMATRIX 50000
done;done;done

#########################################
############# PARSE IN A TABLE ##########
#########################################

for MYTARGET in Psy1 CRTISO CRTISO.49and50bp Psy172h CRTISO0h CRTISO.49and50bp0h PhyB2.1 PhyB2.2 PhyB2.3 PhyB2.1m PhyB2.3m;do
for MYINDUCTION in gRNA RNP;do
for mymodel in modelDSBs1i1_realimprecise modelDSBs1i1_realimpnor11;do
MYTARGETM=$( echo $MYTARGET | sed 's/0h//' | sed 's/m$//' )
for MYDELAY in '_mydelay0.25' '';do
if [ MYDELAY == '' ];then 
	XDELAY=0; 
elif [ MYDELAY == '_mydelay0.25' ];then
	XDELAY=0.25
fi	
ERRORMATRIX=~/workspace/daniela/error_matrices/error_matrix4_${MYTARGETM}_E_errorsfromunbroken.txt
TIMECOURSE=~/workspace/daniela/input_datasets/timecourse_${MYINDUCTION}_${MYTARGET}${MYDELAY}.txt
OUTPUT=~/workspace/daniela/resultsv1/results_${mymodel}_${MYINDUCTION}_${MYTARGET}${MYDELAY}
#bsub -q new-long -J m -e ~/mylogs/m.e%J -o ~/mylogs/m.o%J -R rusage[mem=16000] ./launch_optimization_cas9_v1.sh $TIMECOURSE $ERRORMATRIX $OUTPUT ${mymodel}
INPUT=${OUTPUT}.parsed.CI
if [ -f ${INPUT} ];then
cat ${INPUT} | awk -v v1=$MYTARGET -v v2=$MYINDUCTION -v v3=$mymodel -v v4=$XDELAY -v OFS="\t" '{print v1,v2,v3,v4,$0}' >> ~/workspace/daniela/resultsv1/table_CI.tsv
fi
done;done;done;done


#old calculations
XDELAY=0
for MYINDUCTION in gRNA RNP;do
for mymodel in modelDSBs1i1_realimprecise modelDSBs1i1_realimpnor11;do
for MYTARGET in Psy1 CRTISO CRTISO.49and50bp PhyB2.1 PhyB2.2 PhyB2.3;do
MYTARGETM=$( echo $MYTARGET | sed 's/0h//' | sed 's/m$//' )
ERRORMATRIX=~/workspace/daniela/error_matrices/error_matrix4_${MYTARGETM}_E_errorsfromunbroken.txt
TIMECOURSE=~/workspace/daniela/input_datasets/timecourse_${MYINDUCTION}_${MYTARGET}${MYDELAY}.txt
OUTPUT=~/workspace/daniela/${mymodel}_${MYTARGET}_${MYINDUCTION}_optimization.tsv
#bsub -q new-long -J m -e ~/mylogs/m.e%J -o ~/mylogs/m.o%J -R rusage[mem=16000] ./launch_optimization_cas9_v1.sh $TIMECOURSE $ERRORMATRIX $OUTPUT ${mymodel}
INPUT=${OUTPUT}.parsed.CI
if [ -f ${OUTPUT} ];then
        echo $INPUT
cat ${INPUT} | awk -v v1=$MYTARGET -v v2=$MYINDUCTION -v v3=$mymodel -v v4=$XDELAY -v OFS="\t" '{print v1,v2,v3,v4,$0}' >> ~/workspace/daniela/resultsv1/table_CI.tsv
fi
done;done;done;




##################
#### PLOTTING ####
##################


for MYINDUCTION in gRNA RNP;do
for mymodel in modelDSBs1i1_realimprecise modelDSBs1i1_realimpnor11;do
for MYTARGET in Psy1 CRTISO CRTISO.49and50bp Psy172h CRTISO0h CRTISO.49and50bp0h PhyB2.1 PhyB2.2 PhyB2.3 PhyB2.1m PhyB2.3m;do
MYTARGETM=$( echo $MYTARGET | sed 's/0h//' | sed 's/m$//' | sed 's/72h$//' )
for MYDELAY in '' '_mydelay0.25';do # '';do
ERRORMATRIX=~/workspace/daniela/error_matrices/error_matrix4_${MYTARGETM}_E_errorsfromunbroken.txt
TIMECOURSE=~/workspace/daniela/input_datasets/timecourse_${MYINDUCTION}_${MYTARGET}${MYDELAY}.txt
OUTPUT=~/workspace/daniela/resultsv1/results_${mymodel}_${MYINDUCTION}_${MYTARGET}${MYDELAY}
#bsub -q new-long -J m -e ~/mylogs/m.e%J -o ~/mylogs/m.o%J -R rusage[mem=16000] ./launch_optimization_cas9_v1.sh $TIMECOURSE $ERRORMATRIX $OUTPUT ${mymodel}
INPUT=${OUTPUT}.parsed.tsv
bsub -q new-medium -J plot -e ~/mylogs/plot.e%J -o ~/mylogs/plot.o%J -R rusage[mem=24000] ./plot.v1.sh $INPUT $TIMECOURSE $ERRORMATRIX $mymodel 1000
done;done;done;done




#./plot_CI.v1.R -i ~/workspace/daniela/resultsv1/results_modelDSBs1i1_realimpnor11_gRNA_CRTISO0h.parsed.CI.RData -d ~/workspace/daniela/input_datasets/timecourse_gRNA_CRTISO0h.txt -E ~/workspace/daniela/error_matrices/error_matrix4_CRTISO_E_errorsfromunbroken.txt -m modelDSBs1i1_realimpnor11 -n 10 -o temp


############################# I wish it was final.....30-11-2021##################









