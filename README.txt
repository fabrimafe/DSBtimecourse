#!/bin/bash
#module load R/4.0.4

#README for the analyses of DSB repair for Daniela's project

#Project: writing code so that can be run on any dataset
a) improving script so that it can take an error matrix and a timecourse dataframe and it can run like this.
Now, first I parse data with the script: prepare_data.R. I then modified some of them to remove outlier point in Phyb2 (phyb2.1m),
to set non-intact to 0 in CRTISO.

b) writing script that read control and generate error matrix. This is actually not necessary at the moment.
Rather, I do:
load("RData/error_matrices4_l.RData")
for ( igene in names(error_matrices4_l)){
for ( ierror in c("E_errorsfromintact","E_noerrors","E_errorsfromunbroken")){
write.table(file=paste0("~/workspace/daniela/error_matrices/error_matrix4_",igene,"_",ierror,".txt"),error_matrices4_l[[igene]][[ierror]],quote=FALSE, col.names=FALSE,row.names=FALSE,sep="\t")
}}
ps: now it is an additional argument -E, but I could make -e optional and when not given parsing -E in some default way
c) clean general scripts and files

MYINDUCTION=RNP
for mymodel in modelDSBs1i1_realimprecise modelDSBs1i1_realimpnor11;do
for MYTARGET in Psy1 CRTISO CRTISO.49and50bp CRTISO0h CRTISO.49and50bp0h PhyB2.1 PhyB2.2 PhyB2.3 PhyB2.1m PhyB2.3m;do	 
MYDELAY=_mydelay0.25
ERRORMATRIX=~/workspace/daniela/error_matrices/error_matrix4_${MYTARGET}_E_errorsfromunbroken.txt
TIMECOURSE=~/workspace/daniela/input_datasets/timecourse_${MYINDUCTION}_${MYTARGET}${MYDELAY}.txt
OUTPUT=~/workspace/daniela/resultsv1/results_${mymodel}_${MYINDUCTION}_${MYTARGET}${MYDELAY}
bsub -q new-long -J m -e ~/mylogs/m.e%J -o ~/mylogs/m.o%J -R rusage[mem=16000] ./launch_optimization_cas9_v1.sh $TIMECOURSE $ERRORMATRIX $OUTPUT ${mymodel}
done;done

MYINDUCTION=gRNA
for mymodel in modelDSBs1i1_realimprecise modelDSBs1i1_realimpnor11;do
for MYTARGET in Psy1 CRTISO CRTISO.49and50bp CRTISO0h CRTISO.49and50bp0h;do	 
MYDELAY=_mydelay0.25
ERRORMATRIX=~/workspace/daniela/error_matrices/error_matrix4_${MYTARGET}_E_errorsfromunbroken.txt
TIMECOURSE=~/workspace/daniela/input_datasets/timecourse_${MYINDUCTION}_${MYTARGET}${MYDELAY}.txt
OUTPUT=~/workspace/daniela/resultsv1/results_${mymodel}_${MYINDUCTION}_${MYTARGET}${MYDELAY}
bsub -q new-long -J m -e ~/mylogs/m.e%J -o ~/mylogs/m.o%J -R rusage[mem=16000] ./launch_optimization_cas9_v1.sh $TIMECOURSE $ERRORMATRIX $OUTPUT ${mymodel}
done;done








#OLD RUNS FOR OLDER OPTIMIZATION (BEFORE SUMMER: delay 0) ========================================
#TO DO:
#-normalize k11*curve
#-constrain error rates to be lower than 1
#-plot everything separated (not RNP+gRNA)
#-check if I need to change something in predict_model when using different likelihood function
#
#


for target in CRTISO Psy1 CRTISO_49and50bp PhyB2.1 PhyB2.2 PhyB2.3;do
for induction in gRNA RNP;do
#for mymodel in model5i1 model5i1_nor11;do 
#for mymodel in modelDSBs1i1_fullimpreciseDSB modelDSBs1i1_realimprecise modelDSBs1i1_realimprecise modelDSBs1i1_realimpnor11 modelDSBs1i1_fullimpreciseDSB_nor11;do
for mymodel in modelDSBs1i1_realimprecise modelDSBs1i1_realimpnor11;do
#bsub -q new-long -J m1 -e m1.e -o m1.o -R rusage[mem=16000] ./launch_optimization_cas9.sh ${mymodel} ${target} $induction 1
bsub -q new-long -J m2 -e m2.e -o m2.o -R rusage[mem=16000] ./launch_optimization_cas9.sh ${mymodel} ${target} $induction 2
done
done
done

for target in CRTISO Psy1 CRTISO_49and50bp PhyB2.1 PhyB2.2 PhyB2.3;do
for induction in gRNA RNP;do
for mymodel in PhyB2.1 PhyB2.2 PhyB2.3 modelDSBs1i1_fullimpreciseDSB modelDSBs1i1_realimprecise modelDSBs1i1_realimprecise modelDSBs1i1_realimpnor11 modelDSBs1i1_fullimpreciseDSB_nor11;do
bsub -q new-long -J ${target}${induction}${mymodel} -e ${target}${induction}.e -o ${target}${induction}.o -R rusage[mem=16000] ./launch_optimization_cas9.sh ${mymodel} ${target} $induction
done
done
done


for target in PhyB2.1 PhyB2.2 PhyB2.3;do
induction=RNP
for mymodel in modelDSBs1i1_fullimpreciseDSB modelDSBs1i1_realimprecise modelDSBs1i1_realimpnor11 modelDSBs1i1_fullimpreciseDSB_nor11;do
bsub -q new-long -J ${target}${induction}${mymodel} -e ${target}${induction}.e -o ${target}${induction}.o -R rusage[mem=16000] ./launch_optimization_cas9.sh ${mymodel} ${target} $induction
done
done

myfiles=$( ls modelDSBs1i1_realimprecise_Psy1*tsv )
for i in $myfiles;do
cat ${i} | sed 's/^\t\t\t\t/target\tinduction\tmodel\terror\t/g' > results/${i}
#bsub -q new-long -J CI -e ${target}${induction}${mymodel}confint.e -o ${target}${induction}${mymodel}confint.o -R rusage[mem=24000] 
./calculate_CI_cas9.sh results/${i} 10000
done


myfiles=$( ls *tsv )
for i in $myfiles;do
NLINE=$( echo ${i} | grep model | wc -l )
if [ $NLINE -eq 0 ];then j=modelDSBs1i1_fullimpreciseDSB_${i}; echo $i;else j=${i};fi
#cat ${i} | head -1 | sed 's/^\t\t\t\t/target\tinduction\tmodel\terror\t/g' > results/${j}
#cat ${i} | grep -v k11 | grep -v "NA.*NA.*NA.*NA" >> results/${j}
bsub -q new-medium -J CI -e ${i}.confint.e -o ${i}.confint.o -R rusage[mem=24000] ./calculate_CI_cas9.sh results/${j} 50000
done
