#!/bin/bash
#module load R/4.0.4

#README for the analyses of DSB repair for Daniela's project

#####################rerun everything for final numbers########################

#b) writing script that read control and generate error matrix. This is actually not necessary at the moment.
#Rather, I do:
#load("RData/error_matrices4_l.RData")
#for ( igene in names(error_matrices4_l)){
#for ( ierror in c("E_errorsfromintact","E_noerrors","E_errorsfromunbroken")){
#if (igene=="CRTISO_49and50bp"){igene2<-"CRTISO.49and50bp"} else {igene2<-igene}
#write.table(file=paste0("~/workspace/daniela/error_matrices/error_matrix4_",igene2,"_",ierror,".txt"),error_matrices4_l[[igene]][[ierror]],quote=FALSE, col.names=FALSE,row.names=FALSE,sep="\t")
#}}


### for tomato PhyB2 no controls yet, so only from time 0 
for IGENE in PhyB2.1 PhyB2.2 PhyB2.3; do cat ~/workspace/daniela/input_datasets/timecourse_RNP_${IGENE}.txt | head -3 > ~/workspace/daniela/input_datasets/control_${IGENE}_time0.txt
./calculate_error_matrix.R -i ~/workspace/daniela/input_datasets/control_${IGENE}_time0.txt -o ~/workspace/daniela/error_matrices/error_matrix4_${IGENE}_ -f 1
done

for IGENE in Psy1 CRTISO CRTISO.49and50bp; do 
./calculate_error_matrix.R -i ~/workspace/daniela/input_datasets/control_${IGENE}_Types_MH_df.csv -o ~/workspace/daniela/error_matrices/error_matrix4_${IGENE}_ -f 2
done


for MYDELAY in '' _mydelay0.25;do
for MYINDUCTION in gRNA RNP;do
for mymodel in modelDSBs1i1_realimprecise modelDSBs1i1_realimpnor11;do
for MYTARGET in Psy1 CRTISO CRTISO.49and50bp PhyB2.2 PhyB2.1m PhyB2.3m;do # CRTISO0h;do #CRTISO.49and50bp0h PhyB2.1m PhyB2.3m;do
MYTARGETM=$( echo $MYTARGET | sed 's/0h//' | sed 's/m$//' | sed 's/72h//' )
ERRORMATRIX=~/workspace/daniela/error_matrices/error_matrix4_${MYTARGETM}_errorsfromunbroken.tsv
TIMECOURSE=~/workspace/daniela/input_datasets/timecourse_${MYINDUCTION}_${MYTARGET}${MYDELAY}.txt
OUTPUTA=~/workspace/daniela/resultsv2/results_${mymodel}_${MYINDUCTION}_${MYTARGET}${MYDELAY}.a
OUTPUTB=~/workspace/daniela/resultsv2/results_${mymodel}_${MYINDUCTION}_${MYTARGET}${MYDELAY}.b
OUTPUTC=~/workspace/daniela/resultsv2/results_${mymodel}_${MYINDUCTION}_${MYTARGET}${MYDELAY}.c
OUTPUTD=~/workspace/daniela/resultsv2/results_${mymodel}_${MYINDUCTION}_${MYTARGET}${MYDELAY}.d
bsub -q new-long -J mrep -e ~/mylogs/m1.e%J -o ~/mylogs/m1.o%J -R rusage[mem=16000] ./launch_optimization_cas9_v1.sh $TIMECOURSE $ERRORMATRIX $OUTPUTA ${mymodel} 1250
bsub -q new-long -J mrep -e ~/mylogs/m1.e%J -o ~/mylogs/m1.o%J -R rusage[mem=16000] ./launch_optimization_cas9_v1.sh $TIMECOURSE $ERRORMATRIX $OUTPUTB ${mymodel} 1250
bsub -q new-long -J mrep -e ~/mylogs/m1.e%J -o ~/mylogs/m1.o%J -R rusage[mem=16000] ./launch_optimization_cas9_v1.sh $TIMECOURSE $ERRORMATRIX $OUTPUTC ${mymodel} 1250
bsub -q new-long -J mrep -e ~/mylogs/m1.e%J -o ~/mylogs/m1.o%J -R rusage[mem=16000] ./launch_optimization_cas9_v1.sh $TIMECOURSE $ERRORMATRIX $OUTPUTD ${mymodel} 1250 
done;done;done;done


for MYN in 20000;do #50000 150000 500000 2000000
for MYDELAY in '' _mydelay0.25;do
for MYINDUCTION in gRNA RNP;do
for mymodel in modelDSBs1i1_realimprecise modelDSBs1i1_realimpnor11;do
for MYTARGET in Psy1 CRTISO CRTISO.49and50bp PhyB2.2 PhyB2.1m PhyB2.3m;do # CRTISO0h;do #CRTISO.49and50bp0h PhyB2.1m PhyB2.3m;do
MYTARGETM=$( echo $MYTARGET | sed 's/0h//' | sed 's/m$//' | sed 's/72h//' )
ERRORMATRIX=~/workspace/daniela/error_matrices/error_matrix4_${MYTARGETM}_errorsfromunbroken.tsv
TIMECOURSE=~/workspace/daniela/input_datasets/timecourse_${MYINDUCTION}_${MYTARGET}${MYDELAY}.txt
OUTPUTA=~/workspace/daniela/resultsv2/results_${mymodel}_${MYINDUCTION}_${MYTARGET}${MYDELAY}.a
OUTPUTB=~/workspace/daniela/resultsv2/results_${mymodel}_${MYINDUCTION}_${MYTARGET}${MYDELAY}.b
OUTPUTC=~/workspace/daniela/resultsv2/results_${mymodel}_${MYINDUCTION}_${MYTARGET}${MYDELAY}.c
OUTPUTD=~/workspace/daniela/resultsv2/results_${mymodel}_${MYINDUCTION}_${MYTARGET}${MYDELAY}.d
OUTPUTALL=~/workspace/daniela/resultsv2/results_${mymodel}_${MYINDUCTION}_${MYTARGET}${MYDELAY}.all.tsv
if [ -f $OUTPUTA ];then
echo $MYN $TIMECOURSE
#cat $OUTPUTA | head -1 > $OUTPUTALL 
#cat $OUTPUTA $OUTPUTB $OUTPUTC $OUTPUTD | grep -v k11 | head -5000 >> $OUTPUTALL
bsub -q new-short -J CI -e ~/mylogs/confint.e%J -o ~/mylogs/confint.o%J -R rusage[mem=48000] ./calculate_CI_v1.sh $OUTPUTALL $TIMECOURSE $ERRORMATRIX $MYN $mymodel
fi
done;done;done;done;done

#for MYN in 20000 150000 500000 2000000;do #MYN=50000
echo "aaa" | awk '{printf "max\tCIlow\tCIhigh\trate\ttarget\tinduction\tmodel\terror\ttimecourse\n"}' > ~/workspace/daniela/resultsv2/table.CI
for MYDELAY in '' _mydelay0.25;do
for MYINDUCTION in gRNA RNP;do
for mymodel in modelDSBs1i1_realimprecise modelDSBs1i1_realimpnor11;do
for MYTARGET in Psy1 CRTISO CRTISO.49and50bp PhyB2.2 PhyB2.1m PhyB2.3m;do # CRTISO0h;do #CRTISO.49and50bp0h PhyB2.1m PhyB2.3m;do
MYTARGETM=$( echo $MYTARGET | sed 's/0h//' | sed 's/m$//' | sed 's/72h//' )
ERRORMATRIX=~/workspace/daniela/error_matrices/error_matrix4_${MYTARGETM}_errorsfromunbroken.tsv
TIMECOURSE=~/workspace/daniela/input_datasets/timecourse_${MYINDUCTION}_${MYTARGET}${MYDELAY}.txt
OUTPUT=~/workspace/daniela/resultsv1/results_${mymodel}_${MYINDUCTION}_${MYTARGET}${MYDELAY}
OUTPUTA=~/workspace/daniela/resultsv2/results_${mymodel}_${MYINDUCTION}_${MYTARGET}${MYDELAY}.a
if [ -f $OUTPUTA ];then
trialfile1=~/workspace/daniela/resultsv2/results_${mymodel}_${MYINDUCTION}_${MYTARGET}${MYDELAY}.all_n500000
trialfile2=~/workspace/daniela/resultsv2/results_${mymodel}_${MYINDUCTION}_${MYTARGET}${MYDELAY}.all_n150000
trialfile3=~/workspace/daniela/resultsv2/results_${mymodel}_${MYINDUCTION}_${MYTARGET}${MYDELAY}.all_n50000
trialfile4=~/workspace/daniela/resultsv2/results_${mymodel}_${MYINDUCTION}_${MYTARGET}${MYDELAY}.all_n20000
ls ~/workspace/daniela/resultsv2/results_${mymodel}_${MYINDUCTION}_${MYTARGET}${MYDELAY}.all_n*
if [ -f ${trialfile4}.CI.RData ];then myfile=$trialfile4;fi 
if [ -f ${trialfile3}.CI.RData ];then myfile=$trialfile3;fi 
if [ -f ${trialfile2}.CI.RData ];then myfile=$trialfile2;fi 
if [ -f ${trialfile1}.CI.RData ];then myfile=$trialfile1;fi 
ERRORMATRIXM=$( basename $ERRORMATRIX .tsv )
TIMECOURSEM=$( basename $TIMECOURSE .txt )
myfile=${myfile}
cat ${myfile}.CI | awk -v var1=$MYTARGET -v var2=$MYINDUCTION -v var3=$mymodel -v var4=$ERRORMATRIXM -v var5=$TIMECOURSEM '{print $0,var1,var2,var3,var4,var5}' | grep -v max >> ~/workspace/daniela/resultsv2/table.CI 
#bsub -q new-medium -J plot -e ~/mylogs/plot.e%J -o ~/mylogs/plot.o%J -R rusage[mem=24000] ./plot.v1.sh $myfile $TIMECOURSE $ERRORMATRIX $mymodel 1000 ~/workspace/daniela/resultsv2/plots/plot_${mymodel}_${MYINDUCTION}_${MYTARGET}${MYDELAY}
fi
done;done;done;done


##############################################################################

















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











