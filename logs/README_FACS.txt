#FACS
#ESTIMATE INDUCTION CURVES
cd /home/labs/alevy/fabrizio/workspace/daniela/FACStable
library(tidyverse)
library(optimx)

for ( xFACS in c("FACStable","FACStable.5","FACStable.90sec")){ 
xFACSinputfile<-paste0(xFACS,".txt")
xdata<-read.table(xFACSinputfile,header=TRUE)
induction_curve_3params<-function(x,K,r0,r2,yno1=0) { K1<-(1-K)*2^(-x*r2)/(1+exp(-r0*(x-log(10^6)/(r0+10^(-12)))))-yno1; if (K1<0){K1<-0};K1 }
induction_curve_vectorized_3params<-function(x,params) {induction_curve(x,params[1],params[2],params[3],params[4])}
induction_curve_vectorized_3params_default<-function(x,params) {induction_curve(x,params[1],params[2],params[3],0)}
nparams.ind<-3
if (nparams.ind==3)
        {
        induction_curve<-induction_curve_3params
        induction_curve_vectorized<-induction_curve_vectorized_3params
        induction_curve_vectorized_default<-induction_curve_vectorized_3params_default
        }

loglik<-function(z) {
        py1<-sapply(xdata$time, function(t) induction_curve_vectorized_default(t,z))
        print(py1)
        logls<-sapply(1:length(py1), function(x) dbinom(xdata$y1[x],xdata$y1[x]+xdata$y2[x],prob=py1[x],log=TRUE))
        xpenalty<-sapply(z,function(y) {if (y<0) {penalty<-10^7*y^2} else {penalty<-0}; return(penalty)})
        sum(logls)-sum(xpenalty)
        }

xparms<-c(K=10^(-5),r0=10,r2=0)
res<-optimx(par=xparms,fn=function(z) loglik(z), control=list(maximize=TRUE,hessian=FALSE,gr=NULL)) 
output<-data.frame(max=unlist(res[1,1:3]),rate=names(res[1,1:3]))
write.table(output,file=paste0("params_",xFACS,".txt"),quote=FALSE,row.names=FALSE,sep="\t")
}

for fixedinduction in 5 90sec;do
for MYDELAY in '';do 
for MYINDUCTION in RNP;do 
for MYINDUCTION_CURVE in 3;do
for mymodel in modelDSBs1i1_3x4 modelDSBs1i1_nok12;do 
for MYTARGET in Psy1 CRTISO.nov2021.49and50 PhyB2.2.nov2021 Psy1.20231005.I72h CRTISO.49and50.20231005.I72h PhyB2.20231005.I72h Psy1.20230910 CRTISO.49and50.20230910 PhyB2.20231115 Psy1.20230909 CRTISO.49and50.20230909 PhyB2.20231117;do
 MYTARGETM=$( echo $MYTARGET | sed 's/_allb//' | sed 's/0h//' | sed 's/m$//' | sed 's/.no72h//' | sed 's/_mydelay0.25//' | sed 's/_R2_Feb2022//' | sed 's/.9timepoint12h//' | sed 's/.9timepoint24h//' | sed 's/.9timepoint//' | sed 's/_real0//' ) #| sed 's/.24h//g') # | sed 's/.20230903//' ) .24h should be plaed only for tries on Psy1, since matrix is the same.
ERRORMATRIX=~/workspace/daniela/error_matrices/error_matrix4_${MYTARGETM}_errorsfromunbroken.tsv
TIMECOURSE=~/workspace/daniela/input_datasets/timecourse_${MYINDUCTION}_${MYTARGET}${MYDELAY}.txt
OUTPUTFOLDER=~/workspace/daniela/FACS/${fixedinduction}/results_${mymodel}_${MYINDUCTION}_ind.c${MYINDUCTION_CURVE}_${MYTARGET}
OUTPUT=${OUTPUTFOLDER}/timecourse_rates
mkdir -p $OUTPUTFOLDER
INDUCTIONPARAMS=~/workspace/daniela/FACStable/params_FACStable.${fixedinduction}.txt
likfun=0; 
echo $OUTPUT; NAMESIMS=v4
MEMS=2000; NITER=1250;
for IITER in a b c d;do
rm ${OUTPUT}.${IITER}
if [ ! -f ${OUTPUT}.${IITER} ];then 
bsub -q new-long -J ${NAMESIMS} -e ~/mylogs/${NAMESIMS}-%J.e -o ~/mylogs/${NAMESIMS}-%J.o -R rusage[mem=${MEMS}] scripts/run_optimizationc.sh $TIMECOURSE $ERRORMATRIX ${OUTPUT}.${IITER} ${mymodel} ${NITER} $likfun ${MYINDUCTION_CURVE} ${INDUCTIONPARAMS}; fi
done;done;done;done;done;done


#------------------------------------BOOTSTRAPS----------------------------------
bootstraptype=stratifiedbootstraps
for fixedinduction in 5;do #  90sec;do
for MYDELAY in '';do 
for MYINDUCTION in RNP;do 
for MYINDUCTION_CURVE in 3;do
for mymodel in modelDSBs1i1_3x4 modelDSBs1i1_nok12;do 
for MYTARGET in Psy1 CRTISO.nov2021.49and50 PhyB2.2.nov2021 Psy1.20231005.I72h CRTISO.49and50.20231005.I72h PhyB2.20231005.I72h Psy1.20230910 CRTISO.49and50.20230910 PhyB2.20231115 Psy1.20230909 CRTISO.49and50.20230909 PhyB2.20231117;do
 MYTARGETM=$( echo $MYTARGET | sed 's/_allb//' | sed 's/0h//' | sed 's/m$//' | sed 's/.no72h//' | sed 's/_mydelay0.25//' | sed 's/_R2_Feb2022//' | sed 's/.9timepoint12h//' | sed 's/.9timepoint24h//' | sed 's/.9timepoint//' | sed 's/_real0//' ) #| sed 's/.24h//g') # | sed 's/.20230903//' ) .24h should be plaed only for tries on Psy1, since matrix is the same.
ERRORMATRIX=~/workspace/daniela/error_matrices/error_matrix4_${MYTARGETM}_errorsfromunbroken.tsv
OUTPUTFOLDER=~/workspace/daniela/FACS/${fixedinduction}/results_${mymodel}_${MYINDUCTION}_ind.c${MYINDUCTION_CURVE}_${MYTARGET}
mkdir -p $OUTPUTFOLDER
INDUCTIONPARAMS=~/workspace/daniela/FACStable/params_FACStable.${fixedinduction}.txt
likfun=0; 
echo $OUTPUT; NAMESIMS=boots
MEMS=2000; NITER=1250;
for ip in `seq 1 100`;do 
TIMECOURSEPATH=~/workspace/daniela/input_datasets/${bootstraptype}/timecourse_${MYINDUCTION}_${MYTARGET}${MYDELAY}
TIMECOURSE=${TIMECOURSEPATH}/timecourse.bs${ip}.txt
OUTPUT=${OUTPUTFOLDER}/timecourse_rates_bs${ip}.txt
rm $OUTPUT
bsub -q new-long -J ${NAMESIMS} -e ~/mylogs/${NAMESIMS}-%J.e -o ~/mylogs/${NAMESIMS}-%J.o -R rusage[mem=${MEMS}] scripts/run_optimizationc.sh $TIMECOURSE $ERRORMATRIX ${OUTPUT} ${mymodel} ${NITER} $likfun ${MYINDUCTION_CURVE} ${INDUCTIONPARAMS}
done
done;done;done;done;done;done
done


#########################################################################################
################################### CALCULATE CI ########################################
#########################################################################################

#-------------------------------------SINGLE RUNS----------------------------------------

MYN=50000
for fixedinduction in 5;do # 90sec;do
for MYDELAY in '';do 
for MYINDUCTION in RNP;do 
for MYINDUCTION_CURVE in 3;do
for mymodel in modelDSBs1i1_3x4 modelDSBs1i1_nok12;do 
for MYTARGET in Psy1 CRTISO.nov2021.49and50 PhyB2.2.nov2021 Psy1.20231005.I72h CRTISO.49and50.20231005.I72h PhyB2.20231005.I72h Psy1.20230910 CRTISO.49and50.20230910 PhyB2.20231115 Psy1.20230909 CRTISO.49and50.20230909 PhyB2.20231117;do
 MYTARGETM=$( echo $MYTARGET | sed 's/_allb//' | sed 's/0h//' | sed 's/m$//' | sed 's/.no72h//' | sed 's/_mydelay0.25//' | sed 's/_R2_Feb2022//' | sed 's/.9timepoint12h//' | sed 's/.9timepoint24h//' | sed 's/.9timepoint//' | sed 's/_real0//' ) #| sed 's/.24h//g') # | sed 's/.20230903//' ) .24h should be plaed only for tries on Psy1, since matrix is the same.
ERRORMATRIX=~/workspace/daniela/error_matrices/error_matrix4_${MYTARGETM}_errorsfromunbroken.tsv
TIMECOURSE=~/workspace/daniela/input_datasets/timecourse_${MYINDUCTION}_${MYTARGET}${MYDELAY}.txt
OUTPUTFOLDER=~/workspace/daniela/FACS/${fixedinduction}/results_${mymodel}_${MYINDUCTION}_ind.c${MYINDUCTION_CURVE}_${MYTARGET}
OUTPUT=${OUTPUTFOLDER}/timecourse_rates
OUTPUTALL=${OUTPUT}.all.tsv
mkdir -p $OUTPUTFOLDER
INDUCTIONPARAMS=~/workspace/daniela/FACStable/params_FACStable.${fixedinduction}.txt
likfun=0; if [ $mymodel == "modelDSBs1i1_3x4" ];then likfun=3; fi
echo $OUTPUT; NAMESIMS=v4
MEMS=5000; NITER=1250;
if [ -f ${OUTPUT}.a ];then
	echo $MYN $TIMECOURSE
	cat ${OUTPUT}.a | head -1 > temp.tab
	for IITER in a b c d;do
		cat ${OUTPUT}.${IITER} | grep -v k11 >> temp.tab
	done
	cat temp.tab | head -5001 > $OUTPUTALL
	bsub -q risk -J CIs -e ~/mylogs2/confint.e%J -o ~/mylogs2/confint.o%J -R rusage[mem=1000] ./scripts/calculate_CI_v1c.sh $OUTPUTALL $TIMECOURSE $ERRORMATRIX $MYN $mymodel $likfun $MYINDUCTION_CURVE ${INDUCTIONPARAMS} $INDUCTIONPARAMS
fi
done;done;done;done;done;done

#------------------------------------BOOTSTRAPS----------------------------------
bootstraptype=stratifiedbootstraps
MYN=50000
for fixedinduction in 5;do # 90sec;do
for MYDELAY in '';do 
for MYINDUCTION in RNP;do 
for MYINDUCTION_CURVE in 3;do
for mymodel in modelDSBs1i1_3x4 modelDSBs1i1_nok12;do 
for MYTARGET in Psy1 CRTISO.nov2021.49and50 PhyB2.2.nov2021 Psy1.20231005.I72h CRTISO.49and50.20231005.I72h PhyB2.20231005.I72h Psy1.20230910 CRTISO.49and50.20230910 PhyB2.20231115 Psy1.20230909 CRTISO.49and50.20230909 PhyB2.20231117;do
 MYTARGETM=$( echo $MYTARGET | sed 's/_allb//' | sed 's/0h//' | sed 's/m$//' | sed 's/.no72h//' | sed 's/_mydelay0.25//' | sed 's/_R2_Feb2022//' | sed 's/.9timepoint12h//' | sed 's/.9timepoint24h//' | sed 's/.9timepoint//' | sed 's/_real0//' ) #| sed 's/.24h//g') # | sed 's/.20230903//' ) .24h should be plaed only for tries on Psy1, since matrix is the same.
ERRORMATRIX=~/workspace/daniela/error_matrices/error_matrix4_${MYTARGETM}_errorsfromunbroken.tsv
OUTPUTFOLDER=~/workspace/daniela/FACS/${fixedinduction}/results_${mymodel}_${MYINDUCTION}_ind.c${MYINDUCTION_CURVE}_${MYTARGET}
mkdir -p $OUTPUTFOLDER
INDUCTIONPARAMS=~/workspace/daniela/FACStable/params_FACStable.${fixedinduction}.txt
likfun=0; 
echo $OUTPUT; NAMESIMS=boots
MEMS=5000; NITER=1250;
for ip in `seq 1 100`;do 
TIMECOURSEPATH=~/workspace/daniela/input_datasets/${bootstraptype}/timecourse_${MYINDUCTION}_${MYTARGET}${MYDELAY}
TIMECOURSE=${TIMECOURSEPATH}/timecourse.bs${ip}.txt
OUTPUT=${OUTPUTFOLDER}/timecourse_rates_bs${ip}.txt
cat ${OUTPUT} | head -1 > ${OUTPUT}.parsed.tab 
cat ${OUTPUT} | grep -v k11 >> ${OUTPUT}.parsed.tab
bsub -q risk -J CI -e ~/mylogs3/confint.e%J -o ~/mylogs3/confint.o%J -R rusage[mem=1000] ./scripts/calculate_CI_v1c.sh ${OUTPUT}.parsed.tab $TIMECOURSE $ERRORMATRIX $MYN $mymodel $likfun $MYINDUCTION_CURVE ${INDUCTIONPARAMS}
done
done;done;done;done;done;done
#done


#########################################################################################
################################### PLOTTING / FLOW #####################################
#########################################################################################

#-------------------------------------SINGLE RUNS----------------------------------------

#in plot use z!=1 not to normalize flow on induction curve (with new induction function)
TIMERESOLUTION=0.002
for TIMEFLOW in 24 72;do 
for MYDELAY in '';do #_mydelay0.25 .no72h;do #''
for MYINDUCTION in RNP;do #RNP
for MYINDUCTION_CURVE in 3; do #3 gRNA;do
for mymodel in modelDSBs1i1_nok12 modelDSBs1i1_3x4;do 
for MYTARGET in Psy1 CRTISO.nov2021.49and50 PhyB2.2.nov2021 Psy1.20231005.I72h CRTISO.49and50.20231005.I72h PhyB2.20231005.I72h Psy1.20230910 CRTISO.49and50.20230910 PhyB2.20231115 Psy1.20230909 CRTISO.49and50.20230909 PhyB2.20231117;do
MYTARGETM=$( echo $MYTARGET | sed 's/_allb//' |  sed 's/_all//' | sed 's/_R2_Feb2022//' | sed 's/.no72h//' | sed 's/_R2_Feb2022//' | sed 's/_mydelay0.25//' |  sed 's/.9timepoint24h//'  | sed 's/_real0//g') # | sed 's/.24h//g' )
ERRORMATRIX=~/workspace/daniela/error_matrices/error_matrix4_${MYTARGETM}_errorsfromunbroken.tsv
TIMECOURSE=~/workspace/daniela/input_datasets/timecourse_${MYINDUCTION}_${MYTARGET}${MYDELAY}.txt
likfun=0; if [ $mymodel == "modelDSBs1i1_3x4" ] || [ $mymodel == "modelDSBs1i1_3x4nor11" ];then likfun=3; fi
if [ $mymodel == "model5i1" ];then
ERRORMATRIX=~/workspace/daniela/error_matrices/error_matrix3_${MYTARGETM}_errorsfromunbroken.tsv
TIMECOURSE=~/workspace/daniela/input_datasets/timecourse_${MYINDUCTION}_${MYTARGET}${MYDELAY}_3states.txt
fi
OUTPUTFOLDER=~/workspace/daniela/FACS/5/results_${mymodel}_${MYINDUCTION}_ind.c${MYINDUCTION_CURVE}_${MYTARGET}
OUTPUT=${OUTPUTFOLDER}/timecourse_rates
if [ -f ${OUTPUT}.a ];then
        myfile="notexist"
        if [ -f ${OUTPUT}.all_n50000.CI.RData ];then myfile=${OUTPUT}.all_n50000.CI.RData;fi
        ERRORMATRIXM=$( basename $ERRORMATRIX .tsv )
        TIMECOURSEM=$( basename $TIMECOURSE .txt )
        if [ -f ${myfile} ];then
		echo $myfile
	bsub -q new-short -J plotshort -e ~/mylogs/plot.e%J -o ~/mylogs/plot.o%J -R rusage[mem=2000] ./scripts/plot.v2.sh $myfile $TIMECOURSE $ERRORMATRIX $mymodel 1000 ${OUTPUTFOLDER}/plot_${mymodel}_${MYINDUCTION}_ind.c${MYINDUCTION_CURVE}_${MYTARGET}_${TIMEFLOW} $likfun ${MYINDUCTION_CURVE} ${TIMERESOLUTION} 1 1 ${TIMEFLOW}
	if [ $mymodel == "modelDSBs1i1_3x4" ] || [ $mymodel == "modelDSBs1i1_3x4nor11" ];then
	bsub -q new-short -J plotshort -e ~/mylogs/plot.e%J -o ~/mylogs/plot.o%J -R rusage[mem=2000] ./scripts/plot.v2.sh $myfile $TIMECOURSE $ERRORMATRIX $mymodel 1000 ${OUTPUTFOLDER}/plot_${mymodel}_${MYINDUCTION}_ind.c${MYINDUCTION_CURVE}_${MYTARGET}_${TIMEFLOW}_noprocessed $likfun ${MYINDUCTION_CURVE} ${TIMERESOLUTION} 1 0 ${TIMEFLOW}
	fi
        fi
fi
done;done;done;done;done;done

#-------------------------------------BOOTSTRAPS----------------------------------------
#/home/labs/alevy/fabrizio/workspace/daniela/resultsv4/stratifiedbootstraps/timecourse_RNP_PhyB2.2_allb/modelDSBs1i1_mini/results_modelDSBs1i1_mini_RNP_ind.c2_PhyB2.2_allb/


fixedinduction=5
TIMEFLOWS=$( echo 24 72 ) #24 48 72 )
TIMERESOLUTION=0.002
#stratifiedbootstraps2 goes with all and modelDSBs1i1_mini.bytarget; stratifiedbootstraps with allb
for MYDELAY in '';do #_mydelay0.25 .no72h;do #''
for MYINDUCTION in RNP; do # RNP;do
for MYINDUCTION_CURVE in 3;do # 3 4; do # 2 4
for mymodel in modelDSBs1i1_nok12 modelDSBs1i1_3x4;do 
for MYTARGET in Psy1 CRTISO.nov2021.49and50 PhyB2.2.nov2021 Psy1.20231005.I72h CRTISO.49and50.20231005.I72h PhyB2.20231005.I72h Psy1.20230910 CRTISO.49and50.20230910 PhyB2.20231115 Psy1.20230909 CRTISO.49and50.20230909 PhyB2.20231117;do
for bootstraptype in stratifiedbootstraps;do 
for TIMEFLOW in $TIMEFLOWS;do
for ip in `seq 1 100`;do #2 34 68 82 99
MYTARGETM=$( echo $MYTARGET | sed 's/_allb//' |  sed 's/_all//' | sed 's/_R2_Feb2022//' | sed 's/_R2_Feb2022//' | sed 's/_mydelay0.25//'  | sed 's/.9timepoint24h//'  | sed 's/_real0//g' ) # | sed 's/.24h//g' )
ERRORMATRIX=~/workspace/daniela/error_matrices/error_matrix4_${MYTARGETM}_errorsfromunbroken.tsv
TIMECOURSEPATH=~/workspace/daniela/input_datasets/${bootstraptype}/timecourse_${MYINDUCTION}_${MYTARGET}${MYDELAY}
TIMECOURSE=${TIMECOURSEPATH}/timecourse.bs${ip}.txt
OUTPUTA=~/workspace/daniela/resultsv5/${bootstraptype}/timecourse_${MYINDUCTION}_${MYTARGET}${MYDELAY}/${mymodel}
OUTPUT=~/workspace/daniela/resultsv5/${bootstraptype}/timecourse_${MYINDUCTION}_${MYTARGET}${MYDELAY}/${mymodel}/results_${mymodel}_${MYINDUCTION}_ind.c${MYINDUCTION_CURVE}_${MYTARGET}${MYDELAY}
OUTPUTFOLDER=~/workspace/daniela/FACS/${fixedinduction}/results_${mymodel}_${MYINDUCTION}_ind.c${MYINDUCTION_CURVE}_${MYTARGET}
TIMECOURSEPATH=~/workspace/daniela/input_datasets/${bootstraptype}/timecourse_${MYINDUCTION}_${MYTARGET}${MYDELAY}
TIMECOURSE=${TIMECOURSEPATH}/timecourse.bs${ip}.txt
OUTPUTA=${OUTPUTFOLDER}/timecourse_rates_bs${ip}.txt
mkdir -p $OUTPUT;
OUTPUTTOT=${OUTPUTFOLDER}/results.tot
likfun=0; if [ $mymodel == "modelDSBs1i1_3x4" ];then likfun=3; fi
#echo $OUTPUTA; NAMESIMS=boot2 #MEboot #${MYTARGET};
MEMS=18000; NITER=1250;
echo $ip
        myfile="notexist"
        if [ -f ${OUTPUTA}.parsed.tab_n20000.CI.RData ];then myfile=${OUTPUTA}.parsed.tab_n20000.CI.RData;fi
        if [ -f ${OUTPUTA}.parsed.tab_n50000.CI.RData ];then myfile=${OUTPUTA}.parsed.tab_n50000.CI.RData;fi
        ERRORMATRIXM=$( basename $ERRORMATRIX .tsv )
        TIMECOURSEM=$( basename $TIMECOURSE .txt )
        if [ -f ${myfile} ];then
                myfileCI=$( echo $myfile | sed 's/.CI.RData/.CI/' )
#                cat ${myfileCI} | awk -v var1=$MYTARGET -v var2=$MYINDUCTION -v var3=$mymodel -v var4=$ERRORMATRIXM -v var5=$TIMECOURSEM '{print $0,var1,var2,var3,var4,var5}' | grep -v max >> ~/workspace/daniela/resultsv5/table.CI
#bsub -q new-short -J plotshort -e ~/mylogs4/plot.e%J -o ~/mylogs4/plot.o%J -R rusage[mem=1500] ./scripts/plot.v2.sh $myfile $TIMECOURSE $ERRORMATRIX $mymodel 100 ${OUTPUTA}.${TIMEFLOW}h.plot $likfun ${MYINDUCTION_CURVE} ${TIMERESOLUTION} 1 1 ${TIMEFLOW}
#	if [ $ip -eq 1 ]; then rm ${OUTPUTTOT};fi
#	cat ${OUTPUTA} | awk -v var1=$ip -v OFS='\t' '{if (NR==1){print $0,"bootstrap"} else {print $0,var1}}' >> ${OUTPUTTOT}
        fi
done
#####build flow tables
cat ${OUTPUTFOLDER}/timecourse_rates_bs*.txt.${TIMEFLOW}h.plot_plot.flow.tab > ${OUTPUTFOLDER}/results.totall_${TIMEFLOW}h.plot.flow.tab 
echo ${OUTPUTFOLDER}/results.totall_${TIMEFLOW}h.plot.flow.tab
done #timeflow
done;done;done;done;done;done

for mymodel in modelDSBs1i1_3x4 modelDSBs1i1_nok12;do 
for mytarget in Psy1 CRTISO.nov2021.49and50 PhyB2.2.nov2021 Psy1.20231005.I72h CRTISO.49and50.20231005.I72h PhyB2.20231005.I72h Psy1.20230910 CRTISO.49and50.20230910 PhyB2.20231115 Psy1.20230909 CRTISO.49and50.20230909 PhyB2.20231117;do
myfolder=~/workspace/daniela/FACS/5/results_${mymodel}_${MYINDUCTION}_ind.c${MYINDUCTION_CURVE}_${mytarget}
rm ${myfolder}/results.tot
echo $myfolder
for i in `seq 1 100`;do
cat ${myfolder}/timecourse_rates_bs${i}.txt | awk -v var=${i} -v OFS='\t' '{if (NR==1){print $0,"bootstrap"} else {print $0,var}}' >> ${myfolder}/results.tot
done
done
done



TIMEFLOWS=$( echo 24 72 ) #24 48 72 )
TIMERESOLUTION=0.002
#stratifiedbootstraps2 goes with all and modelDSBs1i1_mini.bytarget; stratifiedbootstraps with allb
for MYDELAY in '';do 
for MYINDUCTION in RNP; do 
for MYINDUCTION_CURVE in 3;do 
for mymodel in modelDSBs1i1_nok12 modelDSBs1i1_3x4;do #modelDSBs1i1_3x4nor11
for MYTARGET in Psy1 CRTISO.nov2021.49and50 PhyB2.2.nov2021 Psy1.20231005.I72h CRTISO.49and50.20231005.I72h PhyB2.20231005.I72h Psy1.20230910 CRTISO.49and50.20230910 PhyB2.20231115 Psy1.20230909 CRTISO.49and50.20230909 PhyB2.20231117;do
for bootstraptype in stratifiedbootstraps;do 
for TIMEFLOW in $TIMEFLOWS;do
MYTARGETM=$( echo $MYTARGET | sed 's/_allb//' |  sed 's/_all//' | sed 's/_R2_Feb2022//' | sed 's/_R2_Feb2022//' | sed 's/_mydelay0.25//'  | sed 's/.9timepoint24h//'  | sed 's/_real0//g' ) # | sed 's/.24h//g' )
ERRORMATRIX=~/workspace/daniela/error_matrices/error_matrix4_${MYTARGETM}_errorsfromunbroken.tsv
TIMECOURSEPATH=~/workspace/daniela/input_datasets/${bootstraptype}/timecourse_${MYINDUCTION}_${MYTARGET}${MYDELAY}
TIMECOURSE=${TIMECOURSEPATH}/timecourse.bs${ip}.txt
OUTPUT=~/workspace/daniela/FACS/5/results_${mymodel}_${MYINDUCTION}_ind.c${MYINDUCTION_CURVE}_${MYTARGET}
#~/workspace/daniela/resultsv5/${bootstraptype}/timecourse_${MYINDUCTION}_${MYTARGET}${MYDELAY}/${mymodel}/results_${mymodel}_${MYINDUCTION}_ind.c${MYINDUCTION_CURVE}_${MYTARGET}${MYDELAY}
mkdir -p $OUTPUT;
OUTPUTTOT=${OUTPUT}/results.tot
OUTPUTA=${OUTPUT}/results_${ip}.txt
likfun=0; if [ $mymodel == "modelDSBs1i1_3x4" ];then likfun=3; fi
echo $OUTPUTA; NAMESIMS=boot2 #MEboot #${MYTARGET};
MEMS=4000; NITER=1250;
#####plots
TIMECOURSE=~/workspace/daniela/input_datasets/timecourse_${MYINDUCTION}_${MYTARGET}.txt
listCIFILE=${OUTPUT}/listCIfiles.txt
myfiles=$( ls ${OUTPUT}/timecourse_rates.all_n50000.CI ${OUTPUT}/timecourse_rates_bs*.txt.parsed.tab_n50000.CI )
for ifile in $myfiles;do
cat $ifile | awk -v OFS='\t' '{if (ar[$4]==1){$4=$4".1"}; ar[$4]=1; print $0}' > ${ifile}.t.CI 
echo $ifile
done
ls ${OUTPUT}/timecourse_rates.all_n50000.CI.t.CI > $listCIFILE 
ls ${OUTPUT}/timecourse_rates_bs*.txt.parsed.tab_n50000.CI.t.CI >> $listCIFILE
bsub -q new-short -J plotshort -e ~/mylogs4/plot.e%J -o ~/mylogs4/plot.o%J -R rusage[mem=16000] ./scripts/plot.v2.sh $listCIFILE $TIMECOURSE $ERRORMATRIX $mymodel 1000 ${OUTPUTTOT}.${TIMEFLOW}h.plot $likfun ${MYINDUCTION_CURVE} ${TIMERESOLUTION} 1 1 ${TIMEFLOW}
        if [ $mymodel == "modelDSBs1i1_3x4" ] || [ $mymodel == "modelDSBs1i1_3x4nor11" ];then
        bsub -q new-short -J plotshort -e ~/mylogs4/plot.e%J -o ~/mylogs4/plot.o%J -R rusage[mem=16000] ./scripts/plot.v2.sh $listCIFILE $TIMECOURSE $ERRORMATRIX $mymodel 1000 ${OUTPUTTOT}.${TIMEFLOW}h.plot.noproc $likfun ${MYINDUCTION_CURVE} ${TIMERESOLUTION} 1 0 ${TIMEFLOW}
        fi
done;done;done;done;done;done;done




