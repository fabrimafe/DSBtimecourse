fixedinduction=5
TIMEFLOWS=$( echo 24 72 ) #$( echo 72 ) #24 48 72 )
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
                cat ${myfileCI} | awk -v var1=$MYTARGET -v var2=$MYINDUCTION -v var3=$mymodel -v var4=$ERRORMATRIXM -v var5=$TIMECOURSEM '{print $0,var1,var2,var3,var4,var5}' | grep -v max >> ~/workspace/daniela/resultsv5/table.CI
bsub -q new-short -J plotshort -e ~/mylogs4/plot.e%J -o ~/mylogs4/plot.o%J -R rusage[mem=4000] ./scripts/plot.v2.sh $myfile $TIMECOURSE $ERRORMATRIX $mymodel 100 ${OUTPUTA}.${TIMEFLOW}h.plot $likfun ${MYINDUCTION_CURVE} ${TIMERESOLUTION} 1 1 ${TIMEFLOW}
#       if [ $ip -eq 1 ]; then rm ${OUTPUTTOT};fi
#       cat ${OUTPUTA} | awk -v var1=$ip -v OFS='\t' '{if (NR==1){print $0,"bootstrap"} else {print $0,var1}}' >> ${OUTPUTTOT}
        fi
done
#####build flow tables
cat ${OUTPUTFOLDER}/results_*.txt.${TIMEFLOW}h.plot_plot.flow.tab > ${OUTPUTFOLDER}/results.totall_${TIMEFLOW}h.plot.flow.tab
echo ${OUTPUTFOLDER}/results.totall_${TIMEFLOW}h.plot.flow.tab
done #timeflow
done;done;done;done;done;done

