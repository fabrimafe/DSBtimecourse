library(data.table)
library(tidyverse)

xdata<-fread("/home/labs/alevy/fabrizio/workspace/daniela/resultsv3/MEbootstraps/timecourse_RNP_CRTISO.nov2021.49and50/modelDSBs1i1_3x4/results.tot",header=TRUE)
xdata<-fread("/home/labs/alevy/fabrizio/workspace/daniela/resultsv3/MEbootstraps/timecourse_RNP_CRTISO.nov2021.49and50/modelDSBs1i1_mini/results.tot",header=TRUE)

xdata<-fread("/home/labs/alevy/fabrizio/workspace/daniela/resultsv3/MEbootstraps/timecourse_RNP_Psy1/modelDSBs1i1_3x4/results.tot",header=TRUE)
xdata<-fread("/home/labs/alevy/fabrizio/workspace/daniela/resultsv3/MEbootstraps/timecourse_RNP_Psy1/modelDSBs1i1_mini/results.tot",header=TRUE)


xdata<-fread("/home/labs/alevy/fabrizio/workspace/daniela/resultsv3/stationarybootstraps/timecourse_RNP_Psy1/modelDSBs1i1_3x4/results.tot",header=TRUE)


xdata<-fread("/home/labs/alevy/fabrizio/workspace/daniela/resultsv3/stationarybootstraps/timecourse_RNP_Psy1/modelDSBs1i1_mini/results.tot",header=TRUE)
xdata<-fread("/home/labs/alevy/fabrizio/workspace/daniela/resultsv3/MEbootstraps/timecourse_RNP_Psy1/modelDSBs1i1_mini/results.tot",header=TRUE)

xdata<-fread("/home/labs/alevy/fabrizio/workspace/daniela/resultsv3/MEbootstraps/timecourse_RNP_Psy1x2/modelDSBs1i1_mini/results.tot",header=TRUE)
myfile="/home/labs/alevy/fabrizio/workspace/daniela/resultsv4/MEbootstraps/timecourse_RNP_CRTISO.cleanedandnov2021/modelDSBs1i1_nok12/results_modelDSBs1i1_nok12_RNP_ind.c3_CRTISO.cleanedandnov2021/results.tot"

xdata<-fread(


myfolder=/home/labs/alevy/fabrizio/workspace/daniela/resultsv3/MEbootstraps/timecourse_RNP_Psy1x2/modelDSBs1i1_3x4
myfolder=/home/labs/alevy/fabrizio/workspace/daniela/resultsv3/MEbootstraps/timecourse_RNP_Psy1.MockEarlyTimePoints/modelDSBs1i1_3x4
myfolder=/home/labs/alevy/fabrizio/workspace/daniela/resultsv3/stationarybootstraps/timecourse_RNP_Psy1.MockEarlyTimePoints/modelDSBs1i1_3x4

myfolder=/home/labs/alevy/fabrizio/workspace/daniela/resultsv3/stationarybootstraps/timecourse_RNP_Psy1x2/modelDSBs1i1_3x4
rm results.tot
for i in `seq 1 100`;do
cat ${myfolder}/results_${i}.txt | awk -v var=${i} -v OFS='\t' '{if (NR==1){print $0,"bootstrap"} else {print $0,var}}' >> results.tot
done


bootstrap_type="stratifiedbootstraps"  #    "stationarybootstraps" #"MEbootstraps"
mymodel="modelDSBs1i1_3x4" # "modelDSBs1i1_3x4" #"modelDSBs1i1_mini" #"modelDSBs1i1_nok12"
ninduction<-2
myfile=paste0("/home/labs/alevy/fabrizio/workspace/daniela/resultsv4/",bootstrap_type,"/timecourse_RNP_CRTISO.cleanedandnov2021/",mymodel,"/results_",mymodel,"_RNP_ind.c",ninduction,"_CRTISO.cleanedandnov2021/results.tot")

myfile="results.tot"
library(data.table)
library(tidyverse)
xdata<-fread(myfile,header=TRUE)
#xdata<-xdata %>% mutate(perfectrepair=as.numeric(r11)+as.numeric(r21)) %>% relocate(perfectrepair, .before = value)
xdata<-xdata %>% filter(er1<0.3)
maxl<- xdata  %>% group_by(bootstrap) %>% summarize(maxl=max(value)) %>% ungroup
maxls<-left_join(xdata,maxl) %>% filter(value==maxl,k11!="k11",as.numeric(maxl)>-100000)
res.all<-apply(maxls[,1:which(names(maxls)=="value")],MARGIN=2,FUN=function(x) quantile(as.numeric(x),c(0.01,0.05,0.10,0.25,0.5,0.95,0.99)))
res.all

res.filtered<-apply(maxls[,1:which(names(maxls)=="value")],MARGIN=2,FUN=function(x) quantile(as.numeric(x),c(0.01,0.05,0.5,0.95,0.99)))

