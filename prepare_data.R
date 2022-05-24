#!/usr/bin/env Rscript
########################################################################################################
##################### LOAD DATA ########################################################################
########################################################################################################
setwd("~/workspace/daniela")
#mydelay<-0
mydelay<-0
ntypes<-4
library(tidyverse)
#errormatrix<-as.matrix(read.table(myerrorE, header=FALSE))
#load("RData/error_matrices4_l.RData")
#error_mat_l<-error_matrices4_l
#if (mymodel=="model5i1" || mymodel=="model5i1_nor11")
#        {
#        ntypes<-3 ; load("RData/error_matrices3_l.RData");
#        error_mat_l<-error_matrices3_l
#        }

#if (optimize_errorDSB2indel==2) { mydelay<-0 }

mydata_l<-list()
mydata.p_l<-list()
#in addition CRTISO0h0 was edited manually by introducing 0 counts at time 0
#
for (mytarget_i in c("PhyB2")) #c("Psy1_R2_Feb2022","PhyB2.2_R2_Feb2022")) #, "CRTISO.cleanedandnov2021.49and50","CRTISO.nov2021.49and50")) #c("CRTISO.cleanedandnov2021","CRTISO.nov2021","PhyB2.2.nov2021")) #c("CRTISOcleaned")) #,"CRTISO.49and50bp","CRTISO","Psy1","PhyB2.1","PhyB2.2","PhyB2.3"))
        {
        for (myinduction in c("gRNA")) #,"RNP"))
                {
#                mytarget<-mytarget_i
                    nophyB2<-(myinduction=="RNP" || (myinduction=="gRNA" && ( mytarget_i!="PhyB2.1" && mytarget_i!="PhyB2.2" && mytarget_i!="PhyB2.3" )))
		    print(paste0("~/workspace/daniela/csv/",mytarget_i,"_",myinduction,"_Types_MH_df.csv"))
                    if (nophyB2)
                            {
                            mydata<-read.csv(paste0("~/workspace/daniela/csv/",mytarget_i,"_",myinduction,"_Types_MH_df.csv"))
                            if (mytarget_i=="CRTISO" && myinduction=="gRNA") { mydata<-mydata[mydata$X!="CRTISO_6h_A",]};
                            if (mytarget_i=="CRTISO" && myinduction=="RNP") { mydata<-mydata[mydata$X!="CRTISO_0h_B",];};
#mydata<-mydata[1:11,]};
                            mydata[is.na(mydata)]<-0; mydata<-cbind(matrix(unlist(sapply(as.character(mydata$X),function(x) strsplit(x,"_"))),ncol=3,byrow=T),mydata); names(mydata)[1:4]<-c("target","time","replicate","ID_name");
                            mydata$time<-as.numeric(sapply(mydata$time, function(x) gsub("h","",x)))
                            mydata$time<-mydata$time+mydelay
                            mydata<-mydata %>% mutate(preciseDSB=Perfect.DSB,impreciseDSB=Extended.DSB+PAM_side.DSB+guide_side.DSB,indels=Ins+Del+MH_Del+MH_Del_2)
                            mydata<- mydata %>% select(time,WT.Sub,indels,preciseDSB,impreciseDSB)
                            names(mydata)<-c("time","y1","y2","y3","y4")
                            mydata.p<-as.data.frame(cbind(mydata$time,t(apply(mydata[,2:5],MARGIN=1,FUN=function(x) x/sum(x)))))
                            names(mydata.p)<-c("time","y1","y2","y3","y4")
                            if (ntypes==3)
                                    {
                                    mydata.p<-mydata.p %>% mutate(y5=y2,y3=y3+y4,y4=NULL) %>% mutate(y2=y3,y3=y5,y5=NULL)
                                    mydata<-mydata %>% mutate(y5=y2,y3=y3+y4,y4=NULL) %>% mutate(y2=y3,y3=y5,y5=NULL)
                                    }
                        mydata_l[[paste0(myinduction,"_",mytarget_i)]]<-mydata
                        mydata.p_l[[paste0(myinduction,"_",mytarget_i)]]<-mydata.p
                        if (mydelay==0){
			write.table(mydata,file=paste0("~/workspace/daniela/input_datasets/timecourse_",myinduction,"_",mytarget_i,".txt"),quote=FALSE, col.names=TRUE,row.names=FALSE,sep="\t")} else 
			{write.table(mydata,file=paste0("~/workspace/daniela/input_datasets/timecourse_",myinduction,"_",mytarget_i,"_mydelay",mydelay,".txt"),quote=FALSE, col.names=TRUE,row.names=FALSE,sep="\t")}
                        }
                }
        }

