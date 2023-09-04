#!/usr/bin/env Rscript
########################################################################################################
##################### LOAD DATA ########################################################################
########################################################################################################
library(argparser, quietly=TRUE,warn.conflicts=FALSE)

argv<- arg_parser("Parse arguments")

argv <- add_argument(argv, "-i", help="input csv file", default="input.csv")
argv <- add_argument(argv, "-o", help="output file", default="output.txt")

args <- parse_args(argv)
input.file<-args$i
output.file<-args$o

#mydelay<-0
mydelay<-0
ntypes<-4
library(tidyverse)

mydata_l<-list()
mydata.p_l<-list()

mydata<-read.csv(input.file) #paste0("~/workspace/daniela/csv/",mytarget_i,"_",myinduction,"_Types_MH_df.csv"))
print("csv imported")
mydata[is.na(mydata)]<-0; mydata<-cbind(matrix(unlist(sapply(as.character(mydata$X),function(x) strsplit(x,"_"))),ncol=3,byrow=T),mydata); names(mydata)[1:4]<-c("target","time","replicate","ID_name");
mydata$time<-as.numeric(sapply(mydata$time, function(x) gsub("h","",x)))
print("add delay")
mydata$time<-mydata$time+mydelay
print("merge cols")
mydata<-mydata %>% mutate(preciseDSB=Perfect.DSB,impreciseDSB=Extended.DSB+PAM_side.DSB+guide_side.DSB,indels=Ins+Del+MH_Del+MH_Del_2)
print("select cols")
mydata<- mydata %>% select(time,WT.Sub,indels,preciseDSB,impreciseDSB)
names(mydata)<-c("time","y1","y2","y3","y4")
print("normalize")
mydata.p<-as.data.frame(cbind(mydata$time,t(apply(mydata[,2:5],MARGIN=1,FUN=function(x) x/sum(x)))))
names(mydata.p)<-c("time","y1","y2","y3","y4")
if (ntypes==3)
	{
	mydata.p<-mydata.p %>% mutate(y5=y2,y3=y3+y4,y4=NULL) %>% mutate(y2=y3,y3=y5,y5=NULL)
	mydata<-mydata %>% mutate(y5=y2,y3=y3+y4,y4=NULL) %>% mutate(y2=y3,y3=y5,y5=NULL)
        }
print("write file")
write.table(mydata,file=output.file,quote=FALSE, col.names=TRUE,row.names=FALSE,sep="\t")

#if (mydelay==0)
#	{
#	write.table(mydata,file=paste0("~/workspace/daniela/input_datasets/timecourse_",myinduction,"_",mytarget_i,".txt"),quote=FALSE, col.names=TRUE,row.names=FALSE,sep="\t")
#        } else 
#	{
#	write.table(mydata,file=paste0("~/workspace/daniela/input_datasets/timecourse_",myinduction,"_",mytarget_i,"_mydelay",mydelay,".txt"),quote=FALSE, col.names=TRUE,row.names=FALSE,sep="\t")
#	}

