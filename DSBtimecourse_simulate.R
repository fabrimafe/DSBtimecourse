#!/usr/bin/env Rscript
#example to run
#./optimize_model_backbone.R  -t CRTISO -i RNP -m modelDSBs1i1_fullimpreciseDSB -e E_errorsfromunbroken -n 10000 -o CRTISO_RNP_optimization.tsv 

library(argparser, quietly=TRUE,warn.conflicts=FALSE)
library(tidyverse)
library(optimx)
library(deSolve)

################################################################################
##################IMPORT FUNCTIONS##############################################
################################################################################

source("likelihood_functions.R")

################################################################################
################## READ ARGUMENTS ##############################################
################################################################################

argv<- arg_parser("Parse arguments")

argv <- add_argument(argv, "-T", help="time course of target site. A dataset with time as 1st column and then the number of molecules")
argv <- add_argument(argv, "-m", help="model. One of modelDSBs1i1_3x4,modelDSBs1i1_nok12,modelDSBs1i1_mini,modelDSBs1i1_realimprecise,modelDSBs1i1_3x4nor11.")
argv <- add_argument(argv, "-o", help="output file", default="output.txt")
#argv <- add_argument(argv, "-e", help="errors")
argv <- add_argument(argv, "-E", help="error matrix. A tab separated matrix specifying in the rows the source type and in the column the observed type, e.g. the second row/first column indicates the proportion of y2 (precise DSB) that are observed as intact molecules (y1)")
argv <- add_argument(argv, "-z", help="n parameters in induction curve.", default=3)
argv <- add_argument(argv, "-p", help="parameters. A dataframe with a column 'max' specifying the rate parameters, and a column 'rate' specifying the name of the parameters") #e.g. ~/workspace/daniela/resultsv5/results_modelDSBs1i1_realimprecise_RNP_ind.c3_Psy1_R2_Feb2022.all_n500000.CI

args <- parse_args(argv)
input.file<-args$T
#myerror<-args$e
myerrorE<-args$E
mymodel<-args$m
nparamsind<-as.numeric(args$z)
output.file<-args$o
n.max<-as.numeric(args$n)
inputparameters<-args$p

#if (is.na(myerror)){ myerror<-"error" }
if (is.na(output.file)){ output.file<-"DSBtimecourse_optimize.tsv" }
if (is.na(nparamsind)){ nparamsind<-2 }
define_ODE.functions(nparamsind)


########################################################################################################
##################### IMPORT DATA ######################################################################
########################################################################################################
#INPUT PARAMETERS
mydelay<-0
ntypes<-4
print("read input parameters")
myinputparameters<-read.table(inputparameters, header=TRUE)

#INPUT TIMECOURSE (used to define the time points and the number of molecules for each time point)
print("read input timecourse")
mydata<-read.table(input.file, header=TRUE)
time_courses_begins<-c(which(mydata$time[-1]-mydata$time[-length(mydata$time)]<0)+1)
nameparms<-model2nameparams(mymodel,nparamsind,0) #model2nameparams set the parameter names for the input model
errormatrix<-as.matrix(read.table(myerrorE, header=FALSE))
ntypes<-model2ntypes(mymodel)
xmodel<-model2xmodel(mymodel)
nrates<-length(nameparms)
print(ntypes)

#################################################################################
####################### PREDICT MODEL ###########################################
#################################################################################
print("input parameters are:")
bestmodels_t<-myinputparameters %>% select(rate,max) %>% pivot_wider(names_from=rate,values_from=max) %>% ungroup
print(bestmodels_t)
print("predicted proportions are:")
res<-predict_models(bestmodels_t,ntypes=ntypes,nparms=nrates,nheaders=0,errormatrix=errormatrix,mymodel=xmodel,timest=mydata$time)
print(res)
print("add multinomial sampling noise")

#Function that sample using a multinomial distribution given an expected proportion of molecules of each given type
prop2nmolecules<-function(xpoint)
{
xtime<-mydata$time[xpoint]
nmolecules<-sum(mydata[xpoint,2:ncol(mydata)])
ptypes<-res %>% filter(time==xtime) %>% head(ntypes) %>% select(p) %>% unlist
#print(ptypes)
tot<-rmultinom(size=nmolecules,n=1,prob=ptypes)
return(tot)
}

restot<-sapply(1:nrow(mydata),function(xpoint) prop2nmolecules(xpoint)) %>% t
restot<-as.data.frame(cbind(mydata$time,restot))
if (ntypes==3) { names(restot)<-c("time","y1","y2","y3") } else if (ntypes==4) { names(restot)<-c("time","y1","y2","y3","y4") }


#################################################################################
################### SAVE THE NUMBER OF MOLECULES FOR EACH TIME POINT ############
#################################################################################
write.table(restot,file=output.file,quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")
