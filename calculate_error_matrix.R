#!/usr/bin/env Rscript
#example to run
#./optimize_model_backbone.R  -t CRTISO -i RNP -m modelDSBs1i1_fullimpreciseDSB -e E_errorsfromunbroken -n 10000 -o CRTISO_RNP_optimization.tsv 

library(argparser, quietly=TRUE,warn.conflicts=FALSE)
library(tidyverse)
library(optimx)
library(ggplot2)
library(deSolve)

################################################################
# DRAW THE MODEL (an ODE function) AND GENERATE SOME RANDOM DATA
################################################################
parms_induction<-c(K=1,x0=1,r0=10,r2=1)
#induction_curve_logistic<-function(x,K,x0,x2,r0,r2) K/(1+exp(-r0*(x-x0)))/(1+exp(+r2*(x-x2)))
induction_curve<-function(x,K,x0,r0,r2) K*2^(-x*r2)/(1+exp(-r0*(x-x0)))
induction_curve_vectorized<-function(x,params) 
  {
  return(induction_curve(x,params[1],params[2],params[3],params[4]))
  }

###################################################################################
################### DEFINE FUNCTIONS ##############################################
###################################################################################


loglik_er_f<-function(parms,my_data=mydata,ODEfunc=model1,E.matrix=error_matrix){
  #ODEfunc: a function of class desolve 
  #my_data: a dataframe with "time" course (which should always include 0) as 1st column, and col-types compatible with ODEfunc
  #parms: vector specifying the set of parameters (compatible with ODEfunc)
  #ncells<-apply(mydata[,-1],MARGIN=1,sum)
  yinit<-c(1,rep(0,ncol(my_data)-2))
  out <- ode (times = my_data$time, y = yinit, func = ODEfunc, parms = parms)
  probs<-(out[,2:ncol(out)]) %*% E.matrix
  probs[probs<10^(-16)]<-10^(-16)
  loglik<-sapply(1:nrow(my_data), function(x) dmultinom(x=my_data[x,-1],prob=probs[x,],log=TRUE))
  sum(loglik)
}

loglik_er_f.pen<-function(parms,my_data=mydata,ODEfunc=model1,E.matrix=error_matrix,induction_curve=induction_curve_vectorized){
  penalty<-induction_curve_vectorized(0,parms[(length(parms)-3):length(parms)])
  if (penalty>0.00001){ penalty<-(10^7)*(penalty-0.00001)^2 } else {penalty<-0}
  loglik_er_f(parms,my_data,ODEfunc,E.matrix)-penalty
  }

loglik_er_f.nopen<-function(parms,my_data=mydata,ODEfunc=model1,E.matrix=error_matrix,induction_curve=induction_curve_vectorized){
  #penalty<-induction_curve_vectorized(0,parms[(length(parms)-3):length(parms)])
  #if (penalty>0.00001){ penalty<-(10^7)*(penalty-0.00001)^2 } else {penalty<-0}
  loglik_er_f(parms,my_data,ODEfunc,E.matrix) #-penalty
  }

loglik_er_f.pen_errorDBS2indel_4states_m1<-function(parms,my_data=mydata,ODEfunc=model1,E.matrix=error_matrix,induction_curve=induction_curve_vectorized){
  erDSB2indel<-parms[length(parms)]
  penalty<-0
  if (erDSB2indel>0.5) { penalty<-100000*(erDSB2indel-0.7)^2 }
  E.matrix_t<-E.matrix
  E.matrix_t[3,1]<-erDSB2indel/2
  E.matrix_t[4,1]<-erDSB2indel/2
  E.matrix_t[3,2]<-erDSB2indel/2
  E.matrix_t[3,2]<-erDSB2indel/2
  E.matrix_t[3,3]<-1-erDSB2indel
  E.matrix_t[4,4]<-1-erDSB2indel
  penalty<-penalty+induction_curve_vectorized(0,parms[(length(parms)-3):length(parms)])
  if (penalty>0.00001){ penalty<-(10^7)*(penalty-0.00001)^2 } else {penalty<-0}
  loglik_er_f(parms,my_data,ODEfunc,E.matrix)-penalty
  }

#function for unconstrained optimization
loglik_er_f.pen.unconstrainedopt<-function(parms,my_data=mydata,ODEfunc=model1,E.matrix=error_matrix,induction_curve=induction_curve_vectorized){
  penalty<-induction_curve_vectorized(0,parms[(length(parms)-3):length(parms)])
  if (penalty>0.00001){ penalty<-(10^7)*(penalty-0.00001)^2 } else {penalty<-0}
  penalty<-penalty+sum(apply(cbind(parms,rep(0,length(parms))),MARGIN=1,min)^2)*10^16
  loglik_er_f(parms,my_data,ODEfunc,E.matrix)-penalty
  }

################################################################################
################## READ ARGUMENTS ##############################################
################################################################################


argv<- arg_parser("Parse arguments")

argv <- add_argument(argv, "-i", help="input data")
argv <- add_argument(argv, "-o", help="output file root", default="error_matrix")
argv <- add_argument(argv, "-e", help="error model")
argv <- add_argument(argv, "-d", help="number of types -> dimension matrix: 3,4", default=4)
argv <- add_argument(argv, "-f", help="input format: 1 is target,time,y1,y2...; 2 is csv with MMEJ from daniela")

args <- parse_args(argv)
myerror<-args$e
input.file<-args$i
output.file<-args$o
ntypes<-args$d
input.format<-args$f

if (is.na(myerror)){ myerror<-"E_noerrors" }
if (is.na(ntypes)){ ntypes<-4 }
if (is.na(input.format)){ ntypes<-1 }
if (is.na(output.file)){ output.file<-"error_matrix" }

########################################################################################################
##################### LOAD DATA ########################################################################
########################################################################################################

#print(mydata)
#################################################################################
##################### START OPTIMIZATION ########################################
#################################################################################

#=====================================4 types=
if (ntypes==4){ 
error_matrices4_l<-list()
if (input.format==2)
	{
	mydata<-read.table(input.file, header=TRUE,sep=",")
	mydata[is.na(mydata)]<-0;
	mydata<-cbind(matrix(unlist(sapply(as.character(mydata$X),function(x) strsplit(x,"_"))),ncol=4,byrow=T),mydata);mydata[,2]<-NULL
	names(mydata)[1:4]<-c("target","time","replicate","ID_name"); mydata$time<-as.numeric(sapply(mydata$time, function(x) gsub("h","",x)))

	mydata<-mydata %>% mutate(preciseDSB=Perfect.DSB,impreciseDSB=Extended.DSB+PAM_side.DSB+guide_side.DSB,indels=Ins+Del+MH_Del+MH_Del_2)
	mydata<- mydata %>% select(time,WT.Sub,indels,preciseDSB,impreciseDSB)
	names(mydata)<-c("time","y1","y2","y3","y4")
	} else if (input.format==1)
	{
	mydata<-read.table(input.file, header=TRUE,sep="\t")
	mydata[is.na(mydata)]<-0;
	}
print(mydata)
mydata.p<-as.data.frame(cbind(mydata$time,t(apply(mydata[,2:5],MARGIN=1,FUN=function(x) x/sum(x)))))
names(mydata.p)<-c("time","y1","y2","y3","y4")

#cor.test(apply(mydata.p[,3:4],MARGIN=1,sum),apply(mydata[,2:4],MARGIN=1,sum))
loglik_errorcontrol_f<-function(parms,mydata=mydata){
  #mydata: a dataframe with "time" course (which should always include 0) as 1st column, and col-types compatible with ODEfunc
  #parms: vector specifying the set of parameters (compatible with ODEfunc)
  loglik<-sapply(1:nrow(mydata), function(y) dmultinom(x=mydata[y,-1],prob=parms,log=TRUE))
  return(sum(loglik))
}


xparms<-c(p11=0.01,p12=0.01,p13=0.01,p14=0.01)
sumerror<-rep(1/(ncol(mydata)-2),ncol(mydata)-2)/10
error_proxy<-optimx(par=xparms,fn=function(z) loglik_errorcontrol_f(z,mydata=mydata),control=list(maximize=TRUE),method="L-BFGS-B",lower=rep(10^(-9),ncol(mydata)-1),hessian=FALSE,gr=NULL)

error_rates<-error_proxy[1:4]/sum(error_proxy[1:4]); error_rates[1]<-1-sum(error_rates[-1]); error_rates<-unlist(error_rates)

xerror_rates<-error_rates
xprobs<-c(1,0,0,0)
#GENERATE ERROR MATRICES
n_types<-ncol(mydata)-1
#--no errors
E_noerrors<-diag(n_types)
E_noerrors<-(E_noerrors+10^(-6))/(1+n_types*(10^(-6)))
#--errors only from intact molecules
E_errorsfromintact<-unname(rbind(error_rates,diag(n_types)));E_errorsfromintact<-E_errorsfromintact[-2,]
#--errors only and equally from unbroken molecules
error_matrix<-E_errorsfromintact
for ( i in 2:(n_types/2)){
error_matrix[i,]<-error_rates
error_matrix[i,1]<-0
error_matrix[i,i]<-error_rates[1]+error_rates[i]
}
for ( i in (n_types/2+1):n_types){
error_matrix[i,]<-rep(0,n_types)
error_matrix[i,i]<-1
}
E_errorsfromunbroken<-error_matrix

#check that all sum to 1
xprobs<-rep(1/4,4)
sum(xprobs %*% E_noerrors)
sum(xprobs %*% E_errorsfromintact)
sum(xprobs %*% E_errorsfromunbroken)

error_matrices4_l[["E_noerrors"]]<-E_noerrors
error_matrices4_l[["E_errorsfromintact"]]<-E_errorsfromintact
error_matrices4_l[["E_errorsfromunbroken"]]<-E_errorsfromunbroken



write.table(error_matrices4_l[["E_noerrors"]],file=paste0(output.file,"noerrors.tsv"),quote=FALSE,row.names=FALSE,sep="\t",col.names=FALSE)
write.table(error_matrices4_l[["E_errorsfromintact"]],file=paste0(output.file,"errorsfromintact.tsv"),quote=FALSE,row.names=FALSE,sep="\t",col.names=FALSE)
write.table(error_matrices4_l[["E_errorsfromunbroken"]],file=paste0(output.file,"errorsfromunbroken.tsv"),quote=FALSE,row.names=FALSE,sep="\t",col.names=FALSE)

} else if (ntypes==3) #########################################=3 types=################################3
{

error_matrices3_l<-list()
if (input.format==2)
	{
	mydata<-read.table(input.file, header=TRUE,sep=",")
	mydata[is.na(mydata)]<-0;
	mydata<-cbind(matrix(unlist(sapply(as.character(mydata$X),function(x) strsplit(x,"_"))),ncol=4,byrow=T),mydata);mydata[,2]<-NULL
	names(mydata)[1:4]<-c("target","time","replicate","ID_name"); mydata$time<-as.numeric(sapply(mydata$time, function(x) gsub("h","",x)))

	mydata<-mydata %>% mutate(DSB=Perfect.DSB+Extended.DSB+PAM_side.DSB+guide_side.DSB,indels=Ins+Del+MH_Del+MH_Del_2)
	mydata<- mydata %>% select(time,WT.Sub,DSB,indels)
	names(mydata)<-c("time","y1","y2","y3")
	} else if (input.format==1)
	{
	mydata<-read.table(input.file, header=TRUE,sep="\t")
	mydata[is.na(mydata)]<-0;
	}
print(mydata)
mydata.p<-as.data.frame(cbind(mydata$time,t(apply(mydata[,2:4],MARGIN=1,FUN=function(x) x/sum(x)))))
names(mydata.p)<-c("time","y1","y2","y3")

loglik_errorcontrol_f<-function(parms,mydata=mydata){
  #mydata: a dataframe with "time" course (which should always include 0) as 1st column, and col-types compatible with ODEfunc
  #parms: vector specifying the set of parameters (compatible with ODEfunc)
  loglik<-sapply(1:nrow(mydata), function(y) dmultinom(x=mydata[y,-1],prob=parms,log=TRUE))
  return(sum(loglik))
}

xparms<-c(p11=0.01,p12=0.01,p13=0.01)
sumerror<-rep(1/(ncol(mydata)-2),ncol(mydata)-2)/10
error_proxy<-optimx(par=xparms,fn=function(z) loglik_errorcontrol_f(z,mydata=mydata),control=list(maximize=TRUE),method="L-BFGS-B",lower=rep(10^(-9),ncol(mydata)-1),hessian=FALSE,gr=NULL)

error_rates<-error_proxy[1:3]/sum(error_proxy[1:3]); error_rates[1]<-1-sum(error_rates[-1]); error_rates<-unlist(error_rates)

xerror_rates<-error_rates
xprobs<-c(1,0,0)
#GENERATE ERROR MATRICES
n_types<-ncol(mydata)-1
#--no errors
E_noerrors<-diag(n_types)
E_noerrors<-(E_noerrors+10^(-6))/(1+n_types*(10^(-6)))
#--errors only from intact molecules
E_errorsfromintact<-unname(rbind(error_rates,diag(n_types)));E_errorsfromintact<-E_errorsfromintact[-2,]
#--errors only and equally from unbroken molecules
error_matrix<-E_errorsfromintact
for ( i in 2:(n_types/2)){
error_matrix[i,]<-error_rates
error_matrix[i,1]<-0
error_matrix[i,i]<-error_rates[1]+error_rates[i]
}
for ( i in (n_types/2+1):n_types){
error_matrix[i,]<-rep(0,n_types)
error_matrix[i,i]<-1
}
E_errorsfromunbroken<-error_matrix

#check that all sum to 1
xprobs<-rep(1/6,3)
sum(xprobs %*% E_noerrors)
sum(xprobs %*% E_errorsfromintact)
sum(xprobs %*% E_errorsfromunbroken)

error_matrices3_l[["E_noerrors"]]<-E_noerrors
error_matrices3_l[["E_errorsfromintact"]]<-E_errorsfromintact
error_matrices3_l[["E_errorsfromunbroken"]]<-E_errorsfromunbroken

write.table(error_matrices3_l[["E_noerrors"]],file=paste0(output.file,"noerrors.tsv"),quote=FALSE,row.names=FALSE,sep="\t",col.names=FALSE)
write.table(error_matrices3_l[["E_errorsfromintact"]],file=paste0(output.file,"errorsfromintact.tsv"),quote=FALSE,row.names=FALSE,sep="\t",col.names=FALSE)
write.table(error_matrices3_l[["E_errorsfromunbroken"]],file=paste0(output.file,"errorsfromunbroken.tsv"),quote=FALSE,row.names=FALSE,sep="\t",col.names=FALSE)

}


