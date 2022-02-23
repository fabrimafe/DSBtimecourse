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

source("likelihood_functions.R")

###################################################################################
################### DEFINE FUNCTIONS ##############################################
###################################################################################

#predict_models <- function(df,nparms=6,ntypes=6,nheaders=3,errormatrix=error_matrices3_l,mymodel=0) {
#dfl<-unlist(df)
#res_parms<-as.numeric(dfl[(nheaders+1):(nheaders+nparms)]);
#names(res_parms)<-names(dfl[(nheaders+1):(nheaders+nparms)])
#times<-seq(0,72,0.1)
#yini<-c(y1 = 1, y2 = 0, y3 = 0,y4 = 0, y5 = 0, y6 = 0)
#if (ntypes==3) { yini<-c(y1 = 1, y2 = 0, y3 = 0) }
#if (ntypes==4) { yini<-c(y1 = 1, y2 = 0, y3 = 0, y4 = 0) }
#if (is(xmodel)[1]!="function"){
#mydata.fitted <- ode (times = times, y = yini, func = get(as.character(df[["model"]])[1]), parms = res_parms)} else {
#mydata.fitted <- ode (times = times, y = yini, func = mymodel, parms = res_parms)
#}
#mydata.fitted.df<-as.data.frame(mydata.fitted)
#mydata.fitted.df[,2:ncol(mydata.fitted.df)]<-t(sapply(1:nrow(mydata.fitted.df), function(x) as.matrix(mydata.fitted.df)[x,2:ncol(mydata.fitted.df)] %*% errormatrix))
#if (ntypes==6) { names(mydata.fitted.df)<-c("time","x0","x+","x-","y0","y+","y-") }
#if (ntypes==3) { names(mydata.fitted.df)<-c("time","intact","DSB","indels") } 
#if (ntypes==4) { names(mydata.fitted.df)<-c("time","intact","indels","preciseDSB","impreciseDSB") } 
#tidy.mydata.fitted.df<-mydata.fitted.df %>% pivot_longer(names(mydata.fitted.df)[2:(ntypes+1)],names_to = "types", values_to = "p")
#return(tidy.mydata.fitted.df)
#}

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


################################################################################
################## READ ARGUMENTS ##############################################
################################################################################


argv<- arg_parser("Parse arguments")

argv <- add_argument(argv, "-T", help="time course of target site. A dataset with time as 1st column and then the number of molecules")
argv <- add_argument(argv, "-m", help="model")
argv <- add_argument(argv, "-o", help="output file", default="output.txt")
argv <- add_argument(argv, "-e", help="errors")
argv <- add_argument(argv, "-E", help="error matrix")
argv <- add_argument(argv, "-n", help="n iterations", default=100)
argv <- add_argument(argv, "-z", help="n parameters in induction curve", default=2)
argv <- add_argument(argv, "-l", help="switch to change likelihood function. To estimate a common error from DSB select 1. Default is no estimate from data, only from controls, to sample 0.2h after induction (0). To set induction curve without delay select 2. To model imprecise DSB as misread precise DSB select 3", default=0)



args <- parse_args(argv)
input.file<-args$T
#myerror<-args$e
myerrorE<-args$E
mymodel<-args$m
nparamsind<-as.numeric(args$z)
output.file<-args$o
n.max<-as.numeric(args$n)
optimize_errorDSB2indel<-as.numeric(args$l)

#if (is.na(myerror)){ myerror<-"error" }
if (is.na(n.max)){ n.max<-100 }
if (is.na(output.file)){ output.file<-"DSBtimecourse_optimize.tsv" }
if (is.na(optimize_errorDSB2indel)){ optimize_errorDSB2indel<-0 }
if (is.na(nparamsind)){ nparamsind<-2 }
define_ODE.functions(nparamsind)


########################################################################################################
##################### LOAD DATA ########################################################################
########################################################################################################
mydelay<-0
ntypes<-4

mydata<-read.table(input.file, header=TRUE)
time_courses_begins<-c(which(mydata$time[-1]-mydata$time[-length(mydata$time)]<0)+1)
nameparms<-model2nameparams(mymodel,nparamsind)
errormatrix<-as.matrix(read.table(myerrorE, header=FALSE))
xmodel<-NA

if (mymodel=="model5i1" || mymodel=="model5i1_nor11")
                {
                ntypes<-3
                };
if ( mymodel=="modelDSBs1i1_realimprecise.inductionx3")
                {
                loglik_er_f.pen<-loglik_er_f.pen_modelinductionx3
                xmodel<-get("modelDSBs1i1_realimprecise")
	        #nameparms <-c("k11","k11","k11","k12","k12","k12","rr12","rr12","rr12","r11","r11","r11","r12","r12","r12","r21","r21","r21","r22","r22","r22","K","x0","r0","r2")
                };
if ( optimize_errorDSB2indel==1 ) 
	{ 
	loglik_er_f.pen<-loglik_er_f.pen_errorDBS2indel_4states_m1
	nameparms<-c(nameparms,"er1")
	} else 
if ( optimize_errorDSB2indel==2 )
	{
	loglik_er_f.pen<-loglik_er_f.nopen
	} else 
if ( optimize_errorDSB2indel==3 )
	{
	loglik_er_f.pen<-loglik_er_f.pen_errorimpDSB2pDSB_4states_m1
	nameparms<-c(nameparms,"er1")
	};

if (optimize_errorDSB2indel==2) { mydelay<-0 } 

if (!"function" %in% is(xmodel)) { xmodel<-get(mymodel) }
print(xmodel)
#################################################################################
##################### START OPTIMIZATION ########################################
#################################################################################
print("prepare optimization")

#if (mytarget_i=="CRTISO_49and50bp") {mytarget<-"CRTISO"} else {mytarget<-mytarget_i};
for ( counter in 1:n.max)
	{
        nrates<-length(nameparms)
	step_exp<-sample(6:8,1);mysteps<-rep(10^(-step_exp),nrates)
	parscales_exp<-sample(1:3,nrates,replace=TRUE);myparscale<-10^(-parscales_exp)
	if (sample(0:1,1)==1) { myparscale<-rep(1,nrates) };
	xparms<-10^runif(nrates,-6,1)
        lower_bound<-c(rep(0,nrates-2),10^(-7),0)
	mylmm<-sample(5:8,1)
	myfactr<-10^(-sample(c(7,8,10),1))
#	print(c(counter,mytarget_i,myinduction,myerror))
	names(xparms)<-nameparms
	res<-optimx(par=xparms,fn=function(z)
        loglik_er_f.pen(z,my_data=mydata,ODEfunc=xmodel,E.matrix=errormatrix),control=list(maximize=TRUE,maxit=1000,parscale=myparscale,ndeps=mysteps,lmm=mylmm,factr=myfactr),method="L-BFGS-B",lower=lower_bound,hessian=FALSE,gr=NULL) #dowarn=FALSE
	names(myparscale)<-paste0("parscale",1:nrates)
        names(xparms)<-paste0("init_",names(xparms))
        names(mysteps)<-paste0("step",1:nrates)
        res_allinfo<-t(unlist(list(res,xparms,myparscale,mysteps)))
        if (counter==1) { print_header=TRUE } else { print_header=FALSE }
        write.table(res_allinfo,file=output.file,quote=FALSE,row.names=FALSE,append=TRUE,col.names=print_header,sep="\t")
        };
print("done") 

