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

#-->Logistic decay -> discard and change to make decay constant as Daniela suggests and Brinkam does
#parms_induction<-c(K=3.5,x0=1,x2=0.7,r0=0,r2=0.02)
#induction_curve<-function(x,K,x0,x2,r0,r2) K/(1+exp(-r0*(x-x0)))/(1+exp(+r2*(x-x2)))
#induction_curve_vectorized<-function(x,params) 
#  {
#  K<-params[1];x0<-params[2];x2<-params[3];r0<-params[4];r2<-params[5];
#  return(K/(1+exp(-r0*(x-x0)))/(1+exp(+r2*(x-x2))))
#  }
#x<-seq(0,2,0.01)
#plot(x,induction_curve(x,parms_induction[1],parms_induction[2],parms_induction[3],parms_induction[4],parms_induction[5]),type="l")
#-->Half-life decay as in Brinkman. r2 is 1/half life
#load("RData/notebook_2errormatrices.RData")
parms_induction<-c(K=1,x0=1,r0=10,r2=1)
#induction_curve_logistic<-function(x,K,x0,x2,r0,r2) K/(1+exp(-r0*(x-x0)))/(1+exp(+r2*(x-x2)))
induction_curve<-function(x,K,x0,r0,r2) K*2^(-x*r2)/(1+exp(-r0*(x-x0)))
induction_curve_vectorized<-function(x,params) 
  {
  return(induction_curve(x,params[1],params[2],params[3],params[4]))
  }

x<-seq(0,72,0.01)

model5i1 <- function(t, y, parms,induction_c=induction_curve) {
  with(as.list(c(y, parms)), 
       {
        dy1 <- -(k11*induction_c(t,K,x0,r0,r2) )*y[1]+r11*y[2]
        dy2 <- -(r11+r12)*y[2]+(k11*induction_c(t,K,x0,r0,r2))*y[1]
        dy3 <- r12*y[2]
        list(c(dy1, dy2, dy3))
      })
}
model5i1_nor11 <- function(t, y, parms,induction_c=induction_curve) {
  with(as.list(c(y, parms)), 
       {
        dy1 <- -(k11*induction_c(t,K,x0,r0,r2) )*y[1]
        dy2 <- -r12*y[2]+(k11*induction_c(t,K,x0,r0,r2))*y[1]
        dy3 <- r12*y[2]
        list(c(dy1, dy2, dy3))
      })
}

modelDSBs1i1_impfromcut <- function(t, y, parms,induction_c=induction_curve) {
  with(as.list(c(y, parms)), 
       {
        dy1 <- -((k11+k12)*induction_c(t,K,x0,r0,r2) )*y[1]+r11*y[3]+r21*y[4]
        dy2 <- r12*y[3]+r22*y[4]
        dy3 <- -(r11+r12)*y[3]+(k11*induction_c(t,K,x0,r0,r2))*y[1]
        dy4 <- -(r21+r22)*y[4]+(k12*induction_c(t,K,x0,r0,r2))*y[1]
        list(c(dy1, dy2, dy3, dy4))
      })
}

modelDSBs1i1_impfrompDSB <- function(t, y, parms,induction_c=induction_curve) {
  with(as.list(c(y, parms)), 
       {
        dy1 <- -(k11*induction_c(t,K,x0,r0,r2) )*y[1]+r11*y[3]+r21*y[4]
        dy2 <- r12*y[3]+r22*y[4]
        dy3 <- -(r11+r12+k12)*y[3]+(k11*induction_c(t,K,x0,r0,r2))*y[1]
        dy4 <- -(r21+r22)*y[4]+k12*y[3]
        list(c(dy1, dy2, dy3, dy4))
      })
}

modelDSBs1i1_imp2DSB <- function(t, y, parms,induction_c=induction_curve) {
  with(as.list(c(y, parms)), 
       {
        dy1 <- -(k11*induction_c(t,K,x0,r0,r2) )*y[1]+r11*y[3]+r21*y[4]
        dy2 <- r12*y[3]+r22*y[4]
        dy3 <- -(r11+r12)*y[3]+k12*y[4]
        dy4 <- -(r21+r22+k12)*y[4]+(k11*induction_c(t,K,x0,r0,r2))*y[1]
        list(c(dy1, dy2, dy3, dy4))
      })
}

modelDSBs1i1_realimprecise <- function(t, y, parms,induction_c=induction_curve) {
  with(as.list(c(y, parms)), 
       {
        dy1 <- -((k11+k12)*induction_c(t,K,x0,r0,r2) )*y[1]+r11*y[3]+r21*y[4]
        dy2 <- r12*y[3]+r22*y[4]
        dy3 <- -(r11+r12+rr12)*y[3]+(k11*induction_c(t,K,x0,r0,r2))*y[1]
        dy4 <- -(r21+r22)*y[4]+(k12*induction_c(t,K,x0,r0,r2))*y[1]+rr12*y[3]
        list(c(dy1, dy2, dy3, dy4))
      })
}

modelDSBs1i1_fullimpreciseDSB <- function(t, y, parms,induction_c=induction_curve) {
  with(as.list(c(y, parms)), 
       {
        dy1 <- -((k11+k12)*induction_c(t,K,x0,r0,r2) )*y[1]+r11*y[3]+r21*y[4]
        dy2 <- r12*y[3]+r22*y[4]
        dy3 <- -(r11+r12+rr12)*y[3]+(k11*induction_c(t,K,x0,r0,r2))*y[1]+rr21*y[4]
        dy4 <- -(r21+r22+rr21)*y[4]+(k12*induction_c(t,K,x0,r0,r2))*y[1]+rr12*y[3]
        list(c(dy1, dy2, dy3, dy4))
      })
}

modelDSBs1i1_fullimpreciseDSB_nor11 <- function(t, y, parms,induction_c=induction_curve) {
  with(as.list(c(y, parms)), 
       {
        dy1 <- -((k11+k12)*induction_c(t,K,x0,r0,r2) )*y[1]
        dy2 <- r12*y[3]+r22*y[4]
        dy3 <- -(r12+rr12)*y[3]+(k11*induction_c(t,K,x0,r0,r2))*y[1]+rr21*y[4]
        dy4 <- -(r22+rr21)*y[4]+(k12*induction_c(t,K,x0,r0,r2))*y[1]+rr12*y[3]
        list(c(dy1, dy2, dy3, dy4))
      })
}

modelDSBs1i1_realimpnor11 <- function(t, y, parms,induction_c=induction_curve) {
  with(as.list(c(y, parms)), 
       {
        dy1 <- -(k11+k12)*induction_c(t,K,x0,r0,r2)*y[1]
        dy2 <- r12*y[3]+r22*y[4]
        dy3 <- -(r12+rr12)*y[3]+(k11*induction_c(t,K,x0,r0,r2))*y[1]
        dy4 <- -r22*y[4]+(k12*induction_c(t,K,x0,r0,r2))*y[1]+rr12*y[3]
        list(c(dy1, dy2, dy3, dy4))
      })
}



###################################################################################
################### DEFINE FUNCTIONS ##############################################
###################################################################################

predict_models <- function(df,nparms=6,ntypes=6,nheaders=3,errormatrix=error_matrices3_l,mymodel=0) {
dfl<-unlist(df)
res_parms<-as.numeric(dfl[(nheaders+1):(nheaders+nparms)]);
names(res_parms)<-names(dfl[(nheaders+1):(nheaders+nparms)])
times<-seq(0,72,0.1)
yini<-c(y1 = 1, y2 = 0, y3 = 0,y4 = 0, y5 = 0, y6 = 0)
if (ntypes==3) { yini<-c(y1 = 1, y2 = 0, y3 = 0) }
if (ntypes==4) { yini<-c(y1 = 1, y2 = 0, y3 = 0, y4 = 0) }
if (is(xmodel)[1]!="function"){
mydata.fitted <- ode (times = times, y = yini, func = get(as.character(df[["model"]])[1]), parms = res_parms)} else {
mydata.fitted <- ode (times = times, y = yini, func = mymodel, parms = res_parms)
}
mydata.fitted.df<-as.data.frame(mydata.fitted)
mydata.fitted.df[,2:ncol(mydata.fitted.df)]<-t(sapply(1:nrow(mydata.fitted.df), function(x) as.matrix(mydata.fitted.df)[x,2:ncol(mydata.fitted.df)] %*% errormatrix))
if (ntypes==6) { names(mydata.fitted.df)<-c("time","x0","x+","x-","y0","y+","y-") }
if (ntypes==3) { names(mydata.fitted.df)<-c("time","intact","DSB","indels") } 
if (ntypes==4) { names(mydata.fitted.df)<-c("time","intact","indels","preciseDSB","impreciseDSB") } 
tidy.mydata.fitted.df<-mydata.fitted.df %>% pivot_longer(names(mydata.fitted.df)[2:(ntypes+1)],names_to = "types", values_to = "p")
return(tidy.mydata.fitted.df)
}


predict_induction <- function(df,nparms=6,nparms_induction=5,nheaders=3,induction_function=induction_curve_vectorized) {
#nparms=3;nparms_induction=5;ntypes=6;induction_function=induction_curve_vectorized
#df<-bestmodels_l[1,]
dfl<-unlist(df)
res_parms<-as.numeric(dfl[(nheaders+nparms+1):(nheaders+nparms+nparms_induction)]);
names(res_parms)<-names(dfl[(nheaders+nparms+1):(nheaders+nparms+nparms_induction)])
times<-seq(0,72,0.1)
mydata.fitted <-data.frame(time=times,curve_fitted=induction_function(times,res_parms))
return(mydata.fitted)
}

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

argv <- add_argument(argv, "-T", help="time course of target site. A dataset with time as 1st column and then the number of molecules")
argv <- add_argument(argv, "-m", help="model")
argv <- add_argument(argv, "-o", help="output file", default="output.txt")
argv <- add_argument(argv, "-e", help="errors")
argv <- add_argument(argv, "-E", help="error matrix")
argv <- add_argument(argv, "-n", help="n iterations", default=100)
argv <- add_argument(argv, "-d", help="switch to change likelihood function. To estimate a common error from DSB select 1. Default is no estimate from data, only from controls, to sample 0.2h after induction (0). To set induction curve without delay select 2.", default=0)



args <- parse_args(argv)
input.file<-args$T
myerror<-args$e
myerrorE<-args$E
mymodel<-args$m
output.file<-args$o
n.max<-as.numeric(args$n)
optimize_errorDSB2indel<-as.numeric(args$d)

if (is.na(myerror)){ myerror<-"error" }
if (is.na(n.max)){ n.max<-100 }
if (is.na(output.file)){ output.file<-"DSBtimecourse_optimize.tsv" }
if (is.na(optimize_errorDSB2indel)){ optimize_errorDSB2indel<-0 }



########################################################################################################
##################### LOAD DATA ########################################################################
########################################################################################################
mydelay<-0
ntypes<-4

mydata<-read.table(input.file, header=TRUE)
errormatrix<-as.matrix(read.table(myerrorE, header=FALSE))
if (mymodel=="model5i1" || mymodel=="model5i1_nor11")
	{ 
	ntypes<-3 ; load("RData/error_matrices3_l.RData");
	error_mat_l<-error_matrices3_l
 	}

if (optimize_errorDSB2indel==2) { mydelay<-0 } 


#################################################################################
##################### START OPTIMIZATION ########################################
#################################################################################
print("prepare optimization")


#if (mytarget_i=="CRTISO_49and50bp") {mytarget<-"CRTISO"} else {mytarget<-mytarget_i};
for ( counter in 1:n.max)
	{
	if (mymodel=="modelDSBs1i1_fullimpreciseDSB")
		{
    		nameparms <-c("k11","k12","rr12","rr21","r11","r12","r21","r22","K","x0","r0","r2")
  		} else 
	if ( mymodel=="modelDSBs1i1_impfromcut" || mymodel=="modelDSBs1i1_impfrompDSB" || mymodel=="modelDSBs1i1_imp2DSB")
		{
		nameparms <-c("k11","k12","r11","r12","r21","r22","K","x0","r0","r2")
		} else
	if ( mymodel=="modelDSBs1i1_realimprecise")
		{
	        nameparms <-c("k11","k12","rr12","r11","r12","r21","r22","K","x0","r0","r2")
		} else
	if ( mymodel=="modelDSBs1i1_realimpnor11")
		{
		nameparms <-c("k11","k12","rr12","r12","r22","K","x0","r0","r2")
		} else
	if ( mymodel=="modelDSBs1i1_fullimpreciseDSB_nor11")
		{
    		nameparms <-c("k11","k12","rr12","rr21","r12","r22","K","x0","r0","r2")
		} else
	if ( mymodel=="model5i1")
		{
    		nameparms <-c("k11","r11","r12","K","x0","r0","r2")
		} else
	if ( mymodel=="model5i1_nor11")
		{
    		nameparms <-c("k11","r12","K","x0","r0","r2")
		};
	if ( optimize_errorDSB2indel==1 ) 
	{ 
	loglik_er_f.pen<-loglik_er_f.pen_errorDBS2indel_4states_m1
	nameparms<-c(nameparms,"er1")
	} else if ( optimize_errorDSB2indel==2 )
	{
	loglik_er_f.pen<-loglik_er_f.nopen
	};
	xmodel=get(mymodel);
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
        #mydata<-mydata_l[[paste0(myinduction,"_",mytarget_i)]]
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

